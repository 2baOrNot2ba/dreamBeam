'''
This module provides some common observational scenarios of interest.
It implements basic Jones chains or Measurement Equations.
'''
import numpy as np
import matplotlib.pyplot as plt
import dreambeam.rime.jones
from dreambeam.rime.conversionUtils import crt2sph


def on_pointing_axis_tracking(telescope, stnID, ObsTimeBeg, duration,
                              ObsTimeStp, CelDir, do_parallactic_rot=True,
                              xtra_results=False):
    """Computes the Jones matrix along pointing axis while tracking a fixed
    celestial source. """  # # FIXME: Doesn't use freq
    #    *Setup Source*
    #celAz, celEl, celRef = CelDir.split(',')
    #celAz = float(celAz)
    #celEl = float(celEl)
    (celAz, celEl, celRef) = CelDir
    srcfld = dreambeam.rime.jones.DualPolFieldPointSrc((celAz, celEl, celRef))

    stnBD = telescope['Station'][stnID]
    stnRot = stnBD.stnRot

    #    *Setup PJones*
    #duration = ObsTimeEnd-ObsTimeBeg
    timespy = []
    nrTimSamps = int((duration.total_seconds()/ObsTimeStp.seconds))+1
    for ti in range(0, nrTimSamps):
        timespy.append(ObsTimeBeg+ti*ObsTimeStp)
    pjones = dreambeam.rime.jones.PJones(timespy, np.transpose(stnRot),
                                         do_parallactic_rot=do_parallactic_rot)

    #    *Setup EJones*
    ejones = stnBD.getEJones(CelDir)
    stnDPolel = stnBD.feed_pat

    #    *Setup MEq*
    pjonesOfSrc = pjones.op(srcfld)
    res = ejones.op(pjonesOfSrc)

    #Get the resulting Jones matrices
    #(structure is Jn[freqIdx, timeIdx, chanIdx, compIdx] )
    Jn = res.getValue()
    freqs = stnDPolel.getfreqs()
    if xtra_results:
        return timespy, freqs, Jn, srcfld, res, pjonesOfSrc
    else:
        return timespy, freqs, Jn


def beamfov(telescope, stnID, ObsTime, CelDir, freq):
    """Computes the Jones matrix over the beam fov for pointing.
    """
    #    *Setup Source*
    (celAz, celEl, celRef) = CelDir
    srcfld = dreambeam.rime.jones.DualPolFieldRegion()

    stnBD = telescope['Station'][stnID]
    stnRot = stnBD.stnRot

    #    *Setup Parallatic Jones*
    pjones = dreambeam.rime.jones.PJones([ObsTime], np.transpose(stnRot))

    #    *Setup EJones*
    stnDPolel = stnBD.feed_pat
    #    **Select frequency
    freqs = stnDPolel.getfreqs()
    frqIdx = np.where(np.isclose(freqs,freq,atol=190e3))[0][0]
    ejones = stnBD.getEJones(CelDir, [freqs[frqIdx]]) #Ejones doesnt use CelDir

    #    *Setup MEq*
    pjonesOfSrc = pjones.op(srcfld)
    res = ejones.op(pjonesOfSrc)

    #Get the resulting Jones matrices
    #(structure is Jn[freqIdx, timeIdx, chanIdx, compIdx] )
    Jn = res.getValue()
    #compute_paral(srcfld, stnRot, res, pjonesOfSrc)
    #res.get_basis()
    return srcfld.azmsh, srcfld.elmsh, Jn, ejones.thisjones


def compute_paral(srcfld, stnRot, res, pjonesOfSrc, ObsTimeBeg):
    """Compute parallactic rotation. Also displays pointings in horizontal
    coordinates as seen on sky from below."""
    def ISO2horz(az=None, aztype='NoE'):
        if aztype == 'NoE':
            azoffset = np.pi/2
            azsign = -1
        if az is None:
            return azoffset
        else:
            return azsign*az+azoffset

    def el2theta(el):
        return np.pi/2-el

    basisITRF_lcl = res.jonesbasis
    ax = plt.subplot(111, projection='polar')
    # Will use horizontal crd sys here with N over E
    ax.set_theta_offset(ISO2horz())
    nrsamps = basisITRF_lcl.shape[0]
    az = np.zeros((nrsamps))
    el = np.zeros((nrsamps))
    Jn = res.getValue()
    Jnf = Jn[256, :, :, :].squeeze()  # Midpoint freq.
    itrf_z_stn = np.matmul(stnRot.T, [[0], [0], [1]])
    itrf_z_az, itrf_z_el = crt2sph(itrf_z_stn)
    itrf_z_az = ISO2horz(itrf_z_az)
    itrf_z_tht = np.rad2deg(el2theta(itrf_z_el))

    antmeasang = []  # Antenna measured angle
    for i in range(nrsamps):
        az[i], el[i] = crt2sph(basisITRF_lcl[i, :, 0])
        az[i] = ISO2horz(az[i])
        # Compute angle of IAU x direction as seen by dualpol ants
        # (angle is from N over E)
        antmeasang.append(-(np.arctan2(
                    np.real(Jnf[i, 1, 0]), np.real(Jnf[i, 0, 0])) + 3*np.pi/4))
    showIAUx = False
    if showIAUx:
        print map(np.rad2deg, antmeasang)
        # Create vectors out of ang from origin
        angsaz = np.stack((antmeasang, np.zeros((nrsamps,))))
        angsr = np.stack((45*np.ones((nrsamps,)), np.zeros((nrsamps,))))
        # Plot IAU x vectors
        ax.plot(angsaz, angsr, '-')

    # Display pointings in horizontal coordinates
    ax.plot(az, np.rad2deg(el2theta(el)), '+', label='Trajectory')
    # Mark out start point
    ax.plot(az[0], np.rad2deg(el2theta(el[0])), 'r8',
            label='Start: '+ObsTimeBeg.isoformat()+'UT')
    ax.plot(itrf_z_az, itrf_z_tht, '*', label='NCP')
    ax.set_rmax(90)
    plt.title('Source trajectory [local coords, ARC projection]')
    ax.legend(numpoints=1, loc='best')
    plt.annotate('N (stn)', xy=(0, 90))
    plt.annotate('E (stn)', xy=(np.pi/2, 90))
    plt.draw()
