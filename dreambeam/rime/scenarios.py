'''
This module provides some common observational scenarios of interest.
It implements basic Jones chains or Measurement Equations.
'''
import numpy as np
import matplotlib.pyplot as plt
from antpat.reps.sphgridfun import pntsonsphere
import dreambeam.rime.jones


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
    #res.getBasis()
    return srcfld.azmsh, srcfld.elmsh, Jn, ejones.thisjones


def compute_paral(srcfld, stnRot, res, pjonesOfSrc, ObsTimeBeg):
    """Compute parallactic rotation. Also displays pointings in horizontal
    coordinates as seen on sky from below."""
    srcbasis = srcfld.jonesbasis
    basisITRF_lcl = res.jonesbasis
    basisJ2000_ITRF = pjonesOfSrc.jonesbasis
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_offset(np.pi/2)
    nrsamps = basisITRF_lcl.shape[0]
    az = np.zeros((nrsamps))
    el = np.zeros((nrsamps))
    Jn = res.getValue()
    Jnf = Jn[256, :, :, :].squeeze()  # Midpoint freq.
    ITRFz_stn = np.matmul(stnRot.T, [[0], [0], [1]])
    ITRFz_stnaz = np.arctan2(ITRFz_stn[0, 0], ITRFz_stn[1, 0])
    ITRFz_stntht = np.rad2deg(np.arccos(ITRFz_stn[2, 0]))
    ang = []
    for i in range(nrsamps):
        basisJ2000_ITRF_to = np.matmul(stnRot, basisJ2000_ITRF[i, :, :])
        paramat = np.matmul(basisITRF_lcl[i, :, :].T, basisJ2000_ITRF_to)
        az[i], el[i] = pntsonsphere.crt2sphHorizontal(
                                            basisITRF_lcl[i, :, 0].squeeze())
        # Compute ang of the IAU x direction (sign here is wrt N over E)
        ang.append(np.arctan2(np.real(Jnf[i, 1, 0]), np.real(Jnf[i, 0, 0]))
                   - 5*np.pi/4)
    showIAUx = False
    if showIAUx:
        # Create vectors out of ang from origin
        angsaz = np.stack((ang, np.zeros((len(ang),))))
        angsr = np.stack((45*np.ones((len(ang),)), np.zeros((len(ang),))))
        # Plot IAU x vectors
        ax.plot(angsaz, angsr, '-')

    # Display pointings in horizontal coordinates
    ax.plot(az, 90-el/np.pi*180, '+', label='Trajectory')
    # Mark out start point
    ax.plot(az[0], 90-el[0]/np.pi*180, 'r8',
            label='Start: '+ObsTimeBeg.isoformat()+'UT')
    ax.plot(ITRFz_stnaz, ITRFz_stntht, '*', label='NCP')
    ax.set_rmax(90)
    plt.title('Source trajectory [local coords, ARC projection]')
    ax.legend(numpoints=1, loc='best')
    plt.annotate('N (stn)', xy=(0, 90))
    plt.annotate('E (stn)', xy=(np.pi/2, 90))
    plt.draw()
