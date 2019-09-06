'''
This module provides some common observational scenarios of interest.
It implements basic Jones chains or Measurement Equations.
'''
import numpy as np
import matplotlib.pyplot as plt
import dreambeam.rime.jones
from dreambeam.telescopes.rt import TelescopesWiz
from dreambeam.rime.conversion_utils import basis2basis_transf, IAUtoC09, \
                                            C09toIAU


def on_pointing_axis_tracking(telescopename, stnid, band, antmodel, obstimebeg,
                              obsdur, obstimestp, pointingdir,
                              do_parallactic_rot=True):
    """Computes the Jones matrix along pointing axis while tracking a fixed
    celestial source.

    This function computes the Jones of an observational setup where a
    telescope station is tracking a source (on-axis) on the sky using one of
    its bands.

    The Jones matrix computed is the matrix which maps the two transverse
    components of the E-field vector (V/m or unitless) at a given frequency
    propagating from the pointing direction to the telescope stations two
    polarization channels. The x,y components of the source are as specified by
    the IAU for polarized emissions; and the output components are the two
    ordered polarization channels. The matrix inverse of the output Jones
    matrix multiplied from the left on the channelized 2D voltage data will
    produce an estimate of the E-field vector of the source in the IAU for all
    the frequencies in the band over the times of the tracking.

    The response of the dual-polarized feed is modeled using the `antmodel`
    specified.

    Parameters
    ----------
    telescopename : str
        Name of telescope, as registered in TelescopesWiz() instance.
    stnid : str
        Name or ID of the station, as registered in TelescopesWiz() instance.
    band : str
        Name of band, as registered in TelescopesWiz() instance.
    antmodel : str
        Name of antenna model, e.g. 'Hamaker', as registered in TelescopesWiz()
        instance.
    obstimebeg : datetime.datetime
        Date-time when the tracking observation begins.
    obsdur : datetime.deltatime
        Duration of the entire tracking observation in seconds. The sample
        at obstimebeg+duration is included.
    obstimestp : datetime.deltatime
        Time step in seconds for which the jones matrix should be sampled at.
    pointingdir : (float, float, str)
        Length 3 tuple encoding the tracking direction on the celestial sphere.
        The last tuple element should usually be 'J2000', in which case the
        the first two tuple elements are the right ascension and declination,
        respectively, in radians.
    do_parallactic_rot : bool (optional)
        Whether of not to perform parallactic rotation (default True).

    Returns
    -------
    timespy : array
        The python datetime of the samples.
    freqs : array
        The frequency [Hz] at which the Jones was computed.
    jones : array_like
        Array over time steps and frequencies of the Jones matrix corresponding
        to RIME for this set of input parameters.

    Notes
    -----
    Specifically, this function only considers an RIME consisting of the
    projection from the celestial (Earth-Centered-Inertial) frame to the
    topocentric station frame, and the polarimetric responce of the
    dual-polarized antenna feed. These are known as the P-Jones and the E-Jones
    respectively.

    Examples
    --------
    Here is an example where Cassiopeia A was tracked with the Swedish
    LOFAR telescope HBA band starting at '2019-08-03T12:00:00' and lasting
    24 hours (86400s) sampled every hour (3600s):

    >>> from dreambeam.rime.scenarios import *
    >>> from datetime import datetime, deltatime
    >>> obstimebeg = datetime.strptime('2019-08-30T12:00:00',
    ... "%Y-%m-%dT%H:%M:%S")
    >>> duration = deltatime(hours=24)
    >>> obstimestp = deltatime(hours=1)
    >>> pointingdir = (6.11, 1.02, 'J2000')
    >>> samptimes, freqs, jones, jonesobj = on_pointing_axis_tracking('LOFAR',
    ... 'HBA', 'Hamaker', 'SE607', obstimebeg, duration, obstimestp,
    ... pointingdir)
    >>> print(jones.shape)
    (1024, 25, 2, 2)
    >>> print(freqs.shape)
    (1024,)
    >>> print(samptimes.shape)
    (25,)
    >>> print(jones[512,0,:,:])
    [[ 0.35038040-0.02403619j  0.46850486+0.01755369j]
     [ 0.39298880-0.02620409j -0.38686167-0.01691861j]]

    """
    # Startup a telescope wizard
    tw = TelescopesWiz()

    # Get the telescopeband instance:
    telescope = tw.getTelescopeBand(telescopename, band, antmodel)

    #    *Setup Source*
    srcfld = dreambeam.rime.jones.DualPolFieldPointSrc(pointingdir)

    stnBD = telescope['Station'][stnid]
    stnRot = stnBD.stnRot

    #    *Setup PJones*
    timespy = []
    nrTimSamps = int((obsdur.total_seconds()/obstimestp.seconds))+1
    for ti in range(0, nrTimSamps):
        timespy.append(obstimebeg+ti*obstimestp)
    pjones = dreambeam.rime.jones.PJones(timespy, np.transpose(stnRot),
                                         do_parallactic_rot=do_parallactic_rot)

    #    *Setup EJones*
    ejones = stnBD.getEJones(pointingdir)
    stnDPolel = stnBD.feed_pat
    freqs = stnDPolel.getfreqs()

    #    *Setup MEq*
    pjonesOfSrc = pjones.op(srcfld)
    jonesobj = ejones.op(pjonesOfSrc)

    # Get the resulting Jones matrices
    # (structure is Jn[freqIdx, timeIdx, chanIdx, compIdx] )
    jones = jonesobj.getValue()
    if not do_parallactic_rot:
        basis_from = jonesobj.get_basis()
        basis_to = pjonesOfSrc.get_basis()
        btransmat2d = basis2basis_transf(basis_from, basis_to)[..., 1:, 1:]
        # Tranformation from IAU2C09 has to be (right) transformed back 1st
        transmat2d = np.matmul(C09toIAU[1:, 1:], btransmat2d)
        jones = np.matmul(jones, transmat2d)
        jonesobj.jones = jones
    return timespy, freqs, jones, jonesobj


def beamfov(telescopename, stnid, band, antmodel, freq,
            pointing=(0., np.pi/2, 'STN'), obstime=None, lmgrid=None):
    """Computes the Jones matrix over the beam fov for pointing.
    """
    # Startup a telescope wizard
    tw = TelescopesWiz()

    # Get the telescopeband instance:
    telescope = tw.getTelescopeBand(telescopename, band, antmodel)
    stnBD = telescope['Station'][stnid]
    stnRot = stnBD.stnRot

    #    *Setup Source*
    (az, el, refframe) = pointing
    srcfld = dreambeam.rime.jones.DualPolFieldRegion(refframe, iaucmp=False,
                                                     lmgrid=lmgrid)

    #    *Setup Parallatic Jones*
    pjones = dreambeam.rime.jones.PJones([obstime], np.transpose(stnRot))

    #    *Setup EJones*
    stnDPolel = stnBD.feed_pat
    #    **Select frequency
    freqs = stnDPolel.getfreqs()
    frqIdx = np.where(np.isclose(freqs, freq, atol=190e3))[0][0]
    # N.B. Ejones doesn't really use pointing
    ejones = stnBD.getEJones(pointing, [freqs[frqIdx]])

    #    *Setup MEq*
    pjones_src = pjones.op(srcfld)
    if refframe == dreambeam.rime.jones.Jones._topo_frame:
        j2000basis = pjones_src.jonesbasis
        # Since refFrame is STN, pjones is inverse of ordinary J2000 to STN.
        # By inverting it, one gets the ordinary conversion back.
        pjones_src = dreambeam.rime.jones.inverse(pjones_src)
    else:
        j2000basis = srcfld.jonesbasis
    res = ejones.op(pjones_src)
    # Because we started off with iaucmp=False, but want IAU components:
    res.convert2iaucmp()

    dreambeam.rime.jones.fix_imaginary_directions(res)
    # NOTE: Not using get_basis() method, in order to get station basis instead
    # of antenna basis:
    stnbasis = res.jonesbasis
    # Get the resulting Jones matrices
    # (structure is Jn[freqIdx, timeIdx, chanIdx, compIdx] )
    res_jones = res.getValue()

    return res_jones, stnbasis, j2000basis


def display_pointings(jones, obsinfo=None, do_3D=False,
                      do_parallactic_rot=None):
    """Display pointings in topocentric station coordinates, the antenna
    basis and the antenna response to IAU x and y. The plots are based on the
    final cumulative jones.
    For informative labelling, include the optional argument obsinfo, which
    contains fields 'stnid', 'band', 'freq', 'starttime', and 'antmodel'.
    """
    def plotsphgrid():
        # Plot grid for the LCL spherical crd sys
        nr_theta_ticks = 10  # 0..90 => 10 deg
        nr_phi_ticks = 4*4+1
        theta_ticks = np.linspace(0., np.pi/2, nr_theta_ticks)
        phi_ticks = np.linspace(0., 2*np.pi, nr_phi_ticks, endpoint=True)
        sg_linewidth = 0.5
        stn_tick_clrmrk = 'y:'
        itrf_tick_clrmrk = 'g:'
        # Compute iso-phi tick lines:
        for phi_tick in phi_ticks:
            xth = np.sin(theta_ticks)*np.sin(phi_tick)
            yth = np.sin(theta_ticks)*np.cos(phi_tick)
            zth = np.cos(theta_ticks)*np.ones((nr_theta_ticks,))
            xyzth_itrf = np.matmul(jones.stnRot.T, [xth, yth, zth])
            (xth_itrf, yth_itrf, zth_itrf) = (xyzth_itrf[0], xyzth_itrf[1],
                                              xyzth_itrf[2])
            if hidebelowhrz:
                abovehrz = zth_itrf > 0
                xth_itrf = xth_itrf[abovehrz]
                yth_itrf = yth_itrf[abovehrz]
                zth_itrf = zth_itrf[abovehrz]
            # plot iso-phi lines
            if do_3D:
                ax.plot(xth, yth, zth, stn_tick_clrmrk, linewidth=sg_linewidth)
                ax.plot(xth_itrf, yth_itrf, zth_itrf, itrf_tick_clrmrk,
                        linewidth=sg_linewidth)
            else:
                ax.plot(xth, yth, stn_tick_clrmrk, linewidth=sg_linewidth)
                ax.plot(xth_itrf, yth_itrf, itrf_tick_clrmrk,
                        linewidth=sg_linewidth)

        # Compute iso-theta tick lines:
        # N.B. the curvaure of iso-theta lines requires more phi points
        phi_xtrticks = np.linspace(0., 2*np.pi, 361)
        nr_phi_xtrticks = len(phi_xtrticks)
        for theta_tick in theta_ticks:
            xph = np.sin(theta_tick)*np.sin(phi_xtrticks)
            yph = np.sin(theta_tick)*np.cos(phi_xtrticks)
            zph = np.cos(theta_tick)*np.ones((nr_phi_xtrticks,))
            xyzph_itrf = np.matmul(jones.stnRot.T, [xph, yph, zph])
            (xph_itrf, yph_itrf, zph_itrf) = (xyzph_itrf[0], xyzph_itrf[1],
                                              xyzph_itrf[2])
            # plot iso-theta lines
            if do_3D:
                ax.plot(xph, yph, zph, stn_tick_clrmrk,
                        linewidth=sg_linewidth)
                ax.plot(xph_itrf, yph_itrf, zph_itrf, itrf_tick_clrmrk,
                        linewidth=sg_linewidth)
            else:
                ax.plot(xph, yph, stn_tick_clrmrk, linewidth=sg_linewidth)
                ax.plot(xph_itrf, yph_itrf, itrf_tick_clrmrk,
                        linewidth=sg_linewidth)

    nrsamps = jones.jonesbasis.shape[0]

    jn = jones.getValue()
    jnf = jn[256, :, :, :].squeeze()  # Midpoint freq.

    # NCP is z-base-vec of ITRF in stn crdsys
    itrf_z_stn = np.matmul(jones.stnRot.T, [[0], [0], [1]])

    # Pointings in Cartesian station crds
    xp = jones.jonesbasis[:, 0, 0]
    yp = jones.jonesbasis[:, 1, 0]
    zp = jones.jonesbasis[:, 2, 0]

    # Cartesian resp. of antenna basis in Ludwig3 w.r.t stn crdsys
    jbant = jones.get_basis()

    # N.B: imag part of ant response not used:
    jbresp = np.real(np.matmul(jbant[:, :, 1:], jnf))

    # Optionally remove data below stn horizon
    hidebelowhrz = True
    if hidebelowhrz:
        abovehrz = zp > 0
        xp = xp[abovehrz]
        yp = yp[abovehrz]
        zp = zp[abovehrz]
        jbant = jbant[abovehrz]
        jbresp = jbresp[abovehrz]
        nrsamps = len(zp)

    # Plot using 3d or 2d. (2d uses orthographic projection)
    fig = plt.figure()
    mplprojection = '3d' if do_3D else None
    ax = fig.add_subplot(111, projection=mplprojection)
    plotsphgrid()

    # Display pointings in horizontal coordinates
    if do_3D:
        ax.scatter(xp, yp, zp, c='c', marker='.')
    ax.plot(xp, yp, 'c.', label='Pointing')

    # Mark out start point
    ax.plot([xp[0]], [yp[0]], 'rP', label='Start')

    # Plot antenna dipole basis
    s = 0.1
    for j in range(1, 3):
        lw = 2 if j == 1 else 1
        ant = 'X' if j == 1 else 'Y'
        for i in range(nrsamps):
            if do_3D:
                ax.plot([xp[i], xp[i]+s*jbant[i, 0, j]],
                        [yp[i], yp[i]+s*jbant[i, 1, j]],
                        [zp[i], zp[i]+s*jbant[i, 2, j]],
                        'm', linewidth=lw)
                # ax.quiver()
            else:
                ax.plot([xp[i], xp[i]+s*jbant[i, 0, j]],
                        [yp[i], yp[i]+s*jbant[i, 1, j]],
                        'm', linewidth=lw)
        # Do it once again with label so that legend only does it once
        if do_3D:
            ax.plot([xp[i], xp[i]+s*jbant[i, 0, j]],
                    [yp[i], yp[i]+s*jbant[i, 1, j]],
                    [zp[i], zp[i]+s*jbant[i, 2, j]],
                    'm', linewidth=lw, label='antdip_'+ant)
            # ax.quiver()
        else:
            ax.plot([xp[i], xp[i]+s*jbant[i, 0, j]],
                    [yp[i], yp[i]+s*jbant[i, 1, j]],
                    'm', linewidth=lw, label='antdip_'+ant)

    # Plot Jones antenna X & Y-channels
    s = 0.2
    for j in range(2):
        lw = 2 if j == 0 else 1
        respiaucmp = 'x' if j == 0 else 'y'
        for i in range(nrsamps):
            if do_3D:
                ax.plot([xp[i], xp[i]+s*jbresp[i, 0, j]],
                        [yp[i], yp[i]+s*jbresp[i, 1, j]],
                        [zp[i], zp[i]+s*jbresp[i, 2, j]],
                        'b', linewidth=lw)
            else:
                ax.plot([xp[i], xp[i]+s*jbresp[i, 0, j]],
                        [yp[i], yp[i]+s*jbresp[i, 1, j]],
                        'b', linewidth=lw)
        # Do it once again with label so that legend only does it once
        if do_3D:
            ax.plot([xp[i], xp[i]+s*jbresp[i, 0, j]],
                    [yp[i], yp[i]+s*jbresp[i, 1, j]],
                    [zp[i], zp[i]+s*jbresp[i, 2, j]],
                    'b', linewidth=lw, label='respSKY_'+respiaucmp)
        else:
            ax.plot([xp[i], xp[i]+s*jbresp[i, 0, j]],
                    [yp[i], yp[i]+s*jbresp[i, 1, j]],
                    'b', linewidth=lw, label='respSKY_'+respiaucmp)

    # Plot NCP (ITRF z-base in STN crdsys)
    if do_3D:
        ax.plot(itrf_z_stn[0], itrf_z_stn[1], itrf_z_stn[2], 'y*', label='NCP')
    else:
        ax.plot(itrf_z_stn[0], itrf_z_stn[1], 'y*', label='NCP')

    # Fix plot settings
    title = "Pointing map"
    if obsinfo:
        title += """ [{} STN crdsys], Band: {}, Freq: {:.2f} MHz
Start @ {}, Model: {}, Pararot: {}"""\
                 .format(obsinfo['stnid'], obsinfo['band'],
                         obsinfo['freq']/1e6,
                         obsinfo['starttime'].isoformat()+' UT',
                         obsinfo['antmodel'],
                         do_parallactic_rot)
    # Plot origin
    ax.plot([0.], [0.], 'k.', label='Origin/Zenith')
    if do_3D:
        ax.set_xlim3d(left=-1.0, right=1.0)
        ax.set_ylim3d(bottom=-1.0, top=1.0)
        ax.set_zlim3d(bottom=0.0, top=1.0)
        ax.text(0, 1, 0, 'N (stn)')
        ax.text(1, 0, 0, 'E (stn)')
        ax.set_zticks([])
    else:
        ax.text(0, 1, 'N (stn)')
        ax.text(1, 0, 'E (stn)')
        ax.axis('equal')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.title(title)
    ax.legend(numpoints=1, loc='lower left')
    plt.draw()
