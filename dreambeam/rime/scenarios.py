'''
This module provides some common observational scenarios of interest.
It implements basic Jones chains or Measurement Equations.
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import dreambeam.rime.jones
from dreambeam.rime.conversionUtils import crt2sph


def on_pointing_axis_tracking(telescope, stnID, ObsTimeBeg, duration,
                              ObsTimeStp, CelDir, do_parallactic_rot=True):
    """Computes the Jones matrix along pointing axis while tracking a fixed
    celestial source. """  # # FIXME: Doesn't use freq
    #    *Setup Source*
    # celAz, celEl, celRef = CelDir.split(',')
    (celAz, celEl, celRef) = CelDir
    srcfld = dreambeam.rime.jones.DualPolFieldPointSrc((celAz, celEl, celRef))

    stnBD = telescope['Station'][stnID]
    stnRot = stnBD.stnRot

    #    *Setup PJones*
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

    # Get the resulting Jones matrices
    # (structure is Jn[freqIdx, timeIdx, chanIdx, compIdx] )
    Jn = res.getValue()
    freqs = stnDPolel.getfreqs()
    return timespy, freqs, Jn, res


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
    frqIdx = np.where(np.isclose(freqs, freq, atol=190e3))[0][0]
    # N.B. Ejones doesn't really use CelDir
    ejones = stnBD.getEJones(CelDir, [freqs[frqIdx]])

    #    *Setup MEq*
    pjonesOfSrc = pjones.op(srcfld)
    res = ejones.op(pjonesOfSrc)

    # Get the resulting Jones matrices
    # (structure is Jn[freqIdx, timeIdx, chanIdx, compIdx] )
    Jn = res.getValue()
    # res.get_basis()  # Get basis set for cumulative Jones
    return srcfld.azmsh, srcfld.elmsh, Jn, ejones.thisjones


def display_pointings(jones, obsinfo=None, do_3D=False,
                      do_parallactic_rot=None):
    """Display pointings in topocentric station coordinates, the antenna
    basis and the antenna response to IAU x and y. The plots are based on the
    final cumulative jones.
    For informative labelling, include the optional arguments obsinfo, The I The time of the first Jones
    matrix obstimebeg is used for labelling.
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
                abovehrz =  zth_itrf > 0
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
                    'b', linewidth=lw, label='respIAU_'+respiaucmp)
        else:
            ax.plot([xp[i], xp[i]+s*jbresp[i, 0, j]],
                    [yp[i], yp[i]+s*jbresp[i, 1, j]],
                    'b', linewidth=lw, label='respIAU_'+respiaucmp)

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
                         obsinfo['starttime'].isoformat()+'UT',
                         obsinfo['antmodel'],
                         do_parallactic_rot)
    # Plot origin
    ax.plot([0.], [0.], 'k.', label='Origin')
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
