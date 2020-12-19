import numpy as np
import matplotlib.pyplot as plt


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

    # Find starting point (Save in case below horizon)
    xp0, yp0, zp0 = xp[0], yp[0], zp[0]

    nrsamps = jones.jonesbasis.shape[0]

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
    label_start = 'Start' if zp0 > 0 else 'Start (below horizon)'
    ax.plot([xp0], [yp0], 'rP', label=label_start)

    # Plot antenna dipole basis
    s = 0.1
    for j in range(1, 3):
        lw = 2 if j == 1 else 1
        ant = 'X' if j == 1 else 'Y'
        for i in range(nrsamps):
            # Label only first samp so that legend only has it once
            label = 'antdip_'+ant if i == 0 else None
            if do_3D:
                ax.plot([xp[i], xp[i]+s*jbant[i, 0, j]],
                        [yp[i], yp[i]+s*jbant[i, 1, j]],
                        [zp[i], zp[i]+s*jbant[i, 2, j]],
                        'm', linewidth=lw, label=label)
                # ax.quiver()
            else:
                ax.plot([xp[i], xp[i]+s*jbant[i, 0, j]],
                        [yp[i], yp[i]+s*jbant[i, 1, j]],
                        'm', linewidth=lw, label=label)

    # Plot Jones antenna X & Y-channels
    s = 0.2
    for j in range(2):
        lw = 2 if j == 0 else 1
        respiaucmp = 'x' if j == 0 else 'y'
        for i in range(nrsamps):
            # Label only first samp so that legend only has it once
            label = 'respSKY_'+respiaucmp if i == 0 else None
            if do_3D:
                ax.plot([xp[i], xp[i]+s*jbresp[i, 0, j]],
                        [yp[i], yp[i]+s*jbresp[i, 1, j]],
                        [zp[i], zp[i]+s*jbresp[i, 2, j]],
                        'b', linewidth=lw, label=label)
            else:
                ax.plot([xp[i], xp[i]+s*jbresp[i, 0, j]],
                        [yp[i], yp[i]+s*jbresp[i, 1, j]],
                        'b', linewidth=lw, label=label)

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
