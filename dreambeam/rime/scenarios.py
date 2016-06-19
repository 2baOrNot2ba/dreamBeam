'''
This module provides some common observational scenarios of interest.
It implements basic Jones chains or Measurement Equations.
'''
import numpy as np
import matplotlib.pyplot as plt
from antpat.reps.sphgridfun import pntsonsphere
import dreambeam.rime.jones

def on_pointing_axis_tracking(telescope, stnID, ObsTimeBeg, duration,
                              ObsTimeStp, CelDir):
    """Computes the Jones matrix along pointing axis while tracking a fixed
    celestial source. """ #Fix: Doesn't use freq
    #    *Setup Source*
    celAz, celEl, celRef = CelDir.split(',')
    celAz = float(celAz)
    celEl = float(celEl)
    srcfld = dreambeam.rime.jones.DualPolFieldPointSrc((celAz, celEl, celRef))
    
    #    *Setup Parallatic Jones*
    #duration = ObsTimeEnd-ObsTimeBeg
    timespy = []
    nrTimSamps = int((duration.total_seconds()/ObsTimeStp.seconds))+1
    for ti in range(0, nrTimSamps):
        timespy.append(ObsTimeBeg+ti*ObsTimeStp)
    pjones = dreambeam.rime.jones.PJones(timespy)

    #    *Setup EJones*
    stnBD = telescope['Station'][stnID]
    ejones = stnBD.getEJones(CelDir)
    stnRot = stnBD.stnRot
    stnDPolel = stnBD.feed_pat

    #    *Setup MEq*
    pjonesOfSrc = pjones.op(srcfld)
    res = ejones.op(pjonesOfSrc)

    #Get the resulting Jones matrices 
    Jn = res.getValue()
    compute_paral(srcfld, stnRot, res, pjonesOfSrc)
    freqs = stnDPolel.getfreqs()
    return timespy, freqs, Jn


def compute_paral(srcfld, stnRot, res, pjonesOfSrc):
    """Compute parallactic rotation. Also display's pointings in horizontal coordinates."""
    #print("Parallactic rotation matrix:")
    srcbasis = srcfld.jonesbasis
    basisITRF_lcl = res.jonesbasis
    basisJ2000_ITRF = pjonesOfSrc.jonesbasis
    ax = plt.subplot(111, projection='polar')
    nrsamps=basisITRF_lcl.shape[0]
    az = np.zeros((nrsamps))
    el = np.zeros((nrsamps))
    for i in range(nrsamps):
        basisJ2000_ITRF_to = np.matmul(stnRot, basisJ2000_ITRF[i,:,:])
        paramat = np.matmul(basisITRF_lcl[i,:,:].T, basisJ2000_ITRF_to)
        #print paramat
        az[i], el[i] = pntsonsphere.crt2sphHorizontal(basisITRF_lcl[i,:,0].squeeze())
    
    #Display pointings in horizontal coordinates
    #print("th, ph", np.rad2deg(np.array([np.pi/2-el, az]).T))
    ax.plot(az, 90-el/np.pi*180,'+')
    #Mark out start point
    ax.plot(az[0], 90-el[0]/np.pi*180,'r8')
    ax.set_rmax(90)
    plt.draw()