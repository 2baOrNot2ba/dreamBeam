#!/usr/bin/env python
"""Model of a LOFAR station. Gets Jones matrix towards a given direction
   and frequency.
"""
import sys
sys.path.append('..')
import os
import optparse
import pickle
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from antpat.dualpolelem import plot_polcomp_dynspec
from antpat.io import NECread
from antpat.reps.sphgridfun import pntsonsphere
import rime.jones
import telescopes

projdir = os.path.dirname(os.path.abspath(__file__))
dataformdir = projdir+'/data_formats/'
NECdir = dataformdir+'/NEC_out/'


def inittelescope(name, patmodel):
    filename = "teldat_"+name+"_"+patmodel+".p"
    teldatdir = os.path.dirname(telescopes.__file__)+"/"+name+"/data/"
    telescope = pickle.load(open(teldatdir+filename,'rb'))
    return telescope


def computeJones(tele_name, ArrBand, stnID, modeltype,
               ObsTimeBeg, duration, ObsTimeStp,
               CelDir, freq):
    #    *Setup Source*
    celAz, celEl, celRef = CelDir.split(',')
    celAz = float(celAz)
    celEl = float(celEl)
    srcfld = rime.jones.DualPolFieldPointSrc((celAz, celEl, celRef))
    
    #    *Setup Parallatic Jones*
    #duration = ObsTimeEnd-ObsTimeBeg
    timespy = []
    nrTimSamps = int((duration.total_seconds()/ObsTimeStp.seconds))+1
    for ti in range(0, nrTimSamps):
       timespy.append(ObsTimeBeg+ti*ObsTimeStp)
    pjones = rime.jones.PJones(timespy)
    #    *Setup EJones*
    telescope = inittelescope(tele_name, modeltype)
    stnBD = telescope['Station'][stnID][ArrBand]
    ejones = stnBD.getEJones()
    stnRot = stnBD.stnRot
    stnDPolel = stnBD.stnDPolel
    #    *Setup MEq*
    pjonesOfSrc = pjones.op(srcfld)
    res = ejones.op(pjones.op(srcfld))
        
    #print np.matmul(basisITRF_lcl, basisJ2000_ITRF)
    #Do something with it
    Jn = res.getValue()
    print_paral(srcfld, stnRot, res, pjonesOfSrc)
    freqs = stnDPolel.getfreqs()
    return timespy, freqs, Jn


def printJonesFreq(timespy, Jnf):
    #Select one frequency
    for ti in range(len(timespy)):
        print freq, timespy[ti], Jnf[ti,0,0], Jnf[ti,0,1], Jnf[ti,1,0], Jnf[ti,1,1]
        #print("{0} {1} {2}".format(ti, np.abs(Jnf[ti,0,0])**2+np.abs(Jnf[ti,0,1])**2, np.abs(Jnf[ti,1,0])**2+np.abs(Jnf[ti,1,1])**2) )


def plotJonesFreq(timespy, Jnf):
    p_ch = np.abs(Jnf[:,0,0].squeeze())**2+np.abs(Jnf[:,0,1].squeeze())**2
    q_ch = np.abs(Jnf[:,1,1].squeeze())**2+np.abs(Jnf[:,1,0].squeeze())**2
    plt.figure()
    plt.subplot(211)
    plt.plot(timespy, 10*np.log10(p_ch))
    plt.title('p-channel')
    plt.subplot(212)
    plt.plot(timespy, 10*np.log10(q_ch))
    plt.title('q-channel')
    plt.xlabel('Time')
    plt.show()


def printAllJones(timespy, freqs, Jn):
    print "Time, Freq, J11, J12, J21, J22"
    #duration.seconds/ObsTimeStp.seconds
    for ti in range(0, len(timespy)):
        for fi,freq in enumerate(freqs):
            print timespy[ti].isoformat(), freq, Jn[fi,ti,0,0], Jn[fi,ti,0,1], Jn[fi,ti,1,0], Jn[fi,ti,1,1]


def plotAllJones(timespy, freqs, Jn):
    plot_polcomp_dynspec(timespy, freqs, Jn)


def print_paral(srcfld, stnRot, res, pjonesOfSrc):
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
    #print("th, ph", np.rad2deg(np.array([np.pi/2-el, az]).T))
    ax.plot(az, 90-el/np.pi*180,'+')
    #Mark out start point
    ax.plot(az[0], 90-el[0]/np.pi*180,'r8')
    ax.set_rmax(90)
    plt.draw()


def args2inpparms(args):
    stnName=args[0]
    bTime = datetime.strptime(args[1], "%Y-%m-%d %H:%M:%S")
    duration =timedelta(0,float(args[2]))
    stepTime =timedelta(0,float(args[3]))
    #ra=args[4]+'rad'
    #dec=args[5]+'rad'
    CelDir=str(args[4])+','+str(args[5])+',J2000'
    return stnName,bTime,duration,stepTime,CelDir

if __name__ == "__main__":
#    stnID="SE607"
#    ArrBand="LBA"
#    ObsTimeBeg=datetime(2015,12,29,0,0,0)
#    ObsTimeEnd=datetime(2015,12,29,22,0,0)
#    ObsTimeStp=timedelta(minutes=60)
#    CelDir="6.11378655886310,1.021919366,J2000"
#    modelLOFARobservation(stnID, ArrBand, ObsTimeBeg, ObsTimeEnd, ObsTimeStp,
#                          CelDir,  modeltype = 'Hamaker_Arts' )
    usage = "usage: %prog act telescope band antmodel stnName beginUTC duration timeStep pointingRA pointingDEC [frequency]"
    #Example: $ pointing_jones.py action LOFAR LBA Hamaker SE607 '2012-04-01 01:02:03' 60 1 0 0 60E6
    opt = optparse.OptionParser(usage=usage)
    options, args = opt.parse_args()
    if len(args) >0:
        action = args.pop(0)
        telescopeName = args.pop(0)
        band = args.pop(0)
        antmodel=args.pop(0)
    if len(args) == 6 or len(args) == 7:
        stnID,bTime,duration,stepTime,CelDir=args2inpparms(args)
        if len(args) == 7:
            freq=float(args[6])
        else:
            freq=0.
        timespy, freqs, Jn = computeJones(telescopeName, band, stnID, antmodel,
                   bTime, duration, stepTime, CelDir, freq)
        if len(args) == 6:
            if action == "plot":
                plotAllJones(timespy, freqs, Jn)
            else:
                printAllJones(timespy, freqs, Jn)
        else:
            frqIdx = np.where(np.isclose(freqs,freq,atol=190e3))[0][0]
            Jnf = Jn[frqIdx,:,:,:].squeeze()
            if action == "plot":
                plotJonesFreq(timespy, Jnf)
            else:
                printJonesFreq(timespy, Jnf)
    else :
        opt.error("incorrect number of arguments")

