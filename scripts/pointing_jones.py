#!/usr/bin/env python
"""Model of a LOFAR station. Gets Jones matrix towards a given direction
   and frequency.
"""
import sys
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from antpat.dualpolelem import plot_polcomp_dynspec
from dreambeam.rime.scenarios import on_pointing_axis_tracking
from dreambeam.telescopes.rt import TelescopesWiz


def printJonesFreq(timespy, Jnf):
    #Select one frequency
    for ti in range(len(timespy)):
        print freq, timespy[ti], Jnf[ti,0,0], Jnf[ti,0,1], Jnf[ti,1,0], Jnf[ti,1,1]
        #Print out data for BST-mode comparison (ie powers of p & q channels):
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


def getnextcmdarg(args, mes):
    try:
        arg = args.pop(0)
    except IndexError:
        print("Specify "+mes)
        print(USAGE)
        exit()
    return arg

SCRIPTNAME = sys.argv[0].split('/')[-1]
USAGE = "Usage:\n  {} print|plot telescope band stnID beammodel beginUTC duration timeStep pointingRA pointingDEC [frequency]".format(SCRIPTNAME)
#Example: 
#$ pointing_jones.py print LOFAR LBA SE607 Hamaker 2012-04-01T01:02:03 60 1 6.11 1.02 60E6
if __name__ == "__main__":
    #Startup a telescope wizard
    TW = TelescopesWiz()
    #Process cmd line arguments
    args = sys.argv[1:] 
    action = getnextcmdarg(args, "output-type:\n  'print' or 'plot'")
    telescopeName = getnextcmdarg(args, "telescope:\n  "+', '.join(TW.get_telescopes()))
    band = getnextcmdarg(args, "band/feed:\n  "+', '.join(TW.get_bands(telescopeName)))
    stnID =  getnextcmdarg(args, "station-ID:\n  "+ ', '.join(TW.get_stations(telescopeName, band)))
    antmodel = getnextcmdarg(args, "beam-model:\n  "+', '.join(TW.get_beammodels(telescopeName, band)))
    try:
        bTime = datetime.strptime(args[0], "%Y-%m-%dT%H:%M:%S")
    except IndexError:
        print("Specify start-time [UTC].")
        print(USAGE)
        exit()
    try:
        duration =timedelta(0,float(args[1]))
    except IndexError:
        print("Specify duration.")
        print(USAGE)
        exit()
    try:
        stepTime =timedelta(0,float(args[2]))
    except IndexError:
        print("Specify step-time.")
        print(USAGE)
        exit()
    try:
        #ra=args[4]+'rad'
        #dec=args[5]+'rad'
        CelDir=str(args[3])+','+str(args[4])+',J2000'
    except IndexError:
        print("Specify pointing direction.")
        print(USAGE)
        exit()
    try:
        freq=float(args[5])
    except IndexError:
        freq=0.

    #Get the telescopeband instance:
    telescope = TW.getTelescopeBand(telescopeName, band, antmodel)
    #Compute the Jones matrices
    timespy, freqs, Jn = on_pointing_axis_tracking(telescope, stnID, bTime, duration,
                                               stepTime, CelDir) #fix: freq not used
    #Do something with resulting Jones according to cmdline args
    if freq == 0.:
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
