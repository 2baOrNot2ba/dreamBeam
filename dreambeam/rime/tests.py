#!/usr/bin/python

import math
from datetime import datetime, timedelta
import numpy as np
from jones import PJones, DualPolFieldPointSrc
from conversionUtils import CEL2TOPOpnts, sph2crt_me, getParallacticRot, \
                            printJones, pyTimes2meTimes, sph2crt, crt2sph, \
                            getSph2CartTransf, setEpoch, convertBasis
import antpat.reps.sphgridfun
from casacore.measures import measures
import dreambeam.telescopes.geometry_ingest as gi


def setupObsInstance():
    # Observation interval (approx vernal equinox)
    beginTime = datetime(2011, 3, 20, 0, 0, 0)
    endTime = datetime(2011, 3, 21, 0, 0, 0)
    stepTime = timedelta(minutes=60)
    td = endTime-beginTime
    Times = []
    nrTimSamps = int(td.total_seconds()/stepTime.seconds)+1
    for ti in range(0, nrTimSamps):
        Times.append(beginTime+ti*stepTime)

    # Source direction
    #   CasA:
    # celSrcTheta_CasA = np.pi/2-1.026515
    # celSrcPhi_CasA = 6.123487
    #   Celestial origin:
    celSrcTheta = 0.001*math.pi/2
    celSrcPhi = 0.
    celSrcDir = celSrcPhi, (math.pi/2-celSrcTheta), 'J2000'

    # Station position and rotation
    #   Alt1 arbitrarily:
    me = measures()
    stnPos_meWGS = measures().position('wgs84', '0deg', '0deg', '0m')
    stnPos_meITRF = me.measure(stnPos_meWGS, 'ITRF')
    stnPos = stnPos_meITRF['m2']['value']*sph2crt_me(stnPos_meITRF)[:,np.newaxis]
    stnRot = antpat.reps.sphgridfun.pntsonsphere.rot3Dmat(0., 0., 1*math.pi/2)
    #   Alt2 set via official LOFAR geodetic data:
    stnPos, stnRot = gi.getArrayBandParams('LOFAR', 'SE607', 'LBA')

    return Times, celSrcDir, stnPos, stnRot


def tgetParallacticRot():
    Times, celSrcDir, stnPos, stnRot = setupObsInstance()
    RotP = getParallacticRot(Times, stnPos, celSrcDir, doPolPrec=False)
    printJones(RotP)


def tPJones():
    Times, celSrcDir, stnPos, stn2ITRFrot = setupObsInstance()
    srcfld = DualPolFieldPointSrc(celSrcDir)
    pjones = PJones(Times, np.transpose(stn2ITRFrot))
    res = pjones.op(srcfld)
    print("Value", res.getValue())
    print("Basis", res.getBasis())


def tModFuncs_CEL2TOPOpnts():
    Times, celSrcDir, stnPos, stnRot = setupObsInstance()
    CelRot, rotang = CEL2TOPOpnts(Times, stnPos, celSrcDir)
    print "CelRot", CelRot
    #print np.rad2deg(rotang)
    #RotP=getParallacticRot(Times, stnPos, stnRot, celSrcDir, doPolPrec=False)
    #printJones(RotP)


def tModFuncs_crt2sph():
    azi = np.array([0.1, 0.2])
    ele = np.array([-0.3, 0.4])
    dir_crt = sph2crt(azi, ele)
    print dir_crt
    dir_sph = crt2sph(dir_crt)
    print dir_sph


def tcomputeSphBasis():
    Times, celSrcDir, stnPos, stnRot = setupObsInstance()
    obsTimesArr, obsTimeUnit = pyTimes2meTimes(Times)
    jonesbasis = np.array(getSph2CartTransf(sph2crt(celSrcDir[0],
                                                    celSrcDir[1])))
    for ti in range(0, len(obsTimesArr)):
        me = setEpoch(obsTimesArr[ti], obsTimeUnit)
        jonesrbasis_to = np.asmatrix(convertBasis(me, jonesbasis, 'J2000',
                                                  'ITRF'))
        jonesbasisMat = getSph2CartTransf(jonesrbasis_to[:, 0])
        print jonesrbasis_to[:, 0]
        #print obsTimesArr[ti], jonesbasisMat[:,1:].H*jonesrbasis_to[:,1:]
        print obsTimesArr[ti], jonesrbasis_to, jonesbasisMat


if __name__ == '__main__':
    pass
    #print IAU_pol_basis(0.*np.pi, 0.0001*np.pi)
    tPJones()
    # tgetParallacticRot()
    # tModFuncs_CEL2TOPOpnts()
    # tModFuncs_crt2sph()
    # tcomputeSphBasis()
