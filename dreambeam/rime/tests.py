#!/usr/bin/python
"""Some tests of `rime` package."""

import math
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from .jones import PJones, DualPolFieldPointSrc
from .conversion_utils import CEL2TOPOpnts, sph2crt_me, getParallacticRot, \
                            printJones, pyTimes2meTimes, sph2crt, crt2sph, \
                            getSph2CartTransf, setEpoch, convertBasis
import antpat.reps.sphgridfun.pntsonsphere
from casacore.measures import measures
import dreambeam.telescopes.geometry_ingest as gi
from dreambeam.rime.scenarios import on_pointing_axis_tracking


def setupObsInstance():
    # Observation interval (approx vernal equinox)
    beginTime = datetime(2011, 3, 20, 6, 0, 0)
    endTime = datetime(2011, 3, 21, 6, 0, 0)
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
    celSrcTheta = 0.4*math.pi/2
    celSrcEl = 0.6  # (math.pi/2-celSrcTheta)
    celSrcPhi = 0.0
    celSrcDir = celSrcPhi, celSrcEl, 'J2000'

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


def test_getParallacticRot():
    Times, celSrcDir, stnPos, stnRot = setupObsInstance()
    RotP = getParallacticRot(Times, stnPos, celSrcDir, doPolPrec=False)
    printJones(RotP)


def test_PJones():
    Times, celSrcDir, stnPos, stn2ITRFrot = setupObsInstance()
    srcfld = DualPolFieldPointSrc(celSrcDir)
    pjones = PJones(Times, np.transpose(stn2ITRFrot))
    res = pjones.op(srcfld)
    print("Value", res.getValue())
    print("Basis", res.get_basis())


def test_CEL2TOPOpnts():
    Times, celSrcDir, stnPos, stnRot = setupObsInstance()
    CelRot, rotang = CEL2TOPOpnts(Times, stnPos, celSrcDir)
    print("CelRot", CelRot)
    # print(np.rad2deg(rotang))
    # RotP=getParallacticRot(Times, stnPos, stnRot, celSrcDir, doPolPrec=False)
    # printJones(RotP)


def test_crt2sph():
    azi = np.array([0.1, 0.2])
    ele = np.array([-0.3, 0.4])
    dir_crt = sph2crt(azi, ele)
    print(dir_crt)
    dir_sph = crt2sph(dir_crt)
    print(dir_sph)


def test_computeSphBasis():
    Times, celSrcDir, stnPos, stnRot = setupObsInstance()
    obsTimesArr, obsTimeUnit = pyTimes2meTimes(Times)
    jonesbasis = np.array(getSph2CartTransf(sph2crt(celSrcDir[0],
                                                    celSrcDir[1])))
    print(jonesbasis)
    for ti in range(0, len(obsTimesArr)):
        me = setEpoch(obsTimesArr[ti], obsTimeUnit)
        jonesrbasis_to = np.asmatrix(convertBasis(me, jonesbasis, 'J2000',
                                                  'ITRF'))
        jonesbasisMat = getSph2CartTransf(jonesrbasis_to[:, 0])
        print(jonesrbasis_to[:, 0])
        # print obsTimesArr[ti], jonesbasisMat[:,1:].H*jonesrbasis_to[:,1:]
        print(obsTimesArr[ti], jonesrbasis_to, jonesbasisMat)


def test_Stokes():
    freq = 80e6
    # cohmatt = 0.5*np.array([[2-1, 0+0.8j], [0-0.8j, 2+1]])  # I,Q,U,V=2,1,0,0.8
    cohmatt = 0.5*np.array([[2-2, 0], [0, 2+2]])
    samptimes, srcdir, stnPos, stnRot = setupObsInstance()
    btime = samptimes[0]
    duration = (samptimes[-1]-samptimes[0])
    steptime = (samptimes[1]-samptimes[0])  # .total_seconds()
    telescopename = 'LOFAR'  # TelescopesWiz().getTelescopeBand('LOFAR', 'LBA', 'Hamaker')
    stnid = 'SE607'
    band = 'LBA'
    antmodel = 'Hamaker-default'
    timespy, freqs, Jn, res = \
        on_pointing_axis_tracking(telescopename, stnid, band, antmodel, btime,
                                  duration, steptime, srcdir)
    frqIdx = np.where(np.isclose(freqs, freq, atol=190e3))[0][0]
    Jnf = Jn[frqIdx, :, :, :].squeeze()
    if True:
        plt.figure()
        plt.plot(np.real(Jnf[:,0,0]))
        plt.plot(np.real(Jnf[:,0,1]))
        plt.show()
    sI = np.zeros(len(timespy))
    sQ = np.zeros(len(timespy))
    sU = np.zeros(len(timespy))
    sV = np.zeros(len(timespy))
    for ti in range(len(timespy)):
        jones = Jnf[ti, :, :].squeeze()
        cohmat = np.matmul(np.matmul(jones, cohmatt), jones.conj().T)
        sI[ti] = +np.real(cohmat[0, 0]+cohmat[1, 1])
        sQ[ti] = +np.real(cohmat[0, 0]-cohmat[1, 1])
        sU[ti] = +np.real(cohmat[0, 1]+cohmat[1, 0])
        sV[ti] = +np.imag(cohmat[0, 1]-cohmat[1, 0])
    plt.figure()
    plt.plot(timespy, sI, label='SI')
    plt.plot(timespy, sQ, label='SQ')
    plt.plot(timespy, sU, label='SU')
    plt.plot(timespy, sV, label='SV')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    pass
    # print(IAU_pol_basis(0.*np.pi, 0.0001*np.pi))
    # test_PJones()
    # test_getParallacticRot()
    # test_CEL2TOPOpnts()
    # test_crt2sph()
    # test_computeSphBasis()
    test_Stokes()
