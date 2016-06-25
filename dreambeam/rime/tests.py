#!/usr/bin/python

import sys
from datetime import datetime, timedelta
from jones import *
from conversionUtils import *
from dreambeam.telescopes.LOFAR.native.parseAntennaField import getArrayBandParams


def setupObsInstance():
   #Observation interval
   beginTime=datetime(2011,10,25,0,0,0)
   endTime=datetime(2011,10,26,0,0,0)
   stepTime=timedelta(minutes=60)
   td=endTime-beginTime
   Times=[]
   nrTimSamps=int(td.total_seconds()/stepTime.seconds)+1
   for ti in range(0,nrTimSamps):
       Times.append(beginTime+ti*stepTime)


   #Source direction
   ##CasA
   celSrcTheta=0.548876961
   celSrcPhi=6.11378655886310
   celSrcDir= celSrcPhi, (math.pi/2-celSrcTheta), 'J2000'
   
   #Station position and rotation
   #stnPos_me=measures().position('wgs84','12deg','58deg','0m')
   stnPos, stnRot, stnRelPos = getArrayBandParams('SE607', 'LBA')

   return Times, celSrcDir, stnPos, stnRot


def tgetParallacticRot():
   Times, celSrcDir, stnPos, stnRot=setupObsInstance()
   RotP=getParallacticRot(Times, stnPos, celSrcDir, doPolPrec=False)
   printJones(RotP)


def tPJones():
   Times, celSrcDir, stnPos, stnRot=setupObsInstance()
   srcfld=DualPolFieldPointSrc(celSrcDir[:2])
   pjones=PJones(Times)
   #snk=DualPolFieldSink()
   res=pjones.op(srcfld)
   print("Value", res.getValue())
   print("Basis", res.getBasis())


def tModFuncs_CEL2TOPOpnts():
   Times, celSrcDir, stnPos, stnRot=setupObsInstance()
   CelRot, rotang=CEL2TOPOpnts(Times, stnPos, celSrcDir)
   print CelRot
   #print np.rad2deg(rotang)
   #RotP=getParallacticRot(Times, stnPos, stnRot, celSrcDir, doPolPrec=False)
   #printJones(RotP)


def tModFuncs_crt2sph():
    azi=np.array([0.1, 0.2])
    ele=np.array([-0.3, 0.4])
    dir_crt=sph2crt(azi,ele)
    print dir_crt
    dir_sph=crt2sph(dir_crt)


def tcomputeSphBasis():
  Times, celSrcDir, stnPos, stnRot=setupObsInstance()
  celSrcDir=(0., 0*1.570796327)
  obsTimesArr, obsTimeUnit= pyTimes2meTimes(Times)
  jonesbasis=np.array(getSph2CartTransf(sph2crt(celSrcDir[0], celSrcDir[1])))
  for ti in range(0,len(obsTimesArr)):
    me=setEpoch(obsTimesArr[ti], obsTimeUnit)
    jonesrbasis_to=np.asmatrix(convertBasis(me, jonesbasis, 'J2000', 'ITRF'))
    jonesbasisMat=getSph2CartTransf(jonesrbasis_to[:,0])
    #print jonesrbasis_to[:,1:]
    print ti, jonesbasisMat[:,1:].H*jonesrbasis_to[:,1:]


if __name__ == '__main__':
  pass
  #tPJones()
  #tgetParallacticRot()
  #tModFuncs_CEL2TOPOpnts()
  #tModFuncs_crt2sph_me()
  #tcomputeSphBasis()
