"""
  This module provides a Jones matrix framework for radio interometric
  measurement equations.
"""
import math
import numpy as np
from casacore.measures import measures
from casacore.quanta import quantity
from conversionUtils import sph2crt, crt2sph, computeSphBasis, convertBasis, \
                            getSph2CartTransf, IAU_pol_basis
from antpat.reps.sphgridfun.tvecfun import getSph2CartTransfMat

class Jones(object):
    """This is the base class for Jones algebra. It contains the Jones matrix
    itself and a basis w.r.t. which the Jones matrix is given.
    The basis is such that:
        self.jonesbasis=array([[r_hat], [phi_hat], [theta_hat]]).
    """
    def __init__(self):
        pass
    
    def op(self,jonesobjright):
        """Operate this Jones on to the Jones passed in the argument."""
        self.jonesr=jonesobjright.getValue()
        self.jonesrbasis=jonesobjright.getBasis()
        self.jonesrmeta=jonesobjright.getMetadata()
        self.computeJonesRes()
        return self
        
    def getValue(self):
        """Return value of the Jones matrix"""
        return self.jones
        
    def getBasis(self):
        """Return basis of the Jones matrix"""
        return self.jonesbasis
    
    def getMetadata(self):
        """Return the data about the Jones matrix."""
        return self.jonesmeta
    
    def computeJonesRes(self):
        pass

class JonesChain(object):
    jonesproducts=[]
    def __init__(self):
        self.jonesproducts=[]


class PJones(Jones):
    """This is a P-Jones or parallactic Jones. This has a temporal dependence
    given by the epoch of observation."""
    def __init__(self, obsTimespy):
        super(PJones, self).__init__()
        obsTimes_lst = []
        for obsTimepy in obsTimespy:
            obsTimes_lst.append(quantity(obsTimepy.isoformat()).get_value())
        obsTimes_me = quantity(obsTimes_lst, 'd')
        self.obsTimes = obsTimes_me.get_value()
        self.obsTimeUnit = obsTimes_me.get_unit()
        self.jonesmeta = {}
        self.jonesmeta['refFrame'] = 'ITRF'
    
    def computeJonesRes(self):
        nrOfTimes = len(self.obsTimes)
        paraRot = np.zeros((nrOfTimes, 2, 2))
        me = measures()
        me.doframe(measures().position('ITRF', '0m', '0m', '0m'))
        self.jonesbasis = np.zeros((nrOfTimes, 3, 3))
        (az_from, el_from) = crt2sph(self.jonesrbasis[:,0])
        r_sph_me = measures().direction(self.jonesrmeta['refFrame'], 
                                        quantity(az_from, 'rad'),
                                        quantity(el_from, 'rad'))
        for ti in range(0, nrOfTimes):
            #Set current time in reference frame 
            timEpoch=me.epoch('UTC', quantity(self.obsTimes[ti], self.obsTimeUnit))
            me.doframe(timEpoch)
            
            #paraRot[ti,:,:]=computeParaMat_me(self.jonesrmeta['refFrame'],
            #                          self.jonesmeta['refFrame'], r_sph_me, me)
            
            jonesrbasis_to = np.asmatrix(convertBasis(me, self.jonesrbasis,
                                                    self.jonesrmeta['refFrame'],
                                                    self.jonesmeta['refFrame']))
            jonesrbasis_to2 = computeSphBasis(self.jonesrmeta['refFrame'], 
                                           self.jonesmeta['refFrame'], 0, r_sph_me, me)
            jonesbasisMat = getSph2CartTransf(jonesrbasis_to[:,0])
            paraRot[ti,:,:] = jonesbasisMat[:,1:].H*jonesrbasis_to[:,1:]
            #paraRot[ti,:,:]=jonesrbasis_to2*jonesbasisMat[:,1:]
            self.jonesbasis[ti,:,:] = jonesbasisMat
        self.jones = np.matmul(paraRot, self.jonesr)
        self.thisjones = paraRot


class DualPolFieldPointSrc(Jones):
    """This is a mock Jones point source. It does not model a real source. It's
    purpose is for testing. It can be seen as a source that first transmits in one
    polarization and then in another, then 2 transmissions given in the 2 columns.
    It may have a spectral dimension. The src_dir should be a tuple with
    (az, el, ref).""" 
    def __init__(self, src_dir, dualPolField=np.identity(2)):
        (src_az, src_el, src_ref)=src_dir
        self.jones = dualPolField
        #self.jonesbasis = np.array(getSph2CartTransf(sph2crt(src_az, src_el)))
        self.jonesbasis = IAU_pol_basis(src_az, src_el)
        self.jonesmeta = {}
        self.jonesmeta['refFrame'] = src_ref


class EJones(Jones):
    """This is the antenna or feed Jones. It is given by a set of complex gain
    patterns for each frequency and polarization channel."""
    def __init__(self, dualPolElem, position, stnRot, freqSel=None):
        self.position = position
        #posrot = getPosRot(position)
        #self.stnRot = npmatmul(posrot, offsetRot)
        self.stnRot = stnRot
        self.dualPolElem = dualPolElem
        self.dualPolElem.rotateframe(np.asarray(stnRot)) #FIX asarray from mat
        if freqSel is None:
            self.freqChan = self.dualPolElem.getfreqs()
        else:
            self.freqChan = freqSel

    def computeJonesRes(self):
        jonesrbasis_to = np.matmul(np.asarray(self.stnRot.T), self.jonesrbasis)
        (az_from, el_from) = crt2sph(self.jonesrbasis[:,:,0].squeeze().T)
        theta_phi_view = (np.pi/2-el_from, az_from)
        ejones = self.dualPolElem.getJonesAlong(self.freqChan, theta_phi_view)
        #(theta_lcl, phi_lcl) = self.dualPolElem.getBuildCoordinates(math.pi/2-r_sph[1], r_sph[0])
        #print theta_lcl, phi_lcl
        r_lcl = crt2sph(jonesrbasis_to[:,:,0].squeeze().T)
        #print np.rad2deg(r_lcl)
        jonesbasisMat = getSph2CartTransfMat(jonesrbasis_to[:,:,0].squeeze())
        #paraRot = np.matmul(np.conjugate(jonesbasisMat), jonesrbasis_to)
        self.jonesbasis = jonesbasisMat
        #This is the actual MEq multiplication:
        self.jones = np.matmul(ejones, self.jonesr)
        self.thisjones = ejones
        
    def getPosRot(self, position):
        """Compute the nominal transformation from the geodetic position to
        ITRF. (Not implemented yet)""" 
        return np.identity(2)


class DualPolFieldSink(Jones):
    def computeJonesRes(self):
        self.jones=self.jonesr
        self.jonesmeta=self.jonesrmeta
