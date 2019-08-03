"""
  This module provides a Jones matrix framework for radio interometric
  measurement equations.
"""
#import math
import numpy as np
import matplotlib.pyplot as plt
from casacore.measures import measures
from casacore.quanta import quantity
from conversionUtils import sph2crt, crt2sph, computeSphBasis, convertBasis, \
                            getSph2CartTransf, getSph2CartTransfArr, \
                            IAU_pol_basis, shiftmat2back
from antpat.reps.sphgridfun.tvecfun import getSph2CartTransfMat
from antpat.reps.sphgridfun.pntsonsphere import sphericalGrid


class Jones(object):
    """This is the base class for Jones algebra. It contains the Jones matrix
    itself and a basis w.r.t. which the Jones matrix is given.
    The basis is such that:
        self.jonesbasis=array([[r_hat], [phi_hat], [theta_hat]]).
    """

    def __init__(self):
        pass

    def op(self, jonesobjright):
        """Operate this Jones on to the Jones passed in the argument."""
        self.jonesr = jonesobjright.getValue()
        self.jonesrbasis = jonesobjright.getBasis()
        self.jonesrmeta = jonesobjright.getMetadata()
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
    jonesproducts = []

    def __init__(self):
        self.jonesproducts = []


class PJones(Jones):
    """This is a P-Jones or parallactic Jones. This has a temporal dependence
    given by the epoch of observation."""

    def __init__(self, obsTimespy, ITRF2stnrot):
        super(PJones, self).__init__()
        obsTimes_lst = []
        for obsTimepy in obsTimespy:
            obsTimes_lst.append(quantity(obsTimepy.isoformat()).get_value())
        obsTimes_me = quantity(obsTimes_lst, 'd')
        self.obsTimes = obsTimes_me.get_value()
        self.obsTimeUnit = obsTimes_me.get_unit()
        self.jonesmeta = {}
        self.jonesmeta['refFrame'] = 'ITRF'
        self.ITRF2stnrot = ITRF2stnrot

    def computeJonesRes(self):
        if type(self.obsTimes) is float:
            self.computeJonesRes_overfield()
        else:
            self.computeJonesRes_overtime()

    def computeJonesRes_overtime(self):
        """Compute the resulting Jones matrix when the parallactic rotation
        matrix is applied to a source brightness. The structure is:

         jones[ti, sphcompIdx, skycompIdx] =
             paraRot[timeIdx,sphcompIdx,compIdx]*jonesr[compIdx,skycompIdx]

        """
        nrOfTimes = len(self.obsTimes)
        paraRot = np.zeros((nrOfTimes, 2, 2))
        me = measures()
        me.doframe(measures().position('ITRF', '0m', '0m', '0m'))
        self.jonesbasis = np.zeros((nrOfTimes, 3, 3))
        (az_from, el_from) = crt2sph(self.jonesrbasis[:, 0])
        # r_sph_me = measures().direction(self.jonesrmeta['refFrame'],
        #                                 quantity(az_from, 'rad'),
        #                                 quantity(el_from, 'rad'))
        for ti in range(0, nrOfTimes):
            # Set current time in reference frame
            timEpoch = me.epoch('UTC', quantity(self.obsTimes[ti],
                                                self.obsTimeUnit))
            me.doframe(timEpoch)
            # paraRot[ti,:,:]=computeParaMat_me(self.jonesrmeta['refFrame'],
            #                         self.jonesmeta['refFrame'], r_sph_me, me)

            jonesrbasis_to = np.asmatrix(convertBasis(me, self.jonesrbasis,
                                                   self.jonesrmeta['refFrame'],
                                                   self.jonesmeta['refFrame']))
            jonesrbasis_to = np.matmul(self.ITRF2stnrot, jonesrbasis_to)
            # jonesrbasis_to2 = computeSphBasis(self.jonesrmeta['refFrame'],
            #                                   self.jonesmeta['refFrame'],
            #                                   0, r_sph_me, me)
            jonesbasisMat = getSph2CartTransf(jonesrbasis_to[:, 0])
            #print("to", jonesrbasis_to)
            paraRot[ti, :, :] = jonesbasisMat[:, 1:].H*jonesrbasis_to[:, 1:]
            # paraRot[ti,:,:]=jonesrbasis_to2*jonesbasisMat[:,1:]
            self.jonesbasis[ti, :, :] = jonesbasisMat
        self.jones = np.matmul(paraRot, self.jonesr)
        self.thisjones = paraRot

    def computeJonesRes_overfield(self):
        paraRot = np.zeros(self.jonesrbasis.shape[0:-2]+(2, 2))
        me = measures()
        me.doframe(measures().position('ITRF', '0m', '0m', '0m'))
        self.jonesbasis = np.zeros(self.jonesrbasis.shape)
        timEpoch = me.epoch('UTC', quantity(self.obsTimes, self.obsTimeUnit))
        me.doframe(timEpoch)
        for idxi in range(self.jonesrbasis.shape[0]):
            for idxj in range(self.jonesrbasis.shape[1]):
                jonesrbasis_to = np.asmatrix(convertBasis(me,
                                            self.jonesrbasis[idxi, idxj, :, :],
                                                  self.jonesrmeta['refFrame'],
                                                  self.jonesmeta['refFrame']))
                jonesrbasis_to = np.matmul(self.ITRF2stnrot, jonesrbasis_to)
                jonesbasisMat = getSph2CartTransf(jonesrbasis_to[..., 0])
                paraRot[idxi, idxj, :, :] = jonesbasisMat[:, 1:].H*jonesrbasis_to[:, 1:]
                self.jonesbasis[idxi, idxj, :, :] = jonesbasisMat
        self.jones = np.matmul(paraRot, self.jonesr)
        self.thisjones = paraRot


class DualPolFieldPointSrc(Jones):
    """This is a mock Jones point source. It does not model a real source. It's
    purpose is for testing. It can be seen as a source that first transmits in
    one polarization and then in another, then 2 transmissions given in the 2
    columns.
    It may have a spectral dimension. The src_dir should be a tuple with
    (az, el, ref)."""

    def __init__(self, src_dir, dualPolField=np.identity(2)):
        (src_az, src_el, src_ref) = src_dir
        self.jones = dualPolField
        #self.jonesbasis = np.array(getSph2CartTransf(sph2crt(src_az, src_el)))
        self.jonesbasis = IAU_pol_basis(src_az, src_el)
        self.jonesmeta = {}
        self.jonesmeta['refFrame'] = src_ref


class DualPolFieldRegion(Jones):
    """This is a Jones unit flux density field."""

    def __init__(self, dualPolField=np.identity(2)):
        #(src_az0, src_el0, src_ref) = src_dir
        thetamsh, phimsh = sphericalGrid()
        self.elmsh = np.pi/2-thetamsh
        self.azmsh = phimsh
        self.jones = np.broadcast_to(dualPolField,
                                     self.elmsh.shape+dualPolField.shape)
        self.jonesbasis = shiftmat2back(getSph2CartTransfArr(sph2crt(self.azmsh,
                                                                  self.elmsh)))
        #self.jonesbasis = IAU_pol_basis(azmsh, elmsh)
        self.jonesmeta = {'refFrame': 'J2000'}


class EJones(Jones):
    """This is the antenna or feed Jones. It is given by a set of complex gain
    patterns for each frequency and polarization channel."""

    def __init__(self, dualPolElem, position, stnRot, freqSel=None):
        self.position = position
        self.stnRot = stnRot
        self.dualPolElem = dualPolElem
        # self.dualPolElem.rotateframe(np.asarray(stnRot))
        if freqSel is None:
            self.freqChan = self.dualPolElem.getfreqs()
        else:
            self.freqChan = freqSel

    def computeJonesRes(self):
        """Compute the Jones that results from applying the E-Jones to the
        right.
        The structure of the jonesrbasis is [timeIdx, sphIdx, skycompIdx].
        """
        idxshape = self.jonesrbasis.shape[0:-2]
        jonesrbasis = np.reshape(self.jonesrbasis, (-1, 3, 3))
        #jonesrbasis_to = np.matmul(np.asarray(self.stnRot.T), jonesrbasis)
        jonesrbasis_to = jonesrbasis
        (az_from, el_from) = crt2sph(jonesrbasis[..., 0].squeeze().T)
        theta_phi_view = (np.pi/2-el_from.flatten(), az_from.flatten())
        ejones = self.dualPolElem.getJonesAlong(self.freqChan, theta_phi_view)
        #(theta_lcl, phi_lcl) = self.dualPolElem.getBuildCoordinates(math.pi/2-r_sph[1], r_sph[0])
        #print theta_lcl, phi_lcl
        r_lcl = crt2sph(jonesrbasis_to[..., 0].squeeze().T)
        #print np.rad2deg(r_lcl)
        jonesbasisMat = getSph2CartTransfMat(jonesrbasis_to[..., 0].squeeze())
        #paraRot = np.matmul(np.conjugate(jonesbasisMat), jonesrbasis_to)
        self.jonesbasis = np.reshape(jonesbasisMat,
                                     idxshape+jonesbasisMat.shape[1:])
        # This is the actual MEq multiplication:
        if ejones.ndim > 3:
            frqdimsz = (ejones.shape[0],)
        else:
            frqdimsz = ()
        self.jones = np.reshape(
                        np.matmul(ejones, np.reshape(self.jonesr, (-1, 2, 2))),
                        frqdimsz+idxshape+(2, 2)
                        )
        self.thisjones = np.reshape(ejones, frqdimsz+idxshape+(2, 2))

    def getPosRot(self, position):
        """Compute the nominal transformation from the geodetic position to
        ITRF. (Not implemented yet)"""
        return np.identity(2)


class DualPolFieldSink(Jones):
    def computeJonesRes(self):
        self.jones = self.jonesr
        self.jonesmeta = self.jonesrmeta


def plotJonesField(az, el, jonesfld, rep='abs-Jones'):
    if rep == 'abs-Jones':
        restitle = 'Beam Jones on sky'
        res00 = np.abs(jonesfld[:, :, 0, 0])
        res00 = np.ma.masked_invalid(res00)
        res00lbl = r'|J_{p\phi}|'
        res01 = np.abs(jonesfld[:, :, 0, 1])
        res01 = np.ma.masked_invalid(res01)
        res01lbl = r'|J_{p\theta}|'
        res10 = np.abs(jonesfld[:, :, 1, 0])
        res10 = np.ma.masked_invalid(res10)
        res10lbl = r'|J_{q\phi}|'
        res11 = np.abs(jonesfld[:, :, 1, 1])
        res11 = np.ma.masked_invalid(res11)
        res11lbl = r'|J_{q\theta}|'
    elif rep == 'Stokes':
        corrmat = np.matmul(jonesfld, np.swapaxes(jonesfld.conj(), -2, -1))
        S0 = np.real(corrmat[:, :, 0, 0]+corrmat[:, :, 1, 1])
        SQ = np.real(corrmat[:, :, 0, 0]-corrmat[:, :, 1, 1])
        SU = np.real(corrmat[:, :, 0, 1]+corrmat[:, :, 1, 0])
        SV = np.imag(corrmat[:, :, 0, 1]-corrmat[:, :, 1, 0])
        restitle = 'Antenna Stokes on sky'
        res00 = S0
        res00lbl = 'I'
        res01 = SQ/S0
        res01lbl = 'q'
        res10 = SU/S0
        res10lbl = 'u'
        res11 = SV/S0
        res11lbl = 'v'
    else:
        print "Unknown Jones representation."
        exit(1)

    fig = plt.figure()
    fig.suptitle(restitle)
    ax = plt.subplot(221, polar=False)
    plt.pcolormesh(az, el, res00, vmin=0., vmax=2.0)
    plt.colorbar()
    ax.set_title(res00lbl)
    plt.ylabel('DEC')

    ax = plt.subplot(222, polar=False)
    plt.pcolormesh(az, el, res01, vmin=-1., vmax=1.)
    plt.colorbar()
    ax.set_title(res01lbl)

    ax = plt.subplot(223, polar=False)
    plt.pcolormesh(az, el, res10, vmin=-1., vmax=1.)
    plt.colorbar()
    ax.set_title(res10lbl)
    plt.xlabel('RA')
    plt.ylabel('DEC')

    ax = plt.subplot(224, polar=False)
    plt.pcolormesh(az, el, res11, vmin=-1., vmax=1.)
    plt.colorbar()
    ax.set_title(res11lbl)
    plt.xlabel('RA')

    plt.show()
