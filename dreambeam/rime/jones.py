"""
  This module provides a Jones matrix framework for radio interometric
  measurement equations.
"""
import copy
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from casacore.measures import measures
from casacore.quanta import quantity
from conversionUtils import sph2crt, crt2sph, convertBasis, \
                            getSph2CartTransf, getSph2CartTransfArr, \
                            IAU_pol_basis, shiftmat2back, IAUtoC09, \
                            sphmeshgrid, dc_hrz2vrt


class Jones(object):
    """This is the base class for Jones algebra. It contains the Jones matrix
    itself and a basis w.r.t. which the Jones matrix is given.
    The basis is such that:
        self.jonesbasis=array([[r_hat], [phi_hat], [theta_hat]]).
    """
    _ecef_frame = 'ITRF'
    _eci_frame = 'J2000'
    _topo_frame = 'STN'

    def __init__(self):
        pass

    def op(self, jonesobjright):
        """Operate this Jones on to the Jones passed in the argument."""
        self.jonesr = jonesobjright.getValue()
        self.jonesrbasis_from = jonesobjright.get_basis()
        self.refframe_r = jonesobjright.get_refframe()
        self.iaucmp = jonesobjright.iaucmp
        self.computeJonesRes()
        return self

    def getValue(self):
        """Return value of the Jones matrix"""
        return self.jones

    def get_basis(self):
        """Return basis of the Jones matrix"""
        return self.jonesbasis

    def get_refframe(self):
        """Return the reference frame of the Jones matrix."""
        return self.refframe

    def computeJonesRes(self):
        pass

    def sph2lud3_basis(self, jonesbasis_sph, alignment=None):
        """Convert sph basis to Ludwig3 frame with an optional rotation
        alignment."""
        # The jonesbasis for the antennas is taken to be the Ludwig3 def.
        # with r,u,v basis expressed wrt the station frame
        r_refframe = jonesbasis_sph[..., 0]
        if alignment is not None:
            r = np.tensordot(r_refframe, alignment, axes=([-1, 1]))
        else:
            r = r_refframe
        (az, el) = crt2sph(r.T)
        lugwig3rot = np.zeros((3, 3, len(az)))
        lugwig3rot[0, 0, :] = 1.
        lugwig3rot[1:, 1:, :] = np.array([[np.cos(az), np.sin(az)],
                                          [-np.sin(az), np.cos(az)]])
        lugwig3rot = np.moveaxis(lugwig3rot, -1, 0)
        jonesbasis_lud3 = np.matmul(jonesbasis_sph, lugwig3rot)
        # ang_u = np.rad2deg(
        #           np.arctan2(jonesbasis_lud3[:,1,1], jonesbasis_lud3[:,0,1]))
        # print ang_u
        return jonesbasis_lud3

    def convert2iaucmp(self):
        if not self.iaucmp:
            self.jones = np.matmul(self.jones, IAUtoC09[1:, 1:])
            self.iaucmp = True


class JonesChain(object):
    jonesproducts = []

    def __init__(self):
        self.jonesproducts = []


class PJones(Jones):
    """This is a P-Jones or parallactic Jones. This has a temporal dependence
    given by the epoch of observation."""

    def __init__(self, obsTimespy, ITRF2stnrot, do_parallactic_rot=True):
        super(PJones, self).__init__()
        obsTimes_lst = []
        for obsTimepy in obsTimespy:
            obsTimes_lst.append(quantity(obsTimepy.isoformat()).get_value())
        obsTimes_me = quantity(obsTimes_lst, 'd')
        self.obsTimes = obsTimes_me.get_value()
        self.obsTimeUnit = obsTimes_me.get_unit()
        self.ITRF2stnrot = ITRF2stnrot
        self.do_parallactic_rot = do_parallactic_rot

    def computeJonesRes(self):
        if type(self.obsTimes) is float:
            self.computeJonesRes_overfield()
        else:
            self.computeJonesRes_overtime()

    def computeJonesRes_overtime(self):
        """Compute and apply the P Jones matrix. The structure is:

          jones[time, sphcomp, skycomp] =
                Pjones[time, sphcomp, comp]*jonesr[comp, skycomp]

        The Pjones matrix is computed as follows: consider a direction
        vector d. Let jonesrbasis be the column concatenation of the 3
        spherical basis vectors corresponding to d in the J2000 reference
        frame, so
          jonesrbasis = [[r_J2000],[phi_J2000],[theta_J2000]].T
        where r_J2000 is along the direction d and theta, phi are the remaining
        two spherical basis vectors. Let jonesbasis be the basis vectors
        corresponding to the same direction d but in the STN reference frame,
        so
          jonesbasis = [[r_STN],[phi_STN],[theta_STN]].T
        where r_STN is along the direction d and theta, phi are the remaining
        two spherical basis vectors in the spherical system associated with the
        STN.

        The algorithm takes r_J2000 from component 0 of the jonesrbasis and
        converts it to STN (i.e. finds r_STN) using casacore measures module,
        along with the other 2 J2000 basis vectors. These converted vectors are
        called jonesrbasis_to. With r_STN, it also computes the corresponding
        jonesbasis. A vector in the cartesian J2000 ref sys converted to STN
        must be equal to the same vector expressed in the cartesian STN ref sys
        via a coversion from spherical, so
          jonesbasis * V_STN^sph = jonesrbasis_to * V_J2000^sph
        which implies that we can convert directly from spherical J2000 to the
        sppherical STN like this
          V_STN^sph = (jonesbasis.H * jonesrbasis_to) * V_J2000^sph
        where the matrix in parentheses is the P Jones matrix.

        The P Jones matrix is then applied to the operand Jones matrix.
        """
        nrOfTimes = len(self.obsTimes)
        pjones = np.zeros((nrOfTimes, 2, 2))
        me = measures()
        me.doframe(measures().position(self._ecef_frame, '0m', '0m', '0m'))
        self.jonesbasis = np.zeros((nrOfTimes, 3, 3))
        if self.refframe_r == self._eci_frame:
            convert2irf = self._ecef_frame
            jonesrbasis_from = self.jonesrbasis_from
            jr_refframe = self.refframe_r
        else:
            convert2irf = self._eci_frame
            jonesrbasis_from = np.matmul(self.ITRF2stnrot.T,
                                         self.jonesrbasis_from)
            jr_refframe = self._ecef_frame
        for ti in range(0, nrOfTimes):
            # Set current time in reference frame
            timEpoch = me.epoch('UTC', quantity(self.obsTimes[ti],
                                                self.obsTimeUnit))
            me.doframe(timEpoch)
            jonesrbasis_to = np.asmatrix(convertBasis(me, jonesrbasis_from,
                                                      jr_refframe,
                                                      convert2irf))
            if convert2irf == self._ecef_frame:
                jonesrbasis_to = np.matmul(self.ITRF2stnrot, jonesrbasis_to)
            jonesbasisMat = getSph2CartTransf(jonesrbasis_to[:, 0])
            if self.do_parallactic_rot:
                pjones[ti, :, :] = jonesbasisMat[:, 1:].H \
                                    * jonesrbasis_to[:, 1:]
            else:
                pjones[ti, :, :] = np.asmatrix(np.identity(2))
            self.jonesbasis[ti, :, :] = jonesbasisMat
        if convert2irf == self._ecef_frame:
            self.refframe = 'STN'  # Final Ref frame is station
        else:
            self.refframe = self._eci_frame
        self.jones = np.matmul(pjones, self.jonesr)
        self.thisjones = pjones

    def computeJonesRes_overfield(self):
        """Compute the PJones over field of directions for one frequency.
        """
        pjones = np.zeros(self.jonesrbasis_from.shape[0:-2]+(2, 2))
        me = measures()
        me.doframe(measures().position(self._ecef_frame, '0m', '0m', '0m'))
        self.jonesbasis = np.zeros(self.jonesrbasis_from.shape)
        if self.refframe_r == self._eci_frame:
            convert2irf = self._ecef_frame
            jonesrbasis_from = self.jonesrbasis_from
            jr_refframe = self.refframe_r
        else:
            convert2irf = self._eci_frame
            jonesrbasis_from = np.matmul(self.ITRF2stnrot.T,
                                         self.jonesrbasis_from)
            jr_refframe = self._ecef_frame
        timEpoch = me.epoch('UTC', quantity(self.obsTimes, self.obsTimeUnit))
        me.doframe(timEpoch)
        for idxi in range(self.jonesrbasis_from.shape[0]):
            for idxj in range(self.jonesrbasis_from.shape[1]):
                jonesrbasis_to = np.asmatrix(convertBasis(
                                            me,
                                            jonesrbasis_from[idxi, idxj, :, :],
                                            jr_refframe, convert2irf))
                if convert2irf == self._ecef_frame:
                    jonesrbasis_to = np.matmul(self.ITRF2stnrot,
                                               jonesrbasis_to)
                jonesbasisMat = getSph2CartTransf(jonesrbasis_to[..., 0])
                pjones[idxi, idxj, :, :] = jonesbasisMat[:, 1:].H \
                    * jonesrbasis_to[:, 1:]
                self.jonesbasis[idxi, idxj, :, :] = jonesbasisMat
        if convert2irf == self._ecef_frame:
            self.refframe = 'STN'  # Final Ref frame is station
        else:
            self.refframe = self._eci_frame
        self.jones = np.matmul(pjones, self.jonesr)
        self.thisjones = pjones


class DualPolFieldPointSrc(Jones):
    """This is a mock Jones point source. It does not model a real source. It's
    purpose is for testing. It can be seen as a source that first transmits in
    one polarization and then in another, then 2 transmissions given in the 2
    columns.
    It may have a spectral dimension. The src_dir should be a tuple with
    (az, el, ref)."""

    def __init__(self, src_dir, dualPolField=np.identity(2), iaucmp=True):
        (src_az, src_el, src_ref) = src_dir
        dualPolField3d = np.asmatrix(np.identity(3))
        dualPolField3d[1:, 1:] = np.asmatrix(dualPolField)
        if iaucmp:
            jones = np.matmul(IAUtoC09, dualPolField3d)[1:, 1:]
            self.iaucmp = True
        else:
            jones = dualPolField3d[1:, 1:]
            self.iaucmp = False
        self.jones = np.asarray(jones)
        self.jonesbasis = np.asarray(IAU_pol_basis(src_az, src_el))
        self.refframe = src_ref


class DualPolFieldRegion(Jones):
    """This is a Jones unit flux density field."""

    def __init__(self, refframe='J2000', dualPolField=np.identity(2),
                 iaucmp=True, lmgrid=None):
        if not lmgrid:
            azimsh, elemsh = sphmeshgrid()
            lmn = sph2crt(azimsh, elemsh)
        else:
            nn = dc_hrz2vrt(*lmgrid)
            lmn = np.array(lmgrid+(nn,))
            azimsh, elemsh = crt2sph(lmn)
        self.azmsh = azimsh
        self.elmsh = elemsh
        dualPolField3d = np.asmatrix(np.identity(3))
        dualPolField3d[1:, 1:] = np.asmatrix(dualPolField)
        if iaucmp:
            jones = np.matmul(IAUtoC09, dualPolField3d)[1:, 1:]
            self.iaucmp = True
        else:
            jones = dualPolField3d[1:, 1:]
            self.iaucmp = False
        self.jones = np.broadcast_to(jones,
                                     elemsh.shape+dualPolField.shape)
        self.jonesbasis = shiftmat2back(getSph2CartTransfArr(lmn))
        self.refframe = refframe


class EJones(Jones):
    """This is the antenna or feed Jones. It is given by a set of complex gain
    patterns for each frequency and polarization channel."""

    def __init__(self, dualPolElem, position, stnRot, freqSel=None):
        self.position = position
        self.stnRot = stnRot
        self.dualPolElem = dualPolElem
        if freqSel is None:
            self.freqChan = self.dualPolElem.getfreqs()
        else:
            self.freqChan = freqSel
        self.refframe = 'STN'

    def computeJonesRes(self):
        """Compute the Jones that results from applying the E-Jones to the
        right.
        The structure of the jonesrbasis is [timeIdx, sphIdx, skycompIdx].
        """
        idxshape = self.jonesrbasis_from.shape[0:-2]
        jonesrbasis = np.reshape(self.jonesrbasis_from, (-1, 3, 3))
        (az_from, el_from) = crt2sph(jonesrbasis[..., 0].squeeze().T)
        theta_phi_view = (np.pi/2-el_from.flatten(), az_from.flatten())
        ejones = self.dualPolElem.getJonesAlong(self.freqChan, theta_phi_view)
        # AntPat has column order theta_hat, phi_hat
        # while C09 has phi_hat, vartheta_hat (where vartheta_hat=-theta_hat).
        # So flip order and change sign of theta_hat.
        ejones = ejones[..., ::-1]
        ejones[..., 1] = -ejones[..., 1]

        self.jonesbasis = self.jonesrbasis_from  # Basis does not change
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
        self.refframe = self.refframe_r


def inverse(jonesobj):
    """Return a Jones object that is the inverse of jonesobj."""
    inv_jones = copy.copy(jonesobj)
    jmat = jonesobj.getValue()
    inv_jones.jones = np.linalg.inv(jmat)
    # Swap basis between left and right:
    inv_jones.jonesbasis = jonesobj.jonesrbasis_from
    inv_jones.jonesrbasis_from = jonesobj.jonesbasis
    jframe = jonesobj.get_refframe()
    jrframe = jonesobj.refframe_r
    inv_jones.refframe = jrframe
    inv_jones.refframe_r = jframe
    return inv_jones


def fix_imaginary_directions(jonesobj, fill=np.identity(2)):
    """Replace jones matrices with imaginary directions in a Jones object.

When specifying 2D Cartesian direction cosines, it is possible that the
corresponding direction is not physical, e.g. when l,m = 1,1. In such cases,
the Jones radius basis will have an imaginary vertical component. This function
will find such 'directions' and replace the corresponding Jones matrix will the
fill matrix specified by the `fill` argument.
    """
    idxs = np.where(np.imag(jonesobj.jonesbasis[..., 0, 2]))
    jonesobj.jones[idxs[0], idxs[1], ...] = fill


def plotJonesField(az, el, jonesfld, jbasis, refframe, rep='abs-Jones'):
    """Plot a Jones field."""
    if rep == 'abs-Jones':
        restitle = 'Beam Jones on sky'
        res00 = np.abs(jonesfld[:, :, 0, 0])
        res00 = ma.masked_invalid(res00)
        res00lbl = r'|J_{p\phi}|'
        res01 = np.abs(jonesfld[:, :, 0, 1])
        res01 = ma.masked_invalid(res01)
        res01lbl = r'|J_{p\theta}|'
        res10 = np.abs(jonesfld[:, :, 1, 0])
        res10 = ma.masked_invalid(res10)
        res10lbl = r'|J_{q\phi}|'
        res11 = np.abs(jonesfld[:, :, 1, 1])
        res11 = ma.masked_invalid(res11)
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
        raise Exception("Unknown Jones representation {}.".format(rep))
    if refframe == 'STN':
        # Directions in Cartesian station crds
        x = jbasis[..., 0, 0]
        y = jbasis[..., 1, 0]
        z = jbasis[..., 2, 0]
        belowhrz = ma.getmask(ma.masked_less(z, .0))
        x = ma.MaskedArray(x, mask=belowhrz)
        y = ma.MaskedArray(y, mask=belowhrz)
        res00 = ma.MaskedArray(res00, mask=belowhrz)
        res01 = ma.MaskedArray(res01, mask=belowhrz)
        res10 = ma.MaskedArray(res10, mask=belowhrz)
        res11 = ma.MaskedArray(res11, mask=belowhrz)
        xlabel = 'STN X'
        ylabel = 'STN Y'
    elif refframe == 'J2000':
        x = az
        y = el
        xlabel = 'RA'
        ylabel = 'DEC'
    fig = plt.figure()
    fig.suptitle(restitle)
    ax = plt.subplot(221, polar=False)
    plt.pcolormesh(x, y, res00, vmin=0., vmax=2.0)
    plt.colorbar()
    ax.set_title(res00lbl)
    plt.ylabel(ylabel)

    ax = plt.subplot(222, polar=False)
    plt.pcolormesh(x, y, res01)
    plt.colorbar()
    ax.set_title(res01lbl)

    ax = plt.subplot(223, polar=False)
    plt.pcolormesh(x, y, res10)
    plt.colorbar()
    ax.set_title(res10lbl)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    ax = plt.subplot(224, polar=False)
    plt.pcolormesh(x, y, res11, vmin=np.nanmin(res11), vmax=np.nanmax(res11))
    plt.colorbar()
    ax.set_title(res11lbl)
    plt.xlabel(xlabel)

    plt.show()
