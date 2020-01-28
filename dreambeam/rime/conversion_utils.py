"""Functions for converting between astronomical reference frames, etc."""

import warnings
import numpy as np
from casacore.measures import measures
from casacore.quanta import quantity
import math


def basis2basis_transf(basis_from, basis_to):
    """Compute tranform matrix between two bases.

    The basis matrices are assumed to be in the last two axes, and that the
    component vector is summed over the last index (i.e. applied to basis
    matrix from the left).

    Parameters
    ----------
    basis_from: array
        The basis matrices that wishes to transform from. Matrices reside
        in last two axes.
    basis_to: array
        The basis matrices that wishes to transform to. Matrices reside
        in last two axes.

    Returns
    -------
    array:
        The transformation matrix.

    Notes
    -----
    Given some coord sys, and two bases, any vector must be equal in these
    two bases:

    basis_to . v_to =  basis_from . v_from

    where '.' stands from matrix vector multiplication. So the transformation
    from the 'from' system to the 'to' is

    v_to = basis_to^H . basis_from . v_from = T . v_from

    where '^H' stands for hermitian transpose and

    T = basis_to^H . basis_from

    is the transformation matrix that is returned.
    """
    return np.matmul(np.conj(np.swapaxes(basis_to, -2, -1)), basis_from)


def convertBasis(me, rbasis, from_refFrame, to_refFrame):
    basis = np.zeros((3, 3))
    for comp in range(3):
        vr = np.squeeze(rbasis[:, comp])
        (az, el) = crt2sph(vr)
        vr_sph_me = measures().direction(from_refFrame,
                                         quantity(az, 'rad'),
                                         quantity(el, 'rad'))
        v_sph_me = me.measure(vr_sph_me, to_refFrame)
        v_me = sph2crt_me(v_sph_me)
        basis[:, comp] = v_me
    return basis


def CEL2TOPOpnts(obsTimes, stnPos, celPnt):
    # Convert python times to pyrap times
    obsTimes_lst = []
    for obsTime in obsTimes:
        obsTimes_lst.append(quantity(obsTime.isoformat()).get_value())
    obsTimes_me = quantity(obsTimes_lst, 'd')
    # Convert source direction to pyrap
    celPntTheta, celPntPhi, celPntRefFrame = celPnt
    celPnt = celPntRefFrame, str(celPntPhi), str(math.pi/2-celPntTheta)
    celPnt_me = measures().direction(celPnt[0], celPnt[1]+'rad',
                                     celPnt[2]+'rad')
    stnPos_me = measures().position('ITRF', str(stnPos[0, 0])+'m',
                                    str(stnPos[1, 0])+'m',
                                    str(stnPos[2, 0])+'m')
    celPntBasis = getSph2CartTransf(sph2crt_me(celPnt_me))
    obsTimesArr = obsTimes_me.get_value()
    obsTimeUnit = obsTimes_me.get_unit()
    rotang = np.zeros(len(obsTimesArr))
    me = measures()
    # Set position of reference frame w.r.t. ITRF
    me.doframe(stnPos_me)

    me.doframe(me.epoch('UTC', quantity(obsTimesArr[0], obsTimeUnit)))
    CelRot0 = getRotbetweenRefFrames(celPntRefFrame, 'ITRF', me)
    rotang[0] = 0.0

    for ti in range(1, len(obsTimesArr)):
        # Set current time in reference frame
        timEpoch = me.epoch('UTC', quantity(obsTimesArr[ti], obsTimeUnit))
        me.doframe(timEpoch)
        # Incomplete
        CelRot = getRotbetweenRefFrames(celPntRefFrame, 'ITRF', me)
        IncRot = CelRot*CelRot0.T
        rotang[ti] = rotzMat2ang(IncRot)
    return CelRot0, rotang


def getParallacticRot(obsTimes, stnPos, srcDir, doPolPrec=True):
    # Convert python times to pyrap times
    obsTimes_lst = []
    for obsTime in obsTimes:
        obsTimes_lst.append(quantity(obsTime.isoformat()).get_value())
    obsTimes_me = quantity(obsTimes_lst, 'd')
    # Convert source direction to pyrap
    srcTheta, srcPhi, srcRefFrame = srcDir
    srcDir = srcRefFrame, str(srcPhi), str(math.pi/2-srcTheta)
    srcDir_me = measures().direction(srcDir[0],
                                     srcDir[1]+'rad',
                                     srcDir[2]+'rad')
    stnPos_me = measures().position('ITRF',
                                    str(stnPos[0, 0])+'m',
                                    str(stnPos[1, 0])+'m',
                                    str(stnPos[2, 0])+'m')

    obsTimesArr = obsTimes_me.get_value()
    obsTimeUnit = obsTimes_me.get_unit()
    paraMat = np.zeros((len(obsTimesArr), 2, 2))

    me = measures()
    # Set position of reference frame w.r.t. ITRF
    me.doframe(stnPos_me)

    if doPolPrec:
        # Get sky precession rotation matrix
        # (Assuming no change over data interval)
        me.doframe(me.epoch('UTC', quantity(obsTimesArr[0], obsTimeUnit)))
        precMat = getSkyPrecessionMat(me, srcDir_me)
    for ti in range(len(obsTimesArr)):
        # Set current time in reference frame
        timEpoch = me.epoch('UTC', quantity(obsTimesArr[ti], obsTimeUnit))
        me.doframe(timEpoch)

        # Compute polariz comps in spherical sys to cartesian Station coord sys
        # paraMtc=computeParaMat_tc('J2000', 'ITRF', srcDir_me, me)

        # Alternatively:
        paraMme = computeParaMat_me('J2000', 'AZEL', srcDir_me, me)

        paraM = paraMme

        if doPolPrec:
            # With precession:
            paraM = paraM * precMat
            # else:
            # Do not apply precession rotation of polarimetric frame.
            # This is then the apparent polarization frame.

        paraMat[ti, :, :] = paraM
    return paraMat


def getSkyPrecessionMat(me, srcDirection):
    """Compute precession matrix. This is the 2D rotation along srcDirection
    between J2000 and current epoch.

    At J2000 epoch:
    compute 2 cartesian vectors orthogonal to srcDirection: alpha & delta.
    alpha is the orthogonal direction on equator."""
    alphaJ2000ra = srcDirection['m0']['value']+math.pi/2
    alphaJ2000dec = 0.0
    alphaJ2000vec = sph2crt(alphaJ2000ra, alphaJ2000dec)
    # delta is the orthogonal direction along the meridian.
    deltaJ2000ra = srcDirection['m0']['value']
    deltaJ2000dec = srcDirection['m1']['value']+math.pi/2
    deltaJ2000vec = sph2crt(deltaJ2000ra, deltaJ2000dec)

    alpha = me.direction('J2000',
                         str(alphaJ2000ra)+'rad',
                         str(alphaJ2000dec)+'rad')
    delta = me.direction('J2000',
                         str(deltaJ2000ra)+'rad',
                         str(deltaJ2000dec)+'rad')
    # Convert alpha & delta to directions in the current epoch
    alphaTru = me.measure(alpha, 'JTRUE')  # 'JTRUE' isn't stable
    deltaTru = me.measure(delta, 'JTRUE')
    raA = alphaTru['m0']['value']
    decA = alphaTru['m1']['value']
    alphaTruvec = sph2crt(raA, decA)
    raD = deltaTru['m0']['value']
    decD = deltaTru['m1']['value']
    deltaTruvec = sph2crt(raD, decD)

    cosPrecAng = ((alphaTruvec*alphaJ2000vec.T)[0, 0]
                  + (deltaTruvec*deltaJ2000vec.T)[0, 0]) / 2.0
    sinPrecAng = ((alphaTruvec*deltaJ2000vec.T)[0, 0]
                  - (deltaTruvec*alphaJ2000vec.T)[0, 0]) / 2.0
    precMat = np.matrix([[+cosPrecAng, sinPrecAng],
                         [-sinPrecAng, cosPrecAng]])
    return precMat


def getRotbetweenRefFrames(rfFrom, rfTo, me):
    x_fr = me.direction(rfFrom,  '0deg',  '0deg')
    y_fr = me.direction(rfFrom, '90deg',  '0deg')
    z_fr = me.direction(rfFrom, '90deg', '90deg')
    x_to = np.asmatrix(sph2crt_me(me.direction(rfTo,    '0deg',  '0deg')))
    y_to = np.asmatrix(sph2crt_me(me.direction(rfTo,   '90deg',  '0deg')))
    z_to = np.asmatrix(sph2crt_me(me.direction(rfTo,   '90deg', '90deg')))
    x_fr_to = np.asmatrix(sph2crt_me(me.measure(x_fr, rfTo)))
    y_fr_to = np.asmatrix(sph2crt_me(me.measure(y_fr, rfTo)))
    z_fr_to = np.asmatrix(sph2crt_me(me.measure(z_fr, rfTo)))
    xyz_fr_to = np.bmat([[x_fr_to], [y_fr_to], [z_fr_to]])
    xyz_to = np.bmat([[x_to], [y_to], [z_to]])
    frameRot = xyz_fr_to*xyz_to.T
    return frameRot


def computeParaMat_me(rfFrom, rfTo, aDir_me, me):
    """Compute parallactic rotation matrix. me contains epoch"""
    NEvFrom = computeSphBasis(rfFrom, rfTo, 0, aDir_me, me)
    NEvTo = computeSphBasis(rfFrom, rfTo, 1, aDir_me, me)
    paraMat = NEvTo * NEvFrom.T
    return paraMat


def computeSphBasis(rfFrom, rfTo, order, aDir_me, me):
    RbaseSph = me.measure(aDir_me, rfTo)
    Rbase = sph2crt_me(RbaseSph)
    if order == 0:
        NbaseFrom = me.direction(rfFrom,
                                 str(aDir_me['m0']['value'])+'rad',
                                 str(aDir_me['m1']['value']+math.pi/2.0)+'rad')
        NbaseSph = me.measure(NbaseFrom, rfTo)
    else:
        NbaseSph = me.direction(rfTo,
                                str(RbaseSph['m0']['value'])+'rad',
                                str(RbaseSph['m1']['value']+math.pi/2.0)+'rad')
    Nbase = sph2crt_me(NbaseSph)
    Ebase = np.cross(Nbase, Rbase)
    NEbasis = np.bmat([[Nbase], [Ebase]])
    return NEbasis


def computeParaMat_tc(rfFrom, rfTo, aDir_me, me):
    framRot = getRotbetweenRefFrames(rfFrom, rfTo, me)
    srcFromcrt = sph2crt_me(me.measure(aDir_me, rfFrom))
    pol2cartFrom = computeSph2CrtMat(srcFromcrt)
    srcTocrt = sph2crt_me(me.measure(aDir_me, rfTo))
    pol2cartTo = computeSph2CrtMat(srcTocrt)
    paraMat = pol2cartFrom.T*framRot*pol2cartTo
    return paraMat


def getSph2CartTransf(r):
    """Compute the transformation matrix from a spherical basis to a Cartesian
    basis at the field point given by the input 'r'.
    The output 'transf_sph2cart' is defined such that:

    [[v_x], [v_y], [v_z]]=transf_sph2cart*matrix([[v_r], [v_phi], [v_theta]]).
    """
    r = np.matrix(r)
    r = r.squeeze()
    ru = r/np.sqrt((r*r.T)[0, 0])
    xu = ru[0, 0]
    yu = ru[0, 1]
    zu = ru[0, 2]
    tang2cart = 1.0/np.sqrt(xu*xu+yu*yu) * np.matrix(
                        [[-yu,        -xu*zu],
                         [ xu,        -yu*zu],
                         [ 0., (xu*xu+yu*yu)]])
    transf_sph2cart = np.bmat([[ru.T, tang2cart]])
    return transf_sph2cart


def getSph2CartTransfArr(r):
    """(Array version) Compute the transformation matrix from a spherical basis
    to a Cartesian basis at the field point given by the input 'r'.
    The output 'transf_sph2cart' is defined such that:

    [[v_x], [v_y], [v_z]]=transf_sph2cart*matrix([[v_r], [v_phi], [v_theta]]).
    """
    r = np.array(r)
    r = r.squeeze()
    # ru = r/np.sqrt(np.tensordot(r,r, axes=([0],[0])))
    ru = r/np.sqrt(r[0, ...]**2+r[1, ...]**2+r[2, ...]**2)
    xu = ru[0, ...]
    yu = ru[1, ...]
    zu = ru[2, ...]
    nrf = 1.0/np.sqrt(xu*xu+yu*yu)
    transf_sph2cart = np.array([
                          [xu,            -yu*nrf,        -xu*zu*nrf],
                          [yu,             xu*nrf,        -yu*zu*nrf],
                          [zu, np.zeros(xu.shape), (xu*xu+yu*yu)*nrf]])
    return transf_sph2cart


def computeSph2CrtMat(lmnMatrix):
    lmn = lmnMatrix.squeeze()
    l = lmn[0, 0]
    m = lmn[0, 1]
    n = lmn[0, 2]
    polz2cart = 1.0/np.sqrt(l*l+m*m) * np.matrix(
                        [[ -m,      -l*n],
                         [  l,      -m*n],
                         [ .0, (l*l+m*m)]])
    return polz2cart


def rotzMat2ang(rotMat):
    ang = np.arctan2(rotMat[0, 1], rotMat[0, 0])
    return ang


def crt2sph(dir_crt):
    """Cartesian to spherical conversion.

N.B. When using direction cosines (Db uses them), regions outside the horizon
(such as l,m=1,1) will have an imaginary vertical component. This function will
neglect the imaginary part of the vertical component, so the elevation will be
0, rather than raise an exception. This is so that further processing can
proceed. The user should check the vertical component of the radial base of the
basis attribute, i.e. jones.basis[...,0,2], to see if it is imaginary and act
accordingly.
    """
    x = np.squeeze(dir_crt[0])
    y = np.squeeze(dir_crt[1])
    z = np.squeeze(dir_crt[2])
    # x,y should be real but since z could be imaginary, dir_crt may have
    # dtype complex, so x,y will seem complex. Therefore take the real parts
    # of x and y here:
    azi = np.arctan2(np.real(y), np.real(x))
    ele = np.arcsin(np.real(z))
    if np.any(np.imag(z)):
         warnings.warn("Some z-components are imaginary.")
    # Radial component (not needed for unit directions)
    # r=np.sqrt(x**2+y**2+z**2)
    # ele=np.arcsin(z/r)
    return azi, ele


def sph2crt_me(sphme):
    return sph2crt(sphme['m0']['value'], sphme['m1']['value'])


def sph2crt(azi, ele):
    """Converts spherical angles to cartesian coordinates. The spherical
    angles are azimuth and elevation.
    """
    x = np.cos(ele)*np.cos(azi)
    y = np.cos(ele)*np.sin(azi)
    z = np.sin(ele)
    return(np.array([x, y, z]))


# Here C09 refers to Carozzi2009 vCZ paper: it has
#   (r_hat, phi_hat, theta_hat).
# While IAU has:
#   +x pointing to +vartheta_hat, +y along +phi_hat, +z along -r_hat.
# HOWEVER, implementation here keeps r_hat in first column, so:
IAUtoC09 = np.array([[1., 0., 0.],
                     [0., 0., 1.],
                     [0., 1., 0.]])
C09toIAU = np.linalg.inv(IAUtoC09)
# C09toIAU = np.array([[1., 0., 0.],
#                      [0., 0., 1.],
#                      [0., 1., 0.]])


def IAU_pol_basis(src_az, src_el):
    """Compute the (x_hat, y_hat, z_hat) basis in IAU polarization system for a
    direction given by (azimuth, elevation) tuple typically Ra, Dec."""

    # IAUtoC09 = np.identity(3)
    basis_C09 = np.array(getSph2CartTransf(sph2crt(src_az, src_el)))
    # basis_IAU = np.matmul(basis_C09, IAUtoC09)
    basis_C09 = np.asmatrix(basis_C09)
    return basis_C09


def sphmeshgrid(nr_ele=128, nr_azi=256):
    """Provides a polar angles grid on a 2-sphere. Azimuthal dimension should
    not wrap around to ensure uniqueness on sphere. The mesh has theta
    increasing along rows and phi along columns."""
    # Include endpoint:
    ele = np.linspace(math.pi/2, -math.pi/2, nr_ele, endpoint=True)
    # Don't overlap 0 and 2*pi
    azi = np.linspace(0.0, 2*math.pi, nr_azi, endpoint=False)
    azimsh, elemsh = np.meshgrid(azi, ele)
    return azimsh, elemsh


def dc_hrz2vrt(ll, mm, hemisph=+1):
    """Compute vertical component of an array of 2D horizontal direction
    cosines."""
    nn = hemisph * np.sqrt((1-ll**2-mm**2).astype(complex))
    return nn


def pyTimes2meTimes(pyTimes):
    obsTimes_lst = []
    for obsTime in pyTimes:
        obsTimes_lst.append(quantity(obsTime.isoformat()).get_value())
    obsTimes_me = quantity(obsTimes_lst, 'd')
    obsTimesArr = obsTimes_me.get_value()
    obsTimeUnit = obsTimes_me.get_unit()
    return obsTimesArr, obsTimeUnit


def setEpoch(obsTimesArr, obsTimeUnit):
    stnPos_me = measures().position('ITRF', '0m', '0m', '0m')
    me = measures()
    # Set position of reference frame w.r.t. ITRF
    me.doframe(stnPos_me)
    timEpoch = me.epoch('UTC', quantity(obsTimesArr, obsTimeUnit))
    me.doframe(timEpoch)
    return me


def printJones(Jn):
    for ti in range(0, Jn.shape[0]):
        print(Jn[ti, :, :])


def shiftmat2back(arr):
    arr = np.rollaxis(arr, 0, arr.ndim)
    arr = np.rollaxis(arr, 0, arr.ndim)
    return arr
