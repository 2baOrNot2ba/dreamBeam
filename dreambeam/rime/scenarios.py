'''
This module provides some common observational scenarios of interest.
It implements basic Jones chains or Measurement Equations.
'''
import numpy as np
import dreambeam.rime.jones
from dreambeam.telescopes.rt import TelescopesWiz
from dreambeam.rime.conversion_utils import basis2basis_transf, C09toIAU


def on_pointing_axis_tracking(telescopename, stnid, band, antmodel, obstimebeg,
                              obsdur, obstimestp, pointingdir,
                              do_parallactic_rot=True):
    """Computes the Jones matrix along pointing axis while tracking a fixed
    celestial source.

    This function computes the Jones of an observational setup where a
    telescope station is tracking a source (on-axis) on the sky using one of
    its bands.

    The Jones matrix computed is the matrix which maps the two transverse
    components of the E-field vector (V/m or unitless) at a given frequency
    propagating from the pointing direction to the telescope stations two
    polarization channels. The x,y components of the source are as specified by
    the IAU for polarized emissions; and the output components are the two
    ordered polarization channels. The matrix inverse of the output Jones
    matrix multiplied from the left on the channelized 2D voltage data will
    produce an estimate of the E-field vector of the source in the IAU for all
    the frequencies in the band over the times of the tracking.

    The response of the dual-polarized feed is modeled using the `antmodel`
    specified.

    Parameters
    ----------
    telescopename : str
        Name of telescope, as registered in TelescopesWiz() instance.
    stnid : str
        Name or ID of the station, as registered in TelescopesWiz() instance.
    band : str
        Name of band, as registered in TelescopesWiz() instance.
    antmodel : str
        Name of antenna model, e.g. 'Hamaker', as registered in TelescopesWiz()
        instance.
    obstimebeg : datetime.datetime
        Date-time when the tracking observation begins.
    obsdur : datetime.deltatime
        Duration of the entire tracking observation in seconds. The sample
        at obstimebeg+duration is included.
    obstimestp : datetime.deltatime
        Time step in seconds for which the jones matrix should be sampled at.
    pointingdir : (float, float, str)
        Length 3 tuple encoding the tracking direction on the celestial sphere.
        The last tuple element should usually be 'J2000', in which case the
        the first two tuple elements are the right ascension and declination,
        respectively, in radians.
    do_parallactic_rot : bool (optional)
        Whether of not to perform parallactic rotation (default True).

    Returns
    -------
    timespy : array
        The python datetime of the samples.
    freqs : array
        The frequency [Hz] at which the Jones was computed.
    jones : array_like
        Array over time steps and frequencies of the Jones matrix corresponding
        to RIME for this set of input parameters.
    jonesobj : Jones(object)
        Resulting instance of Jones(). jonesobj.jones is a copy of the
        other return variable: jones. In addition, it has more info regarding
        this Jones matrix such as its basis.

    Notes
    -----
    Specifically, this function only considers an RIME consisting of the
    projection from the celestial (Earth-Centered-Inertial) frame to the
    topocentric station frame, and the polarimetric responce of the
    dual-polarized antenna feed. These are known as the P-Jones and the E-Jones
    respectively.

    Examples
    --------
    Here is an example where Cassiopeia A was tracked with the Swedish
    LOFAR telescope HBA band starting at '2019-08-03T12:00:00' and lasting
    24 hours (86400s) sampled every hour (3600s):

    >>> from dreambeam.rime.scenarios import *
    >>> from datetime import datetime, deltatime
    >>> obstimebeg = datetime.strptime('2019-08-30T12:00:00',
    ... "%Y-%m-%dT%H:%M:%S")
    >>> duration = deltatime(hours=24)
    >>> obstimestp = deltatime(hours=1)
    >>> pointingdir = (6.11, 1.02, 'J2000')
    >>> samptimes, freqs, jones, jonesobj = on_pointing_axis_tracking('LOFAR',
    ... 'HBA', 'Hamaker', 'SE607', obstimebeg, duration, obstimestp,
    ... pointingdir)
    >>> print(jones.shape)
    (1024, 25, 2, 2)
    >>> print(freqs.shape)
    (1024,)
    >>> print(samptimes.shape)
    (25,)
    >>> print(jones[512,0,:,:])
    [[ 0.35038040-0.02403619j  0.46850486+0.01755369j]
     [ 0.39298880-0.02620409j -0.38686167-0.01691861j]]

    """
    # Startup a telescope wizard
    tw = TelescopesWiz()

    # Get the telescopeband instance:
    telescope = tw.getTelescopeBand(telescopename, band, antmodel)

    #    *Setup Source*
    srcfld = dreambeam.rime.jones.DualPolFieldPointSrc(pointingdir)

    stnBD = telescope['Station'][stnid]
    stnRot = stnBD.stnRot

    #    *Setup PJones*
    timespy = []
    nrTimSamps = int((obsdur.total_seconds()/obstimestp.seconds))+1
    for ti in range(0, nrTimSamps):
        timespy.append(obstimebeg+ti*obstimestp)
    pjones = dreambeam.rime.jones.PJones(timespy, np.transpose(stnRot),
                                         do_parallactic_rot=do_parallactic_rot)

    #    *Setup EJones*
    ejones = stnBD.getEJones(pointingdir)
    stnDPolel = stnBD.feed_pat
    freqs = stnDPolel.getfreqs()

    #    *Setup MEq*
    pjonesOfSrc = pjones.op(srcfld)
    jonesobj = ejones.op(pjonesOfSrc)

    # Get the resulting Jones matrices
    # (structure is Jn[freqIdx, timeIdx, chanIdx, compIdx] )
    jones = jonesobj.getValue()
    if not do_parallactic_rot:
        basis_from = jonesobj.get_basis()
        basis_to = pjonesOfSrc.get_basis()
        btransmat2d = basis2basis_transf(basis_from, basis_to)[..., 1:, 1:]
        # Tranformation from IAU2C09 has to be (right) transformed back 1st
        transmat2d = np.matmul(C09toIAU[1:, 1:], btransmat2d)
        jones = np.matmul(jones, transmat2d)
        jonesobj.jones = jones
    return timespy, freqs, jones, jonesobj


def primarybeampat(telescopename, stnid, band, antmodel, freq,
            pointing=(0., np.pi/2, 'STN'), obstime=None, lmgrid=None):
    """Computes the Jones matrix over the beam fov for pointing.
    """
    # Startup a telescope wizard
    tw = TelescopesWiz()

    # Get the telescopeband instance:
    telescope = tw.getTelescopeBand(telescopename, band, antmodel)
    stnBD = telescope['Station'][stnid]
    stnRot = stnBD.stnRot

    #    *Setup Source*
    (az, el, refframe) = pointing
    srcfld = dreambeam.rime.jones.DualPolFieldRegion(refframe, iaucmp=False,
                                                     lmgrid=lmgrid)

    #    *Setup Parallatic Jones*
    pjones = dreambeam.rime.jones.PJones([obstime], np.transpose(stnRot))

    #    *Setup EJones*
    stnDPolel = stnBD.feed_pat
    #    **Select frequency
    freqs = stnDPolel.getfreqs()
    frqIdx = np.where(np.isclose(freqs, freq, atol=190e3))[0][0]
    # N.B. Ejones doesn't really use pointing
    ejones = stnBD.getEJones(pointing, [freqs[frqIdx]])

    #    *Setup MEq*
    pjones_src = pjones.op(srcfld)
    if refframe == dreambeam.rime.jones.Jones._topo_frame:
        j2000basis = pjones_src.jonesbasis
        # Since refFrame is STN, pjones is inverse of ordinary J2000 to STN.
        # By inverting it, one gets the ordinary conversion back.
        pjones_src = dreambeam.rime.jones.inverse(pjones_src)
    else:
        j2000basis = srcfld.jonesbasis
    res = ejones.op(pjones_src)
    # Because we started off with iaucmp=False, but want IAU components:
    res.convert2iaucmp()

    dreambeam.rime.jones.fix_imaginary_directions(res)
    # NOTE: Not using get_basis() method, in order to get station basis instead
    # of antenna basis:
    stnbasis = res.jonesbasis
    # Get the resulting Jones matrices
    # (structure is Jn[freqIdx, timeIdx, chanIdx, compIdx] )
    res_jones = res.getValue()

    return res_jones, stnbasis, j2000basis
