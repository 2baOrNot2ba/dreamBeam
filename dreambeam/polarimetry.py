import numpy


def convertxy2stokes(cov_xx, cov_xy, cov_yx, cov_yy, iau_convention=True,
                     phase_sign_pos=True):
    r"""Convert hermitian form XX,XY,YX,YY covariance matrix components to
    Stokes parameters.

    A covariance matrix in standard basis takes the form:
        covmat = [[cov(x, x^*), cov(x, y^*)],
                  [cov(y, x^*), cov(y, y^*]]
    where cov(,) is the covariance operation and x^* is the conjugate of x.
    Note the conjugation order adopted here.

    The Stokes representation of such a matrix takes two possible forms:
        covmat = 1/2*[[I+Q,   U-/+i*V],
                      [U+/-i*V, I+Q]]
    where -/+ and +/- represents either - and + respectively, case 1, or + and
    -, case 2. Disregarding the factor 1/2, case 1 corresponds to the standard
    Pauli spin matrices, while case 2 is the standard Pauli matrices except
    that Pauli matrix sigma_3 has a its sign flipped.

    The sign in front of the Stokes V depends on the conventions used w.r.t.
    two different quantities. Firstly, the sign the of Stokes V w.r.t. IAU/IEEE
    right-handed circular polarized signals should be positive. (Note that
    radio polarimetry usual uses the opposite convention) Secondly, assuming a
    positive frequency, the sign of the temporal phase can either be increasing
    (positive sign) or decreasing (negative sign) in time. The parameter
    `iau_convention` set to True means that the output should conform to the
    IAU convention. The parameter `phase_sign_pos` set to True means that the
    quantities all have a temporal dependence proportional to exp(+i*omega*t),
    in the opposite case the dependency is exp(-i*omega*t).

    The returned Stokes parameters are thus determined by the following table:

    +---------------------------------+--------+--------+
    | iau_convention \ phase_sign_pos | True   | False  |
    +=================================+========+========+
    |   True                          | case 2 | case 1 |
    +---------------------------------+--------+--------+
    |   False                         | case 1 | case 2 |
    +---------------------------------+--------+--------+

    Parameters
    ----------
    cov_xx : array
        xy components covariance matrix, i.e. Cov(x, conj(x))
    cov_xy : array
        xy components covariance matrix, i.e. Cov(x, conj(y))
    cov_yx : array
        yx components covariance matrix, i.e. Cov(y, conj(x))
    cov_yy : array
        yy components covariance matrix, i.e. Cov(y, conj(xy))
    iau_convention : bool
        Compute Stokes according to IAU conventions? (Default True)
    phase_sign_pos : bool
        Is phase sign positive, i.e. does phase increase in time?
        (Default True)

    Returns
    -------
    stokes : tuple
        Index refers to Stokes vector components:
        (stokes[0], stokes[1], stokes[2], stokes[3]) = (I,Q,U,V) in accordance
        with the preference for the IAU convention or not, and the sign of
        phase.

    Notes
    -----
    This function returns a Stokes I that is a factor 2 larger than the average
    of the XX and YY components, as is used in radio software such as CASA and
    AIPS. For completeness, note that the covariance matrix is assumed to be in
    a linear, not circular basis. Finally, note the conjugation order in the
    covariance matrix is such that XY represents the expectation of X*Y^*.

    References
    ----------
    .. [1] van Straten et al., PASA, 27, 104, 2010.
    .. [2] Hamaker et al., A&A Sup., 117, 161, 1996.

    Examples
    --------
    >>> convertxy2stokes(1., 1.+1.j, 1.-1.j, 1.,
    ...                  iau_convention=True, phase_sign_pos=True)
    (2.0, 0.0, 2.0, 2.0)
    >>> convertxy2stokes(1., 1.+1.j, 1.-1.j, 1.,
    ...                  iau_convention=True, phase_sign_pos=False)
    (2.0, 0.0, 2.0, -2.0)

    """
    # Compute according to Pauli convention:
    si = +1*numpy.real(cov_xx + cov_yy)
    sq = +1*numpy.real(cov_xx - cov_yy)
    su = +1*numpy.real(cov_xy + cov_yx)
    sv = -1*numpy.imag(cov_xy - cov_yx)
    pauli_std = not(iau_convention == phase_sign_pos)
    if not pauli_std:
        # Standard Pauli is not to be used, so flip sign of Stokes V.
        sv *= -1
    return (si, sq, su, sv)
