#!/usr/bin/env python
"""Show LOFAR element beam pattern.
"""
import sys
from datetime import datetime
import numpy as np
from dreambeam.rime.scenarios import primarybeampat
from dreambeam.telescopes.rt import TelescopesWiz
from dreambeam.rime.jones import plotJonesField


def printJonesField(jnf, jbasis):
    (nr_xs, nr_ys, _, _) = jbasis.shape
    print("x, y, J00, J01, J10, J11")
    for idxi in range(nr_xs):
        for idxj in range(nr_ys):
            x = np.real(jbasis[idxi, idxj, 0, 0])
            y = np.real(jbasis[idxi, idxj, 1, 0])
            J00 = jnf[idxi, idxj, 0, 0]
            J01 = jnf[idxi, idxj, 0, 1]
            J10 = jnf[idxi, idxj, 1, 0]
            J11 = jnf[idxi, idxj, 1, 1]
            jones_1f_outstring = ",".join(map(str, [x, y, J00, J01, J10, J11]))
            print(jones_1f_outstring)


def getnextcmdarg(args, mes):
    try:
        arg = args.pop(0)
    except IndexError:
        print("Specify "+mes)
        print(USAGE)
        exit()
    return arg


SCRIPTNAME = sys.argv[0].split('/')[-1]
USAGE = """Usage:{}
         print|plot telescope band stnID beammodel timeUTC
         pointingRA pointingDEC frequency""".format(SCRIPTNAME)


def main():
    """
    Plot or print the Jones matrices over the field-of-view.

    Example
    -------
    >>> FoV_jones print LOFAR LBA SE607 Hamaker 2012-04-01T01:02:03 \
           6.11 1.02 60E6
    """
    # Startup a telescope wizard
    tw = TelescopesWiz()
    # Process cmd line arguments
    args = sys.argv[1:]
    action = getnextcmdarg(args, "output-type:\n  'print' or 'plot'")
    telescopename = getnextcmdarg(args, "telescope:\n  "
                                  + ', '.join(tw.get_telescopes()))
    band = getnextcmdarg(args, "band/feed:\n  "
                         + ', '.join(tw.get_bands(telescopename)))
    stnid = getnextcmdarg(args, "station-ID:\n  "
                          + ', '.join(tw.get_stations(telescopename, band)))
    antmodel = getnextcmdarg(args, "beam-model:\n  "
                             + ', '.join(tw.get_beammodels(telescopename,
                                                           band)))
    del tw
    try:
        obstime = datetime.strptime(args[0], "%Y-%m-%dT%H:%M:%S")
    except IndexError:
        print("Specify time (UTC in ISO format: yy-mm-ddTHH:MM:SS ).")
        print(USAGE)
        exit()
    try:
        (az, el, refframe) = args[1].split(',')
        az, el = float(az), float(el)
    except ValueError:
        print("""Specify pointing direction (in radians): 'RA,DEC,J200'
                                         or 'AZ,EL,AZEL'""")
        print(USAGE)
        exit()
    try:
        freq = float(args[2])
    except ValueError:
        print("Specify frequency (in Hz).")
        print(USAGE)
        exit()

    if refframe == 'AZEL':
        refframe = 'STN'
        l = np.linspace(-1., 1., 10)
        m = np.linspace(-1., 1., 10)
        ll, mm = np.meshgrid(l, m)
        lmgrid = (ll, mm)
    else:
        lmgrid = None
    pointing = (az, el, refframe)

    # Compute the Jones matrices
    jonesfld, stnbasis, j2000basis = primarybeampat(
                                    telescopename, stnid, band, antmodel, freq,
                                    pointing=pointing, obstime=obstime,
                                    lmgrid=lmgrid)
    if refframe == 'STN':
        jbasis = stnbasis
    else:
        jbasis = j2000basis
    # Do something with resulting Jones according to cmdline args
    if action == "plot":
        plotJonesField(jonesfld, jbasis, refframe, rep='Stokes')
    else:
        printJonesField(jonesfld, jbasis)


if __name__ == "__main__":
    main()
