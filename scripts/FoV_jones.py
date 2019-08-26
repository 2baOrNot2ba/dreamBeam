#!/usr/bin/env python
"""Show LOFAR element beam pattern.
"""
import sys
from datetime import datetime
from dreambeam.rime.scenarios import beamfov
from dreambeam.telescopes.rt import TelescopesWiz
from dreambeam.rime.jones import plotJonesField


def printJonesField(az, el, Jnf):
    # Select one frequency
    for idxi in range(az.shape[0]):
        for idxj in range(az.shape[1]):
            print("az, el:", az[idxi, idxj], el[idxi, idxj])
            print("Jones:", Jnf[idxi, idxj, 0, 0], Jnf[idxi, idxj, 0, 1],
                  Jnf[idxi, idxj, 1, 0], Jnf[idxi, idxj, 1, 1])


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

# Example:
# $ pointing_jones.py print LOFAR LBA SE607 Hamaker 2012-04-01T01:02:03 60 1\
#   6.11 1.02 60E6

if __name__ == "__main__":
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
        celdir = (float(args[1]), float(args[2]), 'J2000')
    except IndexError:
        print("Specify pointing direction (in radians): RA DEC")
        print(USAGE)
        exit()
    try:
        freq = float(args[3])
    except ValueError:
        print("Specify frequency (in Hz).")
        print(USAGE)
        exit()

    # Compute the Jones matrices
    az, el, jonesfld, jbasis = beamfov(telescopename, band, antmodel, stnid,
                                       obstime, celdir, freq)
    # Do something with resulting Jones according to cmdline args
    if action == "plot":
        plotJonesField(az, el, jonesfld, jbasis, rep='Stokes')
    else:
        printJonesField(az, el, jonesfld)
