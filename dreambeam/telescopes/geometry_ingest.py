"""This package retrieves the metadata for a telescope, such as position,
alignment and bands.

The metadata for a telescope should be in the folder '<telescope>/share/'.
Subfolders should be 'simmos' (for array configurations of stations) and
'alignment' for station rotation.
"""
import os
import numpy as np

CASA_CFG_DTYPE = [('X', float), ('Y', float), ('Z', float), ('Diam', float),
                  ('Name', 'S7')]
CASA_CFG_FMT = '%12f %12f %12f %4.1f %s'
SIMOBS_PATH = os.path.dirname(os.path.realpath(__file__))


def readarrcfg(telescope, band):
    """Read a CASA array configuration file. This gives the positions of all
    the stations that make up the telescope.
    """
    arrcfg_folder = os.path.join(SIMOBS_PATH, telescope, 'share/simmos')
    arrcfg_file = os.path.join(arrcfg_folder, '{}_{}.cfg'.format(telescope,
                                                                 band))
    x, y, z, diam, name = np.loadtxt(arrcfg_file, CASA_CFG_DTYPE, unpack=True)
    return x, y, z, diam, name


def readalignment(telescope, stnid, band):
    """Read a alignment matrix file for a telescope station band.
    The alignment is 3x3 matrix that represents the rotation between the
    feed coord sys and the ITRF.
    """
    stnrot_folder = os.path.join(SIMOBS_PATH, telescope, 'share/alignment')
    stnrot_file = os.path.join(stnrot_folder, '{}_{}.txt'.format(stnid, band))
    stnrot = np.loadtxt(stnrot_file)
    return stnrot


def getArrayBandParams(telescope, stnid, band):
    """Get array and band parameters for a telescope, station and band."""
    x, y, z, diam, names = readarrcfg(telescope, band)
    stnid_idx = names.tolist().index(stnid)
    stnPos = [x[stnid_idx], y[stnid_idx], z[stnid_idx]]
    stnRot = readalignment(telescope, stnid, band)
    return stnPos, stnRot


def list_stations(telescope, band=None):
    """List all available stations."""
    x, y, z, diam, names = readarrcfg(telescope, band)
    return names


if __name__ == "__main__":
    print list_stations('LOFAR')
    stnPos, stnRot = getArrayBandParams('LOFAR', 'SE607', 'HBA')
    print stnPos
    print stnRot
