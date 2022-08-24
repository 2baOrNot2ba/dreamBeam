"""Module that creates plugin for the NenuFAR telescope"""
from os.path import abspath, basename, dirname
import numpy
from dreambeam.telescopes.mounts import MountedFeedFixed
from dreambeam.telescopes.rt import TelescopePlugin

OOSR2 = 1./numpy.sqrt(2)
# Rotation of NenuFAR antennas from build frame to station frame.
# It takes x/y and directed dipoles and places them along (+1,-1)/(-1,-1) resp.
# Note that this is flipped w.r.t. LOFAR alignment i.e. X/Y <-> Y/X.
# Reference: Communication with NenuFAR people; unfortunately not documented. 
POLCRDROT = numpy.array([[+OOSR2, -OOSR2,  0.],
                         [-OOSR2, -OOSR2,  0.],
                         [    0.,     0.,  1.]])
TELESCOPE_NAME = basename(dirname(abspath(__file__)))

telhelper = TelescopePlugin(TELESCOPE_NAME, MountedFeedFixed, POLCRDROT)
