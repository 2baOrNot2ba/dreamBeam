"""Script to generate LOFAR antenna response data."""
import os
import numpy
from dreambeam.telescopes.rt import TelWizHelper


TELESCOPE_NAME = 'LOFAR'
MODELTYPE = 'Hamaker'

SAMPFREQ = 100e6
NR_CHANNELS = 512
OOSR2 = 1./numpy.sqrt(2)
# Rotation of LOFAR antennas from build frame to station frame.
# It takes x/y and directed dipoles and places them along (-1,-1)/(+1,-1) resp.
POLCRDROT = numpy.array([[-OOSR2, +OOSR2,  0.],
                         [-OOSR2, -OOSR2,  0.],
                         [    0.,     0.,  1.]])

HAVERSION = 'Default'
# Adding nominal frequency channels. The HA model for LOFAR has two
# bands while data recording S/W has 3 intervals based on sampling
# frequency, namely (0,100), (100,200), (200,300), each with
# 512 channels.
# Here I concatenate the two latter intervals.
BANDCHNS = {'LBA': numpy.linspace(0., SAMPFREQ, NR_CHANNELS, endpoint=False),
            'HBA': numpy.linspace(SAMPFREQ, 3*SAMPFREQ, 2*NR_CHANNELS,
                                  endpoint=False)}


class TelWizHelper_LOFAR(TelWizHelper):
    """Plugin for LOFAR telescope."""
    path_ = os.path.dirname(os.path.abspath(__file__))


telhelper = TelWizHelper_LOFAR(TELESCOPE_NAME, BANDCHNS, HAVERSION,
                               MODELTYPE, POLCRDROT)


if __name__ == "__main__":
    """Run this as:
    $ python telwizhelper.py
    """
    telhelper.initialize()
    print("Completed setup.")
