"""Script to generate LOFAR antenna response data."""
import os
import numpy
from antpat.reps.hamaker import convLOFARcc2DPE
from dreambeam.telescopes import rt
from dreambeam.telescopes.LOFAR.feeds import LOFAR_LHBA_stn

TELESCOPE_NAME = 'LOFAR'
NR_POLS = 2
SAMPFREQ = 100e6
NR_CHANNELS = 512
OOSR2 = 1./numpy.sqrt(2)
# Rotation of LOFAR antennas from build frame to station frame.
# It takes x/y and directed dipoles and places them along (-1,-1)/(+1,-1) resp.
POLCRDROT = numpy.array([[-OOSR2, +OOSR2,  0.],
                         [-OOSR2, -OOSR2,  0.],
                         [    0.,     0.,  1.]])

CCARTS_LBA_DEF = 'DefaultCoeffLBA.cc'
CCARTS_HBA_DEF = 'DefaultCoeffHBA.cc'
DP_LBA_FILE_DEF = 'DP_model_LBA.pkl'
DP_HBA_FILE_DEF = 'DP_model_HBA.pkl'
DP_BAFILES = {'LBA': DP_LBA_FILE_DEF, 'HBA': DP_HBA_FILE_DEF}
BANDS = DP_BAFILES.keys()


def gen_antmodelfiles(inpfileL=CCARTS_LBA_DEF,
                      inpfileH=CCARTS_HBA_DEF,
                      outfileL=DP_LBA_FILE_DEF,
                      outfileH=DP_HBA_FILE_DEF):
    """Reads the 'lofar_elem_resp' packages c++ header files of Hamaker-Arts
    coefficients and writes pickled instances of DualPolElem class.
    Also adds nominal LOFAR frequency channels."""
    inpfileL = os.path.join(rt.SHAREDIR, inpfileL)
    inpfileH = os.path.join(rt.SHAREDIR, inpfileH)
    outfileL = os.path.join(rt.DATADIR, outfileL)
    outfileH = os.path.join(rt.DATADIR, outfileH)
    # Adding nominal frequency channels. The HA model for LOFAR has two bands
    # while data recording S/W has 3 intervals based on sampling frequency,
    # namely (0,100), (100,200), (200,300), each with 512 channels.
    # Here I concatenate the two latter intervals.
    channelsL = numpy.linspace(0., SAMPFREQ, NR_CHANNELS, endpoint=False)
    convLOFARcc2DPE(inpfileL, channelsL, outfileL)
    channelsH = numpy.linspace(SAMPFREQ, 3*SAMPFREQ, 2*NR_CHANNELS,
                               endpoint=False)
    convLOFARcc2DPE(inpfileH, channelsH, outfileH)


if __name__ == "__main__":
    """Use this to produce telescope data files for use in dreamBeam, or when
    configuration data has changed. Run this as:
    $ python telwizhelper.py
    """
    gen_antmodelfiles()
    antmodel = 'Hamaker'
    for band in BANDS:
        rt.save_telescopeband(TELESCOPE_NAME, band, DP_BAFILES[band],
                              LOFAR_LHBA_stn, POLCRDROT, antmodel)
    print("Completed setup.")
