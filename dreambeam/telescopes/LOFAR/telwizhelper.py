"""Script to generate LOFAR antenna response data."""
import numpy
import pickle
from antpat.reps.hamaker import convLOFARcc2HA, convHA2DPE
from dreambeam.telescopes import rt
from dreambeam.telescopes.geometry_ingest import readarrcfg, readalignment
from dreambeam.telescopes.LOFAR.feeds import LOFAR_LHBA_stn

TELESCOPE_NAME = 'LOFAR'
NR_POLS = 2
SAMPFREQ = 100e6
NR_CHANNELS = 512
ANTMODELS = ['Hamaker']
OOSR2 = 1./numpy.sqrt(2)
PICKLE_PROTO = pickle.HIGHEST_PROTOCOL
# Rotation of LOFAR antennas from build frame to station frame.
# It takes x/y and directed dipoles and places them along (-1,-1)/(+1,-1) resp.
POLCRDROT = numpy.array([[-OOSR2, +OOSR2,  0.],
                         [-OOSR2, -OOSR2,  0.],
                         [    0.,     0.,  1.]])
LOFAR_HA_DATADIR = './share/'  # Dir for native telescope project data.
TELEDATADIR = 'data/'          # Dir for telescope data for RIME level work.
HA_LBA_FILE_DEF = TELEDATADIR+'HA_LOFAR_elresp_LBA.p'
HA_HBA_FILE_DEF = TELEDATADIR+'HA_LOFAR_elresp_HBA.p'
DP_LBA_FILE_DEF = TELEDATADIR+'DP_model_LBA.p'
DP_HBA_FILE_DEF = TELEDATADIR+'DP_model_HBA.p'
DP_BAFILES = {'LBA': DP_LBA_FILE_DEF, 'HBA': DP_HBA_FILE_DEF}

BANDS = DP_BAFILES.keys()

# Start up a telescope wizard:
TW = rt.TelescopesWiz()


def gen_antmodelfiles(inpfileL=LOFAR_HA_DATADIR+'DefaultCoeffLBA.cc',
                      inpfileH=LOFAR_HA_DATADIR+'DefaultCoeffHBA.cc',
                      outfileL=HA_LBA_FILE_DEF,
                      outfileH=HA_HBA_FILE_DEF):
    """A convenience function to produce the pickled 'artsdata' for the default
    LOFAR model data stored in the 'lofar_elem_resp' packages c++ header files.
    Also adds nominal LOFAR frequency channels."""

    # Adding nominal frequency channels. The HA model for LOFAR has two bands
    # while data recording S/W has 3 intervals based on sampling frequency,
    # namely (0,100), (100,200), (200,300), each with 512 channels.
    # Here I concatenate the two latter intervals.

    channels = numpy.linspace(0., SAMPFREQ, NR_CHANNELS, endpoint=False)
    convLOFARcc2HA(inpfileL, outfileL, channels)
    inpfileL = outfileL
    convHA2DPE(inpfileL, DP_LBA_FILE_DEF)
    channels = numpy.linspace(SAMPFREQ, 3*SAMPFREQ, 2*NR_CHANNELS,
                              endpoint=False)
    convLOFARcc2HA(inpfileH, outfileH, channels)
    inpfileH = outfileH
    convHA2DPE(inpfileH, DP_HBA_FILE_DEF)


def save_telescopeband(band, antmodel='Hamaker'):
    """Save all the data relevant to the telescope-band beam modeling into
    one file."""
    assert band in BANDS, ("Error: {} is not one of the available bands.\n"
                           "(Available bands are: {})").format(band, BANDS)
    assert antmodel in ANTMODELS, (
        """Error: {} is not one of the available models.
        Available models are: {}""").format(antmodel, ANTMODELS)
    print("""Generating '{}' beam-model for the band {} of the {} telescope
          with stations:""".format(antmodel, band, TELESCOPE_NAME))
    # Create telescope-bands metadata:
    telescope = {'Name': TELESCOPE_NAME, 'Band': band, 'Beam-model': antmodel}
    #   * Create station's antenna model
    DP_BAfile = DP_BAFILES[band]
    if antmodel == 'Hamaker':
        # This is an example of dual-pol element built from a monolithic
        # Jones representation.
        stnDPolel = pickle.load(open(DP_BAfile, 'rb'))
    # Rotate 45 degrees since LOFAR elements are 45 degrees to meridian:
    stnDPolel.rotateframe(POLCRDROT)

    # Create telescope_band_station metadata:
    telescope['Station'] = {}
    LOFAR_BA_stn = LOFAR_LHBA_stn
    x, y, z, diam, stnIds = readarrcfg(TELESCOPE_NAME, band)

    for stnId in stnIds:
        print(stnId)
        #    *Setup station Jones*
        # Get metadata for the LOFAR station. stnRot is transformation matrix
        #  ITRF_crds = stnRot*LOFAR_crds
        stnid_idx = stnIds.tolist().index(stnId)
        stnPos = [x[stnid_idx], y[stnid_idx], z[stnid_idx]]
        stnRot = readalignment(TELESCOPE_NAME, stnId, band)
        # Create a StationBand object for this
        stnbnd = LOFAR_BA_stn(stnPos, stnRot)
        stnbnd.feed_pat = stnDPolel
        telescope['Station'][stnId] = stnbnd
    teldatdir, saveName = TW.telbndmdl2dirfile(TELESCOPE_NAME, band, antmodel)
    pickle.dump(telescope, open('/'.join((teldatdir, saveName)), 'wb'),
                PICKLE_PROTO)
    print("Saved '"+saveName+"' in "+teldatdir)


if __name__ == "__main__":
    """Use this to produce telescope data files for use in dreamBeam, or when
    configuration data has changed. Run this as:
    $ python telwizhelper.py
    """
    gen_antmodelfiles()
    antmodel = 'Hamaker'
    for band in BANDS:
        save_telescopeband(band, antmodel)
    print("Completed setup.")
