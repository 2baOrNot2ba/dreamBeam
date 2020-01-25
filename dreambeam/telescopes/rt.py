"""
rt (i.e. Radio Telescopes) module is for handling real telescope meta-data.
"""
import os
import glob
import pickle
import dreambeam.telescopes as dbtel
import dreambeam.telescopes.geometry_ingest as gi

TELESCOPES_DIR = os.path.dirname(dbtel.__file__)
ANTMODELS = ['Hamaker']
SHAREDIR = 'share'  # Dir for native telescope project data.
DATADIR = 'data'    # Dir for telescope data for RIME level work.


class TelescopeBndStn(object):
    """Model of one station and one band of a telescope."""
    feed_pat = None

    def __init__(self, stnPos, stnRot):
        """Set the station's position and attitude."""
        self.stnPos = stnPos
        self.stnRot = stnRot

    def getEJones(self):
        """Create ejones for station based on antenna patterns."""
        ejones = None
        return ejones


def open_telescopebndmodel(tscopename, band, beammodel):
    """
    Open the TelescopeBndStn object given by tscopename, band and beammodel.
    """
    tbdata_dir = _get_telbnddatadir(tscopename)
    tbdata_fname = _get_teldat_fname(band, beammodel)
    tbdata_path = os.path.join(tbdata_dir, tbdata_fname)
    with open(tbdata_path, 'rb') as f:
        telbnddata = pickle.load(f)
    return telbnddata


class TelescopesWiz():
    """Database over available telescopes patterns."""

    def __init__(self):
        # Register telescope plugin
        # A telescope plugin must exist in the TELESCOPES_DIR and have a
        # directory named DATADIR.
        ls = os.listdir(TELESCOPES_DIR)
        ds = []
        for p in ls:
            if os.path.isdir(os.path.join(TELESCOPES_DIR, p)):
                ds.append(p)
        self.tbdata = {}
        for dd in ds:
            tbdata_dir = _get_telbnddatadir(dd)
            if os.path.isdir(tbdata_dir):
                t = os.path.basename(dd)
                self.tbdata[t] = {}
        # Find bands & models per telescope
        for tel in self.tbdata.keys():
            teldat_path = _get_telbnddatadir(tel)
            tbfiles = glob.glob(teldat_path+'/*_*.teldat.pkl')
            bands = []
            antmodels = []
            for tbfile in tbfiles:
                filename = os.path.basename(tbfile)
                (band, modelsuffix) = filename.split('_')
                bands.append(band)
                antmodel = modelsuffix.split('.', 2)[0]
                antmodels.append(antmodel)
            self.tbdata[tel] = {}
            for band in bands:
                self.tbdata[tel][band] = {}
                for antmodel in antmodels:
                    self.tbdata[tel][band][antmodel] = {}
        # Find stations for telescope band antmodel:
        for tel in self.tbdata.keys():
            for band in self.tbdata[tel].keys():
                for beammodel in self.tbdata[tel][band].keys():
                    telbnddata = open_telescopebndmodel(tel, band, beammodel)
                    self.tbdata[tel][band][antmodel] = \
                        telbnddata['Station'].keys()

    def get_telescopes(self):
        return self.tbdata.keys()

    def get_bands(self, telescope):
        return self.tbdata[telescope].keys()

    def get_stations(self, telescope, band):
        abeammodel = self.tbdata[telescope][band].keys()[0]
        return self.tbdata[telescope][band][abeammodel]

    def get_beammodels(self, telescope, band):
        return self.tbdata[telescope][band].keys()


def _get_teldat_fname(band, beammodel):
    """Get telescope data file name.
    """
    tbdata_fname = band+"_"+beammodel+".teldat.pkl"
    return tbdata_fname


def _get_telbnddatadir(tscopename):
    """Get path to data directory for tscopename.
    Format of the path returned is:
       absolute_directory=/<TELESCOPES_DIR>/<tscopename>/data/
    """
    tbdata_dir = os.path.join(TELESCOPES_DIR, tscopename, DATADIR)
    return tbdata_dir


def save_telescopeband(tscopename, band, DP_BAfile, telbndstn_class,
                       polcrdrot, antmodel='Hamaker'):
    """
    Save all the data relevant to the telescope-band beam modeling into
    one file.
    """
    print("""Generating '{}' beam-model for the band {} of the {} telescope
          with stations:""".format(antmodel, band, tscopename))
    # Create telescope-bands metadata:
    telescope = {'Name': tscopename, 'Band': band, 'Beam-model': antmodel}
    #   * Create station's antenna model
    if antmodel == 'Hamaker':
        # This is an example of dual-pol element built from a monolithic
        # Jones representation.
        dpepath = os.path.join(_get_telbnddatadir(tscopename), DP_BAfile)
        with open(dpepath, 'rb') as fp:
            stnDPolel = pickle.load(fp)
    # Rotate 45 degrees since LOFAR elements are 45 degrees to meridian:
    stnDPolel.rotateframe(polcrdrot)

    # Create telescope_band_station metadata:
    telescope['Station'] = {}
    x, y, z, diam, stnIds = gi.readarrcfg(tscopename, band)

    for stnId in stnIds:
        print(stnId)
        #    *Setup station Jones*
        # Get metadata for the LOFAR station. stnRot is transformation matrix
        #  ITRF_crds = stnRot*LOFAR_crds
        stnid_idx = stnIds.tolist().index(stnId)
        stnPos = [x[stnid_idx], y[stnid_idx], z[stnid_idx]]
        stnRot = gi.readalignment(tscopename, stnId, band)
        # Create a StationBand object for this
        stnbnd = telbndstn_class(stnPos, stnRot)
        stnbnd.feed_pat = stnDPolel
        telescope['Station'][stnId] = stnbnd
    teldatdir = _get_telbnddatadir(tscopename)
    savename = _get_teldat_fname(band, antmodel)
    with open(os.path.join(teldatdir, savename), 'wb') as fp:
        pickle.dump(telescope, fp, pickle.HIGHEST_PROTOCOL)
    print("Saved '"+savename+"' in "+teldatdir)
