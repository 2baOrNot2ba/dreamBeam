"""
rt (i.e. Radio Telescopes) module is for handling real telescope meta-data.
"""
import os
import pickle
import importlib
import pkgutil
from antpat.reps.hamaker import convLOFARcc2DPE
import dreambeam.telescopes as dbtel
import dreambeam.telescopes.geometry_ingest as gi
from dreambeam.telescopes.feeds import FixedMountStn


TELESCOPES_DIR = os.path.dirname(dbtel.__file__)
ANTMODELS = ['Hamaker']
SHAREDIR = 'share'  # Dir for native telescope project data.
DATADIR = 'data'    # Dir for telescope data for RIME level work.


def _get_helpers():
    telescope_plugins = {}
    # print(dbtel.__path__, dbtel.__name__ )
    for _, pluginname, ispkg in pkgutil.iter_modules(dbtel. __path__):
        # print("Found submodule %s (is a package: %s)" % (pluginname, ispkg))
        if ispkg:
            tpi = os.path.join(dbtel.__path__[0], pluginname)
            # print(tpi)
            for _, modnamesub, ispkg in pkgutil.iter_modules([tpi], ""):
                if modnamesub == "telwizhelper":
                    pluginpath = dbtel.__name__+"."+pluginname
                    telwizmod = importlib.import_module(".telwizhelper",
                                                        pluginpath)
                    telescope_plugins[pluginname] = telwizmod.telhelper
    return telescope_plugins


def open_telescopebndmodel(tscopename, band, beammodel):
    """
    Open the TelescopeBndStn object given by tscopename, band and beammodel.
    """
    telescope_plugins = _get_helpers()
    telbnddata = telescope_plugins[tscopename].open_bndmodel(band, beammodel)
    return telbnddata


class TelescopesWiz():
    """Database over available telescopes patterns."""

    def __init__(self):
        telescope_plugins = _get_helpers()
        self.tbdata = {}
        # Find bands & models per telescope
        for tel in telescope_plugins.keys():
            bands = telescope_plugins[tel].get_bands()
            antmodels = telescope_plugins[tel].get_beammodels()
            self.tbdata[tel] = {}
            for band in bands:
                self.tbdata[tel][band] = {}
                for antmodel in antmodels:
                    self.tbdata[tel][band][antmodel] = {}
        # Find stations for telescope band antmodel:
        for tel in self.tbdata.keys():
            for band in self.tbdata[tel].keys():
                for beammodel in self.tbdata[tel][band].keys():
                    bnddata = telescope_plugins[tel].open_bndmodel(band,
                                                                   beammodel)
                    self.tbdata[tel][band][antmodel] = \
                        bnddata['Station'].keys()

    def get_telescopes(self):
        return self.tbdata.keys()

    def get_bands(self, telescope):
        return self.tbdata[telescope].keys()

    def get_stations(self, telescope, band):
        abeammodel = self.tbdata[telescope][band].keys()[0]
        return self.tbdata[telescope][band][abeammodel]

    def get_beammodels(self, telescope, band):
        return self.tbdata[telescope][band].keys()


class TelWizHelper(object):
    """Plugin for a generic telescope."""

    def __init__(self, tscopename, bandchns, haversion, modeltype, polcrdrot):
        self.name = tscopename
        self.bandchns = bandchns
        self.haversion = haversion
        self.modeltype = modeltype
        self.polcrdrot = polcrdrot
        self.bands = self.bandchns.keys()

    def get_bands(self):
        return self.bands

    def get_beammodels(self):
        return [self.modeltype]

    def _get_teldat_fname(self, band, modeltype):
        """Get telescope data file name.
        """
        tbdata_fname = band+"_"+modeltype+".teldat.pkl"
        return tbdata_fname

    def _get_inpfile(self, band):
        inpfile = self.haversion+"Coeff"+band+".cc"
        return inpfile

    def _get_outfile(self, band):
        outfile = "DP_model_"+band+".pkl"
        return outfile

    def get_stndpolel(self, dp_bafile):
        """Get DualPolElem object for this station."""
        dpepath = os.path.join(self.path_, "data", dp_bafile)
        with open(dpepath, 'rb') as fp:
            stnDPolel = pickle.load(fp)
        return stnDPolel

    def open_bndmodel(self, band, beammodel):
        """
        Open the TelescopeBndStn object given by band and beammodel.
        """
        tbdata_fname = self._get_teldat_fname(band, self.modeltype)
        tbdata_path = os.path.join(self.path_, "data", tbdata_fname)
        with open(tbdata_path, 'rb') as f:
            telbnddata = pickle.load(f)
        return telbnddata

    def gen_antmodelfiles(self, band):
        """Reads the 'lofar_elem_resp' packages c++ header files of
        Hamaker-Arts coefficients and writes pickled instances of DualPolElem
        class.
        Also adds nominal LOFAR frequency channels."""
        inpfile = self._get_inpfile(band)
        outfile = self._get_outfile(band)
        inppath = os.path.join(self.path_, SHAREDIR, inpfile)
        outpath = os.path.join(self.path_, DATADIR, outfile)
        channels = self.bandchns[band]
        convLOFARcc2DPE(inppath, channels, outpath)

    def telescope_specific_init(self):
        """Override this method for telescope specific initialization."""
        pass

    def initialize(self):
        """Use this to produce telescope data files for use in dreamBeam, or
        when configuration data has changed."""
        self.telescope_specific_init()
        for band in self.bands:
            outfile = self._get_outfile(band)
            # channels = self.bandchns[band]
            self.gen_antmodelfiles(band)
            self.save_telescopeband(band, outfile, FixedMountStn,
                                    self.polcrdrot, self.modeltype)

    def save_telescopeband(self, band, DP_BAfile, telbndstn_class,
                           polcrdrot, antmodel='Hamaker'):
        """
        Save all the data relevant to the telescope-band beam modeling into
        one file.
        """
        tscopename = self.name
        print("""Generating '{}' beam-model for the band {} of the {} telescope
              with stations:""".format(antmodel, band, tscopename))
        # Create telescope-bands metadata:
        telescope = {'Name': tscopename, 'Band': band, 'Beam-model': antmodel}
        #   * Create station's antenna model
        if antmodel == 'Hamaker':
            stnDPolel = self.get_stndpolel(DP_BAfile)
        # Rotate 45 degrees since LOFAR elements are 45 degrees to meridian:
        stnDPolel.rotateframe(polcrdrot)

        # Create telescope_band_station metadata:
        telescope['Station'] = {}
        x, y, z, diam, stnIds = gi.readarrcfg(tscopename, band)

        for stnId in stnIds:
            print(stnId)
            #    *Setup station Jones*
            # Get metadata for the LOFAR station. stnRot is transformation
            # matrix:
            #   ITRF_crds = stnRot*LOFAR_crds
            stnid_idx = stnIds.tolist().index(stnId)
            stnPos = [x[stnid_idx], y[stnid_idx], z[stnid_idx]]
            stnRot = gi.readalignment(tscopename, stnId, band)
            # Create a StationBand object for this
            stnbnd = telbndstn_class(stnPos, stnRot)
            stnbnd.feed_pat = stnDPolel
            telescope['Station'][stnId] = stnbnd
        teldatdir = self.path_
        savename = self._get_teldat_fname(band, antmodel)
        with open(os.path.join(teldatdir, DATADIR, savename), 'wb') as fp:
            pickle.dump(telescope, fp, pickle.HIGHEST_PROTOCOL)
        print("Saved '"+savename+"' in "+teldatdir)
