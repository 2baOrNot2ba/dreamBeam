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


def get_tel_plugins():
    telescope_plugins = {}
    # print(dbtel.__path__, dbtel.__name__ )
    for _, pluginname, ispkg in pkgutil.iter_modules(dbtel. __path__):
        # print("Found submodule %s (is a package: %s)" % (pluginname, ispkg))
        if ispkg:
            plgospath = os.path.join(dbtel.__path__[0], pluginname)
            # print(tpi)
            for _, modnamesub, ispkg in pkgutil.iter_modules([plgospath], ""):
                if modnamesub == "telwizhelper":
                    pluginpath = dbtel.__name__+"."+pluginname
                    telwizmod = importlib.import_module(".telwizhelper",
                                                        pluginpath)
                    telwizmod.telhelper.path_ = plgospath
                    telescope_plugins[pluginname] = telwizmod.telhelper

    return telescope_plugins


def open_telescopebndmodel(tscopename, band, model):
    """
    Open the TelescopeBndStn object given by tscopename, band and model.
    """
    telescope_plugins = get_tel_plugins()
    telbnddata = telescope_plugins[tscopename].load_teldat(band, model)
    return telbnddata


class TelWizHelper(object):
    """Plugin for a generic telescope."""
    SHAREDIR = 'share'  # Dir for native telescope project data.
    DATADIR = 'data'    # Dir for telescope data for RIME level work.

    def __init__(self, tscopename, bandchns, haversion, modeltype, polcrdrot):
        self.name = tscopename
        self.bandchns = bandchns
        self.haversion = haversion
        self.modeltype = modeltype
        self.polcrdrot = polcrdrot
        self.bands = self.bandchns.keys()
        self.stations = {}
        self.positions = {}
        self.diams = {}
        for band in self.bands:
            xs, ys, zs, diams, stnids = gi.readarrcfg(self.name, band)
            self.stations[band] = stnids.tolist()
            self.positions[band] = zip(xs.tolist(), ys.tolist(), zs.tolist())
            self.diams[band] = diams.tolist()
        self.bandstnrot = {}
        for band in self.bands:
            self.bandstnrot[band] = {}
            for station in self.stations[band]:
                self.bandstnrot[band][station] = gi.readalignment(
                                                    self.name, station, band)

    def get_stations(self, band):
        return self.stations[band]

    def get_bandpositions(self, band):
        bandpositions = {}
        for stnidx, station in enumerate(self.stations[band]):
            bandpositions[station] = self.positions[band][stnidx]
        return bandpositions

    def get_diam(self, station, band):
        stnidx = self.stations[band].index(station)
        diam = self.diams[band][stnidx]
        return diam

    def get_bandstnrot(self):
        """Return a dict of rotation matrices of bands on stations.
        The dict has two keys: stnrot[<band>][<stnid>]. Value is the
        transformation matrix:
           ITRF_crds = stnrot*LOCAL_crds
        """
        return self.bandstnrot

    def get_bands(self):
        return self.bands

    def get_beammodels(self):
        return [self.modeltype]

    def _get_teldat_fname(self, band):
        """Get telescope data file name.
        """
        tbdata_fname = band+"_"+self.modeltype+".teldat.pkl"
        return tbdata_fname

    def _get_ccfile(self, band):
        inpfile = self.haversion+"Coeff"+band+".cc"
        return inpfile

    def _get_dpfile(self, band):
        outfile = "DP_model_"+band+".pkl"
        return outfile

    def load_stndpolel(self, band):
        """Get DualPolElem object for this station."""
        dp_bafile = self._get_dpfile(band)
        dpepath = os.path.join(self.path_, self.DATADIR, dp_bafile)
        try:
            with open(dpepath, 'rb') as fp:
                stndpolel = pickle.load(fp)
        except IOError:
            stndpolel = self.save_stndpolel(band)
        return stndpolel

    def save_stndpolel(self, band):
        """Generate dualPolElem for the given band.
        Reads the 'lofar_elem_resp' packages c++ header files of
        Hamaker-Arts coefficients and writes pickled instances of DualPolElem
        class.
        Also adds nominal LOFAR frequency channels."""
        inpfile = self._get_ccfile(band)
        outfile = self._get_dpfile(band)
        inppath = os.path.join(self.path_, self.SHAREDIR, inpfile)
        outpath = os.path.join(self.path_, self.DATADIR, outfile)
        channels = self.bandchns[band]
        try:
            stndpolel = convLOFARcc2DPE(inppath, channels, outpath)
        except IOError:
            self.telescope_specific_init()
            stndpolel = convLOFARcc2DPE(inppath, channels, outpath)
        return stndpolel

    def load_teldat(self, band, model):
        """
        Load the teldat dict given by band and beammodel for this telescope.
        """
        tbdata_fname = self._get_teldat_fname(band)
        tbdata_path = os.path.join(self.path_, self.DATADIR, tbdata_fname)
        try:
            with open(tbdata_path, 'rb') as f:
                telbnddata = pickle.load(f)
        except IOError:
            telbnddata = self.save_teldat(band, FixedMountStn)  ## FIXME:
        return telbnddata

    def save_teldat(self, band, telbndstn_class):
        """
        Save all the data relevant to the telescope-band beam modeling into
        one file.
        teldat is a dict with structure:
            Name: <tscopename>
            Band: <bandname>
            Beam-model: <modeltype>
            Station:
                <stnId0>: <stnbnd0>
                <stnId1>: <stnbnd1>
                ...
        """
        tscopename = self.name
        print("""Generating '{}' beam-model for the band {} of the {} telescope
              with stations:""".format(self.modeltype, band, tscopename))
        # Create telescope-bands metadata:
        telescope = {'Name': tscopename, 'Band': band,
                     'Beam-model': self.modeltype}
        #   * Create station's antenna model
        if self.modeltype == 'Hamaker':
            stnDPolel = self.load_stndpolel(band)
        # Rotate by polcrdrot:
        stnDPolel.rotateframe(self.polcrdrot)

        # Create telescope_band_station metadata:
        telescope['Station'] = {}
        for station in self.stations[band]:
            print station
            stnpos = self.get_bandpositions(band)[station]
            stnrot = self.get_bandstnrot()[band][station]
            telbndstn = telbndstn_class(stnpos, stnrot)
            telbndstn.feed_pat = stnDPolel
            telescope['Station'][station] = telbndstn
        teldatdir = self.path_
        savename = self._get_teldat_fname(band)
        with open(os.path.join(teldatdir, self.DATADIR, savename), 'wb') as fp:
            pickle.dump(telescope, fp, pickle.HIGHEST_PROTOCOL)
        print("Saved '"+savename+"' in "+teldatdir)
        return telescope

    def telescope_specific_init(self):
        """Override this method for telescope specific initialization."""
        pass

    def initialize(self):
        """Use this to produce telescope data files for use in dreamBeam, or
        when configuration data has changed."""
        self.telescope_specific_init()
        for band in self.bands:
            self.save_stndpolel(band)
            self.save_teldat(band, FixedMountStn)
