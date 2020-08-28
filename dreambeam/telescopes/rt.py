"""
rt (i.e. Radio Telescopes) module is for handling real telescope meta-data.
"""
import os
import pickle
import importlib
import pkg_resources
import dreambeam
from dreambeam.feeds.feedplugins import FeedWiz
import dreambeam.telescopes.geometry_ingest as gi


def get_tel_plugins():
    resource_path = '/'.join(('configs', 'telescope_paths.txt'))
    rootpath = os.path.dirname(os.path.dirname(dreambeam.__file__))
    tele_paths_file = pkg_resources.resource_filename(dreambeam.__name__,
                                                      resource_path)
    tele_paths = []
    with open(tele_paths_file) as fp:
        lines = fp.readlines()
        for line in lines:
            line = line.rstrip()
            if not (line.startswith('#') or line == ''):
                tele_path = line
                tele_paths.append(tele_path)
    teleplugins = {}
    for tele_path in tele_paths:
        telemodpath = os.path.join(rootpath, tele_path, '_telescope.py')
        if os.path.exists(telemodpath):
            pluginpath = tele_path.replace("/", ".")
            telwizmod = importlib.import_module("._telescope",
                                                pluginpath)
            telwizmod.telhelper.path_ = os.path.join(rootpath, tele_path)
            telescopename = os.path.basename(tele_path)
            teleplugins[telescopename] = telwizmod.telhelper
    return teleplugins


def load_mountedfeed(tscopename, station, band, modelname):
    """
    Open the TelescopeBndStn object given by tscopename, band and model.
    """
    telescope_plugins = get_tel_plugins()
    stnfeed = telescope_plugins[tscopename].getstationfeed(station, band,
                                                           modelname)
    return stnfeed


class TelescopePlugin(object):
    """Plugin for a generic telescope."""
    DATADIR = 'data'    # Dir for telescope data for RIME level work.

    def __init__(self, tscopename, mountedfeed_class, feedrot):
        self.name = tscopename
        self.mountedfeed_class = mountedfeed_class
        self.feedrot = feedrot
        self.feedplugin = FeedWiz()[tscopename]
        self.bands = self.feedplugin.get_bands()
        self.stations = {}
        self.positions = {}
        self.diams = {}
        for band in self.bands:
            xs, ys, zs, diams, stnids = gi.readarrcfg(self.name, band)
            self.stations[band] = stnids.tolist()
            self.positions[band] = list(zip(xs.tolist(), ys.tolist(), zs.tolist()))
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

    def get_beammodels(self, band):
        return self.feedplugin.list_models4band(band)

    def _get_teldat_fname(self, band):
        """Get telescope data file name.
        """
        tbdata_fname = band+".teldat.pkl"
        return tbdata_fname

    def getstationfeed(self, station, band, modelname):
        #   * Create station's antenna model
        stndpolel = self.feedplugin.load_dpolel(band, modelname)
        # Rotate by polcrdrot:
        stndpolel.rotateframe(self.feedrot)
        stnpos = self.get_bandpositions(band)[station]
        stnrot = self.get_bandstnrot()[band][station]
        stationfeed = self.mountedfeed_class(stnpos, stnrot)
        stationfeed.mountfeed(stndpolel)
        return stationfeed
