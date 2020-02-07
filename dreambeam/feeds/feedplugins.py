import os
import shutil
import pkg_resources
import pickle
import dreambeam
from antpat import radfarfield
from antpat.io.NECread import readNECout_tvecfuns
from antpat.reps.hamaker import convDPE2LOFARcc, convLOFARcc2DPE


class FeedWiz(object):

    def __init__(self):
        feed_paths = self._get_feed_paths()
        self.feedplugins = {}
        for feed_path in feed_paths:
            feedmodpath = os.path.join(feed_path, '_feeds.py')
            if os.path.exists(feedmodpath):
                feed_category = os.path.basename(feed_path)
                self.feedplugins[feed_category] = feed_path

    def _get_feed_paths(self):
        resource_path = '/'.join(('configs', 'feed_paths.txt'))
        rootpath = os.path.dirname(os.path.dirname(dreambeam.__file__))
        feed_paths_file = pkg_resources.resource_filename(dreambeam.__name__,
                                                          resource_path)
        feed_paths = []
        with open(feed_paths_file) as fp:
            lines = fp.readlines()
            for line in lines:
                line = line.rstrip()
                if not (line.startswith('#') or line == ''):
                    feed_path = os.path.join(rootpath, line)
                    feed_paths.append(feed_path)
        return feed_paths

    def __getitem__(self, feed_category):
        try:
            feed_path = self.feedplugins[feed_category]
        except KeyError:
            return None
        feedplugin = FeedPlugin(feed_path)
        return feedplugin

    def list_names(self):
        names = []
        for feed_plugin in self.feed_plugins:
            names.extend(feed_plugin.get_bands())
        return names

    def list_models4receptor(self, feed_name):
        models = []
        for feed_plugin in self.feed_plugins:
            models.extend(feed_plugin.list_models4receptor(feed_name))
        return models


class FeedPlugin(object):
    dpesuffix = '.dpe.pkl'
    _bndversep = '-'  # Band Version Separator
    _mdlversep = '-'  # Model Version Separator
    SHAREDIR = 'share'  # Dir for native telescope project data.
    DATADIR = 'data'    # Dir for telescope data for RIME level work.

    def __init__(self, path2plugin):
        self.path_ = path2plugin

    def prep_dpefiles(self):
        coeffilemetas = self._search_coeff_files()
        if coeffilemetas != []:
            for coeffilemeta in coeffilemetas:
                self.get_stndpolel(coeffilemeta['band'],
                                   coeffilemeta['modeltype'],
                                   coeffilemeta['version'])

    def find_necoutfiles(self):
        necfiles = []
        necoutdir = os.path.join(self.path_, self.SHAREDIR)
        ls = os.listdir(necoutdir)
        for lsentry in ls:
            if lsentry.endswith('.out'):
                necfiles.append(lsentry)
        return necfiles

    def _parsefilename(self, filename):
        filebase, fileext = filename.split('.', 1)
        band, version, = filebase.split(self._bndversep, 1)
        return band, version

    def _convNEC2antpat(self, necoutfile, scalefac):
        print("Generating antpat from necfile: {}".format(necoutfile))
        necoutfile = os.path.join(self.path_, self.SHAREDIR, necoutfile)
        tvf = readNECout_tvecfuns(necoutfile)
        tvf.scale(scalefac)
        antpat = radfarfield.RadFarField(tvf)
        return antpat

    def compute_ccfile(self, necoutfile, scalefac, freq_center, freq_range,
                       kord, tord, ford, channels):
        """Create a HamakerArts coefficients file from a NEC out file
        and save it in a .cc file.
        """
        antpat = self._convNEC2antpat(necoutfile, scalefac)
        band, haversion = self._parsefilename(necoutfile)
        # freq_center should be close to maximum gain frequency,
        # while freq_range should avoid going outside the frequencies
        # of NEC4 file which is 0 to 100e6 Hz (also avoid edges).
        print("Generating coefficients files...")
        artsdata, filename = convDPE2LOFARcc(antpat, freq_center, freq_range,
                                             HAcoefband=band,
                                             HAcoefversion=haversion,
                                             kord=kord, tord=tord, ford=ford,
                                             channels=channels)
        srcpath = os.path.abspath(filename)
        destfilepath = os.path.join(self.path_, self.SHAREDIR, filename)
        shutil.move(srcpath, destfilepath)

    def _search_coeff_files(self):
        coeffilemetas = []
        modelspecdir = os.path.join(self.path_, self.SHAREDIR)
        ls = os.listdir(modelspecdir)
        for lsentry in ls:
            if lsentry.endswith('.cc'):
                modeltype = 'Hamaker'
                ccfilename = lsentry.rstrip('.cc')
                version, band = ccfilename.split('Coeff')
                coeffilemetas.append({'band': band, 'modeltype': modeltype,
                                      'version': version})
        return coeffilemetas

    def get_bands(self):
        coeffilemetas = self._search_coeff_files()
        bands = []
        for coeffilemeta in coeffilemetas:
            bands.append(coeffilemeta['band'])
        return bands

    def list_models4band(self, band):
        coeffilemetas = self._search_coeff_files()
        models = []
        for coeffilemeta in coeffilemetas:
            if coeffilemeta['band'] == band:
                modelstring = '{}{}{}'.format(coeffilemeta['modeltype'],
                                              self._mdlversep,
                                              coeffilemeta['version'])
                models.append(modelstring)
        return models

    def list_feednames(self):
        coeffilemetas = self._search_coeff_files()
        feednames = []
        for coeffilemeta in coeffilemetas:
            feedname = self._mdlversep.join([coeffilemeta['band'],
                                             coeffilemeta['modeltype'],
                                             coeffilemeta['version']])
            feednames.append(feedname)
        return feednames

    def list_versions4bandmodels(self, band, modeltype):
        coeffilemetas = self._search_coeff_files()
        versions = []
        for coeffilemeta in coeffilemetas:
            if coeffilemeta['band'] == band and coeffilemeta['modeltype'] == modeltype:
                version = coeffilemeta['version']
                versions.append(version)
        return versions

    def _get_ccfile(self, band, haversion):
        inpfile = haversion+"Coeff"+band+".cc"
        return inpfile

    def _get_dpfile(self, band, version):
        outfile = band+self._bndversep+version+self.dpesuffix
        return outfile

    def list_dpes(self):
        dpes = []
        dpedir = os.path.join(self.path_, self.DATADIR)
        ls = os.listdir(dpedir)
        for lsentry in ls:
            if lsentry.endswith(self.dpesuffix):
                modeltype = 'Hamaker'  # FIXME:
                dpefilename = lsentry.rstrip(self.dpesuffix)
                band, version = dpefilename.split(self._bndversep, 1)
                dpes.append({'band': band, 'modeltype': modeltype,
                             'version': version})
        return dpes

    def load_dpolel(self, band, modelstring):
        modeltype, version = modelstring.split(self._mdlversep, 1)
        if version is None:
            versions = self.list_versions4bandmodels(band, modeltype)
            version = versions[0]
        if modeltype == 'Hamaker':
            dpfile = self._get_dpfile(band, version)
            if os.path.exists(dpfile):
                dpolel = pickle.load(open(dpfile, 'rb'))
            else:
                dpolel = self._gen_dpolel_hamaker(band, version)
        return dpolel

    def _gen_dpolel_hamaker(self, band, haversion):
        """Generate dualPolElem for the given band.
        Reads the 'lofar_elem_resp' packages c++ header files of
        Hamaker-Arts coefficients and writes pickled instances of DualPolElem
        class.
        Also adds nominal LOFAR frequency channels."""
        inpfile = self._get_ccfile(band, haversion)
        outfile = self._get_dpfile(band, haversion)
        inppath = os.path.join(self.path_, self.SHAREDIR, inpfile)
        outpath = os.path.join(self.path_, self.DATADIR, outfile)
        stndpolel = convLOFARcc2DPE(inppath, outpath)
        return stndpolel
