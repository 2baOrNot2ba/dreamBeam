"""rt (i.e. Radio Telescopes) module is for handling real telescope meta-data."""
import os
import glob
import pickle
import dreambeam.telescopes


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


class TelescopesWiz():
    """Database over available telescopes patterns."""
    def __init__(self):
        #Register telescope experts
        self.telescopes_dir = os.path.dirname(dreambeam.telescopes.__file__)
        ls = os.listdir(self.telescopes_dir)
        ds = []
        for p in ls:
            if os.path.isdir(self.telescopes_dir+'/'+p):
                ds.append(p)
        self.tbdata = {}
        for dd in ds:
            tbdata_dir = self.telbndmdl2dirfile(dd, '', '')[0]
            if os.path.isdir(tbdata_dir):
                t = os.path.basename(dd)
                self.tbdata[t] = {}
        #Find bands & models per telescope
        for tel in self.tbdata.keys():
            teldat_path = self.telbndmdl2dirfile(tel, '', '')[0]
            tbfiles = glob.glob(teldat_path+'*_*.teldat.p')
            bands = []
            antmodels = []
            for tbfile in tbfiles:
                filename = os.path.basename(tbfile)
                (band, modelsuffix)=filename.split('_')
                bands.append(band)
                antmodel=modelsuffix.split('.',2)[0]
                antmodels.append(antmodel)
            self.tbdata[tel] = {}
            for band in bands:
                self.tbdata[tel][band] = {}
                for antmodel in antmodels:
                    self.tbdata[tel][band][antmodel] = {}
        #Find stations for telescope band antmodel:
        for tel in self.tbdata.keys():
            for band in self.tbdata[tel].keys():
                for beammodel in self.tbdata[tel][band].keys():
                    telbnddata = self.getTelescopeBand(tel, band, beammodel)
                    self.tbdata[tel][band][antmodel] = telbnddata['Station'].keys()
        
    def get_telescopes(self):
        return self.tbdata.keys()

    def get_bands(self, telescope):
        return self.tbdata[telescope].keys()
    
    def get_stations(self, telescope, band):
        abeammodel = self.tbdata[telescope][band].keys()[0]
        return self.tbdata[telescope][band][abeammodel]

    def get_beammodels(self, telescope, band):
        return self.tbdata[telescope][band].keys()
    
    def getTelescopeBand(self, tscopename, band, beammodel):
        tbdata_dir, tbdata_fname  = self.telbndmdl2dirfile(tscopename, band, beammodel)
        tbdata_path = tbdata_dir+tbdata_fname
        with open(tbdata_path,'rb') as f:
            telbnddata = pickle.load(f)
        return telbnddata

    def telbndmdl2dirfile(self, tscopename, band, beammodel):
        """Map tscopename, band, beammodel tuple to file-path. file-path is a tuple
        of (absolute_directory, filename), where
           absolute_directory=/path-to-telescopes/TELESCOPENAME/data/
        and
           filename BAND_MODEL.teldat.p"""
        metadata_dir = "data/" #subdir within telescope dir with telbnd metadata. 
        #Currently it only maps requests to filename
        tbdata_fname = band+"_"+beammodel+".teldat.p"
        tbdata_dir = self.telescopes_dir+"/"+tscopename+"/"+metadata_dir
        return tbdata_dir, tbdata_fname
    
