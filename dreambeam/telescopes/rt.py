"""rt (i.e. Radio Telescopes) module is for handling real telescope meta-data."""
import os
import pickle
import dreambeam.telescopes


class TelescopeStnBnd(object):
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


class TelescopeWiz():
    """Database over available telescopes patterns."""
    def __init__(self):
        self.telescopes_dir = os.path.dirname(dreambeam.telescopes.__file__)

    def list_telescopes(self):
        ls = os.listdir(self.telescopes_dir)
        ds = []
        for p in ls:
            if os.path.isdir(self.telescopes_dir+'/'+p):
                ds.append(p)
        telescope_list = []
        for dd in ds:
            if os.path.isdir(self.telescopes_dir+'/'+dd+"/data/"):
                telescope_list.append(os.path.basename(dd))
        return telescope_list

    def getTelescopeband(self, tscopename, band, beammodel):
        teldatapath = self.telbnd2path(tscopename, band, beammodel)
        with open(teldatapath,'rb') as f:
            telescope = pickle.load(f)
        return telescope
    
    def telbnd2path(self, tscopename, band, beammodel):
        #Currently it only maps requests to filename
        filename = "teldat_"+tscopename+"_"+band+"_"+beammodel+".p"
        teldatdir = os.path.dirname(dreambeam.telescopes.__file__)+"/"+tscopename+"/data/"
        teldatapath = teldatdir+filename
        return teldatapath

