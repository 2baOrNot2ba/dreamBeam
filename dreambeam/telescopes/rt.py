import dreambeam.rime.jones

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
    telescopeDir = "./"
    def __init__(self, tscopename, band, beammodel):
        mapname2filen(tscopename, band, beammodel)
        telescopestnbnd = pickle.load(open(filename,'rb'))
        return telescopestnbnd
    
    def mapname2filen(tscopename, band, beammodel):
        #Currently it only maps requests to filename
        filename = telescopeDir+tscopename+'/tel_'+band+beammodel+'.p'
        return filename

