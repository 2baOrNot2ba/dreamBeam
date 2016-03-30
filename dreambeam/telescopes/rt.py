import dreambeam.rime.jones

class TelescopeStnBnd(object):
    """Model of one station and one band of a telescope."""
    def __init__(self, stnDPolel, stnPos, stnRot):
        self.stnDPolel = stnDPolel
        self.stnPos = stnPos
        self.stnRot = stnRot
    
    def getEJones(self):
        """Create ejones for station based on antenna model."""
        ejones = dreambeam.rime.jones.EJones(self.stnDPolel, self.stnPos, self.stnRot)
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

