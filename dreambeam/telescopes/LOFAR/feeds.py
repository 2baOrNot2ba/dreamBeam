from dreambeam.telescopes.rt import TelescopeStnBnd
import dreambeam.rime.jones

class LOFAR_LBA_stn(TelescopeStnBnd):
    """Class for LOFAR LBA station."""
    def __init__(self, stnPos, stnRot):
        super(LOFAR_LBA_stn, self).__init__(stnPos, stnRot)
    def getEJones(self, pointing):
        """Get e-jones for this pointing for LOFAR LBA station.
        The basic model used here is that the antenna elements are all the
        same on all the stations, i.e. they have the same far-field pattern.
        Furthermore LOFAR are fixed mounted so 'pointing' LOFAR means
        electronic pointing. The simple pointing model used here for the LBA
        is to take the the response along the pointing-direction to be
        proportional to the response along the corresponding direction
        for the far-field pattern of the single antenna element."""
        ejones = dreambeam.rime.jones.EJones(self.feed_pat, self.stnPos, self.stnRot)
        return ejones
   # def __getstate__(self):
   #     print self.__dict__
   #     return self.__dict__
   # def __setstate__(self, pdict):
   #     print pdict
   #     self.__dict__ = pdict

class LOFAR_HBA_stn(TelescopeStnBnd):
    """Class for LOFAR HBA station."""
    def __init__(self, stnPos, stnRot):
        super(LOFAR_HBA_stn, self).__init__(stnPos, stnRot)
    def getEJones(self, pointing):
        """Get e-jones for this pointing for LOFAR HBA station.
        The basic model used here is that the antenna elements are all the
        same on all the stations, i.e. they have the same far-field pattern.
        Furthermore LOFAR are fixed mounted so 'pointing' LOFAR means
        electronic pointing. The simple pointing model used here for the HBA
        is to take the the response along the pointing-direction to be
        proportional to the response along the corresponding direction
        for the far-field pattern of the single antenna element."""
        ejones = dreambeam.rime.jones.EJones(self.feed_pat, self.stnPos, self.stnRot)
        return ejones

