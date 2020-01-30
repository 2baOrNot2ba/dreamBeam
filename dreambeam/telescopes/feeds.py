import dreambeam.rime.jones


class FixedMountEJones(dreambeam.rime.jones.EJones):
    def get_basis(self):
        jb_sph = self.jonesbasis
        feedbasis_stn = self.dualPolElem.basis
        jb_lud3 = self.sph2lud3_basis(jb_sph, feedbasis_stn)
        return jb_lud3


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


class FixedMountStn(TelescopeBndStn):
    """Class for LOFAR LBA or HBA station."""

    def __init__(self, stnPos, stnRot):
        super(FixedMountStn, self).__init__(stnPos, stnRot)

    def getEJones(self, pointing, freqsel=None):
        """Get e-jones for this pointing for LOFAR LBA or HBA station.
        The basic model used here is that the antenna elements are all the
        same on all the stations, i.e. they have the same far-field pattern.
        Furthermore LOFAR are fixed mounted so 'pointing' LOFAR means
        electronic pointing. The simple pointing model used here for the LOFAR
        bands is to take the the response along the pointing-direction to be
        proportional to the response along the corresponding direction
        for the far-field pattern of the single antenna element."""
        ejones = FixedMountEJones(self.feed_pat, self.stnPos, self.stnRot,
                                  freqsel)
        return ejones
