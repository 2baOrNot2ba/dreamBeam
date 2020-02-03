import dreambeam.rime.jones


class FixedMountEJones(dreambeam.rime.jones.EJones):
    def get_basis(self):
        jb_sph = self.jonesbasis
        feedbasis_stn = self.dualPolElem.basis
        jb_lud3 = self.sph2lud3_basis(jb_sph, feedbasis_stn)
        return jb_lud3


# PosedFeed
class TelescopeBndStn(object):
    """Base model of a band feed on a telescope station."""
    feed_pat = None

    def __init__(self, stnPos, stnRot):
        """Set the station's position and attitude."""
        self.stnPos = stnPos
        self.stnRot = stnRot

    def getEJones(self):
        """Create ejones for station based on antenna patterns.
        Meant to be overriden.
        """
        ejones = None
        return ejones

    def getfreqs(self):
        freqs = self.feed_pat.getfreqs()
        return freqs


# PosedFeed_FixMnt
class FixedMountStn(TelescopeBndStn):
    """Class for fixed mount station."""

    def __init__(self, stnPos, stnRot):
        super(FixedMountStn, self).__init__(stnPos, stnRot)

    def getEJones(self, pointing, freqsel=None):
        """Get e-jones for this pointing for a fixed mount station.
        The basic model used here is that the antenna elements are all the
        same on all the stations, i.e. they have the same far-field pattern.
        For fixed mounted antennas 'pointing' means array based pointing, i.e.
        no mechanical slewing. The simple pointing model used here for
        is to take the the response along the pointing-direction to be
        proportional to the response along the corresponding direction
        for the far-field pattern of the single antenna element."""
        ejones = FixedMountEJones(self.feed_pat, self.stnPos, self.stnRot,
                                  freqsel)
        return ejones
