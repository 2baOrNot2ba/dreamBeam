"""Module to indicate that a feed plugin for NenuFAR can be generated from
here but first time it needs to initialize.
"""


if __name__ == "__main__":
    """Script to generate Hamaker model data files for NEC simulations of the
    NenuFAR antenna element.

    Only need to run this if there is no HamakerArts .cc files in ``share``
    folder. These should be provided in repo. If it needs to be recreated you
    will need to put an approriate NEC .out file in the ``share`` folder.
    In that case do the following in the folder containing this file:
    ::

        $ python _feeds.py
    """
    import numpy
    from dreambeam.feeds.feedplugins import FeedPlugin
    feedplugin = FeedPlugin(".")
    necfiles = feedplugin.find_necoutfiles()
    necfile = necfiles[0]
    SAMPFREQ = 100e6
    NR_CHANNELS = 512
    # Overriding channels for NEC simulation:
    channels = numpy.linspace(0., SAMPFREQ, NR_CHANNELS, endpoint=False)
    feedplugin.compute_ccfile(necfile, scalefac=-1/240.0, freq_center=67e6,
                              freq_range=32e6, kord=1, tord=1, ford=2,
                              channels=channels)
    # Default with same ords as LOFAR default:
    # feedplugin.compute_ccfile(necfile, scalefac=-1/240.0, freq_center=67e6,
    #                          freq_range=32e6, kord=2, tord=5, ford=5,
    #                          channels=channels)
    print("Completed setup.")
