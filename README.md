# dreamBeam #
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18433630.svg)](https://doi.org/10.5281/zenodo.18433630)

`dreamBeam` is a python package that implements radio astronomical measurement equations. It is based on the 3D Jones formalism introduced in:
> T. D. Carozzi, G. Woan, A generalized measurement equation and van Cittert-Zernike theorem for wide-field radio astronomical interferometry, Monthly Notices of the Royal Astronomical Society, Volume 395, Issue 3, May 2009, Pages 1558–1568, [(https://doi.org/10.1111/j.1365-2966.2009.14642.x)](https://doi.org/10.1111/j.1365-2966.2009.14642.x)


### What is this repository for? ###

`dreamBeam` let's you create generic measurement equations, i.e. multiplicative
chains of Jones matrices. It also comes with some predetermined, basic beam
models of real telescopes, currently:

* LOFAR (LBA & HBA)
* NenuFAR
* ALMA (contact me)

### How do I get set up? ###

* Use setup.py script to install
* Requires: python-casacore, numpy and [AntPat](https://github.com/2baOrNot2ba/AntPat)

### Status ###

`dreamBeam` is now fairly mature.
