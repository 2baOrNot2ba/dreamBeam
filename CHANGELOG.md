# Changelog
All notable changes to this project will be documented in this file.


## [0.5] - 2020-01-30
### Added
- New base class `TelWizHelper` that is meant to be inherited for each
  telescope plug-in.
- New module `polarimetry` with function `convertxy2stokes()`.

### Changed
- Moved `TelescopeBndStn` to `feeds` module.
- `TelescopesWiz` now uses `pkgutil` to find telescope plug-ins.
- Moved `_get_teldat_fname` and `_get_telbnddatadir` from `rt` to telescope
 `telwizhelper`s.
- Renamed `LofarFeedJones` to `FixedMountEJones`
  and `LOFAR_LHBA_stn` to `FixedMountStn`.
- Moved `Feeds` module out of `LOFAR` and put it directly under `telescopes`.
- Renamed `beamfov()` in scenarios modules to `primarybeampat()`.
- Moved `read_LOFAR_HAcc()`, `convLOFARccHA()`, `convHA2DPE()` functions to
  `AntPat` package.
- Pickled instances are now suffixed with `.pkl`.
- Moved `getTelescopeBand` method out `TelescopeWizard` to become a function in
  `rt` and renamed it `open_telescopebndmodel`.
- `scenarios` does not use `TelescopeWizard` to get a feed model.
- Moved `save_telescopeband` to 'rt' and generalized it.
- Major refactor of `telescope` package, mainly to facilitate ingestion of new
  telescope plugins via `telwizhelper.py` scripts. `telwizhelper.py` scripts
  can now be quite minimal.
- Made `scripts`, under project root, a package so that modules there can be
  used as console_scripts entry points.


## [0.4] - 2019-09-06
### Added
- New method in `Jones` class called `sph2lud3_basis` that converts spherical
  basis vectors to the Ludwig3 basis. It has an optional `alignment` rotation
  argument that allows the p/q directions of the dualpol antennas to be
  rotated w.r.t. to the reference spherical system.
- New class `LofarFeedJones` (subclass of `EJones`) to handle LOFAR
  specific details of its `EJones`. Basis of final Jones is now the unit
  vectors of the X/Y dipole vectors projected onto the plane transverse to the
  pointing direction, using a conversion to Ludwig3.
- New class `LOFAR_LHBA_stn` (subclass of `TelescopeBndStn`) which is a merger
  of `LOFAR_LBA_stn` and `LOFAR_HBA_stn` classes (which were identical/redundant
  to each other).
- `CHANGELOG.md` file.
- Another way of specifying pointing in FoV_jones.py script: AZEL or STN frame,
  in addition to previous J2000.
- Direction cosine map argument to `beamfov()` in `scenarios` module. This means
  that primary beam corrections to interferometric images can be done in drB.
- New method `convert2iau()` in `Jones` class, to allow enforcing IAU coord sys
  of a Jones matrix.
- `inverse()` function for `Jones` objects, which computes the inverse of the
  Jones matrix. It is used e.g. to compute sky coordinates from pointings in
  station coordinates in `PJones` objects.
- New `diagnostics` module, for analyzing results.

### Changed
- New interpretation of the no-pararot algorithm: instead of computing in terms
  of an alt-az basis, the algorithm uses the actual antenna basis. This becomes
  the usual no-pararot for alt-az mounts but is different for a fixed mount
  such as LOFAR.
- Moved `display_pointings()` to `diagnostics` module. It has been heavily
  reworked to plot sph grid sys and bases, and support 3D option.
- Refactored `Jones` class to make it simpler (`jonesmeta` and `jonesmetar`
  dict data, now full attributes) and more python naming style (`get_basis()`
  instead of `getBasis()`).


## [0.3] - 2019-08-19
### Added
- New argument `--fmt` for `pointing_jones.py` script that controls print
  output. Can have values: `csv` (default) or `pac` (compatible)
- New argument `--no-pararot` and conversely `--pararot` (default) for
  `pointing_jones.py` script that either does not or does (resp.) apply the
  parallactic rotation.

### Changed
- Flipped X/Y signs in the LOFAR antenna build to station rotation.
  This results in a overall sign change of final Jones matrix (all components
  flip sign).


## [0.2] - 2019-04-14
### Added
- Flag in `PJones()` to turn off parallactic rotation.

### Changed
- Algorithm for parallactic rotation. Previously `PJones()` returned Jones
  in geocentric ITRF, and the ITRF2stn rotation was applied to the feed_pat
  (using its `rotateframe()`) in `EJones()`. Since this calculation was done in
  ITRF, it was difficult to access the parallactic rotation which is normally
  viewed as rotation from celestial frame to topocentric frame.
  Now `PJones` returns Jones in local, geocentric STN frame (topocentric) using
  the ITRF2stn rotation and the `EJones` is also computed in the STN frame.
- Start using C09 convention consistently for spherical basis vectors
  [Carozzi2009](https://doi.org/10.1111/j.1365-2966.2009.14642.x). This implied
  two major changes: the `DualPolFieldPointSrc` class is now defined in the C09
  basis and converts the Jones components accordingly, rather than as
  previously, keep the IAU components and convert the basis C09 to IAU.
- Mapping of the Jones columns (`AntPat` sph basis comps.) returned from
  `dualPolElem.getJonesAlong()` in `EJones()` to corresponding C09 sph basis
  components, rather than identity mapping.


## [0.1] - 2019-04-11
### Added
- Basic framework for project, consisting of a `rime` package that implements
  generic Radio Interferometric Measurement Equations, and a `telescopes`
  package which contains runtime pluggable modules for specific telescopes.
- `LOFAR` module under `telescopes`.
- print output from `pointing_jones.py` is now CSV.


## [0.0] - 2016-03-26
### Added
- Start using git as VCS for project.
