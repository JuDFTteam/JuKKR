---
title: Changelog
---

# Changelog

## Conventions used for this changelog

 - keep it concise and human readable
 - keep the *UNRELEASED* section up to date with the `develop` branch
 - create a new subsection for each release version
 - each version should have the following information:
   - a release date in the format `YYYY-MM-DD`
   - a list of added new feature
   - a list of changed functionnality of existing features
   - a list of deprecated features (features that will be deleted in a future release)
   - a list of removed feature (previously marked deprecated)
   - a list of bug fixes

----

## *UNRELEASED* (last updated: 2018-10-16)

**Here we collect the list of *added*, *changed*, *deprecated*, *removed* and *fixed* features in preparation for the next release.**

Major code refactoring getting rid of the `inc.p` files which eliminates the need to recompile the code for different system sizes.

### Added
- `README`, `CHANGELOG`, `CONTRIBUTING` files and better documentation
- use cmake to auto-generate dependencies and makefile
- use FORD for automatic source code documentation
- readin of things for Bogoliubov-de-Gennes formalism
- test cases for ASA, Dirac, CPA, SOC, ...
- test option `noserial` to omit writing serial number to `kkrflex_*` files (for backwards compatibility)

### Changed
- convert all files to Fortran 90, including putting everything into modules
- polish code using NAG compiler, gfortran and ifort with debug options
- change 'NOSOC' test option such that now is explicitly does a spin-loop with decoupled (i.e. a factor 2 smaller) matrices. The old behavior can be triggered by setting the '<SOCSCL>' values to 0.

### Deprecated
- `inc.p` dependecy for array dimensions
-  old versions of the source code (`*.f` and `*.F` files) moved to old_src subdirectory

### Removed
- makefile (use cmake instead)

### Fixed
- gfortran compilation
- interfaces to Dirac routines (now also work in debug mode)

----

## v2.4 (2018-10-16)

Include Godfrin slab inversion, rhoq-writeout and more digits in alat writeout.

### Added
- output for QPI tool (activated with `rhoqtest` option)
- merged `godfrin` branch for sparse-matrix inversion in 2D mode
- added Options in `makefile` for `hazelhen` cluster in Stuttgart

### Changed
- increased number of digits in `alat` writeout of potential file

### Fixed
- bugfix in Dirac solver, introduced earlier due to code modernization
- writeout of energy-intrgrated DOS in last line of `dos.atom` files
- array size of `DEN` in `wrldos` routine corrected
- bugfix MPI communication in `rhovalnew` for `den_out` array

----

## v2.3 (2018-05-30)

Major changes merging branches `qdos_parallel`, and `develop_auto_tests` with some code cleanup.

### Added
- `RLL-SLL` option which speeds up computation of wavefunctions
- PBEsol GGA functional
- MPI parallel version of qdos
- `WRTGREEN` and `GREENIMP` run options 
- automatic tests (continuous integration with gitlab-ci)

### Changed
- nonco_angles_imp now has default values (theta=phi=0)

### Fixed
- GGA bugfix
- `NCLS=1` bug in `WRTGREEN` option
- case sensitivity of `CPP_HYBRID` preprocessor flag

----

## v2.2 (2017-04-21)

Some corrections, in particular change of sign in DMI vector calculation (with SOC).

### Added
- writeout files for FS code (run options `FERMIOUT`)

### Changed
- try to use either `md5` or `md5sum` commands to generate MD% checksum for potential and shapefun files
- use 3D Madelung sum also for 2D systems (only use 2D sum with option `ewald2d`). Needed since 2D sum seems to be buggy.
- sign change is changed to be consistent with KKRdoku, i.e. -Dij.(SixSj) is used

### Fixed
- MPI communication in XCPL + SOC
- MPI bug that set angles to 90 degree
- kkrflex writeout for DOS contour
- corrected readin for shapefunction if `LMAX` is smaller than `LMAX` of `shapefun` file
- idreals had a hard-coded system size
- orbital moments were not reinitialized in rhoval routines

----

## v2.1 (2016-09-08)

Some bugfixes and full Jij-tensor calculation (using new solver). 

### Added
- `XCPL` option for full Jij tensor (newsolver with SOC)
- GGA (PBE) and LDA+U for both versions with and without SOC by Long
- storing of wavefunctions
- writeout of MD5 sums for `potential` and `shapefun` files
- automatically guess if `MPIatom` or `MPIenerg` parallelization scheme should be used

### Changed
- removed common blocks for `optc` and `testc`, now in derived data type
- Jij energy dependent files only written for `NPOL=0` or `Jijenerg` test option

### Fixed
- writeout for `green` file of `XCPL` corrected
- MPI communication of `nonco_angles` fixed
- initialization in LLY for more than one iteration

----

## v2.0 (2016-03-24)

Major improvement to MPI parallelization.

### Added
- 2-level MPI parallelization
- MPI parallelization for Lloyd, Dirac, and exchange coupling constants
- hybrid version of rllsll
- writeout verbosity levels

### Changed
- `bzirr3d` can now use up to 500x500x100 kpoints
- change accuracy for kpoints (now more than single precision is kept)
- include compile flags and libs to version info

### Fixed
- t-matrix writeout for KKRimp was in done in local frame (needed to be in global frame)

----

## v1.4 (2015-11-13)

Improvement to MPI parallelization and new version of Lloyd.

### Added
- ouput moved to it's own file (`output.myrank.txt`, unit 1337)

### Changed
- change some static to allocatable arrays (makes -mcmodel=medium flag obsolete)
- Updated version of Long's implementation for LLOYD (based on `codeLLY_10_08_2015`)

### Fixed
- bugfix for serial run
- LLOYD with MPI 
- energy MPI parallelization for qdos
- MPI parallelization CPA

----

## v1.3 (2015-04-29)

First version with MPI parallelization and Long's implementation of Lloyd's formula

### Added
- Long's version of LLOYD implementation
- MPI parallelization 

----

## v1.2 (2015-04-21)

Got rid of unformatted files for data transfer between 1a, 1b, 1c parts of the code.

### Added
- storing of gmat, gref, tmat, in t_tgmat

### Changed
- renamed `fort.37` to `tmat.qdos` for qdos option

----

## v1.1 (2015-04-17)

First version with a single executable.

### Added
- Single executable with common makefile

----

## v1.0 (2015-04-15)

Started version control with git based on KKR code `JMcode_2015_03_10`
This version is the JM code containing Long's addition of David's routine for FP+SOC
and Philipp's addition of qdos with SOC. All changes mentioned below are included.

----

## Old versions (before version control with *git*)

### 2015-03-10

Fixed a bug in rhovalnew.f that was causing erroneous results when calculating
qdos.

----

### 2015-01-22

Again: 
Fixed a bug in rinput13.f90 that did not allow the runopt XCPL to work
properly. (Was not working with NPOL=0)

----

### 2015-01-11

Fixed a bug in rinput13.f90 that did not allow the runopt XCPL to work
properly.

----

### 2014-11-13

Fixed a bug that would not allow using the REFPOT information for the
left/right region in 2D if ATOMINFO was used (the routine clsgen_tb was
finding its own rmtref).

Changed the calculation of rmtrefat in clsgen_tb to round up to 2 digits after
the decimal point.

----

### 2014-11-06

inc.cls not needed any more.
Parameters NCLSD, NACLSD transfered to inc.p, with default NACLSD = NAEZD + NEMBD
Cluster info found automatically.

----

### 2014-10-13

`IQAT(NAEZD,NATYPD)` is changed to `IQAT(NATYPD)` because only the `IQAT(1,*)` was
ever used.

Array `NAT(NATYPD)` removed (was=1 always)

Changes in startb1.f to incorporate array fpradius.

ICC and IGREENFUN are automatically set to 1 if `OPT('KKRFLEX ')` is used.

----

### 2014-07-25

Philipp's OMP parallelization is included. 
Lloyd's formula is included for non-relativistic non-cpa calculations.

----

### 2013-09-27

rinput99 is changed to rinput13.
- everything useless was thrown out
- for most variables, default values are used (given in rinput13 routine)
- file inputcard_generated is created contains most read-in values

ioinput now can read in 5000 lines x 200 columns and symbols "<" and ">" 
in addition to capitals, numerals and "-".

Keyword LAMBDA_XC (default=1.) is introduced to mix the magnetic part of the
xc potential. LAMBDA_XC=0 corresponds to non-magn. calculation, LAMBDA_XC=1
to magn. calculation, and 0<LAMBDA_XC<1 to suppression of moments.
Result of xc-energy difference writter out as EXCDIFF.
