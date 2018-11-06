# Changelog

## Conventions used for this changelog

 - keep it concise but human readable
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

## *UNRELEASED* (last updated: 2018-11-06)

**Here we collect the list of *added*, *changed*, *deprecated*, *removed* and *fixed* features in preparation for the next release.**

Start of large KKR repository holding *voronoi*, *KKRhost*, *KKRimp*, *KKRsusc*, and *PKKprime* with major refactoring of code structure.


### Added
- None

### Changed
- `install.py` script now deals with *KKRhost*, *KKRimp*, *PKKprime*, and *voronoi*
- updated documentation
- tests structure
- directory structure following `source`, `docs`, `utils`, `tests`, etc. 
- change factor (2-korbit) to (nspin-korbit) for NOSOC with NSPIN=1

### Deprecated
- makefiles of *PKKprime*, *voronoi*, *KKRimp*, *rhoq*
- duplicated files from *KKRhost* in *rhoq*

### Removed
- None

### Fixed
- None

----

## kkrimp-v1.3 (2018-10-31)

Some bugfixed and som new funcitonalities.

### Added
- simple version of Lloyd's formula: `LLYsimple` option
- use other XC functionals than VWN (other LDAs and PW91 GGA, no PBE yet!)
- auto tests with *gitlab-ci*
- ford documentation of source code

### Changed
- default behavior to use SRA-trick automatically
- relaxed comparison between old and new mesh in `cheb2oldgrid`

### Fixed
- doubling of allocations in rllsll
- bugfix Jijsymmetries (string of wavefunctions did not work properly)

----

## kkrhost-v3.0 (2018-10-30)

Major code refactoring getting rid of the `inc.p` files which eliminates the need to recompile the code for different system sizes.

### Added
- MIT license
- `README`, `CHANGELOG`, `CONTRIBUTING` files and better documentation
- use cmake to auto-generate dependencies and makefile
- use FORD for automatic source code documentation
- read-in of things for Bogoliubov-de-Gennes formalism
- test cases for ASA, Dirac, CPA, SOC, godfrin, LDA+U, XCs
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
- wrong symmetrization of shells for Jij calculation without inversino symmetry (issue #64)

----

## pkkprime-v1.1 (2018-10-24)

Major update of Fermi surface and scattering code.

### Added
- ford source code annotations
- masked integration for lifetime
- feast library
- optimized memory footprint
- MIT license
- different weights of multiple impurities
- torkance calculation

### Changed
- kfixstart==kfixstop option
- signchange spinflux

### Fixed
- some small bugs corrected

----

## kkrhost-v2.4 (2018-10-16)

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

## kkrhost-v2.3 (2018-05-30)

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

## rhoq-v1.1 (2018-01-30)

First working (beta) version of the rhoq module (KKR-QPI) functionality.
Needs writeout of files from host code (use `rhoqtest` test option).

### Added
- hybrid parallelization (MPI+OpenMP)
- timing and version info
- kmask selective integration
- exclude region (C_M-term)

----

## kkrhost-v2.2 (2017-04-21)

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

## kkrhost-v2.1 (2016-09-08)

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

## kkrhost-v2.0 (2016-03-24)

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

## kkrimp-v1.2 (2016-02-22)

New version tag, now everything is faster because wavefunctions are not
always recalculated and hybrid parallelisation (explicit OpenMP
parallelism in rllsll) is implemented.

### Added
- OpenMP parallel version of rllsll by Sachin

### Fixed
- bugfix usespinorbit

----

## rhoq-v1.0 (2016-02-01)

Initialized rhoq module and started development of KKR-QPI feature.

----

## kkrimp-v1.1 (2015-11-20)

Version which prints build information and has minor improvements

### Added
- version information

### Changed
- case sensitivity of run/test options removed
- default behavior of storing more wavefunctions

### Fixed
- bug in spinorbitperatom

----

## kkrhost-v1.4 (2015-11-13)

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

## kkrimp-v1.0 (2015-10-19)

First version of KKRimp that is tracked with git. Started from *kkrimp_source_2014_10_01*.

----

## kkrhost-v1.3 (2015-04-29)

First version with MPI parallelization and Long's implementation of Lloyd's formula

### Added
- Long's version of LLOYD implementation
- MPI parallelization 

----

## kkrhost-v1.2 (2015-04-21)

Got rid of unformatted files for data transfer between 1a, 1b, 1c parts of the code.

### Added
- storing of gmat, gref, tmat, in t_tgmat

### Changed
- renamed `fort.37` to `tmat.qdos` for qdos option

----

## kkrhost-v1.1 (2015-04-17)

First version with a single executable.

### Added
- Single executable with common makefile

----

## kkrhost-v1.0 (2015-04-15)

Started version control with git based on KKR code `JMcode_2015_03_10`
This version is the JM code containing Long's addition of David's routine for FP+SOC
and Philipp's addition of qdos with SOC. All changes mentioned below are included.

----

## Old versions of kkrhost (before version control with *git*)

### kkrhost: 2015-03-10

Fixed a bug in rhovalnew.f that was causing erroneous results when calculating
qdos.

----

### kkrhost: 2015-01-22

Again: 
Fixed a bug in rinput13.f90 that did not allow the runopt XCPL to work
properly. (Was not working with NPOL=0)

----

### kkrhost: 2015-01-11

Fixed a bug in rinput13.f90 that did not allow the runopt XCPL to work
properly.

----

## pkkprime-v1.0 (2015-01-07)

Initial version of git tracking.

--

## pkkprime-*untracked* (2014-11-19)

The output of the atom resolved torque for bulk systems is now implemented (but not tested though).

----

### kkrhost: 2014-11-13

Fixed a bug that would not allow using the REFPOT information for the
left/right region in 2D if ATOMINFO was used (the routine clsgen_tb was
finding its own rmtref).

Changed the calculation of rmtrefat in clsgen_tb to round up to 2 digits after
the decimal point.

----

## pkkprime-*untracked* (2014-11-10)

Implemented the atom resolved torque (set LTORQATOM=1 in the input file). The output for bulk systems is not yet implemented though.

----

### kkrhost: 2014-11-06

inc.cls not needed any more.
Parameters NCLSD, NACLSD transfered to inc.p, with default NACLSD = NAEZD + NEMBD
Cluster info found automatically.

----

## pkkprime-*untracked* (2014-10-24)

- Implemented writeout of Pkk' for option SCATTFIX and LLIFETIME.
- For option SCATTFIX, changed the implementation such that the INCOMMING k-vector is fixed (before it was outgoing k-vector).
- To solve memory-issues with very large slab calculations, the paramter 'NROOTMAX' can be set in the inputcard. If NROOTMAX= 0 assumes the maximal possible value.
- Included checks whether actually TBKKR_torq orr TBKKR_rhod-files are present when they shall be used.

----

### kkrhost: 2014-10-13

`IQAT(NAEZD,NATYPD)` is changed to `IQAT(NATYPD)` because only the `IQAT(1,*)` was
ever used.

Array `NAT(NATYPD)` removed (was=1 always)

Changes in startb1.f to incorporate array fpradius.

ICC and IGREENFUN are automatically set to 1 if `OPT('KKRFLEX ')` is used.

----

### kkrimp: 2014-10-01

Further corrections for running with lower lmax than given in
the host GF. Changes mainly in _dysonviratom.f90_ (again) and
in rhooutnew. Now the lmax-per-atom (lmaxatom) is given
to rhooutnew (before it was only lmaxd) and some according
changes were made in the routine. Works when all atoms
have the same lmax (lower-or-equal than the host), but crashes
in rhooutnew when lmax is different per atom.
However, different-lmax-per-atom runs in case of no-spin-orbit
(old radial solver). 


----

### kkrimp: 2014-07-31

Corrected a bug in _dysonviratom.f90_ that was causing problems 
when lmax of impurity was smaller than lmax of host and 
spin-orbit of host was on.

Small change in _cheb2oldgridc.f90_ for improved numerical accuracy.

Change in _rllsll.f90_ so that the subr. _inverse_ defined within this 
file is commented out, because the subr. _inverse_  of the module
_rllslltools.f90_ is used when _inverse_ is called from _rllsll_. 

Change of record length definition (recl) in direct-access files.
Defined parameter wlength=1 (or 4 depending on compiler) in _nrtype.f90_ 
and included this in all places where a direct-access file is opened.
The changed files are: 
gdyson.f90, utrafo.f90, wavefunctodisc.f90, energyloop.F90, preconditioning.F90

----

### kkrhost: 2014-07-25

Philipp's OMP parallelization is included. 
Lloyd's formula is included for non-relativistic non-cpa calculations.

----

### kkrhost: 2013-09-27

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
