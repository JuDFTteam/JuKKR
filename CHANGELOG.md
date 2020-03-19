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

## *UNRELEASED* (last updated: 2020-03-19)

**Here we collect the list of *added*, *changed*, *deprecated*, *removed* and *fixed* features in preparation for the next release.**

### Added
- `EFSET` option for voronoi: read in whished value of Fermi level from inputcard (core state energies of starting potential is shifted accordingly)
- `<AFAC_RAD>` input for voronoi: read in a factor of radial mesh (r(i)=b*(exp(a*(i-1))-1)) where b is determined from a and rmt. Default value (a=0.025) is obtained with setting AFAC_RAD to a negative number. If a is increased (decreased) the start of the radial mesh is shifted to a lower (higher) radius.
- `POT_NS_CUTOFF` option for KKRhost and KKRimp (see issue #124)

### Changed
- cutoff of potential components which are (almost) zero (see issue #124 and `POT_NS_CUTOFF`)

### Deprecated
- None

### Removed
- None

### Fixed
- None


## v3.4 (2019-08-30)

### Added
- allow up to LMAX=8 in voronoi code and KKRhost code

### Changed
- Renamed keywords that have `files` in the name to avoid clash with `FILES` keyword in inputcard
- if GGA is used, set Vxc of empty cells to LDA (VWN)
- remove a lot of unnessecary files written out by each rank running the KKRimp code. Now by default only the master writes the files (old behavior can be reactivated using the `write_all_ranks` test flag)

### Deprecated
- `IMPURITY` option of voronoi code, not working properly and thus commented out

### Fixed
- Fix auto tests
- Small bugfixes

----

## v3.3 (2019-04-15)

Removed a lot of code duplicates among different codes.

### Added
- script to check the wronskian (see PhD Bauer, p.48)
- calc_wronskian option in KKRhost code to check single-site wavefunctions
- auto tests for OpenMP and MPI functionality
- xsf maker utility to create vesta input file
- improved error handling of inputcard reading errors
- generate code coverage report

### Changed
- more routines and modules in source/common directory
- use static library for pkkprime file among different apps
- use static library for kkrhost routines
- use routines from common/radial_solver_Chebychev in KKRimp as well
- rename lmsize to lmmaxd0 whenever the value is fixed (orbitalmoment.f90, rhoqtools.F90, rhovalnew.F90, tmat_newsolver.F90, rllsllsourceterms.f90)
- replace lmmaxso with lmmaxd and old lmmaxd with lmmmax0d wherever no spin-doubling occurs. Now lmmax0d=(lmax+1)**2 and lmmaxd=(1+krel+korbit)*(lmax+1)**2
- rename lmmax with lmmax0d or lmsize (in cases where the subroutine is called with different matrix sizes)
- use mpiatom parallelization scheme by default if NATYP>=IELAST
- moved wronskian from KKRimp to common to be used in KKRhost as well
- added lmdos writeout mode for qdos
- find NPRINCD to lowest possible divisor of NAEZ

### Removed
- removed code duplicates

### Fixed
- header handling in complexdos3 tool
- position of chebint in rhooutnew of KKRimp (probably) corrected
- bug in SRA-trick usage (issue #108)
- bug FERMIOUT option (issue #109)
- fix for issue #113
- fix for issue #114
- fix for issue #115
- fix for issue #116
- gfortran compilation fixed for kkrhost, kkrimp, voronoi
- fix debug compilation of kkrimp

----

## v3.2 (2018-12-11)

Improved memory consumption for XCPL calculation and small bugfixes.

### Changed
- greatly reduced memory consumption for XCPL calculation and improved performance of loops in shellgen2k

### Fixed
- issue #104 (voronoi-genpot mode broken for ASA)
- issue #105 (kkrflex_green wirteout for single-atom cluster)

----

## v3.1 (2018-11-21)

Start of large KKR repository holding *voronoi*, *KKRhost*, *KKRimp*, *KKRsusc*, and *PKKprime* with major refactoring of code structure.

### Added
- cmake installation scripts
- add rhoq repo (KKR-QPI)
- introduced new keywords in *KKRhost*: <INVMODE>, <VERBOSITY> and <MPI_SCHEME>, as well as new style of runoptions.

### Changed
- `install.py` script now deals with *KKRhost*, *KKRimp*, *PKKprime*, *voronoi*, and *rhoq*
- updated documentation
- restructured tests (see `gitlab-ci.yml` and `tests/gitlab-ci/*.yml` files)
- directory structure following `source`, `docs`, `utils`, `tests`, etc.
- mayor refactoring of the treatment of runoptions for *KKRhost* (new, descriptive keywords; no fixed format but style `runoption= T/F`; backwards compatibility is mostly ensured)
- escaped keywords (like `<ZATOM>` as opposed to `LMAX`) are now allowed to be case-insensitive (in both, `inputcard` and source code).
- refactoring of `source/common/ioinput.f90` for more simplicity, readability and flexibility
- total energy is also calculated and written out in case of non-scf calculation
- change factor (2-korbit) to (nspin-korbit) for NOSOC with NSPIN=1

### Deprecated
- makefiles of *PKKprime*, *voronoi*, *KKRimp*, *rhoq*
- duplicated files from *KKRhost* in *rhoq*
- common/test.f90 and common/opt.f90

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
- change `NOSOC` test option such that now is explicitly does a spin-loop with decoupled (i.e. a factor 2 smaller) matrices. The old behavior can be triggered by setting the `<SOCSCL>` values to 0.

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

### voronoi: 2016-08-18

Corrected a bug in the convolution of the starting potentials (concerning the
convolution of the 2Z/r term).

Added the possibility of CPA preparation (NATYP>NAEZ).
Then, the <SITE> keyword & column must be added next to the <ZATOM>.

Added input parameter <NFACELIM>= (integer) in order to be able to limit the
number of neighbours of the polyhedron that are examined for faster processing.

Moved all inputcard reads to readinput12 and readimpatoms12.
Other small points improved in the standard output.

----

### voronoi: 2016-06-20

The starting potentials are now convoluted fully with shapes (also the l>0
components). 
Gives improved results in some cases.
Does not apply to the case of interpolated potentials from previous KKR runs.

----

### voronoi: 2016-06-16

The starting potentials (l=0 component) are now convoluted with shapes.
Changes on this were done mainly in Subr. JELLSTART.
Does not apply to the case of interpolated potentials from previous KKR runs.

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

### voronoi: 2015-05-16

Added a print-out of the nearest-neighbor distance, dist(NN), per site and the ratio
Rout/dist(NN). Written out in file radii.dat. 

Attention, in case of shifted atoms, dist(NN) refers to the unshifted
positions. In order to find it for the shifted positions, re-run with the
parameter NBSHIFT=0 and also declaring "0" in the proper position under the
IMPINFO inplying there are no shifted positions. This will make the program
treat the shifted positions as if they were unshifted. 
The problem arises because the dist(NN) is calculated in the clsgen_voronoi
and clsgenimp routines that use the unshifted positions to calculate the 
voronoi cells, and the shift is introduced only later for the shape calculation.
The Rout on the other hand is defined from the shifted center.

----

## kkrhost-v1.3 (2015-04-29)

First version with MPI parallelization and Long's implementation of Lloyd's formula

### Added
- Long's version of LLOYD implementation
- MPI parallelization 

----

### voronoi: 2015-04-29

Change in subr. suggestpts.f fixing a problem that occured
if all panels had NMIN points. Then (NPAN-1)*NMIN points
were suggested, which was creating a problem in subr. mesh0.
Now at least (NPAN-1)*NMIN+1 points are suggested fixing the 
problem.
Added a write-out in subr. mesh0.

Added in subr. ritesone.f the option to write out in the potential
file a number of radial points larger than 999 (was formatted to 
i3, now format i4 is used for NR>999).

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

### voronoi: 2014-12-12

Fixed a bug with KAOEZ by setting NATYP --> NAEZ in readinput
(for old-type input of surf. calculations).

----

### voronoi: 2014-12-11

Introduced running option (RUNOPT) called "findsize" or "FINDSIZE". When this
is used, the weights of the atoms are automatically calculated from their
touching muffin-tin radii rmt (weight=rmt**2), where rmt is defined as the
half-distance of the atom to the nearest neighbour. Tested fot graphene on Co. 
Works also for impurity (not tested yet).
Option "findsize" overrides all other weight definitions, e.g. <MTWAL>.

Introduced keywords <LFMTWAL>, <RTMTWAL> and <LFMTWAU>, <RTMTWAU> for the
weights of the embedded atoms (outside the physical region) in slab
calculations.

Introduced keywords for reading in tolerance from inputcard:
<TOLHS> for subroutine halfplane;
<TOLVD> for minimum acceptable vertex distance;
<TOLAREA> for minimum acceptable face area.

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

### voronoi: 2014-10-29

readinput changed to comply with new inputcard keywords of JM code
(concerning the left/right basis in 2D mode).

----

## pkkprime-*untracked* (2014-10-24)

- Implemented writeout of Pkk' for option SCATTFIX and LLIFETIME.
- For option SCATTFIX, changed the implementation such that the INCOMMING k-vector is fixed (before it was outgoing k-vector).
- To solve memory-issues with very large slab calculations, the paramter `NROOTMAX` can be set in the inputcard. If `NROOTMAX= 0` assumes the maximal possible value.
- Included checks whether actually TBKKR_torq or TBKKR_rhod-files are present when they shall be used.

----

### kkrhost: 2014-10-13

`IQAT(NAEZD,NATYPD)` is changed to `IQAT(NATYPD)` because only the `IQAT(1,*)` was
ever used.

Array `NAT(NATYPD)` removed (was=1 always)

Changes in startb1.f to incorporate array fpradius.

ICC and IGREENFUN are automatically set to 1 if `OPT('KKRFLEX ')` is used.

----

### voronoi: 2014-10-10

Fixed a bug that was giving FPRADIUS not scaled with latt.parameter.

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

### voronoi: 2014-07-31

Corrected a bug arising from DSORT when sorting the polyhedron faces
at the end of VORONOI12 that resulted in giving slightly different shape
function mesh for symmetry-equivalent cells. 
Introduced new routine DSORT_COMP (in file
dsort.f) that sorts arrays according to increasing value of 1st component ,
then 2nd, then 3rd etc.

----

### kkrhost: 2014-07-25

Philipp's OMP parallelization is included. 
Lloyd's formula is included for non-relativistic non-cpa calculations.

----

### voronoi: 2014-07-25

Info on clusters and RMT of embedded region is now calculated and written
out. File atominfo.txt augmented with this information. With this, the
reference-system atom type in the embedded region is automatically found.

Proposed full-potential radius is written out in atominfo.txt 
(under keyword "<FPRADIUS>"), correspinding to IRNS of each atom. 

Corrected some bugs.

----

### voronoi: 2014-06-11

If the potential filename given in inputcard (under the FILES keyword, position
I13) is not found, then the program switches automatically to database starting
potentials (should be like this since 2 Oct 2013 but was not working). Fixed
by Benedikt Schweflinghaus in routine maindriver12.

Keyword FILES not needed any more. Then database starting potential is assumed. 
Change was done in routine rinput12.

----

### voronoi: 2014-06-05

Added a header at the beginning of each shape-function in the shape file.
E.g.
    4  135  NPAN,MESHN;  Shape number     1
instead of just
    4  135

----

### voronoi: 2014-03-18

Added parameter <ZATOM> for nuclear number. (Default=29)
Added parameter <RMTCORE> for (wished non-rouching) RMT. (Default=1.E10, reduced
automatically to touching-2% later)
These override ATOMINFO. 
Using these, ATOMINFO is made obsolete.

Also, a file "atominfo.dat" is created containing the ATOMINFO
information for the KKR code.

----

### voronoi: 2014-03-05

Most parameters are now default.
Only the following are needed:
ALATBASIS,BRAVAIS,NAEZ,RBASIS,RCLUSTZ

----

### voronoi: 2013-10-02

If potential filename is given in inputcard but not found in run-directory
then program does not stop but switches to jellium potentials automatically.

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

----

### voronoi: 2013-09-04

Added option to choose the muffin-tin radii of the atoms.
Place a keyword "<MTWAL>" (do not forget < and >) meaning 
"Muffin tin weight in alat units "
or "<MTWAU>" meaning
"Muffin tin weight in atomic units "
that represents the wished muffin-tin radius.
Under it should be all muffin tin radii (one under the other).
Most convenient place the whole column next to RBASIS.
Then these rmt's are squared and set as weights, w(i)=rmt(i)**2, 
and the weight(i) in atominfo is overridden. This results in a 
voronoi construction with these rmt radii.

For impurity atom weights: Here, in the place reserved for the weight 
place the desired MT-radius in the same units as in the normal atoms.

----

### voronoi: 2013-03-31

Routine clsgen2000 was polished, now in file clsgen_voronoi.f
Calling of clsgen list changed (also in maindriver12).
In maindriver12: Dimension of RMTCL fixed to NSHAPED

New routines clsgen_tb, clustcomp_tb introduced.
These distinguish the clusters also based on the MT-radius
of the reference potential. For Voronoi this is irrelevant,
but it gives useful info for setting up the TB calculation.
For this new (optional) variables are introduced in inputcard:
RMTREF, LEFTMTREF, RIGHMTREF (the last two for 2D-systems).
Result on clusters with clsgen_tb written in standard output.

----

### voronoi: 2012-10-22

Bug fixed in maindriver12.f where in case of 2D
the TLEFT and TRIGHT basis were not rationalized
creating a problem with clusters.

Added keyword "CARTESIMP= " (T or F) in scalevecimp
(see inputcard_example)

Added keyword "TOLIMP= ..." in clsgenimp12
(see inputcard_example)

----

### voronoi: 2012-09-12

Two bugs fixed in writing out file "radii.dat".
No new version issued since the bugs did not affect
shape-function or potential calculation.

----

### voronoi: 2012-03-25

A subroutine "dividepanels" was added that divides large panels into
smaller ones keeping the exact same radial points. The purpose
is to reduce the panel size for tests with the new integral solver.
The size of the new panels is determined by "NSMALL" which one can place
as a keyword in the inputcard: e.g., NSMALL=10 makes panels of 10-points each.
If the original panels have left-over points after division, these are 
placed in the last of the new panels.


Also:
Parameters NMIN and NRAD (minimum number of points per panel and
muffintinization points) can now be read in from the inputcard.
If they are not, then default values are used.

----

### voronoi: 2012-03-02

Furhter simplifications:

Now the code recognises the necessary number of radial points in the 
shape-region before calculating the shapes and pre-sets the value NMESH
accordingly. This is done by finding the critical points once before
entering the subroutine "shape". A parameter is used for this, "DENPT"
(set as a DATA statement in maindriver12) which declares the wished density
of radial points in the shape-function region. The panels are then assumed
to have this density of points, unless they are too small, when they have
NMIN points.

Also the code suggests a full-potential radius per atom, controlled by the
DATA-statement parameter STARTFP (again in maindriver12). This is the ratio
of the starting radius for full-potential to the muffin-tin radius. From this,
the parameter IRNS is calculated.

The number of radial points in the muffin-tin region is also set as a
DATA-statement (NMT).

Other changes: INS=1 and KSHAPE=0 means calculation of geometrical
information, shifting the centers, finding the panels etc, but avoiding
the actual calculation of shape functions ("spherical" shape functions are
calculated instead, same as option "SIMULASA").

A few bugs were removed.

The inc.geometry file was rationalized, unnecessary parameters removed.

----

### voronoi: 2012-02-01

Additional option "SIMULASA" was added.
Fake shape functions are generated in order to
simulate an ASA calculation if full-potential
mode. Only the (0,0) component of the shape
function is written out in a panel that extends
from Rmt until Rws (ASA radius).


----

### voronoi: 2012-01-30

The new version has a couple of new features:

1. It can give the shape function for impurity clusters.
2. It can expand the shape functions around "shifted" positions,
    different from the cell center, both for impurities and for host.

When in impurity mode, the code does not write out the host shape
functions, only the impurity ones.

In the directory there is a file "inputcard_example" where
the new features are explained (look for comments preceded by
the symbols #####). From the inputcard_example I have
omitted all input that would be redundant in a voronoi
calculation. In an actual calculation you can use directly
the inputcard that you prepared for the KKR code, just as before.

The impurity atoms are read in from the inputcard, just as the
host atoms. Older inputcards should be compatible with this.


A short description of further changes in the code follows:

*** Changes concerning the end-user ***

- The read-in of the inputcard was rationalised, omitting all redundant 
parameters. In some cases default values are used if a keyword is not
found (e.g., cartesian=.f. or kshape=2 are default).

- If you want to perform only a cell construction but not a calculation of
the shape functions, set KSHAPE=2. (kshape=1 is equivalent to =2 for
the voronoi prog. but not for the host kkr code). For full calculations put
INS=1 and KSHAPE=1 or 2.

- You can now write out the shapes for all atoms, as required by David's
code, or the shapes of only the representative atoms, as required by
the KKR host-code. Just use running option "WRITEALL" (with obvious meaning).

-The shape header is now included automatically by the program in the shape-file.

- Two new test-options were introduced, "verb0" and "verb1" (to be placed
in the TESTOPT part of the inputcard).
If you use none, then only little information is written out in the standard output,
compared to what was happening so-far. "verb0" gives more, and "verb0" together
with "verb1" even more.

- New files are written out to help some checks, plots, or post-processing:
verices.dat     for plotting the cells e.g. with gnuplot, 
radii.dat         for comparing the ratio of inner to outer radius, shows if the cell is too asymmetric,
cellinfo.dat     is the full set of face-equations and vertices, to be read in e.g. by a future
                      version of the KKR program so that the shapes are constructed on-the-fly.


*** Changes concerning the developer ***

- A routine creating the impurity clusters was added (clsgenimp12).
It accounts for host positions, impurity positions, and "killed" sites,
i.e. host sites that are to be completely removed in the impurity
calculation (e.g. when substituting more than one host atoms by only one
impurity).

- The basis vectors are rationalised by a new routine. If a basis atom
is too far away from the (0,0,0) lattice position, and close to another lattice
position, it is translated by a lattice vector back close to (0,0,0).
This was important for considering the impurity clusters.

- The write-out of the shape functions was moved from subr. "mtmesh"
into routine "writeshape" which is called from the main program.

- All shapes are kept into the main memory, so that the code is easier
to handle internally. As a result there is a ceiling of approx. 2000 different
shape functions, then the required memory exceeds 2GB. This can be 
improved and perhaps doubled by a reduction of array-size parameters,
in particular IBDMAXD (beginning of maindriver12) and NFACED,NVERTD.

- The subroutine shape.f (now called shape12.f) has only minor modifications:
the NFACED and NVERTD parameters have been moved to the inc.geometry
file so that they are uniquely defined in the code, and some write-out statements
are only given with "verb0" or "verb1" into the standard output (but I saw to it
that error messages are not restrained by "verb0" or "verb1").

