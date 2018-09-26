This is the JM code containing Long's addition of David's routine for FP+SOC
and Philipp's addition of qdos qith SOC.


10.03.15

Fixed a bug in rhovalnew.f that was causing erroneous results when calculating
qdos.

22.01.15

Again: 
Fixed a bug in rinput13.f90 that did not allow the runopt XCPL to work
properly. (Was not working with NPOL=0)

11.01.15

Fixed a bug in rinput13.f90 that did not allow the runopt XCPL to work
properly.


13.11.14

Fixed a bug that would not allow using the REFPOT information for the
left/right region in 2D if ATOMINFO was used (the routine clsgen_tb was
finding its own rmtref).

Changed the calculation of rmtrefat in clsgen_tb to round up to 2 digits after
the decimal point.

6.11.14

inc.cls not needed any more.
Parameters NCLSD, NACLSD transfered to inc.p, with default NACLSD = NAEZD + NEMBD
Cluster info found automatically.


13.10.14

IQAT(NAEZD,NATYPD) is changed to IQAT(NATYPD) because only the IQAT(1,*) was
ever used.

Array NAT(NATYPD) removed (was=1 always)

Changes in startb1.f to incorporate array fpradius.

ICC and IGREENFUN are automatically set to 1 if OPT('KKRFLEX ') is used.


25.7.14

Philipp's OMP parallelization is included. 
Lloyd's formula is included for non-relativistic non-cpa calculations.

27.9.13

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

