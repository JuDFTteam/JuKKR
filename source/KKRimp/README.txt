2014.10.01

Further corrections for running with lower lmax than given in
the host GF. Changes mainly in _dysonviratom.f90_ (again) and
in rhooutnew. Now the lmax-per-atom (lmaxatom) is given
to rhooutnew (before it was only lmaxd) and some according
changes were made in the routine. Works when all atoms
have the same lmax (lower-or-equal than the host), but crashes
in rhooutnew when lmax is different per atom.
However, different-lmax-per-atom runs in case of no-spin-orbit
(old radial solver). 

2014.07.31

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
