##################################################################
# for iff371c... cluster use
#FC = ifort
#FFLAGS = -O2 -c -fomit-frame-pointer -fno-alias -r8
#LLIBS = -L/work/iff_th3/fmavrop/LIB/ifort_ATHLONSSE1_2 -Wl,-rpath=/usr/local/intel/fc/current/lib -lpthread -llapack -lptcblas -lptf77blas -latlas -lfxdropt
#LD = ifort
#
 FC = ifort
# LD = $(FC)
# FFLAGS = -O2 -check all -r8 
# FFLAGS = -O2 -r8 
#FFLAGS = -O2 -r8 -mcmodel=medium -i-dynamic -check all -traceback
FFLAGS = -O2 -r8 -mcmodel=medium -i-dynamic
# LLIBS = -L /usr/local/intel/mkl/current/lib/32 -lmkl -lguide -lpthread -i_dynamic -L /work/iff_th3/fmavrop/LIB/fxdr -lfxdropt -Wl,-rpath=/usr/local/intel/mkl/current/lib/32 -Wl,-rpath=/usr/local/intel/fc/current/lib
 LLIBS = -L /usr/local/intel/mkl/lib/ia32 -lmkl -lguide -lpthread -i_dynamic -L /work/iff_th3/fmavrop/LIB/fxdr -lfxdropt -Wl,-rpath=/usr/local/intel/mkl/lib/ia32 -Wl,-rpath=/usr/local/intel/fc/current/lib
##################################################################
 LLIBS = -L/usr/local/nag/lib -lnag -llapack_ifort -lblas_ifort 
# LLIBS = -L /usr/users/iff_th1/heers/Programme -lfxdr -L/usr/local/nag/lib -lnag -llapack_ifort -lblas_ifort 
##################################################################
#LLIBS= -L /usr/users/iff_th1/heers/Programme -lfxdr -L/usr/local/intel/current/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
#LLIBS= -L /usr/users/iff_th1/heers/Programme -lfxdr -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
LLIBS= -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
# for Swantjes macbook:
#FC = ifort
#FFLAGS = -O3 -c -fomit-frame-pointer -fno-alias -r8 -heap-arrays
#LLIBS = -L/usr/local/lib -lpthread -llapack -lblas -latlas -lfxdr -L /Users/swantjeheers/Documents/PhD/fxdr/library_fxdr/
#LD = ifort

PROGRAM = kkrscatter.exe

dir=.


host = `hostname`


#matrixESSL = dgetrf.o dgetrs.o zgetrf.o zgetrs.o zgeevy.o zhpevy.o

#matrixLAPACK = zgeevy.o zhpevy.o

#laplib1= dgetf2.o  dgetrs.o  ieeeck.o  zgetf2.o \
#dgetrf.o  dlaswp.o  ilaenv.o  zgetrf.o  zgetrs.o  zlaswp.o zgeev_plus.o

## for Fermi surface link bandinput_FS3.o kkrmateigen_FS.o

HPfiles = derfc.o

#ggafiles= vxcgga.o gradr.o gradrl.o mkxcpe.o gxcpt.o \
exch91.o corlsd.o cpw91.o gcor91.o spherenew.o

files    =  rinput99.o \
rcstop.o \
ioinput.o \
opt.o \
integ.o \
csimpk1.o \
csimp3.o \
testdim.o \
lattix99.o \
ddet33.o \
cross.o \
vadd.o \
veq.o \
crospr.o \
spatpr.o \
rrgen.o \
vmul.o \
scalpr.o \
scalevec.o \
test.o \
clsgen99.o \
dsort.o \
clustcomp.o \
gaunt.o \
cnwf01.o \
gauss.o \
corelb.o \
intcor.o \
intout.o \
intin.o \
breit.o \
spherenew.o \
lebedev.o \
inout.o \
rinit.o \
cinit.o \
calrmt.o \
calctref.o \
potcut.o \
emesht.o \
deltamat.o \
expidl.o \
tmatrx01_SO.o \
cplxwb01.o \
regsol.o \
irwsol.o \
regns_volterra.o \
regns_SO.o \
regns.o \
wfint.o \
wfint_SO.o \
wfint0.o \
wfint0_SO.o \
csinwd.o \
csinwd_SO.o \
csout.o \
csout_spline.o \
vllns.o \
vllmat.o \
wftsca.o \
wftsca_SO.o \
irwns.o \
gll95.o \
gfree.o \
grefsy.o \
ham_so.o \
spinorbit_ham.o \
spinorbit_ham_test.o \
spin_orbit.o \
spin_orbit_alle_l.o \
rot_rll.o \
normcoeff_SO.o \
normcoeff_SO_torq.o \
normcoeff_SO_spinflux.o \
calc_rho_ll_ss.o \
calc_torq_ll_ss.o \
calc_rho_ll_ss_lmdos.o \
sigma_drehung.o \
rotbrill.o \
pointgrp.o \
findgroup.o \
latvec.o \
length.o \
symlat.o \
setgroup.o \
sname.o \
bzirr3d.o \
crtstar.o \
gfshells.o \
rotmat.o \
dmpy.o \
dinv33.o \
shellgen2k.o \
dinit.o \
errmsg.o \
nrmliz.o \
ddot1.o \
dswap1.o \
dscal1.o \
gfmask.o \
inversion_SO.o \
sp2.o \
sp2_SO.o \
gllcopy1.o \
tautog1_SO.o \
invslab_SO.o \
invsupercell.o \
btom.o \
bofm.o \
main.o \
matrix_oper.o \
check_WF_real.o \
changerep.o \
cmatstr.o \
bastrmat.o \
drvbastrans.o \
calccgc.o \
transform_tmat.o \
transform_LSM.o \
transform_PZM_back.o \
potenshift1.o \
convol.o \
bandinput.o \
bandstr_SO.o \
calc_wapr_bstr.o \
rotatecoeff_SO.o \
calculate_lifetime_SO.o \
lifetime_2D_SO.o \
scattering_matrix_2D_SO.o \
scattering_matrix_2D_NOSOC.o \
scattering_matrix_2D_SO_ani.o \
read_coeff_2D_SO.o \
read_coeff_2D_SO_ani.o \
read_delta_tmat.o \
read_green_LL.o \

dyson= dlke0.o \
       dlke1.o \
       calc_wapr.o \
       calc_wapr_fv.o \
       calc_wapr_sof.o \
       tetraeder.o \
       volume_tetr.o \
       dgdk.o \
       gllcopya.o \
       zsort.o \
       zsort_complex.o \
       zsort_complex_abs.o \
       zgeev_plus.o \

newsol = tmat_newsolver.o \
	 create_newmesh.o \
	 interpolspline.o \
         beshank.o \
         chebint.o \
         readimppot.o \
         tmatimp_newsolver.o \
         writegreen.o \
	 writetbkkrdata.o \
         greenimp.o \
         getvolbz.o \
         gfshells1.o \

mode90 = rotatespinframe.o \
         ../common/DataTypes.f90 \
         ../common/splint.o \
         ../common/spline.o \
         vllmatsra.o \
         rllsllsourceterms.o \
         rllsll.o \
         calcsph.o \
         chebyshev.o \
         cheb2oldgrid.o 

SRC77 = ${files} ${dyson} ${newsol}
SRC90 = ${mode90}

SRC = ${SRC90} ${SRC77}

OBJS = ${SRC:.f=.o}
.SUFFIXES: .f
.f.o: 
	$(FC) $(FFLAGS) -c $<

OBJS = ${SRC:.f90=.o}
.SUFFIXES: .f90
.f90.o:
	$(FC) $(FFLAGS) -c $<

$(PROGRAM): $(OBJS) 
	$(FC) -o $@ $(OBJS) $(LLIBS)
${files}:	inc.p inc.cls
${dyson}:	inc.p

# Clean up
clean:
	rm -f *.mod *.o ../$(files) ../$(dyson) ../$(newsol)
