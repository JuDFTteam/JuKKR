####################################################
####################################################
#         Julich KKR code makefile
####################################################
#
# to compile use the following options:
#
# make           -    standart compilation using
#                     ifort on a single core
# make mpi       -    for using MPI parallelization 
#                     with ifort
# make clean     -    removes all *.o files
# make openmp    -    uses openmp
#
default: install
####################################################
# specify the compiler which is going to be used   #
####################################################
###  default value  ###
FC =  ifort
# change default value if a flag is matched
mpi:       FC =  mpiifort
hybrid:    FC =  mpiifort
hybriddebug:    FC =  mpiifort
###  mac  ###
macmpi:     FC =  /usr/local/impi/bin/mpif90  
mpidebug:   FC = mpiifort
###  rwth cluster  ###
rwthmpi: FC=$(MPIFC)
rwthmpidebug: FC=$(MPIFC)
# marmotf90
rwthhybrid: FC=$(MPIFC)
###  jureca  ###
jureca: FC = mpif90
# $(MPIFC)


# information for version-tracking in executable:
# default
COMPVER ='serial'
# special cases
debug:     COMPVER ='debug'
mpi:       COMPVER ='mpi'
mpidebug:  COMPVER ='mpidebug'
hybrid:    COMPVER ='hybrid'
hybriddebug:    COMPVER ='hybriddebug'
openmp:    COMPVER ='openmp'
ompdebug:  COMPVER ='openmpdebug'
mac:       COMPVER ='serial-mac'
macmpi:    COMPVER ='mpi-mac'
macomp:    COMPVER ='openmp-mac'
macompprof:    COMPVER ='openmp-mac-prof'
rwth:      COMPVER ='serial-rwth'
rwthmpi:   COMPVER ='mpi-rwth'
rwthomp:   COMPVER ='openmp-rwth'
rwthmpidebug: COMPVER ='mpidebug-rwth'
rwthhybrid:   COMPVER ='hybrid-rwth'
jureca   : COMPVER='jureca'



####################################################
# Objects to compile                               #
####################################################


OBJ = ./OBJ
SRC = ./ALL_SOURCE_FILES

INCP = ./inc.p
INCR = ./sprkkr_rmesh.dim

FILENAME0 = kkr.x

MODULES = $(OBJ)/mod_version.o $(OBJ)/mod_version_info.o  $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJ)/mod_save_wavefun.o $(OBJ)/mod_wunfiles.o $(OBJ)/mod_md5sums.o $(OBJ)/mod_jijhelp.o $(OBJ)/initabjij.o $(OBJ)/tbxccpljij.o $(OBJ)/tbxccpljijdij.o $(OBJ)/create_newmesh.o

OBJECTS0 = $(OBJ)/getbr3.o $(OBJ)/write_tbkkr_files.o $(OBJ)/gaunt.o $(OBJ)/main0.o	$(OBJ)/rinput13.o $(OBJ)/symtaumat.o $(OBJ)/taustruc.o $(OBJ)/bzkint0.o $(OBJ)/findgroup.o $(OBJ)/calcrotmat.o $(OBJ)/errortrap.o $(OBJ)/bzkmesh.o $(OBJ)/bzirr3d.o $(OBJ)/madelung3d.o $(OBJ)/madelung2d.o $(OBJ)/madelout.o $(OBJ)/idreals.o $(OBJ)/deciopt.o $(OBJ)/impcheck.o $(OBJ)/impcoefs.o $(OBJ)/setgijtab.o $(OBJ)/gijxcpl.o $(OBJ)/gijcond.o $(OBJ)/calctref.o $(OBJ)/getclusnxyz.o $(OBJ)/decitset.o    $(OBJ)/decitmat.o $(OBJ)/outpothost.o  $(OBJ)/decipothead.o $(OBJ)/getscratch.o $(OBJ)/decipotbas.o  $(OBJ)/outtmathost.o $(OBJ)/lngstring.o $(OBJ)/startldau.o   $(OBJ)/readldaupot.o $(OBJ)/epathtb.o $(OBJ)/addvirtual14.o   $(OBJ)/interpolspline.o $(OBJ)/spline_real.o $(OBJ)/splint_real.o $(OBJ)/calrmt.o     $(OBJ)/clsgen_tb.o   $(OBJ)/clustcomp.o $(OBJ)/cmomsread.o $(OBJ)/crospr.o $(OBJ)/cross.o $(OBJ)/crtstar.o $(OBJ)/ddet33.o $(OBJ)/ddot1.o $(OBJ)/dscal1.o $(OBJ)/dsort.o $(OBJ)/dswap1.o $(OBJ)/emesht.o $(OBJ)/erfcex.o $(OBJ)/errmsg.o $(OBJ)/gamfc.o $(OBJ)/gaufd.o  $(OBJ)/gauleg.o $(OBJ)/gaunt2.o $(OBJ)/gfmask.o $(OBJ)/gfshells.o $(OBJ)/grule.o $(OBJ)/ioinput.o $(OBJ)/lattice2d.o $(OBJ)/lattice3d.o $(OBJ)/lattix99.o $(OBJ)/latvec.o $(OBJ)/length.o $(OBJ)/nrmliz.o $(OBJ)/opt.o $(OBJ)/pointgrp.o $(OBJ)/potcut.o $(OBJ)/rcstop.o $(OBJ)/rinit.o $(OBJ)/decimaread.o $(OBJ)/rrgen.o $(OBJ)/scalevec.o $(OBJ)/cinit.o $(OBJ)/scalpr.o $(OBJ)/shellgen2k.o $(OBJ)/sname.o $(OBJ)/spatpr.o $(OBJ)/strmat.o $(OBJ)/test.o $(OBJ)/generalpot.o $(OBJ)/ymy.o $(OBJ)/veq.o $(OBJ)/vadd.o $(OBJ)/vmul.o $(OBJ)/testdim.o $(OBJ)/startb1.o $(OBJ)/madelcoef.o $(OBJ)/madelgaunt.o $(OBJ)/ewald2d.o $(OBJ)/fplaneg.o $(OBJ)/fplaner.o $(OBJ)/shape.o $(OBJ)/wfmesh.o $(OBJ)/regsol.o $(OBJ)/beshan.o $(OBJ)/drvbastrans.o $(OBJ)/relpotcvt.o $(OBJ)/drvreltmat.o $(OBJ)/calccgc.o $(OBJ)/ikmlin.o $(OBJ)/strsmat.o $(OBJ)/bastrmat.o $(OBJ)/ikapmue.o $(OBJ)/cmatstr.o $(OBJ)/checkrmat.o $(OBJ)/rinvgj.o $(OBJ)/ssite.o $(OBJ)/cintabr.o $(OBJ)/sumupint.o $(OBJ)/dirac_bs.o $(OBJ)/dirbslag.o $(OBJ)/dirbsmid.o $(OBJ)/dirbsrad.o $(OBJ)/dirbsrze.o $(OBJ)/dirbsstp.o $(OBJ)/cjlz.o $(OBJ)/cnlz.o $(OBJ)/cdjlzdz.o $(OBJ)/cdnlzdz.o $(OBJ)/cint4pts.o $(OBJ)/cinthff.o $(OBJ)/ylag.o $(OBJ)/rnuctab.o $(OBJ)/dirac_bi.o $(OBJ)/dirac_op.o $(OBJ)/dirac_soc.o $(OBJ)/dirac_soc2.o $(OBJ)/changerep.o $(OBJ)/writehoststructure.o $(OBJ)/readimppot.o $(OBJ)/tmatimp_newsolver.o $(OBJ)/greenimp.o

OBJECTS1a= $(OBJ)/main1a.o $(OBJ)/calctmat.o $(OBJ)/initldau.o $(OBJ)/phicalc.o $(OBJ)/tmat_newsolver.o $(OBJ)/cheb.o $(OBJ)/interpolate_poten.o $(OBJ)/vllmat.o $(OBJ)/rotatespinframe.o $(OBJ)/spinorbit_ham.o $(OBJ)/spin_orbit.o $(OBJ)/spin_orbit_compl.o $(OBJ)/vllmatsra.o $(OBJ)/beshank.o $(OBJ)/rllsllsourceterms.o $(OBJ)/calcsph.o $(OBJ)/rllsll.o $(OBJ)/chebint.o $(OBJ)/csum.o $(OBJ)/gfree.o $(OBJ)/gll95.o $(OBJ)/grefsy.o $(OBJ)/ioben.o $(OBJ)/mapblock.o $(OBJ)/cradwf.o $(OBJ)/csinwd.o $(OBJ)/csout.o $(OBJ)/hankel.o $(OBJ)/irwsol.o $(OBJ)/regns.o $(OBJ)/tbref.o $(OBJ)/wfint.o $(OBJ)/wfint0.o $(OBJ)/wftsca.o $(OBJ)/pnstmat.o $(OBJ)/vllns.o $(OBJ)/simpk.o $(OBJ)/soutk.o $(OBJ)/ssum.o $(OBJ)/rll_local_solutions.o $(OBJ)/sll_local_solutions.o $(OBJ)/rll_global_solutions.o $(OBJ)/sll_global_solutions.o

OBJECTS1b= $(OBJ)/cmatmul.o $(OBJ)/main1b.o $(OBJ)/kloopz1.o  $(OBJ)/mssinit.o $(OBJ)/kkrmat01.o $(OBJ)/decimate.o $(OBJ)/rotgll.o $(OBJ)/setfactl.o $(OBJ)/gijdmat.o $(OBJ)/symetrmat.o $(OBJ)/projtau.o $(OBJ)/cpamillsx.o $(OBJ)/rotate.o $(OBJ)/getdmat.o $(OBJ)/bofm.o $(OBJ)/btom.o $(OBJ)/dlke0.o $(OBJ)/dlke1.o $(OBJ)/inversion.o $(OBJ)/invslab.o $(OBJ)/invsupercell.o $(OBJ)/surfgf.o

OBJECTS1c= $(OBJ)/main1c.o $(OBJ)/rhoval.o $(OBJ)/rhoval0.o $(OBJ)/rhocore.o $(OBJ)/wrldos.o $(OBJ)/wrmoms.o $(OBJ)/wrmomssoc.o $(OBJ)/densitymat.o $(OBJ)/wrldaupot.o $(OBJ)/overlap.o $(OBJ)/wmatldau.o $(OBJ)/mixldau.o $(OBJ)/rclm.o $(OBJ)/rhovalnew.o $(OBJ)/rhooutnew.o $(OBJ)/intcheb_cell.o $(OBJ)/cheb2oldgrid.o $(OBJ)/orbitalmoment.o $(OBJ)/wmatldausoc.o $(OBJ)/corel.o $(OBJ)/csimpk.o $(OBJ)/intcor.o $(OBJ)/intin.o $(OBJ)/intout.o $(OBJ)/irwns.o $(OBJ)/pnsqns.o $(OBJ)/rhoin.o $(OBJ)/rholm.o $(OBJ)/rhons.o $(OBJ)/rhoout.o $(OBJ)/renorm_lly.o $(OBJ)/simp3.o $(OBJ)/drvcore.o $(OBJ)/drvrho.o $(OBJ)/scfchrdns.o $(OBJ)/calcgf.o $(OBJ)/rintsimp.o $(OBJ)/core.o $(OBJ)/coredir.o $(OBJ)/hffcore.o $(OBJ)/corehff.o $(OBJ)/coreerr.o $(OBJ)/rint4pts.o $(OBJ)/calcmvec.o $(OBJ)/amemagvec.o $(OBJ)/rmatstr.o $(OBJ)/mvecglobal.o

OBJECTS2 = $(OBJ)/main2.o $(OBJ)/vmadelblk.o $(OBJ)/vinterface.o $(OBJ)/etotb1.o $(OBJ)/rhototb.o $(OBJ)/brydbm.o $(OBJ)/brysh1.o $(OBJ)/brysh2.o $(OBJ)/brysh3.o $(OBJ)/convol.o $(OBJ)/ecoub.o $(OBJ)/epotinb.o $(OBJ)/espcb.o $(OBJ)/force.o $(OBJ)/forceh.o $(OBJ)/forcxc.o $(OBJ)/mixstr.o $(OBJ)/mtzero.o $(OBJ)/reg2.o $(OBJ)/rhosymm.o $(OBJ)/rites.o $(OBJ)/sinwk.o $(OBJ)/vintras.o $(OBJ)/corlsd.o $(OBJ)/cpw91.o $(OBJ)/cylm02.o $(OBJ)/exch91.o $(OBJ)/gcor91.o $(OBJ)/gradr.o $(OBJ)/gradrl.o $(OBJ)/gxcpt.o $(OBJ)/mkxcpe.o $(OBJ)/spher.o $(OBJ)/sphere_gga.o $(OBJ)/sphere_nogga.o $(OBJ)/trarea.o $(OBJ)/vosko.o $(OBJ)/vxcdrv.o $(OBJ)/vxcgga.o $(OBJ)/vxclm.o $(OBJ)/vxcspo.o $(OBJ)/lebedev.o $(OBJ)/mkxcpe2.o $(OBJ)/mdirnewang.o $(OBJ)/scfiterang.o $(OBJ)/writekkrflex.o

# only these conatin explicit !$omp statements:
$(OBJ)/main0.o: $(SRC)/main0.F
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/rhovalnew.o: $(SRC)/rhovalnew.F
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/strmat.o: $(SRC)/strmat.F
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/rllsll.o: $(SRC)/rllsll.F90
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/rll_local_solutions.o: $(SRC)/rll_local_solutions.F90
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/sll_local_solutions.o: $(SRC)/sll_local_solutions.F90
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/rll_global_solutions.o: $(SRC)/rll_global_solutions.F90
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/sll_global_solutions.o: $(SRC)/sll_global_solutions.F90
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/tmat_newsolver.o: $(SRC)/tmat_newsolver.F
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
$(OBJ)/gll95.o: $(SRC)/gll95.F
	$(FC) $(FFLAGS) $(OPENMPFLAGS) $(CPP) -c -o $@ $<
# these are called in the above routines and thus have to be thread safe
$(OBJ)/gamfc.o: $(SRC)/gamfc.f
	$(FC) $(FFLAGS) $(OPENMPFLAGS) -c -o $@ $<
$(OBJ)/ymy.o: $(SRC)/ymy.f
	$(FC) $(FFLAGS) $(OPENMPFLAGS) -c -o $@ $<
$(OBJ)/gfree.o: $(SRC)/gfree.f
	$(FC) $(FFLAGS) $(OPENMPFLAGS) -c -o $@ $<

# wildcards for module files
$(OBJ)/mod_%.o: $(SRC)/%.f
	$(FC) $(FFLAGS)  -c -o $@ $<
$(OBJ)/mod_%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS)  -c -o $@ $<
$(OBJ)/mod_%.o: $(SRC)/%.F
	$(FC) $(FFLAGS)  $(CPP) -c -o $@ $<
$(OBJ)/mod_%.o: $(SRC)/%.F90
	$(FC) $(FFLAGS)  $(CPP) -c -o $@ $<

# wildcards for normal .o files
$(OBJ)/%.o: $(SRC)/%.f
	$(FC) $(FFLAGS)   -c -o $@ $<
$(OBJ)/%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS)   -c -o $@ $<
$(OBJ)/%.o: $(SRC)/%.F
	$(FC) $(FFLAGS)   $(CPP) -c -o $@ $<
$(OBJ)/%.o: $(SRC)/%.F90
	$(FC) $(FFLAGS)   $(CPP) -c -o $@ $<

####################################################
# Libraries                                        #
####################################################

# standart libs
# serial:
LLIBS  = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
# mac libs:
MACLLIBS = -mkl -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000
# RWTH cluster:
RWTHLLIBS = $(FLAGS_MKL_LINKER) #-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
# JURECA:
JURECALIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread

mac:    LLIBS=$(MACLLIBS)
macmpi: LLIBS=$(MACLLIBS)
macomp: LLIBS=$(MACLLIBS)
macompprof: LLIBS=$(MACLLIBS)

rwth:          LLIBS=$(RWTHLLIBS)
rwthmpi:       LLIBS=$(RWTHLLIBS)
rwthmpidebug:  LLIBS=$(RWTHLLIBS)

jureca: LLIBS=$(JURECALIBS)

# openmp libs
LLIBSOPENMP = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
# mac:
#LLIBSOPENMP = -L${libfolder} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
MACLLIBSOPENMP = -mkl -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000 -lmkl_intel_thread -liomp5
# RWTH cluster:
RWTHLLIBSOMP = $(FLAGS_MKL_LINKER) -lmkl_intel_lp64 -lmkl_intel_thread -liomp5


mac:    LLIBSOPENMP=$(MACLLIBSOPENMP)
macmpi: LLIBSOPENMP=$(MACLLIBSOPENMP)
macomp: LLIBSOPENMP=$(MACLLIBSOPENMP)
macompprof: LLIBSOPENMP=$(MACLLIBSOPENMP)

rwth:          LLIBSOPENMP=$(RWTHLLIBSOMP)
rwthmpi:       LLIBSOPENMP=$(RWTHLLIBSOMP)
rwthmpidebug:  LLIBSOPENMP=$(RWTHLLIBSOMP)
rwthomp:       LLIBSOPENMP=$(RWTHLLIBSOMP)
rwthompdebug:  LLIBSOPENMP=$(RWTHLLIBSOMP)

openmp: LLIBS=$(LLIBSOPENMP)
hybrid: LLIBS=$(LLIBSOPENMP)
hybriddebug: LLIBS=$(LLIBSOPENMP)
ompdebug: LLIBS=$(LLIBSOPENMP)
rwthomp: LLIBS=$(RWTHLLIBSOMP)
rwthompdebug: LLIBS=$(RWTHLLIBSOMP)
rwthhybrid: LLIBS=$(RWTHLLIBSOMP)

####################################################
# C preprocessor options                           #
####################################################
# Default option
CPP= -D CPP_TIMING
# change to if flag is matched
openmp:   CPP = -D CPP_OMP
ompdebug: CPP = -D CPP_OMP
mpidebug: CPP = -D CPP_MPI -D CPP_TIMING
mpi:      CPP = -D CPP_MPI -D CPP_TIMING
macmpi:   CPP = -D CPP_MPI -D CPP_TIMING
macomp:   CPP = -D CPP_OMP -D CPP_TIMING
macompprof:   CPP = -D CPP_OMP -D CPP_TIMING
rwthmpi:  CPP = -D CPP_MPI -D CPP_TIMING
rwthmpidebug:  CPP = -D CPP_MPI -D CPP_TIMING
rwthomp:  CPP = -D CPP_OMP -D CPP_TIMING
rwthhybrid:  CPP = -D CPP_MPI -D CPP_TIMING -D CPP_HYBRID
hybrid:   CPP = -D CPP_MPI -D CPP_TIMING -D CPP_HYBRID
hybriddebug:   CPP = -D CPP_MPI -D CPP_TIMING -D CPP_HYBRID -D CPP_NOMD5
jureca:   CPP = -D CPP_MPI -D CPP_TIMING

####################################################
# Compiler flags                                   #
####################################################
# Default:
FFLAGS = -O2 -r8 -module $(OBJ) -sox -traceback
# -xHost 
# flags for ifort: optimization, default 8 bytes, trackeback to line, compiled module files in $(OBJ), save compiler version info in executable

# Change to if flag is matched:

# on rwth cluster (CLAIX and bull)
rwthomp: FFLAGS = -O2 -r8 -traceback   -module $(OBJ) -sox
rwthhybrid: FFLAGS = -O2 -r8 -traceback   -module $(OBJ) -sox

# profiling on mac
macompprof: FFLAGS = -O2 -r8 -module $(OBJ) -sox -traceback -prof-gen:srcpos

# openmp flags
OPENMPFLAGS = 
hybrid: OPENMPFLAGS = -openmp
openmp: OPENMPFLAGS = -openmp
rwthomp: OPENMPFLAGS = -qopenmp
rwthhybrid: OPENMPFLAGS = -qopenmp


# standard debug flags
DEBUGFLAGS = -O0 -r8 -CB -check all -check uninit -ftrapuv -gen-interfaces -warn all -warn notruncated_source -fpe0 -debug extended -traceback -module $(OBJ) -g -sox

debug:       FFLAGS = $(DEBUGFLAGS)
hybriddebug: FFLAGS = $(DEBUGFLAGS)
mpidebug:    FFLAGS = $(DEBUGFLAGS)
ompdebug:    FFLAGS = $(DEBUGFLAGS)
macmpidebug: FFLAGS = $(DEBUGFLAGS)
macompdebug: FFLAGS = $(DEBUGFLAGS)
rwthmpidebug: FFLAGS = $(DEBUGFLAGS)
rwthompdebug: FFLAGS = $(DEBUGFLAGS)

hybriddebug: OPENMPFLAGS = -openmp
ompdebug: OPENMPFLAGS = -openmp
rwthompdebug: OPENMPFLAGS = -openmp



####################################################
# helper                                           #
####################################################
.PHONY: default mkdirobj isincpthere copyinc clean clear serial debug openmp mpi mpidebug install getversion compile hybrid hybriddebug mac macmpi macomp macompprof rwth rwthmpi rwthomp rwthmpidebug rwthompdebug rwthhybrid jureca

getversion:
	@echo "Wrinting Version file: "; git describe;\
	echo $(COMPVER)>compver; \
	echo $(FFLAGS) $(OPENMPFLAGS)>compflag; \
	echo $(LLIBS) >complib; \
	$(SRC)/getversion.sh
mkdirobj:
	@if test ! -d $(OBJ); \
	then mkdir $(OBJ); fi
isincpthere:
	@if test ! -e inc.p; \
	then echo "inc.p not found. Copy default values from inc.p_0"; cp inc.p_0 inc.p; \
	fi;
copyinc: isincpthere
	cp -p inc.p   $(SRC)/inc.p; \
	cp -p sprkkr_rmesh.dim   $(SRC)/sprkkr_rmesh.dim
clean:
	rm -rf *.mod *.o OBJ/ *__genmod.*
clear:
	rm -rf *.mod *.o OBJ/ *__genmod.*


# define three steps of compilation

# first compile modules (has to be done in serial)
compile0: mkdirobj getversion copyinc $(MODULES)

# then compile everything exept main_all (this can be done in parallel if -j is given)
compile1:  $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2)

# finally compile main_all and do linking in serial
compile2:  compile0 compile1 $(OBJ)/main_all.o
	  $(FC) -o $(FILENAME0) $(MODULES) $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)

# set compile steps 0,1,2 together:
compile: compile0 compile1 compile2

###################################################################




serial:   compile
debug:    compile
openmp:   compile
ompdebug: compile
mpi:      compile
hybrid:   compile
hybriddebug:   compile
mpidebug: compile
rwthmpidebug: compile

install: clean mpi

mac:    compile
macmpi: compile
macomp: compile
macompprof: compile
rwth:    compile
rwthmpi: compile
rwthomp: compile
rwthmpidebug: compile
rwthompdebug: compile
rwthhybrid:   compile
jureca: compile
