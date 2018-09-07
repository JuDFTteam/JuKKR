cmake_minimum_required(VERSION 2.6.4)

project(JuKKR)
enable_language(Fortran)

#call cmake like this:
#
# mkdir build
# cd build
# FC=<compiler-you-want-to-use> cmake -D<options> ..
#
# where with FC you coose the compiler, e.g. 
#  FC=gfortran, FC=mpif90 (gfortran with MPI), FC=ifort, FC=mpiifort
#
# List of -D<options>=default:
#  -DENABLE_MPI=ON
#  -DENABLE_OMP=OFF
#  -DENABLE_COV=OFF
#  -DCMAKE_BUILD_TYPE=Release
#  -DENABLE_BdG=OFF

# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# here the executables are written to
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ../ )

# we define some compiler flags
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    message("Intel compiler found")
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -g -traceback")
    if(ENABLE_OMP)
        set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -mkl=parallel")
    else()
        set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -mkl=sequential")
    endif()
    #set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -r8 -CB -check all -check uninit -ftrapuv -gen-interfaces -warn all -warn notruncated_source -fpe0 -debug extended -traceback -module $(OBJ) -g -sox")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -r8 -CB -check all -check uninit -gen-interfaces -warn all -warn notruncated_source -debug extended -traceback -module $(OBJ) -g -sox")
    if(ENABLE_COV)
        set(CMAKE_Fortran_FLAGS     "${CMAKE_Fortran_FLAGS} -prof-gen:srcpos")
    endif()
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -r8 -ip -xHOST")
    set(compver "intel")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    message("GNU compiler found")
    set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -lblas -llapack -ffree-line-length-none -fdefault-real-8 -fimplicit-none")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3 -fcheck=all -fbacktrace -Wall -Wextra")
    set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -march=native -fbacktrace")
    if(ENABLE_COV)
        set(CMAKE_Fortran_FLAGS     "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage -fPIC")
    endif()
    # disable MD5 sums for gnu compiler
    add_definitions(-DCPP_NOMD5)
    set(compver "gnu")
endif()

if(CMAKE_BUILD_TYPE MATCHES Debug)
    set(compver "${compver}-debug")
endif()

# common preprocessing options for all compilers
add_definitions(-DCPP_TIMING)

# default is Release build type (use '-DCMAKE_BUILD_TYPE=Debug' for Debug mode)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    message("Release build")
endif()

# MPI option (enable with '-DENABLE_MPI=ON')
# and OpenMP option (enable with '-DENABLE_OMP=ON')
# hybrid version with -DENABLE_MPI=ON and -DENABLE_OMP=ON
option(ENABLE_MPI "Enable MPI parallelization" ON)
option(ENABLE_OMP "Enable OpenMP parallelization" OFF)
set(ompflags "") # default
if(ENABLE_MPI)
    message("MPI is enabled")
    add_definitions(-DCPP_MPI)
    if(ENABLE_OMP)
        message("OpenMP (hybrid MPI/OpenMP) is enabled")
        add_definitions(-DCPP_HYBRID)
        if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
            set(ompflags "-qopenmp")
        elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            set(ompflags "-fopenmp")
        endif()
        set(compflag "hybrid")
    else()
        set(compflag "mpi")
    endif()
elseif(ENABLE_OMP)
    message("OpenMP is enabled")
    add_definitions(-DCPP_OMP)
    if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        set(ompflags "-qopenmp")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(ompflags "-fopenmp")
    endif()
    set(compflag "openmp")
else()
    message("both MPI and OpenMP is disabled")
endif()

# maybe for future librbary dependecies
set(complib "")


# extract version info
# generate a Git version string
execute_process(COMMAND git describe OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

# write compiler flags
add_definitions(-DCOMPVER1="${GIT_VERSION}")
add_definitions(-DCOMPVER2="${compver}")
add_definitions(-DCOMPVER3="${compflag}")
add_definitions(-DCOMPVER4="${complib}")


# deal with OpenMP flags:
# only these files contain OpenMP statements:
set_source_files_properties(ALL_SOURCE_FILES/gll95.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
set_source_files_properties(ALL_SOURCE_FILES/main0.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
set_source_files_properties(ALL_SOURCE_FILES/rhoqtools.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
set_source_files_properties(ALL_SOURCE_FILES/strmat.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
if(ENABLE_MPI) # this works only in hybrid mode at the moment and not with OpenMP only
    set_source_files_properties(ALL_SOURCE_FILES/tmat_newsolver.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
    set_source_files_properties(ALL_SOURCE_FILES/rhovalnew.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
    set_source_files_properties(ALL_SOURCE_FILES/rllsll.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
endif()

# parallelization of global/local solutions broken, probably because local solutions are called within global solutions
#set_source_files_properties(ALL_SOURCE_FILES/rll_global_solutions.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
#set_source_files_properties(ALL_SOURCE_FILES/sll_global_solutons.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
#set_source_files_properties(ALL_SOURCE_FILES/rll_local_solutions.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
#set_source_files_properties(ALL_SOURCE_FILES/sll_local_solutons.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
# here parallelization can be optimized (calls to subroutines dlke0 etc are serialized at the moment)
set_source_files_properties(ALL_SOURCE_FILES/kkrmat01.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")
# omp statenments commented out:
#set_source_files_properties(ALL_SOURCE_FILES/rhoval.F90 PROPERTIES COMPILE_FLAGS "${ompflags}")


# the executable is built from this list of files
add_executable(
    kkr.x
    ALL_SOURCE_FILES/DataTypes.f90
    ALL_SOURCE_FILES/addvirtual14.f90
    ALL_SOURCE_FILES/amemagvec.f90
    ALL_SOURCE_FILES/bastrmat.f90
    ALL_SOURCE_FILES/beshan.f90
    ALL_SOURCE_FILES/beshank.f90
    ALL_SOURCE_FILES/bofm.f90
    ALL_SOURCE_FILES/brydbm.f90
    ALL_SOURCE_FILES/brysh1.f90
    ALL_SOURCE_FILES/brysh2.f90
    ALL_SOURCE_FILES/brysh3.f90
    ALL_SOURCE_FILES/btom.f90
    ALL_SOURCE_FILES/bzirr3d.f90
    ALL_SOURCE_FILES/bzkint0.f90
    ALL_SOURCE_FILES/bzkmesh.f90
    ALL_SOURCE_FILES/calc_rho_ll_ss.f90
    ALL_SOURCE_FILES/calc_rho_ll_ss_lmdos.f90
    ALL_SOURCE_FILES/calc_torq_ll_ss.f90
    ALL_SOURCE_FILES/calccgc.f90
    ALL_SOURCE_FILES/calcgf.f90
    ALL_SOURCE_FILES/calcmvec.f90
    ALL_SOURCE_FILES/calcrotmat.f90
    ALL_SOURCE_FILES/calcsph.f90
    ALL_SOURCE_FILES/calctmat.F90
    ALL_SOURCE_FILES/calctref.f90
    ALL_SOURCE_FILES/calrmt.f90
    ALL_SOURCE_FILES/cdjlzdz.f90
    ALL_SOURCE_FILES/cdnlzdz.f90
    ALL_SOURCE_FILES/changerep.f90
    ALL_SOURCE_FILES/cheb.f90
    ALL_SOURCE_FILES/cheb2oldgrid.f90
    ALL_SOURCE_FILES/chebint.f90
    ALL_SOURCE_FILES/checkrmat.f90
    ALL_SOURCE_FILES/cinit.f90
    ALL_SOURCE_FILES/cint4pts.f90
    ALL_SOURCE_FILES/cintabr.f90
    ALL_SOURCE_FILES/cinthff.f90
    ALL_SOURCE_FILES/cjlz.f90
    ALL_SOURCE_FILES/clsgen_tb.f90
    ALL_SOURCE_FILES/cmatmul.f90
    ALL_SOURCE_FILES/cmatstr.f90
    ALL_SOURCE_FILES/cmomsread.f90
    ALL_SOURCE_FILES/cnlz.f90
    ALL_SOURCE_FILES/constants.f90
    ALL_SOURCE_FILES/convol.f90
    ALL_SOURCE_FILES/core.f90
    ALL_SOURCE_FILES/coredir.f90
    ALL_SOURCE_FILES/coreerr.f90
    ALL_SOURCE_FILES/corehff.f90
    ALL_SOURCE_FILES/corel.f90
    ALL_SOURCE_FILES/corlsd.f90
    ALL_SOURCE_FILES/cpamillsx.f90
    ALL_SOURCE_FILES/cpw91.f90
    ALL_SOURCE_FILES/cradwf.f90
    ALL_SOURCE_FILES/create_newmesh.f90
    ALL_SOURCE_FILES/crospr.f90
    ALL_SOURCE_FILES/cross.f90
    ALL_SOURCE_FILES/crtstar.f90
    ALL_SOURCE_FILES/csimpk.f90
    ALL_SOURCE_FILES/csinwd.f90
    ALL_SOURCE_FILES/csout.f90
    ALL_SOURCE_FILES/csum.f90
    ALL_SOURCE_FILES/cylm02.f90
    ALL_SOURCE_FILES/ddet33.f90
    ALL_SOURCE_FILES/ddot1.f90
    ALL_SOURCE_FILES/decimaread.f90
    ALL_SOURCE_FILES/decimate.f90
    ALL_SOURCE_FILES/deciopt.f90
    ALL_SOURCE_FILES/decipotbas.f90
    ALL_SOURCE_FILES/decipothead.f90
    ALL_SOURCE_FILES/decitmat.f90
    ALL_SOURCE_FILES/decitset.f90
    ALL_SOURCE_FILES/densitymat.f90
    ALL_SOURCE_FILES/dinv33.f90
    ALL_SOURCE_FILES/dirac_bi.f90
    ALL_SOURCE_FILES/dirac_bs.f90
    ALL_SOURCE_FILES/dirac_op.f90
    ALL_SOURCE_FILES/dirac_soc.f90
    ALL_SOURCE_FILES/dirac_soc2.f90
    ALL_SOURCE_FILES/dirbslag.f90
    ALL_SOURCE_FILES/dirbsmid.f90
    ALL_SOURCE_FILES/dirbsrad.f90
    ALL_SOURCE_FILES/dirbsrze.f90
    ALL_SOURCE_FILES/dirbsstp.f90
    ALL_SOURCE_FILES/dlke0.f90
    ALL_SOURCE_FILES/dlke1.f90
    ALL_SOURCE_FILES/dmpy.f90
    ALL_SOURCE_FILES/drvbastrans.f90
    ALL_SOURCE_FILES/drvcore.f90
    ALL_SOURCE_FILES/drvreltmat.f90
    ALL_SOURCE_FILES/drvrho.f90
    ALL_SOURCE_FILES/dscal1.f90
    ALL_SOURCE_FILES/dsort.f90
    ALL_SOURCE_FILES/dswap1.f90
    ALL_SOURCE_FILES/ecoub.f90
    ALL_SOURCE_FILES/emesht.f90
    ALL_SOURCE_FILES/epathtb.f90
    ALL_SOURCE_FILES/epotinb.f90
    ALL_SOURCE_FILES/erfcex.f90
    ALL_SOURCE_FILES/errmsg.f90
    ALL_SOURCE_FILES/errortrap.f90
    ALL_SOURCE_FILES/espcb.f90
    ALL_SOURCE_FILES/etotb1.f90
    ALL_SOURCE_FILES/ewald2d.f90
    ALL_SOURCE_FILES/exch91.f90
    ALL_SOURCE_FILES/findgroup.f90
    ALL_SOURCE_FILES/force.f90
    ALL_SOURCE_FILES/forceh.f90
    ALL_SOURCE_FILES/forcxc.f90
    ALL_SOURCE_FILES/fplaneg.f90
    ALL_SOURCE_FILES/fplaner.f90
    ALL_SOURCE_FILES/gamfc.f90
    ALL_SOURCE_FILES/gaufd.f90
    ALL_SOURCE_FILES/gauleg.f90
    ALL_SOURCE_FILES/gaunt.f90
    ALL_SOURCE_FILES/gaunt2.f90
    ALL_SOURCE_FILES/gcor91.f90
    ALL_SOURCE_FILES/generalpot.f90
    ALL_SOURCE_FILES/getbr3.f90
    ALL_SOURCE_FILES/getclusnxyz.f90
    ALL_SOURCE_FILES/getdmat.f90
    ALL_SOURCE_FILES/getscratch.f90
    ALL_SOURCE_FILES/gfmask.f90
    ALL_SOURCE_FILES/gfree.f90
    ALL_SOURCE_FILES/gfshells.f90
    ALL_SOURCE_FILES/gijcond.f90
    ALL_SOURCE_FILES/gijdmat.f90
    ALL_SOURCE_FILES/gijxcpl.f90
    ALL_SOURCE_FILES/gll95.F90
    ALL_SOURCE_FILES/global_variables.F90
    ALL_SOURCE_FILES/godfrin.F90
    ALL_SOURCE_FILES/gradr.f90
    ALL_SOURCE_FILES/gradr_d0.f90
    ALL_SOURCE_FILES/gradrl.f90
    ALL_SOURCE_FILES/greenimp.F90
    ALL_SOURCE_FILES/grefsy.f90
    ALL_SOURCE_FILES/grule.f90
    ALL_SOURCE_FILES/gxcpt.f90
    ALL_SOURCE_FILES/hankel.f90
    ALL_SOURCE_FILES/hffcore.f90
    ALL_SOURCE_FILES/idreals.f90
    ALL_SOURCE_FILES/ikapmue.f90
    ALL_SOURCE_FILES/ikmlin.f90
    ALL_SOURCE_FILES/impcheck.f90
    ALL_SOURCE_FILES/impcoefs.f90
    ALL_SOURCE_FILES/initabjij.f90
    ALL_SOURCE_FILES/initldau.f90
    ALL_SOURCE_FILES/intcheb_cell.f90
    ALL_SOURCE_FILES/intcor.f90
    ALL_SOURCE_FILES/interpolate_poten.f90
    ALL_SOURCE_FILES/interpolspline.f90
    ALL_SOURCE_FILES/intin.f90
    ALL_SOURCE_FILES/intout.f90
    ALL_SOURCE_FILES/inversion.f90
    ALL_SOURCE_FILES/invslab.f90
    ALL_SOURCE_FILES/invsupercell.f90
    ALL_SOURCE_FILES/ioben.f90
    ALL_SOURCE_FILES/ioinput.f90
    ALL_SOURCE_FILES/irwns.f90
    ALL_SOURCE_FILES/irwsol.f90
    ALL_SOURCE_FILES/jijhelp.f90
    ALL_SOURCE_FILES/kkrmat01.F90
    ALL_SOURCE_FILES/kloopz1.F90
    ALL_SOURCE_FILES/lattice2d.f90
    ALL_SOURCE_FILES/lattice3d.f90
    ALL_SOURCE_FILES/lattix99.f90
    ALL_SOURCE_FILES/latvec.f90
    ALL_SOURCE_FILES/lebedev.f90
    ALL_SOURCE_FILES/length.f90
    ALL_SOURCE_FILES/lngstring.f90
    ALL_SOURCE_FILES/madelcoef.f90
    ALL_SOURCE_FILES/madelgaunt.f90
    ALL_SOURCE_FILES/madelout.f90
    ALL_SOURCE_FILES/madelung2d.f90
    ALL_SOURCE_FILES/madelung3d.f90
    ALL_SOURCE_FILES/main0.F90
    ALL_SOURCE_FILES/main1a.F90
    ALL_SOURCE_FILES/main1b.F90
    ALL_SOURCE_FILES/main1c.F90
    ALL_SOURCE_FILES/main2.F90
    ALL_SOURCE_FILES/main_all.F90
    ALL_SOURCE_FILES/mapblock.f90
    ALL_SOURCE_FILES/md5sums.F90
    ALL_SOURCE_FILES/mdirnewang.f90
    ALL_SOURCE_FILES/memoryhandling.f90
    ALL_SOURCE_FILES/mixldau.f90
    ALL_SOURCE_FILES/mixstr.f90
    ALL_SOURCE_FILES/mkxcpe.f90
    ALL_SOURCE_FILES/mkxcpe2.f90
    ALL_SOURCE_FILES/mssinit.f90
    ALL_SOURCE_FILES/mtzero.f90
    ALL_SOURCE_FILES/mvecglobal.f90
    ALL_SOURCE_FILES/mympi.F90
    ALL_SOURCE_FILES/normcoeff_SO.F90
    ALL_SOURCE_FILES/normcoeff_SO_spinflux.F90
    ALL_SOURCE_FILES/normcoeff_SO_torq.F90
    ALL_SOURCE_FILES/nrmliz.f90
    ALL_SOURCE_FILES/operators_for_FScode.F90
    ALL_SOURCE_FILES/opt.f90
    ALL_SOURCE_FILES/orbitalmoment.f90
    ALL_SOURCE_FILES/outpothost.f90
    ALL_SOURCE_FILES/outtmathost.f90
    ALL_SOURCE_FILES/overlap.f90
    ALL_SOURCE_FILES/phicalc.f90
    ALL_SOURCE_FILES/pnsqns.f90
    ALL_SOURCE_FILES/pnstmat.f90
    ALL_SOURCE_FILES/pointgrp.f90
    ALL_SOURCE_FILES/potcut.f90
    ALL_SOURCE_FILES/profiling.f90
    ALL_SOURCE_FILES/projtau.f90
    ALL_SOURCE_FILES/rclm.f90
    ALL_SOURCE_FILES/rcstop.f90
    ALL_SOURCE_FILES/readimppot.f90
    ALL_SOURCE_FILES/readldaupot.f90
    ALL_SOURCE_FILES/reg2.f90
    ALL_SOURCE_FILES/regns.f90
    ALL_SOURCE_FILES/regsol.f90
    ALL_SOURCE_FILES/relpotcvt.f90
    ALL_SOURCE_FILES/renorm_lly.f90
    ALL_SOURCE_FILES/rhocore.f90
    ALL_SOURCE_FILES/rhoin.f90
    ALL_SOURCE_FILES/rholm.f90
    ALL_SOURCE_FILES/rhons.F90
    ALL_SOURCE_FILES/rhoout.f90
    ALL_SOURCE_FILES/rhooutnew.f90
    ALL_SOURCE_FILES/rhoqtools.F90
    ALL_SOURCE_FILES/rhosymm.f90
    ALL_SOURCE_FILES/rhototb.f90
    ALL_SOURCE_FILES/rhoval.F90
    ALL_SOURCE_FILES/rhoval0.f90
    ALL_SOURCE_FILES/rhovalnew.F90
    ALL_SOURCE_FILES/rinit.f90
    ALL_SOURCE_FILES/rinput13.F90
    ALL_SOURCE_FILES/rint4pts.f90
    ALL_SOURCE_FILES/rintsimp.f90
    ALL_SOURCE_FILES/rinvgj.f90
    ALL_SOURCE_FILES/rites.f90
    ALL_SOURCE_FILES/rll_global_solutions.F90
    ALL_SOURCE_FILES/rll_local_solutions.F90
    ALL_SOURCE_FILES/rllsll.F90
    ALL_SOURCE_FILES/rllsllsourceterms.f90
    ALL_SOURCE_FILES/rmatstr.f90
    ALL_SOURCE_FILES/rnuctab.f90
    ALL_SOURCE_FILES/rotate.f90
    ALL_SOURCE_FILES/rotatespinframe.f90
    ALL_SOURCE_FILES/rotgll.F90
    ALL_SOURCE_FILES/rotmat.f90
    ALL_SOURCE_FILES/rrgen.f90
    ALL_SOURCE_FILES/save_wavefun.F90
    ALL_SOURCE_FILES/scalevec.f90
    ALL_SOURCE_FILES/scalpr.f90
    ALL_SOURCE_FILES/scfchrdns.f90
    ALL_SOURCE_FILES/scfiterang.f90
    ALL_SOURCE_FILES/setfactl.f90
    ALL_SOURCE_FILES/setgijtab.f90
    ALL_SOURCE_FILES/shape_corr.f90
    ALL_SOURCE_FILES/shellgen2k.f90
    ALL_SOURCE_FILES/simp3.f90
    ALL_SOURCE_FILES/simpk.f90
    ALL_SOURCE_FILES/sinwk.f90
    ALL_SOURCE_FILES/sll_global_solutions.F90
    ALL_SOURCE_FILES/sll_local_solutions.F90
    ALL_SOURCE_FILES/sname.f90
    ALL_SOURCE_FILES/soutk.f90
    ALL_SOURCE_FILES/spatpr.f90
    ALL_SOURCE_FILES/spher.f90
    ALL_SOURCE_FILES/sphere_gga.f90
    ALL_SOURCE_FILES/sphere_nogga.f90
    ALL_SOURCE_FILES/spin_orbit.f90
    ALL_SOURCE_FILES/spin_orbit_compl.f90
    ALL_SOURCE_FILES/spinorbit_ham.f90
    ALL_SOURCE_FILES/spline_real.f90
    ALL_SOURCE_FILES/splint_real.f90
    ALL_SOURCE_FILES/ssite.f90
    ALL_SOURCE_FILES/ssum.f90
    ALL_SOURCE_FILES/startb1.f90
    ALL_SOURCE_FILES/startldau.f90
    ALL_SOURCE_FILES/strmat.F90
    ALL_SOURCE_FILES/strsmat.f90
    ALL_SOURCE_FILES/sumupint.f90
    ALL_SOURCE_FILES/surfgf.f90
    ALL_SOURCE_FILES/symetrmat.f90
    ALL_SOURCE_FILES/symlat.f90
    ALL_SOURCE_FILES/symtaumat.f90
    ALL_SOURCE_FILES/taustruc.f90
    ALL_SOURCE_FILES/tbref.F90
    ALL_SOURCE_FILES/tbxccpljij.F90
    ALL_SOURCE_FILES/tbxccpljijdij.F90
    ALL_SOURCE_FILES/test.f90
    ALL_SOURCE_FILES/testdim.f90
    ALL_SOURCE_FILES/timing.F90
    ALL_SOURCE_FILES/tmat_newsolver.F90
    ALL_SOURCE_FILES/tmatimp_newsolver.F90
    ALL_SOURCE_FILES/trarea.f90
    ALL_SOURCE_FILES/types.F90
    ALL_SOURCE_FILES/vadd.f90
    ALL_SOURCE_FILES/veq.f90
    ALL_SOURCE_FILES/version.F90
    ALL_SOURCE_FILES/version_info.F90
    ALL_SOURCE_FILES/vinterface.f90
    ALL_SOURCE_FILES/vintras.f90
    ALL_SOURCE_FILES/vllmat.f90
    ALL_SOURCE_FILES/vllmatsra.f90
    ALL_SOURCE_FILES/vllns.f90
    ALL_SOURCE_FILES/vmadelblk.f90
    ALL_SOURCE_FILES/vmul.f90
    ALL_SOURCE_FILES/vosko.f90
    ALL_SOURCE_FILES/vxcdrv.f90
    ALL_SOURCE_FILES/vxcgga.f90
    ALL_SOURCE_FILES/vxclm.f90
    ALL_SOURCE_FILES/vxcspo.f90
    ALL_SOURCE_FILES/wfint.f90
    ALL_SOURCE_FILES/wfmesh.f90
    ALL_SOURCE_FILES/wftsca.f90
    ALL_SOURCE_FILES/wmatldau.f90
    ALL_SOURCE_FILES/wmatldausoc.f90
    ALL_SOURCE_FILES/write_tbkkr_files.f90
    ALL_SOURCE_FILES/writehoststructure.f90
    ALL_SOURCE_FILES/writekkrflex.f90
    ALL_SOURCE_FILES/wrldaupot.f90
    ALL_SOURCE_FILES/wrldos.f90
    ALL_SOURCE_FILES/wrmoms.f90
    ALL_SOURCE_FILES/wunfiles.F90
    ALL_SOURCE_FILES/ylag.f90
    ALL_SOURCE_FILES/ymy.f90
    )


if(ENABLE_BdG)
    # for wf testing (BdG)
    add_executable(
        kkr_BdG_radial.x
        ALL_SOURCE_FILES/main_tmat_runner.F90
        ALL_SOURCE_FILES/DataTypes.f90
        ALL_SOURCE_FILES/version.F90
        ALL_SOURCE_FILES/version_info.F90 
        ALL_SOURCE_FILES/constants.f90 
        ALL_SOURCE_FILES/mympi.F90 
        ALL_SOURCE_FILES/global_variables.F90 
        ALL_SOURCE_FILES/types.F90 
        ALL_SOURCE_FILES/timing.F90 
        ALL_SOURCE_FILES/profiling.f90 
        ALL_SOURCE_FILES/save_wavefun.F90 
        ALL_SOURCE_FILES/wunfiles.F90 
        ALL_SOURCE_FILES/md5sums.F90 
        ALL_SOURCE_FILES/jijhelp.f90 
        ALL_SOURCE_FILES/initabjij.f90 
        ALL_SOURCE_FILES/tbxccpljij.F90 
        ALL_SOURCE_FILES/tbxccpljijdij.F90 
        ALL_SOURCE_FILES/create_newmesh.f90 
        ALL_SOURCE_FILES/cheb.f90 
        ALL_SOURCE_FILES/rhoqtools.F90 
        ALL_SOURCE_FILES/godfrin.F90 
        ALL_SOURCE_FILES/memoryhandling.f90 
        ALL_SOURCE_FILES/rinput13.F90
        ALL_SOURCE_FILES/getbr3.f90 
        ALL_SOURCE_FILES/write_tbkkr_files.f90 
        ALL_SOURCE_FILES/gaunt.f90 
        ALL_SOURCE_FILES/symtaumat.f90 
        ALL_SOURCE_FILES/taustruc.f90 
        ALL_SOURCE_FILES/bzkint0.f90 
        ALL_SOURCE_FILES/findgroup.f90 
        ALL_SOURCE_FILES/calcrotmat.f90 
        ALL_SOURCE_FILES/errortrap.f90 
        ALL_SOURCE_FILES/bzkmesh.f90 
        ALL_SOURCE_FILES/bzirr3d.f90 
        ALL_SOURCE_FILES/madelung3d.f90 
        ALL_SOURCE_FILES/madelung2d.f90 
        ALL_SOURCE_FILES/madelout.f90 
        ALL_SOURCE_FILES/idreals.f90 
        ALL_SOURCE_FILES/deciopt.f90 
        ALL_SOURCE_FILES/impcheck.f90 
        ALL_SOURCE_FILES/impcoefs.f90 
        ALL_SOURCE_FILES/setgijtab.f90 
        ALL_SOURCE_FILES/gijxcpl.f90 
        ALL_SOURCE_FILES/gijcond.f90 
        ALL_SOURCE_FILES/getclusnxyz.f90 
        ALL_SOURCE_FILES/decitset.f90    
        ALL_SOURCE_FILES/decitmat.f90 
        ALL_SOURCE_FILES/outpothost.f90  
        ALL_SOURCE_FILES/decipothead.f90 
        ALL_SOURCE_FILES/getscratch.f90 
        ALL_SOURCE_FILES/decipotbas.f90  
        ALL_SOURCE_FILES/outtmathost.f90 
        ALL_SOURCE_FILES/lngstring.f90 
        ALL_SOURCE_FILES/startldau.f90   
        ALL_SOURCE_FILES/readldaupot.f90 
        ALL_SOURCE_FILES/epathtb.f90 
        ALL_SOURCE_FILES/addvirtual14.f90   
        ALL_SOURCE_FILES/interpolspline.f90 
        ALL_SOURCE_FILES/spline_real.f90 
        ALL_SOURCE_FILES/splint_real.f90 
        ALL_SOURCE_FILES/calrmt.f90     
        ALL_SOURCE_FILES/clsgen_tb.f90   
        ALL_SOURCE_FILES/clustcomp.f90 
        ALL_SOURCE_FILES/cmomsread.f90 
        ALL_SOURCE_FILES/crospr.f90 
        ALL_SOURCE_FILES/cross.f90 
        ALL_SOURCE_FILES/crtstar.f90 
        ALL_SOURCE_FILES/ddet33.f90 
        ALL_SOURCE_FILES/ddot1.f90 
        ALL_SOURCE_FILES/dscal1.f90 
        ALL_SOURCE_FILES/dsort.f90 
        ALL_SOURCE_FILES/dswap1.f90 
        ALL_SOURCE_FILES/emesht.f90 
        ALL_SOURCE_FILES/erfcex.f90 
        ALL_SOURCE_FILES/errmsg.f90 
        ALL_SOURCE_FILES/gamfc.f90 
        ALL_SOURCE_FILES/gaufd.f90  
        ALL_SOURCE_FILES/gauleg.f90 
        ALL_SOURCE_FILES/gaunt2.f90 
        ALL_SOURCE_FILES/gfmask.f90 
        ALL_SOURCE_FILES/gfshells.f90 
        ALL_SOURCE_FILES/grule.f90 
        ALL_SOURCE_FILES/ioinput.f90 
        ALL_SOURCE_FILES/lattice2d.f90 
        ALL_SOURCE_FILES/lattice3d.f90 
        ALL_SOURCE_FILES/lattix99.f90 
        ALL_SOURCE_FILES/latvec.f90 
        ALL_SOURCE_FILES/length.f90 
        ALL_SOURCE_FILES/nrmliz.f90 
        ALL_SOURCE_FILES/opt.f90 
        ALL_SOURCE_FILES/pointgrp.f90 
        ALL_SOURCE_FILES/potcut.f90 
        ALL_SOURCE_FILES/rcstop.f90 
        ALL_SOURCE_FILES/rinit.f90 
        ALL_SOURCE_FILES/decimaread.f90 
        ALL_SOURCE_FILES/rrgen.f90 
        ALL_SOURCE_FILES/scalevec.f90 
        ALL_SOURCE_FILES/cinit.f90 
        ALL_SOURCE_FILES/scalpr.f90 
        ALL_SOURCE_FILES/shellgen2k.f90 
        ALL_SOURCE_FILES/sname.f90 
        ALL_SOURCE_FILES/spatpr.f90 
        ALL_SOURCE_FILES/strmat.F90 
        ALL_SOURCE_FILES/test.f90 
        ALL_SOURCE_FILES/generalpot.f90 
        ALL_SOURCE_FILES/ymy.f90 
        ALL_SOURCE_FILES/veq.f90 
        ALL_SOURCE_FILES/vadd.f90 
        ALL_SOURCE_FILES/vmul.f90 
        ALL_SOURCE_FILES/testdim.f90 
        ALL_SOURCE_FILES/startb1.f90 
        ALL_SOURCE_FILES/madelcoef.f90 
        ALL_SOURCE_FILES/madelgaunt.f90 
        ALL_SOURCE_FILES/ewald2d.f90 
        ALL_SOURCE_FILES/fplaneg.f90 
        ALL_SOURCE_FILES/fplaner.f90 
        ALL_SOURCE_FILES/shape_corr.f90 
        ALL_SOURCE_FILES/wfmesh.f90 
        ALL_SOURCE_FILES/regsol.f90 
        ALL_SOURCE_FILES/beshan.f90 
        ALL_SOURCE_FILES/drvbastrans.f90 
        ALL_SOURCE_FILES/relpotcvt.f90 
        ALL_SOURCE_FILES/drvreltmat.f90 
        ALL_SOURCE_FILES/calccgc.f90 
        ALL_SOURCE_FILES/ikmlin.f90 
        ALL_SOURCE_FILES/strsmat.f90 
        ALL_SOURCE_FILES/bastrmat.f90 
        ALL_SOURCE_FILES/ikapmue.f90 
        ALL_SOURCE_FILES/cmatstr.f90 
        ALL_SOURCE_FILES/checkrmat.f90 
        ALL_SOURCE_FILES/rinvgj.f90 
        ALL_SOURCE_FILES/ssite.f90 
        ALL_SOURCE_FILES/cintabr.f90 
        ALL_SOURCE_FILES/sumupint.f90 
        ALL_SOURCE_FILES/dirac_bs.f90 
        ALL_SOURCE_FILES/dirbslag.f90 
        ALL_SOURCE_FILES/dirbsmid.f90 
        ALL_SOURCE_FILES/dirbsrad.f90 
        ALL_SOURCE_FILES/dirbsrze.f90 
        ALL_SOURCE_FILES/dirbsstp.f90 
        ALL_SOURCE_FILES/cjlz.f90 
        ALL_SOURCE_FILES/cnlz.f90 
        ALL_SOURCE_FILES/cdjlzdz.f90 
        ALL_SOURCE_FILES/cdnlzdz.f90 
        ALL_SOURCE_FILES/cint4pts.f90 
        ALL_SOURCE_FILES/cinthff.f90 
        ALL_SOURCE_FILES/ylag.f90 
        ALL_SOURCE_FILES/rnuctab.f90 
        ALL_SOURCE_FILES/dirac_op.f90 
        ALL_SOURCE_FILES/dirac_soc.f90 
        ALL_SOURCE_FILES/dirac_soc2.f90 
        ALL_SOURCE_FILES/changerep.f90 
        ALL_SOURCE_FILES/writehoststructure.f90 
        ALL_SOURCE_FILES/readimppot.f90 
        ALL_SOURCE_FILES/tmatimp_newsolver.F90 
        ALL_SOURCE_FILES/greenimp.F90 
        ALL_SOURCE_FILES/calctref.f90 
        ALL_SOURCE_FILES/main0.F90
        ALL_SOURCE_FILES/tmat_newsolver.F90 
        ALL_SOURCE_FILES/main1a_dummy_BdG.f90 
        ALL_SOURCE_FILES/calctmat.F90 
        ALL_SOURCE_FILES/initldau.f90 
        ALL_SOURCE_FILES/phicalc.f90 
        ALL_SOURCE_FILES/interpolate_poten.f90 
        ALL_SOURCE_FILES/vllmat.f90 
        ALL_SOURCE_FILES/rotatespinframe.f90 
        ALL_SOURCE_FILES/spinorbit_ham.f90 
        ALL_SOURCE_FILES/spin_orbit.f90 
        ALL_SOURCE_FILES/spin_orbit_compl.f90 
        ALL_SOURCE_FILES/vllmatsra.f90 
        ALL_SOURCE_FILES/beshank.f90 
        ALL_SOURCE_FILES/rllsllsourceterms.f90 
        ALL_SOURCE_FILES/calcsph.f90 
        ALL_SOURCE_FILES/rllsll.F90 
        ALL_SOURCE_FILES/chebint.f90 
        ALL_SOURCE_FILES/csum.f90 
        ALL_SOURCE_FILES/gfree.f90 
        ALL_SOURCE_FILES/gll95.F90 
        ALL_SOURCE_FILES/grefsy.f90 
        ALL_SOURCE_FILES/ioben.f90 
        ALL_SOURCE_FILES/mapblock.f90 
        ALL_SOURCE_FILES/cradwf.f90 
        ALL_SOURCE_FILES/csinwd.f90 
        ALL_SOURCE_FILES/csout.f90 
        ALL_SOURCE_FILES/hankel.f90 
        ALL_SOURCE_FILES/irwsol.f90 
        ALL_SOURCE_FILES/regns.f90 
        ALL_SOURCE_FILES/trarea.f90
        ALL_SOURCE_FILES/tbref.F90 
        ALL_SOURCE_FILES/wfint.f90 
        ALL_SOURCE_FILES/wfint0.f90 
        ALL_SOURCE_FILES/wftsca.f90 
        ALL_SOURCE_FILES/pnstmat.f90 
        ALL_SOURCE_FILES/vllns.f90 
        ALL_SOURCE_FILES/simpk.f90 
        ALL_SOURCE_FILES/soutk.f90 
        ALL_SOURCE_FILES/ssum.f90 
        ALL_SOURCE_FILES/rll_local_solutions.F90 
        ALL_SOURCE_FILES/sll_local_solutions.F90 
        ALL_SOURCE_FILES/rll_global_solutions.F90 
        ALL_SOURCE_FILES/sll_global_solutions.F90 
        ALL_SOURCE_FILES/intcheb_cell.f90 
        ALL_SOURCE_FILES/cmatmul.f90 
        ALL_SOURCE_FILES/convol.f90
        )
endif()
