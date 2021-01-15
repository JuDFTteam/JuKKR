!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_scattering

  use mod_ioformat, only: MODE_INT, MODE_VIS
  implicit none

  private
  public :: impcls_TYPE, read_scattmat, calc_scattering_inputcard



  type :: sca_TYPE
    integer :: N1 = 17
    integer :: lscatfixk=-1
    integer :: llifetime=-1
    integer :: lboltzmann=-1
    integer :: mode=-1
    integer :: nsteps=-1
    integer :: niter=-1
    integer :: roottake=-1
    integer :: savepkk=-1
    integer :: naverage=1
    integer :: gammamode=0
    double precision :: rooteps   = -1d0
    double precision :: kfix(3,2) =  1d38
    double precision :: gammaval  =  0.001837465441 ! 25 meV in Rydbergs
    double precision :: impconc   =  1d0
    integer :: subarr_inp(2) = -1
    integer :: maskint = 0
    !allocatable arrays
    integer :: N2 = 2
    double precision, allocatable :: weight_imp(:)
!   integer,          allocatable :: ispincomb(:)
!   double precision, allocatable :: nvect(:,:)
  end type sca_TYPE



  type :: impcls_TYPE
    integer          :: N1 = 3
    integer          :: nCluster = -1
    integer          :: clmso = -1
    !allocatable arrays
    integer          :: N2 = 3
    double precision,  allocatable :: RCluster(:,:)
    integer,           allocatable :: ihosttype(:)
  end type impcls_TYPE



  logical, save :: sca_read=.false.
  type(sca_type), save :: sca


contains



  subroutine calc_scattering_inputcard(inc, lattice, cluster, tgmatrx)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx

    type(impcls_type), allocatable :: impcls(:)
    double complex, allocatable :: Amat(:,:,:)

    if(.not.sca_read) then
      call read_sca()
    end if

    if(sca%lscatfixk==1 .or. sca%llifetime==1 .or. sca%lboltzmann==1) call read_scattmat(inc, impcls, Amat)
    if(sca%lscatfixk==1)  call calc_scattering_fixk(inc, lattice, cluster, tgmatrx, impcls(1), Amat(:,:,1) )
    if(sca%llifetime==1)  call calc_lifetime(inc, lattice, impcls(1), Amat(:,:,1) )
    if(sca%lboltzmann==1) call calc_boltzmann(inc, lattice, impcls, Amat)

  end subroutine calc_scattering_inputcard



  subroutine calc_boltzmann(inc, lattice, impcls, Amat)
    use type_inc,       only: inc_type
    use type_data,      only: lattice_type
    use mod_symmetries, only: symmetries_type, set_symmetries, rotate_kpoints, expand_areas, expand_spinvalues, expand_torqvalues
    use mod_read,       only: read_kpointsfile_int, read_weights, read_fermivelocity, read_spinvalue, read_torqvalue, read_torqvalue_atom,&
                            & read_spinvec_atom, read_spinflux_atom
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_iohelp,     only: open_mpifile_setview, close_mpifile, getBZvolume, file_present
    use mod_ioformat,   only: filemode_int, filename_eigvect, filename_scattmat, fmt_fn_ext, fmt_fn_atom_sub_ext, fmt_fn_sub_ext, ext_vtkxml, ext_mpiio
    use mod_ioformat,   only: filename_torq, filename_spinvec, filename_spinflux, ext_formatted
    use mod_mympi,      only: myrank, nranks, master
    use mod_mathtools,  only: pi
#ifdef CPP_TIMING
    use mod_timing,     only: timing_start, timing_stop
#endif
    use mpi
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(impcls_type),  intent(in) :: impcls(sca%naverage)
    double complex,     intent(in) :: Amat(impcls(1)%clmso,impcls(1)%clmso,sca%naverage)

    double precision :: BZVol

    !symmetry arrays
    integer :: nsym
    integer, allocatable :: isym(:)
    type(symmetries_type) :: symmetries

    !local k-point arrays
    integer :: nkpts
    double precision, allocatable :: kpoints(:,:), areas(:), weights(:), fermivel(:,:)

    !spinvalue locals
    integer :: nsqa, ndegen1
    integer, allocatable :: ispincomb(:)
    double precision, allocatable :: nvect(:,:), spinval(:,:,:), spinval1(:,:,:)

    !torque values locals
    integer                       :: itmp
    double precision, allocatable :: torqval(:,:,:), torqval_atom(:,:,:,:), spinvec_atom(:,:,:,:), spinflux_atom(:,:,:,:)
    double precision, allocatable :: arrtmp1(:,:), arrtmp2(:,:)

    !subarray locals
    integer :: myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col
    integer :: nkpt1, nkpt2, ioff1, ioff2, subarr_dim(2)
    integer :: dataarr_lb(0:nranks-1,2),  &
             & dataarr_ub(0:nranks-1,2),  &
             & dataarr_nkpt(0:nranks-1,2)
    double complex, allocatable :: rveig_dim1(:,:,:,:,:), rveig_dim2(:,:,:,:,:)
    double precision, allocatable :: Pkksub(:,:,:,:,:)

    !result locals
    double precision, allocatable :: tau(:,:,:), tau2(:,:,:), tau_avg(:,:), meanfreepath(:,:,:,:), chcond(:,:,:), spcond(:,:,:)

    !temp k-point arrays
    integer :: nkpts1, nkpts2
    double precision, allocatable :: areas1(:), weights1(:), kpoints1(:,:), fermivel1(:,:), temparr(:,:)

    !file-io
    logical :: compute_scattmat=.true.
    character(len=256) :: filename

    integer :: ierr, ikp, ii3, iat

    logical :: l_torqfile, l_torqfile_atom, l_spinvecfile_atom, l_spinfluxfile_atom

    !parameter
    double precision, parameter :: RyToinvfs = 20.67068667282055d0

    if(.not.sca_read) then
      call read_sca()
    end if
    call set_symmetries(inc, lattice, symmetries)
    BZVol = getBZvolume(lattice%recbv)

#ifdef CPP_TIMING
    call timing_start('Read in of data')
#endif

    !======================================!
    != determine whether to calculate the =!
    !=   scattering matrix or read it in  =!
    if(sca%savepkk==2)then
      write(filename,fmt_fn_ext) filename_scattmat, ext_mpiio
      if(file_present(filename)) compute_scattmat=.false.
    end if

    !=============================!
    != Read in the k-point files =!
    call read_kpointsfile_int(nkpts1, kpoints1, areas1, nsym, isym)
    call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
    call expand_areas(nsym,nkpts1,areas1,areas)
    deallocate(isym, areas1, kpoints1)

    !===================================!
    != Read in the integration weights =!
    call read_weights(nkpts1, weights1, nsym, isym)
    call expand_areas(nsym,nkpts1,weights1,weights)
    deallocate(isym, weights1)

    !===========================================!
    != Check for which observables files exist =!
    write(filename,fmt_fn_sub_ext) filename_torq, trim(filemode_int), ext_formatted
    inquire(file=filename, exist=l_torqfile)

    ! l_torqfile_atom is true only if files are found for all atoms
    l_torqfile_atom = .true.
    do iat=1,inc%natypd
      write(filename,fmt_fn_atom_sub_ext) filename_torq, iat, trim(filemode_int), ext_formatted
      inquire(file=filename, exist=l_torqfile_atom)
      if(.not.l_torqfile_atom) exit
    end do!iat=1,inc%natypd

    ! l_spinvecfile_atom is true only if files are found for all atoms
    l_spinvecfile_atom = .true.
    do iat=1,inc%natypd
      write(filename,fmt_fn_atom_sub_ext) filename_spinvec, iat, trim(filemode_int), ext_formatted
      inquire(file=filename, exist=l_spinvecfile_atom)
      if(.not.l_spinvecfile_atom) exit
    end do!iat=1,inc%natypd

    ! l_spinfluxfile_atom is true only if files are found for all atoms
    l_spinfluxfile_atom = .true.
    do iat=1,inc%natypd
      write(filename,fmt_fn_atom_sub_ext) filename_spinflux, iat, trim(filemode_int), ext_formatted
      inquire(file=filename, exist=l_spinfluxfile_atom)
      if(.not.l_spinfluxfile_atom) exit
    end do!iat=1,inc%natypd

    !==============================!
    != Read in the fermi velocity =!
    call read_fermivelocity(filemode_int, nkpts1, fermivel1, nsym, isym)
    call rotate_kpoints(symmetries%rotmat, nkpts1, fermivel1, nsym, isym, nkpts2, fermivel)
!   call project_fermivel_newaxis(nkpts,fermivel)
    if(nkpts2/=nkpts) stop 'inconsistency in number of k-points'
    deallocate(isym, fermivel1)
    if(myrank==master) then
      write(*,'(A,3ES25.16)') 'kpoints-sum:  ', sum(kpoints,  dim=2)
      write(*,'(A,3ES25.16)') 'fermivel-sum: ', sum(fermivel, dim=2)

      !compute weighted sums
      allocate(temparr(3,nkpts), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating temparr'
      do ii3=1,3
        temparr(ii3,:) = kpoints(ii3,:)*weights(:)
      end do!ii3
      write(*,'(A,3ES25.16)') 'weighted kpoints-sum:  ', sum(temparr,  dim=2)
      do ii3=1,3
        temparr(ii3,:) = fermivel(ii3,:)*weights(:)
      end do!ii3
      write(*,'(A,3ES25.16)') 'weighted fermivel-sum: ', sum(temparr,  dim=2)

    end if

    !=========================!
    != Read in the spinvalue =!
    call read_spinvalue(filemode_int, nkpts1, nsqa, ndegen1, ispincomb, nvect, spinval1, nsym, isym)
    if(nkpts1*nsym/=nkpts) stop 'inconsistency in number of k-points'
    if(ndegen1/=inc%ndegen) stop 'inconsistency in ndegen'
    call expand_spinvalues(nsym,ndegen1,nsqa,nkpts1,spinval1,spinval)
    deallocate(isym, spinval1)

    !=============================!
    != Read in the torque values =!
    if(l_torqfile) then
      call read_torqvalue(filemode_int, nkpts1, ndegen1, torqval, nsym, isym)

      !next apply symmetries to torqval; as a first step, flatten the array
      itmp = size(torqval)/3
      allocate(arrtmp1(3,itmp), STAT=ierr)
      if(ierr/=0) stop 'problem alloaction arrtmp1'
      arrtmp1 = reshape(torqval, (/3, itmp/))

      deallocate(torqval)

      !now apply the symmetries
      call rotate_kpoints(symmetries%rotmat, nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
      if(nkpts2 /= nkpts*ndegen1) then
        write(*,*) 'nkpts2', nkpts2, 'nkpts', nkpts, 'ndegen1', ndegen1
        stop 'nkpts2 /= nkpts'
      endif

      !transform back to old shape
      allocate(torqval(3,ndegen1,nkpts1*nsym), STAT=ierr)
      if(ierr/=0) stop 'problem allocating torqval'
      torqval = reshape(arrtmp2,(/3,ndegen1,nkpts1*nsym/))

      deallocate(arrtmp1, arrtmp2, isym)
    endif

    !======================================!
    != Read in the torque values per atom =!
    if(l_torqfile_atom) then
      call read_torqvalue_atom(filemode_int, inc%natypd, nkpts1, ndegen1, torqval_atom, nsym, isym)

      !next apply symmetries to torqval_atom; as a first step, flatten the array
      itmp = size(torqval_atom)/3
      allocate(arrtmp1(3,itmp), STAT=ierr)
      if(ierr/=0) stop 'problem alloaction arrtmp1'
      arrtmp1 = reshape(torqval_atom, (/3, itmp/))

      deallocate(torqval_atom)

      !now apply the symmetries
      call rotate_kpoints(symmetries%rotmat, inc%natypd*nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
      if(nkpts2 /= nkpts*inc%natypd*ndegen1) stop 'nkpts2 /= nkpts*natpyd'

      !transform back to old shape
      allocate(torqval_atom(3,inc%natypd,ndegen1,nkpts1*nsym), STAT=ierr)
      if(ierr/=0) stop 'problem allocating torqval_atom'
      torqval_atom = reshape(arrtmp2,(/3,inc%natypd,ndegen1,nkpts1*nsym/))

      deallocate(arrtmp1, arrtmp2, isym)
    endif

    !====================================!
    != Read in the spin values per atom =!
    if(l_spinvecfile_atom) then
      call read_spinvec_atom(filemode_int, inc%natypd, nkpts1, ndegen1, spinvec_atom, nsym, isym)

      !next apply symmetries to spinvec_atom; as a first step, flatten the array
      itmp = size(spinvec_atom)/3
      allocate(arrtmp1(3,itmp), STAT=ierr)
      if(ierr/=0) stop 'problem alloaction arrtmp1'
      arrtmp1 = reshape(spinvec_atom, (/3, itmp/))

      deallocate(spinvec_atom)

      !now apply the symmetries
      call rotate_kpoints(symmetries%rotmat, inc%natypd*nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
      if(nkpts2 /= nkpts*inc%natypd*ndegen1) stop 'nkpts2 /= nkpts*natpyd'

      !transform back to old shape
      allocate(spinvec_atom(3,inc%natypd,ndegen1,nkpts1*nsym), STAT=ierr)
      if(ierr/=0) stop 'problem allocating spinvec_atom'
      spinvec_atom = reshape(arrtmp2,(/3,inc%natypd,ndegen1,nkpts1*nsym/))

      deallocate(arrtmp1, arrtmp2, isym)
    endif

    !====================================!
    != Read in the spin fluxes per atom =!
    if(l_spinfluxfile_atom) then
      call read_spinflux_atom(filemode_int, inc%natypd, nkpts1, ndegen1, spinflux_atom, nsym, isym)

      !next apply symmetries to spinflux_atom; as a first step, flatten the array
      itmp = size(spinflux_atom)/3
      allocate(arrtmp1(3,itmp), STAT=ierr)
      if(ierr/=0) stop 'problem alloaction arrtmp1'
      arrtmp1 = reshape(spinflux_atom, (/3, itmp/))

      deallocate(spinflux_atom)

      !now apply the symmetries
      call rotate_kpoints(symmetries%rotmat, inc%natypd*nkpts1*ndegen1, arrtmp1, nsym, isym, nkpts2, arrtmp2)
      if(nkpts2 /= nkpts*inc%natypd) stop 'nkpts2 /= nkpts*natpyd'

      !transform back to old shape
      allocate(spinflux_atom(3,inc%natypd,ndegen1,nkpts1*nsym), STAT=ierr)
      if(ierr/=0) stop 'problem allocating spinflux_atom'
      spinflux_atom = reshape(arrtmp2,(/3,inc%natypd,ndegen1,nkpts1*nsym/))

      deallocate(arrtmp1, arrtmp2, isym)
    endif

    !=====================================!
    !======== redefine symmetries ========!
    != (set to unit transformation only) =!
    nsym = 1
    allocate(isym(nsym))
    isym = (/ 1 /)

#ifdef CPP_TIMING
    call timing_stop('Read in of data')
#endif

    !============================================!
    !=== Create a subarray communication grid ===!
    != The big Pkk' array is subdivided into MxN submatrices.
    !=  Here, M is the number of rows and N the number of columns.
    !=   For example, M=2 and N=3 yields the following arrangement of processes.
    !=
    !=     <----------k_1 ------------>
    !=  <  (  id=0  |  id=2  |  id=4  )
    !=  |  (        |        |        )
    !=  |  (        |        |        )
    != k_2 (--------|--------|--------)
    !=  |  (  id=1  |  id=3  |  id=5  )
    !=  |  (        |        |        )
    !=  >  (        |        |        )
    !=
    !=  Then, the communicatior comm_row forms the two subgroups (0, 2 and 4) and (1, 3 and 5).
    !=  Therefore, the k2 axis is split in two parts, the eigenvectors are read in on
    !=    process 0 (master of group #1)
    !=    process 1 (master of group #2)
    !=  and then brodcastet to the other members of the group.

    ! determine size of subarray
    if(all(sca%subarr_inp>0))then
      subarr_dim = sca%subarr_inp
    else
      subarr_dim(1) = max( 1, int( sqrt( nranks + .01 ) ) )
      subarr_dim(2) = max( 1, nranks / subarr_dim(1) )
    end if
    if(myrank==master) write(*,*) 'subarr_dim=', subarr_dim
    ! create a cartesian communicator for dim(1) x dim(2) processes
    call create_subarr_comm( subarr_dim, myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col )
    call create_subarr( subarr_dim, nkpts, dataarr_lb, dataarr_ub, dataarr_nkpt)
    nkpt1 = dataarr_nkpt(myrank_grid,1)
    nkpt2 = dataarr_nkpt(myrank_grid,2)
    ioff1 = dataarr_lb(myrank_grid,1)
    ioff2 = dataarr_lb(myrank_grid,2)
!   write(*,'(5(A,I0))') 'myrank_grid=', myrank_grid, ', ioff1=',  ioff1, ', ioff2=', ioff2, ', myrank_row=', myrank_row, ', myrank_col=', myrank_col

    !***************************************************************************************
    if(compute_scattmat)then!*** Switches for calculation or read in of scattering matrix **
    !***************************************************************************************

      !=============================!
      != Read the eigenvector file =!
#ifdef CPP_TIMING
      call timing_start('Read in of eigenvectors')
#endif
      call read_eigv_part(inc, nsqa, ioff1, nkpt1, myMPI_comm_row, myrank_col, myMPI_comm_col, rveig_dim1)
      call read_eigv_part(inc, nsqa, ioff2, nkpt2, myMPI_comm_col, myrank_row, myMPI_comm_row, rveig_dim2)
#ifdef CPP_TIMING
      call timing_stop('Read in of eigenvectors')
#endif

      !===================================!
      != calculate the scattering matrix =!
#ifdef CPP_TIMING
      call timing_start('Calculation of Pkksub')
#endif
      call calculate_Pkksub( inc, lattice, impcls, nsqa, myrank, master, nkpt1, nkpt2, ioff1,     &
                           & ioff2, nkpts, rveig_dim1, rveig_dim2, kpoints, weights, Amat, Pkksub )
#ifdef CPP_TIMING
      call timing_stop('Calculation of Pkksub')
#endif

      !==============================!
      != save the scattering matrix =!
      if(sca%savepkk>0) then
#ifdef CPP_TIMING
        call timing_start('Writing of Pkksub to disk')
#endif
        call Pkkmpifile_write(myMPI_comm_grid, nkpts, nkpt1, nkpt2, ioff1, ioff2, inc%ndegen, nsqa, Pkksub)
#ifdef CPP_TIMING
        call timing_stop('Writing of Pkksub to disk')
#endif
        if(sca%savepkk==1)then
          call MPI_Finalize(ierr)
          stop 'Saved scattering matrix to file. Stopping.'
        end if!sca%savepkk==1
      end if!sca%savepkk>0

    !************************************************************************************
    else!compute_scattmat!*** Switches for calculation or read in of scattering matrix **
    !************************************************************************************

      !read the scattering matrix from file
      if(myrank==master) write(*,*) 'reading scatteringmatrix from file...'
      call Pkkmpifile_read(myMPI_comm_grid, nkpts, nkpt1, nkpt2, ioff1, ioff2, inc%ndegen, nsqa, Pkksub)

    !**************************************************************************************
    end if!compute_scattmat!*** Switches for calculation or read in of scattering matrix **
    !**************************************************************************************

#ifdef CPP_TIMING
    call timing_start('Calculating the lifetime')
#endif

    call calculate_lifetime_minmem( myrank_grid, myMPI_comm_grid, inc%ndegen, nsqa, nkpts, nkpt1,   &
                                  & nkpt2, ioff1, ioff2, BZVol, weights, Pkksub, tau, tau2, tau_avg )
#ifdef CPP_TIMING
    call timing_stop('Calculating the lifetime')
    call timing_start('Converging the meanfreepath[total]')
#endif

!   call meanfreepath_RTA(nkpts, nsqa, inc%ndegen, fermivel, tau, tau_avg, meanfreepath )
!   write(5000+myrank,*) 'check for NaNs in Pkksub:', sum(Pkksub)
    call converge_meanfreepath( myrank_grid, myMPI_comm_grid, nkpts, nkpt1, nkpt2,        &
                              & ioff1, ioff2, nsqa, inc%ndegen, BZVol, weights,           &
                              & fermivel, tau, tau_avg, Pkksub, meanfreepath,             &
                              & sca%gammamode, gammaval=sca%gammaval, impconc=sca%impconc )

    if(myrank==master .and. inc%ndegen==1) then
      open(unit=1325,file='output_lifetime',form='formatted',action='write')
      write(1325,'(A,I0)') '#nkpts= ', nkpts
      write(1325,'(A,ES25.16)') '# conversion factor Rydberg to inv.femtoseconds: ', RyToinvfs
      write(1325,'(A,6ES25.16)') '# averaged tau:', tau_avg(1,1)
      write(1325,'(A)') '# k[xyz], v[xyz], tau, tau2, lamba[xyz]' 
!     write(1325,'(A)') '# kx, ky, kz, vx, vy, vz, weight, tau[a.u.]' 
!     write(1325,'(A)') '# kx, ky, kz, vx, vy, vz, weight, tau[fs]' 
      do ikp=1,nkpts
        write(1325,'(14ES18.9)') kpoints(:,ikp), fermivel(:,ikp), tau(1,1,ikp), tau2(1,1,ikp), meanfreepath(:,:,1,ikp)
!       write(1325,'(8ES18.9)') kpoints(:,ikp), fermivel(:,ikp), weights(ikp), tau(1,1,ikp)
!       write(1325,'(8ES18.9)') kpoints(:,ikp), fermivel(:,ikp), weights(ikp), tau(1,1,ikp)/RyToinvfs
      end do!ikp
      close(1325)
    end if!myrank==master

    if(myrank==master) then
      open(unit=1326,file='lifetime_boltzmann.int.txt',form='formatted',action='write')
      write(1326,'(3I8)') nkpts, nsym, inc%ndegen
      write(1326,'(12I8)') isym
      do ikp=1,nkpts
        write(1326,'(10ES25.16)') tau(:,1,ikp), tau2(:,1,ikp), meanfreepath(:,:,1,ikp)
      end do!ikp
      close(1326)
    end if!myrank==master


#ifdef CPP_TIMING
    call timing_stop('Converging the meanfreepath[total]')
#endif

!   if(myrank==master) then
!     open(unit=1325,file='meanfreepath',form='formatted',action='write')
!     write(1325,'(3ES18.9)') meanfreepath
!     close(1325)
!   end if!myrank==master

    if(sca%lboltzmann==1) then
      call calc_condtensor( nkpts, nsqa, inc%ndegen, fermivel, spinval,         &
                           & meanfreepath, weights, lattice%alat, inc%nBZdim, chcond, spcond )

      if(allocated(torqval).or.allocated(torqval_atom).or.allocated(spinvec_atom).or.allocated(spinflux_atom)) then
        if(nsqa>1) write(*,*)"nsqa>1 : response functions will be computed only for sqa=1 !"
        call calc_response_functions_tensors( inc%natypd, nkpts, nsqa, inc%ndegen, fermivel, torqval,   &
                          & torqval_atom, spinvec_atom, spinflux_atom, meanfreepath, weights, lattice%alat, BZVol)
      end if!allocated(torqval).or.allocated(torqval_atom).or.allocated(spinvec_atom).or.allocated(spinflux_atom)
    end if!sca%lboltzmann==1

  end subroutine calc_boltzmann



  subroutine calc_condtensor( nkpts, nsqa, ndegen, fermivel, spinvalue,   &
                            & meanfreepath, weights, alat, nBZdim, chcond, spcond )
    use mod_mathtools, only: tpi
    use mod_mympi, only: myrank, master
    use mpi
    implicit none

    integer,          intent(in) :: nkpts, nsqa, ndegen, nBZdim
    double precision, intent(in) :: fermivel(3,nkpts), spinvalue(ndegen,nsqa,nkpts), meanfreepath(3,ndegen,nsqa,nkpts), weights(nkpts), alat


    double precision, allocatable, intent(out) ::  chcond(:,:,:), spcond(:,:,:)

    integer :: ihelp, ierr, ixyz1, ixyz2, ikp, isqa, ispin
    double precision :: chcond_tmp(3,3,nsqa), spcond_tmp(3,3,nsqa), dtmp, fac
    character(len=256) :: testfile1, testfile2, unitstr

    logical, parameter :: SIunits=.true.    
    double precision, parameter :: e2byhbar = 2.434134807664281d-4, abohr = 0.52917721092d-10


    !the equation for the [charge- and spin-] conductivity reads:
    !  \sigma = e^2/hbar * 1/(2pi)**3 * \int{ dS/|v_F| v_F * \lambda }
    !  The integration over the dimension of the Fermi-surface element yields a factor (2pi/a)**2,
    !   and the dimension of the mean free path is (a/2pi), thus this equation reduces to
    !  For a.u., this yields a factor of 2/(2pi)**2 / alat
    !  For SI, this yields [e^2/hbar = 2.4341E-4] * 1/(2pi)**2 / alat
    if(SIunits) then
      fac = e2byhbar/(tpi**2)/alat/abohr
      unitstr = 'SI units'
    else
      fac = 2d0/(tpi**2)/alat
      unitstr = 'atomic units'
    end if

    allocate( chcond(3,3,nsqa), &
            & spcond(3,3,nsqa), &
            & STAT=ierr  )
    if(ierr/=0) stop 'Problem allocating conductivity'

   !testwrite
!  write(testfile1,'(A,I0,A)') 'testfile1_',myrank,'.txt'
!  write(testfile2,'(A,I0,A)') 'testfile2_',myrank,'.txt'
!  open(unit=1235,file=testfile1,action='write',form='formatted')
!  open(unit=1236,file=testfile2,action='write',form='formatted')
!  do ikp=1,nkpts
!    write(1235,'(I0,4ES25.16)') ikp, fermivel(:,ikp), weights(ikp)
!    write(1236,'(I0,6ES25.16)') ikp, spinvalue(:,:,ikp)
!  end do
!  close(1235)
!  close(1236)

    chcond_tmp = 0d0
    chcond = 0d0
    do ikp=1,nkpts
      do ixyz1=1,3
        dtmp = fermivel(ixyz1,ikp)*weights(ikp)
        do isqa=1,nsqa
          do ispin=1,ndegen
            do ixyz2=1,3
              chcond_tmp(ixyz2,ixyz1,isqa) = chcond_tmp(ixyz2,ixyz1,isqa) + dtmp*meanfreepath(ixyz2,ispin,isqa,ikp)
!  TEST       chcond_tmp(ixyz2,ixyz1,isqa) = chcond_tmp(ixyz2,ixyz1,isqa) + dtmp*(1/(2*25*0.001/13.60569253))*fermivel(ixyz2,ikp)
            end do!ixyz2
          end do!ispin
        end do!isqa
      end do!ixyz1
    end do!ikp


!   if(ndegen==2)then
      spcond_tmp = 0d0
      spcond = 0d0
      do ikp=1,nkpts
        do ixyz1=1,3
          dtmp = fermivel(ixyz1,ikp)*weights(ikp)
          do isqa=1,nsqa
            do ispin=1,ndegen
              do ixyz2=1,3
                spcond_tmp(ixyz2,ixyz1,isqa) = spcond_tmp(ixyz2,ixyz1,isqa) + dtmp*meanfreepath(ixyz2,ispin,isqa,ikp)*spinvalue(ispin,isqa,ikp)
!  TEST         spcond_tmp(ixyz2,ixyz1,isqa) = spcond_tmp(ixyz2,ixyz1,isqa) + dtmp*(1/(2*25*0.001/13.60569253))*fermivel(ixyz2,ikp)*spinvalue(ispin,isqa,ikp)
               end do!ixyz2
            end do!ispin
          end do!isqa
        end do!ixyz1
      end do!ikp

!     do ikp=1,nkpts
!       do ixyz1=1,3
!         dtmp = fermivel(ixyz1,ikp)*weights(ikp)
!         do isqa=1,nsqa
!           do ixyz2=1,3
!             spcond_tmp(ixyz2,ixyz1,isqa) = spcond_tmp(ixyz2,ixyz1,isqa) + spinvalue(1,isqa,ikp)*dtmp*( meanfreepath(ixyz2,1,isqa,ikp) - meanfreepath(ixyz2,2,isqa,ikp) )
!           end do!ixyz2
!         end do!isqa
!       end do!ixyz1
!     end do!ikp
!   end if!ndegen==2

    !=== multiply in factors
    chcond = chcond_tmp*fac
    spcond = spcond_tmp*fac

    if(myrank==master)then
      write(*,'("Charge conductivity / 1 at.% in ",A,":")') trim(unitstr)
      if(nBZdim==2) write(*,'("Attention: these values have to be divided by the thickness of the film in units of the lattice constant to get proper units")')
      do isqa=1,nsqa
        write(*,'(2X,"isqa= ",I0)') isqa
        write(*,'(4X,"(",3ES18.9,")")') chcond(:,:,isqa)
      end do!isqa

      write(*,'("Spin conductivity / 1 at.% in ",A,":")') trim(unitstr)
      if(nBZdim==2) write(*,'("Attention: these values have to be divided by the thickness of the film in units of the lattice constant to get proper units")')
      do isqa=1,nsqa
        write(*,'(2X,"isqa= ",I0)') isqa
        write(*,'(4X,"(",3ES18.9,")")') spcond(:,:,isqa)
      end do!isqa

      write(*,'("AHE angle:")')
      do isqa=1,nsqa
        write(*,'(4X,"isqa=",I0," : ",ES18.9)') isqa, chcond(2,1,isqa)/chcond(1,1,isqa)
      end do!isqa

      write(*,'("SHE angle:")')
      do isqa=1,nsqa
        write(*,'(4X,"isqa=",I0," : ",ES18.9)') isqa, spcond(2,1,isqa)/chcond(1,1,isqa)
      end do!isqa
    end if!myrank==master

  end subroutine calc_condtensor




  subroutine calc_response_functions_tensors( natyp, nkpts, nsqa, ndegen, fermivel, torqval,   &
                            & torqval_atom, spinvec_atom, spinflux_atom, meanfreepath, weights, alat, BZVol)
    use mod_mathtools, only: tpi
    use mod_mympi, only: myrank, master
    use mpi
    implicit none

    integer,          intent(in) :: natyp, nkpts, nsqa, ndegen
    double precision, intent(in) :: fermivel(3,nkpts), meanfreepath(3,ndegen,nsqa,nkpts), weights(nkpts), alat, BZVol
    double precision, allocatable, intent(inout) :: torqval(:,:,:), torqval_atom(:,:,:,:), spinvec_atom(:,:,:,:), spinflux_atom(:,:,:,:) !torqval(3,ndegen,nkpts), torqval_atom(3,natyp,ndegen,nkpts), spinvec_atom(3,natyp,ndegen,nkpts), spinflux_atom(3,natyp,ndegen,nkpts)

    double precision, allocatable ::  torkance(:,:), torkance_atom(:,:,:), spinacc_resp_atom(:,:,:), spinflux_resp_atom(:,:,:)

    integer :: ihelp, ierr, ixyz1, ixyz2, ikp, isqa, ispin, iat
    double precision :: dtmp
    double precision :: torkance_tmp(3,3), torkance_atom_tmp(3,3,natyp), spinacc_resp_atom_tmp(3,3,natyp), spinflux_resp_atom_tmp(3,3,natyp)
    double precision :: fac_torkance, fac_spinacc, fac_spinflux
    character(len=256) :: testfile1, testfile2, unitstr_torkance, unitstr_spinacc, unitstr_spinflux

    logical, parameter :: SIunits=.false.    
    double precision, parameter :: e = 1.602176565d-19, abohr = 0.52917721092d-10, Rd=2.1798741d-18, mu_B=9.27400968d-24

    integer, parameter :: iounit=13521

    ! the equation for the torkance reads:
    !  \torkance = e/hbar * 1/(BZVol) * \int{ dS/|v_F| torq * \lambda }
    ! The weights of k-pts have the unit of (2pi/a)**2
    ! The fermi velocity has the unit of a/2pi*Ry/hbar
    ! The vector mean free path has the unit of a/pi/second
    ! The volume of the BZ has the volume of (a/2pi)**3
    ! The torkance has the unit of Ry
    ! In total, this yields an additional factor of hbar*a/2pi in the previous equation for the torkance, i.e. :
    !  \torkance = e * 1/(BZVol) * a/2pi * \sum{ area/|v_F| torq * \lambda }
    if(SIunits) then
      if(allocated(torqval)) fac_torkance = e*alat/(tpi)*abohr/BZVol
      if(allocated(spinvec_atom)) fac_spinacc = e*alat/(tpi)*abohr/BZVol*mu_B/Rd
      if(allocated(spinflux_atom)) fac_spinflux = e*alat/(tpi)*abohr/BZVol
      unitstr_torkance = 'SI units'
      unitstr_spinacc  = 'SI units'
      unitstr_spinflux = 'SI units'
    else
      if(allocated(torqval)) fac_torkance = alat/(tpi)/BZVol
      if(allocated(spinvec_atom)) fac_spinacc = alat/(tpi)/BZVol
      if(allocated(spinflux_atom)) fac_spinflux = alat/(tpi)/BZVol
      unitstr_torkance = 'unit of (e a_0)'
      unitstr_spinacc  = 'unit of (e a_0 mu_B)/Rd'
      unitstr_spinflux = 'unit of (e a_0)'
    end if!SIunits

    if(allocated(torqval))then
      allocate( torkance(3,3), STAT=ierr  )
      if(ierr/=0) stop 'Problem allocating torkance'
      torkance_tmp = 0d0
      torkance     = 0d0
    end if!allocated(torqval)

    if(allocated(torqval_atom))then
      allocate( torkance_atom(3,3,natyp), STAT=ierr  )
      if(ierr/=0) stop 'Problem allocating torkance_atom'
      torkance_atom_tmp = 0d0
      torkance_atom     = 0d0
    end if!allocated(torqval_atom)

    if(allocated(spinvec_atom))then
      allocate( spinacc_resp_atom(3,3,natyp), STAT=ierr  )
      if(ierr/=0) stop 'Problem allocating spinacc_resp_atom'
      spinacc_resp_atom_tmp = 0d0
      spinacc_resp_atom     = 0d0
    end if!allocated(spinvec_atom)

    if(allocated(spinflux_atom))then
      allocate( spinflux_resp_atom(3,3,natyp), STAT=ierr  )
      if(ierr/=0) stop 'Problem allocating spinflux_resp_atom'
      spinflux_resp_atom_tmp = 0d0
      spinflux_resp_atom     = 0d0
    end if!allocated(spinflux_atom)

    do ikp=1,nkpts
      do ixyz2=1,3
          do ispin=1,ndegen
            dtmp=weights(ikp)*meanfreepath(ixyz2,ispin,1,ikp)
            do ixyz1=1,3
              if(allocated(torqval))       torkance_tmp(ixyz2,ixyz1)                = torkance_tmp(ixyz2,ixyz1)               &
                                                                                  & + torqval(ixyz1,ispin,ikp)*dtmp
!  TEST       torkance_tmp(ixyz2,ixyz1,isqa) = torkance_tmp(ixyz2,ixyz1,isqa) + dtmp*(1/(2*25*0.001/13.60569253))*fermivel(ixyz2,ikp)
             do iat=1,natyp
              if(allocated(torqval_atom))  torkance_atom_tmp(ixyz2,ixyz1,iat)       = torkance_atom_tmp(ixyz2,ixyz1,iat)      &
                                                                                  & + torqval_atom(ixyz1,iat,ispin,ikp)*dtmp
              if(allocated(spinvec_atom))  spinacc_resp_atom_tmp(ixyz2,ixyz1,iat)   = spinacc_resp_atom_tmp(ixyz2,ixyz1,iat)  &
                                                                                  & + spinvec_atom(ixyz1,iat,ispin,ikp)*dtmp
              if(allocated(spinflux_atom)) spinflux_resp_atom_tmp(ixyz2,ixyz1,iat)  = spinflux_resp_atom_tmp(ixyz2,ixyz1,iat) &
                                                                                  & + spinflux_atom(ixyz1,iat,ispin,ikp)*dtmp
             end do!iat
            end do!ixyz2
          end do!ispin
      end do!ixyz1
    end do!ikp

    if(allocated(torqval))then
      torkance = torkance_tmp * fac_torkance

      if(myrank==master)then
        write(*,'("Torkance / 1 at.% in ",A,":")') trim(unitstr_torkance)
          write(*,'(4X,"(",3ES18.9,")")') torkance
      end if!myrank==master
    endif!allocated(torqval)

    if(allocated(torqval_atom))then
      torkance_atom = torkance_atom_tmp*fac_torkance

      if(myrank==master)then
        open(unit=iounit,file="torkance.SCAT.int.dat",form='formatted',action='write')
        write(iounit,'("#Torkance / 1 at.% in ",A,":")') trim(unitstr_torkance)
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') torkance_atom
        close(iounit)
      end if!myrank==master
    end if!allocated(torqval_atom)

    if(allocated(spinvec_atom))then
      spinacc_resp_atom = spinacc_resp_atom_tmp*fac_spinacc

      if(myrank==master)then
        open(unit=iounit,file="spinacc_resp.SCAT.int.dat",form='formatted',action='write')
        write(iounit,'("#Response of the spin accumulation to the electric field / 1 at.% in ",A,":")') trim(unitstr_spinacc)
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') spinacc_resp_atom
        close(iounit)
      end if!myrank==master
    end if!allocated(spinvec_atom)

    if(allocated(spinflux_atom))then
      spinflux_resp_atom = spinflux_resp_atom_tmp*fac_spinflux

      if(myrank==master)then
        open(unit=iounit,file="spinflux_resp.SCAT.int.dat",form='formatted',action='write')
        write(iounit,'("#Response of the spin fluxes to the electric field / 1 at.% in ",A,":")') trim(unitstr_spinflux)
        write(iounit,'(1I8)') natyp
        write(iounit,'(9(ES25.16))') spinflux_resp_atom
        close(iounit)
      end if!myrank==master
    end if!allocated(spinflux_atom)

  end subroutine calc_response_functions_tensors



  subroutine meanfreepath_RTA(nkpts, nsqa, ndegen, fermivel, tau, tau_avg, meanfreepath_new )
    implicit none

    integer,          intent(in) :: nkpts, nsqa, ndegen
    double precision, intent(in) :: fermivel(3,nkpts), tau(ndegen,nsqa,nkpts), tau_avg(ndegen,nsqa)
    double precision, allocatable, intent(out) :: meanfreepath_new(:,:,:,:)
 
    integer :: ikp, isqa, ispin, ixyz, ierr

    !=== init
    allocate( meanfreepath_new(3,ndegen,nsqa,nkpts), &
            & STAT=ierr                              )
    if(ierr/=0) stop 'Problem allocating meanfreepath'
    meanfreepath_new = 0d0

    do ikp=1,nkpts
      do isqa=1,nsqa
        do ispin=1,ndegen
          do ixyz=1,3
            meanfreepath_new(ixyz,ispin,isqa,ikp) = tau(ispin,isqa,ikp)*fermivel(ixyz,ikp)
!           meanfreepath_new(ixyz,ispin,isqa,ikp) = tau_avg(ispin,isqa)*fermivel(ixyz,ikp)
          end do!ixyz
        end do!ispin
      end do!isqa
    end do!ikp


  end subroutine meanfreepath_RTA




  subroutine converge_meanfreepath( myrank_grid, comm_grid, nkpts, nkpt1, nkpt2, &
                                  & ioff1, ioff2, nsqa, ndegen, BZVol, weights,  &
                                  & fermivel, tau, tau_avg, Pkksub,              &
                                  & meanfreepath_new, add_gamma_mode, gammaval,  &
                                  & impconc                                     )

#ifdef CPP_TIMING
    use mod_timing, only: timing_start, timing_stop, timing_pause
#endif
    use mod_mympi, only: master, myrank
    use mpi
    implicit none

    integer,          intent(in) :: myrank_grid, comm_grid, nkpts, nkpt1, nkpt2, ioff1, ioff2, nsqa, ndegen
    double precision, intent(in) :: BZVol, weights(nkpts), fermivel(3,nkpts), &
                                  & tau(ndegen,nsqa,nkpts), tau_avg(ndegen,nsqa), Pkksub(ndegen,nkpt1,ndegen,nsqa,nkpt2)
    integer,                       intent(in)  :: add_gamma_mode
    double precision, optional,    intent(in)  :: gammaval, impconc

    double precision, allocatable, intent(out) :: meanfreepath_new(:,:,:,:)

    double precision, allocatable :: meanfreepath_old_t(:,:,:,:), &
                                   & meanfreepath_old(:,:,:,:),   &
                                   & meanfreepath_tmp(:,:,:,:),   &
                                   & fveltau_sum(:,:,:),          &
                                   & fveltau_sum_w(:,:,:)

    double precision :: tau_tilde(ndegen,nsqa,nkpts), Pkksub_tilde(ndegen,nkpt1,ndegen,nsqa,nkpt2)
    double precision :: dos_Ef(nsqa)

    integer :: ikp, ikp1, ikp2, isqa, ispin, ispin1, ispin2, ixyz, iter, ierr, ihelp
    double precision :: dist, distxyz(3)

    !parameter
    integer, parameter :: nitermax=100
    double precision, parameter :: alpha=1.0d0, eps=1d-8

    !==========================!
    !====     S T A R T     ===!
    !==========================!
    !=== BOLTZMANN EQUATION ===!
    !=== ITERATIVE SOLUTION ===!
    !==========================!

    !=== init
    allocate( meanfreepath_old_t(ndegen,nsqa,nkpts,3), &
            & meanfreepath_old(3,ndegen,nsqa,nkpts),   &
            & meanfreepath_new(3,ndegen,nsqa,nkpts),   &
            & meanfreepath_tmp(3,ndegen,nsqa,nkpts),   &
            & fveltau_sum(3,ndegen,nsqa),              &
            & fveltau_sum_w(3,ndegen,nsqa),            &
            & STAT=ierr )
    if(ierr/=0) stop 'Problem allocating meanfreepath'
    meanfreepath_old = 0d0
    meanfreepath_new = 0d0

    !=== compute weighted sums
    fveltau_sum   = 0d0
    fveltau_sum_w = 0d0
    do ikp=1,nkpts
      do isqa=1,nsqa
        do ispin=1,ndegen
          do ixyz=1,3
            fveltau_sum(ixyz,ispin,isqa)   = fveltau_sum(ixyz,ispin,isqa)   + tau(ispin,isqa,ikp)*fermivel(ixyz,ikp)
            fveltau_sum_w(ixyz,ispin,isqa) = fveltau_sum_w(ixyz,ispin,isqa) + tau(ispin,isqa,ikp)*fermivel(ixyz,ikp)*weights(ikp)
          end do!ixyz
        end do!ispin
      end do!isqa
    end do!ikp
    if(myrank==master) write(*,'("         fermivel-sum with taus: ",3ES25.16)') fveltau_sum
    if(myrank==master) write(*,'("weighted fermivel-sum with taus: ",3ES25.16)') fveltau_sum_w

    ! if add_gamma_mode==1 or 2, a constant contribution hbar/(2*gamma) is added to the relaxation times
    if(add_gamma_mode==1.or.add_gamma_mode==2)then
      if(.not.present(gammaval)) stop 'gammaval should be present if add_gamma_mode==1 or 2'
      tau_tilde(:,:,:) = tau(:,:,:)/impconc        ! rescale the tau for the specified impurity concentration
      tau_tilde(:,:,:) = 1/(1/tau_tilde(:,:,:) + 2*gammaval)   ! 1/tau_tilde = 1/tau + 2*gamma/hbar
    else
      tau_tilde = tau
    end if!add_gamma_mode==1.or.add_gamma_mode==2


    !=== first guess for the mean free path
#ifdef CPP_TIMING
    call timing_start('  Boltzmann eq first guess')
#endif
!$omp parallel private(ikp, isqa, ispin,ixyz)
!$omp do collapse (4)
    do ikp=1,nkpts
      do isqa=1,nsqa
        do ispin=1,ndegen
          do ixyz=1,3
!           meanfreepath_new(ixyz,ispin,isqa,ikp) = tau(ispin,isqa,ikp)*fermivel(ixyz,ikp)/fveltau_sum(ixyz,ispin,isqa)
            if(add_gamma_mode==1.or.add_gamma_mode==2)then
              meanfreepath_new(ixyz,ispin,isqa,ikp) = tau_tilde(ispin,isqa,ikp)*fermivel(ixyz,ikp)
            else
              meanfreepath_new(ixyz,ispin,isqa,ikp) = tau(ispin,isqa,ikp)*fermivel(ixyz,ikp)
            end if!present(add_gamma_mode).and.present(gammaval)
          end do!ixyz
        end do!ispin
      end do!isqa
    end do!ikp
!$omp end do
!$omp barrier
!$omp end parallel
#ifdef CPP_TIMING
    call timing_stop('  Boltzmann eq first guess')
#endif

!   write(*,'("rank ",I0," finds ",I0," NaNs in ",A)') myrank, count(IsNan(Pkksub)), 'Pkksub'
!   write(*,'("rank ",I0," finds ",I0," NaNs in ",A)') myrank, count(IsNan(weights)), 'weights'
!   write(*,'("rank ",I0," finds ",I0," NaNs in ",A)') myrank, count(IsNan(fermivel)), 'fermivel'
!   write(*,'("rank ",I0," finds ",I0," NaNs in ",A)') myrank, count(IsNan(tau)), 'tau'

    !=== compute density of states at Ef (needed if add_gamma_mode==2 only)
    if(add_gamma_mode==2)then
#ifdef CPP_TIMING
        call timing_start('  Density of states at Ef')
#endif

      dos_Ef = 0d0 
      do ikp=1,nkpts
        do isqa=1,nsqa
          dos_Ef(isqa) = dos_Ef(isqa) + weights(ikp)*ndegen ! ndegen accounts for spin degeneracy
        end do!isqa
      end do!ikp

#ifdef CPP_TIMING
    call timing_stop('  Density of states at Ef')
#endif
    end if!add_gamma_mode==2

    ! if add_gamma_mode==2, a constant contribution 2*Gamma/(hbar*dos(Ef)) is added to the P matrix (for consistency with tau_tilde)
    ! and the original (scattering) Pkksub is multiplied my the impurity concentration (in %) which may now differ from 1
    if(add_gamma_mode==2)then
      do isqa=1,nsqa
        Pkksub_tilde(:,:,:,isqa,:) = Pkksub(:,:,:,isqa,:)*impconc + 2*gammaval/dos_Ef(isqa)  !Pkk'_tilde = Pkk' + 2*gamma/(hbar*dos)
      end do!isqa
    end if!add_gamma_mode==2


    !======== BEGIN =========!
    !=== Convergency loop ===!
    !========================!
    iter = 0
    dist = 1d38
    distxyz = 1d38
    do while (dist > eps .and. iter<nitermax)

      iter = iter+1

      meanfreepath_old = meanfreepath_new

      do ikp=1,nkpts
        do isqa=1,nsqa
          do ispin=1,ndegen
            do ixyz=1,3
              meanfreepath_old_t(ispin,isqa,ikp,ixyz) = meanfreepath_new(ixyz,ispin,isqa,ikp)
            end do!ixyz
          end do!ispin
        end do!isqa
      end do!ikp

      !========= BEGIN ==========!
      !==== compute r.h.s. of ===!
      !=== boltzmann equation ===!

#ifdef CPP_TIMING
      call timing_start('  Boltzmann eq. big loop')
#endif
      !calculate sum over k1 on this processor, \sum_{k1}{ P(k1,k2) \lambda(k1)}
      meanfreepath_tmp = 0d0
!$omp parallel default(shared) private(ikp2,isqa,ispin2,ikp1,ispin1,ixyz)
!$omp do collapse (4)
      do ikp2=1,nkpt2
        do isqa=1,nsqa
          do ispin2=1,ndegen
            do ixyz=1,3
              do ikp1=1,nkpt1
                do ispin1=1,ndegen

                  if(add_gamma_mode==2) then
                    meanfreepath_tmp(ixyz,ispin2,isqa,ikp2+ioff2) = meanfreepath_tmp(ixyz,ispin2,isqa,ikp2+ioff2) + &
                    & Pkksub_tilde(ispin1,ikp1,ispin2,isqa,ikp2)*weights(ikp1+ioff1)*meanfreepath_old_t(ispin1,isqa,ikp1+ioff1,ixyz)
                  else 
                    meanfreepath_tmp(ixyz,ispin2,isqa,ikp2+ioff2) = meanfreepath_tmp(ixyz,ispin2,isqa,ikp2+ioff2) + &
                    & Pkksub(ispin1,ikp1,ispin2,isqa,ikp2)*weights(ikp1+ioff1)*meanfreepath_old_t(ispin1,isqa,ikp1+ioff1,ixyz)
                  endif!add_gamma_mode==2
                end do!ispin1
              end do!ikp1
            end do!ixyz
          end do!ispin2
        end do!isqa
      end do!ikp2
!$omp end do
!$omp end parallel
#ifdef CPP_TIMING
      call timing_pause('  Boltzmann eq. big loop')
#endif

      !sum up the sum over k1 from different processors
      meanfreepath_new = 0d0
      ihelp = 3*ndegen*nsqa
      call MPI_Allreduce( meanfreepath_tmp, meanfreepath_new, ihelp*nkpts, MPI_DOUBLE_PRECISION, MPI_SUM, comm_grid, ierr )
      if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allreduce( meanfreepath_tmp )'

      !add fermivelocity and lifetime
#ifdef CPP_TIMING
      call timing_start('  Boltzmann eq. small loop')
#endif
      do ikp=1,nkpts
        do isqa=1,nsqa
          do ispin=1,ndegen
            if(add_gamma_mode==1.or.add_gamma_mode==2)then
              meanfreepath_new(:,ispin,isqa,ikp) = ( fermivel(:,ikp) + meanfreepath_new(:,ispin,isqa,ikp)/BZVol ) * tau_tilde(ispin,isqa,ikp)
            else
              meanfreepath_new(:,ispin,isqa,ikp) = ( fermivel(:,ikp) + meanfreepath_new(:,ispin,isqa,ikp)/BZVol ) * tau(ispin,isqa,ikp)
            end if!add_gamma_mode==1.or.add_gamma_mode==2
          end do!ispin
        end do!isqa
      end do!ikp
#ifdef CPP_TIMING
      call timing_pause('  Boltzmann eq. small loop')
#endif
      !==== compute r.h.s. of ===!
      !=== boltzmann equation ===!
      !=========  END  ==========!

!     !perform mixing 
!     meanfreepath_new = alpha*meanfreepath_new + (1d0-alpha)*meanfreepath_old

      !calculate distance
      dist = sqrt( sum( (meanfreepath_new - meanfreepath_old)**2 ) )/(3*ndegen*nsqa*nkpts)
      do ixyz=1,3
        distxyz(ixyz) = sqrt( sum( (meanfreepath_new(ixyz,:,:,:) - meanfreepath_old(ixyz,:,:,:))**2 ) )/(ndegen*nsqa*nkpts)
      end do!ixyz
      if(myrank_grid==master) then
         write(*,'("distance=",ES18.9," after iteration ",I0," - distxyz=",3ES18.9," - random lambda=",3ES18.9)') dist, iter, distxyz, meanfreepath_new(:,1,1,361)
      end if

    end do! while(dist > eps .and. iter<nitermax)
    !========================!
    !=== Convergency loop ===!
    !========  END  =========!
#ifdef CPP_TIMING
    call timing_stop('  Boltzmann eq. big loop')
    call timing_stop('  Boltzmann eq. small loop')
#endif

    deallocate( meanfreepath_old, meanfreepath_tmp)

  end subroutine converge_meanfreepath





  subroutine calculate_lifetime_minmem( myrank_grid, comm_grid, ndegen, nsqa, nkpts, nkpt1, nkpt2, ioff1, ioff2, BZVol, weights, Pkksub, tau, tau2, tau_avg )
    use mod_mympi, only: master
    use mpi
    implicit none

    integer,          intent(in) :: myrank_grid, comm_grid
    integer,          intent(in) :: ndegen, nsqa, nkpts, nkpt1, nkpt2, ioff1, ioff2
    double precision, intent(in) :: BZVol, weights(nkpts), Pkksub(ndegen,nkpt1,ndegen,nsqa,nkpt2)

    double precision, allocatable, intent(out) :: tau(:,:,:), tau2(:,:,:), tau_avg(:,:)

    integer :: ikp1, ikp2, ispin1, ispin2, isqa, ierr
    double precision :: tau_tmp(ndegen,nsqa,nkpts),  tau_inv(ndegen,nsqa,nkpts), &
                      & tau_tmp2(ndegen,nsqa,nkpts), tau_inv2(ndegen,nsqa,nkpts), &
                      & tau_avg_tmp(ndegen,nsqa),    tau_avg_inv(ndegen,nsqa)
    double precision :: tempsum1, tempsum2,   &
                      & tauinv_flip(nsqa),    &
                      & tauinv_cons(nsqa),    &
                      & tauinv_flip_sum(nsqa),&
                      & tauinv_cons_sum(nsqa)

    character(len=256) :: filename, formatstring

    double precision, parameter :: RyToinvfs = 20.67068667282055d0

    allocate( tau(ndegen,nsqa,nkpts), tau2(ndegen,nsqa,nkpts), tau_avg(ndegen,nsqa), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating tau'

    !initialize the tempprary arrays (sum on each processor)
    tau_tmp  = 0d0
    tau_tmp2 = 0d0
    tau_avg_tmp = 0d0

    !initialize the tempprary arrays (combined sum from all processor)
    tau_inv  = 0d0
    tau_inv2 = 0d0
    tau_avg_inv = 0d0

!!!!$omp parallel private(ikp1, ispin1)
    do ikp2=1,nkpt2
      do isqa=1,nsqa
        do ispin2=1,ndegen
!!!!$omp do collapse (2)
          do ikp1=1,nkpt1
            do ispin1=1,ndegen
              ! tau_k  = sum_k' P(k,k')
              tau_tmp(ispin1,isqa,ikp1+ioff1) = tau_tmp(ispin1,isqa,ikp1+ioff1) + & 
                                             & Pkksub(ispin1,ikp1,ispin2,isqa,ikp2)*weights(ikp2+ioff2)
              ! tau_k' = sum_k  P(k,k')
              tau_tmp2(ispin2,isqa,ikp2+ioff2) = tau_tmp2(ispin2,isqa,ikp2+ioff2) + & 
                                             & Pkksub(ispin1,ikp1,ispin2,isqa,ikp2)*weights(ikp1+ioff1)
              ! tau_av = sum_{kk'} P(k,k')
              tau_avg_tmp(ispin1,isqa) = tau_avg_tmp(ispin1,isqa) + &
                                             & Pkksub(ispin1,ikp1,ispin2,isqa,ikp2)*weights(ikp2+ioff2)*weights(ikp1+ioff1)
            end do!ispin1
          end do!ikp1
!!!!$omp end do
        end do!ispin2
      end do!isqa
    end do!ikp2

!!!!$omp end parallel
  
    call MPI_Allreduce( tau_tmp,  tau_inv,  nkpts*nsqa*ndegen, MPI_DOUBLE_PRECISION, MPI_SUM, comm_grid, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allreduce( tau_tmp )'

    call MPI_Allreduce( tau_tmp2, tau_inv2, nkpts*nsqa*ndegen, MPI_DOUBLE_PRECISION, MPI_SUM, comm_grid, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allreduce( tau_tmp2 )'

    call MPI_Allreduce( tau_avg_tmp, tau_avg_inv, nsqa*ndegen, MPI_DOUBLE_PRECISION, MPI_SUM, comm_grid, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allreduce( tau_tmp2 )'

    tau  = 1d0/tau_inv*BZVol
    tau2 = 1d0/tau_inv2*BZVol
    tau_avg = sum(weights)*BZvol/tau_avg_inv

!   write(*,'(A,I0,A,ES18.9)') 'myrank=', myrank_grid, ', tau_avg [fs] =', tau_avg/RyToinvfs

    !=============
    !=== START ===
    ! perform additional calulations

    !init
    tauinv_cons     = 0d0
    tauinv_cons_sum = 0d0
    tauinv_flip     = 0d0
    tauinv_flip_sum = 0d0

    !calculate
    select case (ndegen)
    case (2)
      do ikp2=1,nkpt2
        do isqa=1,nsqa
          tempsum1 = 0d0
          tempsum2 = 0d0
          do ikp1=1,nkpt1
            tempsum1 = tempsum1 + ( Pkksub(1,ikp1,1,isqa,ikp2) + Pkksub(2,ikp1,2,isqa,ikp2) )*weights(ikp1+ioff1)
            tempsum2 = tempsum2 + ( Pkksub(1,ikp1,2,isqa,ikp2) + Pkksub(2,ikp1,1,isqa,ikp2) )*weights(ikp1+ioff1)
          end do!ikp1
          tauinv_cons(isqa) = tauinv_cons(isqa) + tempsum1*weights(ikp2+ioff2)
          tauinv_flip(isqa) = tauinv_flip(isqa) + tempsum2*weights(ikp2+ioff2)
        end do!isqa
      end do!ikp2
    case (1)
      do ikp2=1,nkpt2
        do isqa=1,nsqa
          tempsum1 = 0d0
          do ikp1=1,nkpt1
            tempsum1 = tempsum1 + Pkksub(1,ikp1,1,isqa,ikp2)*weights(ikp1+ioff1)
          end do!ikp1
          tauinv_cons(isqa) = tauinv_cons(isqa) + tempsum1*weights(ikp2+ioff2)
        end do!isqa
      end do!ikp2
    case default; stop 'case(ndegen) not known'
    end select


    !collect results
    call MPI_Reduce( tauinv_cons, tauinv_cons_sum, nsqa, MPI_DOUBLE_PRECISION, MPI_SUM, master, comm_grid, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allreduce( tauinv_cons )'
    if(ndegen==2) then
      call MPI_Reduce( tauinv_flip, tauinv_flip_sum, nsqa, MPI_DOUBLE_PRECISION, MPI_SUM, master, comm_grid, ierr )
      if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allreduce( tauinv_flip )'
    end if!ndegen==2

    !output
    if(myrank_grid==master)then

      write(*,'(A)') 'momentum relaxation times tau [fs] / 1 at.%:'
      do isqa=1,nsqa
        write(*,'(2X,"isqa=",I0,", ",ES18.9)') isqa, 2d0/tauinv_cons_sum(isqa)*(BZVol*sum(weights))/RyToinvfs
      end do!isqa

      if(ndegen==2)then
        write(*,'(A)') 'spin relaxation times tau [fs] / 1 at.%:'
        do isqa=1,nsqa
          write(*,'(2X,"isqa=",I0,", ",ES18.9)') isqa, 1d0/tauinv_flip_sum(isqa)*(BZVol*sum(weights))/RyToinvfs
        end do!isqa
      end if!ndegen==2

      do isqa=1,nsqa
        write(filename,'(A,I0,A)') 'test_lifetimes_sqa=',isqa,'.txt'
        open(unit=1358,file=trim(filename),form='formatted',action='write')
        do ikp2=1,nkpts
          write(formatstring,'(A,I0,A,I0,A)') '(',2*ndegen,'ES25.16,',2*ndegen,'F10.3)'
          write(1358,formatstring) tau_inv(:,isqa,ikp2), tau_inv2(:,isqa,ikp2), 2d0*( tau_inv(:,isqa,ikp2)-tau_inv2(:,isqa,ikp2) )/(tau_inv(:,isqa,ikp2)+tau_inv2(:,isqa,ikp2))*100
        end do!ikp2
        close(1358)
      end do!isqa

    end if!myrank_grid==master

  end subroutine calculate_lifetime_minmem




  subroutine calculate_Pkksub(inc, lattice, impcls, nsqa, myrank, master, nkpt1, nkpt2, ioff1, ioff2, ntot, rveig1, rveig2, kpoints, weights, Amat, Pkksub)

    use type_inc,      only: inc_type
    use type_data,     only: lattice_type
    use mod_mathtools, only: pi, tpi
    use mod_iohelp,    only: getBZvolume
    use mpi
    implicit none

    type(inc_TYPE),    intent(in) :: inc
    type(lattice_TYPE),intent(in) :: lattice
    type(impcls_TYPE), intent(in) :: impcls(sca%naverage)
    integer,           intent(in) :: nsqa, myrank, master, nkpt1, nkpt2, ioff1, ioff2, ntot
    double precision,  intent(in) :: kpoints(3,ntot), weights(ntot)
    double complex,    intent(in) :: Amat(impcls(1)%clmso,impcls(1)%clmso,sca%naverage), &
                                   & rveig1(inc%lmmaxso, inc%natypd, inc%ndegen, nsqa, nkpt1),& 
                                   & rveig2(inc%lmmaxso, inc%natypd, inc%ndegen, nsqa, nkpt2)
    double precision, allocatable, intent(out) :: Pkksub(:,:,:,:,:)

    double precision, allocatable :: optical_left(:,:,:),&
                                   & optical_left_all(:,:,:),&
                                   & optical_right(:,:,:),&
                                   & optical_right_all(:,:,:),&
                                   & optical_right_thread(:,:,:)

    integer          :: clmso, ierr, ikp1, ikp2, ispin1, ispin2, isqa, ihelp, iaverage
    double precision :: Tkk_abs(nkpt1,nkpt2), fac(sca%naverage), righttmp, BZvol
    double complex   :: Tkk_tmp(nkpt1,nkpt2), &
                      & rvcls1(impcls(1)%clmso,nkpt1),&
                      & rvcls2(impcls(1)%clmso,nkpt2),&
                      & Ctmp(nkpt2,impcls(1)%clmso)
    character(len=256) :: filename
    logical, parameter :: l_optical=.true.

    ! parameter
    double complex, parameter :: CZERO=(0d0, 0d0), CONE=(1d0,0d0), CI=(0d0,1d0)

    ! weights are rescaled such that sum(weights) = 1
    fac     = 2d0*pi*0.01d0*(sca%weight_imp/sum(sca%weight_imp)) !concentration = 1 atom percent
    clmso   = impcls(1)%clmso
    BZVol   = getBZvolume(lattice%recbv)


    allocate( Pkksub(inc%ndegen,nkpt1,inc%ndegen,nsqa,nkpt2),&
            & STAT=ierr )
    if(ierr/=0) stop 'Problem allocating Pkksub'
    Pkksub = 0d0

    if(l_optical) then
      allocate( optical_left(inc%ndegen,nsqa,ntot),&
              & optical_left_all(inc%ndegen,nsqa,ntot),&
              & optical_right(inc%ndegen,nsqa,ntot),&
              & optical_right_all(inc%ndegen,nsqa,ntot),&
              & optical_right_thread(inc%ndegen,nsqa,ntot),&
              & STAT=ierr )
      if(ierr/=0) stop 'Problem allocating optical arrays'
      optical_left      = 0d0
      optical_right     = 0d0
      optical_left_all  = 0d0
      optical_right_all = 0d0
    end if!l_optical

    !======= BEGIN ======!
    !=== Calculations ===!

    do iaverage=1,sca%naverage

      if(l_optical) then
        optical_left      = 0d0
        optical_right     = 0d0
        optical_left_all  = 0d0
        optical_right_all = 0d0
      end if!l_optical

      if(myrank==master) write(*,'("Start calculating the Pkk'' submatrix")')
!      if(myrank==master) write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
!      if(myrank==master) write(*,FMT=190) !beginning of statusbar

      !calculate Pkk'
      do isqa=1,nsqa
        do ispin2=1,inc%ndegen

          call extend_coeffvector2cluster_allkpt( inc, lattice, impcls(iaverage), nkpt2, kpoints(:,ioff2+1:ioff2+nkpt2), rveig2(:,:,ispin2,isqa,:), rvcls2)

          ! Ctmp = conjg(c_k') * Amat
          call ZGEMM( 'C','N', nkpt2, clmso, clmso, CONE, rvcls2, clmso, Amat(:,:,iaverage), clmso, CZERO, Ctmp, nkpt2)

          do ispin1=1,inc%ndegen

            call extend_coeffvector2cluster_allkpt( inc, lattice, impcls(iaverage), nkpt1, kpoints(:,ioff1+1:ioff1+nkpt1), rveig1(:,:,ispin1,isqa,:), rvcls1 )

            ! Tkk' = conjg(c_k') * Amat * c_k
            call ZGEMM( 'T','T', nkpt1, nkpt2, clmso, CONE, rvcls1, clmso, Ctmp, nkpt2, CZERO, Tkk_tmp, nkpt1)

            ! Pkk' = |Tkk'|^2
            Tkk_abs(:,:) = dble(Tkk_tmp(:,:))**2 + aimag(Tkk_tmp(:,:))**2
            Pkksub(ispin1,:,ispin2,isqa,:) = Pkksub(ispin1,:,ispin2,isqa,:) + Tkk_abs(:,:)*fac(iaverage)
          end do!ispin1
        end do!ispin2
      end do!isqa

      if(myrank==master) write(*,'("Checking optical theorem")')
!$omp parallel private(ikp1, ispin1, optical_right_thread)
      if(l_optical) optical_right_thread = 0d0
      do ikp2=1,nkpt2
        do isqa=1,nsqa
          do ispin2=1,inc%ndegen
!$omp do collapse (2)
            do ikp1=1,nkpt1
              do ispin1=1,inc%ndegen
                if(l_optical)then
                  if((ikp1+ioff1)==(ikp2+ioff2) .and. ispin1==ispin2) optical_left(ispin2,isqa,ikp2+ioff2) =  dimag(Tkk_tmp(ikp1,ikp2))
                  optical_right_thread(ispin2,isqa,ikp2+ioff2) = optical_right_thread(ispin2,isqa,ikp2+ioff2) + weights(ikp1+ioff1)*Tkk_abs(ikp1,ikp2)
                end if!l_optical
              end do!ispin1
            end do!ikp1
!$omp end do 
          end do!ispin2
        end do!isqa
      end do!ikp2

      if(l_optical)then
!$omp critical
        optical_right = optical_right+optical_right_thread
!$omp end critical
      end if!l_optical

!$omp end parallel
      if(myrank==master) write(*,*) ''

      !=== Calculations ===!
      !=======  END  ======!



      !===============  BEGIN  ===============!
      !=== Gather info for optical theorem ===!
      if(l_optical)then
        ihelp = inc%ndegen*nsqa*ntot
        optical_right_all = 0d0
        call MPI_Reduce(optical_right,optical_right_all,ihelp,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
        if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Reduce(optical_right)'

        call MPI_Reduce(optical_left,optical_left_all,ihelp,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,ierr)
        if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Reduce(optical_left)'

        if(myrank==master)then
          do isqa=1,nsqa
            do ispin1=1,inc%ndegen

              if(sca%naverage==1)then
                write(filename,'("optical.",I0,".",I0)') isqa, ispin1
              else!sca%naverage==1
                write(filename,'("optical.",I0,".",I0,".",I0)') iaverage,isqa, ispin1
              end if!sca%naverage==1

              open(unit=13512,file=trim(filename),form='formatted',action='write')
              do ikp1=1,ntot
                righttmp = pi*optical_right_all(ispin1,isqa,ikp1)/BZVol
                write(13512,'(4ES18.9,F10.4)')  -optical_left_all(ispin1,isqa,ikp1), &
                                       &  righttmp,&
                                       &  optical_left_all(ispin1,isqa,ikp1) + righttmp,&
                                       &  optical_left_all(ispin1,isqa,ikp1)/righttmp,&
                                       &  (optical_left_all(ispin1,isqa,ikp1)/righttmp+1)*100d0
              end do!ikp1
              close(13512)

            end do!ispin1
          end do!isqa
        end if!myrank==master
      end if!l_optical

    end do!iaverage

  end subroutine calculate_Pkksub





  subroutine calc_lifetime(inc, lattice, impcls, Amat)
    use type_inc,       only: inc_type
    use type_data,      only: lattice_type
    use mod_calconfs,   only: get_nsqa
    use mod_symmetries, only: symmetries_type, set_symmetries, rotate_kpoints, expand_visarrays, expand_areas
    use mod_read,       only: read_kpointsfile_vis, read_kpointsfile_int, read_weights, read_fermivelocity
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_iohelp,     only: open_mpifile_setview, close_mpifile, getBZvolume
    use mod_ioformat,   only: filemode_vis, filemode_int, filename_eigvect, fmt_fn_ext, filename_lifetime, ext_vtkxml, filename_intmask, ext_mpiio
    use mod_vtkxml,     only: write_pointdata_rot
    use mod_mympi,      only: myrank, nranks, master
    use mod_mathtools,  only: pi
    use mpi
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(impcls_type),  intent(in) :: impcls
    double complex,     intent(in) :: Amat(impcls%clmso,impcls%clmso)

    integer :: nsqa
    double precision, allocatable :: taukinv(:,:)

    !for masked integration, e.g. to get surface-to-bulk scattering rates etc.
    double precision, allocatable :: taukinv_masked(:,:)
    logical :: l_maskfile
    integer, allocatable :: intmask(:)

    !symmetry arrays
    integer :: nsym
    integer, allocatable :: isym(:)
    type(symmetries_type) :: symmetries

    !local k-point arrays
    integer :: nkpts_in, nkpts_all_in, nkpts_out
    integer, allocatable :: kpt2irr_in(:), irr2kpt_in(:)
    double precision, allocatable :: kpoints_in(:,:), kpoints_out(:,:)
    double precision, allocatable :: areas_in(:), areas_out(:), weights_out(:), fermivel_in(:,:)

    !parallelization
    integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1),   &
             & recvcounts(0:nranks-1), displs(0:nranks-1), &
             & nkpt, ioff

    !eigenvector-file
    integer :: dimens(4), fh_eigv_in, fh_eigv_out
    character(len=256) :: filemode, filename
    integer :: mpistatus(MPI_STATUS_SIZE)
    logical :: l_exist ! logical to check if eigenvector file exists
    double complex, allocatable :: rveig_out(:,:,:,:,:),&
                                 & rveig_in(:,:,:,:),   &
                                 & rveig_out_impcls(:), &
                                 & rveig_in_impcls(:,:,:)

    !temp k-point arrays
    integer :: nkpts1, nkpts2, nkpts_all1
    double precision :: fac, BZVol
    double complex :: Ctmp(impcls%clmso), Tkk_tmp
    integer,          allocatable :: kpt2irr1(:), irr2kpt1(:)
    double precision, allocatable :: areas1(:), weights1(:), kpoints1(:,:), fermivel1(:,:)

    !visualization
    integer                       :: nscalar, nvector, iscalar, ivector
    double precision, allocatable :: scalardata(:,:), &
                                   & vectordata(:,:,:)
    character(len=256), allocatable :: scalarstring(:), vectorstring(:)

    !loop counter and temp
    integer :: ierr, irveigel, ikpi, ikpo, ispini, ispino, isqa, istore, ihelp, clmso, printstep
    integer :: itmp1(1)
    character(len=64) :: fsRyunit

    !parameter
    double precision, parameter :: RyToinvfs = 20.67068667282055d0
    double complex, parameter :: CZERO=(0d0,0d0), CONE=(1d0,0d0), CI=(0d0,1d0)
    character(len=1) :: cspinlabel(2) = (/ 'u','d' /)

    if(.not.sca_read) then
      call read_sca()
    end if
    call set_symmetries(inc, lattice, symmetries)
    nsqa = get_nsqa()
    BZVol = getBZvolume(lattice%recbv)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+                 2 pi c*N      1                                   +
    !+    Tau^(-1) = ------------ * ---- * int_FS{ |Tkk'|^2  / v_F }dS   +
    !+                   hbar       V_BZ                                 +
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   fac = 2d0*pi*0.01d0/BZVol
!   if(myrank==master) write(*,*) 'Output in Rydbergs!'
!   fsRyunit = '(Ry^{-1}/at%)'
    fac = 2d0*pi*0.01d0/BZVol*RyToinvfs
    fsRyunit = '(fs/at%)'
    if(myrank==master) write(*,*) 'Output in femtoseconds!'



    !====================================================!
    !=== read in the (irreducible) k-points and apply ===!
    !===   symmetry operations to expand to full BZ   ===!
    !====================================================!

    !==========================!
    != for incomming k-vector =!
    select case (sca%mode)

    case (MODE_VIS)
      filemode = filemode_vis
      call read_kpointsfile_vis(nkpts_all1, nkpts1, kpoints1, nsym, isym, kpt2irr1, irr2kpt1)
      call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts_in, kpoints_in)
      call expand_visarrays(nsym, nkpts_all1, nkpts1, kpt2irr1, irr2kpt1, kpt2irr_in, irr2kpt_in)
      nkpts_all_in = nsym*nkpts_all1
      deallocate(kpt2irr1, irr2kpt1, kpoints1, isym)
      call read_fermivelocity(filemode_vis, nkpts1, fermivel1, nsym, isym)
      call rotate_kpoints(symmetries%rotmat, nkpts1, fermivel1, nsym, isym, nkpts2, fermivel_in)
      if(nkpts2/=nkpts_in) stop 'inconsistency in number of k-points'
      deallocate(fermivel1, isym)

    case (MODE_INT)
      filemode = filemode_int
      call read_kpointsfile_int(nkpts1, kpoints1, areas1, nsym, isym)
      call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts_in, kpoints_in)
      call expand_areas(nsym,nkpts1,areas1,areas_in)
      deallocate(areas1, kpoints1, isym)
!     call read_fermivelocity(filemode_int, nkpts1, fermivel1, nsym, isym)
!     call rotate_kpoints(symmetries%rotmat, nkpts1, fermivel1, nsym, isym, nkpts2, fermivel_in)
!     if(nkpts2/=nkpts_in) stop 'inconsistency in number of k-points'
!     deallocate(fermivel1, isym)

    case default; stop 'case not known in select case (cfg%mode)'
    end select
    != for incomming k-vector =!
    !==========================!


    !=========================!
    != for outgoing k-vector =!
    call read_kpointsfile_int(nkpts1, kpoints1, areas1, nsym, isym)
    call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts_out, kpoints_out)
    call expand_areas(nsym,nkpts1,areas1,areas_out)
    deallocate(areas1, kpoints1, isym)
    call read_weights(nkpts1, weights1, nsym, isym)
    call expand_areas(nsym,nkpts1,weights1,weights_out)
    deallocate(isym)
    != for outgoing k-vector =!
    !=========================!

    !redefine symmetries (set to unit transformation only)
    nsym = 1
    allocate(isym(nsym))
    isym = (/ 1 /)

    !allocate eigenvector files
    allocate( rveig_out(inc%lmmaxso,inc%natypd,inc%ndegen,nsqa,nkpts_out),&
            & rveig_in(inc%lmmaxso,inc%natypd,inc%ndegen,nsqa), &
            & rveig_out_impcls(impcls%clmso),                   &
            & rveig_in_impcls(impcls%clmso,inc%ndegen,nsqa),    &
            & STAT=ierr                                         )
    if(ierr/=0) stop 'Problem allocating rveig_out etc'

    !read in outgoing eigenvectors
    write(filename,'(A,A)') filename_eigvect, filemode_int
    dimens = (/ inc%lmmaxso,inc%natypd,inc%ndegen, nsqa /)
    itmp1(1) = nkpts_out
    inquire(file=trim(filename)//ext_mpiio, exist=l_exist)
    if(.not.l_exist) stop 'Error: eigenvector file (int) not present but needed for scattering calculations'
    call open_mpifile_setview( trim(filename), 'read', 4, dimens, itmp1, MPI_DOUBLE_COMPLEX, &
                             & 0, 1 , MPI_COMM_WORLD, fh_eigv_out                            )
    ihelp = product(dimens)*nkpts_out
    call MPI_File_read_all(fh_eigv_out, rveig_out, ihelp, MPI_DOUBLE_COMPLEX, mpistatus, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Problem reading outgoing eigenvectors in calc_lifetime'
    call close_mpifile(fh_eigv_out)

    !parallelize over incomming k-vectors
    call distribute_linear_on_tasks(nranks, myrank, master, nkpts_in, ntot_pT, ioff_pT, .true.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    write(filename,'(A,A)') filename_eigvect, trim(filemode)
    dimens = (/ inc%lmmaxso,inc%natypd,inc%ndegen, nsqa /)
    inquire(file=trim(filename)//ext_mpiio, exist=l_exist)
    if(.not.l_exist) stop 'Error: eigenvector file (vis) not present but needed for scattering calculations'
    call open_mpifile_setview( trim(filename), 'read', 4, dimens, ntot_pT, MPI_DOUBLE_COMPLEX, &
                             & myrank, nranks, MPI_COMM_WORLD, fh_eigv_in                      )



    !==============================!
    !=== Calculate the Lifetime ===!
    !==============================!
    irveigel = product(dimens)
    clmso = impcls%clmso
    allocate( taukinv(inc%ndegen*inc%ndegen*nsqa,nkpts_in), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating taukinv'
    taukinv  = 0d0

    !masked integration:
    if(sca%maskint==1) then
       ! allocate second array
       allocate( taukinv_masked(inc%ndegen*inc%ndegen*nsqa,nkpts_in), STAT=ierr )
       if(ierr/=0) stop 'Problem allocating taukinv_masked'
       taukinv_masked  = 0d0
       ! allocate and read in integration mask from file
       allocate( intmask(nkpts_out), STAT=ierr )
       if(ierr/=0) stop 'Problem allocating intmask'
       intmask  = 0
       if(myrank==master) then
          ! check if file exists
          inquire(file=filename_intmask, exist=l_maskfile)
          if(.not.l_maskfile) then
             write(*,'(3A)') 'Error: file "',filename_intmask,'" not found! Please provide file for masked integration.'
             stop
          endif
          ! read file
          open(998877, file=filename_intmask, form='formatted')
          do istore=1,nkpts_out
             read(998877, *) intmask(istore)
          end do
          close(998877)
       endif!myrank==master
       !Bcast intmask to all processors
       call MPI_Bcast(intmask, nkpts_out, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
       if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Bcast for intmask'
    endif

    !calculate Pkk' and sum it up
    printstep = nkpt/50
    if(printstep==0) printstep=1
    !print header of statusbar
    if(myrank==master) write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
    if(myrank==master) write(*,FMT=190) !beginning of statusbar
    do ikpi=1,nkpt
      !update statusbar
      if(mod(ikpi,printstep)==0 .and. myrank==master) write(*,FMT=200)

      !read incomming eigenvector
      call MPI_File_read( fh_eigv_in, rveig_in, irveigel, MPI_DOUBLE_COMPLEX, mpistatus, ierr )
      if(ierr/=MPI_SUCCESS) stop 'error reading fh_eigv_in'
      call extend_coeffvector2cluster_allsqa( inc, lattice, impcls, nsqa, &
                                            & kpoints_in(:,ioff+ikpi), rveig_in, rveig_in_impcls )

      do isqa=1,nsqa
        do ispini=1,inc%ndegen

          ! Ctmp = Amat * c_k
          call ZGEMM( 'N','N', clmso, 1, clmso, CONE, Amat, clmso,              &
                    & rveig_in_impcls(:,ispini,isqa), clmso, CZERO, Ctmp, clmso )

          do ikpo=1,nkpts_out
            do ispino=1,inc%ndegen

              call extend_coeffvector2cluster( inc, lattice, impcls, kpoints_out(:,ikpo),        &
                                             & rveig_out(:,:,ispino,isqa,ikpo), rveig_out_impcls )

              istore = (isqa-1)*inc%ndegen**2 + (ispini-1)*inc%ndegen + ispino
              ! Tkk' = conjg(c_k') * Amat * c_k
              call ZGEMM( 'C','N', 1, 1, clmso, CONE, rveig_out_impcls, &
                        & clmso, Ctmp, clmso, CZERO, Tkk_tmp, 1         )

              ! Pkk' = |Tkk'|^2
              taukinv(istore,ikpi) = taukinv(istore,ikpi) + (dble(Tkk_tmp)**2 + aimag(Tkk_tmp)**2)*weights_out(ikpo)

              !masked integration:
              if(sca%maskint==1) then
                 ! integrate only those kpts where intmask==1
                 if(intmask(ikpo)==1) taukinv_masked(istore,ikpi) = taukinv_masked(istore,ikpi) + (dble(Tkk_tmp)**2 + aimag(Tkk_tmp)**2)*weights_out(ikpo)
              endif

            end do!ispino
          end do!ikpo

        end do!ispini
      end do!isqa
    end do!ikpi
    if(myrank==master) write(*,*) ''
    !==============================!
    !=== Calculate the Lifetime ===!
    !==============================!



    !=======================!
    !=== Collect results ===!
    !=======================!
    ihelp      = inc%ndegen*inc%ndegen*nsqa
    recvcounts = ntot_pT*ihelp
    displs     = ioff_pT*ihelp
    call MPI_Allgatherv( taukinv, ihelp*nkpt, MPI_DOUBLE_PRECISION,         &
                       & taukinv, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                       & MPI_COMM_WORLD, ierr )

    if(sca%maskint==1) then
       call MPI_Allgatherv( taukinv_masked, ihelp*nkpt, MPI_DOUBLE_PRECISION,         &
                          & taukinv_masked, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                          & MPI_COMM_WORLD, ierr )
    endif



    !==========================!
    !=== Calculate averages ===!
    !==========================!
    if(myrank==master) then
      select case(sca%mode)
        case(MODE_VIS); call calculate_lifetimeaverage_vis(inc,nsqa,nkpts_in,nkpts_all_in,kpt2irr_in,kpoints_in,fermivel_in,taukinv,.true.,fac,fsRyunit)
        case(MODE_INT); call calculate_lifetimeaverage_int(inc,nsqa,nkpts_in,weights_out,taukinv,.true.,fac,fsRyunit)
        case default; stop 'sca%mode not known in calculating lifetimeaverages'
      end select
    end if

    !======================!
    !=== Save results   ===!
    !======================!
    if(myrank==master) call save_lifetime(trim(filemode), nkpts_in, nsqa, inc%ndegen, taukinv, nsym, isym, fac, fsRyunit)
    if(sca%maskint==1) then
       write(filename,'(A,A)') trim(filemode), '.masked'
       if(myrank==master) call save_lifetime(trim(filename), nkpts_in, nsqa, inc%ndegen, taukinv_masked, nsym, isym, fac, fsRyunit)
    endif

    !======================!
    !=== Visualize data ===!
    !======================!
    if(myrank==master .and. sca%mode==MODE_VIS .and. inc%nBZdim==3)then


      if(inc%ndegen==2)then

        nscalar = 2*nsqa
        nvector=0
        allocate(scalardata(nkpts_in,nscalar), scalarstring(nscalar), STAT=ierr)
        if(ierr/=0) stop 'Problem allocating scalardata'

        iscalar = 0
        do isqa=1,nsqa
          iscalar = iscalar+1
          scalardata(:,iscalar) = 2d0/(taukinv(1+4*(isqa-1),:) + taukinv(4+4*(isqa-1),:))/fac
          write(scalarstring(iscalar),'(A,I0,A,A)') 'sqa=',isqa,'_tau',trim(fsRyunit)

          iscalar = iscalar+1
          scalardata(:,iscalar) = 1d0/(taukinv(2+4*(isqa-1),:) + taukinv(3+4*(isqa-1),:))/fac
          write(scalarstring(iscalar),'(A,I0,A,A)') 'sqa=',isqa,'_T1',trim(fsRyunit)
        end do!isqa


      else!inc%ndegen==2


        nscalar=nsqa
        nvector=0
        allocate(scalardata(nkpts_in,nscalar), scalarstring(nscalar), STAT=ierr)
        if(ierr/=0) stop 'Problem allocating scalardata'

        iscalar = 0
        do isqa=1,nsqa
          iscalar = iscalar+1
          scalardata(:,iscalar) = 1d0/(taukinv(1+4*(isqa-1),:))/fac 
          write(scalarstring(iscalar),'(A,I0,A,A)') 'sqa=',isqa,'_tau',trim(fsRyunit)
        end do!isqa


      end if!inc%ndegen==2


      write(filename,fmt_fn_ext) filename_lifetime, ext_vtkxml
      call write_pointdata_rot( trim(filename),nkpts_in,kpoints_in,&
                              & nscalar,scalardata,scalarstring,   &
                              & nvector,vectordata,vectorstring,   &
                              & nsym,symmetries%rotmat,isym,       &
                              & nkpts_all_in, kpt2irr_in           )

    end if!sca%mode==MODE_VIS

190 FORMAT('                 |'$)
200 FORMAT('|'$)
  end subroutine calc_lifetime




  subroutine calc_scattering_fixk(inc, lattice, cluster, tgmatrx, impcls, Amat)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_calconfs,   only: get_nsqa
    use mod_symmetries, only: symmetries_type, set_symmetries, rotate_kpoints, expand_visarrays, expand_areas
    use mod_read,       only: read_kpointsfile_vis, read_kpointsfile_int
    use mod_ioformat,   only: filemode_vis, filemode_int, filename_eigvect, fmt_fn_ext, filename_scattfixk, ext_vtkxml
    use mod_mympi,      only: myrank, nranks, master
    use mod_parutils,   only: distribute_linear_on_tasks
    use mod_iohelp,     only: open_mpifile_setview, close_mpifile
    use mod_vtkxml,     only: write_pointdata_rot
    use mpi
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx
    type(impcls_type),  intent(in) :: impcls
    double complex,     intent(in) :: Amat(impcls%clmso,impcls%clmso)

    type(symmetries_type) :: symmetries
    integer :: nsqa

    !arrays for the incomming k-vector
    double precision :: kpoint_fix(3)
    double precision, allocatable :: spinval_fix(:,:)
    double complex,   allocatable :: rveig_fix(:,:,:,:), rveig_fix_impcls(:,:,:), rveig_in(:,:,:,:), rveig_in_impcls(:,:,:)

    !arrays for the outgoing k-vector
    integer :: nkpts, nkpts_all, nsym
    integer, allocatable :: isym(:), kpt2irr(:), irr2kpt(:)
    double precision, allocatable :: kpoints(:,:)
    double precision, allocatable :: areas(:)

    !scattering
    double precision, allocatable :: Pkkfix(:,:,:,:), Pkkfix_outfix(:,:,:,:)!, Pkkfix_loc(:,:,:,:)
    double precision, allocatable :: tmp_global(:,:,:,:) ! for communication

    !parallelization
    integer :: ntot_pT(0:nranks-1),    ioff_pT(0:nranks-1), &
             & recvcounts(0:nranks-1), displs(0:nranks-1),  &
             & nkpt, ioff

    !eigenvector-file
    integer :: dimens(4), fh_eigv
    character(len=256) :: filemode, filename
    integer :: mpistatus(MPI_STATUS_SIZE)

    !visualization
    integer                       :: nscalar, nvector, iscalar, ivector
    double precision, allocatable :: scalardata(:,:), &
                                   & vectordata(:,:,:)
    character(len=256), allocatable :: scalarstring(:), vectorstring(:)

    !counter and temp arrays
    integer :: ierr, isqa, ikp, ispin1, ispin2, clmso, ihelp, irveigel, printstep
    integer :: nkpts1, nkpts_all1
    double complex :: Tkk_tmp,Tkk_tmp_outfix
    integer,          allocatable :: kpt2irr1(:), irr2kpt1(:)
    double precision, allocatable :: areas1(:), kpoints1(:,:)
    double complex,   allocatable :: Ctmp(:,:,:),Ctmp_outfix(:,:,:)

    !parameter
    double complex, parameter :: CZERO=(0d0,0d0), CONE=(1d0,0d0), CI=(0d0,1d0)
    character(len=1) :: cspinlabel(2) = (/ 'u','d' /)

    if(.not.sca_read) then
      call read_sca()
    end if
    call set_symmetries(inc, lattice, symmetries)

    if(myrank==master) write(*,*) 'get_incident_kvector'
    call get_incident_kvector(inc, lattice, cluster, tgmatrx, kpoint_fix)

    nsqa = get_nsqa()
    allocate( rveig_fix(inc%lmmaxso,inc%natypd,inc%ndegen,nsqa),&
            & rveig_in(inc%lmmaxso,inc%natypd,inc%ndegen,nsqa), &
            & rveig_fix_impcls(impcls%clmso,inc%ndegen,nsqa),   &
            & rveig_in_impcls(impcls%clmso,inc%ndegen,nsqa),    &
            & spinval_fix(inc%ndegen,nsqa),                     &
            & STAT=ierr                                         )
    if(ierr/=0) stop 'Problem allocating rveig_fix etc'

    if(myrank==master) write(*,*) 'calc_wavefun_kpoint'
    call calc_wavefun_kpoint(inc, lattice, cluster, tgmatrx, nsqa, kpoint_fix, spinval_fix, rveig_fix)
    if(myrank==master) write(*,*) 'extend_coeffvector2cluster_allsqa'
    call extend_coeffvector2cluster_allsqa(inc, lattice, impcls, nsqa, kpoint_fix, rveig_fix, rveig_fix_impcls)


    !read in the (irreducible) k-points and apply symmetry operations to expand to full BZ
    select case (sca%mode)
    case (MODE_VIS)
      filemode = filemode_vis
      call read_kpointsfile_vis(nkpts_all1, nkpts1, kpoints1, nsym, isym, kpt2irr1, irr2kpt1)
      call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
      call expand_visarrays(nsym, nkpts_all1, nkpts1, kpt2irr1, irr2kpt1, kpt2irr, irr2kpt)
      nkpts_all = nsym*nkpts_all1
      deallocate(kpt2irr1, irr2kpt1)
    case (MODE_INT)
      filemode = filemode_int
      call read_kpointsfile_int(nkpts1, kpoints1, areas1, nsym, isym)
      call rotate_kpoints(symmetries%rotmat, nkpts1, kpoints1, nsym, isym, nkpts, kpoints)
      call expand_areas(nsym,nkpts1,areas1,areas)
      deallocate(areas1)
    case default; stop 'case not known in select case (cfg%mode)'
    end select

    !redefine symmetries (set to unit transformation only)
    nsym = 1
    deallocate(isym, kpoints1)
    allocate(isym(nsym))
    isym = (/ 1 /)

    !allocate(kpoints(3,2))
    !kpoints(:,1) = kpoint_fix(:)
    !kpoints(:,2) = -kpoint_fix(:)
    !if(myrank==master) write(*,*) 'kpoint fix, kpoint',kpoint_fix,kpoints

    !nkpts = 2
    call distribute_linear_on_tasks(nranks, myrank, master, nkpts, ntot_pT, ioff_pT, .true.)
    nkpt = ntot_pT(myrank)
    ioff = ioff_pT(myrank)

    write(filename,'(A,A)') filename_eigvect, trim(filemode)
    dimens = (/ inc%lmmaxso,inc%natypd,inc%ndegen, nsqa /)
    irveigel = product(dimens)
    call open_mpifile_setview( trim(filename), 'read', 4, dimens, ntot_pT, MPI_DOUBLE_COMPLEX, &
                             & myrank, nranks, MPI_COMM_WORLD, fh_eigv                         )

    clmso = impcls%clmso
    allocate( Pkkfix(inc%ndegen,inc%ndegen,nsqa,nkpts), &
            & Pkkfix_outfix(inc%ndegen,inc%ndegen,nsqa,nkpts), &
            & Ctmp(clmso,inc%ndegen,nsqa), &
            & Ctmp_outfix(clmso,inc%ndegen,nsqa), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating Pkkfix, Ctmp in calc_scattering_fixk'
    Ctmp = CZERO

    !================ BEGIN =============!
    !=== Prepare the incomming vector ===!
    do isqa=1,nsqa
      do ispin1=1,inc%ndegen
!OLD--------------------------------------------------------------------------------------------!OLD
!OLD    ! Ctmp = conjg(c_kfix) * Amat                                                           !OLD
! uncommented for k0==kfinal writeout as well
        call ZGEMM( 'C','N', 1, clmso, clmso, CONE, rveig_fix_impcls(:,ispin1,isqa), &          !OLD
                  & clmso, Amat, clmso, CZERO, Ctmp_outfix(:,ispin1,isqa), 1                )          !OLD
!OLD--------------------------------------------------------------------------------------------!OLD
!NEW    in the old implementation, the fixed k was actually an outgoing k-vector                !NEW
!NEW    ! Ctmp =  Amat*c_kfix                                                                   !NEW
!NEW--------------------------------------------------------------------------------------------!NEW
        call ZGEMM( 'N','N', clmso, 1, clmso, CONE, Amat, clmso,                              & !NEW
                  & rveig_fix_impcls(:,ispin1,isqa), clmso, CZERO, Ctmp(:,ispin1,isqa), clmso ) !NEW
      end do!ispin1
    end do!isqa
    !=== Prepare the incomming vector ===!
    !================  END  =============!


    !================ BEGIN ================!
    !=== Calculate the scattering matrix ===!
    !===      for fixed incomming k      ===!
    printstep = nkpt/50
    if(printstep==0) printstep=1
    !print header of statusbar
    if(myrank==master) write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
    if(myrank==master) write(*,FMT=190) !beginning of statusbar
    do ikp=1,nkpt
      !update statusbar
      if(mod(ikp,printstep)==0 .and. myrank==master) write(*,FMT=200)

      call MPI_File_read(fh_eigv, rveig_in, irveigel, MPI_DOUBLE_COMPLEX, mpistatus, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error reading fh_eigv'
      call extend_coeffvector2cluster_allsqa(inc, lattice, impcls, nsqa, kpoints(:,ioff+ikp), rveig_in, rveig_in_impcls)

      do isqa=1,nsqa
        do ispin1=1,inc%ndegen
          do ispin2=1,inc%ndegen

!OLD---------------------------------------------------------------------------------------------!OLD
!OLD        ! Tkk' = conjg(c_kfix) * Amat * c_k                                                  !OLD
! uncommented for k0==kfinal writeout as well
            call ZGEMM( 'N','N', 1, 1, clmso, CONE, Ctmp_outfix(:,ispin1,isqa), 1,      &               !OLD
                      & rveig_in_impcls(:,ispin2,isqa), clmso, CZERO, Tkk_tmp_outfix, 1 )               !OLD
!OLD---------------------------------------------------------------------------------------------!OLD
!NEW         in the old implementation, the fixed k was actually an outgoing k-vector         ---!NEW
!NEW        ! Tkk' = conjg(c_k) * Amat * c_kfix                                                  !NEW
            call ZGEMM( 'C','N', 1, 1, clmso, CONE, rveig_in_impcls(:,ispin2,isqa), clmso,&      !NEW
                      & Ctmp(:,ispin1,isqa), clmso, CZERO, Tkk_tmp, 1                     )      !NEW

            ! Pkk = |Tkk|^2
            Pkkfix(ispin2,ispin1,isqa,ikp) = dble(Tkk_tmp)**2 + aimag(Tkk_tmp)**2
            Pkkfix_outfix(ispin2,ispin1,isqa,ikp) = dble(Tkk_tmp_outfix)**2 + aimag(Tkk_tmp_outfix)**2

          end do!ispin2
        end do!ispin1
      end do!isqa
    end do!ikp
    if(myrank==master) write(*,*) ''
    !=== Calculate the scattering matrix ===!
    !===      for fixed incomming k      ===!
    !================  END  ================! 

    call close_mpifile(fh_eigv)


    ihelp      = inc%ndegen*inc%ndegen*nsqa
    recvcounts = ntot_pT*ihelp
    displs     = ioff_pT*ihelp

    if(myrank==master) write(*,*) 'before mpiallgather Pkkfix:',ihelp,recvcounts,displs,nkpt

    allocate( tmp_global(inc%ndegen,inc%ndegen,nsqa,nkpts), stat=ierr)
    if(ierr/=0) stop 'Error allocating tmp_global for Pkkfix communication in mod_scattering'

    tmp_global(:,:,:,:) = 0.0d0
    call MPI_Allgatherv( Pkkfix, ihelp*nkpt, MPI_DOUBLE_PRECISION,         &
                       & tmp_global, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                       & MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allgatherv for PKKfix'
    Pkkfix = tmp_global

! gather Pkk' for fixed outgoing k
    tmp_global(:,:,:,:) = 0.0d0
    call MPI_Allgatherv( Pkkfix_outfix, ihelp*nkpt, MPI_DOUBLE_PRECISION,         &
                       & tmp_global, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                       & MPI_COMM_WORLD, ierr )
    if(ierr/=MPI_SUCCESS) stop 'Error in MPI_Allgatherv for PKKfix_outfix'
    Pkkfix_outfix = tmp_global

    deallocate( tmp_global, stat=ierr)
    if(ierr/=0) stop 'Error deallocating tmp_global for Pkkfix communication in mod_scattering'


    !save scattfix to file
    if(myrank==master) call save_scattfix(trim(filemode), nkpts, nsqa, inc%ndegen, kpoint_fix, Pkkfix, nsym, isym, 'start')
    if(myrank==master) call save_scattfix(trim(filemode), nkpts, nsqa, inc%ndegen, kpoint_fix, Pkkfix_outfix, nsym, isym, 'outfix')


    if(myrank==master .and. sca%mode==MODE_VIS .and. inc%nBZdim==3) then
      nscalar = inc%ndegen**2 * nsqa
      allocate(scalardata(nkpts,nscalar), scalarstring(nscalar), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating scalardata'

      iscalar = 0
      nvector=0
      do isqa=1,nsqa
        do ispin1=1,inc%ndegen
          do ispin2=1,inc%ndegen
            iscalar = iscalar+1
            scalardata(:,iscalar) = Pkkfix(ispin2,ispin1,isqa,:)
            write(scalarstring(iscalar),'(A,I0,A,A,A)') 'sqa=',isqa,'_',cspinlabel(ispin2),cspinlabel(ispin1)
          end do!ispin2
        end do!ispin1
      end do!isqa


      write(filename,fmt_fn_ext) filename_scattfixk, ext_vtkxml
      call write_pointdata_rot( trim(filename),nkpts,kpoints,   &
                              & nscalar,scalardata,scalarstring,&
                              & nvector,vectordata,vectorstring,&
                              & nsym,symmetries%rotmat,isym,    &
                              & nkpts_all, kpt2irr              )

    end if!sca%mode==MODE_VIS


190 FORMAT('                 |'$)
200 FORMAT('|'$)
  end subroutine calc_scattering_fixk





  subroutine get_incident_kvector(inc, lattice, cluster, tgmatrx, kpoint_fix)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_fermisurf_basic,  only: roots_along_edge, ROOT_IMAG
    use mod_mympi,      only: myrank, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    type(inc_type),     intent(in)  :: inc
    type(lattice_type), intent(in)  :: lattice
    type(cluster_type), intent(in)  :: cluster
    type(tgmatrx_type), intent(in)  :: tgmatrx
    double precision,   intent(out) :: kpoint_fix(3)

    integer :: nroots, iroot, nsqa, ierr
    integer, allocatable :: lmroot(:)
    double precision, allocatable :: kroot(:,:)
    double complex, allocatable :: LVroot(:,:,:), &
                                 & RVroot(:,:,:), &
                                 & eigwroot(:,:)

   double precision, parameter :: eps=1e-6!threshold which determines if kstart==kstop so that line-search is not done

    if(myrank==master)then

      if(sqrt((sca%kfix(1,1)-sca%kfix(1,2))**2+(sca%kfix(2,1)-sca%kfix(2,2))**2+(sca%kfix(3,1)-sca%kfix(3,2))**2)<eps) then
          write(*,*) 'take input as k-point'
          ! if start==stop for fixed interval take only this point without finding the intersection with the FS
          kpoint_fix = sca%kfix(:,1)
      else

      allocate( lmroot(inc%nrootmax),&
              & kroot(3,inc%nrootmax),&
              & LVroot(inc%almso,inc%almso,inc%nrootmax),&
              & RVroot(inc%almso,inc%almso,inc%nrootmax),&
              & eigwroot(inc%almso,inc%nrootmax),&
              & STAT=ierr)
      if(ierr/=0) stop 'Problem allocating lmroot etc. in calc_scattering_fix'

      !find the FS-point on a given path in the BZ
      call roots_along_edge( inc, lattice, cluster, tgmatrx, sca%nsteps, sca%kfix, sca%niter, ROOT_IMAG, &
                           & sca%rooteps, nroots, lmroot, kroot, LVroot, RVroot, eigwroot, -1            )

      write(*,'(A,I0)') 'Number of roots found on path: ', nroots
      do iroot=1,nroots
        write(*,'(4X,A,I0,A,3ES25.16,A,I0,A,2ES25.16,A)') 'root #', iroot,&
                                                   & ' kvect=[ ', kroot(:,iroot),&
                                                   & ' ],  LM= ', lmroot(iroot),&
                                                   & ', eigw=(', eigwroot(lmroot(iroot),iroot), ')'
      end do!iroot
      write(*,'(A,I0)') 'Take k-point #', sca%roottake
      kpoint_fix = kroot(:,sca%roottake)

      deallocate(lmroot, LVroot, RVroot, eigwroot)

      end if !kfix_start==kfix_stop

    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(kpoint_fix, 3, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting kpoint_fix'
#endif

  end subroutine get_incident_kvector





  subroutine calc_wavefun_kpoint(inc, lattice, cluster, tgmatrx, nsqa, kpoint, spinval, eigvect_rot)

    use type_inc,       only: inc_type
    use type_data,      only: lattice_type, cluster_type, tgmatrx_type
    use mod_calconfs,   only: calc_spinvalue_state
    use mod_kkrmat,     only: compute_kkr_eigenvectors2
    use mod_mathtools,  only: findminindex, bubblesort
    use mod_mympi,      only: myrank, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    type(tgmatrx_type), intent(in) :: tgmatrx
    integer,            intent(in) :: nsqa
    double precision,   intent(in) :: kpoint(3)
    double precision,   intent(out) :: spinval(inc%ndegen,nsqa)
    double complex,     intent(out) :: eigvect_rot(inc%lmmaxso,inc%natypd,inc%ndegen,nsqa)


    integer        :: lm_fs, lm_fs2, ihelp, ierr
    double complex :: eigwEF(inc%almso),&
                    & LVeigEF(inc%almso,inc%almso),&
                    & RVeigEF(inc%almso,inc%almso),&
                    & rveig(inc%almso,inc%ndegen)

    integer          :: indsort(inc%almso)
    double precision :: deigsort(inc%almso)

    if(myrank==master)then
      !===============================================================!
      !=== find the eigenvector and eigenvalue at the fermi energy ===!
      !===============================================================!
      call compute_kkr_eigenvectors2( inc, lattice, cluster, tgmatrx%ginp(:,:,:,1),&
                                    & tgmatrx%tmat(:,:,1), kpoint,                 &
                                    & eigwEF, LVeigEF, RVeigEF                     )
      deigsort = abs(eigwEF)
      call findminindex(inc%almso, deigsort, lm_fs)

      !copy the band
      rveig(:,1) = RVeigEF(:,lm_fs)

      !find degenerate band for degnerate cases
      if(inc%ndegen==2)then
        deigsort = abs(eigwEF-eigwEF(lm_fs))
        call bubblesort(inc%almso, deigsort, indsort)
        lm_fs2 = indsort(2)
        rveig(:,2) = RVeigEF(:,lm_fs2)
      end if

      call calc_spinvalue_state( inc, tgmatrx, rveig, spinval, eigvect_rot )

    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(spinval, inc%ndegen*nsqa, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting spinval in calc_wavefun_kpoint'

    ihelp = inc%lmmaxso*inc%natypd*inc%ndegen*nsqa
    call MPI_Bcast(eigvect_rot, ihelp, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting spinval in calc_wavefun_kpoint'
#endif

  end subroutine calc_wavefun_kpoint





  subroutine read_sca()

    use mod_ioinput, only: IoInput
    use mod_mympi,   only: myrank, nranks, master
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer           :: ierr
    integer           :: conducti ! used to ensure compatibility of old format input files
    double precision  :: dtmpin(3), theta, phi
    character(len=80) :: uio

    if(myrank==master) then


      call IoInput('SCATTFIX  ',uio,1,7,ierr)
      read(unit=uio,fmt=*) sca%lscatfixk

      if(sca%lscatfixk==1) then
        call IoInput('SCAMODE   ',uio,1,7,ierr)
        read(unit=uio,fmt=*) sca%mode
        call IoInput('KFIXSTART ',uio,1,7,ierr)
        read(unit=uio,fmt=*) sca%kfix(:,1)
        call IoInput('KFIXSTOP  ',uio,1,7,ierr)
        read(unit=uio,fmt=*) sca%kfix(:,2)
        call IoInput('SXNSTEPS  ',uio,1,7,ierr)
        read(unit=uio,fmt=*) sca%nsteps
        call IoInput('SXNITER   ',uio,1,7,ierr)
        read(unit=uio,fmt=*) sca%niter
        call IoInput('SXROOTNUM ',uio,1,7,ierr)
        read(unit=uio,fmt=*) sca%roottake
        call IoInput('SXROOTEPS ',uio,1,7,ierr)
        read(unit=uio,fmt=*) sca%rooteps
      end if!sca%lscatfixk==1



      call IoInput('LLIFETIME ',uio,1,7,ierr)
      read(unit=uio,fmt=*) sca%llifetime

      if(sca%llifetime==1) then
        call IoInput('SCAMODE   ',uio,1,7,ierr)
        read(unit=uio,fmt=*) sca%mode
      end if!sca%llifetime==1



      call IoInput('LBOLTZ    ',uio,1,7,ierr)
      if(ierr==0) then ! new input file format
        read(unit=uio,fmt=*) sca%lboltzmann

        if(sca%lboltzmann==1) then
          call IoInput('SAVEPKK   ',uio,1,7,ierr)
          read(unit=uio,fmt=*) sca%savepkk

          call IoInput('PROCMAP   ',uio,1,7,ierr)
          read(unit=uio,fmt=*) sca%subarr_inp

          call IoInput('NAVERAGE  ',uio,1,7,ierr)
          read(unit=uio,fmt=*) sca%naverage

          allocate(sca%weight_imp(sca%naverage), STAT=ierr)
          if(ierr/=0) stop 'Problem allocating weight_imp'

          call IoInput('WEIGHTIMP ',uio,1,7,ierr)
          if(ierr==0) then
            read(unit=uio,fmt=*) sca%weight_imp
          else ! ensures compatibility of old format input files
            write(*,*) "Warning : WEIGHTIMP set to 1/naverage by default !"
            sca%weight_imp(:)=1d0/sca%naverage
          end if!ierr/=0

          call IoInput('GAMMAMODE ',uio,1,7,ierr)
          read(unit=uio,fmt=*) sca%gammamode
          if(sca%gammamode==1.or.sca%gammamode==2)then
            call IoInput('GAMMAVAL  ',uio,1,7,ierr)
            read(unit=uio,fmt=*) sca%gammaval
            call IoInput('IMPCONC   ',uio,1,7,ierr)
            read(unit=uio,fmt=*) sca%impconc
          end if!sca%gammamode==1.or.sca%gammamode==2
        end if!sca%lboltzmann==1

      else ! old input file format
        call IoInput('LCONDUCTI ',uio,1,7,ierr)
        if(ierr/=0) stop 'Either LCONDUCTI or LBOLTZ should be present in input file'
        read(unit=uio,fmt=*) conducti

        if(conducti==1) then
          call IoInput('SAVEPKK   ',uio,1,7,ierr)
          read(unit=uio,fmt=*) sca%savepkk
          call IoInput('PROCMAP   ',uio,1,7,ierr)
          read(unit=uio,fmt=*) sca%subarr_inp
          call IoInput('NAVERAGE  ',uio,1,7,ierr)
          read(unit=uio,fmt=*) sca%naverage

          ! uses default values for keywords that are absent in old format input file
          write(*,*) "Warning : the old format LCONDUCTI keyword has been found !"
          write(*,*) "The following new format keywords will be assigned default values:"
          write(*,*) "LBOLTZ    = 1"
          write(*,*) "GAMMAMODE = 0"
          write(*,*) "GAMMAVAL  = 0"
          write(*,*) "IMPCONC   = 1"
          sca%lboltzmann = 1
          sca%gammamode  = 0
          sca%gammaval   = 0
          sca%impconc    = 1
        endif
      end if!ierr


      if(sca%llifetime==1) then
         call IoInput('MASKINT   ',uio,1,7,ierr)
         if(ierr==0) then
            read(unit=uio,fmt=*) sca%maskint
         else
            sca%maskint = 0
            write(unit=*,fmt=*) 'Warning: MASKINT not found in inputcard. Take default value:', sca%maskint
         end if
      end if

    end if!myrank==master

#ifdef CPP_MPI
    call myBcast_sca(sca)
#endif

    sca_read = .true.

  end subroutine read_sca





#ifdef CPP_MPI
  subroutine myBcast_sca(sca)

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(sca_type), intent(inout)  :: sca

    integer :: blocklen1(sca%N1), etype1(sca%N1), myMPItype1
    integer :: blocklen2(sca%N2), etype2(sca%N2), myMPItype2
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp1(sca%N1), disp2(sca%N2), base

    call MPI_Get_address(sca%N1,            disp1(1), ierr)
    call MPI_Get_address(sca%lscatfixk,     disp1(2), ierr)
    call MPI_Get_address(sca%llifetime,     disp1(3), ierr)
    call MPI_Get_address(sca%lboltzmann,    disp1(4), ierr)
    call MPI_Get_address(sca%mode,          disp1(5), ierr)
    call MPI_Get_address(sca%nsteps,        disp1(6), ierr)
    call MPI_Get_address(sca%niter,         disp1(7), ierr)
    call MPI_Get_address(sca%roottake,      disp1(8), ierr)
    call MPI_Get_address(sca%savepkk,       disp1(9), ierr)
    call MPI_Get_address(sca%naverage,      disp1(10), ierr)
    call MPI_Get_address(sca%gammamode,     disp1(11), ierr)
    call MPI_Get_address(sca%rooteps,       disp1(12), ierr)
    call MPI_Get_address(sca%kfix,          disp1(13), ierr)
    call MPI_Get_address(sca%gammaval,      disp1(14), ierr)
    call MPI_Get_address(sca%impconc,       disp1(15), ierr)
    call MPI_Get_address(sca%subarr_inp,    disp1(16), ierr)
    call MPI_Get_address(sca%maskint,       disp1(17), ierr)

    base  = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:12) =1
    blocklen1(13)   =6
    blocklen1(14:15)=1
    blocklen1(16)   =2
    blocklen1(17)   =1

    etype1(1:12) = MPI_INTEGER
    etype1(13:15)   = MPI_DOUBLE_PRECISION
    etype1(16)   = MPI_INTEGER
    etype1(17)   = MPI_INTEGER ! maskint

    call MPI_Type_create_struct(sca%N1, blocklen1, disp1, etype1, myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_sca_1'

    call MPI_Type_commit(myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_sca_1'

    call MPI_Bcast(sca%N1, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting sca_1'

    call MPI_Type_free(myMPItype1, ierr)

    if(sca%lboltzmann==1)then
      if(myrank/=master) then
        allocate( sca%weight_imp(sca%naverage), STAT=ierr )
        if(ierr/=0) stop 'Problem allocating cfg_arrays on slaves'
      end if!if(myrank/=master)

      call MPI_Get_address(sca%N2,         disp2(1), ierr)
      call MPI_Get_address(sca%weight_imp, disp2(2), ierr)

      base  = disp2(1)
      disp2 = disp2 - base

      blocklen2(1)=1
      blocklen2(2)=size(sca%weight_imp)

      etype2(1)   = MPI_INTEGER
      etype2(2)   = MPI_DOUBLE_PRECISION

      call MPI_Type_create_struct(sca%N2, blocklen2, disp2, etype2, myMPItype2, ierr)
      if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_sca_2'

      call MPI_Type_commit(myMPItype2, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_sca_2'

      call MPI_Bcast(sca%N2, 1, myMPItype2, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error brodcasting sca_2'

      call MPI_Type_free(myMPItype2, ierr)

    end if!sca%lboltzmann==1

  end subroutine myBcast_sca
#endif





#ifdef CPP_MPI
  subroutine myBcast_impcls(impcls)

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(impcls_type), intent(inout)  :: impcls

    integer :: blocklen1(impcls%N1), etype1(impcls%N1), myMPItype1
    integer :: blocklen2(impcls%N2), etype2(impcls%N2), myMPItype2
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp1(impcls%N1), disp2(impcls%N2), base

    call MPI_Get_address(impcls%N1,       disp1(1), ierr)
    call MPI_Get_address(impcls%nCluster, disp1(2), ierr)
    call MPI_Get_address(impcls%clmso,    disp1(3), ierr)

    base  = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:3)=1

    etype1(1:3) = MPI_INTEGER

    call MPI_Type_create_struct(impcls%N1, blocklen1, disp1, etype1, myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_impcls_1'

    call MPI_Type_commit(myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_impcls_1'

    call MPI_Bcast(impcls%N1, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting impcls_1'

    call MPI_Type_free(myMPItype1, ierr)


    !brodcast allocatable arrays
    if(myrank/=master) then
      allocate( impcls%RCluster(3,impcls%nCluster), impcls%ihosttype(impcls%nCluster), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating impcls_arrays on slaves'
    end if

    call MPI_Get_address(impcls%N2,        disp2(1), ierr)
    call MPI_Get_address(impcls%RCluster,  disp2(2), ierr)
    call MPI_Get_address(impcls%ihosttype, disp2(3), ierr)

    base  = disp2(1)
    disp2 = disp2 - base

    blocklen2(1)=1
    blocklen2(2)=size(impcls%RCluster)
    blocklen2(3)=size(impcls%ihosttype)

    etype2(1) = MPI_INTEGER
    etype2(2) = MPI_DOUBLE_PRECISION
    etype2(3) = MPI_INTEGER

    call MPI_Type_create_struct(impcls%N2, blocklen2, disp2, etype2, myMPItype2, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_impcls_2'

    call MPI_Type_commit(myMPItype2, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_impcls_2'

    call MPI_Bcast(impcls%N2, 1, myMPItype2, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting impcls_2'

    call MPI_Type_free(myMPItype2, ierr)

  end subroutine myBcast_impcls
#endif





  subroutine read_scattmat(inc, impcls, Amat, storeAmatin)

    use type_inc,  only: inc_TYPE
    use mod_mympi, only: myrank, master
    use mod_ioformat, only: fmt_fn_ext, filename_scattAmat, ext_formatted
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    type(inc_type),                  intent(in)  :: inc
    type(impcls_type), allocatable,  intent(out) :: impcls(:)
    double complex, allocatable,     intent(out) :: Amat(:,:,:)
    logical, optional,               intent(in)  :: storeAmatin

    integer :: clmso, ierr, lm1, lm2, ihelp, iaverage
    double complex,   allocatable :: tmat(:,:), deltamat(:,:), Gll0(:,:), Atmp1(:,:)
    character(len=256) :: filename
    logical :: file_exist, storeAmat=.false.

    integer, parameter :: iofile = 1658
    double complex, parameter :: CZERO=(0d0,0d0), CONE=(1d0,0d0)

    if(present(storeAmatin)) storeAmat = storeAmatin

    if(.not.sca_read) then
      call read_sca()
    end if

    allocate(impcls(sca%naverage), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating impcls'

    ! construct Amat on master
    if(myrank==master)then

      !check if a file contains already AMAT, otherwise construct it
      write(filename,fmt_fn_ext) filename_scattAmat, ext_formatted
      inquire(file=filename, exist=file_exist)

      !**********************************************
      if(file_exist)then !read in data from file
      !**********************************************
        write(*,*) 'Amat-file exists. Read in Amat from here.'

        do iaverage=1,sca%naverage
          call read_DTMTRX(iaverage, inc, impcls(iaverage), tmat, deltamat, .true. )
        end do!iaverage
        clmso = impcls(1)%clmso

        allocate(Amat(clmso,clmso,sca%naverage), STAT=ierr )
        if(ierr/=0) stop 'Problem allocating Amat on master'

        write(filename,fmt_fn_ext) filename_scattAmat, ext_formatted
        open(iofile,file=trim(filename),form='formatted',action='read')
        read(iofile,'(2ES25.16)') Amat
        close(iofile)

      !***************************************************
      else!file_exist  !construct AMAT from GMAT and TMAT
      !***************************************************

        do iaverage=1,sca%naverage


          !read in tmat, deltamat and Gll0
          call read_DTMTRX(iaverage, inc, impcls(iaverage), tmat, deltamat )
          if(iaverage==1)then
            clmso = impcls(1)%clmso
            allocate(Amat(clmso,clmso,sca%naverage), STAT=ierr )
            if(ierr/=0) stop 'Problem allocating Amat on master'
            Amat  = CZERO
          else!iaverage==1
            if(impcls(iaverage)%clmso/=clmso) stop "Error for averaging Pkk: clmso's are different"
          end if!iaverage==1

          call read_green_ll(iaverage, clmso, Gll0)


          !==================== BEGIN ==================!
          !=== calculate Amat = \Delta * ( 1 + G*t ) ===!

          !allocate arrays for Amat
          allocate( Atmp1(clmso,clmso), STAT=ierr )
          if(ierr/=0) stop 'Problem allocating Atmp1 on master'
          Atmp1 = CZERO
      
          ! calculate G*t
          call ZGEMM( 'N','N', clmso, clmso, clmso, CONE, Gll0, &
                    & clmso, tmat, clmso, CZERO, Atmp1, clmso   )

          ! calculate 1 + (G*t)
          do lm1=1,clmso
            Atmp1(lm1,lm1) = CONE + Atmp1(lm1,lm1)
          end do

          !calculate  \Delta * (1 + G*t)
          call ZGEMM( 'N','N', clmso, clmso, clmso, CONE, Deltamat,         &
                    & clmso, Atmp1, clmso, CZERO, Amat(:,:,iaverage), clmso )

          !=== calculate Amat = \Delta * ( 1 + G*t ) ===!
          !====================  END  ==================!

          deallocate(Atmp1, tmat, deltamat, Gll0 )
        end do!iaverage

        if(storeAmat)then
          write(filename,fmt_fn_ext) filename_scattAmat, ext_formatted
          open(iofile,file=trim(filename),form='formatted',action='write')
          write(iofile,'(2ES25.16)') Amat
          close(iofile)
        end if

      !**********************************************
      end if!file_exist
      !**********************************************

    end if!myrank==master

#ifdef CPP_MPI
    ! Brodcast impcls
    do iaverage=1,sca%naverage
      call myBcast_impcls(impcls(iaverage))
    end do
    clmso = impcls(1)%clmso

    ! Allocate Amat on the slaves
    if(myrank/=master) then
      allocate( Amat(clmso,clmso,sca%naverage), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating Amat on slaves'
    end if

    ! Brodcast Amat
    ihelp = clmso**2*sca%naverage
    call MPI_Bcast(Amat, ihelp, MPI_DOUBLE_COMPLEX, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_Success) stop 'error in Bcast(Amat)'
#endif

#ifdef CPP_DEBUG
     do iaverage=1,sca%naverage
      write(filename,'(A,I0,A,I0)') 'checkAmat_myrank=', myrank, '_iavg=', iaverage
      open(6842,file=trim(filename),form='formatted',action='write')
      do lm1=1,clmso
       do lm2=1,clmso
        write(6842,'(2I8,2ES25.16)') lm1, lm2, Amat(lm2,lm1,iaverage)
       end do!lm2
      end do!lm1
      close(6842)
     end do!iaverage
#endif

  end subroutine read_scattmat





  subroutine read_DTMTRX( iaverage, inc, impcls, tmat, deltamat, readtmatin )

    use type_inc, only: inc_type
    implicit none

    ! .... Arguments ....
    integer,           intent(in)  :: iaverage
    type(inc_type),    intent(in)  :: inc
    type(impcls_type), intent(out) :: impcls
    double complex,    allocatable, intent(out) :: tmat(:,:), deltamat(:,:)
    logical, optional,              intent(in)  :: readtmatin

    integer :: ierr, icls, lm1, lm2, idummy1, idummy2, nClustertest
    double precision :: Rclstest(3), Zatom, Rdist
    character(len=32) :: filename
    logical :: readtmat = .true.

    if(present(readtmatin)) readtmat = readtmatin

    !============ BEGIN ===========!
    !=== read the file 'DTMTRX' ===!
    if(sca%naverage==1)then
      open(unit=562, file='DTMTRX', form='formatted', action='read')
      write(*,*) 'reading file DTMTRX'
    else
      write(filename,'(A,I0)') 'DTMTRX.', iaverage
      write(*,*) 'reading file ', trim(filename)
      open(unit=562, file=trim(filename), form='formatted', action='read')
    end if

    !read the number of atoms in the impurity-cluster
    read(562, fmt=*) impcls%nCluster
    impcls%clmso = inc%lmmaxso*impcls%nCluster

    !allocate arrays
    if(readtmat)then
      allocate( tmat(impcls%clmso,impcls%clmso), deltamat(impcls%clmso,impcls%clmso), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating tmat, deltamat etc.'
    end if!readtmat

    allocate(impcls%RCluster(3,impcls%nCluster), impcls%ihosttype(impcls%nCluster), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating impcls%RCluster etc.'

    !read in the position of the atoms in the impurity cluster
    do icls=1,impcls%nCluster
      read(562, fmt=*) impcls%RCluster(:,icls)
    end do!icls

    if(readtmat)then
      !read in the t-matrix and delta-matrix
      do lm1=1,impcls%clmso
        do lm2=1,impcls%clmso
          read(562,"(2I5,4e17.9)") idummy1, idummy2, tmat(lm2,lm1), deltamat(lm2,lm1)
        end do!lm2
      end do!lm1
    end if!readtmat

    close(562)
    !=== read the file 'DTMTRX' ===!
    !============  END  ===========!

    !============ BEGIN ==========!
    !=== read the file 'scoef' ===!
    if(sca%naverage==1)then
      open(unit=565, file='scoef', form='formatted', action='read')
    else
      write(filename,'(A,I0)') 'scoef.', iaverage
      open(unit=565, file=trim(filename), form='formatted', action='read')
    end if

    read(565, fmt=*) nClustertest
    if(nClustertest/=impcls%nCluster) stop 'number of impurity atoms inconsistent'

    do icls=1,impcls%nCluster
      read(565, fmt=*) Rclstest, impcls%ihosttype(icls), Zatom, Rdist
      if( any( abs(Rclstest-impcls%RCluster(:,icls))>1d-6 ) ) stop 'Rcluster inconsistent'
    end do!icls
    close(565)
    !=== read the file 'scoef' ===!
    !============  END  ==========!

  end subroutine read_DTMTRX





  subroutine read_green_ll(iaverage,clmso, Gll0)

    implicit none
    integer,                     intent(in)  :: iaverage, clmso
    double complex, allocatable, intent(out) :: Gll0(:,:)

    double precision :: dEimag
    double complex :: energy(3), dE1, dE2, ctmp1, ctmp2
    double complex :: Gll_3(clmso, clmso ,3), &
                    &  dGll(clmso, clmso),    &
                    & d2Gll(clmso, clmso)

    integer :: ienergy, lm1, lm2, id1, id2, ierr
    character(len=32) :: filename

    double complex, parameter :: CI=(0d0,1d0)

    if(.not.allocated(Gll0)) then
      allocate( Gll0(clmso,clmso), STAT=ierr )
      if(ierr/=0) stop 'In read_green_ll: Problem allocating Gll0'
    end if!.not.allocated

    ! Read the gll_3 (for the three energies)
    if(sca%naverage==1)then
      open(unit=1283, file='GMATLL_GES', form='formatted', action='read')
    else
      write(filename,'(A,I0)') 'GMATLL_GES.', iaverage
      open(unit=1283, file=trim(filename), form='formatted', action='read')
    end if

    do ienergy=1,3

      read(1283,"(2(e17.9,X))") energy(ienergy)

      do lm1 = 1,clmso
        do lm2 = 1,clmso
          read(1283,"((2I5),(2e17.9))") id1, id2, Gll_3(lm2,lm1,ienergy)
        end do!lm2
      end do!lm1

    end do!ienergy

    ! Checks
    dE1 = energy(2)-energy(1)
    dE2 = energy(3)-energy(2)
    if(abs(dE1-dE2)>1d-8) stop '3 Energy points not equidistant'

    ! Construct first and second derivative of G(E)
    ctmp1 = 0.5d0/dE1
    ctmp2 = 1d0/dE1**2
    dGll  = ctmp1*( Gll_3(:,:,3)                    - Gll_3(:,:,1) )
    d2Gll = ctmp2*( GLL_3(:,:,3) - 2d0*Gll_3(:,:,2) + Gll_3(:,:,1) )

    ! extrapolate to energy with imag(E)=0
    dEimag = aimag(energy(2))
    Gll0 = Gll_3(:,:,2) - CI*dEimag*dGll(:,:) - 0.5d0*dEimag**2 * d2Gll(:,:)

  end subroutine read_green_ll





  subroutine extend_coeffvector2cluster_allsqa(inc,lattice,impcls,nsqa,kpoint,rveig_in,rveig_big)

    use type_inc,      only: inc_type
    use type_data,     only: lattice_type
    use mod_mathtools, only: tpiimag
    implicit none

    type(inc_type),     intent(in)  :: inc
    type(lattice_type), intent(in)  :: lattice
    type(impcls_type),  intent(in)  :: impcls
    integer,            intent(in)  :: nsqa
    double precision,   intent(in)  :: kpoint(3)
    double complex,     intent(in)  :: rveig_in(inc%lmmaxso,inc%natypd,inc%ndegen,nsqa)
    double complex,     intent(out) :: rveig_big(impcls%clmso,inc%ndegen,nsqa)

    integer :: isqa, ideg, icls, lb, ub
    double precision :: Rshift(3)
    double complex   :: expfac

    ! extent coefficient vector to size of impurity cluster
    do isqa=1,nsqa
      do ideg=1,inc%ndegen
        do icls=1,impcls%nCluster
          !calculate bloch factor
          Rshift = impcls%RCluster(:,icls) - lattice%rbasis(:,impcls%ihosttype(icls))
          expfac = exp( tpiimag*dot_product(kpoint,Rshift) )

          !copy array
          lb = (icls-1)*inc%lmmaxso+1
          ub =     icls*inc%lmmaxso
          rveig_big(lb:ub,ideg,isqa) = expfac*rveig_in(:,impcls%ihosttype(icls),ideg,isqa)
        end do!icls
      end do!ispin1
    end do!isqa

  end subroutine extend_coeffvector2cluster_allsqa





  subroutine extend_coeffvector2cluster_allkpt(inc,lattice,impcls,nkpts,kpoints,rveig_in,rveig_big)

    use type_inc,      only: inc_type
    use type_data,     only: lattice_type
    use mod_mathtools, only: tpiimag
    implicit none

    type(inc_type),     intent(in)  :: inc
    type(lattice_type), intent(in)  :: lattice
    type(impcls_type),  intent(in)  :: impcls
    integer,            intent(in)  :: nkpts
    double precision,   intent(in)  :: kpoints(3,nkpts)
    double complex,     intent(in)  :: rveig_in(inc%lmmaxso,inc%natypd,nkpts)
    double complex,     intent(out) :: rveig_big(impcls%clmso,nkpts)

    integer :: icls, lb, ub, k
    double precision :: Rshift(3)
    double complex   :: expfac

    ! extent coefficient vector to size of impurity cluster
    do icls=1,impcls%nCluster
      !calculate bloch factor
      Rshift = impcls%RCluster(:,icls) - lattice%rbasis(:,impcls%ihosttype(icls))
      do k=1,nkpts
        expfac = exp( tpiimag*dot_product(reshape(kpoints(1:3,k),[3]),Rshift) )

        !copy array
        lb = (icls-1)*inc%lmmaxso+1
        ub =     icls*inc%lmmaxso
        rveig_big(lb:ub,k) = expfac*rveig_in(:,impcls%ihosttype(icls),k)
      end do!k
    end do!icls

  end subroutine extend_coeffvector2cluster_allkpt




  subroutine extend_coeffvector2cluster(inc,lattice,impcls,kpoint,rveig_in,rveig_big)

    use type_inc,      only: inc_type
    use type_data,     only: lattice_type
    use mod_mathtools, only: tpiimag
    implicit none

    type(inc_type),     intent(in)  :: inc
    type(lattice_type), intent(in)  :: lattice
    type(impcls_type),  intent(in)  :: impcls
    double precision,   intent(in)  :: kpoint(3)
    double complex,     intent(in)  :: rveig_in(inc%lmmaxso,inc%natypd)
    double complex,     intent(out) :: rveig_big(impcls%clmso)

    integer :: icls, lb, ub
    double precision :: Rshift(3)
    double complex   :: expfac

    ! extent coefficient vector to size of impurity cluster
    do icls=1,impcls%nCluster
      !calculate bloch factor
      Rshift = impcls%RCluster(:,icls) - lattice%rbasis(:,impcls%ihosttype(icls))
      expfac = exp( tpiimag*dot_product(kpoint,Rshift) )

      !copy array
      lb = (icls-1)*inc%lmmaxso+1
      ub =     icls*inc%lmmaxso
      rveig_big(lb:ub) = expfac*rveig_in(:,impcls%ihosttype(icls))
    end do!icls

  end subroutine extend_coeffvector2cluster



  subroutine calculate_lifetimeaverage_vis(inc,nsqa,nkpts,nkpts_all,kpt2irr,kpoints,fermivel,taukinv,printout,fac,fsRyunit)
    use type_inc,  only: inc_type
    use mod_mympi, only: myrank, master
    use mod_mathtools, only: crossprod, simple_integration
    use mpi
    implicit none

    type(inc_type),   intent(in) :: inc
    integer,          intent(in) :: nsqa, nkpts, nkpts_all, kpt2irr(nkpts_all)
    double precision, intent(in) :: kpoints(3,nkpts), fermivel(3,nkpts), taukinv(inc%ndegen**2*nsqa,nkpts)
    logical,          intent(in) :: printout
    double precision, intent(in) :: fac
    character(len=*), intent(in) :: fsRyunit

    integer :: iloop, isqa, itri, pointspick(3)
    double precision :: integ(inc%ndegen**2*nsqa), tauinvss(inc%ndegen**2*nsqa), kpoints_triangle(3,3), fermivel_triangle(3,3), taukinv_triangle(3), k21(3), k31(3), kcross(3), area, d_integ, d_dos, dos

    write(*,*) 'in calculate_lifetimeaverage_vis'

    dos   = 0d0
    integ = 0d0

    do itri=1,nkpts_all/3

      pointspick(1) = kpt2irr((itri-1)*3+1)
      pointspick(2) = kpt2irr((itri-1)*3+2)
      pointspick(3) = kpt2irr((itri-1)*3+3)

      !rewrite corner point data
      kpoints_triangle(:,:)  = kpoints(:,pointspick(:))
      fermivel_triangle(:,:) = fermivel(:,pointspick(:))

      !compute area
      k21 = kpoints_triangle(:,2) - kpoints_triangle(:,1)
      k31 = kpoints_triangle(:,3) - kpoints_triangle(:,1)
      call crossprod(k21, k31, kcross)
      area = 0.5d0*sqrt(sum(kcross**2))

      do iloop=1,inc%ndegen**2*nsqa
        taukinv_triangle(:) = taukinv(iloop,pointspick(:))
        call simple_integration(area, fermivel_triangle, taukinv_triangle, d_integ, d_dos)
        integ(iloop) = integ(iloop)+d_integ
      end do!iloop

      dos = dos+d_dos

    end do!ikp

    tauinvss = (integ*fac)/dos

    if(printout .and. myrank==master)then
      if(inc%ndegen==2) then
        do isqa=1,nsqa
          write(*,'(A,I3,A,ES25.16,1X,A)') "spin conserving lifetime: isqa= ", isqa, ", tau= ", 2d0/(tauinvss(1+4*(isqa-1)) + tauinvss(4+4*(isqa-1))), trim(fsRyunit)
          write(*,'(A,I3,A,ES25.16,1X,A)') "spin flip       lifetime: isqa= ", isqa, ", T_1= ", 1d0/(tauinvss(2+4*(isqa-1)) + tauinvss(3+4*(isqa-1))), trim(fsRyunit)
        end do!isqa
      elseif(inc%ndegen==1)then
        do isqa=1,nsqa
          write(*,'(A,I3,A,ES25.16,1X,A)') "lifetime: isqa= ", isqa, ", tau= ", 1d0/tauinvss(isqa), trim(fsRyunit)
        end do
      else
        write(*,*) 'Calculation of lifetimes not implemented for inc%ndegen /= 2'
      end if!inc%ndegen==2
    end if!printout

  end subroutine calculate_lifetimeaverage_vis


  subroutine calculate_lifetimeaverage_int(inc,nsqa,nkpts,weights,taukinv,printout,fac,fsRyunit)
    use type_inc,  only: inc_type
    use mod_mympi, only: myrank, master
    use mod_mathtools, only: crossprod, simple_integration
    use mpi
    implicit none

    type(inc_type),   intent(in) :: inc
    integer,          intent(in) :: nsqa, nkpts
    double precision, intent(in) :: weights(nkpts), taukinv(inc%ndegen**2*nsqa,nkpts)
    logical,          intent(in) :: printout
    double precision, intent(in) :: fac
    character(len=*), intent(in) :: fsRyunit

    integer :: isqa, ikp
    double precision :: integ(inc%ndegen**2*nsqa), tauinvss(inc%ndegen**2*nsqa)

    write(*,*) 'in calculate_lifetimeaverage_int'

    integ = 0d0
    do ikp=1,nkpts
      integ = integ + weights(ikp)*taukinv(:,ikp)
    end do!ikp

    tauinvss = (integ*fac)/sum(weights)

    if(printout .and. myrank==master)then
      if(inc%ndegen==2) then
        do isqa=1,nsqa
          write(*,'(A,I3,A,ES25.16,1X,A)') "spin conserving lifetime: isqa= ", isqa, ", tau= ", 2d0/(tauinvss(1+4*(isqa-1)) + tauinvss(4+4*(isqa-1))), trim(fsRyunit)
          write(*,'(A,I3,A,ES25.16,1X,A)') "spin flip       lifetime: isqa= ", isqa, ", T_1= ", 1d0/(tauinvss(2+4*(isqa-1)) + tauinvss(3+4*(isqa-1))), trim(fsRyunit)
        end do!isqa
      elseif(inc%ndegen==1)then
        do isqa=1,nsqa
          write(*,'(A,I3,A,ES25.16,1X,A)') "lifetime: isqa= ", isqa, ", tau= ", 1d0/tauinvss(isqa), trim(fsRyunit)
        end do
      else
        write(*,*) 'Calculation of lifetimes not implemented for inc%ndegen /= 2'
      end if!inc%ndegen==2
    end if!printout

  end subroutine calculate_lifetimeaverage_int





  subroutine create_subarr_comm( subarr_dim, myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col )
    use mpi
    implicit none
    integer, intent(in) :: subarr_dim(2)
    integer, intent(out) :: myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col, myrank_grid, myrank_row, myrank_col

    integer :: ierr
    logical :: logic2(2)
    logical, parameter :: periodic(2) = .false., reorder(2) = .false.


    call MPI_Cart_create( MPI_COMM_WORLD, 2, subarr_dim, periodic, reorder, myMPI_comm_grid, ierr )
    call MPI_Comm_rank( myMPI_comm_grid, myrank_grid, ierr )

    logic2 = (/ .true., .false. /)
    call MPI_Cart_sub( myMPI_comm_grid, logic2, myMPI_comm_row, ierr ) ! row communicator
    logic2 = (/ .false., .true. /)
    call MPI_Cart_sub( myMPI_comm_grid, logic2, myMPI_comm_col, ierr ) ! col communicator

    call MPI_Comm_rank( myMPI_comm_row, myrank_row, ierr )
    call MPI_Comm_rank( myMPI_comm_col, myrank_col, ierr )

  end subroutine create_subarr_comm




  subroutine create_subarr(subarr_dim, ntot, dataarr_lb, dataarr_ub, dataarr_nkpt)
    use mod_mympi, only: myrank, nranks, master
    implicit none
    ! This subroutine reads in the required data for computation of the NxN-matrix Pkk' with minimum amount of memory
    ! It is achieved by dividing the NxN-matrix into AxB blocks, where each MPI-rank treats one subblock and needs only
    ! the corresponding data (eigenvectors).
    !

    integer, intent(in)  :: ntot, subarr_dim(2)
    integer, intent(out) :: dataarr_lb(0:nranks-1,2), &
                          & dataarr_ub(0:nranks-1,2), &
                          & dataarr_nkpt(0:nranks-1,2)

    integer :: ierr, idimen, irank, irank1, irank2, nranks2, irest
    integer, allocatable :: subblocks_nkpt(:,:), subblocks_ioff(:,:)
    character(len=80) :: fmtstr

    if(myrank==master .and. product(subarr_dim)<nranks) write(*,'("WARNING: ",I0," ranks stay idle")') nranks-product(subarr_dim)
    if(product(subarr_dim)>nranks) stop 'Error: number of blocks is larger than number of ranks'

    allocate( subblocks_nkpt(0:maxval(subarr_dim)-1,2), subblocks_ioff(0:maxval(subarr_dim)-1,2), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating kpoint-dimensions in subblock'

    subblocks_nkpt= -1
    subblocks_ioff= -1
    do idimen=1,2

      nranks2 = subarr_dim(idimen)

      subblocks_nkpt(0:nranks2-1,idimen) = int(ntot/nranks2)
      subblocks_ioff(0:nranks2-1,idimen) = int(ntot/nranks2)*(/ (irank, irank=0,nranks2-1) /)
      irest   = ntot-int(ntot/nranks2)*nranks2

      if(irest>0) then

        do irank=0,irest-1
          subblocks_nkpt(irank,idimen) = subblocks_nkpt(irank,idimen) + 1
          subblocks_ioff(irank,idimen) = subblocks_ioff(irank,idimen) + irank
        end do!irank

        do irank=irest,nranks2-1
          subblocks_ioff(irank,idimen) = subblocks_ioff(irank,idimen) + irest
        end do!irank

      end if!irest>0

    end do!idimen
    

    if(myrank==master)then
      write(*,'("=== DISTRIBUTION OF K-POINTS ON A GRID: ===")')
      irank=0
      do irank1=0,subarr_dim(1)-1
       do irank2=0,subarr_dim(2)-1
        write(*,'(2X,"Processor ",I0," treats ",I0," x ",I0," submatrix.")') &
                 & irank1*subarr_dim(2) + irank2, subblocks_nkpt(irank1,1), subblocks_nkpt(irank2,2)
       end do!irank2
      end do!irank1
!     write(fmtstr,'(A,I0,A)') '("--> ",',maxval(subarr_dim),'I8," <--")'
!     write(*,'("=== subblocks_nkpt ===")')
!     write(*,fmtstr) subblocks_nkpt
!     write(*,'("=== subblocks_ioff ===")')
!     write(*,fmtstr) subblocks_ioff
    end if

    !determine the lower and upper bounds of receive-buffers
    do irank1=0,subarr_dim(1)-1
     do irank2=0,subarr_dim(2)-1
      irank = irank1*subarr_dim(2) + irank2
      dataarr_lb(irank,1) = subblocks_ioff(irank1,1)
      dataarr_lb(irank,2) = subblocks_ioff(irank2,2)
      dataarr_ub(irank,1) = subblocks_ioff(irank1,1)+subblocks_nkpt(irank1,1)
      dataarr_ub(irank,2) = subblocks_ioff(irank2,2)+subblocks_nkpt(irank2,2)
      dataarr_nkpt(irank,1) = subblocks_nkpt(irank1,1)
      dataarr_nkpt(irank,2) = subblocks_nkpt(irank2,2)
      end do!irank2
    end do!irank1

    deallocate(subblocks_nkpt, subblocks_ioff)

  end subroutine create_subarr





  subroutine read_eigv_part(inc, nsqa, ioff, nkpt, file_comm, subrank, subcomm, rveig)

    use type_inc,     only: inc_type
    use mod_ioformat, only: filename_eigvect, filemode_int, ext_mpiio
    use mpi

    type(inc_TYPE),              intent(in)  :: inc
    integer,                     intent(in)  :: nsqa, ioff, nkpt, file_comm, subrank, subcomm
    double complex, allocatable, intent(out) :: rveig(:,:,:,:,:)

    ! MPI variables
    integer :: myMPI_iotype, filehandle
    integer :: mpistatus(MPI_STATUS_SIZE)
    integer(kind=MPI_OFFSET_KIND) :: disp=0
    integer, parameter :: mo = kind(disp)

    integer :: ierr, ihelp
    integer, parameter :: submaster=0
    character(len=256) :: filename

    ihelp=inc%lmmaxso*inc%natypd*inc%ndegen*nsqa

    !allocate the result-arrays
    allocate( rveig(inc%lmmaxso, inc%natypd, inc%ndegen, nsqa, nkpt), STAT=ierr )
    if(ierr/=0) stop 'problem allocating rveig_dim'

    if(subrank==submaster)then

      !open the eigenvector file
      write(filename,'(A,A,A)') filename_eigvect, filemode_int, ext_mpiio
      call MPI_File_open( file_comm, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, filehandle, ierr )
      if(ierr/=MPI_SUCCESS) stop 'MPI_File_open filename_eigv_rot'

      !prepare the file view (see only the part of this task)
      call MPI_Type_Contiguous( ihelp*nkpt, MPI_DOUBLE_COMPLEX, myMPI_iotype,ierr )
      if(ierr/=MPI_SUCCESS) stop 'MPI_Type_Vector myMPI_iotype'

      !commit the datatype
      call MPI_Type_commit(myMPI_iotype,ierr)
      if(ierr/=MPI_SUCCESS) stop 'MPI_Type_commit myMPI_iotype'

      !set the file view
      disp = ioff*(ihelp*16_mo)!byte
      call MPI_File_set_view( filehandle, disp, MPI_DOUBLE_COMPLEX, myMPI_iotype, 'native', MPI_INFO_NULL, ierr )
      if(ierr/=MPI_SUCCESS) stop 'MPI_File_set_view filehandle'

      !read
      call MPI_File_Read( filehandle, rveig, ihelp*nkpt, MPI_DOUBLE_COMPLEX, mpistatus, ierr )

      !free+close
      call MPI_Type_free(myMPI_iotype,ierr)
      call MPI_File_close( filehandle, ierr )
      if(ierr/=MPI_SUCCESS) stop 'MPI_File_close'

    end if!subrank==submaster

    !Brodcast to other processes in subcomm
    call MPI_Bcast(rveig, ihelp*nkpt, MPI_DOUBLE_COMPLEX, submaster, subcomm, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Error in Bcast(rveig)'

  end subroutine read_eigv_part




  subroutine Pkkmpifile_setview(rwmode, my_mpi_comm, nkpts, nkpt1, nkpt2, ioff1, ioff2, ndegen, nsqa, filehandle)
    use mod_ioformat, only: filename_scattmat, ext_mpiio, fmt_fn_ext
    use mod_mympi, only: myrank
    use mpi
    implicit none

    character(len=*), intent(in) :: rwmode
    integer, intent(in) :: my_mpi_comm, nkpts, nkpt1, nkpt2, ioff1, ioff2, ndegen, nsqa
    integer, intent(out) :: filehandle

    integer :: myMPI_iotype, my_mpi_datatype
    integer(kind=MPI_OFFSET_KIND) :: disp=0, datasize=0, nkpts_big=0
    integer, parameter :: mo = kind(disp)

    integer :: ierr, ihelp0, ihelp1, ihelp2
    character(len=256) :: filename

    my_mpi_datatype = MPI_DOUBLE_PRECISION
    datasize=8_mo

    write(filename,fmt_fn_ext) filename_scattmat, ext_mpiio
    select case (rwmode)
    case ('read')
      call MPI_File_open( my_mpi_comm, trim(filename), MPI_MODE_RDONLY, &
                        & MPI_INFO_NULL, filehandle, ierr            )
      if(ierr/=MPI_SUCCESS) stop 'error: Pkkmpifile_setview readmode'
    case ('write')
      call MPI_File_open( my_mpi_comm, trim(filename), MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                        & MPI_INFO_NULL, filehandle, ierr                            )
      if(ierr/=MPI_SUCCESS) stop 'error: Pkkmpifile_setview writemode'
    case default; stop 'dont know how to open the mpi-file. only mode = read/write'
    end select

    !prepare the file view (see only the part of this task)
    ihelp0 = nkpts*ndegen
    ihelp1 = nkpt1*ndegen
    ihelp2 = nkpt2*ndegen*nsqa
    call MPI_Type_Vector( ihelp2, ihelp1, ihelp0, my_mpi_datatype, myMPI_iotype, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error: MPI_Type in Pkkmpifile_setview'

    call MPI_Type_commit(myMPI_iotype,ierr)
    if(ierr/=MPI_SUCCESS) stop 'error: MPI_Type_commit in Pkkmpifile_setview'

    nkpts_big = nkpts*1_mo
    disp = (ioff1+ioff2*ndegen*nkpts_big*nsqa)*(ndegen*datasize)!byte

    call MPI_File_set_view( filehandle, disp, my_mpi_datatype,             &
                          & myMPI_iotype, 'native', MPI_INFO_NULL, ierr )
    if(ierr/=MPI_SUCCESS) stop 'error: MPI_File_set_view in Pkkmpifile_setview'

  end subroutine Pkkmpifile_setview




  subroutine Pkkmpifile_write(my_mpi_comm, nkpts, nkpt1, nkpt2, ioff1, ioff2, ndegen, nsqa, Pkksub)
    use mpi
    implicit none

    integer,          intent(in) :: my_mpi_comm, nkpts, nkpt1, nkpt2, ioff1, ioff2, ndegen, nsqa
    double precision, intent(in) :: Pkksub(ndegen,nkpt1,ndegen,nsqa,nkpt2)

    integer :: ierr, filehandle, ihelp, mpistatus(MPI_STATUS_SIZE)

    call Pkkmpifile_setview('write', my_mpi_comm, nkpts, nkpt1, nkpt2, ioff1, ioff2, ndegen, nsqa, filehandle)
    ihelp = (ndegen**2)*nkpt1*nkpt2*nsqa
    call MPI_File_write_all(filehandle, Pkksub, ihelp, MPI_DOUBLE_PRECISION, mpistatus, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error writing in Pkkmpifile_write'

    call MPI_File_close(filehandle, ierr)

  end subroutine Pkkmpifile_write


  subroutine Pkkmpifile_read(my_mpi_comm, nkpts, nkpt1, nkpt2, ioff1, ioff2, ndegen, nsqa, Pkksub)
    use mpi
    implicit none

    integer,                        intent(in) :: my_mpi_comm, nkpts, nkpt1, nkpt2, ioff1, ioff2, ndegen, nsqa
    double precision, allocatable, intent(out) :: Pkksub(:,:,:,:,:)

    integer :: ierr, filehandle, ihelp, mpistatus(MPI_STATUS_SIZE)

    allocate( Pkksub(ndegen,nkpt1,ndegen,nsqa,nkpt2),&
            & STAT=ierr )
    if(ierr/=0) stop 'Problem allocating Pkksub in Pkkmpifile_read'
    Pkksub = 0d0

    call Pkkmpifile_setview('read', my_mpi_comm, nkpts, nkpt1, nkpt2, ioff1, ioff2, ndegen, nsqa, filehandle)
    ihelp = (ndegen**2)*nkpt1*nkpt2*nsqa
    call MPI_File_read_all(filehandle, Pkksub, ihelp, MPI_DOUBLE_PRECISION, mpistatus, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error reading in Pkkmpifile_read'

    call MPI_File_close(filehandle, ierr)

  end subroutine Pkkmpifile_read


  subroutine project_fermivel_newaxis(nkpts,fermivel)
    implicit none
    integer,          intent(in)    :: nkpts
    double precision, intent(inout) :: fermivel(3,nkpts)

    double precision :: fermivel_tmp(3,nkpts), n1(3), n2(3), n3(3), sqrt13
    integer :: ikp


    n1=(/ 1d0, -1d0,  0d0 /)/sqrt(2d0)
    n2=(/ 1d0,  1d0, -2d0 /)/sqrt(6d0)
    n3=(/ 1d0,  1d0,  1d0 /)/sqrt(3d0)

    fermivel_tmp = fermivel

    do ikp=1,nkpts
      fermivel(1,ikp) = sum(fermivel_tmp(:,ikp)*n1)
      fermivel(2,ikp) = sum(fermivel_tmp(:,ikp)*n2)
      fermivel(3,ikp) = sum(fermivel_tmp(:,ikp)*n3)
    end do!ikp

  end subroutine project_fermivel_newaxis




  subroutine save_lifetime(filemode, nkpts, nsqa, ndegen, taukinv, nsym, isym, fac, fsRyunit)

    use mod_ioformat,   only: fmt_fn_sub_ext, fmt_fn_isqa_sub_ext, ext_formatted, filename_lifetime
    implicit none

    character(len=*), intent(in) :: filemode
    integer,          intent(in) :: nkpts, nsqa, ndegen, nsym, isym(nsym)
    double precision, intent(in) :: taukinv(ndegen**2*nsqa,nkpts), fac
    character(len=*), intent(in) :: fsRyunit

    integer :: ii, ikp, isqa, imin, imax
    double precision :: taukinv_tmp(4), T1_inv, taup_inv, tauk_sum(2)
    character(len=256) :: filename, fmtstr

    if(ndegen==2)then
      write(fmtstr,'(A,I0,A)') '(', 8,'ES25.16)'

      do isqa=1,nsqa
        if(nsqa>1)then
          write(filename,fmt_fn_isqa_sub_ext) filename_lifetime, isqa, filemode, ext_formatted
        else!nsqa>1
          write(filename,fmt_fn_sub_ext) filename_lifetime, filemode, ext_formatted
        end if!nsqa>1

        open(unit=326529, file=trim(filename), form='formatted', action='write')
        write(326529,'(A,A)') '# output in ', trim(fsRyunit)
        write(326529,'(A)')   '# tau_p, T_1, tau^u, tau^d, {tau^uu, ^du, ^ud, ^dd}'
        write(326529,'(A)')   '#                           The second spin index refers to the incoming k'
        write(326529,'(A)')   '#                           (e.g. 1/tau^u = 1/tau^uu + 1/tau^du)'

        write(326529,'(2I8)') nkpts, nsym, nsqa, ndegen
        write(326529,'(12I8)') isym

        imin = (isqa-1)*4+1
        imax = isqa*4
        do ikp=1,nkpts
          taukinv_tmp =  taukinv(imin:imax,ikp)*fac

          taup_inv    = (taukinv_tmp(1) + taukinv_tmp(4))/2  !spin-conserving relaxation time
          T1_inv      =  taukinv_tmp(2) + taukinv_tmp(3)     !Spin relaxation time
          tauk_sum(1) =  taukinv_tmp(1) + taukinv_tmp(3)     !Momentum relaxation time for spin-up
          tauk_sum(2) =  taukinv_tmp(2) + taukinv_tmp(4)     !Momentum relaxation time for spin-down
          write(326529,fmtstr) 1d0/taup_inv, 1d0/T1_inv, 1d0/tauk_sum, 1d0/taukinv_tmp
        end do!ikp

        close(326529)
      end do!isqa

    else!ndegen==2

      write(fmtstr,'(A,I0,A)') '(', nsqa,'ES25.16)'

      write(filename,fmt_fn_sub_ext) filename_lifetime, filemode, ext_formatted
      open(unit=326529, file=trim(filename), form='formatted', action='write')
      write(326529,'(A,A)') '# output in ', trim(fsRyunit)
      write(326529,'(A)')   '# tau_k, one column for each SQA'
      write(326529,'(2I8)') nkpts, nsym, nsqa, ndegen
      write(326529,'(12I8)') isym
      write(326529,fmtstr) fac/taukinv
      close(326529)

    endif!ndegen==2

  end subroutine save_lifetime



  subroutine save_scattfix(filemode, nkpts, nsqa, ndegen, kfix, Pkkfix, nsym, isym, mode)

    use mod_ioformat,   only: fmt_fn_sub_ext, ext_formatted, filename_scattfixk, filename_scattfixk_out
    implicit none

    character(len=*), intent(in) :: filemode
    character(len=*), intent(in) :: mode
    integer,          intent(in) :: nkpts, nsqa, ndegen, nsym, isym(nsym)
    double precision, intent(in) :: kfix(3), Pkkfix(ndegen,ndegen,nsqa,nkpts)

    integer :: ii
    !character(len=256) :: filename, fmtstr
    character(len=512) :: filename, fmtstr

    write(fmtstr,'(A,I0,A)') '(', nsqa*ndegen**2,'ES25.16)'
    write(*,*) 'in save_sattfix:',mode,mode=='start'
    write(*,*) 'filename:',len(filename_scattfixk),filename_scattfixk
    write(*,*) 'filename out:',len(filename_scattfixk_out),filename_scattfixk_out
    write(*,*) 'filemode:',len(filemode),filemode

    if (mode=='start') then
       write(filename,fmt_fn_sub_ext) filename_scattfixk, filemode, ext_formatted
    else
       write(filename,fmt_fn_sub_ext) filename_scattfixk_out,filemode, ext_formatted
    endif
    open(unit=326529, file=trim(filename), form='formatted', action='write')

    if (mode=='start') then
       write(326529,'(A,3ES25.16)') '# incoming k-vector = ', kfix
    else
       write(326529,'(A,3ES25.16)') '# outgoing k-vector = ', kfix
    endif
    if(ndegen==2)then
      write(326529,'(A)') '# for each SQA: P^uu, P^du, P^ud, P^dd, where the second spin index refers to the incoming k'
    else!ndegen==2
      write(326529,'(A)') '# one column for each SQA'
    end if!ndegen==2

    write(326529,'(2I8)') nkpts, nsym, nsqa, ndegen
    write(326529,'(12I8)') isym
    write(326529,fmtstr) Pkkfix
    close(326529)

  end subroutine save_scattfix

end module mod_scattering
