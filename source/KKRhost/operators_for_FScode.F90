!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Interface routine to normcoeff routines that prepare operators for use in `FScode` (compuation of spin expectation value etc.)
!> Author:
!> Interface routine to normcoeff routines that prepare operators for
!> use in `FScode` (compuation of spin expectation value etc.)
!> First wavefuncitons are read in and converted to old mesh and then
!> `normcoeff_*` routines are called
!------------------------------------------------------------------------------------
module mod_operators_for_fscode

contains

  !-------------------------------------------------------------------------------
  !> Summary: Interface routine to normcoeff routines that prepare operators for use in `FScode` (compuation of spin expectation value etc.)
  !> Author:
  !> Category: physical-observables, KKRhost
  !> Deprecated: False
  !> Interface routine to normcoeff routines that prepare operators for
  !> use in `FScode` (compuation of spin expectation value etc.)
  !> First wavefuncitons are read in and converted to old mesh and then
  !> `normcoeff_*` routines are called
  !-------------------------------------------------------------------------------
  subroutine operators_for_fscode(korbit, operator_imp)

#ifdef CPP_MPI
    use :: mod_types, only: t_inc, t_mpi_c_grid
    use :: mpi
    use :: mod_mympi, only: nranks, master, myrank, distribute_linear_on_tasks
#else
    use :: mod_types, only: t_inc
    use :: mod_mympi, only: nranks, master, myrank
#endif
    use :: mod_wunfiles, only: t_params
    use :: mod_runoptions, only: impurity_operator_only
    use :: mod_save_wavefun, only: t_wavefunctions, read_wavefunc
    use :: mod_types, only: t_imp
    use :: mod_datatypes, only: dp
    use :: mod_cheb2oldgrid
    use :: mod_normcoeff_so
    use :: mod_normcoeff_so_spinflux
    use :: mod_normcoeff_so_torq
    use :: mod_rotatespinframe, only: rotatematrix
    use :: mod_constants, only: czero
    use :: global_variables, only: lmmaxd, nspind

    implicit none

    integer, intent (in) :: korbit
    logical, intent (in) :: operator_imp !! logical that determines if second part of computing operators with impurity wavefunctions is done or not

    ! read in wavefunctions
    logical :: rll_was_read_in, sll_was_read_in, rllleft_was_read_in, sllleft_was_read_in
    complex (kind=dp), dimension(:, :), allocatable :: rlltemp
    complex (kind=dp), dimension(:, :), allocatable :: pnstemp
    complex (kind=dp), dimension(:, :, :, :), allocatable :: rll
    complex (kind=dp), dimension(:, :, :, :), allocatable :: sll
    complex (kind=dp), dimension(:, :, :, :), allocatable :: pns_so
    complex (kind=dp), dimension(:, :, :, :), allocatable :: rllleft
    complex (kind=dp), dimension(:, :, :, :), allocatable :: sllleft
    complex (kind=dp), dimension(:, :, :, :, :), allocatable :: pns_so_all

    ! loop counter etc.
    integer :: ie, ie_start, ie_end, ie_num, lm1, lm2, ir, i1, i1_start, i1_end, ierr

    ! array dimensions
    integer :: irmd, natyp, nsra, ncheb, ntotd

    ! arrays for rmeshes (old and new), and nonco_angles
    integer, dimension(:), allocatable :: irws, npan_tot
    integer, dimension(:,:), allocatable :: ipan_intervall
    real (kind=dp), dimension(:), allocatable :: theta, phi
    real (kind=dp), dimension(:,:), allocatable :: rmesh, rpan_intervall

    ! for impurity-wavefunction related stuff
    complex (kind=dp), dimension(:, :, :, :, :), allocatable :: pns_so_imp
    integer :: i1_imp, natomimp

#ifdef CPP_MPI
    ! communcate PNS_SO_ALL for OPERATOR option
    integer :: ihelp
    integer, dimension(0:nranks-1) :: ntot_pt, ioff_pt
    complex (kind=dp), dimension(:, :, :, :, :), allocatable :: work
#endif

    ! for TEST options

    if (t_inc%i_write>0) write (1337, *) 'start computing Operators'
    if (t_inc%i_write>0) write (*, *) 'start computing Operators'

    !--------------------------------------------------------------------------------
    ! Part 1: operators for host wavefunctions
    !--------------------------------------------------------------------------------

    ! first fill scalar and array parameters that are used here
    ! call get_params_operators(lmmaxd, irmd, natyp, nsra, ncheb, ntot, irws,
    ! scalars
    irmd    = t_params%irm
    natyp   = t_params%natyp
    nsra    = t_params%nsra
    ncheb   = t_params%ncheb
    ntotd   = t_params%ntotd
    ! arrays
    allocate (irws(natyp))
    irws = t_params%irws
    allocate (rmesh(irmd,natyp))
    rmesh = t_params%rmesh
    allocate (npan_tot(natyp))
    npan_tot = t_params%npan_tot
    allocate (rpan_intervall(0:ntotd,natyp))
    rpan_intervall = t_params%rpan_intervall
    allocate (ipan_intervall(0:ntotd,natyp))
    ipan_intervall = t_params%ipan_intervall
    allocate (theta(natyp), phi(natyp))
    theta = t_params%theta
    phi = t_params%phi

    if (.not. impurity_operator_only) then ! test option to disable costly recalculation of host operators

      if (t_inc%i_write>0) write (1337, *) 'Operators using host wavefunctions'
      if (t_inc%i_write>0) write (*, *) 'Operators using host wavefunctions'

      ! now get the radial wavefunctions in the correct (i.e. old) radial mesh
      ! called PNS_SO_ALL

      allocate (pns_so_all(lmmaxd,lmmaxd,irmd,2,natyp))
      pns_so_all = czero

#ifdef CPP_MPI
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie,                       &
        t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at, master, natyp, ntot_pt,      &
        ioff_pt, .true.)

      i1_start  = ioff_pt(t_mpi_c_grid%myrank_ie) + 1
      i1_end    = ioff_pt(t_mpi_c_grid%myrank_ie) + ntot_pt(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1    = ntot_pt(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot_pt1 = ntot_pt
      t_mpi_c_grid%ioff_pt1 = ioff_pt
#else
      i1_start = 1
      i1_end = natyp
#endif

      do i1 = i1_start, i1_end

        allocate (rll(nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:0))
        allocate (sll(nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:0))
        allocate (rllleft(nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:0))
        allocate (sllleft(nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:0))
        allocate (pns_so(lmmaxd,lmmaxd,irmd,2))
        rll     = czero
        sll     = czero
        rllleft = czero
        sllleft = czero
        pns_so  = czero

#ifdef CPP_MPI
        ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
        ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
        ie_start = 0               ! offset
        ie_end = t_params%ielast
#endif

        do ie_num = 1, ie_end
          ie = ie_start + ie_num
          ! make sure only calculated at the Fermi level
          ! if(ie_end==1 .or. ie==1) then
          if (ie==1) then
            if (t_wavefunctions%nwfsavemax>0) then ! read wavefunctions?
              ! read in wavefunction from memory
              call read_wavefunc(t_wavefunctions,rll,rllleft,sll,sllleft,i1,ie,nsra,&
                lmmaxd, t_inc%irmdnew,0,1,rll_was_read_in,sll_was_read_in,          &
                rllleft_was_read_in,sllleft_was_read_in)
            end if                 ! t_wavefunctions%Nwfsavemax
          end if                   ! ie==1
        end do                     ! ie_num=1,ie_end

        ! transform radial wavefunction back to old mesh
        allocate (rlltemp(t_inc%irmdnew,lmmaxd))
        allocate (pnstemp(irws(i1),lmmaxd))
        do lm1 = 1, lmmaxd
          rlltemp = czero
          pnstemp = czero
          ir = 0
          do ir = 1, t_inc%irmdnew
            do lm2 = 1, lmmaxd
              rlltemp(ir, lm2) = rll(lm1, lm2, ir, 0)
            end do
          end do
          call cheb2oldgrid(irws(i1),t_inc%irmdnew,lmmaxd,rmesh(:,i1),ncheb,        &
            npan_tot(i1),rpan_intervall(:,i1),ipan_intervall(:,i1),rlltemp,pnstemp, &
            irmd)
          do ir = 1, irws(i1)
            do lm2 = 1, lmmaxd
              pns_so(lm1, lm2, ir, 1) = pnstemp(ir, lm2)
            end do
          end do
        end do                     ! LM1
        ! for small component
        if (nsra==2) then
          do lm1 = 1, lmmaxd
            rlltemp = czero
            pnstemp = czero
            do ir = 1, t_inc%irmdnew
              do lm2 = 1, lmmaxd
                rlltemp(ir, lm2) = rll(lm1+lmmaxd, lm2, ir, 0)
              end do
            end do
            call cheb2oldgrid(irws(i1),t_inc%irmdnew,lmmaxd,rmesh(:,i1),ncheb,      &
              npan_tot(i1),rpan_intervall(:,i1),ipan_intervall(:,i1),rlltemp,       &
              pnstemp,irmd)
            do ir = 1, irws(i1)
              do lm2 = 1, lmmaxd
                pns_so(lm1, lm2, ir, 2) = pnstemp(ir, lm2)
              end do
            end do
          end do                   ! LM1
        end if                     ! NSRA.EQ.2

        ! rotate radial wavefunction to global frame
        do ir = 1, irmd
          call rotatematrix(pns_so(1,1,ir,1), theta(i1), phi(i1), lmmaxd/2, 0)
          call rotatematrix(pns_so(1,1,ir,2), theta(i1), phi(i1), lmmaxd/2, 0)
        end do

        ! finally collect wavefuncitons in global frame and old mesh for all atoms to be used in normcoeff-routines below
        pns_so_all(:, :, :, :, i1) = pns_so(:, :, :, :)

        deallocate (rll, sll, rllleft, sllleft, pns_so)
        deallocate (rlltemp)
        deallocate (pnstemp)

      end do                       ! I1=i1_start, i1_end

#ifdef CPP_MPI
      ! finally gather PNS_SO_ALL on master in case of MPI run
      allocate (work(lmmaxd,lmmaxd,irmd,2,natyp), stat=ierr)
      if (ierr/=0) stop 'Error allocating work for MPI comm of PNS_SO_ALL in main1a'
      ihelp = lmmaxd*lmmaxd*irmd*2*natyp
      call mpi_allreduce(pns_so_all,work, ihelp,mpi_double_complex,mpi_sum,         &
        t_mpi_c_grid%mympi_comm_ie,ierr)

      if (ierr/=mpi_success) stop 'Error in MPI comm of PNS_SO_ALL in main1a'
      pns_so_all(:, :, :, :, :) = work(:, :, :, :, :)
      deallocate (work, stat=ierr)
      if (ierr/=0) stop 'Error deallocating work for MPI comm of PNS_SO_ALL in main1a'
#endif

      ! done with preparations, call normcoeff routines that construct operators
      if (myrank==master) write (*, *) 'Computing spin operator'
      call normcoeff_so(natyp, t_params%ircut, lmmaxd/(1+korbit),pns_so_all,&
        t_params%thetas,t_params%ntcell,t_params%ifunm,t_params%ipan,t_params%lmsp, &
        t_inc%kvrel,t_params%cleb,t_params%icleb,t_params%iend,t_params%drdi,       &
        t_params%irws,1+korbit,0)

      if (myrank==master) write (*, *) 'Computing torq operator'
      call normcoeff_so_torq(natyp,t_params%ircut, lmmaxd/(1+korbit),       &
        pns_so_all,t_params%ntcell,t_params%ifunm,t_params%ipan,t_params%lmsp,      &
        t_inc%kvrel,t_params%cleb,t_params%icleb,t_params%iend,t_params%drdi,       &
        t_params%irws,t_params%visp,nspind,t_params%vins,t_params%irmin,0)

      if (myrank==master) write (*, *) 'Computing spinflux operator'
      call normcoeff_so_spinflux(natyp,t_params%ircut, lmmaxd/(1+korbit),   &
        pns_so_all,t_inc%kvrel,t_params%drdi,0)

    end if                         ! .not. impurity_operator_only

    !--------------------------------------------------------------------------------
    ! Part 2: operators for imp. wavefunctions
    !--------------------------------------------------------------------------------
    if (operator_imp) then

      if (t_inc%i_write>0) write (1337, *) 'Operators using impurity wavefunctions'
      if (t_inc%i_write>0) write (*, *) 'Operators using impurity wavefunctions'

      ! interpolate impurity wavefunctions to old radial mesh and global spin frame

      natomimp = t_imp%natomimp

      allocate (pns_so_imp(lmmaxd,lmmaxd,irmd,2,natomimp))
      pns_so_imp = czero

#ifdef CPP_MPI
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie,                       &
        t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,master,natomimp,ntot_pt,      &
        ioff_pt, .true.)

      i1_start  = ioff_pt(t_mpi_c_grid%myrank_ie) + 1
      i1_end    = ioff_pt(t_mpi_c_grid%myrank_ie) + ntot_pt(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1    = ntot_pt(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot_pt1 = ntot_pt
      t_mpi_c_grid%ioff_pt1 = ioff_pt
#else
      i1_start = 1
      i1_end = natomimp
#endif

      do i1_imp = i1_start, i1_end
        ! use I1 to point to host atom corresponding to position in imp. cluster
        ! (used to map to correct mesh)
        i1 = t_params%atomimp(i1_imp)

        allocate (rll(nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:0))
        allocate (pns_so(lmmaxd,lmmaxd,irmd,2))
        rll = czero
        pns_so = czero

        ! read wavefunctions ...
        rll(:, :, :, 0) = t_imp%rllimp(:, :, :, i1_imp)

        ! transform radial wavefunction back to old mesh
        allocate (rlltemp(t_inc%irmdnew,lmmaxd))
        allocate (pnstemp(irws(i1),lmmaxd))
        do lm1 = 1, lmmaxd
          rlltemp = czero
          pnstemp = czero
          ir = 0
          do ir = 1, t_inc%irmdnew
            do lm2 = 1, lmmaxd
              rlltemp(ir, lm2) = rll(lm1, lm2, ir, 0)
            end do
          end do
          call cheb2oldgrid(irws(i1),t_inc%irmdnew,lmmaxd,rmesh(:,i1),ncheb,        &
            npan_tot(i1),rpan_intervall(:,i1),ipan_intervall(:,i1),rlltemp,pnstemp, &
            irmd)
          do ir = 1, irws(i1)
            do lm2 = 1, lmmaxd
              pns_so(lm1, lm2, ir, 1) = pnstemp(ir, lm2)
            end do
          end do
        end do                     ! LM1
        ! for small component
        if (nsra==2) then
          do lm1 = 1, lmmaxd
            rlltemp = czero
            pnstemp = czero
            do ir = 1, t_inc%irmdnew
              do lm2 = 1, lmmaxd
                rlltemp(ir, lm2) = rll(lm1+lmmaxd, lm2, ir, 0)
              end do
            end do
            call cheb2oldgrid(irws(i1),t_inc%irmdnew,lmmaxd,rmesh(:,i1),ncheb,      &
              npan_tot(i1),rpan_intervall(:,i1),ipan_intervall(:,i1),rlltemp,       &
              pnstemp,irmd)
            do ir = 1, irws(i1)
              do lm2 = 1, lmmaxd
                pns_so(lm1, lm2, ir, 2) = pnstemp(ir, lm2)
              end do
            end do
          end do                   ! LM1
        end if                     ! NSRA.EQ.2
        ! rotate radial wavefunction to global frame
        do ir = 1, irmd
          call rotatematrix(pns_so(1,1,ir,1), t_imp%thetaimp(i1_imp), t_imp%phiimp(i1_imp), lmmaxd/2, 0)
          call rotatematrix(pns_so(1,1,ir,2), t_imp%thetaimp(i1_imp), t_imp%phiimp(i1_imp), lmmaxd/2, 0)
        end do
        ! finally collect wavefuncitons in global frame and old mesh for all atoms to be used in normcoeff-routines below
        pns_so_imp(:, :, :, :, i1_imp) = pns_so(:, :, :, :)
        deallocate (rll, pns_so)
        deallocate (rlltemp)
        deallocate (pnstemp)
      end do                       ! I1_imp=i1_start, i1_end

      ! deallocate temporary arrays
      deallocate (irws, rmesh, npan_tot, rpan_intervall, ipan_intervall, theta, phi)

#ifdef CPP_MPI
      ! finally gather PNS_SO_IMP on master in case of MPI run
      allocate (work(lmmaxd,lmmaxd,irmd,2,natomimp), stat=ierr)
      if (ierr/=0) stop 'Error allocating work for MPI comm of PNS_SO_ALL in main1a'
      ihelp = lmmaxd*lmmaxd*irmd*2*natomimp
      call mpi_allreduce(pns_so_imp,work,ihelp,mpi_double_complex,mpi_sum,          &
        t_mpi_c_grid%mympi_comm_ie,ierr)

      if (ierr/=mpi_success) stop 'Error in MPI comm of PNS_SO_ALL in main1a'
      pns_so_imp(:, :, :, :, :) = work(:, :, :, :, :)
      deallocate (work, stat=ierr)
      if (ierr/=0) stop 'Error deallocating work for MPI comm of PNS_SO_ALL in main1a'
#endif

      ! construct impurity operators using impurity wavefunctions
      if (myrank==master) write (*, *) 'Computing impurity spin operator'
      call normcoeff_so(natomimp,t_params%ircut,t_params%lmmaxd/2,pns_so_imp,       &
        t_params%thetas,t_params%ntcell,t_params%ifunm,t_params%ipan,t_params%lmsp, &
        t_inc%kvrel,t_params%cleb,t_params%icleb,t_params%iend,t_params%drdi,       &
        t_params%irws,1+korbit,1)

      if (myrank==master) write (*, *) 'Computing impurity torq operator'
      call normcoeff_so_torq(natomimp,t_params%ircut,t_params%lmmaxd/2,pns_so_imp,  &
        t_params%ntcell,t_params%ifunm,t_params%ipan,t_params%lmsp,t_inc%kvrel,     &
        t_params%cleb,t_params%icleb,t_params%iend,t_params%drdi,t_params%irws,     &
        t_imp%vispimp,nspind,t_imp%vinsimp,t_params%irmin,1)

      if (myrank==master) write (*, *) 'Computing impurity spinflux operator'
      call normcoeff_so_spinflux(natomimp,t_params%ircut,t_params%lmmaxd/2,         &
        pns_so_imp,t_inc%kvrel,t_params%drdi,1)

    end if                         ! operator_imp

    if (t_inc%i_write>0) write (1337, *) 'Done with Operators'
    if (t_inc%i_write>0) write (*, *) 'Done with Operators'

  end subroutine operators_for_fscode

end module mod_operators_for_fscode
