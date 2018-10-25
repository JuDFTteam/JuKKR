!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculate the reference system for the decimation case
!> Author: 
!> Calculate the reference system for the decimation case. 
!------------------------------------------------------------------------------------
module mod_tbref

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the reference system for the decimation case
  !> Author: 
  !> Category: reference-system, structural-greensfunction, single-site, KKRhost
  !> Deprecated: False 
  !> Calculate the reference system for the decimation case. 
  !-------------------------------------------------------------------------------
  subroutine tbref(ez,ielast,alatc,vref,iend,lmax,ncls,nineq,nref,cleb,rcls,atom,   &
    cls,icleb,loflm,nacls,refpot,rmtref,tolrdif,tmpdir,itmpdir,iltmp,naez,lly)

    use :: mod_mympi, only: myrank, nranks, master, distribute_work_atoms
    use :: mod_types, only: t_tgmat, t_lloyd, t_inc
    use :: mod_types, only: t_mpi_c_grid, init_tgmat, init_tlloyd
#ifdef CPP_MPI
    use :: mpi
    use :: mod_mympi, only: find_dims_2d
#endif
    use :: mod_constants
    use :: mod_profiling
    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_calctref13
    use :: mod_getscratch, only: opendafile
    use :: mod_gll13, only: gll13

    implicit none

    ! .. Input variables
    integer, intent (in) :: lly    !! LLY <> 0 : alpply Lloyds formula
    integer, intent (in) :: iend
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: ncls   !! Number of reference clusters
    integer, intent (in) :: nref   !! Number of diff. ref. potentials
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: nineq  !! Number of ineq. positions in unit cell
    integer, intent (in) :: ielast
    real (kind=dp), intent (in) :: alatc
    real (kind=dp), intent (in) :: tolrdif !! For distance between scattering-centers smaller than [<TOLRDIF>], free GF is set to zero. Units are Bohr radii.
    integer, dimension (naezd+nembd), intent (in) :: cls !! Cluster around atomic sites
    integer, dimension (lm2d), intent (in) :: loflm !! l of lm=(l,m) (GAUNT)
    integer, dimension (nclsd), intent (in) :: nacls !! Number of atoms in cluster
    integer, dimension (naezd+nembd), intent (in) :: refpot !! Ref. pot. card  at position
    integer, dimension (naclsd, naezd+nembd), intent (in) :: atom !! Atom at site in cluster
    integer, dimension (ncleb, 4), intent (in) :: icleb !! Pointer array
    real (kind=dp), dimension (nref), intent (in) :: vref
    real (kind=dp), dimension (nref), intent (in) :: rmtref !! Muffin-tin radius of reference system
    real (kind=dp), dimension (ncleb, 2), intent (in) :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (3, naclsd, nclsd), intent (in) :: rcls !! Real space position of atom in cluster
    ! .. In/Out variables
    integer, intent (inout) :: iltmp
    integer, intent (inout) :: itmpdir
    character (len=80), intent (inout) :: tmpdir
    complex (kind=dp), dimension (iemxd), intent (inout) :: ez
    ! .. Local variables
    integer :: i1, ic, icls, ie, lm1, naclsmax, lrecgrf1, i_stat, i_all
    complex (kind=dp) :: eryd
    complex (kind=dp) :: lly_g0tr_ie ! LLY
    complex (kind=dp) :: lly_g0tr_dum ! LLY dummy variable used if no LLY is chosen to save memory
    ! .. Parameters
    integer :: lrecgrf
    ! .. Local Arrays
    complex (kind=dp), dimension (0:lmax, nref) :: alpharef ! LLY Lloyd Alpha matrix
    complex (kind=dp), dimension (0:lmax, nref) :: dalpharef ! LLY Derivative of the Lloyd Alpha matrix
    complex (kind=dp), dimension (lmgf0d, lmgf0d, nref) :: trefll ! LLY
    complex (kind=dp), dimension (lmgf0d, lmgf0d, nref) :: dtrefll ! LLY

    ! .. Local allocatable arrays
    complex (kind=dp), dimension (:, :), allocatable :: lly_g0tr ! LLY
    complex (kind=dp), dimension (:, :), allocatable :: dginp_dum ! LLY dummy variable used if no LLY is chosen to save memory
    complex (kind=dp), dimension (:, :, :), allocatable :: ginp
    complex (kind=dp), dimension (:, :, :), allocatable :: dginp
#ifdef CPP_MPI
    ! ..
    ! .. MPI variables
    integer :: ntot1, idim
    integer, dimension (0:nranks-1) :: ntot_pt, ioff_pt
    complex (kind=dp), dimension (:, :, :), allocatable :: work
#endif
    integer :: ie_start, ie_end
    integer :: i1_start, i1_end
    integer :: ie_num, ierr
    ! ..
    ! .. External Functions ..
    logical :: test, opt
    external :: test, opt

    naclsmax = 1
    do ic = 1, ncls
      if (nacls(ic)>naclsmax) naclsmax = nacls(ic)
    end do
    lrecgrf1 = wlength*4*naclsmax*lmgf0d*lmgf0d*ncls
    lrecgrf = wlength*4*naclsd*lmgf0d*lmgf0d*nclsd

    ! allocate and initialize ginp
    allocate (ginp(naclsmax*lmgf0d,lmgf0d,ncls), stat=i_stat)
    call memocc(i_stat, product(shape(ginp))*kind(ginp), 'GINP', 'TBREF')
    ginp = czero
    allocate (dginp_dum(naclsmax*lmgf0d,lmgf0d), stat=i_stat)
    call memocc(i_stat, product(shape(dginp_dum))*kind(dginp_dum), 'dginp_dum', 'TBREF')
    dginp_dum = czero

    if (lly/=0) then
      ! allocate and initialize dginp and lly_g0tr
      allocate (dginp(naclsmax*lmgf0d,lmgf0d,ncls), stat=i_stat)
      call memocc(i_stat, product(shape(dginp))*kind(dginp), 'DGINP', 'TBREF')
      dginp = czero
      allocate (lly_g0tr(ielast,nclsd), stat=i_stat)
      call memocc(i_stat, product(shape(lly_g0tr))*kind(lly_g0tr), 'LLY_G0TR', 'TBREF')
      lly_g0tr = czero
    end if

    if (t_tgmat%gref_to_file) then
      call opendafile(68, 'gref', 4, lrecgrf1, tmpdir, itmpdir, iltmp)
    end if
    if (lly/=0) then
      if (t_lloyd%dgref_to_file) then
        call opendafile(681, 'dgrefde', 7, lrecgrf1, tmpdir, itmpdir, iltmp)
      end if
      if (t_lloyd%g0tr_to_file) then
        open (682, file='lly_g0tr_ie.ascii', form='FORMATTED')
      end if
    end if

    ! ----------------------------------------------------------------------------
    if (t_inc%i_write>0) write (1337, *) myrank, 'start tbref e-loop'
    ! get share of energy loop done by each rank
#ifdef CPP_MPI
    ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
    ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
    ie_start = 0
    ie_end = ielast
#endif

    ! now initialize arrays for tmat, gmat, and gref
    call init_tgmat(t_inc, t_tgmat, t_mpi_c_grid)
    if (lly/=0) call init_tlloyd(t_inc, t_lloyd, t_mpi_c_grid)

    ! start energy loop
    do ie_num = 1, ie_end
      ie = ie_start + ie_num
      if (t_inc%i_write>0) write (1337, fmt='(A23,I5,2F14.8)') 'TBREF: GREF for energy:', ie, ez(ie)

      ! reset arrays to zero (needed in case of MPIatom energy/atom distribution
      ginp(:, :, :) = czero
      if (lly>0) then
        dginp(:, :, :) = czero
        lly_g0tr(ie, :) = czero
        dginp_dum = czero
        lly_g0tr_dum = czero
      end if

      ! skip this part with GREENIMP option
      if (opt('GREENIMP')) then
        if (myrank==master) write (*, *) 'Skipping atom loop in tbref for energy point', ie, 'of', ielast
        i1_start = 1
        i1_end = 0
      else
        i1_start = 1
        i1_end = nref
      end if

      eryd = ez(ie)

      do i1 = i1_start, i1_end     ! 1,NREF
        call calctref13(eryd, vref(i1), rmtref(i1), lmax, lm1 & ! LLY Lloyd
          , trefll(1,1,i1), dtrefll(1,1,i1), & ! LLY Lloyd
          alpharef(0,i1), dalpharef(0,i1), lmax+1, lmgf0d) ! LLY Lloyd

        if (test('flow    ') .and. (t_inc%i_write>0)) write (1337, fmt=*) 'tll(ref),  i1 = ', i1
      end do

      if (test('flow    ') .and. (t_inc%i_write>0)) write (1337, fmt=*) 't-mat(Ref) o.k.', ie
      ! ----------------------------------------------------------------------
      ! distribute NCLS loop over atom ranks
      call distribute_work_atoms(ncls, i1_start, i1_end, distribute_rest=.false.)

      ! overwrite  to skip this part with GREENIMP option
      if (opt('GREENIMP')) then
        i1_start = 1
        i1_end = 0
      end if

      do icls = i1_start, i1_end
        i1 = 1
        ic = 0
        do while (ic==0 .and. i1<=nineq)
          if (cls(i1)==icls) ic = i1
          i1 = i1 + 1
        end do

        if (ic==0) stop 'Error in CLS(*) array in tbref'
        if (test('flow    ') .and. (t_inc%i_write>0)) then
          write (1337, fmt=*) 'CLUSTER ', icls, ' at ATOM ', ic
        end if

        call gll13(eryd,cleb(1,2),icleb,loflm,iend,trefll,dtrefll,atom(1,ic),refpot,&
          rcls(1,1,icls),nacls(icls),tolrdif,alatc,0,ginp(1,1,icls),dginp_dum,      &
          naclsmax,lly_g0tr_dum,lly)

        ! copy from dummy variable to arrays
        if (lly/=0) then
          dginp(:, :, icls) = dginp_dum(:, :)
          lly_g0tr(ie, icls) = lly_g0tr_dum
        end if

      end do                       ! icls


#ifdef CPP_MPI
      if (t_mpi_c_grid%nranks_ie>1) then
        ! gather results of ginp, dginp and lly_g0tr from above parallel loop
        idim = naclsmax*lmgf0d*lmgf0d*ncls
        allocate (work(naclsmax*lmgf0d,lmgf0d,ncls), stat=i_stat)
        call memocc(i_stat, product(shape(work))*kind(work), 'work', 'TBREF')

        call mpi_allreduce(ginp,work,idim,mpi_double_complex,mpi_sum,               &
          t_mpi_c_grid%mympi_comm_ie,ierr)
        call zcopy(idim, work, 1, ginp, 1)

        i_all = -product(shape(work))*kind(work)
        deallocate (work, stat=i_stat)
        call memocc(i_stat, i_all, 'work', 'TBREF')
        if (lly/=0) then
          idim = naclsmax*lmgf0d*lmgf0d*ncls
          allocate (work(naclsmax*lmgf0d,lmgf0d,ncls), stat=i_stat)
          call memocc(i_stat, product(shape(work))*kind(work), 'work', 'TBREF')

          call mpi_allreduce(dginp,work,idim,mpi_double_complex,mpi_sum,            &
            t_mpi_c_grid%mympi_comm_ie,ierr)
          call zcopy(idim, work, 1, dginp, 1)
          i_all = -product(shape(work))*kind(work)
          deallocate (work, stat=i_stat)
          call memocc(i_stat, i_all, 'work', 'TBREF')
          idim = ncls
          allocate (work(nclsd,1,1), stat=i_stat)
          call memocc(i_stat, product(shape(work))*kind(work), 'work', 'TBREF')
          call mpi_allreduce(lly_g0tr(ie,:),work,idim,mpi_double_complex,mpi_sum,   &
            t_mpi_c_grid%mympi_comm_ie,ierr)
          call zcopy(idim, work, 1, lly_g0tr(ie,:), 1)
          i_all = -product(shape(work))*kind(work)
          deallocate (work, stat=i_stat)
          call memocc(i_stat, i_all, 'work', 'TBREF')
        end if
      end if                       ! t_mpi_c_grid%nranks_ie>1
#endif
      if (lly/=0) then             ! LLY Lloyd
        lly_g0tr_ie = czero        ! LLY Lloyd
        do i1 = 1, naez            ! LLY Lloyd
          icls = cls(i1)           ! LLY Lloyd
          lly_g0tr_ie = lly_g0tr_ie + lly_g0tr(ie, icls) ! LLY Lloyd
        end do                     ! LLY Lloyd
      end if                       ! LLY Lloyd
      ! ----------------------------------------------------------------------
      if (t_tgmat%gref_to_file) then
#ifdef CPP_MPI
        ! make sure only one processor writes for one energy point
        if (t_mpi_c_grid%myrank_ie==0) then
          write (68, rec=ie) ginp  ! Gref
        end if                     ! if(t_mpi_c_grid%myrank_ie==0) then

        ! human readable writeout if test option is hit
        if (test('fileverb')) then
          ! test writeout
          do i1 = 0, nranks - 1
            if (myrank==i1) then
              write (686868+myrank, *) myrank, ie
              write (686868+myrank, '(2ES21.9)') ginp
            end if
            call mpi_barrier(t_mpi_c_grid%mympi_comm_ie, ierr)
          end do
          ! end test writeout
        end if
#else
        write (68, rec=ie) ginp    ! Gref
#endif
      else                         ! (t_tgmat%gref_to_file)
        t_tgmat%gref(:, :, :, ie_num) = ginp(:, :, :)
      end if                       ! (t_tgmat%gref_to_file)
      ! store result either to file or in memory
      if (lly/=0) then             ! LLY Lloyd
        if (t_lloyd%dgref_to_file) then ! LLY Lloyd
          write (681, rec=ie) dginp ! dGref/dE         ! LLY Lloyd
          if (test('fileverb')) then
            write (681681, '(i9,200000ES15.7)') ie, dginp
          end if
        else                       ! LLY Lloyd
          t_lloyd%dgref(:, :, :, ie_num) = dginp(:, :, :) ! LLY Lloyd
        end if                     ! LLY Lloyd
        if (t_lloyd%g0tr_to_file) then ! LLY Lloyd
          write (682, fmt='(2E24.16)') lly_g0tr_ie ! LLY Lloyd
          if (test('fileverb')) then
            write (682682, '(i9,200000ES15.7)') ie_num, lly_g0tr_ie
          end if
        else                       ! LLY Lloyd
          t_lloyd%g0tr(ie_num) = lly_g0tr_ie ! LLY Lloyd
        end if                     ! LLY Lloyd
      end if                       ! LLY.NE.0                                 ! LLY Lloyd

      if (test('flow    ') .and. (t_inc%i_write>0)) write (1337, fmt=*) 'G(n,lm,n,lm) (Ref) o.k.'
    end do                         ! IE
    ! ----------------------------------------------------------------------------
    if (t_tgmat%gref_to_file) then
      close (68)
    end if

    if (lly/=0 .and. t_lloyd%dgref_to_file) close (681)
    if (lly/=0 .and. t_lloyd%g0tr_to_file) close (682)
    i_all = -product(shape(ginp))*kind(ginp)
    deallocate (ginp, stat=i_stat)
    call memocc(i_stat, i_all, 'GINP', 'TBREF')
    if (lly/=0) then
      i_all = -product(shape(dginp))*kind(dginp)
      deallocate (dginp, stat=i_stat)
      call memocc(i_stat, i_all, 'DGINP', 'TBREF')

      i_all = -product(shape(lly_g0tr))*kind(lly_g0tr)
      deallocate (lly_g0tr, stat=i_stat)
      call memocc(i_stat, i_all, 'LLY_G0TR', 'TBREF')
    end if

  end subroutine tbref

end module mod_tbref
