!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper module for the calculation of the T-matrix for the JM-KKR package
!> Author: Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
!> and many others ...
!> The code uses the information obtained in the main0 module, this is
!> mostly done via the `get_params_1a()` call, that obtains parameters of the type
!> `t_params` and passes them to local variables
!------------------------------------------------------------------------------------
module mod_main1a

  private
  public :: main1a

contains

  ! ----------------------------------------------------------------------------
  !> Summary: Main subroutine regarding the calculation of the t-matrix
  !> Author: Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
  !> and many others ...
  !> Category: single-site, potential, KKRhost
  !> Deprecated: False 
  !> Main subroutine for the calculation of the t-matrix
  !>
  !> Calls routines that compute singe-site wavefunctions and t-matrices.
  !> Two modes are impleneted for old (Born-iteration) solver without SOC or
  !> non-collinear magnetism and new solver (Chebychev mesh) working with larger
  !> spin-coupled matrices for SOC.
  !>
  !> @note
  !> PR: The BdG solver will only be impleneted with the newsolver for the moment.
  !> @endnote
  ! ----------------------------------------------------------------------------
  subroutine main1a()

#ifdef CPP_MPI
    use :: mpi
    use :: mod_types, only: gather_tmat, gather_lly_dtmat, save_t_mpi_c_grid, get_ntot_pt_ioff_pt_2d
    use :: mod_mympi, only: find_dims_2d
#endif
#ifdef CPP_TIMING
    use :: mod_timing, only: timing_start, timing_stop
#endif

    use :: mod_datatypes, only: dp
    use :: mod_constants, only: czero
    use :: mod_profiling, only: memocc
    use :: mod_tmatnewsolver, only: tmat_newsolver
    use :: mod_tbref, only: tbref
    use :: mod_getscratch, only: opendafile
    use :: mod_interpolate_poten, only: interpolate_poten
    use :: mod_initldau, only: initldau
    use :: mod_calctmat, only: calctmat
    use :: mod_types, only: t_tgmat, t_inc, t_lloyd, t_dtmatjij, init_t_dtmatjij, init_t_dtmatjij_at, t_mpi_c_grid
    use :: mod_mympi, only: nranks, master, myrank, distribute_work_atoms, distribute_work_energies
    use :: mod_wunfiles, only: get_params_1a, t_params, read_angles
    use :: mod_jijhelp, only: set_jijcalc_flags
    ! array dimensions
    use :: global_variables, only: natypd, wlength, lmmaxd, nrmaxd, lmpotd, nspotd, irmd, naclsd, nclsd, nrefd, ncleb, nembd, &
      naezd, lm2d, krel, nspind, iemxd, ntotd, nrmaxd, irmind, lmpotd, nspotd, npotd, natomimpd, ipand, knosph, lpotd, irnsd, korbit
#ifdef CPP_BdG
    use :: global_variables, only: mmaxd
#endif
    ! stuff defined in main0 already
    use :: mod_main0, only: ielast, nspin, icst, ipan, ircut, lmax, ncls, nineq, idoldau, lly, atom, cls, icleb, loflm, nacls, &
      refpot, irws, iend, ez, vins, irmin, alat, drdi, rmesh, zat, rcls, visp, rmtref, vref, cleb, cscl, socscale, socscl, erefldau, &
      ueff, jeff, solver, deltae, tolrdif, npan_log_at, npan_eq_at, ncheb, npan_tot, ipan_intervall, rpan_intervall, rnew, r_log, &
      ntldau, jwsrel, zrel, itscf, natomimp, atomimp, iqat, naez, natyp, nref, nsra, ins, itldau, lopt, vtrel, btrel, drdirel, &
      r2drdirel, rmrel, itrunldau, wldau, uldau, phildau

    implicit none

    ! .. Local variables
    integer :: i1
    integer :: ipot
    integer :: iltmp
    integer :: ispin
    integer :: itmpdir
    character (len=80) :: tmpdir
    logical :: lrefsys
    integer :: lrectmt
    integer :: lrectra
    ! .. Local arrays
    real (kind=dp), dimension (natypd) :: phi
    real (kind=dp), dimension (natypd) :: theta
    real (kind=dp), dimension (:, :, :), allocatable :: vinsnew

#ifdef CPP_MPI
    integer :: ntot1, mytot, ii
    integer, dimension (0:nranks-1) :: ntot_pt, ioff_pt, ntot_all, ioff_all
    ! communication of dtmat in case of lloyd
    integer :: iwork
    complex (kind=dp), dimension (:, :, :, :), allocatable :: work_jij
#endif
    integer :: i1_start, i1_end, ierr, i_stat, i_all

    logical, external :: opt, test
    ! ..
    ! data TOLRDIF /1.5D0/ ! Set free GF to zero if R<TOLRDIF in case of virtual atoms
    ! data LLY /0/

    lly = 0
    tolrdif = 1.5d0
    ! LRECTMT=WLENGTH*kind(czero)*LMMAXD*LMMAXD
    ! LRECTRA=WLENGTH*kind(czero)
    lrectmt = wlength*4*lmmaxd*lmmaxd
    lrectra = wlength*4

    allocate (vinsnew(nrmaxd,lmpotd,nspotd), stat=i_stat)
    call memocc(i_stat, product(shape(vinsnew))*kind(vinsnew), 'VINSNEW', 'main1a')
    vinsnew = 0.0d0

    ! Consistency check
    if ((krel<0) .or. (krel>1)) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
    if ((krel==1) .and. (nspind==2)) stop ' set NSPIN = 1 for KREL = 1 in the inputcard'
    ! -------------------------------------------------------------------------
    ! This routine previously used to read from unformatted files created by
    ! the main0 module, now  instead of unformatted files take parameters from
    ! types defined in wunfiles.F90
    ! -------------------------------------------------------------------------
    call get_params_1a(t_params,ipand,natypd,irmd,naclsd,ielast,nclsd,nrefd,ncleb,  &
      nembd,naezd,lm2d,nsra,ins,nspin,icst,ipan,ircut,lmax,ncls,nineq,idoldau,lly,  &
      krel,atom,cls,icleb,loflm,nacls,refpot,irws,iend,ez,vins,irmin,itmpdir,iltmp, &
      alat,drdi,rmesh,zat,rcls,iemxd,visp,rmtref,vref,cleb,cscl,socscale,socscl,    &
      erefldau,ueff,jeff,solver,tmpdir,deltae,tolrdif,npan_log_at,npan_eq_at,ncheb, &
      npan_tot,ipan_intervall,rpan_intervall,rnew,ntotd,nrmaxd,r_log,ntldau,itldau, &
      lopt,vtrel,btrel,drdirel,r2drdirel,rmrel,irmind,lmpotd,nspotd,npotd,jwsrel,   &
      zrel,itscf,natomimpd,natomimp,atomimp,iqat,naez,natyp,nref)

    if (test('Vspher  ')) vins(irmind:irmd, 2:lmpotd, 1:nspotd) = 0.d0

    ! -------------------------------------------------------------------------
    ! End read in variables
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    ! LDA+U treatment
    ! -------------------------------------------------------------------------
    if (idoldau==1) then
      open (67, file='ldau.unformatted', form='unformatted')
      read (67) itrunldau, wldau, uldau, phildau
      close (67)
      ! !---------------------------------------------------------------------
      ! Calculate Coulomb matrix ULDAU it calculates U matrix only once.
      ! Remove the next IF statement to have U calculated for each iteration anew.
      ! !---------------------------------------------------------------------

      ! IF ( ITRUNLDAU.LE.0 ) THEN
      call initldau(nsra, ntldau, itldau, lopt, ueff, jeff, erefldau, visp, nspin, rmesh, drdi, zat, ipan, ircut, phildau, uldau)
      ! END IF
    end if
    ! -------------------------------------------------------------------------
    ! End of LDA+U setup
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    ! No need to recalculate the reference system in SCF decimation case
    ! -------------------------------------------------------------------------
    ! ITSCF is initialised to 0 in main0
    lrefsys = .true.
    if (opt('DECIMATE') .and. (itscf>0)) lrefsys = .false.
    if (opt('rigid-ef') .and. (itscf>0)) lrefsys = .false.
    if (test('no-neutr') .and. (itscf>0)) lrefsys = .false.
    if (opt('no-neutr') .and. (itscf>0)) lrefsys = .false.
    if (test('lrefsysf') .or. opt('lrefsysf')) lrefsys = .false.


    if (t_tgmat%tmat_to_file) then
      call opendafile(69, 'tmat', 4, lrectmt, tmpdir, itmpdir, iltmp)
    end if

    if (lly/=0) then
      if (t_lloyd%dtmat_to_file) then
        call opendafile(691, 'dtmatde', 7, lrectmt, tmpdir, itmpdir, iltmp) ! LLY
      end if
      if (t_lloyd%tralpha_to_file) then
        call opendafile(692, 'tralpha', 7, lrectra, tmpdir, itmpdir, iltmp) ! LLY
      end if
    end if

    ! distribute atoms over ranks
    call distribute_work_atoms(natyp, i1_start, i1_end)

#ifdef CPP_MPI
    call mpi_barrier(mpi_comm_world, ierr)
#endif

    ! skip this part with GREENIMP option
    if (opt('GREENIMP') .or. test('IMP_ONLY')) then
      if (myrank==master) write (*, *) 'Skipping atom loop in main1a'
      i1_start = 1
      i1_end = 0
      ! distribute IE dimension here if atom loop is skipped
      ! otherwise this would be done in tmat_newsolver/calctmat
      call distribute_work_energies(ielast)
    end if

    if (.not. opt('NEWSOSOL')) then
      do i1 = i1_start, i1_end
        do ispin = 1, nspin
          ipot = nspin*(i1-1) + ispin

          call calctmat(icst,ins,ielast,nsra,ispin,nspin,i1,ez,drdi(1,i1),          &
            rmesh(1,i1),vins(irmind,1,knosph*ipot+(1-knosph)),visp(1,ipot),zat(i1), &
            irmin(i1),ipan(i1),ircut(0,i1),cleb,loflm,icleb,iend,solver,            &
            socscl(1,krel*i1+(1-krel)),cscl(1,krel*i1+(1-krel)),vtrel(1,i1),        &
            btrel(1,i1),rmrel(1,i1),drdirel(1,i1),r2drdirel(1,i1),zrel(i1),         &
            jwsrel(i1),idoldau,lopt(i1),wldau(1,1,1,i1),lly,deltae) ! LLY

        end do
      end do

    else ! opt('NEWSOSOL')

      ! !---------------------------------------------------------------------
      ! For calculation of Jij-tensor: create array for additional t-matrices and
      ! set atom-dependent flags which indicate if t-matrix is needed
      ! !---------------------------------------------------------------------
      call init_t_dtmatjij(t_inc, t_dtmatjij)
      if (opt('XCPL    ')) then
        call set_jijcalc_flags(t_dtmatjij, natyp, natomimpd, natomimp, atomimp, iqat)
      end if                       ! OPT('XCPL')

      ! nonco angles: defined in mod_wunfiles
      call read_angles(t_params, natyp, theta, phi)

      ! Interpolate potential
      call interpolate_poten(lpotd,irmd,irnsd,natyp,ipand,lmpotd,nspotd,ntotd,      &
        ntotd*(ncheb+1),nspin,rmesh,irmin,irws,ircut,vins,visp,npan_log_at,         &
        npan_eq_at,npan_tot,rnew,ipan_intervall,vinsnew)

      do i1 = i1_start, i1_end
        do ispin = 1, nspin/(1+korbit) ! run spin-loop only if 'NOSOC' test option is not used

          ipot = nspin*(i1-1) + ispin

#ifdef CPP_BdG
          if (test('BdG_dev ')) then
            call BdG_write_tmatnewsolver_inputs(nranks, i1, i1_start, ielast, &
              nspin, lmax, nsra, iend, lmpotd, lly, deltae, idoldau, ncleb, &
              ncheb, ntotd, mmaxd, nspind, iemxd, nrmaxd, nspotd, cleb, icleb, &
              ez, ipot, npan_tot, ipan_intervall, zat, phi, theta, &
              socscale, rnew, rpan_intervall, wldau, vinsnew, i1_end, natyp, lopt)
          end if
#endif

          call tmat_newsolver(ielast, nspin, lmax, zat(i1), socscale(i1), ez, nsra, cleb(:,1), icleb, iend, ncheb, npan_tot(i1), rpan_intervall(:,i1), ipan_intervall(:,i1), &
            rnew(:,i1), vinsnew, theta(i1), phi(i1), i1, ipot, lmpotd, lly, deltae, idoldau, lopt(i1), wldau(:,:,:,i1), t_dtmatjij(i1), ispin)

        end do ! ispin
      end do ! i1_start, i1_end atom loop

    end if ! NEWSOSOL

    if (idoldau==1) then
      open (67, file='ldau.unformatted', form='unformatted')
      write (67) itrunldau, wldau, uldau, phildau
      close (67)
    end if

    close (69)

    if (lly/=0) then
      if (t_lloyd%dtmat_to_file) close (691)
      if (t_lloyd%tralpha_to_file) close (692)
    end if


#ifdef CPP_MPI
    ! skip this part with GREENIMP option
    if (.not. (opt('GREENIMP') .or. test('IMP_ONLY'))) then

      if (.not. t_tgmat%tmat_to_file) then
        do ii = 0, t_mpi_c_grid%nranks_ie - 1
          ntot_all(ii) = t_mpi_c_grid%ntot_pt1(ii)
          ioff_all(ii) = t_mpi_c_grid%ioff_pt1(ii)
        end do
        mytot = t_mpi_c_grid%ntot_pt1(t_mpi_c_grid%myrank_ie)
        call gather_tmat(t_inc,t_tgmat,t_mpi_c_grid,ntot_all,ioff_all,mytot,        &
          t_mpi_c_grid%mympi_comm_ie,t_mpi_c_grid%nranks_ie)
      end if

      if (lly/=0 .and. .not. t_lloyd%dtmat_to_file) then
        if (t_mpi_c_grid%myrank_ie>(t_mpi_c_grid%dims(1)-1)) then
          ! reset tralpha and dtmat to zero for rest-ranks, otherwise
          ! these contributions are counted twice
          t_lloyd%tralpha = czero
          t_lloyd%dtmat = czero
        end if
        call gather_lly_dtmat(t_mpi_c_grid,t_lloyd,lmmaxd,t_mpi_c_grid%mympi_comm_ie)
      end if

      ! -------------------------------------------------------------------------
      ! for calculation of Jij-tensor
      ! -------------------------------------------------------------------------
      if (opt('XCPL    ') .and. opt('NEWSOSOL')) then
        do i1 = 1, t_inc%natyp
          ! initialize t_dtmatJij on other tasks
          ! t_dtmatJij was already allocated for certain atoms within the atom loop
          ! (in tmat_newsolver). This initialization cannot be made before tmat_newsolver,
          ! because the division of the enegry loop (done in there) influences t_dtmatJij.
          call init_t_dtmatjij_at(t_inc, t_mpi_c_grid, t_dtmatjij(i1))

          ! communicate
          if (t_dtmatjij(i1)%calculate) then
            iwork = product(shape(t_dtmatjij(i1)%dtmat_xyz))
            allocate (work_jij(iwork,1,1,1), stat=i_stat)
            call memocc(i_stat, product(shape(work_jij))*kind(work_jij), 'work_jij', 'main1a')

            call mpi_allreduce(t_dtmatjij(i1)%dtmat_xyz,work_jij,iwork,             &
              mpi_double_complex, mpi_sum, t_mpi_c_grid%mympi_comm_ie, ierr)
            if (ierr/=mpi_success) stop 'error communicating t_dtmatJij'
            call zcopy(iwork, work_jij, 1, t_dtmatjij(i1)%dtmat_xyz, 1)
            i_all = -product(shape(work_jij))*kind(work_jij)
            deallocate (work_jij, stat=i_stat)
            call memocc(i_stat, i_all, 'work_jij', 'main1a')
          end if                   ! t_dtmatJij(I1)%calculate

        end do                     ! I1=1,t_inc%NATYP

      end if                       ! OPT('XCPL    ').and.OPT('NEWSOSOL')

    end if                         ! .not.opt('GREENIMP')
    ! end skip this part with GREENIMP option
#endif
    ! -------------------------------------------------------------------------
    ! End of calculation of Jij-tensor
    ! -------------------------------------------------------------------------


#ifdef CPP_TIMING
    call timing_start('main1a - tbref')
#endif
    if (lrefsys) then
      call tbref(ez, ielast,alat,vref,iend,lmax,ncls,nineq,nref,cleb,rcls,atom,cls, &
        icleb,loflm,nacls,refpot,rmtref,tolrdif,tmpdir,itmpdir,iltmp,naez,lly) ! LLY Lloyd
    end if
#ifdef CPP_TIMING
    call timing_stop('main1a - tbref')
#endif

    if (t_inc%i_write>0) write (1337, '(79("="),/,30X,"< KKR1a finished >",/,79("="),/)')

    ! Deallocate leftover arrays
    if (allocated(vinsnew)) then
      i_all = -product(shape(vinsnew))*kind(vinsnew)
      deallocate (vinsnew, stat=i_stat)
      call memocc(i_stat, i_all, 'VINSNEW', 'main1a')
    end if

  end subroutine main1a


#ifdef CPP_BdG
  !-------------------------------------------------------------------------------
  !> Summary: Write out inputs for tmat_newsolver (BdG develop) 
  !> Author: Philipp Ruessmann  
  !> Deprecated: False 
  !> Category: input-output, unit-test, KKRhost
  !>
  !> @note JC: not sure how this exactly works bur variables do not seem to need 
  !> declaration before being run. Maybe it is a good idea to add it.                    
  !> @endnote                                                                       
  !-------------------------------------------------------------------------------
  subroutine BdG_write_tmatnewsolver_inputs(nranks, i1, i1_start, ielast, &
    nspin, lmax, nsra, iend, lmpotd, lly, deltae, idoldau, ncleb, &
    ncheb, ntotd, mmaxd, nspind, iemxd, nrmaxd, nspotd, cleb, icleb, &
    ez, ipot, npan_tot, ipan_intervall, zat, phi, theta, &
    socscale, rnew, rpan_intervall, wldau, vinsnew, i1_end, natyp, lopt)
    ! write out inputs for tmat_newsolver to extract first BdG
    use mod_datatypes, only: dp
    use :: global_variables, only: nchebd
    implicit none
    integer, intent(in) :: nranks, i1, i1_start, ielast, nspin ,lmax, ncleb, nsra, natyp, iend, lmpotd, lly, idoldau, ncheb, mmaxd, nspind, iemxd, nrmaxd, nspotd, ipot, i1_end
    integer, intent(in) :: ntotd, icleb(ncleb,4), lopt(natyp), ipan_intervall(0:ntotd,natyp), npan_tot(natyp)
    real(kind=dp), intent(in) :: cleb(ncleb,2), zat(natyp), phi(natyp), theta(natyp), socscale(natyp), rnew(ntotd*(nchebd+1),natyp), rpan_intervall(0:ntotd,natyp)
    real(kind=dp), intent(in) :: wldau(:,:,:,:), vinsnew(:,:,:)
    complex(kind=dp), intent(in) :: deltae, ez(iemxd)
    if (nranks>1) stop 'test option BdG_dev can only be used in serial!'
    if (i1==i1_start) open (887766, file='BdG_tmat_inputs.txt', form='formatted')
    if (i1==1) then
      write (887766, '(A25)') 'global parameters:'
      write (887766, *)
      write (887766, '(A25,I9)') 'IELAST= ', ielast
      write (887766, '(A25,I9)') 'NSPIN= ', nspin
      write (887766, '(A25,I9)') 'LMAX= ', lmax
      write (887766, '(A25,I9)') 'NSRA= ', nsra
      write (887766, '(A25,I9)') 'IEND= ', iend
      write (887766, '(A25,I9)') 'LMPOTD= ', lmpotd
      write (887766, '(A25,I9)') 'LLY= ', lly
      write (887766, '(A25,2ES21.9)') 'DELTAE= ', deltae
      write (887766, '(A25,I9)') 'IDOLDAU= ', idoldau
      write (887766, '(A25,I9)') 'NCLEB= ', ncleb
      write (887766, '(A25,I9)') 'NCHEB= ', ncheb
      write (887766, '(A25,I9)') 'NTOTD= ', ntotd
      write (887766, '(A25,I9)') 'MMAXD= ', mmaxd
      write (887766, '(A25,I9)') 'NSPIND= ', nspind
      write (887766, '(A25,I9)') 'IEMXD= ', iemxd
      write (887766, '(A25,I9)') 'NRMAXD= ', nrmaxd
      write (887766, '(A25,I9)') 'NSPOTD= ', nspotd
      write (887766, '(A25)') 'CLEB= '
      write (887766, '(999999999ES21.9)') cleb(:, 1)
      write (887766, '(A25)') 'ICLEB= '
      write (887766, '(999999999I9)') icleb(:, :)
      write (887766, '(A25)') 'EZ= '
      write (887766, '(2ES21.9)') ez
      write (887766, *)
      write (887766, '(A25)') 'atom-dependent input:'
      write (887766, *)
    end if
    write (887766, '(A25,I9)') 'I1= ', i1
    write (887766, '(A25,I9)') 'IPOT= ', ipot
    write (887766, '(A25,I9)') 'NPAN_TOT= ', npan_tot(i1)
    write (887766, '(A25,I9)') 'LPOT= ', lopt(i1)
    write (887766, '(A25,999999999I9)') 'IPAN_INTERVALL= ', ipan_intervall(:, i1)
    write (887766, '(A25,ES21.9)') 'ZAT= ', zat(i1)
    write (887766, '(A25,ES21.9)') 'PHI= ', phi(i1)
    write (887766, '(A25,ES21.9)') 'THETA= ', theta(i1)
    write (887766, '(A25,ES21.9)') 'SOCSCALE= ', socscale(i1)
    write (887766, '(A25,999999999ES21.9)') 'RNEW= ', rnew(:, i1)
    write (887766, '(A25,999999999ES21.9)') 'RPAN_INTERVALL= ', rpan_intervall(:, i1)
    write (887766, '(A25,999999999ES21.9)') 'WLDAU= ', wldau(:, :, :, i1)
    write (887766, '(A25,999999999ES21.9)') 'VINSNEW= ', vinsnew
    write (887766, *)
    if (i1==i1_end) then
      close (887766)
      stop 'done writing tmat_newsolver input of test option BdG_dev'
    end if
  end subroutine bdg_write_tmatnewsolver_inputs
#endif

end module mod_main1a
