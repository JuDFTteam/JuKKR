module mod_rhoqtools
  use :: mod_datatypes, only: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Write the k-mesh to file with its respective weights
  !> Author: 
  !> Category: input-output, k-points, KKRhost 
  !> Deprecated: False 
  !> Write the k-mesh to file with its respective weights
  !-------------------------------------------------------------------------------
  subroutine rhoq_write_kmesh(nofks, nxyz, volbz, bzkp, volcub, recbv, bravais)

    implicit none

    integer, intent (in) :: nofks             !! number of points in irreducible BZ
    real (kind=dp), intent (in) :: volbz      !! volume of the BZ
    integer, dimension(3), intent (in) :: nxyz !! original k-mesh net in the 3 directions of the reciprocal lattice vectors (not xyz directions)
    real (kind=dp), dimension(nofks), intent (in) :: volcub       !! Weight of the k-points
    real (kind=dp), dimension(3, nofks), intent (in)  :: bzkp     !! k-point mesh
    real (kind=dp), dimension(3, 3), intent (in)      :: recbv    !! Reciprocal basis vectors
    real (kind=dp), dimension(3, 3), intent (in)      :: bravais  !! Bravais lattice vectors

    ! .. Local variables
    integer :: ks, i
    ! write out kpoints
    open (8888, file='kpts.txt', form='formatted')
    write (8888, '(3I9)') nofks, nxyz(1), nxyz(2)
    write (8888, '(E16.7)') volbz
    do ks = 1, nofks
      write (8888, '(4E16.7)')(bzkp(i,ks), i=1, 3), volcub(ks)
    end do
    write (8888, '(100E16.7)') recbv(1:3, 1:3), bravais(1:3, 1:3)
    close (8888)

  end subroutine rhoq_write_kmesh

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rhoq_read_mu0_scoef(iatomimp, mu, nscoef, imin)

#ifdef CPP_MPI
    use :: mpi
#endif
    use :: mod_mympi, only: myrank, master
    implicit none

    integer, intent (out) :: mu, nscoef, imin
    integer, allocatable, intent (inout) :: iatomimp(:)
    ! local
    integer :: i1
#ifdef CPP_MPI
    integer :: ierr
#endif

    ! read in mu0
    if (myrank==master) then
      open (8888, file='mu0', form='formatted')
      read (8888, *) mu, nscoef
      allocate (iatomimp(nscoef))
      do i1 = 1, nscoef
        read (8888, *) iatomimp(i1)
      end do
      close (8888)
    end if

#ifdef CPP_MPI
    ! communicate mu and nscoef
    call mpi_bcast(mu, 1, mpi_integer, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'Error Bcast mu0'
    call mpi_bcast(nscoef, 1, mpi_integer, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'Error Bcast nscoef'
    if (.not. allocated(iatomimp)) allocate (iatomimp(nscoef))
    call mpi_bcast(iatomimp, nscoef, mpi_integer, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'Error Bcast iatomimp'
#endif

    ! find imin
    imin = 100000
    do i1 = 1, nscoef
      if (iatomimp(i1)<imin) imin = iatomimp(i1)
    end do
    nscoef = nscoef - 1

  end subroutine rhoq_read_mu0_scoef

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rhoq_find_kmask(nofks, k_end, bzkp, kmask, rhoq_kmask)

#ifdef CPP_MPI
    use :: mpi
#endif
    use :: mod_mympi, only: myrank, master
    use :: mod_datatypes

    implicit none

    integer, intent (in) :: nofks
    integer, intent (out) :: k_end
    integer, allocatable, intent (out) :: kmask(:)
    real (kind=dp), intent (in) :: bzkp(3, nofks)
    real (kind=dp), allocatable, intent (out) :: rhoq_kmask(:, :)
    ! local
    integer :: i, j, kpt, kmask_mode, k_start
    logical :: kmask_info
    real (kind=dp) :: k_mask_bounds(4), recbv(3, 3), kp(3)
#ifdef CPP_MPI
    integer :: ierr
#endif

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (myrank==master) then
      ! read recbv
      open (8888, file='kpts.txt', form='formatted')
      read (8888, *)
      read (8888, *)
      do kpt = 1, nofks
        read (8888, *)
      end do
      read (8888, '(100E16.7)') recbv(1:3, 1:3)
      close (8888)

      allocate (kmask(nofks))

      ! determine kmask parameters
      inquire (file='kmask_info.txt', exist=kmask_info)
      if (kmask_info) then
        write (*, *) 'found ''kmask_info.txt'' file, start reading...'
        open (8888, file='kmask_info.txt', form='formatted')
        read (8888, *) kmask_mode
        if (kmask_mode==1) then    ! read R0=(x,y), then R1 and R2 (outer and inner radius around R0)
          write (*, *) 'kmask_mode is 1: spherical region'
          read (8888, *) k_mask_bounds(1), k_mask_bounds(2)
          write (*, *) 'R0=', k_mask_bounds(1), k_mask_bounds(2)
          read (8888, *) k_mask_bounds(3)
          write (*, *) 'R1=', k_mask_bounds(3)
          read (8888, *) k_mask_bounds(4)
          write (*, *) 'R2=', k_mask_bounds(4)
        else if (kmask_mode==2) then ! read xmin, xmax, ymin, ymax of kmask_box
          write (*, *) 'kmask_mode is 2: box'
          read (8888, *) k_mask_bounds(1)
          read (8888, *) k_mask_bounds(2)
          read (8888, *) k_mask_bounds(3)
          read (8888, *) k_mask_bounds(4)
          write (*, *) 'xmin=', k_mask_bounds(1)
          write (*, *) 'xmax=', k_mask_bounds(2)
          write (*, *) 'ymin=', k_mask_bounds(3)
          write (*, *) 'ymax=', k_mask_bounds(4)
        end if                     ! kmask_mode 1 or 2
        ! close kmask_info.txt
        close (8888)
      else
        kmask_mode = 0
      end if                       ! kmask_info.txt file found

      if (kmask_mode==3) then      ! read kmask from file
        write (*, *) 'kmask_mode is 3: read ''kpts_mask.txt'' file'
        open (8888, file='kpts_mask.txt', form='formatted')
      end if

      ! use these as counters
      k_start = 1
      k_end = 0

      ! find kmask and number points in box
      do kpt = 1, nofks
        ! findig kmask
        ! default is take all
        kmask(kpt) = 1
        if (kmask_mode==1) then    ! sph mode
          kmask(kpt) = 0
          do i = -1, 1, 1
            do j = -1, 1, 1
              kp(1:3) = bzkp(1:3, kpt) + i*recbv(1:3, 1) + j*recbv(1:3, 2)
              ! first shift kpt to be centered around R0
              kp(1) = kp(1) - k_mask_bounds(1)
              kp(2) = kp(2) - k_mask_bounds(2)
              ! then apply rules concerning inner and outer radius
              if (sqrt(kp(1)**2+kp(2)**2)<k_mask_bounds(3)) kmask(kpt) = 1
              if (sqrt(kp(1)**2+kp(2)**2)<k_mask_bounds(4)) kmask(kpt) = 0
            end do                 ! j
          end do                   ! i
        else if (kmask_mode==2) then ! box mode
          do i = -1, 1, 1
            do j = -1, 1, 1
              kp(1:3) = bzkp(1:3, kpt) + i*recbv(1:3, 1) + j*recbv(1:3, 2)
              if (kp(1)<k_mask_bounds(1)) kmask(kpt) = 0
              if (kp(1)>k_mask_bounds(2)) kmask(kpt) = 0
              if (kp(2)<k_mask_bounds(3)) kmask(kpt) = 0
              if (kp(2)>k_mask_bounds(4)) kmask(kpt) = 0
            end do                 ! j
          end do                   ! i
        else if (kmask_mode==3) then ! read kmask from file
          read (8888, *) kmask(kpt)
        end if                     ! kmask_mode
        ! count number of kpts in reduced part
        if (kmask(kpt)>0) k_end = k_end + 1
      end do                       ! kpt loop

      ! close kmask file
      if (kmask_mode==3) then      ! read kmask from file
        close (8888)
      end if

      ! fill rhoq_kmask (on reduced set of kpts)
      allocate (rhoq_kmask(1:5,k_end))
      do kpt = 1, nofks
        if (kmask(kpt)>0) then
          rhoq_kmask(1:3, k_start) = bzkp(1:3, kpt)
          rhoq_kmask(4, k_start) = real(kpt, kind=dp)
          rhoq_kmask(5, k_start) = real(kmask(kpt), kind=dp)
          k_start = k_start + 1
        end if
      end do

      write (*, *) 'found ', k_end, 'kpoints'

      open (8888, file='rhoq_kmask.test', form='formatted')
      do kpt = 1, k_end
        write (8888, '(5F14.7)') rhoq_kmask(1:3, kpt), rhoq_kmask(4, kpt), rhoq_kmask(5, kpt)
      end do
      close (8888)

    end if                         ! (myrank==master)

#ifdef CPP_MPI
    ! communicate kmask stuff from master to all others
    call mpi_bcast(k_end, 1, mpi_integer, master, mpi_comm_world, ierr)
    if (myrank/=master) then
      allocate (kmask(nofks))
      allocate (rhoq_kmask(1:5,k_end))
    end if
    call mpi_bcast(kmask, nofks, mpi_integer, master, mpi_comm_world, ierr)
    call mpi_bcast(rhoq_kmask, 5*k_end, mpi_double_precision, master, mpi_comm_world, ierr)
#endif

  end subroutine rhoq_find_kmask

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rhoq_saveg(nscoef, rhoq_kmask, kpt, k_end, kp, i, j, mu, imin, iatomimp, lmmaxd, g)

#ifdef CPP_HYBRID
    use :: omp_lib
#endif

    implicit none

    integer, intent (in) :: i, j, mu, imin, lmmaxd, nscoef, k_end, kpt
    integer, intent (in) :: iatomimp(nscoef)
    real (kind=dp), intent (in) :: rhoq_kmask(5, k_end), kp(3)
    complex (kind=dp), intent (in) :: g(lmmaxd, lmmaxd)
    ! local
    integer :: ix, jx, lm1, irec

#ifdef CPP_HYBRID
    ! $omp critical
#endif
    irec = (nscoef*2)*(int(rhoq_kmask(4,kpt))-1)
    if (((i==mu) .and. any(j==iatomimp(1:nscoef)))) then
      ix = 0
      jx = 0
      lm1 = 1
      do while (ix==0 .and. lm1<=nscoef)
        if (iatomimp(lm1)==j) ix = j - imin + 1
        lm1 = lm1 + 1
      end do
      irec = irec + nscoef + ix
      write (9889, rec=irec) kp(1:2), g(1:lmmaxd, 1:lmmaxd) ! *
      ! rhoq_kmask(5,kpt)
    end if

    irec = (nscoef*2)*(int(rhoq_kmask(4,kpt))-1)
    if (((j==mu) .and. any(i==iatomimp(1:nscoef)))) then
      ix = 0
      jx = 0
      lm1 = 1
      do while (jx==0 .and. lm1<=nscoef+1)
        if (iatomimp(lm1)==i) jx = i - imin + 1
        lm1 = lm1 + 1
      end do
      irec = irec + jx
      write (9889, rec=irec) kp(1:2), g(1:lmmaxd, 1:lmmaxd) ! * !*ETAIKR(ISYM,NS) not this factor because it is dealt with explicitly in rhoq module
      ! rhoq_kmask(5,kpt)
    end if                         ! i==mu ...
#ifdef CPP_HYBRID
    ! $omp end critical
#endif

  end subroutine rhoq_saveg

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rhoq_write_tau0(nofks, nshell, nsh1, nsh2, nsymat, nscoef, mu, iatomimp, kmask, lmmaxd, bzkp, imin)

    use :: mod_mympi, only: myrank, master
    use :: mod_datatypes
    use :: constants, only: czero

    implicit none

    integer, intent (in) :: nofks, nshell, nsymat, nscoef, mu, lmmaxd, imin
    integer, intent (in) :: nsh1(nshell), nsh2(nshell), iatomimp(nscoef)
    real (kind=dp), intent (in) :: bzkp(3, nofks)
    integer, allocatable, intent (inout) :: kmask(:)
    ! local
    integer :: kpt, ns, i, j, isym, irec, ix, jx, lm1
    real (kind=dp) :: kp(3)
    complex (kind=dp) :: g(lmmaxd, lmmaxd)

    if (myrank==master) then
      write (*, *)                 ! status bar
      write (*, *) 'rhoq: write-out loop'
      write (*, '("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
      write (*, fmt=100, advance='no') ! beginning of statusbar

      ! write out fort.998899, fort.998888
      ! ======================================================================
      do kpt = 1, nofks
        do ns = 1, nshell
          i = nsh1(ns)
          j = nsh2(ns)
          do isym = 1, nsymat
            irec = (nscoef*2)*(kpt-1)
            if (((i==mu) .and. any(j==iatomimp(1:nscoef)))) then
              ix = 0
              jx = 0
              lm1 = 1
              do while (ix==0 .and. lm1<=nscoef)
                if (iatomimp(lm1)==j) ix = j - imin + 1
                lm1 = lm1 + 1
              end do
              irec = irec + nscoef + ix
              if (kmask(kpt)>0) then
                read (9889, rec=irec) kp(1:2), g(1:lmmaxd, 1:lmmaxd)
              else
                kp(1:3) = bzkp(1:3, kpt)
                g(1:lmmaxd, 1:lmmaxd) = czero
              end if
              write (998899, '(10000ES15.7)') kp(1:2), g(1:lmmaxd, 1:lmmaxd)*real(kmask(kpt), kind=dp)
            end if
            irec = (nscoef*2)*(kpt-1)
            if (((j==mu) .and. any(i==iatomimp(1:nscoef)))) then
              ix = 0
              jx = 0
              lm1 = 1
              do while (jx==0 .and. lm1<=nscoef+1)
                if (iatomimp(lm1)==i) jx = i - imin + 1
                lm1 = lm1 + 1
              end do
              irec = irec + jx
              if (kmask(kpt)>0) then
                read (9889, rec=irec) kp(1:2), g(1:lmmaxd, 1:lmmaxd)
              else
                kp(1:3) = bzkp(1:3, kpt)
                g(1:lmmaxd, 1:lmmaxd) = czero
              end if
              write (998888, '(10000ES15.7)') kp(1:2), g(1:lmmaxd, 1:lmmaxd)*real(kmask(kpt), kind=dp)
            end if                 ! iii==mu ...
          end do                   ! ISYM = 1,NSYMAT
        end do                     ! NS = 1,NSHELL

        if (nofks>=50) then
          if (mod(kpt-0,nofks/50)==0) write (6, fmt=110, advance='no')
        else
          write (6, fmt=110, advance='no')
        end if

      end do                       ! kpt=1,nofks
      ! ======================================================================
      write (6, *)                 ! status bar
    end if                         ! myrank==master
    deallocate (kmask)

100 format ('                 |')  ! status bar
110 format ('|')                   ! status bar

  end subroutine rhoq_write_tau0

  !-------------------------------------------------------------------------------
  !> Summary: Write the radial mesh information to file
  !> Author: 
  !> Category: radial-grid, input-output, KKRhost
  !> Deprecated: False
  !> Write the radial mesh information to file
  !-------------------------------------------------------------------------------
  subroutine rhoq_save_rmesh(natyp,irmd,ipand,irmin,irws,ipan,rmesh,ntcell,ircut,   &
    r_log,npan_log,npan_eq)

    implicit none

    integer, intent (in) :: natyp, irmd, ipand, npan_log, npan_eq
    integer, intent (in) :: irmin(natyp), irws(natyp), ipan(natyp), ntcell(natyp), ircut(0:ipand, natyp)
    real (kind=dp), intent (in) :: r_log
    real (kind=dp), intent (in) :: rmesh(irmd, natyp)
    ! local
    integer :: i1

    ! read mu_0
    open (9999, file='mu0', form='formatted')
    read (9999, *) i1
    close (9999)
    ! write out corresponding mesh information
    open (9999, file='rmesh_mu0.txt', form='formatted')
    write (9999, '(A)') '# mu_0, IRMIN, IRWS, IPAN'
    write (9999, '(4I9)') i1, irmin(i1), irws(i1), ipan(i1)
    write (9999, '(A)') '# Rmesh(1:IRWS)'
    write (9999, '(1000E22.15)') rmesh(1:irws(i1), i1)
    write (9999, '(A)') '# NTCELL(1:NATYP)'
    write (9999, '(1000I9)') ntcell(1:natyp)
    write (9999, '(A)') '# IRCUT(1:IPAN)'
    write (9999, '(1000I9)') ircut(1:ipan(i1), i1)
    write (9999, '(A)') '# R_LOG, NPAN_LOG, NPAN_EQ'
    write (9999, '(E22.15,2I9)') r_log, npan_log, npan_eq
    close (9999)

  end subroutine rhoq_save_rmesh

  !-------------------------------------------------------------------------------
  !> Summary: Save the reference potentials to an unformatted file
  !> Author: 
  !> Category: input-output, reference-system, KKRhost
  !> Deprecated: False 
  !> Save the reference potentials to an unformatted file
  !-------------------------------------------------------------------------------
  subroutine rhoq_save_refpot(ielast,i1,nref,natyp,refpot,wlength,lmmaxd,ie,trefll)

    implicit none

    integer, intent (in) :: ielast, i1, nref, natyp, wlength, lmmaxd, ie
    integer, intent (in) :: refpot(natyp)
    complex (kind=dp), intent (in) :: trefll(lmmaxd, lmmaxd, natyp)
    ! local
    integer :: irec

    if (i1==1) then
      open (99991, file='refinfo')
      write (99991, '(2I9)') nref, natyp
      write (99991, '(1000I9)') refpot(1:natyp)
      close (99991)
      open (99992, file='tref', access='direct', recl=wlength*4*lmmaxd**2, form='unformatted')
    end if
    irec = ie + ielast*(i1-1)
    write (99992, rec=irec) trefll(:, :, i1)

  end subroutine rhoq_save_refpot

end module mod_rhoqtools
