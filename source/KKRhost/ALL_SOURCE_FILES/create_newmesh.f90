! -------------------------------------------------------------------------------
! MODULE: mod_create_newmesh
! > @brief Module for the creation of the integration grid for the new solver
! -------------------------------------------------------------------------------
module mod_create_newmesh

  use :: constants
  use :: profiling
  use :: mod_datatypes, only: dp

  implicit none
  private :: dp

contains

  ! ----------------------------------------------------------------------------
  ! SUBROUTINE: CREATE_NEWMESH
  ! > @brief Creation of the integration grid for the new solver
  ! > @note changed interface to get rid of inc.p and to be able to use i
  ! > create_newmesh in tmatimp routine for GREENIMP option
  ! > this is the list of  array dimensions previously importted from inc.p
  ! ----------------------------------------------------------------------------
  subroutine create_newmesh(natyp, lmax, lpot, irm, irnsd, ipand, irid, ntotd, &
    nfund, ncheb, irmdnew, nspin, r, irmin, ipan, ircut, r_log, npan_log, &
    npan_eq, npan_log_at, npan_eq_at, npan_tot, rnew, rpan_intervall, &
    ipan_intervall, ncelld, ntcell, thetas, thetasnew) ! < optional arguments

    implicit none

    ! .. Input variables
    integer, intent (in) :: irm    ! < Maximum number of radial points
    integer, intent (in) :: irid   ! < Shape functions parameters in
                                   ! non-spherical part
    integer, intent (in) :: lmax   ! < Maximum l component in wave function
                                   ! expansion
    integer, intent (in) :: lpot   ! < Maximum l component in potential
                                   ! expansion
    integer, intent (in) :: nspin  ! < Counter for spin directions
    integer, intent (in) :: natyp  ! < Number of kinds of atoms in unit cell
    integer, intent (in) :: irnsd
    integer, intent (in) :: ipand  ! < Number of panels in non-spherical part
    integer, intent (in) :: ntotd
    integer, intent (in) :: nfund  ! < Shape functions parameters in
                                   ! non-spherical part
    integer, intent (in) :: ncheb  ! < Number of Chebychev pannels for the new
                                   ! solver
    integer, intent (in) :: ncelld ! < Number of cells (shapes) in
                                   ! non-spherical part
    integer, intent (in) :: irmdnew
    integer, intent (in) :: npan_eq ! < Variables for the pannels for the new
                                    ! solver
    integer, intent (in) :: npan_log ! < Variables for the pannels for the new
                                     ! solver
    real (kind=dp), intent (in) :: r_log
    integer, dimension (natyp), intent (in) :: ipan ! < Number of panels in
                                                    ! non-MT-region
    integer, dimension (natyp), intent (in) :: irmin ! < Max R for spherical
                                                     ! treatment
    integer, dimension (0:ipand, natyp), intent (in) :: ircut ! < R points of
                                                              ! panel borders
    real (kind=dp), dimension (irm, natyp), intent (in) :: r ! < Radial mesh (
                                                             ! in units a
                                                             ! Bohr)
    ! .. Input/Output variables
    integer, dimension (natyp), intent (inout) :: npan_tot
    integer, dimension (natyp), intent (inout) :: npan_eq_at
    integer, dimension (natyp), intent (inout) :: npan_log_at
    integer, dimension (0:ntotd, natyp), intent (inout) :: ipan_intervall
    real (kind=dp), dimension (irmdnew, natyp), intent (inout) :: rnew
    real (kind=dp), dimension (0:ntotd, natyp), intent (inout) :: &
      rpan_intervall
    ! .. Optional arguments, do interpolation when given
    integer, dimension (natyp), intent (in), optional :: ntcell ! < Index for
                                                                ! WS cell
    real (kind=dp), dimension (irid, nfund, ncelld), intent (in), &
      optional :: thetas           ! < shape function THETA=0 outer space
                                   ! THETA =1 inside WS cell in spherical
                                   ! harmonics expansion
    real (kind=dp), dimension (ntotd*(ncheb+1), nfund, ncelld), &
      intent (inout), optional :: thetasnew

    ! .. Local variables
    integer :: npan_inst, i_stat, i_all
    integer :: i1, ipot, ipotm, ir2, ip, icell
    integer :: imin, imax, iminnew, imaxnew, lm1
    integer :: ishift, ilogpanshift, ilinpanshift, npan_logtemp
    real (kind=dp) :: fac
    real (kind=dp) :: rmin, rmax, rval
    ! .. Allocatable arrays
    real (kind=dp), dimension (:, :, :), allocatable :: thetasin
    ! .. Assignment of values to parameters
    parameter (fac=2e0_dp)


    ! checks for optional arguments
    if (present(ntcell) .and. .not. present(thetas) .or. &
      .not. present(thetasnew)) then
      write (*, *) 'Error in create_newmesh:'
      write (*, *) 'List of optional arguments not complete'
      stop
    end if

    ! allocations
    if (present(ntcell)) then
      allocate (thetasin(irid,nfund,ncelld), stat=i_stat)
      call memocc(i_stat, product(shape(thetasin))*kind(thetasin), 'THETASIN', &
        'CREATE_NEWMESH')
      thetasin = 0.0e0_dp
    end if
    if (present(ntcell)) thetasnew = 0.0e0_dp

    ipotm = 0

    do i1 = 1, natyp

      ipot = nspin*(i1-1) + 1
      npan_inst = ipan(i1) - 1
      npan_tot(i1) = npan_log + npan_eq + npan_inst

      ! log panel
      rmin = r(2, i1)
      rmax = r_log
      rval = 0e0_dp
      ishift = 0
      if (r_log>r(irmin(i1),i1)) then
        ilogpanshift = 1
        ilinpanshift = 0
      else
        ilogpanshift = 0
        ilinpanshift = 1
      end if

      if (ilinpanshift==1) then
        stop 'non-spherical part of the potential needs to be inside &
          &the log panel'
      end if

      do ip = 0, npan_log - ilogpanshift
        rval = (fac**ip-1e0_dp)/(fac**(npan_log-ilogpanshift)-1e0_dp)
        rpan_intervall(ip+ishift, i1) = rmin + rval*(rmax-rmin)
        ipan_intervall(ip+ishift, i1) = (ip+ishift)*(ncheb+1)
        if (ishift==0 .and. rpan_intervall(ip,i1)>r(irmin(i1),i1)) then
          ishift = 1
          npan_logtemp = ip
          rpan_intervall(ip+1, i1) = rpan_intervall(ip, i1)
          ipan_intervall(ip+1, i1) = (ip+ishift)*(ncheb+1)
          rpan_intervall(ip, i1) = r(irmin(i1), i1)
          ipan_intervall(ip, i1) = ip*(ncheb+1)
        end if
      end do                       ! NPAN_LOG

      ! equivalent panel
      ishift = 0
      rmin = r_log
      rmax = r(ircut(1,i1), i1)
      do ip = 0, npan_eq - ilinpanshift
        rpan_intervall(ip+ishift+npan_log, i1) = rmin + &
          ip*(rmax-rmin)/(npan_eq-ilinpanshift)
        ipan_intervall(ip+ishift+npan_log, i1) = (npan_log+ip+ishift)* &
          (ncheb+1)
      end do                       ! NPAN_EQ

      ! intersection zone
      do ip = 1, npan_inst
        rpan_intervall(npan_log+npan_eq+ip, i1) = r(ircut(ip+1,i1), i1)
        ipan_intervall(npan_log+npan_eq+ip, i1) = (npan_log+npan_eq+ip)* &
          (ncheb+1)
      end do                       ! NPAN_INST

      npan_eq_at(i1) = npan_eq + npan_log - npan_logtemp
      npan_log_at(i1) = npan_logtemp

      call chebmesh(npan_tot(i1), ncheb, rpan_intervall(0:,i1), rnew(1,i1))

      ! do interpolation only when optional arguments are given
      if (present(ntcell)) then
        ! interpolate shape function THETAS to new shape function THETASNEW
        ! save THETAS to THETASIN
        icell = ntcell(i1)
        do lm1 = 1, nfund
          thetasin(:, lm1, icell) = thetas(:, lm1, icell)
          ir2 = 0
          do ip = npan_log_at(i1) + npan_eq_at(i1) + 1, npan_tot(i1)
            ir2 = ir2 + 1
            imin = ircut(ir2, i1) + 1
            imax = ircut(ir2+1, i1)
            iminnew = ipan_intervall(ip-1, i1) + 1
            imaxnew = ipan_intervall(ip, i1)
            call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
              thetasin(imin-ircut(1,i1):imax-ircut(1, &
              i1),lm1,icell), thetasnew(iminnew:imaxnew,lm1,icell), &
              imax-imin+1, imaxnew-iminnew+1)
          end do
        end do
      end if                       ! present(ntcell)
    end do                         ! i1

    if (present(ntcell)) then
      i_all = -product(shape(thetasin))*kind(thetasin)
      deallocate (thetasin, stat=i_stat)
      call memocc(i_stat, i_all, 'THETASIN', 'CREATE_NEWMESH')
    end if

  end subroutine create_newmesh

  ! ----------------------------------------------------------------------------
  ! SUBROUTINE: CHEBMESH
  ! ----------------------------------------------------------------------------
  subroutine chebmesh(npan, ncheb, ri, ro)

    implicit none

    integer, intent (in) :: npan
    integer, intent (in) :: ncheb  ! < Number of Chebychev pannels for the new
                                   ! solver
    real (kind=dp), dimension (0:npan), intent (in) :: ri
    real (kind=dp), dimension (npan*(ncheb+1)), intent (inout) :: ro
    integer :: i, k, ik
    real (kind=dp) :: tau

    do i = 1, npan
      do k = 0, ncheb
        ik = i*ncheb + i - k
        tau = cos(((2*k+1)*pi)/(2*(ncheb+1)))
        tau = 0.5e0_dp*((ri(i)-ri(i-1))*tau+ri(i)+ri(i-1))
        ro(ik) = tau
      end do
    end do
  end subroutine chebmesh

end module mod_create_newmesh
