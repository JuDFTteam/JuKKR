!-------------------------------------------------------------------------------
! MODULE: mod_create_newmesh
!> @brief Module for the creation of the integration grid for the new solver
!-------------------------------------------------------------------------------
    Module mod_create_newmesh

      Use constants
      Use profiling
      Use mod_datatypes, Only: dp

      Implicit None
      Private :: dp

    Contains

!----------------------------------------------------------------------------
! SUBROUTINE: CREATE_NEWMESH
!> @brief Creation of the integration grid for the new solver
!> @note changed interface to get rid of inc.p and to be able to use i
!> create_newmesh in tmatimp routine for GREENIMP option
!> this is the list of  array dimensions previously importted from inc.p
!----------------------------------------------------------------------------
      Subroutine create_newmesh(natyp, lmax, lpot, irm, irnsd, ipand, irid, &
        ntotd, nfund, ncheb, irmdnew, nspin, r, irmin, ipan, ircut, r_log, &
        npan_log, npan_eq, npan_log_at, npan_eq_at, npan_tot, rnew, &
        rpan_intervall, ipan_intervall, ncelld, ntcell, thetas, thetasnew) !< optional arguments

        Implicit None

! .. Input variables
        Integer, Intent (In) :: irm !< Maximum number of radial points
        Integer, Intent (In) :: irid !< Shape functions parameters in non-spherical part
        Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
        Integer, Intent (In) :: lpot !< Maximum l component in potential expansion
        Integer, Intent (In) :: nspin !< Counter for spin directions
        Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
        Integer, Intent (In) :: irnsd
        Integer, Intent (In) :: ipand !< Number of panels in non-spherical part
        Integer, Intent (In) :: ntotd
        Integer, Intent (In) :: nfund !< Shape functions parameters in non-spherical part
        Integer, Intent (In) :: ncheb !< Number of Chebychev pannels for the new solver
        Integer, Intent (In) :: ncelld !< Number of cells (shapes) in non-spherical part
        Integer, Intent (In) :: irmdnew
        Integer, Intent (In) :: npan_eq !< Variables for the pannels for the new solver
        Integer, Intent (In) :: npan_log !< Variables for the pannels for the new solver
        Real (Kind=dp), Intent (In) :: r_log
        Integer, Dimension (natyp), Intent (In) :: ipan !< Number of panels in non-MT-region
        Integer, Dimension (natyp), Intent (In) :: irmin !< Max R for spherical treatment
        Integer, Dimension (0:ipand, natyp), Intent (In) :: ircut !< R points of panel borders
        Real (Kind=dp), Dimension (irm, natyp), Intent (In) :: r !< Radial mesh ( in units a Bohr)
! .. Input/Output variables
        Integer, Dimension (natyp), Intent (Inout) :: npan_tot
        Integer, Dimension (natyp), Intent (Inout) :: npan_eq_at
        Integer, Dimension (natyp), Intent (Inout) :: npan_log_at
        Integer, Dimension (0:ntotd, natyp), Intent (Inout) :: ipan_intervall
        Real (Kind=dp), Dimension (irmdnew, natyp), Intent (Inout) :: rnew
        Real (Kind=dp), Dimension (0:ntotd, natyp), Intent (Inout) :: &
          rpan_intervall
! .. Optional arguments, do interpolation when given
        Integer, Dimension (natyp), Intent (In), Optional :: ntcell !< Index for WS cell
        Real (Kind=dp), Dimension (irid, nfund, ncelld), Intent (In), &
          Optional :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
        Real (Kind=dp), Dimension (ntotd*(ncheb+1), nfund, ncelld), &
          Intent (Inout), Optional :: thetasnew

! .. Local variables
        Integer :: npan_inst, i_stat, i_all
        Integer :: i1, ipot, ipotm, ir2, ip, icell
        Integer :: imin, imax, iminnew, imaxnew, lm1
        Integer :: ishift, ilogpanshift, ilinpanshift, npan_logtemp
        Real (Kind=dp) :: fac
        Real (Kind=dp) :: rmin, rmax, rval
! .. Allocatable arrays
        Real (Kind=dp), Dimension (:, :, :), Allocatable :: thetasin
! .. Assignment of values to parameters
        Parameter (fac=2E0_dp)


! checks for optional arguments
        If (present(ntcell) .And. .Not. present(thetas) .Or. &
          .Not. present(thetasnew)) Then
          Write (*, *) 'Error in create_newmesh:'
          Write (*, *) 'List of optional arguments not complete'
          Stop
        End If

! allocations
        If (present(ntcell)) Then
          Allocate (thetasin(irid,nfund,ncelld), Stat=i_stat)
          Call memocc(i_stat, product(shape(thetasin))*kind(thetasin), &
            'THETASIN', 'CREATE_NEWMESH')
          thetasin = 0.0E0_dp
        End If
        If (present(ntcell)) thetasnew = 0.0E0_dp

        ipotm = 0

        Do i1 = 1, natyp

          ipot = nspin*(i1-1) + 1
          npan_inst = ipan(i1) - 1
          npan_tot(i1) = npan_log + npan_eq + npan_inst

! log panel
          rmin = r(2, i1)
          rmax = r_log
          rval = 0E0_dp
          ishift = 0
          If (r_log>r(irmin(i1),i1)) Then
            ilogpanshift = 1
            ilinpanshift = 0
          Else
            ilogpanshift = 0
            ilinpanshift = 1
          End If

          If (ilinpanshift==1) Then
            Stop 'non-spherical part of the potential needs to &
              &be inside the log panel'
          End If

          Do ip = 0, npan_log - ilogpanshift
            rval = (fac**ip-1E0_dp)/(fac**(npan_log-ilogpanshift)-1E0_dp)
            rpan_intervall(ip+ishift, i1) = rmin + rval*(rmax-rmin)
            ipan_intervall(ip+ishift, i1) = (ip+ishift)*(ncheb+1)
            If (ishift==0 .And. rpan_intervall(ip,i1)>r(irmin(i1),i1)) Then
              ishift = 1
              npan_logtemp = ip
              rpan_intervall(ip+1, i1) = rpan_intervall(ip, i1)
              ipan_intervall(ip+1, i1) = (ip+ishift)*(ncheb+1)
              rpan_intervall(ip, i1) = r(irmin(i1), i1)
              ipan_intervall(ip, i1) = ip*(ncheb+1)
            End If
          End Do ! NPAN_LOG

! equivalent panel
          ishift = 0
          rmin = r_log
          rmax = r(ircut(1,i1), i1)
          Do ip = 0, npan_eq - ilinpanshift
            rpan_intervall(ip+ishift+npan_log, i1) = rmin + &
              ip*(rmax-rmin)/(npan_eq-ilinpanshift)
            ipan_intervall(ip+ishift+npan_log, i1) = (npan_log+ip+ishift)* &
              (ncheb+1)
          End Do ! NPAN_EQ

! intersection zone
          Do ip = 1, npan_inst
            rpan_intervall(npan_log+npan_eq+ip, i1) = r(ircut(ip+1,i1), i1)
            ipan_intervall(npan_log+npan_eq+ip, i1) = (npan_log+npan_eq+ip)* &
              (ncheb+1)
          End Do ! NPAN_INST

          npan_eq_at(i1) = npan_eq + npan_log - npan_logtemp
          npan_log_at(i1) = npan_logtemp

          Call chebmesh(npan_tot(i1), ncheb, rpan_intervall(0:,i1), &
            rnew(1,i1))

! do interpolation only when optional arguments are given
          If (present(ntcell)) Then
! interpolate shape function THETAS to new shape function THETASNEW
! save THETAS to THETASIN
            icell = ntcell(i1)
            Do lm1 = 1, nfund
              thetasin(:, lm1, icell) = thetas(:, lm1, icell)
              ir2 = 0
              Do ip = npan_log_at(i1) + npan_eq_at(i1) + 1, npan_tot(i1)
                ir2 = ir2 + 1
                imin = ircut(ir2, i1) + 1
                imax = ircut(ir2+1, i1)
                iminnew = ipan_intervall(ip-1, i1) + 1
                imaxnew = ipan_intervall(ip, i1)
                Call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
                  thetasin(imin-ircut(1,i1):imax-ircut(1, &
                  i1),lm1,icell), thetasnew(iminnew:imaxnew,lm1,icell), &
                  imax-imin+1, imaxnew-iminnew+1)
              End Do
            End Do
          End If ! present(ntcell)
        End Do ! i1

        If (present(ntcell)) Then
          i_all = -product(shape(thetasin))*kind(thetasin)
          Deallocate (thetasin, Stat=i_stat)
          Call memocc(i_stat, i_all, 'THETASIN', 'CREATE_NEWMESH')
        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: CHEBMESH
!----------------------------------------------------------------------------
      Subroutine chebmesh(npan, ncheb, ri, ro)

        Implicit None

        Integer, Intent (In) :: npan
        Integer, Intent (In) :: ncheb !< Number of Chebychev pannels for the new solver
        Real (Kind=dp), Dimension (0:npan), Intent (In) :: ri
        Real (Kind=dp), Dimension (npan*(ncheb+1)), Intent (Inout) :: ro
        Integer :: i, k, ik
        Real (Kind=dp) :: tau

        Do i = 1, npan
          Do k = 0, ncheb
            ik = i*ncheb + i - k
            tau = cos(((2*k+1)*pi)/(2*(ncheb+1)))
            tau = 0.5E0_dp*((ri(i)-ri(i-1))*tau+ri(i)+ri(i-1))
            ro(ik) = tau
          End Do
        End Do
      End Subroutine

    End Module
