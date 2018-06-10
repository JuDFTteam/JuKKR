!-------------------------------------------------------------------------------
! SUBROUTINE: INTERPOLATE_POTEN
!> @brief Routine for the interpolation of the potential in the integration grid.
!> @note Jonathan Chico: Include array dimensions in interface explicitly to get rid of
!> inc.p import and to be able to use routine for different number of atoms
!-------------------------------------------------------------------------------
subroutine interpolate_poten(lpot, irm, irnsd, natyp, ipand, lmpot, nspotd, &
  ntotd, ncheb, irmdnew, nspin, r, irmind, irmin, irws, ircut, vins, visp, &
  npan_log, npan_eq, npan_tot, rnew, ipan_intervall, vinsnew)

  implicit none

  integer, intent (in) :: irm !< Maximum number of radial points
  integer, intent (in) :: lpot !< Maximum l component in potential expansion
  integer, intent (in) :: irnsd
  integer, intent (in) :: ntotd
  integer, intent (in) :: ncheb !< Number of Chebychev pannels for the new solver
  integer, intent (in) :: natyp !< Number of kinds of atoms in unit cell
  integer, intent (in) :: ipand !< Number of panels in non-spherical part
  integer, intent (in) :: lmpot !< (LPOT+1)**2
  integer, intent (in) :: nspin !< Counter for spin directions
  integer, intent (in) :: nspotd
  integer, intent (in) :: irmind !< IRM-IRNSD
  integer, intent (in) :: irmdnew
  integer, dimension (natyp), intent (in) :: irws !< R point at WS radius
  integer, dimension (natyp), intent (in) :: irmin !< Max R for spherical treatment
  integer, dimension (natyp), intent (in) :: npan_eq !< Variables for the pannels for the new solver
  integer, dimension (natyp), intent (in) :: npan_log !< Variables for the pannels for the new solver
  integer, dimension (natyp), intent (in) :: npan_tot
  integer, dimension (0:ipand, natyp), intent (in) :: ircut !< R points of panel borders
  integer, dimension (0:ntotd, natyp), intent (in) :: ipan_intervall
  double precision, dimension (irm, natyp), intent (in) :: r
  double precision, dimension (irm, nspotd), intent (in) :: visp
  double precision, dimension ((irm-irnsd):irm, (lpot+1)**2, nspotd), &
    intent (in) :: vins
! .. Input/Output variables
  double precision, dimension (irmdnew, (lpot+1)**2, nspotd), &
    intent (inout) :: vinsnew

! .. Local variables
  integer :: i1, ipot, ipotm, imin, imax, ip, ir, lm1, ispin, iminnew, &
    imaxnew, ir2
  double precision, dimension (irmdnew, natyp) :: rnew
  double precision, dimension (irm, (lpot+1)**2, nspin) :: vinsin

  vinsnew = 0d0
  ipotm = 0

! interpolate potential to new mesh
  do i1 = 1, natyp

    ipot = nspin*(i1-1) + 1

! save input potential to VINSIN
    vinsin = 0d0
    do ir = 1, irws(i1)
      if (ir<irmin(i1)) then
        vinsin(ir, 1, 1) = visp(ir, ipot)
        vinsin(ir, 1, nspin) = visp(ir, ipot+nspin-1)
      else
        do lm1 = 1, lmpot
          if (lm1==1) then
            vinsin(ir, lm1, 1) = visp(ir, ipot)
            vinsin(ir, lm1, nspin) = visp(ir, ipot+nspin-1)
          else
            vinsin(ir, lm1, 1) = vins(ir, lm1, ipot)
            vinsin(ir, lm1, nspin) = vins(ir, lm1, ipot+nspin-1)
          end if
        end do
      end if
    end do

    do ispin = 1, nspin
      ipotm = ipotm + 1
      do lm1 = 1, lmpot

        imin = 1
        imax = irmin(i1)
        do ip = 1, npan_log(i1)
          iminnew = ipan_intervall(ip-1, i1) + 1
          imaxnew = ipan_intervall(ip, i1)
          call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm), &
            imax-imin+1, imaxnew-iminnew+1)
        end do

        imin = irmin(i1)
        imax = ircut(1, i1)
        do ip = npan_log(i1) + 1, npan_log(i1) + npan_eq(i1)
          iminnew = ipan_intervall(ip-1, i1) + 1
          imaxnew = ipan_intervall(ip, i1)
          call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm), &
            imax-imin+1, imaxnew-iminnew+1)
        end do

        ir2 = 0
        do ip = npan_log(i1) + npan_eq(i1) + 1, npan_tot(i1)
          ir2 = ir2 + 1
          imin = ircut(ir2, i1) + 1
          imax = ircut(ir2+1, i1)
          iminnew = ipan_intervall(ip-1, i1) + 1
          imaxnew = ipan_intervall(ip, i1)
          call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm), &
            imax-imin+1, imaxnew-iminnew+1)
        end do
      end do ! lm1
    end do ! ispin
  end do ! i1
end subroutine
