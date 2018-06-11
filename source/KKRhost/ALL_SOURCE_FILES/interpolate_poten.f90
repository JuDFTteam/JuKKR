!-------------------------------------------------------------------------------
! SUBROUTINE: INTERPOLATE_POTEN
!> @brief Routine for the interpolation of the potential in the integration grid.
!> @note Jonathan Chico: Include array dimensions in interface explicitly to get rid of
!> inc.p import and to be able to use routine for different number of atoms
!-------------------------------------------------------------------------------
    Subroutine interpolate_poten(lpot, irm, irnsd, natyp, ipand, lmpot, &
      nspotd, ntotd, ncheb, irmdnew, nspin, r, irmind, irmin, irws, ircut, &
      vins, visp, npan_log, npan_eq, npan_tot, rnew, ipan_intervall, vinsnew)
      Use mod_datatypes, Only: dp

      Implicit None

      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: lpot !< Maximum l component in potential expansion
      Integer, Intent (In) :: irnsd
      Integer, Intent (In) :: ntotd
      Integer, Intent (In) :: ncheb !< Number of Chebychev pannels for the new solver
      Integer, Intent (In) :: natyp !< Number of kinds of atoms in unit cell
      Integer, Intent (In) :: ipand !< Number of panels in non-spherical part
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: nspotd
      Integer, Intent (In) :: irmind !< IRM-IRNSD
      Integer, Intent (In) :: irmdnew
      Integer, Dimension (natyp), Intent (In) :: irws !< R point at WS radius
      Integer, Dimension (natyp), Intent (In) :: irmin !< Max R for spherical treatment
      Integer, Dimension (natyp), Intent (In) :: npan_eq !< Variables for the pannels for the new solver
      Integer, Dimension (natyp), Intent (In) :: npan_log !< Variables for the pannels for the new solver
      Integer, Dimension (natyp), Intent (In) :: npan_tot
      Integer, Dimension (0:ipand, natyp), Intent (In) :: ircut !< R points of panel borders
      Integer, Dimension (0:ntotd, natyp), Intent (In) :: ipan_intervall
      Real (Kind=dp), Dimension (irm, natyp), Intent (In) :: r
      Real (Kind=dp), Dimension (irm, nspotd), Intent (In) :: visp
      Real (Kind=dp), Dimension ((irm-irnsd):irm, (lpot+1)**2, nspotd), &
        Intent (In) :: vins
! .. Input/Output variables
      Real (Kind=dp), Dimension (irmdnew, (lpot+1)**2, nspotd), &
        Intent (Inout) :: vinsnew

! .. Local variables
      Integer :: i1, ipot, ipotm, imin, imax, ip, ir, lm1, ispin, iminnew, &
        imaxnew, ir2
      Real (Kind=dp), Dimension (irmdnew, natyp) :: rnew
      Real (Kind=dp), Dimension (irm, (lpot+1)**2, nspin) :: vinsin

      vinsnew = 0E0_dp
      ipotm = 0

! interpolate potential to new mesh
      Do i1 = 1, natyp

        ipot = nspin*(i1-1) + 1

! save input potential to VINSIN
        vinsin = 0E0_dp
        Do ir = 1, irws(i1)
          If (ir<irmin(i1)) Then
            vinsin(ir, 1, 1) = visp(ir, ipot)
            vinsin(ir, 1, nspin) = visp(ir, ipot+nspin-1)
          Else
            Do lm1 = 1, lmpot
              If (lm1==1) Then
                vinsin(ir, lm1, 1) = visp(ir, ipot)
                vinsin(ir, lm1, nspin) = visp(ir, ipot+nspin-1)
              Else
                vinsin(ir, lm1, 1) = vins(ir, lm1, ipot)
                vinsin(ir, lm1, nspin) = vins(ir, lm1, ipot+nspin-1)
              End If
            End Do
          End If
        End Do

        Do ispin = 1, nspin
          ipotm = ipotm + 1
          Do lm1 = 1, lmpot

            imin = 1
            imax = irmin(i1)
            Do ip = 1, npan_log(i1)
              iminnew = ipan_intervall(ip-1, i1) + 1
              imaxnew = ipan_intervall(ip, i1)
              Call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
                vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm &
                ), imax-imin+1, imaxnew-iminnew+1)
            End Do

            imin = irmin(i1)
            imax = ircut(1, i1)
            Do ip = npan_log(i1) + 1, npan_log(i1) + npan_eq(i1)
              iminnew = ipan_intervall(ip-1, i1) + 1
              imaxnew = ipan_intervall(ip, i1)
              Call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
                vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm &
                ), imax-imin+1, imaxnew-iminnew+1)
            End Do

            ir2 = 0
            Do ip = npan_log(i1) + npan_eq(i1) + 1, npan_tot(i1)
              ir2 = ir2 + 1
              imin = ircut(ir2, i1) + 1
              imax = ircut(ir2+1, i1)
              iminnew = ipan_intervall(ip-1, i1) + 1
              imaxnew = ipan_intervall(ip, i1)
              Call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
                vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm &
                ), imax-imin+1, imaxnew-iminnew+1)
            End Do
          End Do ! lm1
        End Do ! ispin
      End Do ! i1
    End Subroutine
