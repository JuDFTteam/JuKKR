!-------------------------------------------------------------------------------
! SUBROUTINE: RHOVAL0
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine rhoval0(ez, drdi, rmesh, ipan, ircut, irws, thetas, dos0, dos1, &
      irm, lmax)
!
      Use constants
      Use global_variables
      Use mod_datatypes, Only: dp

      Implicit None

! .. Input variables
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
      Integer, Intent (In) :: ipan !< Number of panels in non-MT-region
      Integer, Intent (In) :: irws !< R point at WS radius
      Complex (Kind=dp), Intent (In) :: ez
      Integer, Dimension (0:ipand), Intent (In) :: ircut !< R points of panel borders
      Real (Kind=dp), Dimension (irm), Intent (In) :: drdi !< Derivative dr/di
      Real (Kind=dp), Dimension (irm), Intent (In) :: rmesh
      Real (Kind=dp), Dimension (irid, nfund), Intent (In) :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
! .. Output variables
      Complex (Kind=dp), Intent (Out) :: dos0
      Complex (Kind=dp), Intent (Out) :: dos1
! .. Local Scalars
      Integer :: ir, l, l1, imt1
      Integer :: lmaxd1
      Real (Kind=dp) :: c0ll
      Complex (Kind=dp) :: ek, ciek, denl
! .. Local Arrays ..
      Complex (Kind=dp), Dimension (0:lmax+1) :: bessjw
      Complex (Kind=dp), Dimension (0:lmax+1) :: bessyw
      Complex (Kind=dp), Dimension (0:lmax+1) :: hankws
      Complex (Kind=dp), Dimension (irm, 0:lmax) :: pz
      Complex (Kind=dp), Dimension (irm, 0:lmax) :: qz
      Complex (Kind=dp), Dimension (irm, 0:lmax+1) :: cden0
      Complex (Kind=dp), Dimension (irm, 0:lmax+1) :: cden1

! .. Intrinsic Functions
      Intrinsic :: atan, sqrt
!
      lmaxd1 = lmax + 1
      ek = sqrt(ez)
      c0ll = 1.0E0_dp/sqrt(16.0E0_dp*atan(1.0E0_dp))
      ciek = ci*ek
!
!----------------------------------------------------------------------------
      Do ir = 2, irws
        Call beshan(hankws, bessjw, bessyw, rmesh(ir)*ek, lmaxd1)
        Do l = 0, lmax
          pz(ir, l) = bessjw(l)*rmesh(ir)
          qz(ir, l) = (bessyw(l)-ci*bessjw(l))*rmesh(ir)
        End Do
      End Do
      imt1 = ircut(1)
      Do l1 = 0, lmaxd1
        cden0(1, l1) = czero
        cden1(1, l1) = czero
      End Do
      Do ir = 2, irws
        cden0(ir, 0) = ek*pz(ir, 0)*qz(ir, 0)
        cden1(ir, 0) = ek*pz(ir, 0)**2*(0.E0_dp, -1.E0_dp)
        cden1(ir, lmaxd1) = ciek*rmesh(ir)**2
      End Do
      Do l1 = 1, lmax
        Do ir = 2, irws
          cden0(ir, l1) = ek*pz(ir, l1)*qz(ir, l1)*(l1+l1+1)
          cden1(ir, l1) = ek*pz(ir, l1)**2*(0.E0_dp, -1.E0_dp)*(l1+l1+1)
        End Do
      End Do
!
      Do l1 = 0, lmaxd1 !LMAXD1
        If (ipan>1) Then
          Do ir = imt1 + 1, irws
            cden0(ir, l1) = cden0(ir, l1)*thetas(ir-imt1, 1)*c0ll
            cden1(ir, l1) = cden1(ir, l1)*thetas(ir-imt1, 1)*c0ll
          End Do
        End If
        Call csimpk(cden0(1,l1), denl, ipan, ircut, drdi)
        dos0 = dos0 + denl
        Call csimpk(cden1(1,l1), denl, ipan, ircut, drdi)
        dos1 = dos1 + denl
      End Do
!
    End Subroutine
