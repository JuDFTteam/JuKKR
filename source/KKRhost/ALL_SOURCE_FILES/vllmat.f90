!-------------------------------------------------------------------------------
! SUBROUTINE: VLLMAT
!> @brief
!-------------------------------------------------------------------------------
    Subroutine vllmat(irmin, nrmaxd, irc, lmmax, lmmaxso, vnspll0, vins, &
      lmpot, cleb, icleb, iend, nspin, z, rnew, use_sratrick, ncleb)
      Use mod_datatypes, Only: dp

      Implicit None

!.. Input variables
      Integer, Intent (In) :: irc !< r point for potential cutting
      Integer, Intent (In) :: iend
      Integer, Intent (In) :: ncleb !< Number of Clebsch-Gordon coefficients
      Integer, Intent (In) :: irmin !< max r for spherical treatment
      Integer, Intent (In) :: lmmax !< (LMAX+1)^2
      Integer, Intent (In) :: nspin !< Counter for spin directions
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: nrmaxd !< NTOTD*(NCHEBD+1)
      Integer, Intent (In) :: lmmaxso
      Integer, Intent (In) :: use_sratrick
      Real (Kind=dp), Intent (In) :: z

      Integer, Dimension (ncleb, 4), Intent (In) :: icleb
      Real (Kind=dp), Dimension (*), Intent (In) :: cleb !< GAUNT coefficients (GAUNT)
      Real (Kind=dp), Dimension (irmin:irc, lmpot, nspin), Intent (In) :: vins !< Non-spherical part of the potential
      Real (Kind=dp), Dimension (irmin:nrmaxd), Intent (In) :: rnew
      Complex (Kind=dp), Dimension (lmmaxso, lmmaxso, irmin:irc), &
        Intent (Out) :: vnspll0

!.. Local variables
      Integer :: isp
      Integer :: i, ir, j, lm1, lm2, lm3
      Real (Kind=dp), Dimension (lmmax, lmmax, irmin:irc, nspin) :: vnspll

      Do isp = 1, nspin
        Do lm1 = 1, lmmax
          Do lm2 = 1, lm1
            Do ir = irmin, irc
              vnspll(lm1, lm2, ir, isp) = 0.0E0_dp
            End Do ! IR
          End Do ! LM2
        End Do ! LM11

        Do j = 1, iend
          lm1 = icleb(j, 1)
          lm2 = icleb(j, 2)
          lm3 = icleb(j, 3)
          Do i = irmin, irc
            vnspll(lm1, lm2, i, isp) = vnspll(lm1, lm2, i, isp) + &
              cleb(j)*vins(i, lm3, isp)
          End Do ! I
        End Do ! J
!-------------------------------------------------------------------------
! Use symmetry of the gaunt coef.
!-------------------------------------------------------------------------
        Do lm1 = 1, lmmax
          Do lm2 = 1, lm1 - 1
            Do i = irmin, irc
              vnspll(lm2, lm1, i, isp) = vnspll(lm1, lm2, i, isp)
            End Do ! I
          End Do ! LM2
        End Do ! LM1

        If (use_sratrick==0) Then
          Do lm1 = 1, lmmax
            Do i = irmin, irc
              vnspll(lm1, lm1, i, isp) = vnspll(lm1, lm1, i, isp) + &
                vins(i, 1, isp) - 2E0_dp*z/rnew(i)
            End Do
          End Do
        End If

      End Do !NSPIN

! Set vnspll as twice as large
      vnspll0(1:lmmax, 1:lmmax, irmin:irc) = cmplx(vnspll(1:lmmax,1:lmmax, &
        irmin:irc,1), 0E0_dp, kind=dp)

      If (nspin==2) Then ! hack to make routine work for Bxc-field
        vnspll0(lmmax+1:lmmaxso, lmmax+1:lmmaxso, irmin:irc) &
          = cmplx(vnspll(1:lmmax,1:lmmax,irmin:irc,nspin), 0E0_dp, kind=dp)
      End If

    End Subroutine
