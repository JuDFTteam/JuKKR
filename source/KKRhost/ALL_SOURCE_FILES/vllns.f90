!-------------------------------------------------------------------------------
! SUBROUTINE: VLLNS
!> @brief Transformation of the wavefunctions for non spherical potentials.
!> @details To determine the non - spherical wavefunctions the potential
!> has to be lm1 and lm2 dependent . the potential is stored
!> only as lm dependent , therefore a transformation in the
!> following way has to be done :
!> \f$ vnsll(r,lm1,lm2)   =  \sum_{lm3} \left\{  c(lm1,lm2,lm3) *vins(r,lm3)  \right\}\f$
!> where c(lm1,lm2,lm3) are the gaunt coeffients. (see notes by B. Drittler)
!> @author B. Drittler
!> @date July 1988
!> @note attention : The gaunt coeffients are stored in an index array only for lm1.gt.lm2
!> (see subroutine gaunt)
!> - R. Zeller Sep. 2000: modified
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine vllns(vnspll, vins, cleb, icleb, iend, irm, ncleb, lmpot, &
      irmind, lmmaxd)
      Use mod_datatypes, Only: dp

      Implicit None

! .. Input variables
      Integer, Intent (In) :: irm !< Maximum number of radial points
      Integer, Intent (In) :: iend
      Integer, Intent (In) :: ncleb !< Number of Clebsch-Gordon coefficients
      Integer, Intent (In) :: lmpot !< (LPOT+1)**2
      Integer, Intent (In) :: irmind !< IRM-IRNSD
      Integer, Intent (In) :: lmmaxd !< (KREL+KORBIT+1)(LMAX+1)^2
! .. Array Arguments
      Integer, Dimension (ncleb, 4), Intent (In) :: icleb !< Pointer array
      Real (Kind=dp), Dimension (ncleb, 2), Intent (In) :: cleb !< GAUNT coefficients (GAUNT)
      Real (Kind=dp), Dimension (irmind:irm, lmpot), Intent (In) :: vins !< Non-spherical part of the potential
! .. Output variables
      Real (Kind=dp), Dimension (lmmaxd, lmmaxd, irmind:irm), &
        Intent (Out) :: vnspll
! .. Local Scalars
      Integer :: ir, j, lm1, lm2, lm3
! ..
      Do lm1 = 1, lmmaxd
        Do lm2 = 1, lm1
          Do ir = irmind, irm
            vnspll(lm1, lm2, ir) = 0.0E0_dp
          End Do
        End Do
      End Do
!
      Do j = 1, iend
        lm1 = icleb(j, 1)
        lm2 = icleb(j, 2)
        lm3 = icleb(j, 3)
        Do ir = irmind, irm
          vnspll(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + &
            cleb(j, 1)*vins(ir, lm3)
        End Do
      End Do
!----------------------------------------------------------------------------
! Use symmetry of the gaunt coef.
!----------------------------------------------------------------------------
      Do lm1 = 1, lmmaxd
        Do lm2 = 1, lm1 - 1
          Do ir = irmind, irm
            vnspll(lm2, lm1, ir) = vnspll(lm1, lm2, ir)
          End Do
        End Do
      End Do

      Do lm1 = 1, lmmaxd
        Do ir = irmind, irm
          vnspll(lm1, lm1, ir) = vnspll(lm1, lm1, ir) + vins(ir, 1)
        End Do
      End Do

    End Subroutine
