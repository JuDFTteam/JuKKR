!     >> Input parameters
!>    @param     ALAT    lattice constant a
!>    @param     NACLS   number of atoms in cluster
!>    @param     RR      array of real space vectors
!>    @param     EZOA
!>    @param     BZKP
!>    @param     nrd     There are nrd+1 real space vectors in RR
!>    @param     naclsd. maximal number of atoms in cluster

!     << Output parameters
!>    @param     EIKRM   Fourier exponential factor with minus sign
!>    @param     EIKRP   Fourier exponential factor with plus sign

      subroutine dlke1(alat, nacls, rr, ezoa, bzkp, eikrm, eikrp, nrd, naclsd)
      implicit none
! ----------------------------------------------------------------------
!     Fourier transformation of the cluster Greens function
!     Prepares the calculation (calculates Fourier factors) for dlke0
! ----------------------------------------------------------------------
      integer, intent(in) :: nrd, naclsd ! todo: remove these and change arrays to deferred shape
      double precision, intent(in) :: alat
      integer, intent(in) :: nacls !< number of vectors in the cluster
      integer, intent(in) :: ezoa(*) !< index list of ...
      double complex, intent(out) :: eikrp(naclsd), eikrm(naclsd)
      double precision, intent(in) :: bzkp(1:3) !< k-point (vector in the Brillouin zone)
      double precision, intent(in) :: rr(3,0:nrd) !< real space cluster vectors
      
      double complex, parameter :: ci=(0.d0,1.d0)
      double precision :: convpuh, tpi
      integer :: m
      double complex :: tt, exparg

      tpi = 8.d0*atan(1.d0)         
      convpuh = alat/tpi * 0.5d0

      do m = 1, nacls
!     
!     Here we do   --                  nn'
!                  \                   ii'          ii'
!                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
!                  --          n'  n   LL'          LL'
!                  n'
!  Be careful a minus sign must be included here. RR is not
!  symmetric around each atom. The minus comes from the fact that
!  the repulsive potential GF is calculated for 0n and not n0!                   
!  and that is why we need a minus sign extra!
!  
        tt = -ci*tpi*dot_product(bzkp(1:3), rr(1:3,ezoa(m))) ! purely imaginary number
!
!  convert to p.u. and multiply with 1/2 (done above)
        exparg = exp(tt)
        eikrp(m) =       exparg  * convpuh
        eikrm(m) = conjg(exparg) * convpuh ! we can re-use exparg here instead of exp(-tt) since tt is purely imaginary
      enddo ! m
      
      endsubroutine ! dlke1
