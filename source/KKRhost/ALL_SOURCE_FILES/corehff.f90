subroutine corehff(kap1, kap2, mj, s, nsol, bhf, gck, fck, rc, drdic, rnuc, &
  nzero, nrc)
!   ********************************************************************
!   *                                                                  *
!   *   CALCULATE THE RELATIVISTIC HYPERFINEFIELDS FOR THE             *
!   *                  CURRENT  CORE STATE S                           *
!   *                                                                  *
!   *   THE WAVE FUNCTION  {G(K,S),F(K,S)}  IS NORMALIZED TO 1         *
!   *                                                                  *
!   ********************************************************************

  implicit none


! PARAMETER definitions
  real *8 :: e0, a0, cautog

!CONVERSION FACTOR FOR HYPERFINE FIELDS FROM A.U. TO GAUSS
!                                 ELECTRON CHARGE     IN ESU
!                                 BOHR-RADIUS         IN CM


  parameter (e0=1.6021892d-19*2.997930d+09, a0=0.52917706d-08, &
    cautog=e0/(a0*a0))

! Dummy arguments
  integer :: kap1, kap2, nrc, nsol, nzero, s
  real *8 :: mj, rnuc
  real *8 :: bhf(2, 2), drdic(nrc), fck(2, 2, nrc), gck(2, 2, nrc), rc(nrc)

! Local variables
  real *8 :: ame(2, 2), xx(5), yi(nrc), yy(5), zi(nrc)
  double precision :: dble, dsqrt
  integer :: i, k1, k2, n
  real *8 :: ylag

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ANGULAR HYPERFINE MATRIX ELEMENTS   SEE E.G.  E.M.ROSE
!        THE FACTOR  I  HAS BEEN OMITTED

  ame(1, 1) = 4.0d0*kap1*mj/(4.0d0*kap1*kap1-1.0d0)
  if (nsol==2) then
    ame(2, 2) = 4.0d0*kap2*mj/(4.0d0*kap2*kap2-1.0d0)
    ame(1, 2) = dsqrt(0.25d0-(mj/dble(kap1-kap2))**2)
    ame(2, 1) = ame(1, 2)
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  do k1 = 1, nsol
    do k2 = 1, nsol
      do n = 1, nzero
        yi(n) = drdic(n)*(gck(k1,s,n)*fck(k2,s,n)+fck(k1,s,n)*gck(k2,s,n))
      end do
      call rint4pts(yi, nzero, zi)
      if (rnuc/=0.0d0) then
        do i = 1, 5
          xx(i) = rc(nzero-5+i)
          yy(i) = zi(nzero-5+i)
        end do
        zi(nzero) = ylag(rnuc, xx, yy, 0, 4, 5)
      end if
      xx(1) = 1.0d0
!                      !RC( 1)
      xx(2) = 6.0d0
!                      !RC( 6)
      xx(3) = 11.0d0
!                      !RC(11)
      yy(1) = zi(nzero) - zi(1)
      yy(2) = zi(nzero) - zi(6)
      yy(3) = zi(nzero) - zi(11)
      bhf(k1, k2) = cautog*ame(k1, k2)*ylag(0.0d0, xx, yy, 0, 2, 3)
    end do
  end do

end subroutine
