  subroutine free_gf(sqrten,rij,g0)
! Free space GF

  use global, only: i4b, r8b, c8b, rgaunt, lmmax, nlmax2, lmmax2, i2lm
  use bessel_new

  implicit none

! Energy
  complex(kind=c8b), intent(in)  :: sqrten
! Vector connecting i to j
  real(kind=r8b),    intent(in)  :: rij(3)
! Free space GF
  complex(kind=c8b), intent(out) :: g0(lmmax,lmmax)
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: iu = (0.d0,1.d0)
  real(kind=r8b),    parameter :: pi = 4.d0*atan(1.d0)
  real(kind=r8b)    :: rijlen, rylm(lmmax2)
  integer(kind=i4b) :: il, ilm, jl, jlm, kl, klm
  complex(kind=c8b) :: h0, gf, fac

  g0(:,:) = 0.d0
  rijlen = sqrt(dot_product(rij,rij))
  if (rijlen < 1.d-8) return
  fac = -iu*sqrten  !/(4.d0*pi)  ! due to normalization of Ylm
! spherical harmonics
  call rymy(rij,nlmax2,lmmax2,rylm)
! put together the free space GF
  do jlm=5,9!1,lmmax
    jl = i2lm(2,jlm)
    do ilm=5,9!1,lmmax
      il = i2lm(2,ilm)
!     internal angular momentum sum
      gf = 0.d0
      do klm=1,lmmax2
!       spherical Hankel function of the first kind
        kl = i2lm(2,klm)
!        if (kl < 4 .or. kl > 4) cycle
        h0 = bessh1(kl,sqrten*rijlen)
        gf = gf + h0*rylm(klm)*rgaunt(ilm,jlm,klm)
      end do
!     put the result in the right place
      g0(ilm,jlm) = gf*fac*(-1)**(il-jl)
      if (il /= jl) g0(ilm,jlm) = 0.d0
    end do
  end do
! All done!
  end subroutine free_gf
