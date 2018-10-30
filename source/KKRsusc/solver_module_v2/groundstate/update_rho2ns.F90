  subroutine update_rho2ns(lmaxpot,irmkd,lmpotd,natypd,rho2ns)
! Save non-spherical charge and magnetization densities
  use global

  implicit none

! lmax for potential; dimensions of rho2ns array
  integer(kind=i4b), intent(in)    :: lmaxpot, irmkd, lmpotd, natypd
! Radial mesh, powers of r and radial integration weights
  real(kind=r8b),    intent(inout) :: rho2ns(irmkd,lmpotd,natypd,nsmax)
! -----------------------------------------------------------------
  real(kind=r8b), parameter :: sqrt4pi = 4.d0*sqrt(atan(1.d0))
  real(kind=r8b), parameter :: tol = 1.d-8
  integer(kind=i4b) :: ia, ih, nr0, nr1, nr, ir, lmpot, ilm, lm
  integer(kind=i4b) :: nonzero(lmpotd)
  real(kind=r8b)    :: qlm(lmpotd), mlm(lmpotd), q, m
  complex(kind=c8b) :: work(nrmax), norm
  character*14      :: filename

  lmpot = (lmaxpot + 1)**2
  if (lmpot > lmpotd) stop 'update_rho2ns: lmpot > lmpotd'
  write(*,*) nasusc
!  write(*,'("update_rho2ns: killing spin density for ia > 3 !!!")')
! **************
  do ia=1,nasusc
! **************
    if (normesh(ia))  stop 'update_rho2ns: save rmesh first!'
    ih  = iasusc(ia)
    nr  = nrpts(ia)
    nr0 = irmkd - nr + 1 !nrpts0(ia)
    nr1 = irmkd!nrpts1(ia)
    do lm=1,lmpot
      rho2ns(nr0:nr1,lm,ih,1) = new_rho2ns(1:nr,lm,1,ia)/sqrt4pi
      rho2ns(nr0:nr1,lm,ih,2) = new_rho2ns(1:nr,lm,2,ia)/sqrt4pi
    end do
! ******
  end do
! ******
! All done!
  end subroutine update_rho2ns
