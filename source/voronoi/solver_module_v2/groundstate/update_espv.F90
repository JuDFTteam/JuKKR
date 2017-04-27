  subroutine update_espv(lmx,natypd,nspind,espv)
! Save non-spherical charge and magnetization densities
  use global

  implicit none

! dimensions of espv array
  integer(kind=i4b), intent(in)    :: lmx, natypd, nspind
! Radial mesh, powers of r and radial integration weights
  real(kind=r8b),    intent(inout) :: espv(0:lmx-1,natypd,nspind)
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, il, is

  if (nlmax > lmx)    stop 'update_espv: nlmax > lmx'
  if (nsmax > nspind) stop 'update_espv: nsmax > nspind'
! **************
  do ia=1,nasusc
! **************
    ih = iasusc(ia)
    do is=1,nsmax
      do il=0,nlmax
        espv(il,ih,is) = ebandv(il,is,ia)
      end do
    end do
! ******
  end do
! ******
! All done!
  end subroutine update_espv

