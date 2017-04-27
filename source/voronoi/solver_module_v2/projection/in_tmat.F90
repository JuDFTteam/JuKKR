  subroutine in_tmat(tmat)
! Get t-matrix from outsusc.dat
! Reads the whole thing for a given spin channel and energy point
  use global

  implicit none

! --> the t-matrix for a given spin channel
  complex(kind=c8b), intent(out) :: tmat(lmmax,lmmax,nasusc)
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ja, lmi, lmj, i3(3)
  real(kind=r8b)    :: x(1000)
  character*60      :: header

!  write(iodb,*) "lmmax=", lmmax
  read(iomain,'(a)') header
!  write(iodb,*) header
  do ia=1,nasusc
    read(iomain,*) i3
!    write(iodb,'(3i4)') i3
    do lmj=1,lmmax
      read(iomain,*) x(1:2*lmmax)
      do lmi=1,lmmax
        tmat(lmi,lmj,ia) = cmplx(x(2*lmi-1),x(2*lmi))
!        write(iodb,'(2i4,6es12.4)') lmi, lmj, tmat(lmi,lmj,ia)
      end do
    end do
  end do
! All done
  end subroutine in_tmat
