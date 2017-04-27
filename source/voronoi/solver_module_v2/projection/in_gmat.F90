  subroutine in_gmat(gmat)
! Get structural GF from outsusc.dat
! Reads the whole thing for a given spin channel and energy point
  use global

  implicit none

! --> the structural GF matrix for a given spin channel
  complex(kind=c8b), intent(out) :: gmat(lmmax,lmmax,nasusc,nasusc)
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ja, lmi, lmj, i6(6)
  real(kind=r8b)    :: x(1000)
  character*60      :: header

!  write(iodb,*) "lmmax=", lmmax
  read(iomain,'(a)') header
!  write(iodb,*) header
  do ja=1,nasusc
    do ia=1,nasusc
      read(iomain,*) i6
!      write(iodb,'(6i4)') i6
      do lmj=1,lmmax
        read(iomain,*) x(1:2*lmmax)
        do lmi=1,lmmax
          gmat(lmi,lmj,ia,ja) = cmplx(x(2*lmi-1),x(2*lmi))
!          write(iodb,'(2i4,6es12.4)') lmi, lmj, gmat(lmi,lmj,ia,ja)
        end do
      end do
    end do
  end do
! All done
  end subroutine in_gmat
