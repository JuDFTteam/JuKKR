  subroutine sort_energies(n0,n,z0,zin,iout)
! Finds the closest n0 values in zin to z0, returns locations in iout

  implicit none

! How many points to output
  integer(kind=i4b), intent(in)  :: n0
! How many points to sort
  integer(kind=i4b), intent(in)  :: n
! Target point
  complex(kind=c8b), intent(in)  :: z0
! Points to sort
  complex(kind=c8b), intent(in)  :: zin(n)
! Location of sorted points
  integer(kind=i4b), intent(out) :: iout(n0)
! -----------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-8
  integer(kind=i4b) :: k(n), i, j, info
  real(kind=r8b)    :: din(n), dout(n)

!  write(iodb,'("sort_energies: n,n0,z0=",2i4,2f16.8)') n, n0, z0
! Prepare data
  din = abs(zin - z0)
  dout = din
! MKL sorting routine (can be replaced by quicksort or heapsort)
!  call dlasrt2('I',n,d,k,info)  ! needs SCALAPACK apparently
  call dlasrt('I',n,dout,info)
  do j=1,n0
    do i=1,n
      if (abs(din(i) - dout(j)) < tol) then
        iout(j) = i
        exit
      end if
    end do
  end do
!  write(iodb,'("sort_energies: zout=",1000f16.8)') zin(iout)
! All done!
  end subroutine sort_energies

