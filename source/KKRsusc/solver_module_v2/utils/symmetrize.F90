  subroutine symmetrize(n,a,tol)
! Compare all elements of a with each other

  implicit none

  integer(kind=i4b), intent(in)    :: n
  complex(kind=c8b), intent(inout) :: a(n,n)
  real(kind=r8b),    intent(in)    :: tol
! ---------------------------------------
  integer(kind=i4b) :: i1, j1, i2, j2
  real(kind=r8b)    :: re1, im1, re2, im2, avg

  do j2=1,n
  do i2=1,n
    do j1=1,n
    do i1=1,n
      re2 = real(a(i2,j2)); im2 = aimag(a(i2,j2))
      re1 = real(a(i1,j1)); im1 = aimag(a(i1,j1))
      if (abs(re1) < tol) re1 = 0.d0
      if (abs(im1) < tol) im1 = 0.d0
      if (abs(re2) < tol) re2 = 0.d0
      if (abs(im2) < tol) im2 = 0.d0
!      if (abs(abs(re1) - abs(re2)) < tol) then
!        avg = 0.5d0*(abs(re1) + abs(re2))
!        re1 = sign(avg,re1)
!        re2 = sign(avg,re2)
!      end if
!      if (abs(abs(im1) - abs(im2)) < tol) then
!        avg = 0.5d0*(abs(im1) + abs(im2))
!        im1 = sign(avg,im1)
!        im2 = sign(avg,im2)
!      end if
!      if (abs(abs(re1) - abs(im2)) < tol) then
!        avg = 0.5d0*(abs(re1) + abs(im2))
!        re1 = sign(avg,re1)
!        im2 = sign(avg,im2)
!      end if
!      if (abs(abs(im1) - abs(re2)) < tol) then
!        avg = 0.5d0*(abs(im1) + abs(re2))
!        im1 = sign(avg,im1)
!        re2 = sign(avg,re2)
!      end if
      a(i1,j1) = cmplx(re1,im1)
      a(i2,j2) = cmplx(re2,im2)
    end do
    end do
  end do
  end do
! All done!
  end subroutine symmetrize

