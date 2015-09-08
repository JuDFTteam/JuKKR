!> tests routine YMY - calculation of spherical harmonics
program test_ymy
  implicit none

  integer, parameter :: LMAX = 32

  double precision :: phi = 2.3
  double precision :: theta = 1.1

  double precision :: ylm( (LMAX+1)**2 )

  double precision :: vec(3)
  double precision :: radius = 1.0

  integer :: LL, MM, count

  vec(1) = radius * sin(theta) * cos(phi)
  vec(2) = radius * sin(theta) * sin(phi)
  vec(3) = radius * cos(theta)

  write(*,*) "cos(theta) = ", cos(theta)
  write(*,*) "cos(phi) = ", cos(phi)
  write(*,*) "sin(phi) = ", sin(phi)

  call YMY(vec(1),vec(2),vec(3),radius,YLM,LMAX) 

  count = 1
  do LL = 0, LMAX
    do MM = -LL, LL
       write(*,*) ylm(count)
       count = count + 1
    end do
  write(*,*) "------------------------"
  end do

end program
