! ************************************************************************
subroutine spin_orbit_one_l(lmax, l_s)
! ************************************************************************
!      in this subroutine the matrix L*S is calculated for the basis of
!      real spherical harmonics

!      schematically it has the form
!      (  -L_z    L_+  )
!      (  L_-     L_z  )
  implicit none

  integer, intent (in) :: lmax
  double complex, intent (out) :: l_s((2*lmax+1)*2, (2*lmax+1)*2)

!  local variables 
  integer :: i1, i2, i1l
  double complex :: icompl
  double complex, allocatable :: l_min(:, :)
  double complex, allocatable :: l_up(:, :)
  double precision :: lfac


  icompl = (0d0, 1d0)


  allocate (l_min(-lmax:lmax,-lmax:lmax))
  allocate (l_up(-lmax:lmax,-lmax:lmax))

!  initialize the matrix

  do i1 = 1, (2*lmax+1)*2
    do i2 = 1, (2*lmax+1)*2
      l_s(i2, i1) = 0d0
    end do
  end do

  do i1 = -lmax, lmax
    do i2 = -lmax, lmax
      l_min(i2, i1) = 0d0
      l_up(i2, i1) = 0d0
    end do
  end do

!  fill the second and the forth quadrant with L_z
! (-L_z,respectively)


  do i1 = 1, 2*lmax + 1
    i1l = i1 - lmax - 1 ! the value of m (varies from -l to +l)
    i2 = 2*lmax + 1 - (i1-1)

!         L_S(i2,i1)=icompl*i1l
    l_s(i2, i1) = -icompl*i1l

  end do

  do i1 = 2*lmax + 2, (2*lmax+1)*2
    i1l = i1 - lmax - 1 - (2*lmax+1) ! the value of m (varies from -l to +l)
    i2 = (2*lmax+1)*2 - (i1-(2*lmax+2))

!         L_S(i2,i1)=-icompl*i1l
    l_s(i2, i1) = icompl*i1l

  end do


!  implement now L_- in the third quadrant

  if (lmax>0) then

    lfac = sqrt(lmax*(lmax+1d0))/sqrt(2d0)
    l_min(0, -1) = -icompl*lfac
!         l_min(0,-1)=icompl*lfac
    l_min(0, 1) = lfac
    l_min(-1, 0) = icompl*lfac
    l_min(1, 0) = -lfac

    if (lmax>1) then

      do i1 = 2, lmax

        lfac = 0.5d0*sqrt(lmax*(lmax+1d0)-i1*(i1-1d0))
        l_min(-i1, -i1+1) = -lfac
        l_min(-i1, i1-1) = icompl*lfac
        l_min(i1, -i1+1) = -icompl*lfac
        l_min(i1, i1-1) = -lfac

        lfac = 0.5d0*sqrt(lmax*(lmax+1d0)-(i1-1)*(i1))
        l_min(-i1+1, -i1) = lfac
        l_min(-i1+1, i1) = icompl*lfac
        l_min(i1-1, -i1) = -icompl*lfac
        l_min(i1-1, i1) = lfac

      end do

    end if
  end if


  do i1 = -lmax, lmax
    do i2 = -lmax, lmax
      l_s(i2+3*lmax+2, i1+lmax+1) = l_min(i1, i2)
    end do
  end do


!  implement now L_+ in the   quadrant

  if (lmax>0) then

    lfac = sqrt(lmax*(lmax+1d0))/sqrt(2d0)
    l_up(0, -1) = -icompl*lfac
    l_up(0, 1) = -lfac
    l_up(-1, 0) = icompl*lfac
    l_up(1, 0) = lfac

    if (lmax>1) then

      do i1 = 2, lmax

        lfac = 0.5d0*sqrt(lmax*(lmax+1d0)-i1*(i1-1d0))
        l_up(-i1, -i1+1) = lfac
        l_up(-i1, i1-1) = icompl*lfac
        l_up(i1, -i1+1) = -icompl*lfac
        l_up(i1, i1-1) = lfac

        lfac = 0.5d0*sqrt(lmax*(lmax+1d0)-(i1-1)*(i1))
        l_up(-i1+1, -i1) = -lfac
        l_up(-i1+1, i1) = icompl*lfac
        l_up(i1-1, -i1) = -icompl*lfac
        l_up(i1-1, i1) = -lfac

      end do

    end if
  end if


  do i1 = -lmax, lmax
    do i2 = -lmax, lmax
      l_s(i2+lmax+1, i1+3*lmax+2) = l_up(i1, i2)
    end do
  end do



  deallocate (l_min)
  deallocate (l_up)


end subroutine
