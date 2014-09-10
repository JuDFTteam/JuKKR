!> Electrostatics test with Morgan charge distribution.
!>
!> Morgan charge distribution:
!> \rho(\vec{r}) = \sum_G \cos( \vec{G} \vec{r} ) where G is sum over first shell of reciprocal lattice vectors
!>               = Re \sum_G \exp( i \vec{G} \vec{r} )
!>
!> This charge distribution is lattice periodic and expansions + potential are known analytically
!>
!> @author Elias Rabel

module debug_morgan_mod
  implicit none
  
  double precision, parameter, private :: PI = 4.0d0 * atan(1.0d0)
  
  contains
  
  !------------------------------------------------------------------------------
  !> Evaluate Morgan charge density
  function eval_morgan_rho(reciprocal, vec) result(val)
    double precision, intent(in) :: reciprocal(:,:) !< reciprocal lattice vectors of first shell
    double precision, intent(in) :: vec(3) !< position to evaluate potential
    
    integer ii
    double precision val
    
    val = 0.0d0
    do ii = 1, size(reciprocal, 2)
      val = val + cos(dot_product(reciprocal(:,ii), vec))
    end do
  end function
  
  !> Calculate the spherical harmonic expansion of Morgan charge density
  subroutine calc_morgan_rho_expansion(coeffs, reciprocal, radius, lmax)
    double precision, intent(inout) :: coeffs(:)
    double precision, intent(in) :: reciprocal(:,:) !< reciprocal lattice vectors of first shell
    double precision, intent(in) :: radius !< position to evaluate potential
    integer, intent(in) :: lmax
    
    double precision norm_G, arg, temp
    double precision, allocatable :: bessels(:,:)
    double precision, allocatable :: ylm(:), sum_ylm(:)
    integer nm, L, M, cnt, ii
    
    norm_G = norm2(reciprocal(:,1))
    arg = norm_G * radius
    
    allocate(bessels(0:lmax,2))
    allocate(ylm((lmax+1)**2))
    allocate(sum_ylm((lmax+1)**2))
    
    ! spherical bessel functions up to lmax needed
    call SPHJ(lmax,arg,nm, bessels(:,1), bessels(:,2))
    
    sum_ylm = 0.0d0
    do ii = 1, size(reciprocal, 2)
      call YMY(reciprocal(1,ii),reciprocal(2,ii),reciprocal(3,ii), temp, ylm, lmax)
      sum_ylm = sum_ylm + ylm
    end do
    
    cnt = 1
    do L = 0, lmax
      do M = -L, L
        if (mod(L,2) == 0) then
          coeffs(cnt) = 4 * PI * (-1)**(L/2) * bessels(L, 1) * sum_ylm(cnt)
        else
          coeffs(cnt) = 0.0d0 ! expansion for odd L is 0 since charge distribution has even parity 
        endif
        cnt = cnt + 1
      enddo
    enddo
  end subroutine
  
  !----------------------------------------------------------------------------
  !> evaluate spherical harmonic expansion at angles given by 'vec'.
  double precision function eval_expansion(coeffs, vec)
    implicit none
    double precision :: coeffs(:)
    double precision :: vec(3)

    integer lmmaxd
    integer lmax
    double precision, allocatable :: ylm(:)
    double precision :: vnorm

    lmmaxd = size(coeffs)
    lmax = int(sqrt(dble(lmmaxd) + 0.1) - 1)

    allocate(ylm(lmmaxd))
    call YMY(vec(1), vec(2), vec(3), vnorm, ylm, LMAX)

    eval_expansion = dot_product(coeffs, ylm)
  end function
  
  !> get first shell reciprocal lattice vectors of fcc lattice (the reciprocal lattice is bcc)
  subroutine reciprocal_fcc(reciprocals)
    double precision, intent(out) :: reciprocals(3, 8)
    
    integer :: sx, sy, sz, cnt
    
    cnt = 1
    do sx = -1, 1, 2
      do sy = -1, 1, 2
        do sz = -1, 1, 2
          reciprocals(1, cnt) = 2 * PI * sx
          reciprocals(2, cnt) = 2 * PI * sy
          reciprocals(3, cnt) = 2 * PI * sz
          cnt = cnt + 1
        enddo
      enddo
    enddo
    
  end subroutine
  
  !> Morgan charge density in direction 'dir' in a fcc lattice
  subroutine print_morgan_rho_dir(dir, npoints)
    double precision :: dir(3), vec(3), reciprocals(3,8), val
    integer npoints
    
    integer :: ii

    call reciprocal_fcc(reciprocals)
    
    do ii = 1, NPOINTS
      vec = dir / (NPOINTS-1) * (ii-1)
      val = eval_morgan_rho(reciprocals, vec)
      print *, norm2(vec), val
    enddo
  end subroutine
  
  subroutine print_morgan_rho_exp_dir(dir, lmax, npoints)
    
    double precision :: dir(3), vec(3), reciprocals(3,8), val
    integer :: lmax
    integer :: npoints
    
    integer :: ii

    double precision :: coeffs((lmax+1)**2)

    call reciprocal_fcc(reciprocals)
    
    do ii = 1, NPOINTS
      vec = dir / (NPOINTS-1) * (ii-1)
      if (norm2(vec) == 0.0d0) vec(1) = 1e-6
      call calc_morgan_rho_expansion(coeffs, reciprocals, norm2(vec), lmax)
      val = eval_expansion(coeffs, vec)
      print *, norm2(vec), val
    enddo
  end subroutine
  
end module

