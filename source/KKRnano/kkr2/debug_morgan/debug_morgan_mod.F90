!> Electrostatics test with Morgan charge distribution.
!>
!> Morgan charge distribution:
!> \rho(\vec{r}) = \sum_G \cos( \vec{G} \vec{r} ) where G is sum over first shell of reciprocal lattice vectors
!>               = Re \sum_G \exp( i \vec{G} \vec{r} )
!>
!> This charge distribution is lattice periodic and expansions + potential are known analytically
!> The potential is 8*PI/G**2 * rho + integration constant
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
  
  !------------------------------------------------------------------------------
  !> Evaluate Morgan potential at 'vec'.
  function eval_morgan_potential(reciprocal, vec) result(val)
    double precision, intent(in) :: reciprocal(:,:) !< reciprocal lattice vectors of first shell
    double precision, intent(in) :: vec(3) !< position to evaluate potential

    double precision val

    val = 8 * PI / norm2(reciprocal(:,1))**2 * eval_morgan_rho(reciprocal, vec)
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
  
!================== generalised Morgan test charge for lattice+basis ==========

  !------------------------------------------------------------------------------
  !> Evaluate generalised Morgan charge density
  function eval_gen_morgan_rho(reciprocal, vec, prefactors, rbasis, center) result(val)
    double precision, intent(in) :: reciprocal(:,:) !< reciprocal lattice vectors of first shell
    double precision, intent(in) :: vec(3) !< position to evaluate potential
    double complex  , intent(in) :: prefactors(:)
    double precision, intent(in) :: rbasis(:,:)
    double precision, intent(in) :: center(3)

    integer g_ind, ii
    double complex val
    double complex, parameter :: IMAGINARY = (0.0d0, 1.0d0)

    val = dcmplx(0.0d0, 0.0d0)

    do g_ind = 1, size(reciprocal, 2)
      do ii = 1, size(rbasis, 2)
        val = val + prefactors(ii) * exp(IMAGINARY * dot_product(reciprocal(:,g_ind), (vec + center - rbasis(:,ii))))
      enddo
    enddo

  end function

  !------------------------------------------------------------------------------
  !> Evaluate generalised Morgan potential at 'vec'.
  function eval_gen_morgan_potential(reciprocal, vec, prefactors, rbasis, center) result(val)
    double precision, intent(in) :: reciprocal(:,:) !< reciprocal lattice vectors of first shell
    double precision, intent(in) :: vec(3) !< position to evaluate potential
    double complex  , intent(in) :: prefactors(:)
    double precision, intent(in) :: rbasis(:,:)
    double precision, intent(in) :: center(3)

    double complex val

    val = 8 * PI / norm2(reciprocal(:,1))**2 * eval_gen_morgan_rho(reciprocal, vec, prefactors, rbasis, center)
  end function

  !----------------------------------------------------------------------------
  !> Calculate the spherical harmonic expansion of \exp(i \vec{G} \vec{r}).
  subroutine calc_exponential_expansion(coeffs, g_vector, radius, lmax)
    double complex, intent(inout) :: coeffs(:)
    double precision, intent(in) :: g_vector(3) !< reciprocal lattice vectors of first shell
    double precision, intent(in) :: radius !< position to evaluate potential
    integer, intent(in) :: lmax

    double complex, parameter :: IMAGINARY = (0.0d0, 1.0d0)
    double precision norm_G, arg, temp
    double precision, allocatable :: bessels(:,:)
    double precision, allocatable :: ylm(:)
    integer nm, L, M, cnt

    norm_G = norm2(g_vector)
    arg = norm_G * radius

    allocate(bessels(0:lmax,2))
    allocate(ylm((lmax+1)**2))

    ! spherical bessel functions up to lmax needed
    call SPHJ(lmax,arg,nm, bessels(:,1), bessels(:,2))

    call YMY(g_vector(1), g_vector(2), g_vector(3), temp, ylm, lmax)

    cnt = 1
    do L = 0, lmax
      do M = -L, L
        coeffs(cnt) = 4 * PI * IMAGINARY**L * bessels(L, 1) * ylm(cnt)
        cnt = cnt + 1
      enddo
    enddo

    deallocate(bessels)
    deallocate(ylm)

  end subroutine

  !----------------------------------------------------------------------------
  !> calculate structure factor for basis 'rbasis' for cell centered at center.
  !>
  !> S_j(\vec{G}) = \sum_i c_i \exp(i \vec{G} (\vec{R_j} - \vec{R_i})
  !> c_i ... prefactors
  !> R_j ... center
  !> R_i ... rbasis
  double complex function structure_factor(g_vector, prefactors, rbasis, center)
    double precision, intent(in) :: g_vector(3)
    double complex  , intent(in) :: prefactors(:)
    double precision, intent(in) :: rbasis(:,:)
    double precision, intent(in) :: center(3)

    double complex, parameter :: IMAGINARY = (0.0d0, 1.0d0)
    integer :: ii

    structure_factor = dcmplx(0.0d0, 0.0d0)
    do ii = 1, size(rbasis, 2)
      structure_factor = structure_factor + prefactors(ii) * &
                                            exp(IMAGINARY * dot_product(g_vector, (center - rbasis(:, ii))))
    enddo
  end function

  !----------------------------------------------------------------------------
  !> Calculate the spherical harmonic expansion of generalised Morgan charge density.
  !>
  !> \rho_(r_j) = \sum_{G first reciprocal shell} \sum_i c_i \exp(i G (r_j + R_j - R_i)
  !> Cell centered coordinates: r_j centered at R_j ('center')
  subroutine calc_gen_morgan_rho_expansion(coeffs, reciprocal, prefactors, rbasis, center, radius, lmax)
    double complex, intent(inout) :: coeffs(:)
    double precision, intent(in) :: reciprocal(:,:) !< reciprocal lattice vectors of first shell
    double complex  , intent(in) :: prefactors(:)
    double precision, intent(in) :: rbasis(:,:)
    double precision, intent(in) :: center(3)
    double precision, intent(in) :: radius !< position to evaluate potential
    integer, intent(in) :: lmax

    integer :: ii
    double complex, allocatable :: temp_coeffs(:)
    double complex :: factor

    coeffs = dcmplx(0.0d0, 0.0d0)
    allocate(temp_coeffs, source = coeffs)

    do ii = 1, size(reciprocal, 2)
      call calc_exponential_expansion(temp_coeffs, reciprocal(:,ii), radius, lmax)
      factor = structure_factor(reciprocal(:,ii), prefactors, rbasis, center)
      coeffs = coeffs + temp_coeffs * factor
    end do

    deallocate(temp_coeffs)
  end subroutine

  !----------------------------------------------------------------------------
  !> get first shell reciprocal lattice vectors of simple cubic Bravais lattice
  subroutine reciprocal_cubic(reciprocals)
    double precision, intent(out) :: reciprocals(3, 6)

    reciprocals(:, 1) = (/  2*PI, 0.0d0, 0.0d0 /)
    reciprocals(:, 2) = (/ -2*PI, 0.0d0, 0.0d0 /)
    reciprocals(:, 3) = (/  0.0d0,  2*PI, 0.0d0 /)
    reciprocals(:, 4) = (/  0.0d0, -2*PI, 0.0d0 /)
    reciprocals(:, 5) = (/  0.0d0, 0.0d0,  2*PI /)
    reciprocals(:, 6) = (/  0.0d0, 0.0d0, -2*PI /)

  end subroutine

end module

