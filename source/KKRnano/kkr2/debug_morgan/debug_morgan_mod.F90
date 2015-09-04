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

  ! On BlueGene/Q norm2 does not exist as an intrinsic - use BLAS
#ifdef __bgq__
  double precision function norm2(vec)
    double precision, intent(in) :: vec(:)
    double precision, external :: dnrm2
    norm2 = dnrm2 (size(vec), vec, 1)
  end function
#endif
    
  !----------------------------------------------------------------------------
  !> evaluate spherical harmonic expansion at angles given by 'vec'.
  double precision function eval_expansion(coeffs, vec)
    use Harmonics_mod, only: ymy
  
    double precision, intent(in) :: coeffs(:)
    double precision, intent(in) :: vec(3)

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
    use Harmonics_mod, only: ymy
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
  !> Calculates cross product.
  function cross_product(vec_a, vec_b)
    double precision :: cross_product(3)
    double precision, intent(in) :: vec_a(3)
    double precision, intent(in) :: vec_b(3)

    cross_product(1) = vec_a(2) * vec_b(3) - vec_a(3) * vec_b(2)
    cross_product(2) = vec_a(3) * vec_b(1) - vec_a(1) * vec_b(3)
    cross_product(3) = vec_a(1) * vec_b(2) - vec_a(2) * vec_b(1)
  end function

  !----------------------------------------------------------------------------
  !> Calculate basis vectors of reciprocal lattice.
  subroutine calc_reciprocal_basis(rec_basis, bravais)
    double precision, intent(out) :: rec_basis(3,3)
    double precision, intent(in)  :: bravais(3,3)

    double precision :: volume

    volume = dot_product(bravais(:,1), cross_product(bravais(:,2), bravais(:,3)))

    rec_basis(:,1) = 2*PI/volume * cross_product(bravais(:,2), bravais(:,3))
    rec_basis(:,2) = 2*PI/volume * cross_product(bravais(:,3), bravais(:,1))
    rec_basis(:,3) = 2*PI/volume * cross_product(bravais(:,1), bravais(:,2))

  end subroutine

  !---------------------------------------------------------------------------
  !> Get first shell of reciprocal lattice vectors.
  subroutine calc_reciprocal_first_shell(reciprocals, rec_basis)
    double precision, allocatable, intent(out) :: reciprocals(:,:)
    double precision, intent(in) :: rec_basis(3,3)

    integer :: na, nb, nc
    double precision :: length, norm, vec(3)
    double precision, parameter :: TOL = 1.d-12
    integer counter

    length = min(norm2(rec_basis(:,1)), &
                 norm2(rec_basis(:,2)), &
                 norm2(rec_basis(:,3)))    ! all vectors in 1st shell have this length

    ! determine first shell rec. vectors - in an inefficient brute force way
    ! count vectors in first shell
    counter = 0
    do na = -1, 1
      do nb = -1, 1
        do nc = -1, 1
          norm = norm2(na * rec_basis(:,1) + nb * rec_basis(:,2) + nc * rec_basis(:,3))
          if (abs(norm - length) < TOL) counter = counter + 1
        enddo
      enddo
    enddo

    allocate(reciprocals(3,counter))

    counter = 0
    do na = -1, 1
      do nb = -1, 1
        do nc = -1, 1
          vec = na * rec_basis(:,1) + nb * rec_basis(:,2) + nc * rec_basis(:,3)
          norm = norm2(vec)
          if (abs(norm - length) < TOL) then
            counter = counter + 1
            reciprocals(:,counter) = vec
          endif
        enddo
      enddo
    enddo
  end subroutine

end module

