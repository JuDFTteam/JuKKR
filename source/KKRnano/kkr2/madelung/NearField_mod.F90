#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

module NearField_mod
  use MadelungCalculator_mod
  implicit none

  type IntracellPot
    integer :: dummy
  end type

  contains

  !----------------------------------------------------------------------------
  !> V_(lm)(r) = ac_wrong(lm) * (-r)**l
  subroutine calc_wrong_contribution_coeff(ac_wrong, dist_vec, dfac, &
                                           charge_mom_total, gaunt)
    implicit none
    double precision, intent(inout) :: ac_wrong(:)
    double precision, intent(in)    :: dist_vec(3)
    double precision, intent(in)    :: dfac(:,:)
    double precision, intent(in)    :: charge_mom_total(:)
    type (MadelungClebschData), intent(in) :: gaunt

    integer :: lmx, L, m, ii, lmmaxd, L1, L2, lm1, lm2, lm3
    double precision :: r, rfac
    double precision :: pi
    double precision, allocatable :: ylm(:)
    double precision, allocatable :: smat(:)
    double precision, allocatable :: avmad(:,:)

    pi = 4.0d0*atan(1.0d0)

    lmx = size(dfac, 1)
    lmmaxd = (lmx + 1) ** 2
    allocate(ylm(lmx))
    allocate(smat(lmx))
    allocate(avmad(lmmaxd, lmmaxd))

    call ymy(dist_vec(1),dist_vec(2),dist_vec(3),r,ylm,lmx)

    do L = 0,lmx
       rfac = (1. / (r**(L+1))) / sqrt(pi)
       do m = -L,L
          lm1 = L*(L+1) + m + 1
          smat(lm1) = ylm(lm1)*rfac
       end do
    end do

    avmad = 0.0d0

    do ii = 1,gaunt%IEND
      LM1 = gaunt%ICLEB(ii,1)
      LM2 = gaunt%ICLEB(ii,2)
      LM3 = gaunt%ICLEB(ii,3)
      L1 = gaunt%LOFLM(LM1)
      L2 = gaunt%LOFLM(LM2)

      ! --> this loop has to be calculated only for l1+l2=l3
      ! that means it contains only 1 term !!!!
      ! L1, L2 are given ===> L3 = L1+L2

      AVMAD(LM1,LM2) = AVMAD(LM1,LM2) + &
                       2.0D0*DFAC(L1,L2)*SMAT(LM3)*gaunt%CLEB(ii)
    end do

    ac_wrong = matmul(avmad, charge_mom_total)

  end subroutine

!  !----------------------------------------------------------------------------
!  subroutine test_calc_wrong_contribution_coeff(madelung_calc, dist_vec, cmom, cmom_inst)
!    implicit none
!    type (MadelungCalculator) :: madelung_calc
!    double precision :: dist_vec(3)
!    double precision :: cmom(:)
!    double precision :: cmom_inst(:)
!
!    integer :: lmmaxd
!    double precision :: ac_wrong(:)
!    double precision :: cmom_total(:)
!
!    lmmaxd = size(cmom)
!    allocate(ac_wrong(lmmaxd))
!    allocate(cmom_total(lmmaxd))
!
!    cmom_total = cmom + cmom_inst
!    call calc_wrong_contribution_coeff(ac_wrong, dist_vec, madelung_calc%dfac, &
!                                          cmom_total, madelung_calc%clebsch)
!
!    write(*,*) ac_wrong
!  end subroutine

  !----------------------------------------------------------------------------
  !> TODO: interpolation/calculation of v_intra
  subroutine get_intracell(v_intra, radius, pot)
    implicit none
    double precision, intent(out) :: v_intra(:)
    double precision, intent(in) :: radius
    type(IntracellPot), intent(in) :: pot

    integer :: ii

    do ii = 1, size(v_intra)
      v_intra(ii) = radius * ii
    end do
  end subroutine

  !----------------------------------------------------------------------------
  !> TODO: this is a general routine for shifting sph. harm. expansions
  subroutine calc_near_field(v_near, radius, dist_vec, pot, lmax_prime)
    implicit none
    double precision, intent(out) :: v_near(:)  !indices (lm)
    double precision, intent(in) :: radius
    double precision, intent(in) :: dist_vec(3)
    type(IntracellPot), intent(in) :: pot
    integer, intent(in) :: lmax_prime

    integer :: lmmaxd, lmax, lmmaxd_prime
    integer, parameter :: NUM_LEBEDEV = 434
    double precision :: v_leb(3)
    double precision :: vec(3)
    double precision :: weight_leb(NUM_LEBEDEV)
    double precision, allocatable :: sph_harm_leb(:)
    double precision, allocatable :: sph_harm(:, :)
    double precision, allocatable :: v_intra_leb(:)
    double precision, allocatable :: integrand(:,:)
    double precision, allocatable :: temp(:)
    double precision :: norm_vec
    double precision :: dummy
    double precision :: FOUR_PI ! prefactor for Lebedev: 4*pi

    integer :: ij, lm

    FOUR_PI = 16.0d0 * atan(1.0d0)

    lmmaxd = size(v_near)
    lmax = int(sqrt(dble(lmmaxd) + 0.1) - 1)
    lmmaxd_prime = (lmax_prime+1)**2

    CHECKASSERT( (lmax + 1)**2 == lmmaxd )

    allocate(sph_harm_leb(lmmaxd_prime))
    allocate(v_intra_leb(lmmaxd_prime))
    allocate(sph_harm(NUM_LEBEDEV, lmmaxd))
    allocate(integrand(NUM_LEBEDEV, lmmaxd))
    allocate(temp(lmmaxd))

    do ij = 1, NUM_LEBEDEV

      call LEBEDEV(ij, v_leb(1), v_leb(2), v_leb(3), weight_leb(ij))

      call ymy(v_leb(1), v_leb(2), v_leb(3), dummy, sph_harm_leb, lmax_prime)

      vec = radius * v_leb + dist_vec
      call ymy(vec(1), vec(2), vec(3), norm_vec, temp, lmax)
      do lm = 1, lmmaxd
        sph_harm(ij, lm) = temp(lm)
      end do

      call get_intracell(v_intra_leb, norm_vec, pot)

      ! perform summation over L'
      integrand(ij,1) = dot_product(sph_harm_leb, v_intra_leb)
    end do

    integrand(:,1) = integrand(:,1) * weight_leb * FOUR_PI

    do lm = 1, lmmaxd
      integrand(:,lm) = integrand(:,1)
    end do

    integrand = integrand * sph_harm

    v_near = sum(integrand, 1)

  end subroutine

  !----------------------------------------------------------------------------
  ! Tests orthgonality of real spherical harmonics up to LMAX
  ! The first column of 'integrand' must be 1.0 and the others = 0.0
  subroutine test_lebedev()
    implicit none
    integer, parameter :: LMAX = 2
    integer, parameter :: LMMAXD = (LMAX+1)**2

    real(8) :: integrand(LMMAXD, LMMAXD)
    real(8) :: ylm(LMMAXD)
    real(8) :: vnorm
    real(8) :: vec(3)
    real(8) :: weight
    real(8) :: FOUR_PI
    integer ij, lm1, lm2

    FOUR_PI = 16.0d0 * atan(1.0d0)

    integrand = 0.0d0
    do ij = 1, 434
      call LEBEDEV(ij, vec(1), vec(2), vec(3), weight)
      call YMY(vec(1), vec(2), vec(3), vnorm, ylm, LMAX)
      do lm1 = 1, LMMAXD
        do lm2 = 1, LMMAXD
        integrand(lm2,lm1) = integrand(lm2, lm1) + ylm(lm2) * ylm(mod((lm2 + lm1 - 2), LMMAXD)+1) * weight * FOUR_PI
        end do
      end do
    end do

    write(*,*) integrand
  end subroutine

end module NearField_mod

!program test_it
!  use NearField_mod
!  implicit none
!  double precision v_near(9)
!  type (IntracellPot) :: pot
!  double precision :: d(3) = (/0.00d0, 0.00d0, -0.05d0 /)
!  call test_lebedev()
!  call calc_near_field(v_near, 1.0d0, d , pot, 2)
!  write(*,*) v_near
!end program
