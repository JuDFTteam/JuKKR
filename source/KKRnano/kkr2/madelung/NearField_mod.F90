
!> @author Elias Rabel

module NearField_mod
#include "macros.h"
  use Constants_mod, only: pi
  implicit none
  private
  public :: calculate
  
  interface calculate
    module procedure calc_near_field, calc_wrong_contribution_coeff
  endinterface

  contains

  !----------------------------------------------------------------------------
  !> V_(lm)(r) = ac_wrong(lm) * (-r)**l
  subroutine calc_wrong_contribution_coeff(ac_wrong, dist_vec, charge_mom_total, Gaunt)
    use MadelungCalculator_mod, only: MadelungClebschData, calculate!createDfac
    use Harmonics_mod, only: ymy
    double precision, intent(inout) :: ac_wrong(:)
    double precision, intent(in)    :: dist_vec(3)

    double precision, intent(in)    :: charge_mom_total(:)
    type(MadelungClebschData), intent(in) :: Gaunt

    integer :: lmx, l, m, ii, lmmaxd, l1, l2, lm1, lm2, lm3, lmx_prime, lmmaxd_prime
    double precision :: r, rfac
    double precision, allocatable :: ylm(:), smat(:), avmad(:,:), dfac(:,:)

    lmmaxd = size(ac_wrong)
    lmx = int(sqrt(dble(lmmaxd) + 0.1) - 1)

    ! the inner summation should run up to 2*lmx
    ! this follows from Gaunt properties
    lmx_prime = 2*lmx
    lmmaxd_prime = (2*lmx+1)**2

    allocate(ylm(lmmaxd_prime), smat(lmmaxd_prime), avmad(lmmaxd,lmmaxd))

    call calculate(dfac, lmx) ! createDfac: calculate complicated prefactor, will be allocated here

    call ymy(dist_vec(1), dist_vec(2), dist_vec(3), r, ylm, lmx_prime)

    do l = 0, lmx_prime
      rfac = r**(-l-1)
      do m = -l, l
        lm1 = l*l + l + m + 1
        smat(lm1) = ylm(lm1)*rfac
      enddo ! m
    enddo ! l
!     or equivalent    
!     do lm1 = 1, lmmaxd_prime
!       smat(lm1) = ylm(lm1) * r**(-Gaunt%loflm(lm1)-1)
!     enddo ! lm1

    avmad = 0.d0

    do ii = 1, Gaunt%iend
      lm1 = Gaunt%icleb(ii,1)
      lm2 = Gaunt%icleb(ii,2)
      lm3 = Gaunt%icleb(ii,3)
      l1 = Gaunt%loflm(lm1)
      l2 = Gaunt%loflm(lm2)
      ! --> Gaunt property: l1 + l2 <= l3

      if (lm1 <= lmmaxd .and. lm2 <= lmmaxd .and. lm3 <= lmmaxd_prime) &
        avmad(lm1,lm2) = avmad(lm1,lm2) + 2.d0*dfac(l1,l2)*smat(lm3)*Gaunt%cleb(ii) ! factor 2 comes from the use of Rydberg units
    enddo ! ii

    ac_wrong = matmul(avmad, charge_mom_total)

    deallocate(ylm, smat, avmad, dfac, stat=ii)
  endsubroutine ! calc

  !----------------------------------------------------------------------------
  !> this is a general routine for shifting sph. harm. expansions
  subroutine calc_near_field(v_near, radius, dist_vec, pot, lmax_prime)
    use NearField_kkr_mod, only: IntracellPotential, get
    use Harmonics_mod, only: ymy
    double precision, intent(out) :: v_near(:)  !indices (lm)
    double precision, intent(in) :: radius
    double precision, intent(in) :: dist_vec(3)
    type(IntracellPotential), intent(inout) :: pot
    integer, intent(in), optional :: lmax_prime

    integer :: lmmaxd, lmax, lmmaxd_prime
    integer, parameter :: NUM_LEBEDEV = 434
    double precision :: v_leb(3), vec(3)
    double precision :: weight_leb(NUM_LEBEDEV)
    double precision, allocatable :: sph_harm_leb(:,:)
    double precision, allocatable :: sph_harm(:)
    double precision, allocatable :: v_intra(:)
    double precision, allocatable :: integrand(:,:)
    double precision, allocatable :: temp(:)
    double precision :: norm_vec, dummy, four_pi ! prefactor for Lebedev: 4*pi

    integer :: ij, lm
    integer :: lmax_p

    four_pi = 4.d0*pi

    lmmaxd = size(v_near)
    lmax = int(sqrt(dble(lmmaxd) + 0.1) - 1)
    
    if (.not. present(lmax_prime)) lmax_p = lmax
    
    lmmaxd_prime = (lmax_p+1)**2

    CHECKASSERT( (lmax + 1)**2 == lmmaxd )

    allocate(sph_harm(lmmaxd_prime), v_intra(lmmaxd_prime), temp(lmmaxd))
    allocate(sph_harm_leb(NUM_LEBEDEV,lmmaxd), integrand(NUM_LEBEDEV,lmmaxd))

    do ij = 1, NUM_LEBEDEV

      call LEBEDEV(ij, v_leb(1), v_leb(2), v_leb(3), weight_leb(ij))

      call ymy(v_leb(1), v_leb(2), v_leb(3), dummy, temp, lmax)
      do lm = 1, lmmaxd
        sph_harm_leb(ij,lm) = temp(lm)
      enddo ! lm

      vec(1:3) = radius*v_leb(1:3) + dist_vec(1:3)

      if (vec(1)**2 + vec(2)**2 + vec(3)**2 > 0.d0) then  ! can be zero
        call ymy(vec(1), vec(2), vec(3), norm_vec, sph_harm, lmax_p)
      else
        norm_vec    = 0.d0
        sph_harm(:) = 0.d0
        sph_harm(1) = 1.d0/sqrt(four_pi)
      endif

#ifdef TEST_POTENTIALS
!     call get_const_monopole(v_intra, norm_vec)
      call get_const_multipole(v_intra, norm_vec)
#else
      call get(pot, v_intra, norm_vec) ! get intracell potential at radius 'norm_vec'
#endif

      integrand(ij,1) = dot_product(sph_harm, v_intra) ! perform summation over L'
    enddo ! ij

    integrand(:,1) = integrand(:,1) * weight_leb * four_pi

    do lm = 1, lmmaxd
      integrand(:,lm) = integrand(:,1)
    enddo ! lm

    integrand = integrand * sph_harm_leb

    v_near = sum(integrand, dim=1)

  endsubroutine ! calc

  !----------------------------------------------------------------------------
  ! Tests orthgonality of real spherical harmonics up to LMAX
  ! The first column of 'integrand' must be 1.0 and the others = 0.0
  subroutine test_lebedev()
    use Harmonics_mod, only: ymy
    integer, parameter :: lmax = 2, lmmaxd = (lmax+1)**2, NUM_LEBEDEV = 434

    double precision :: integrand(lmmaxd,lmmaxd)
    double precision :: ylm(lmmaxd)
    double precision :: vnorm, vec(3), weight, four_pi
    integer ij, lm1, lm2

    four_pi = 4.d0*pi
    integrand = 0.d0
    do ij = 1, NUM_LEBEDEV
      call lebedev(ij, vec(1), vec(2), vec(3), weight)
      call ymy(vec(1), vec(2), vec(3), vnorm, ylm, lmax)
      do lm1 = 1, lmmaxd
        do lm2 = 1, lmmaxd
          integrand(lm2,lm1) = integrand(lm2,lm1) + ylm(lm2) * ylm(mod((lm2 + lm1 - 2), lmmaxd)+1) * weight * four_pi
        enddo ! lm2
      enddo ! lm1
    enddo ! ij

    write(*,*) integrand
  endsubroutine ! test

  !----------------------------------------------------------------------------
  ! evaluate spherical harmonic expansion at angles given by 'vec'.
  double precision function eval_expansion(coeffs, vec)
    use Harmonics_mod, only: ymy
    double precision :: coeffs(:)
    double precision :: vec(3)

    integer :: lmmaxd, lmax
    double precision, allocatable :: ylm(:)
    double precision :: vnorm

    lmmaxd = size(coeffs)
    lmax = int(sqrt(dble(lmmaxd) + 0.1) - 1)

    allocate(ylm(lmmaxd))
    call YMY(vec(1), vec(2), vec(3), vnorm, ylm, LMAX)

    eval_expansion = dot_product(coeffs, ylm)
  endfunction ! eval


#ifdef TEST_POTENTIALS
!+never

!------------------------------------------------------------------------------
! Some test potentials that replace 
!------------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  !> A test potential: constant multipole moments
  subroutine get_const_multipole(v_intra, radius)
    ! Test potential: assume multipoles Q_L = 1.d0
    double precision, intent(out) :: v_intra(:)
    double precision, intent(in) :: radius

    integer :: lm, L, M

    lm = 1
    L = 0
    M = 0
    do while (.true.)
      v_intra(lm) = 4.d0*sqrt(pi)/(radius**(L+1)*(2*L+1))
      lm = lm + 1
      M = M + 1
      if (M > L) then
        L = L + 1
        M = -L
      endif
      if (lm > size(v_intra)) exit
    enddo ! while

  endsubroutine ! get

  !----------------------------------------------------------------------------
  !> A test potential: potential of a (unit) monopole
  subroutine get_const_monopole(v_intra, radius)
    double precision, intent(out) :: v_intra(:)
    double precision, intent(in) :: radius

    v_intra = 0.d0
    v_intra(1) = 4.d0*sqrt(pi)/radius
  endsubroutine ! get

!-never
#endif
endmodule ! NearField_mod

!  program test_it
!    use NearField_mod, only: ...
!    use NearField_kkr_mod, only: ...
!    implicit none
!    integer, parameter :: LMAX = 4
!    integer, parameter :: LMMAXD = (LMAX+1)**2
!    integer, parameter :: NPOINTS = 100
!    double precision :: v_near(LMMAXD)
!    
!    type(TestPotentialConstMulti) :: pot
!    type(IntracellPotential) :: intra_pot
!  
!    double precision :: d(3) = [0.d0, 0.d0, 2.d0]
!    double precision :: vec(3)
!    double precision :: ac_wrong(LMMAXD)
!    double precision :: v_mad_wrong(LMMAXD)
!    double precision :: cmom(LMMAXD)
!  
!    double precision :: radius = 1.4d0
!    double precision :: r_temp
!  
!    type(MadelungClebschData) :: clebsch
!    integer :: L, M, ii
!  
!    !call test_lebedev()
!    call calc_near_field(v_near, radius, d , pot, LMAX)
!    write(*,*) v_near
!  
!    ! konfus
!    call createMadelungClebschData(clebsch, (4*LMAX+1)**2, (2*LMAX+1)**2)
!    call initMadelungClebschData(clebsch, LMAX)
!  
!    !cmom = 0.d0
!    cmom = 1.d0 / sqrt(16.d0 * atan(1.d0)) 
!    !cmom(3) = 1.d0 / sqrt(16.d0 * atan(1.d0)) 
!  
!    call calc_wrong_contribution_coeff(ac_wrong, d, cmom, clebsch)
!  
!    ii = 1
!    do L = 0, LMAX
!      do M = -L, L
!        v_mad_wrong(ii) = (ac_wrong(ii) * (-radius)**L)
!        write(*,*) L, M, v_near(ii), v_mad_wrong(ii) 
!        ii = ii + 1
!      enddo
!    enddo
!  
!    vec = 0.d0
!    vec(3) = -radius
!    
!    write(*,*) eval_expansion(v_near, vec)
!    write(*,*) eval_expansion(v_mad_wrong, vec)
!    
!    !-------- another test ---------------------------------------------
!    write(*,*) "=========================================================="
!    
!    call intra_pot%create(LMMAXD, NPOINTS)
!    
!    intra_pot%charge_moments = cmom
!    
!    do ii = 1, NPOINTS
!      intra_pot%radial_points(ii) = (ii) * 1.d0 / (NPOINTS)
!      call pot%get(v_near, intra_pot%radial_points(ii))
!      intra_pot%v_intra_values(ii, :) = v_near
!    enddo
!    
!    call intra_pot%init()
! 
!    call calc_near_field(v_near, radius, d , intra_pot, LMAX)
!    write(*,*) v_near
!    write(*,*) eval_expansion(v_near, vec)
! 
! !    do ii = 1, NPOINTS
! !      call pot%get(v_near, intra_pot%radial_points(ii))
! !      write(*,*) v_near(1)
! !      call intra_pot%get(v_near, intra_pot%radial_points(ii))
! !      write(*,*) v_near(1)
! !    enddo
!    
!    !do ii = 1, NPOINTS * 4
!    !  r_temp = (ii) * 2.d0 / (NPOINTS*4) + 1.d0
!    !  call intra_pot%get(v_near, r_temp)
!    !  write(*,*) r_temp, v_near(7), v_near(15)
!    !enddo
!    !call calc_near_field(v_near, radius, d , intra_pot, LMAX)
!    !write(*,*) v_near
!  endprogram
