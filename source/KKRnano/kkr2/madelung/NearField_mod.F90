#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!> @author Elias Rabel

module NearField_mod
  implicit none
  private
  public :: Potential
  public :: calc_near_field, calc_wrong_contribution_coeff
  
  !> Objects of types derived from 'Potential' have to be passed to 'calc_near_field'.
  !>
  !> This allows to pass arbitrary (intra-cell) potentials
  !> Classes that inherit from/extend this type have to override 'get_pot' using a subroutine
  !> with an interface as specified by 'potential_func'
  type, abstract :: Potential
    contains
    procedure (potential_func), deferred :: get_pot
  end type

  !> The methods 'get_pot' of any type derived from 'Potential' have to conform to this interface
  abstract interface
    subroutine potential_func(self, v_intra, radius)
    import Potential
    class (Potential), intent(inout) :: self
    double precision, intent(out) :: v_intra(:)
    double precision, intent(in) :: radius
    end subroutine
  end interface

  !> Some test potentials - educational and testing purposes
  type, extends(Potential) :: TestPotentialConstMulti
    contains
    procedure :: get_pot => get_const_multipole
  end type
  
  !> Potential of a monopole (= point charge)
  type, extends(Potential) :: TestPotentialMonopole
    contains
    procedure :: get_pot => get_const_monopole
  end type  

  contains

  !----------------------------------------------------------------------------
  !> V_(lm)(r) = ac_wrong(lm) * (-r)**l
  subroutine calc_wrong_contribution_coeff(ac_wrong, dist_vec, charge_mom_total, gaunt)
    use MadelungCalculator_mod, only: MadelungClebschData, calc_dfac
    double precision, intent(inout) :: ac_wrong(:)
    double precision, intent(in)    :: dist_vec(3)

    double precision, intent(in)    :: charge_mom_total(:)
    type (MadelungClebschData), intent(in) :: gaunt

    integer :: lmx, L, m, ii, lmmaxd, L1, L2, lm1, lm2, lm3
    double precision :: r, rfac
    double precision :: pi
    double precision, allocatable :: ylm(:)
    double precision, allocatable :: smat(:)
    double precision, allocatable :: avmad(:,:)
    double precision, allocatable :: dfac(:,:)

    integer :: lmx_prime, lmmaxd_prime

    pi = 4.0d0*atan(1.0d0)

    lmmaxd = size(ac_wrong)
    lmx = int(sqrt(dble(lmmaxd) + 0.1) - 1)

    ! the inner summation should run up to 2*lmx
    ! this follows from Gaunt properties
    lmx_prime = 2*lmx
    lmmaxd_prime = (2*lmx + 1)**2

    allocate(ylm(lmmaxd_prime))
    allocate(smat(lmmaxd_prime))
    allocate(avmad(lmmaxd, lmmaxd))
    allocate(dfac(0:lmx, 0:lmx))

    ! calculate complicated prefactor
    call calc_dfac(dfac, lmx)

    call ymy(dist_vec(1),dist_vec(2),dist_vec(3),r,ylm,lmx_prime)

    do L = 0,lmx_prime
       rfac = (1. / (r**(L+1)))
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

      ! --> Gaunt property: l1+l2<=l3

      if (lm1 <= lmmaxd .and. lm2 <= lmmaxd .and. lm3 <= lmmaxd_prime) then
          AVMAD(LM1,LM2) = AVMAD(LM1,LM2) + &
                           2.0D0*DFAC(L1,L2)*SMAT(LM3)*gaunt%CLEB(ii)  ! factor 2 comes from the use of Rydberg units
      end if
    end do

    ac_wrong = matmul(avmad, charge_mom_total)

  end subroutine calc_wrong_contribution_coeff

  !----------------------------------------------------------------------------
  !> this is a general routine for shifting sph. harm. expansions
  subroutine calc_near_field(v_near, radius, dist_vec, pot, lmax_prime)
    double precision, intent(out) :: v_near(:)  !indices (lm)
    double precision, intent(in) :: radius
    double precision, intent(in) :: dist_vec(3)
    class(Potential), intent(inout) :: pot
    integer, intent(in), optional :: lmax_prime

    integer :: lmmaxd, lmax, lmmaxd_prime
    integer, parameter :: NUM_LEBEDEV = 434
    double precision :: v_leb(3)
    double precision :: vec(3)
    double precision :: weight_leb(NUM_LEBEDEV)
    double precision, allocatable :: sph_harm_leb(:,:)
    double precision, allocatable :: sph_harm(:)
    double precision, allocatable :: v_intra(:)
    double precision, allocatable :: integrand(:,:)
    double precision, allocatable :: temp(:)
    double precision :: norm_vec
    double precision :: dummy
    double precision :: FOUR_PI ! prefactor for Lebedev: 4*pi

    integer :: ij, lm
    integer :: lmax_p


    FOUR_PI = 16.0d0 * atan(1.0d0)

    lmmaxd = size(v_near)
    lmax = int(sqrt(dble(lmmaxd) + 0.1) - 1)
    
    if (.not. present(lmax_prime)) then
      lmax_p = lmax
    end if
    
    lmmaxd_prime = (lmax_p+1)**2

    CHECKASSERT( (lmax + 1)**2 == lmmaxd )

    allocate(sph_harm(lmmaxd_prime))
    allocate(v_intra(lmmaxd_prime))
    allocate(sph_harm_leb(NUM_LEBEDEV, lmmaxd))
    allocate(integrand(NUM_LEBEDEV, lmmaxd))
    allocate(temp(lmmaxd))

    do ij = 1, NUM_LEBEDEV

      call LEBEDEV(ij, v_leb(1), v_leb(2), v_leb(3), weight_leb(ij))

      call ymy(v_leb(1), v_leb(2), v_leb(3), dummy, temp, lmax)
      do lm = 1, lmmaxd
        sph_harm_leb(ij, lm) = temp(lm)
      end do

      vec = radius * v_leb + dist_vec

      if (vec(1)**2 + vec(2)**2 + vec(3)**2 > 0.0d0) then  ! can be zero
        call ymy(vec(1), vec(2), vec(3), norm_vec, sph_harm, lmax_p)
      else
        norm_vec = 0.0d0
        sph_harm = 0.0d0
        sph_harm(1) = 1.d0 / sqrt(FOUR_PI)
      end if

      ! get intracell potential at radius 'norm_vec'
      call pot%get_pot(v_intra, norm_vec)

      ! perform summation over L'
      integrand(ij,1) = dot_product(sph_harm, v_intra)
    end do

    integrand(:,1) = integrand(:,1) * weight_leb * FOUR_PI

    do lm = 1, lmmaxd
      integrand(:,lm) = integrand(:,1)
    end do

    integrand = integrand * sph_harm_leb

    v_near = sum(integrand, 1)

  end subroutine calc_near_field

  !----------------------------------------------------------------------------
  ! Tests orthgonality of real spherical harmonics up to LMAX
  ! The first column of 'integrand' must be 1.0 and the others = 0.0
  subroutine test_lebedev()
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
  end subroutine test_lebedev

  !----------------------------------------------------------------------------
  ! evaluate spherical harmonic expansion at angles given by 'vec'.
  double precision function eval_expansion(coeffs, vec)
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
  end function eval_expansion

!------------------------------------------------------------------------------
! Some test potentials follow
!------------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  !> A test potential: constant multipole moments
  subroutine get_const_multipole(self, v_intra, radius)
    ! Test potential: assume multipoles Q_L = 1.0d0
    class (TestPotentialConstMulti), intent(inout) :: self
    double precision, intent(out) :: v_intra(:)
    double precision, intent(in) :: radius

    integer :: lm, L, M

    lm = 1
    L = 0
    M = 0
    do while (.true.)
      v_intra(lm) = 4.0d0 * sqrt(4.0d0 * atan(1.0d0)) / (radius**(L+1)) / (2*L + 1)
      lm = lm + 1
      M = M + 1
      if (M > L) then
        L = L + 1
        M = -L
      end if
      if (lm > size(v_intra)) exit
    end do

  end subroutine get_const_multipole

  !----------------------------------------------------------------------------
  !> A test potential: potential of a (unit) monopole
  subroutine get_const_monopole(self, v_intra, radius)
    class (TestPotentialMonopole), intent(inout) :: self
    double precision, intent(out) :: v_intra(:)
    double precision, intent(in) :: radius

    integer :: lm, L, M

    v_intra = 0.0d0
    v_intra(1) = 4.0d0 * sqrt(4.0d0 * atan(1.0d0)) / radius

  end subroutine get_const_monopole

end module NearField_mod

!  program test_it
!    use NearField_mod, only: ...
!    use NearField_kkr_mod, only: ...
!    implicit none
!    integer, parameter :: LMAX = 4
!    integer, parameter :: LMMAXD = (LMAX+1)**2
!    integer, parameter :: NPOINTS = 100
!    double precision v_near(LMMAXD)
!    
!    type(TestPotentialConstMulti) :: pot
!    type(IntracellPotential) :: intra_pot
!  
!    double precision :: d(3) = (/0.0d0, 0.0d0, 2.0000d0 /)
!    double precision :: vec(3)
!    double precision :: ac_wrong(LMMAXD)
!    double precision :: v_mad_wrong(LMMAXD)
!    double precision :: cmom(LMMAXD)
!  
!    double precision :: radius =  1.400000d0
!    double precision :: r_temp
!  
!    type (MadelungClebschData) :: clebsch
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
!    !cmom = 0.0d0
!    cmom = 1.0d0 / sqrt(16.0d0 * atan(1.0d0)) 
!    !cmom(3) = 1.0d0 / sqrt(16.0d0 * atan(1.0d0)) 
!  
!    call calc_wrong_contribution_coeff(ac_wrong, d, cmom, clebsch)
!  
!    ii = 1
!    do L = 0, LMAX
!      do M = -L, L
!        v_mad_wrong(ii) = (ac_wrong(ii) * (-radius)**L)
!        write(*,*) L, M, v_near(ii), v_mad_wrong(ii) 
!        ii = ii + 1
!      end do
!    end do
!  
!    vec = 0.0d0
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
!      intra_pot%radial_points(ii) = (ii) * 1.0d0 / (NPOINTS)
!      call pot%get_pot(v_near, intra_pot%radial_points(ii))
!      intra_pot%v_intra_values(ii, :) = v_near
!    end do
!    
!    call intra_pot%init()
! 
!    call calc_near_field(v_near, radius, d , intra_pot, LMAX)
!    write(*,*) v_near
!    write(*,*) eval_expansion(v_near, vec)
! 
! !    do ii = 1, NPOINTS
! !      call pot%get_pot(v_near, intra_pot%radial_points(ii))
! !      write(*,*) v_near(1)
! !      call intra_pot%get_pot(v_near, intra_pot%radial_points(ii))
! !      write(*,*) v_near(1)
! !    end do
!    
!    !do ii = 1, NPOINTS * 4
!    !  r_temp = (ii) * 2.0d0 / (NPOINTS*4) + 1.0d0
!    !  call intra_pot%get_pot(v_near, r_temp)
!    !  write(*,*) r_temp, v_near(7), v_near(15)
!    !end do
!    !call calc_near_field(v_near, radius, d , intra_pot, LMAX)
!    !write(*,*) v_near
!  end program
