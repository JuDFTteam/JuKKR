!> Near-field electrostatics: Interpolation of intracell potential.
!>
!> Routines to calculate intracell potential values for radii
!> not specified in the radial mesh.
!>
!> When radius is outside of Voronoi cell then use multipole moments
!> to calculate potential value (exactly).
!> When radius is inside Voronoi cell, use cubic splines to interpolate
!> from values on radial mesh.
!>
!> @author Elias Rabel
module NearField_kkr_mod
! use NearField_mod, only: Potential
  implicit none
  private
  public :: IntracellPotential, create, destroy, init, get
  
  !----------------------------------------------------------------------------
  !> Usage:
  !> type(IntracellPotential) :: pot
  !> call create(pot, lmpot, irmd)
  !> pot%charge_moments = ... ! set charge moments
  !> pot%radial_points = ...
  !> pot%v_intra_values = ...
  !> call init(pot)
  !> ...  ! now use it
  !> call destroy(pot)
  type :: IntracellPotential
    double precision, allocatable :: charge_moments(:)
    double precision, allocatable :: radial_points(:)
    double precision, allocatable :: v_intra_values(:, :)
    
    ! The following objects are private - do not modify
    double precision, allocatable :: xarray(:)
    double precision, allocatable :: yarray(:,:)
    double precision, allocatable :: y2ndder(:,:)
    double precision :: rws_bounding
    integer :: counter_spline
  endtype

  interface create
    module procedure createIntracellPotential
  endinterface

  interface get
    module procedure getIntracellPotential
  endinterface

  interface init
    module procedure initIntracellPotential
  endinterface
  
  interface destroy
    module procedure destroyIntracellPotential
  endinterface
  
  contains
  
  !----------------------------------------------------------------------------
  subroutine getIntracellPotential(self, v_intra, radius)
    use Constants_mod, only: pi
    type(IntracellPotential), intent(inout) :: self
    double precision, intent(out) :: v_intra(:)
    double precision, intent(in) :: radius
    
    double precision :: yderiv
    integer :: lmpotd, lm, ell, emm
    
    lmpotd = size(v_intra)
    
    ell = 0
    emm = 0
    
    if (radius > self%rws_bounding) then
    
      do lm = 1, lmpotd

        ! calculate potential from charge moments
        v_intra(lm) = 8.0d0 * pi * self%charge_moments(lm) / (radius**(ell+1) * (2*ell + 1))

        emm = emm + 1
        if (emm > ell) then
          ell = ell + 1
          emm = -ell
        endif
      enddo ! lm

    else

      ! calculate potential from known values within the bounding sphere
      do lm = 1, lmpotd
        call splint(self%xarray, self%yarray(:,lm), self%y2ndder(:,lm), self%counter_spline, radius, v_intra(lm), yderiv)
      enddo ! lm

    endif

  endsubroutine ! get
  
  !----------------------------------------------------------------------------
  subroutine createIntracellPotential(self, lmpotd, irmd)
#ifndef __GFORTRAN__   
  use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_SIGNALING_NAN
#endif
    type(IntracellPotential), intent(inout) :: self
    integer, intent(in) :: lmpotd
    integer, intent(in) :: irmd
#ifndef __GFORTRAN__   
    double precision :: nan
    nan = ieee_value(nan, IEEE_SIGNALING_NAN)
#else
    double precision, parameter :: nan = -999.9d9
#endif

    allocate(self%charge_moments(lmpotd))
    allocate(self%radial_points(irmd))
    allocate(self%v_intra_values(irmd,lmpotd))
    allocate(self%xarray(irmd))
    allocate(self%yarray(irmd,lmpotd))
    allocate(self%y2ndder(irmd,lmpotd))
    
    self%charge_moments = nan
    self%radial_points = nan
    self%v_intra_values = nan
    self%xarray = nan
    self%yarray = nan
    self%y2ndder = nan
  endsubroutine ! create
  
  !----------------------------------------------------------------------------
  subroutine initIntracellPotential(self)
    use Constants_mod, only: pi
    type(IntracellPotential), intent(inout) :: self
    
    integer :: lmpotd, irmd, ii, ell, emm, counter
    double precision :: derivative
    double precision, parameter :: TOL = 1.d-10
    
    lmpotd = size(self%v_intra_values, 2)
    irmd = size(self%v_intra_values, 1)
    
    self%rws_bounding = self%radial_points(size(self%radial_points))
    
    counter = 1
    self%xarray(1) = self%radial_points(1)
    self%yarray(1, :) = self%v_intra_values(1, :)
    
    ! remove duplicate x - values (panels!)
    ! assumption: no discontinuities (jumps)!
    do ii = 2, irmd
      if (abs(self%radial_points(ii) - self%radial_points(ii-1)) > TOL) then
        counter = counter + 1
        self%xarray(counter) = self%radial_points(ii)
        self%yarray(counter,:) = self%v_intra_values(ii,:)
      endif
    enddo ! ii
          
    self%counter_spline = counter
    
    ell = 0
    emm = 0
    do ii = 1, lmpotd

      ! match derivative of spline interpolation at bounding sphere
      derivative = -8.d0 * pi * (ell + 1.d0) * self%charge_moments(ii) / ( self%rws_bounding**(ell+2) * (2*ell + 1) )

      call spline(irmd, self%xarray, self%yarray(:,ii), counter, 1.d35, derivative, self%y2ndder(:,ii))
      
      emm = emm + 1
      if (emm > ell) then
        ell = ell + 1
        emm = -ell
      endif
    enddo ! ii
    
  endsubroutine ! init
  
  !----------------------------------------------------------------------------
  elemental subroutine destroyIntracellPotential(self)
    type(IntracellPotential), intent(inout) :: self
    integer :: ist ! ignore status
    deallocate(self%charge_moments, self%radial_points, &
               self%v_intra_values, self%xarray, &
               self%yarray, self%y2ndder, stat=ist)
  endsubroutine ! destroy

endmodule ! NearField_kkr_mod
