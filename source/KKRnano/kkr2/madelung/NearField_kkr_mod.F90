module NearField_kkr_mod
  use NearField_mod, only: Potential
  use, intrinsic :: ieee_arithmetic, only: ieee_value, IEEE_SIGNALING_NAN
  
  integer, parameter, private :: PI = 3.1415926535897932d0
  
  type, extends(Potential) :: IntracellPotential

    double precision, allocatable :: charge_moments(:)
    double precision, allocatable :: radial_points(:)
    double precision, allocatable :: v_intra_values(:, :)
    
    ! The following objects are private - do not modify
    double precision, allocatable :: xarray(:)
    double precision, allocatable :: yarray(:,:)
    double precision, allocatable :: y2ndder(:,:)
    double precision :: rws_bounding
    integer :: counter_spline
    
    contains
    
    procedure create => createIntracellPot
    procedure destroy => destroyIntracellPotential
    procedure get_pot => get_intracell
  end type
  
  contains
  
  !----------------------------------------------------------------------------
  subroutine get_intracell(self, v_intra, radius)
    implicit none
    class (IntracellPotential), intent(inout) :: self
    double precision, intent(out) :: v_intra(:)
    double precision, intent(in) :: radius
    
    double precision :: yderiv
    integer :: lmpotd
    integer :: lm, L, M
    
    lmpotd = size(v_intra)
    
    L = 0
    M = 0
    
    if (radius > self%rws_bounding) then
    
      do lm = 1, lmpotd

        ! calculate potential from charge moments
        v_intra(lm) = 8 * PI / (2*dble(L) + 1) * self%charge_moments(lm) / (radius**(L+1))

        M = M + 1
        if (M > L) then
          L = L + 1
          M = -L
        end if
      end do

    else

      ! calculate potential from known values within the bounding sphere
      do lm = 1, lmpotd
        call splint(self%xarray, self%yarray(:, lm), self%y2ndder(:, lm), &
                    self%counter_spline, radius, v_intra(lm), yderiv)
      end do

    endif

  end subroutine
  
  !----------------------------------------------------------------------------
  subroutine createIntracellPot(self, lmpotd, irmd)
    implicit none
    class (IntracellPotential), intent(inout) :: self
    integer, intent(in) :: lmpotd
    integer, intent(in) :: irmd
    double precision :: nan
    
    nan = ieee_value(nan, IEEE_SIGNALING_NAN)
    
    allocate(self%charge_moments(lmpotd))
    allocate(self%radial_points(irmd))
    allocate(self%v_intra_values(irmd, lmpotd))
    allocate(self%xarray(irmd))
    allocate(self%yarray(irmd, lmpotd))
    allocate(self%y2ndder(irmd, lmpotd))
    
    self%charge_moments = nan
    self%radial_points = nan
    self%v_intra_values = nan
    self%xarray = nan
    self%yarray = nan
    self%y2ndder = nan
  end subroutine
  
    !----------------------------------------------------------------------------
  subroutine createIntracellPotential(self, lmpotd, irmd)
    implicit none
    class (IntracellPotential), intent(inout) :: self
    integer, intent(in) :: lmpotd
    integer, intent(in) :: irmd
    allocate(self%charge_moments(lmpotd))
    allocate(self%radial_points(irmd))
    allocate(self%v_intra_values(irmd, lmpotd))
    
    ! allocate arrays for spline interpolation
    allocate(self%xarray(irmd))
    allocate(self%yarray(irmd, lmpotd))
    allocate(self%y2ndder(irmd, lmpotd))
  end subroutine
  
  !----------------------------------------------------------------------------
  subroutine initIntracellPotential(self)
    implicit none
    class (IntracellPotential), intent(inout) :: self
    
    integer :: lmpotd
    integer :: irmd
    integer :: ii
    integer :: L, M
    integer :: counter
    double precision :: derivative
    double precision :: pi
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
        self%yarray(counter, :) = self%v_intra_values(ii, :)
      end if
    end do
          
    self%counter_spline = counter
    
    ! TODO: remove duplicate x-values
    L = 0
    M = 0
    do ii = 1, lmpotd

      ! match derivative of spline interpolation at bounding sphere
      derivative = -8 * PI / (2 * dble(L) + 1) * (dble(L) + 1) & 
                   * self%charge_moments(ii) / self%rws_bounding**(L+2)

      call spline(irmd, self%xarray, self%yarray(:, ii), counter, 1.d35, derivative, self%y2ndder(:, ii))
      
      M = M + 1
      if (M > L) then
        L = L + 1
        M = -L
      end if
    end do
    
  end subroutine
  
  !----------------------------------------------------------------------------
  subroutine destroyIntracellPotential(self)
    implicit none
    class (IntracellPotential), intent(inout) :: self
    deallocate(self%charge_moments)
    deallocate(self%radial_points)
    deallocate(self%v_intra_values)
    deallocate(self%xarray)
    deallocate(self%yarray)
    deallocate(self%y2ndder)
  end subroutine
end module