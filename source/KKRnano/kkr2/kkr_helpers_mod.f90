module kkr_helpers_mod
  implicit none
  private
  public :: lmmaxToLmax
  
  contains

  !> Returns the l_max value which corresponds to a given lmmax-value
  !! for non-relativistic case.
  !! Terminates program for illegal value.
  integer function lmmaxToLmax(lmmax)
    implicit none
    integer, intent(in) :: lmmax
    integer :: ii
    integer :: lm

    ii = 0

    do while (.true.)
      lm = (ii + 1)**2
      if (lm == lmmax) then
        exit
      end if
      if (lm > lmmax) then
        write(*,*) "Illegal value passed to lmmaxToLmax. ", lmmax
        stop
      end if
      ii = ii + 1
    end do

    lmmaxToLmax = ii
  end function
end module kkr_helpers_mod
