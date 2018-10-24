!-------------------------------------------------------------------------------
!> Summary: Check if a string is contained in the `opt` array
!> Author: 
!> Category: undefined, KKRhost
!> Deprecated: False 
!> Check if a string is contained in the `opt` array
!-------------------------------------------------------------------------------
!> @note Jonathan Chico: This function is not contained in a module, it should be 
!> included in a tool like module where auxiliary functions and subrotuines are 
!> located.
!> @endnote
!-------------------------------------------------------------------------------
logical function opt(string)
  ! ***********************************************************************

  ! OPT = 'STRING  ' IS CONTAINED IN /OPTC/.

  ! ------------------------------------------------------------------------
  use :: mod_wunfiles, only: t_params

  implicit none
  character (len=8), intent (in) :: string
  integer :: i
  character (len=8), dimension(32) :: optc

  optc = t_params%optc

  opt = .false.
  do i = 1, 32
    if (string==optc(i)) opt = .true.
  end do
  return
end function opt
