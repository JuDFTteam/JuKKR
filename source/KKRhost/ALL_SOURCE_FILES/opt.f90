module mod_opt

contains

! ***********************************************************************
logical function opt(string)
  ! ***********************************************************************

  ! OPT = 'STRING  ' IS CONTAINED IN /OPTC/.

  ! ------------------------------------------------------------------------
  use :: mod_wunfiles, only: t_params

  implicit none
  character (len=8), intent (in) :: string
  integer :: i
  character (len=8) :: optc(32)

  optc = t_params%optc

  opt = .false.
  do i = 1, 32
    if (string==optc(i)) opt = .true.
  end do
  return
end function opt

end module mod_opt
