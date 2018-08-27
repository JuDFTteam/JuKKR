! *********************************************************** 17.05.91 **
logical function test(string)
  ! ***********************************************************************

  ! TEST = 'STRING  ' IS CONTAINED IN /TESTC/.

  ! ------------------------------------------------------------------------
  use :: mod_wunfiles, only: t_params

  implicit none
  character (len=8), intent (in) :: string
  integer :: i
  character (len=8) :: testc(32)

  testc = t_params%testc

  test = .false.
  do i = 1, 32
    if (string==testc(i)) test = .true.
  end do
  return
end function test
