module mod_sname

contains

! ************************************************************************
subroutine sname(name, new, band)
  ! ************************************************************************
  use :: mod_datatypes, only: dp
   use mod_length
  implicit none
  ! .. scalar arguments
  integer :: band
  character (len=40) :: name, new

  ! .. locals
  integer :: i, l, lo
  character (len=1) :: ch(50), poi
  character (len=10) :: s
  ! ------------------------------------------------------------------------
  poi = '.'
  if (band<0) then
    lo = nint(log(real(-band,kind=dp))/log(10.0e0_dp)) + 1
  else if (band==0) then
    lo = 0
  else
    lo = nint(log(real(band,kind=dp))/log(10.0e0_dp))
  end if

  ! write(6,*) 'LO ',lo

  read (name, fmt='(255a1)')(ch(i), i=1, 40)
  l = length(ch, 40)
  ! write(6,*) 'L  ',l

  ! write(6,*) 'CH ',(CH(I),I=1,25)

  write (s, fmt='(I10)') band
  ! write(6,*) 'S  ',s

  read (s, fmt='(255A1)')(ch(i), i=l+1, l+10)
  ! write(6,*) 'CH ',(CH(I),I=L+1,L+10)

  write (new, fmt='(255A1)')(ch(i), i=1, l), poi, (ch(i), i=l+10-lo, l+10)

  return
end subroutine sname

end module mod_sname
