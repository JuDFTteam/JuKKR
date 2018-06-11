subroutine writehoststructure(bravais, nrbasis, rbasis, naezd, nembd)
  use mod_version_info
  use mod_md5sums
  use mod_DataTypes
      Use mod_datatypes, Only: dp
  implicit none
!interface
  real (kind=dp), intent(in) :: bravais(3, 3)
  integer, intent(in) :: nrbasis
  integer, intent(in) :: naezd
  integer, intent(in) :: nembd
  real (kind=dp), intent(in) :: rbasis(3, naezd+nembd)

!local
  integer :: iatom
  real (kind=dp) :: wght
  character (len=256) :: uio
  integer :: ier
  integer :: itemp1(12), it

  open (unit=3463453, file='kkrflex_hoststructure.dat')
  call version_print_header(3463453, '; '//md5sum_potential//'; '// &
    md5sum_shapefun)

  write (3463453, '(100A)') '[bravais]'
  write (3463453, '(100A)') '#   x   y   z - component'
  do iatom = 1, 3
    write (3463453, '(3F21.16)') bravais(:, iatom)
  end do !nbasis
  write (3463453, '(100A)') '[basis]'
  write (3463453, *) nrbasis
  write (3463453, '(100A)') '#   x   y   z    weight'
  do iatom = 1, nrbasis
    call ioinput('ATOMINFO        ', uio, iatom+1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*)(itemp1(it), it=1, 12), wght
    else
      wght = 1.d0
    end if
    write (3463453, '(4F21.16)') rbasis(:, iatom), wght
  end do !nbasis

  close (3463453)



end subroutine
