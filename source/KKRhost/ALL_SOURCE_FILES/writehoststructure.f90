!------------------------------------------------------------------------------------
!> Summary: Printing to file the `kkrflex_hoststructure` file
!> Author:
!> Printing to file the `kkrflex_hoststructure` file containing geometrical information
!> of the host structure to be used by the KKRimp program
!------------------------------------------------------------------------------------
module mod_writehoststructure

contains
  !-------------------------------------------------------------------------------
  !> Summary: Printing to file the `kkrflex_hoststructure` file
  !> Author: 
  !> Category: input-output, KKRhost
  !> Deprecated: False 
  !> Printing to file the `kkrflex_hoststructure` file containing geometrical information
  !> of the host structure to be used by the KKRimp program
  !-------------------------------------------------------------------------------
  subroutine writehoststructure(bravais, nrbasis, rbasis, naezd, nembd)
    use :: mod_version_info
    use :: mod_md5sums
    use :: mod_datatypes, only: dp
    use :: mod_ioinput
    implicit none
    ! .. Input variables
    integer, intent (in) :: nrbasis !! Number of basis layers of right host (repeated units)
    integer, intent (in) :: naezd !! Number of atoms in unit cell
    integer, intent (in) :: nembd !! Number of 'embedding' positions
    real (kind=dp), dimension(3,3), intent (in) :: bravais !! Bravais lattice vectors
    real (kind=dp), dimension(3, naezd+nembd), intent (in) :: rbasis !! Position of atoms in the unit cell in units of bravais vectors

    ! .. Local variables
    integer :: iatom
    real (kind=dp) :: wght
    character (len=256) :: uio
    integer :: ier, it
    integer, dimension(8) :: itemp1
    real (kind=dp), dimension(4) :: ftemp

    open (unit=3463453, file='kkrflex_hoststructure.dat')
    call version_print_header(3463453, '; '//md5sum_potential//'; '//md5sum_shapefun)

    write (3463453, '(100A)') '[bravais]'
    write (3463453, '(100A)') '#   x   y   z - component'
    do iatom = 1, 3
      write (3463453, '(3F21.16)') bravais(:, iatom)
    end do                         ! nbasis
    write (3463453, '(100A)') '[basis]'
    write (3463453, *) nrbasis
    write (3463453, '(100A)') '#   x   y   z    weight'
    do iatom = 1, nrbasis
      call ioinput('ATOMINFO        ', uio, iatom+1, 7, ier)
      if (ier==0) then
        ! read weight from last column in atominfo (ftemp and itemp1 are used to
        ! skip the first 12 entries in that line)
        read (unit=uio, fmt=*) ftemp(1), (itemp1(it), it=1, 8), ftemp(2:4), wght
      else
        wght = 1.d0
      end if
      write (3463453, '(4F21.16)') rbasis(:, iatom), wght
    end do                         ! nbasis

    close (3463453)

  end subroutine writehoststructure

end module mod_writehoststructure
