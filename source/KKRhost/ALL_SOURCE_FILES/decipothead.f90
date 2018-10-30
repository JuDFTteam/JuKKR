!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_decipothead

  private
  public :: decipothead

contains

  !-------------------------------------------------------------------------------
  !> Summary: Read header of decimate potential file
  !> Author: 
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Reads in the header of the host potential-file 'decimate.pot'
  !> checking for consistencies with the actual 2D system
  !>                                         V. Popescu - Munich, Dec 04
  !-------------------------------------------------------------------------------
  subroutine decipothead(ihost, filehost, ilhost, nathost, vacflag, alat, bravsys, nq, nt, bravais, efermi, insh, krelh, nspinh, ins, krel, nspin, kmrot)

    use :: mod_version_info, only: version_check_header
    use :: mod_datatypes, only: dp
    implicit none
    ! ..
    ! .. Arguments
    real (kind=dp) :: alat, efermi
    character (len=40) :: filehost
    integer :: ihost, ilhost, ins, insh, kmrot, krel, krelh, nathost, nq, nspin, nspinh, nt
    real (kind=dp) :: bravais(3, 3), bravsys(3, 3)
    logical :: vacflag(2)
    ! ..
    ! .. Locals
    real (kind=dp) :: alath
    integer :: i, ih, ios, kmroth
    ! ----------------------------------------------------------------------
    if (.not. (filehost(1:7)=='vacuum')) then
      open (36+ihost, file=filehost, status='OLD', iostat=ios)
      call version_check_header(36+ihost)
      ! ......................................................................
      if (ios>0) then
        write (6, '(/,5X,2A)') 'ERROR: Can not open host file ', filehost(1:ilhost)
        write (6, '(12X,A,A)') 'vacuum    entry should be used to', ' set a vacuum-host on this side'
        stop
      end if

      do ih = 1, 2
        read (36+ihost, *)
      end do
      read (36+ihost, 110) krelh, insh, nspinh, kmroth
      read (36+ihost, '(2(7X,I3),7X,F12.8)') nq, nt, alath

      ! --> non-spherical and/or CPA cases not implemented yet

      if (insh/=0) stop ' INS<>0 not implemented '
      if (nq/=nt) stop ' CPA-host not implemented'
      ! ......................................................................
      if ((nq/=nathost) .or. (kmroth/=kmrot) .or. (abs(alath-alat)>1e-6_dp)) then
        write (6, 150) filehost(1:ilhost)
        write (6, 120) '  NAEZ KMROT ALAT '
        write (6, 130) 'syst: ', nathost, kmrot, alat
        write (6, 130) 'host: ', nq, kmroth, alath
        write (6, *)
        stop
      end if
      if ((insh/=ins) .or. (nspinh/=nspin) .or. (krelh/=krel)) then
        write (6, 150) filehost(1:ilhost)
        write (6, 120) '  KREL   INS NSPIN'
        write (6, 140) 'syst: ', krel, ins, nspin
        write (6, 140) 'host: ', krelh, insh, nspinh
        write (6, *)
        stop
      end if
      ! ......................................................................
      read (36+ihost, '(10X,F10.6)') efermi
      read (36+ihost, *)
      do ih = 1, 3
        read (36+ihost, 100)(bravais(i,ih), i=1, 3)
      end do
      alath = 0e0_dp
      do ih = 1, 2
        do i = 1, 2
          alath = alath + bravais(i, ih) - bravsys(i, ih)
        end do
      end do
      if (abs(alath)>1e-6_dp) then
        write (6, 150) filehost(1:ilhost)
        write (6, 160)
        do ih = 1, 2
          write (6, 170) ih, (bravais(i,ih), i=1, 2), (bravsys(i,ih), i=1, 2)
        end do
        write (6, *)
        stop
      end if
      ! ......................................................................
      read (36+ihost, *)
    else
      vacflag(ihost) = .true.
    end if
    ! ----------------------------------------------------------------------
100 format (3f8.4)
110 format (8(7x,i3))
120 format (14x, a, /, 8x, 28('-'))
130 format (8x, a6, 2i6, f10.6)
140 format (8x, a6, 3i6)
150 format (6x, 'ERROR: host-potential ( ', a, ' )', /, 13x, 'not compatible with your input/sytem')
160 format (14x, '2D lattice', 4x, 'host', 14x, 'system', /, 14x, 38('-'))
170 format (10x, 'a_', i1, 1x, 2f9.5, 2x, 2f9.5)
  end subroutine decipothead

end module mod_decipothead
