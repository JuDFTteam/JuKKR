!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Writes out the header of the t-matrices decimation file
!> Author: 
!> Writes out the header of the t-matrices decimation file
!------------------------------------------------------------------------------------
module mod_outtmathost
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Writes out the header of the t-matrices decimation file
  !> Author: 
  !> Category: input-output, single-site, KKRhost
  !> Deprecated: False
  !> Writes out the header of the t-matrices decimation file
  !-------------------------------------------------------------------------------
  subroutine outtmathost(alat,ins,krel,kmrot,nspin,naez,lmmax0d,bravais,rbasis,qmtet, &
    qmphi,e2in,tk,npol,npnt1,npnt2,npnt3)

    use :: mod_version_info, only: version_print_header
    use :: mod_runoptions, only: disable_print_serialnumber
    implicit none
    ! ..
    ! .. Input variables
    integer, intent(in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent(in) :: krel   !! Switch for non- (or scalar-) relativistic/relativistic (Dirac) program (0/1). Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent(in) :: naez   !! Number of atoms in unit cell
    integer, intent(in) :: npol   !! Number of Matsubara Poles (EMESHT)
    integer, intent(in) :: kmrot  !! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
    integer, intent(in) :: nspin  !! Counter for spin directions
    integer, intent(in) :: lmmax0d !! Maximum l component in wave function expansion
    integer, intent(in) :: npnt1  !! number of E points (EMESHT) for the contour integration
    integer, intent(in) :: npnt2  !! number of E points (EMESHT) for the contour integration
    integer, intent(in) :: npnt3  !! number of E points (EMESHT) for the contour integration
    real (kind=dp), intent(in) :: tk    !! Temperature
    real (kind=dp), intent(in) :: alat  !! Lattice constant in a.u.
    real (kind=dp), intent(in) :: e2in
    real (kind=dp), dimension(*), intent(in) :: qmtet !! $$ \theta $$ angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension(*), intent(in) :: qmphi !! $$ \phi $$ angle of the agnetization with respect to the z-axis
    real (kind=dp), dimension(3,*), intent(in) :: rbasis  !! Position of atoms in the unit cell in units of bravais vectors
    real (kind=dp), dimension(3,3), intent(in) :: bravais !! Bravais lattice vectors
    ! .. Local variables
    integer :: i, ih
    ! ----------------------------------------------------------------------
    write (1337, '(5X,A,/)') '< DECIOPT > : writing header of decimation file'

    open (37, file='decifile', status='unknown')
    call version_print_header(37, disable_print=disable_print_serialnumber)
    write (37, fmt=*) 'INVERSE T-MATRIX AND CMOMS'
    write (37, fmt=100)
    write (37, fmt=110) alat, nspin, naez, lmmax0d, ins, krel, kmrot
    write (37, fmt=120) bravais
    if (krel==0) then
      write (37, fmt=130)
      do ih = 1, naez
        write (37, fmt=150)(rbasis(i,ih), i=1, 3)
      end do
    else
      write (37, fmt=140)
      do ih = 1, naez
        write (37, fmt=160)(rbasis(i,ih), i=1, 3), qmtet(ih), qmphi(ih)
      end do
    end if
    write (37, fmt=170) e2in, tk
    write (37, fmt=180) npnt1, npnt2, npnt3, npol
    close (37)
    ! ----------------------------------------------------------------------
100 format (' Vectors in lattice constant units', /, '                                 ')
110 format ('ALAT=', f9.6, ' NSPIN=', i2, '  NAEZ=', i3, ' LMMAX=', i3, ' INS=', i1, ' KREL=', i1, ' KMROT=', i1)
120 format ('BRAVAIS ', /, 3f8.4, /, 3f8.4, /, 3f8.4)
130 format ('RBASIS')
140 format ('RBASIS', 20x, 'MAGNETISATION ANGLES THETA/PHI')
150 format (3f8.4)
160 format (3f8.4, 2f9.4)
170 format ('EF=', f10.6, ' TEMP=', f10.4, ' Kelvin')
180 format ('N1=', i3, ' N2=', i3, ' N3=', i3, ' NPOL=', i3)
  end subroutine outtmathost

end module mod_outtmathost
