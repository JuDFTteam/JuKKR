!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Printing to file the TBKKR files, containing key information of the input parameters
!> Author:
!> Printing to file the TBKKR files, containing key information of the input
!> parameters
!------------------------------------------------------------------------------------
module mod_write_tbkkr_files
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Printing to file the TBKKR files, containing key information of the input parameters
  !> Author: 
  !> Category: input-output, KKRhost
  !> Deprecated: False 
  !> Printing to file the TBKKR files, containing key information of the input
  !> parameters
  !-------------------------------------------------------------------------------
  subroutine write_tbkkr_files(lmax,nemb,ncls,natyp,naez,ielast,ins,alat,bravais,   &
    recbv,rbasis,cls,nacls,rcls,ezoa,atom,rr,nspin,nr,korbit,nclsd,naclsd)

    use :: mod_version_info, only: serialnr

    implicit none
    ! interface
    integer, intent (in) :: nr     !! Number of real space vectors rr
    integer, intent (in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: nemb   !! Number of 'embedding' positions
    integer, intent (in) :: ncls   !! Number of reference clusters
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: nclsd  !! Maximum number of different TB-clusters
    integer, intent (in) :: naclsd !! Maximum number of atoms in a TB-cluster
    integer, intent (in) :: ielast
    integer, intent (in) :: korbit !! Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
    real (kind=dp), intent (in) :: alat !! Lattice constant in a.u.
    real (kind=dp), dimension (3, 3), intent (in) :: recbv !! Reciprocal basis vectors
    real (kind=dp), dimension (3, 3), intent (in) :: bravais !! Bravais lattice vectors
    real (kind=dp), dimension (3, naez+nemb), intent (in) :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
    real (kind=dp), dimension (3, naclsd, nclsd), intent (in) :: rcls !! Real space position of atom in cluster
    real (kind=dp), dimension (3, 0:nr), intent (in) :: rr !! Set of real space vectors (in a.u.)
    integer, dimension (naez), intent (in) :: cls !! Cluster around atomic sites
    integer, dimension (nclsd), intent (in) :: nacls !! Number of atoms in cluster
    integer, dimension (naclsd, naez), intent (in) :: ezoa !! EZ of atom at site in cluster
    integer, dimension (naclsd, naez), intent (in) :: atom !! Atom at site in cluster
    ! .. Local variables
    integer :: i1, i2, j, naclsmax

    naclsmax = 1
    do i1 = 1, ncls
      if (nacls(i1)>naclsmax) naclsmax = nacls(i1)
    end do

    open (934, file='TBkkr_params.txt', form='formatted')
    write (934, '(A,A)') '#FILEVERSION= 2' // '   # serial: ', serialnr
    ! write nspin instead of npsind for
    ! program to work with nspin==1 case
    ! ncls instead of nclsd
    ! for smaller files (see
    ! kloopz writeout)
    ! naclsmax instead of naclsd for smaller
    ! files
    write (934, '(I8,4X,A)') lmax, 'lmaxd', lmax, 'lmax', korbit, 'korbit', nspin, 'nspin, used as nspind', nr, 'nrd', nemb, 'nembd', nemb, 'nemb', ncls, &
      'ncls once again, used as nclsd', ncls, 'ncls', natyp, 'natypd', natyp, 'natyp', naez, 'naezd', naez, 'naez', naclsmax, 'naclsmax, used as naclsd', ielast, 'ielast', ins, &
      'ins'
    close (934)

    open (935, file='TBkkr_container.txt', form='formatted')
    write (935, '(A,A)') '#FILEVERSION= 2' // '   # serial: ', serialnr
    ! write out lattice information
    write (935, '(A)') 'alat:'
    write (935, '(ES25.16)') alat
    write (935, '(A)') 'bravais:'
    write (935, '(3ES25.16)')((bravais(i1,i2),i1=1,3), i2=1, 3)
    write (935, '(A)') 'recbv:'
    write (935, '(3ES25.16)')((recbv(i1,i2),i1=1,3), i2=1, 3)
    write (935, '(A)') 'RBASIS:'
    write (935, '(3ES25.16)')((rbasis(j,i1),j=1,3), i1=1, naez+nemb)

    ! write out cluster information
    write (935, '(A)') 'CLS:'
    write (935, '(1I8)')(cls(i1), i1=1, natyp)
    write (935, '(A)') 'NACLS:'
    write (935, '(1I8)')(nacls(i1), i1=1, ncls)
    write (935, '(A)') 'RCLS:'
    do i2 = 1, ncls
      do i1 = 1, naclsmax
        write (935, '(3ES25.16)') rcls(:, i1, i2)
      end do
    end do
    write (935, '(A)') 'EZOA:'
    write (935, '(1I8)')((ezoa(i1,i2),i1=1,naclsmax), i2=1, naez)
    write (935, '(A)') 'ATOM:'
    write (935, '(1I8)')((atom(i1,i2),i1=1,naclsmax), i2=1, naez)
    write (935, '(A)') 'RR:'
    do i1 = 0, nr
      write (935, '(3ES25.16)') rr(:, i1)
    end do

    close (935)

  end subroutine write_tbkkr_files

end module mod_write_tbkkr_files
