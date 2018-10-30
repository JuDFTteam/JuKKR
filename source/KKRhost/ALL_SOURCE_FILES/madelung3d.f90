!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculation of the Madelung potential coefficients for a 3D structure 
!> Author:
!> Calculation of the Madelung potential coefficients for a 3D structure, the coefficients
!> are then stored in an unformatted file.
!------------------------------------------------------------------------------------
module mod_madelung3d

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the Madelung potential coefficients for a 3D structure 
  !> Author: 
  !> Category: electrostatics, KKRhost, geometry 
  !> Deprecated: False 
  !> Calculation of the Madelung potential coefficients for a 3D structure, the coefficients
  !> are then stored in an unformatted file.
  !-------------------------------------------------------------------------------
  !> @note All positions must be scaled with `ALAT` to get them correct
  !> The record index is simply `(IQ1-1)*NAEZ + IQ2` for record `(IQ1,IQ2)`
  !> @endnote
  !> @todo This routine uses both naez and naezd which should be the same number, 
  !> one should replace this to eliminate redundant variables.
  !> @endtodo
  !-------------------------------------------------------------------------------
  subroutine madelung3d(lpot,yrg,wg,naez,alat,volume0,bravais,recbv,rbasis,rmax,    &
    gmax,naezd,lmxspd,lassld,lpotd, lmpotd,nmaxd,ishld,nembd,wlength)

    use :: mod_datatypes, only: dp
    use :: mod_madelgaunt, only: madelgaunt
    use :: mod_madelcoef, only: madelcoef
    use :: mod_madelout, only: madel3out
    use :: mod_lattice3d, only: lattice3d
    use :: mod_strmat, only: strmat
    use :: mod_types, only: t_madel

    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer, intent(in) :: lpot     !! Maximum l component in potential expansion
    integer, intent(in) :: naez     !! Number of atoms in unit cell
    integer, intent(in) :: naezd    !! Number of atoms in unit cell
    integer, intent(in) :: nmaxd    !! Paremeters for the Ewald summations
    integer, intent(in) :: ishld    !! Paremeters for the Ewald summations
    integer, intent(in) :: lpotd    !! Maximum l component in potential expansion
    integer, intent(in) :: lassld   !! 4*lmax
    integer, intent(in) :: lmpotd   !! (lpot+1)**2
    integer, intent(in) :: lmxspd   !! (2*lpot+1)**2
    integer, intent(in) :: nembd    !! Number of 'embedding' positions
    integer, intent(in) :: wlength  !! Word length for direct access files, compiler dependent ifort/others (1/4)
    real (kind=dp), intent(in) :: alat  !! Lattice constant in a.u.
    real (kind=dp), intent(inout) :: rmax  !! Ewald summation cutoff parameter for real space summation
    real (kind=dp), intent(inout) :: gmax  !! Ewald summation cutoff parameter for reciprocal space summation
    real (kind=dp), intent(in) :: volume0
    ! ..
    ! .. Array Arguments ..
    real (kind=dp), dimension(lassld), intent(in) :: wg !! Integr. weights for Legendre polynomials
    real (kind=dp), dimension(lassld, 0:lassld, 0:lassld), intent(in) :: yrg  !! Spherical harmonics (GAUNT2)
    real (kind=dp), dimension(3,3), intent(in) :: recbv   !! Reciprocal basis vectors
    real (kind=dp), dimension(3,3), intent(in) :: bravais !! Bravais lattice vectors
    real (kind=dp), dimension(3,naezd+nembd), intent(in) :: rbasis  !! Position of atoms in the unit cell in units of bravais vectors
    ! ..
    ! .. Local Scalars ..
    integer :: iend, iprint, iq1, iq2, nclebd
    integer :: ngmax, nrmax, nshlg, nshlr
    integer :: lrecabmad, irec
    integer :: ierr
    ! ..
    ! .. Local Arrays ..
    ! .. Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
    integer, dimension(ishld)                       :: nsg, nsr
    integer, dimension(lmxspd*lmpotd, 3)            :: icleb
    real (kind=dp), dimension(lmpotd)               :: bvmad
    real (kind=dp), dimension(lmxspd*lmpotd)        :: cleb
    real (kind=dp), dimension(lmpotd, lmpotd)       :: avmad
    real (kind=dp), dimension(6,6)                  :: smat1, smat2
    real (kind=dp), dimension(3,nmaxd)              :: gn, rm
    real (kind=dp), dimension(lmxspd, naezd, naezd) :: madelsmat

    logical, external :: test
    ! ......................................................................
    iprint = 0
    nclebd = lmxspd*lmpotd

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, '(79("="))')
    write (1337, '(18X,A)') 'MADELUNG3D: setting bulk Madelung coefficients'
    write (1337, '(79("="))')
    write (1337, *)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    ! ======================================================================
    call lattice3d(alat,bravais,recbv,ngmax,nrmax,nshlg,nshlr,nsg,nsr,gn,rm,rmax,   &
      gmax,iprint,nmaxd,ishld)

    call strmat(alat,lpot,naez,ngmax,nrmax,nsg,nsr,nshlg,nshlr,gn,rm,rbasis,        &
      madelsmat,volume0,iprint,lassld,lmxspd,naezd)
    ! ======================================================================

    lrecabmad = wlength*2*lmpotd*lmpotd + wlength*2*lmpotd
    ! lrecabmad = wlength*kind(0.0_dp)*lmpotd*lmpotd + wlength*kind(0.0_dp)*lmpotd
    if (test('madelfil')) then
      open (69, access='direct', recl=lrecabmad, file='abvmad.unformatted', form='unformatted')
    else
      allocate(t_madel%avmad(naez*naez, lmpotd, lmpotd), stat=ierr)
      allocate(t_madel%bvmad(naez*naez, lmpotd), stat=ierr)
    end if

    ! --> calculate the gaunt coefficients

    call madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)

    ! --> calculate the madelung coefficients to be used for VMAD
    ! call MADELCOEF with first arg. .FALSE. = 3D case

    do iq1 = 1, naez
      do iq2 = 1, naez
        call madelcoef(.false.,lpot,avmad,bvmad,madelsmat(1,iq1,iq2),cleb,icleb,    &
          iend,lpotd,lmpotd,lmxspd,nclebd)

        irec = iq2 + naez*(iq1-1)
        if (test('madelfil')) then
          write (69, rec=irec) avmad, bvmad
        else
          t_madel%avmad(irec,:,:) = avmad(:,:)
          t_madel%bvmad(irec,:) = bvmad(:)
        end if
        ! -----------------------------------------------------------------------
        if ((iq1<=6) .and. (iq2<=6)) then
          smat1(iq1, iq2) = avmad(1, 1)
          smat2(iq1, iq2) = bvmad(1)
        end if
        ! -----------------------------------------------------------------------
      end do
    end do
    if (test('madelfil')) close (69)

    if (iprint<1) return
    ! ======================================================================

    call madel3out(iprint, naez, lrecabmad, smat1, smat2, lmpotd)

  end subroutine madelung3d

end module mod_madelung3d
