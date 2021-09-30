!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Testing the dimension of several arrays
!> Author: 
!> Testing the dimension of several arrays
!------------------------------------------------------------------------------------
!> @note Jonathan Chico: Some of these tests seem unnecessary with the removal of
!> the `inc.p`
!> @endnote
!------------------------------------------------------------------------------------
module mod_testdim

contains

  !-------------------------------------------------------------------------------
  !> Summary: Testing the dimension of several arrays
  !> Author: 
  !> Category: sanity-check, KKRhost 
  !> Deprecated: False 
  !> Testing the dimension of several arrays
  !-------------------------------------------------------------------------------
  !> @note Jonathan Chico: Some of these tests seem unnecessary with the removal of
  !> the `inc.p`
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine testdim(nspin,naez,nemb,natyp,ins,insref,nref,irns,nlayer,krel,nspind, &
    nprincd,knosph,irnsd,korbit,invmod)

    use :: mod_runoptions, only: calc_complex_bandstructure, use_Chebychev_solver, use_cont, use_virtual_atoms, decouple_spin_cheby
    implicit none

    integer, intent (in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: nref   !! Number of diff. ref. potentials
    integer, intent (in) :: krel   !! Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: irnsd
    integer, intent (in) :: insref !! INS for reference pot. (usual 0)
    integer, intent (in) :: knosph !! switch for spherical/non-spherical (0/1) program.
    integer, intent (in) :: korbit !! Spin-orbit/non-spin-orbit (1/0) added to the Schroedinger or SRA equations. Works with FP. KREL and KORBIT cannot be both non-zero.
    integer, intent (in) :: nspind !! KREL+(1-KREL)*(NSPIN+1)
    integer, intent (in) :: nprincd !! Number of principle layers, set to a number >= NRPINC in output of main0
    integer, intent (in) :: invmod
    ! .. In/Out variables
    integer, intent (inout) :: nemb !! Number of 'embedding' positions
    integer, intent (inout) :: nlayer !! Number of principal layer
    integer, dimension (natyp), intent (inout) :: irns !! Position of atoms in the unit cell in units of bravais vectors

    integer :: stop_mark
    integer :: i, j

    ! ---> dimension tests

    write (1337, 120)

    stop_mark = 0
    if ((nspin>nspind) .and. (krel==0)) then
      write (6, *) 'There is an inconsistenciy between spin polarised calculation and relativistic options'
      stop_mark = 1
    end if
    if (max(ins,insref)>knosph) then
      write (6, *) 'Please, change the parameter insd in', ' the inputcard to', max(ins, insref)
      stop_mark = 1
    end if
    j = 1
    do i = 1, natyp
      j = max(j, irns(i))
    end do
    if (ins==0 .and. j>1) then
      write (6, *) 'IRNS(*) is set to 1 in case of ', 'spherical potential treatment.'
      do i = 1, natyp
        irns(i) = 1
      end do
      j = 1
    end if

    if (j>irnsd) then
      write (6, *) 'Please, change the parameter irnsd in', ' the inputcard to', j
      stop_mark = 1
    end if

    if (.not. use_virtual_atoms) then
      if (nref>natyp) then
        write (6, *) 'There are some inconsistencies in the input file./', ' nref(=', nref, ') is greater than natyp (=', natyp, ').'
        stop_mark = 1
      end if
    end if

    if ((krel==1) .and. (korbit==1)) then
      write (6, *) 'Full relativistic for ASA and new SO solver', 'KREL', krel, 'KORBIT', korbit
      stop_mark = 1
    end if

    if (.not. use_Chebychev_solver .and. korbit==1) then
      write (6, *) 'Option NEWSOSOL not found, change KORBIT in the inputcard from', korbit, 'to 0'
      stop_mark = 1
    end if

    if (use_Chebychev_solver .and. korbit==0) then
      if (decouple_spin_cheby) then
        write (6, *) 'Using NEWSOSOL for decoupled spin channels.'
      else
        write (6, *) 'Using option NEWSOSOL, change KORBIT in the inputcard from', korbit, 'to 1'
        stop_mark = 1
      end if
    end if

    if (calc_complex_bandstructure) then
      write (6, *) 'Use option < calc_complex_bandstructure > only for eigenvalue determination.'
      stop_mark = 1
    end if

    if (use_cont) then
      nemb = 0
      write (6, *) 'No usage of embedding points. NEMB is set to ', nemb, '.'
    end if

    if (.not. ( (invmod==0) .or. (invmod==3) ) ) then
      ! -------------------------------------------------------------------------
      ! Constants for O(N) algorithm for matrix inversion
      ! -------------------------------------------------------------------------
      nlayer = naez/nprincd
      write (1337, 100) nprincd, nlayer
      write (1337, 110)
      ! ignore this test if full inversion is done
      if (.not. (invmod==0)) then
        if (nlayer*nprincd/=naez) then
          write (6, *) 'NLAYER*NPRINCD ( = ', nlayer*nprincd, ').NE.NAEZ ( = ', naez, ')'
          stop_mark = 1
        end if
      end if

    end if
    ! ----------------------------------------------------------------------------
    ! STOP IF A DIMENSION ERROR OCCURED
    ! ----------------------------------------------------------------------------
    if (stop_mark>0) stop 'STOP : Dimension Error.'
100 format (' NPRINCD  NLAYER', /, 2i8)
110 format (2(7('-'),'+'), 63('-'))
120 format (' Dimension and Input Data CHECK')
    return
  end subroutine testdim

end module mod_testdim
