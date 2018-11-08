!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Subroutine dealing with the printing of the needed kkrflex files
!> Author: 
!> Subroutine dealing with the printing of the needed kkrflex files for the
!> realization of an impurity calculation with the KKRimp software package.
!> It specifically prints the following files:
!>
!> - `kkrflex_tmat`
!> - `kkrflex_intercell_ref`
!> - `kkrflex_intercell_cmoms`
!------------------------------------------------------------------------------------
module mod_writekkrflex

contains

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine dealing with the printing of the needed kkrflex files
  !> Author: 
  !> Category: input-output, KKRhost
  !> Deprecated: False 
  !> Subroutine dealing with the printing of the needed kkrflex files for the
  !> realization of an impurity calculation with the KKRimp software package.
  !> It specifically prints the following files:
  !> 
  !> - `kkrflex_tmat`
  !> - `kkrflex_intercell_ref`
  !> - `kkrflex_intercell_cmoms`
  !-------------------------------------------------------------------------------
  subroutine writekkrflex(natomimp,nspin,ielast,lmpot,alat,natyp,kshape,vbc,atomimp,&
    hostimp,noq,zat,kaoez,conc,cmom,cminst,vinters,nemb,naez)

    use :: mod_types, only: t_tgmat
    use :: mod_wunfiles, only: t_params, read_angles
    use :: mod_version_info
    use :: mod_md5sums
    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_rotatespinframe, only: rotatematrix

    implicit none

    ! .. Input variables
    integer, intent (in) :: nemb   !! Number of 'embedding' positions
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: kshape !! Exact treatment of WS cell
    integer, intent (in) :: ielast
    integer, intent (in) :: natomimp !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    real (kind=dp), intent (in) :: alat !! Lattice constant in a.u.
    integer, dimension (naez), intent (in) :: noq !! Number of diff. atom types located
    integer, dimension (natomimp), intent (in) :: atomimp
    integer, dimension (0:natyp), intent (in) :: hostimp
    integer, dimension (natyp, naez+nemb), intent (in) :: kaoez !! Kind of atom at site in elem. cell
    real (kind=dp), dimension (2), intent (in) :: vbc !! Potential constants
    real (kind=dp), dimension (natyp), intent (in) :: zat !! Nuclear charge
    real (kind=dp), dimension (natyp), intent (in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension (lmpot, natyp), intent (in) :: cmom
    real (kind=dp), dimension (lmpot, natyp), intent (in) :: cminst
    real (kind=dp), dimension (lmpot, naez), intent (in) :: vinters
    ! .. Local variables
    integer :: ispin, ie, i1, iatom, irec, i, lm
    real (kind=dp), dimension (natyp) :: theta, phi
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: tmat0
    ! .. External Functions
    logical :: opt
    external :: opt

    write (1337, *) 'KKRFLEX WRITEOUT'
    write (1337, *) write_kkrimp_input

    if (write_kkrimp_input) then
      open (6699, file='kkrflex_tmat', status='unknown')
      call version_print_header(6699, '; '//md5sum_potential//'; '//md5sum_shapefun)
      write (6699, *) '#', natomimp, nspin, ielast, lmmaxd, korbit
      if (t_tgmat%tmat_to_file) then
        open (69, access='direct', recl=wlength*4*lmmaxd*lmmaxd, file='tmat', form='unformatted')
      end if

      ! read in non-collinear angles
      call read_angles(t_params, natyp, theta, phi)

      do iatom = 1, natomimp
        i1 = atomimp(iatom)
        do ispin = 1, nspin/(1+korbit)
          do ie = 1, ielast
            if (i1<=natyp) then
              irec = ie + ielast*(ispin-1) + ielast*nspin/(1+korbit)*(i1-1)
              if (t_tgmat%tmat_to_file) then
                read (69, rec=irec) tmat0
              else
                stop 'WRONG tmat_to_file for KKRFLEX writeout!!'
                ! not correctly read in with this option, to fix this communication is needed
                ! tmat0(:,:) = t_tgmat%tmat(:,:,irec)
              end if
              if (korbit==1) then
                ! perform here a transformation from the local to the global
                ! spin-frame of reference
                call rotatematrix(tmat0, theta(i1), phi(i1), lmgf0d, 0)
              end if
            else
              tmat0 = (0d0, 0d0)
            end if
            write (6699, '(4I12,50000E25.16)') iatom, ispin, ie, 0, tmat0
          end do
        end do
      end do
      close (69)
      close (6699)

      open (91, file='kkrflex_intercell_ref', status='unknown')
      call version_print_header(91, '; '//md5sum_potential//'; '//md5sum_shapefun)
      write (91, *) '# Intercell potential of each atom'
      write (91, *) '# '
      write (91, *) '# NATOMIMP', natomimp
      write (91, *) '# lmpot', lmpot
      write (91, *) '# KSHAPE', kshape
      write (91, *) '# NATOMIMP, lmpot, ALAT VBC(1), VBC(2)'
      write (91, '(2I12,10F25.16)') natomimp, lmpot, alat, vbc(1), vbc(2)
      do iatom = 1, natomimp       ! Bauer 2011-10-11
        i = atomimp(iatom)
        write (1337, *) 'ac2', i, hostimp(i), lmpot, (vinters(lm,i), lm=1, lmpot)
        write (91, '(5000F25.16)')(vinters(lm,i), lm=1, lmpot)
      end do
      close (91)

      open (91, file='kkrflex_intercell_cmoms', status='unknown')
      call version_print_header(91, '; '//md5sum_potential//'; '//md5sum_shapefun)
      write (91, *) '# Charge moments of each atom in the unit cell'
      write (91, *) '# Values given are CMOM + CMINST'
      write (91, *) '# First colums is the core charge other'
      write (91, *) '# other colums the charge moments'
      write (91, *) '# NATOMIMP', natomimp
      write (91, *) '# lmpot', lmpot
      write (91, *) '# KSHAPE', kshape

      do iatom = 1, natomimp
        i = atomimp(iatom)
        i1 = kaoez(1, i)
        write (1337, *) 'NOQ', i, noq(i)
        if (noq(i)/=1 .and. noq(i)/=0) stop '[vmadelblk] VIRATOMS: NOQ/=1'
        if (noq(i)==0) then
          write (91, '(5000F25.16)') 0.0d0, (0.0d0, lm=1, lmpot)
        else
          if (kshape/=0) then
            write (91, '(5000F25.16)') zat(i1), ((cmom(lm,i1)+cminst(lm,i1))*conc(i1), lm=1, lmpot)
          else
            write (91, '(5000F25.16)') zat(i1), (cmom(lm,i1)*conc(i1), lm=1, lmpot)
          end if
        end if
      end do
      close (91)
    end if

  end subroutine writekkrflex

end module mod_writekkrflex
