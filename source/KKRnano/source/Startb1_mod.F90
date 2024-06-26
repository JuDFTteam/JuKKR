!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module Startb1_mod
!-------------------------------------------------------------------------------
!> Summary: A wrapper for the subroutine STARTB1
!> Author: Elias Rabel, Alexander R Thiess, Rudolf Zeller, Paul F Baumeister, Marcel Bornemann
!> Category: KKRnano, initialization, input-output, radial-grid, shape-functions
!-------------------------------------------------------------------------------
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: startb1_wrapper_new
  
  contains

! determine the record length
  subroutine startb1_wrapper_new(alat, nspin, EFERMI, Zat, naez, nowrite)
    use read_formatted_shapefun_mod, only: ShapefunFile, create, destroy
    use InputParams_mod, only: InputParams
    double precision, intent(in) :: alat
    integer, intent(in) :: naez
    double precision, intent(out) :: EFERMI
    integer, intent(in) :: nspin
    double precision, intent(out):: Zat(:)
    logical, intent(in) :: nowrite

    integer :: ntcell(naez)
    integer :: max_reclen, max_reclen_mesh
    type(ShapefunFile) :: sfile

!    write(*,*) 'read shapefun file'
    ! read the complete shapefun file to -> sfile
    open(91, file='shapefun', form='formatted', status='old', action='read')
    call create(sfile, 91) ! create_read_ShapefunFile
    close(91)
!    write(*,*) 'finished reading shapefun file'

    ! write atoms file and get maximum record lengths for vpotnew file and meshes file
    call write_atoms_file(alat, nspin, ntcell, Zat, naez, sfile, EFERMI, max_reclen, &
                          max_reclen_mesh, nowrite)

    ! routine to write binary potential and binary meshes
    call write_binary_potential(alat, nspin, ntcell, naez, sfile, max_reclen, &
                                max_reclen_mesh, nowrite)

    call destroy(sfile)

  endsubroutine ! startb1_wrapper_new

!------------------------------------------------------------------------------
!> Write the 'atoms' file (= binary analogue to atominfo - used for parallel reading in kkr2)
!> and determine the record length ('max_reclen') for the binary direct access potential file.
  subroutine write_atoms_file(alat, nspin, ntcell, Zat, naez, sfile, EFERMI, max_reclen, max_reclen_mesh, nowrite)

    use read_formatted_mod, only: PotentialEntry, create, destroy
    use read_formatted_shapefun_mod, only: ShapefunFile
    use BasisAtom_mod, only: BasisAtom, createBasisAtom, openBasisAtomDAFile, writeBasisAtomDA, closeBasisAtomDAFile, destroy
    use RadialMeshData_mod, only: RadialMeshData, getMinReclenMesh, create, destroy
    double precision, intent(in) :: alat
    integer, intent(in) :: nspin
    double precision, intent(out) :: Zat(:)
    integer, intent(in) :: naez
    integer, intent(out) :: ntcell(:)
    type(ShapefunFile), intent(in) :: sfile
    double precision, intent(out) :: EFERMI
    integer, intent(out) :: max_reclen
    integer, intent(out) :: max_reclen_mesh
    logical, intent(in) :: nowrite

    integer :: iatom, ispin, n_warn_backfolding
    integer :: reclen, lpot_atom, ib
    logical :: bf ! backfold
    type(PotentialEntry), allocatable :: pe(:,:)
    type(BasisAtom) :: atom
    type(RadialMeshData) :: mesh
    integer, parameter :: fu = 13, fa = 37
    integer :: cell_index, n_warn_alat_differs
    double precision :: radius_muffin_tin
    
    max_reclen = 0
    max_reclen_mesh = 0
    reclen = 0
    n_warn_backfolding = 0
    n_warn_alat_differs = 0
    
    bf = (sfile%ncell < naez)
    ib = 1 ; if (bf) ib = sfile%ncell
    allocate(pe(2,ib))

    if (.not. nowrite) call openBasisAtomDAFile(atom, fa, 'bin.atoms', action='write')
    
    open(unit=fu, file='potential', status='old', form='formatted', action='read')
    do iatom = 1, naez

      ib = 1 ; if (bf) ib = modulo(iatom - 1, sfile%ncell) + 1
      do ispin = 1, nspin
        if (iatom <= sfile%ncell) &
        call create(pe(ispin,ib), fu, iatom) ! create_read_PotentialEntry

        if (iatom == 1 .and. ispin == 1) EFERMI = pe(ispin,ib)%header%EFERMI ! take approximate Fermi energy from 1st potential entry

  !        The nuclear charge Z is now solely determined by the 'potential' file, 'atominfo' is no longer needed

        Zat(iatom) = pe(ispin,ib)%header%Z_nuclear

  !        The parameter ntcell is now automatically set to the index of the atom (1st atom -> 1st pe in shapefun file,
  !        2nd atom -> 2nd pe in shapefun file), 'atominfo' is no longer needed

        ntcell(iatom) = iatom

  !        The muffin-tin radius (formerly RMT in atominfo) is set to the value of RMT in the 'potential' file,
  !        'atominfo' is no longer needed

        radius_muffin_tin = pe(ispin,ib)%header%RMT

  !        ! do some consistency checks
  !        if (abs(Zat(iatom) - pe(ispin,ib)%header%Z_nuclear) > 1.d-8) then
  !          write(*,*) "ERROR: Mismatch of nuclear charge between atominfo and potential file for pe: ", iatom
  !          write(*,*) Zat(iatom), pe(ispin,ib)%header%Z_nuclear
  !          STOP
  !        endif

        if (abs(alat - pe(ispin,ib)%header%alat) > 1.d-8) n_warn_alat_differs = n_warn_alat_differs + 1
      enddo ! ispin

      lpot_atom = lmpot_to_lpot(pe(1,ib)%sblock%lmpot)
      CHECKASSERT( (lpot_atom+1)**2 == pe(1,ib)%sblock%lmpot)

      ! number of radial points must be the same for every spin direction
      CHECKASSERT( pe(1,ib)%sblock%irt1p == pe(nspin,ib)%sblock%irt1p )
      CHECKASSERT( pe(1,ib)%sblock%irns == pe(nspin,ib)%sblock%irns )
      CHECKASSERT( pe(1,ib)%sblock%lmpot == pe(nspin,ib)%sblock%lmpot )

      call createBasisAtom(atom, 1, lpot_atom, nspin, pe(1,ib)%sblock%IRT1P - pe(1,ib)%sblock%IRNS, pe(1,ib)%sblock%IRT1P)  ! create dummy basis atom

      inquire(iolength=reclen) atom%potential%vins, &
                               atom%potential%visp, &
                               atom%core%ecore

      max_reclen = max(max_reclen, reclen)

      ! write file 'bin.atoms'
      atom%atom_index = iatom
      atom%cell_index = ntcell(iatom)
      atom%Z_nuclear = Zat(iatom)

      atom%radius_muffin_tin = radius_muffin_tin

      atom%core%NCORE = 0
      atom%core%LCORE = 0
      atom%core%ITITLE = 0

      do ispin = 1, nspin
        atom%core%NCORE(ispin) = pe(ispin,ib)%csblock%NCORE
        atom%core%LCORE(:,ispin) = pe(ispin,ib)%csblock%LCORE
        atom%core%ITITLE(:,ispin) = pe(ispin,ib)%header%ITITLE
      enddo

      if (.not. nowrite) call writeBasisAtomDA(atom, fa, iatom)
      
      call destroy(atom)

      cell_index = modulo(ntcell(iatom) - 1, sfile%ncell) + 1
      if (cell_index < ntcell(iatom)) n_warn_backfolding = n_warn_backfolding + 1
      
      call create(mesh, pe(1,ib)%sblock%IRT1P, sfile%mesh(cell_index)%npan+1, sfile%mesh(cell_index)%meshn, sfile%shapes(cell_index)%nfu) ! createRadialMeshData

      if(.not. bf) call destroy(pe(:,1))
      
      ! determine maximal record length for meshes.0 file
      ! this is a bit of a hack
      max_reclen_mesh = max(getMinReclenMesh(mesh), max_reclen_mesh)
     
      ! cleanup
      call destroy(mesh)

    enddo ! iatom ! end loop over atoms
    close(fu)
    
    if (.not. nowrite) call closeBasisAtomDAFile(fa)
    
    call destroy(pe)

    if (n_warn_alat_differs > 0) &
      warn(6, "In"+n_warn_alat_differs/nspin+"potential files ALAT is not the same as in the input!")

    if (n_warn_backfolding > 0) &
      warn(6, "In"+n_warn_backfolding+"cases the shapefun and potential index has been backfolded into [1,"+sfile%ncell-"]!")

  endsubroutine ! write_atoms_file

  integer function lmpot_to_lpot(lmpot)
    integer, intent(in) :: lmpot

    lmpot_to_lpot = int( sqrt(dble(lmpot)) - 1.d0 + 0.01d0)
  endfunction ! lmpot_to_lpot


  subroutine write_binary_potential(alat, nspin, ntcell, naez, sfile, max_reclen, max_reclen_mesh, nowrite)
    use read_formatted_mod, only: PotentialEntry, create, destroy
    use read_formatted_shapefun_mod, only: ShapefunFile
    use BasisAtom_mod, only: BasisAtom, create, destroy
    use BasisAtom_mod, only: openBasisAtomPotentialDAFile, writeBasisAtomPotentialDA, closeBasisAtomPotentialDAFile
#ifndef TASKLOCAL_FILES
    use BasisAtom_mod, only: openBasisAtomPotentialIndexDAFile, writeBasisAtomPotentialIndexDA, closeBasisAtomPotentialIndexDAFile
    use RadialMeshData_mod, only: openRadialMeshDataIndexDAFile, writeRadialMeshDataIndexDA, closeRadialMeshDataIndexDAFile
#endif
    use RadialMeshData_mod, only: RadialMeshData, initRadialMesh, create, destroy
    use RadialMeshData_mod, only: openRadialMeshDataDAFile, writeRadialMeshDataDA, closeRadialMeshDataDAFile
    double precision, intent(in) :: alat
    integer, intent(in) :: nspin
    integer, intent(in) :: naez
    integer, intent(in) :: ntcell(*)
    type(ShapefunFile), intent(in) :: sfile
    integer, intent(in) :: max_reclen
    integer, intent(in) :: max_reclen_mesh
    logical, intent(in) :: nowrite

    type(BasisAtom) :: atom
    integer :: ib
    logical :: bf ! backfold
    type(PotentialEntry), allocatable :: pe(:,:)
    type(RadialMeshData) :: mesh
    integer, parameter :: fu=13
    integer :: cell_index, lpot, irmd, irmind, ipand, irns, irid, iatom, ispin

    bf = (sfile%ncell < naez)
    ib = 1 ; if (bf) ib = sfile%ncell
    allocate(pe(2,ib))
    
    open(unit=fu, file='potential', form='formatted', status='old', action='read')
    do iatom = 1, naez

      ib = 1 ; if (bf) ib = modulo(iatom - 1, sfile%ncell) + 1
      do ispin = 1, nspin
        if (iatom <= sfile%ncell) &
        call create(pe(ispin,ib), fu, iatom) ! create_read_PotentialEntry
      enddo ! spin

      cell_index = modulo(ntcell(iatom) - 1, sfile%ncell) + 1 ! backfold

      lpot = lmpot_to_lpot(pe(1,ib)%sblock%lmpot)
      irmd = pe(1,ib)%sblock%irt1p
      irns = pe(1,ib)%sblock%irns
      irmind = irmd - irns
      ipand = sfile%mesh(cell_index)%npan+1
      irid = sfile%mesh(cell_index)%meshn

      call create(atom, iatom, lpot, nspin, irmind, irmd) ! createBasisAtom
      call create(mesh, irmd, ipand, irid, sfile%shapes(cell_index)%nfu) ! createRadialMeshData

      ! set potential
      do ispin = 1, nspin
        atom%potential%vins(:,:,ispin) = pe(ispin,ib)%nsblocks%vins
        atom%potential%visp(:,ispin) = pe(ispin,ib)%sblock%visp
        atom%core%ecore(:,ispin) = pe(ispin,ib)%csblock%ecore
      enddo ! ispin
      
      if(.not. bf) call destroy(pe(:,1))

      ! initialise radial mesh
      call initRadialMesh(mesh, alat, sfile%mesh(cell_index)%xrn, sfile%mesh(cell_index)%drn, &
                          sfile%mesh(cell_index)%nm, irmd-irid, irns, &
                          irid, sfile%shapes(cell_index)%nfu, sfile%shapes(cell_index)%llmsp, &
                          sfile%shapes(cell_index)%thetas, pe(1,ib)%header%a_log_mesh, &
                          pe(1,ib)%header%b_log_mesh)
      
      if (.not. nowrite) then
      
        if (iatom == 1) then
#ifndef TASKLOCAL_FILES
          call openBasisAtomPotentialIndexDAFile(atom, 38, 'bin.vpotnew.0.idx', action='write')
          call openRadialMeshDataIndexDAFile(mesh, 93, "bin.meshes.0.idx")
#endif
          call openBasisAtomPotentialDAFile(atom, 39, 'bin.vpotnew.0', max_reclen, action='write')
          call openRadialMeshDataDAFile(mesh, 94, "bin.meshes.0", max_reclen_mesh)
        endif ! iatom == 1

        call writeBasisAtomPotentialDA(atom, 39, iatom)
        call writeRadialMeshDataDA(mesh, 94, iatom)

#ifndef TASKLOCAL_FILES
        call writeBasisAtomPotentialIndexDA(atom, 38, iatom, max_reclen)
        call writeRadialMeshDataIndexDA(mesh, 93, iatom, max_reclen_mesh)
#endif
      endif ! not nowrite

      ! cleanup
      call destroy(atom)
      call destroy(mesh)

    enddo ! endloop over atoms

    close(fu)

    call destroy(pe)
    deallocate(pe)

    if (.not. nowrite) then
    
      call closeBasisAtomPotentialDAFile(39)
      call closeRadialMeshDataDAFile(94)

#ifndef TASKLOCAL_FILES
      call closeBasisAtomPotentialIndexDAFile(38)
      call closeRadialMeshDataIndexDAFile(93)
#endif
    endif ! not nowrite

  endsubroutine ! write_binary_potential

endmodule ! Startb1_mod
