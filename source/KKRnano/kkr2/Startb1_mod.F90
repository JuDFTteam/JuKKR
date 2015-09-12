#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!-------------------------------------------------------------------------------
!> A wrapper for the subroutine STARTB1
!> Returns: EFERMI (from first entry in potential file)

module Startb1_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  
  public :: startb1_wrapper_new
  
  contains

! determine the record length
  subroutine startb1_wrapper_new(alat, nspin, ntcell, EFERMI, ZAT, radius_muffin_tin, naezd, nowrite)
    use read_formatted_shapefun_mod, only: ShapefunFile, create_read_ShapefunFile, destroy_ShapefunFile
    double precision, intent(in) :: alat
    integer, intent(in) :: naezd
    double precision, intent(out) :: EFERMI
    integer, intent(in) :: nspin
    double precision, intent(out):: ZAT(*)
    double precision, intent(inout) :: radius_muffin_tin(naezd)
    integer, intent(out) :: ntcell(*)
    logical, intent(in) :: nowrite

    integer :: max_reclen, max_reclen_mesh
    type(ShapefunFile) :: sfile

    ! read the complete shapefun file to -> sfile
    open (91, file='shapefun', form='formatted', status='old', action='read')
    call create_read_ShapefunFile(sfile, 91)
    close(91)

    ! write atoms file and get maximum record lengths for vpotnew file and meshes file
    call write_atoms_file(alat, nspin, ntcell, ZAT, radius_muffin_tin, naezd, sfile, EFERMI, max_reclen, max_reclen_mesh, nowrite)

    ! routine to write binary potential and binary meshes
    call write_binary_potential(alat, nspin, ntcell, naezd, sfile, max_reclen, max_reclen_mesh, nowrite)

    call destroy_ShapefunFile(sfile)

  endsubroutine startb1_wrapper_new

!------------------------------------------------------------------------------
!> Write the 'atoms' file (= binary analogue to atominfo - used for parallel reading in kkr2)
!> and determine the record length ('max_reclen') for the binary direct access potential file.
  subroutine write_atoms_file(alat, nspin, ntcell, ZAT, radius_muffin_tin, naezd, sfile, EFERMI, max_reclen, max_reclen_mesh, nowrite)

    use read_formatted_mod, only: PotentialEntry, create_read_PotentialEntry, destroy_PotentialEntry
    use read_formatted_shapefun_mod, only: ShapefunFile
    use BasisAtom_mod, only: BasisAtom, createBasisAtom, openBasisAtomDAFile, writeBasisAtomDA, closeBasisAtomDAFile, destroyBasisAtom
    use RadialMeshData_mod, only: RadialMeshData, getMinReclenMesh, createRadialMeshData, destroyRadialMeshData
    double precision, intent(in) :: alat
    integer, intent(in) :: nspin
    double precision :: ZAT(*)
    double precision, intent(inout) :: radius_muffin_tin(naezd)
    integer, intent(in) :: naezd
    integer, intent(out) :: ntcell(*)
    type(ShapefunFile), intent(in) :: sfile
    double precision, intent(out) :: EFERMI
    integer, intent(out) :: max_reclen
    integer, intent(out) :: max_reclen_mesh
    logical, intent(in) :: nowrite

    integer :: iatom
    integer :: ispin
    integer :: lpot_atom
    integer :: reclen
    type(PotentialEntry) :: pe(2)
    type(BasisAtom) :: atom
    type(RadialMeshData) :: mesh
    integer, parameter :: fu = 13
    integer :: cell_index
    integer :: n_warn_alat_differs

    max_reclen = 0
    max_reclen_mesh = 0
    reclen = 0
    
    n_warn_alat_differs = 0

    call openBasisAtomDAFile(atom, 37, 'atoms')

    open(unit=fu, file='potential', status='old', form='formatted', action='read')
    do iatom = 1, naezd

      do ispin = 1, nspin
        call create_read_PotentialEntry(pe(ispin), fu, iatom)

        if (iatom == 1 .and. ispin == 1) EFERMI = pe(ispin)%header%EFERMI ! take approximate Fermi energy from 1st potential entry

!        The nuclear charge Z is now solely determined by the 'potential' file, 'atominfo' is no longer needed

         ZAT(iatom) = pe(ispin)%header%Z_nuclear

!        The parameter ntcell is now automatically set to the index of the atom (1st atom -> 1st pe in shapefun file,
!        2nd atom -> 2nd pe in shapefun file), 'atominfo' is no longer needed

         ntcell(iatom) = iatom

!        The muffin-tin radius (formerly RMT in atominfo) is set to the value of RMT in the 'potential' file,
!        'atominfo' is no longer needed

         radius_muffin_tin(iatom) = pe(ispin)%header%RMT

!        ! do some consistency checks
!        if (abs(ZAT(iatom) - pe(ispin)%header%Z_nuclear) > 1.d-8) then
!          write(*,*) "ERROR: Mismatch of nuclear charge between atominfo and potential file for pe: ", iatom
!          write(*,*) ZAT(iatom), pe(ispin)%header%Z_nuclear
!          STOP
!        endif

        if (abs(alat - pe(ispin)%header%alat) > 1.d-8) n_warn_alat_differs = n_warn_alat_differs + 1
      enddo ! ispin

      lpot_atom = lmpot_to_lpot(pe(1)%sblock%lmpot)
      CHECKASSERT( (lpot_atom+1)**2 == pe(1)%sblock%lmpot)

      ! number of radial points must be the same for every spin direction
      CHECKASSERT( pe(1)%sblock%irt1p == pe(nspin)%sblock%irt1p )
      CHECKASSERT( pe(1)%sblock%irns == pe(nspin)%sblock%irns )
      CHECKASSERT( pe(1)%sblock%lmpot == pe(nspin)%sblock%lmpot )

      call createBasisAtom(atom, 1, lpot_atom, nspin, pe(1)%sblock%IRT1P - pe(1)%sblock%IRNS, pe(1)%sblock%IRT1P)  ! create dummy basis atom

      inquire (iolength = reclen) atom%potential%VINS, &
                                  atom%potential%VISP, &
                                  atom%core%ECORE

      max_reclen = max(max_reclen, reclen)

      ! write file 'atoms'
      atom%atom_index = iatom
      atom%cell_index = ntcell(iatom)
      atom%Z_nuclear = ZAT(iatom)

      atom%radius_muffin_tin = radius_muffin_tin(iatom)

      atom%core%NCORE = 0
      atom%core%LCORE = 0
      atom%core%ITITLE = 0

      do ispin = 1, nspin
        atom%core%NCORE(ispin) = pe(ispin)%csblock%NCORE
        atom%core%LCORE(:,ispin) = pe(ispin)%csblock%LCORE
        atom%core%ITITLE(:,ispin) = pe(ispin)%header%ITITLE
      enddo

      if (.not. nowrite) call writeBasisAtomDA(atom, 37, iatom)

      ! determine maximal record length for meshes.0 file
      ! this is a bit of a hack
      cell_index = ntcell(iatom)
      CHECKASSERT (1 <= cell_index .and. cell_index <= sfile%ncell)
      call createRadialMeshData(mesh, pe(1)%sblock%IRT1P, sfile%mesh(cell_index)%npan+1)
      max_reclen_mesh = max(getMinReclenMesh(mesh), max_reclen_mesh)
      
      ! cleanup
      call destroyRadialMeshData(mesh)
      call destroyBasisAtom(atom)
      do ispin = 1, nspin
        call destroy_PotentialEntry(pe(ispin))
      enddo ! ispin

    enddo ! iatom ! end loop over atoms
    close(fu)
    
    call closeBasisAtomDAFile(37)
    
    if (n_warn_alat_differs > 0) &
      warn(6, "In"+n_warn_alat_differs/nspin+"potential files ALAT is not the same as in the input!")

  endsubroutine ! write_atoms_file

  integer function lmpot_to_lpot(lmpot)
    integer, intent(in) :: lmpot

    lmpot_to_lpot = int( sqrt(dble(lmpot)) - 1.d0 + 0.01d0)
  endfunction


  subroutine write_binary_potential(alat, nspin, ntcell, naezd, sfile, max_reclen, max_reclen_mesh, nowrite)
    use read_formatted_mod, only: PotentialEntry, create_read_PotentialEntry, destroy_PotentialEntry
    use read_formatted_shapefun_mod, only: ShapefunFile
    use BasisAtom_mod, only: BasisAtom, createBasisAtom, destroyBasisAtom
    use BasisAtom_mod, only: openBasisAtomPotentialDAFile, writeBasisAtomPotentialDA, closeBasisAtomPotentialDAFile
#ifndef TASKLOCAL_FILES
    use BasisAtom_mod, only: openBasisAtomPotentialIndexDAFile, writeBasisAtomPotentialIndexDA, closeBasisAtomPotentialIndexDAFile
    use RadialMeshData_mod, only: openRadialMeshDataIndexDAFile, writeRadialMeshDataIndexDA, closeRadialMeshDataIndexDAFile
#endif
    use RadialMeshData_mod, only: RadialMeshData, initRadialMesh, createRadialMeshData, destroyRadialMeshData
    use RadialMeshData_mod, only: openRadialMeshDataDAFile, writeRadialMeshDataDA, closeRadialMeshDataDAFile
    double precision, intent(in) :: alat
    integer, intent(in) :: nspin
    integer, intent(in) :: naezd
    integer, intent(in) :: ntcell(*)
    type(ShapefunFile), intent(in) :: sfile
    integer, intent(in) :: max_reclen
    integer, intent(in) :: max_reclen_mesh
    logical, intent(in) :: nowrite

    type(BasisAtom) :: atom
    type(PotentialEntry) :: pe(2)
    type(RadialMeshData) :: meshdata
    integer, parameter :: fu=13
    integer :: cell_index, lpot, irmd, irmind, ipand, irns, irid, iatom, ispin

    open(unit=fu, file='potential', form='formatted', status='old', action='read')
    
    do iatom = 1, naezd

      do ispin = 1, nspin
        call create_read_PotentialEntry(pe(ispin), fu, iatom)
      enddo ! spin

      cell_index = ntcell(iatom)

      lpot = lmpot_to_lpot(pe(1)%sblock%lmpot)
      irmd = pe(1)%sblock%irt1p
      irns = pe(1)%sblock%irns
      irmind = irmd - irns
      ipand = sfile%mesh(cell_index)%npan+1
      irid = sfile%mesh(cell_index)%meshn

      call createBasisAtom(atom, iatom, lpot, nspin, irmind, irmd)
      call createRadialMeshData(meshdata, irmd, ipand)

      ! set potential
      do ispin = 1, nspin
        atom%potential%VINS(:,:,ispin) = pe(ispin)%nsblocks%VINS
        atom%potential%VISP(:,ispin) = pe(ispin)%sblock%VISP
        atom%core%ECORE(:,ispin) = pe(ispin)%csblock%ECORE
      enddo ! ispin

      ! initialise radial mesh
      call initRadialMesh(meshdata, alat, sfile%mesh(cell_index)%xrn, sfile%mesh(cell_index)%drn, sfile%mesh(cell_index)%nm, irmd-irid, irns)

      if (.not. nowrite) then
      
        if (iatom == 1) then
#ifndef TASKLOCAL_FILES
          call openBasisAtomPotentialIndexDAFile(atom, 38, 'vpotnew.0.idx')
          call openRadialMeshDataIndexDAFile(meshdata, 93, "meshes.0.idx")
#endif
          call openBasisAtomPotentialDAFile(atom, 39, 'vpotnew.0', max_reclen)
          call openRadialMeshDataDAFile(meshdata, 94 , "meshes.0", max_reclen_mesh)
        endif ! iatom == 1
      
        call writeBasisAtomPotentialDA(atom, 39, iatom)
        call writeRadialMeshDataDA(meshdata, 94, iatom)

#ifndef TASKLOCAL_FILES
        call writeBasisAtomPotentialIndexDA(atom, 38, iatom, max_reclen)
        call writeRadialMeshDataIndexDA(meshdata, 93, iatom, max_reclen_mesh)
#endif
      endif ! not nowrite

      ! cleanup
      call destroyBasisAtom(atom)
      call destroyRadialMeshData(meshdata)
      do ispin = 1, nspin
        call destroy_PotentialEntry(pe(ispin))
      enddo

    enddo ! endloop over atoms

    close(fu)
    
    if (.not. nowrite) then
    
      call closeBasisAtomPotentialDAFile(39)
      call closeRadialMeshDataDAFile(94)

#ifndef TASKLOCAL_FILES
      call closeBasisAtomPotentialIndexDAFile(38)
      call closeRadialMeshDataIndexDAFile(93)
#endif
    endif ! not nowrite

  endsubroutine write_binary_potential

endmodule ! 
