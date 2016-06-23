#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!-------------------------------------------------------------------------------
!> A wrapper for the subroutine STARTB1
!> Returns: EFERMI (from first entry in potential file)

! determine the record length
subroutine STARTB1_wrapper_new(alat,NSPIN, &
                           NTCELL, &
                           EFERMI,ZAT, radius_muffin_tin, &
                           NAEZD)

  use read_formatted_shapefun_mod
  implicit none

  ! Arguments
  double precision, intent(in) :: alat
  INTEGER, INTENT(IN) :: NAEZD
  DOUBLE PRECISION :: EFERMI
  INTEGER :: NSPIN
  DOUBLE PRECISION, dimension(*) :: ZAT
  double precision, dimension(naezd), intent(inout) :: radius_muffin_tin
  INTEGER, dimension(*) :: NTCELL

  integer :: max_reclen, max_reclen_mesh
  type (ShapefunFile) :: sfile

  ! read the complete shapefun file to -> sfile
  open (91, file='shapefun', status='old', form='formatted')
  call create_read_shapefunfile(sfile, 91)
  close (91)

  ! write atoms file and get maximum record lengths for vpotnew file and meshes file
  call write_atoms_file(alat, NSPIN, &
                            NTCELL, &
                            ZAT, radius_muffin_tin, &
                            NAEZD, sfile, EFERMI, max_reclen, max_reclen_mesh)

  ! routine to write binary potential and binary meshes
  call write_binary_potential(alat, NSPIN, &
                                  NTCELL, &
                                  NAEZD, sfile, max_reclen, max_reclen_mesh)

  call destroy_shapefunfile(sfile)

end subroutine

!------------------------------------------------------------------------------
!> Write the 'atoms' file (= binary analogue to atominfo - used for parallel reading in kkr2)
!> and determine the record length ('max_reclen') for the binary direct access potential file.
subroutine write_atoms_file(alat, NSPIN, &
                            NTCELL, &
                            ZAT, radius_muffin_tin, &
                            NAEZD, sfile, EFERMI, max_reclen, max_reclen_mesh)
  use read_formatted_mod
  use read_formatted_shapefun_mod
  use BasisAtom_mod
  use RadialMeshData_mod
  implicit none

  ! Arguments
  double precision, intent(in) :: alat
  INTEGER :: NSPIN
  DOUBLE PRECISION, dimension(*) :: ZAT
  double precision, dimension(naezd), intent(inout) :: radius_muffin_tin
  integer, intent(in) :: NAEZD
  INTEGER, dimension(*) :: NTCELL
  type (Shapefunfile), intent(in) :: sfile
  double precision :: EFERMI
  integer, intent(out) :: max_reclen
  integer, intent(out) :: max_reclen_mesh

  integer :: iatom
  integer :: ispin
  integer :: lpot_atom
  integer :: reclen

  type (PotentialEntry) :: entry(2)
  type (BasisAtom) :: atom
  type (RadialMeshData) :: mesh
  integer, parameter :: UNIT = 13
  integer :: cell_index
  integer, external :: lmpot_to_lpot

  max_reclen = 0
  max_reclen_mesh = 0
  reclen = 0

  call openBasisAtomDAFile(atom, 37, 'atoms')

  open (UNIT, file='potential', status='old', form='formatted')
    do iatom = 1, naezd

      do ispin = 1, nspin
        call create_read_PotentialEntry(entry(ispin), UNIT)

        ! take approximate Fermi energy from 1st potential entry
        if (iatom == 1 .and. ispin == 1) then
          EFERMI = entry(ispin)%header%EFERMI
        endif

!        The nuclear charge Z is now solely determined by the 'potential' file, 'atominfo' is no longer needed

         ZAT(iatom)=entry(ispin)%header%Z_nuclear

!        The parameter NTCELL is now automatically set to the index of the atom (1st atom -> 1st entry in shapefun file,
!        2nd atom -> 2nd entry in shapefun file), 'atominfo' is no longer needed

         NTCELL(iatom)=iatom

!        The muffin-tin radius (formerly RMT in atominfo) is set to the value of RMT in the 'potential' file,
!        'atominfo' is no longer needed

         radius_muffin_tin(iatom)=entry(ispin)%header%RMT

!        ! do some consistency checks
!        if (abs(ZAT(iatom) - entry(ispin)%header%Z_nuclear) > 1.d-8) then
!          write(*,*) "ERROR: Mismatch of nuclear charge between atominfo and potential file for entry: ", iatom
!          write(*,*) ZAT(iatom), entry(ispin)%header%Z_nuclear
!          STOP
!        endif


        if (abs(alat - entry(ispin)%header%alat) > 1.d-8) then
          write(*,*) "WARNING: alat in potential file is not the same as in input.conf."
        endif
      enddo

      lpot_atom = lmpot_to_lpot(entry(1)%sblock%lmpot)
      CHECKASSERT( (lpot_atom + 1)**2 == entry(1)%sblock%lmpot)

      ! number of radial points must be the same for every spin direction
      CHECKASSERT( entry(1)%sblock%irt1p == entry(nspin)%sblock%irt1p )
      CHECKASSERT( entry(1)%sblock%irns == entry(nspin)%sblock%irns )
      CHECKASSERT( entry(1)%sblock%lmpot == entry(nspin)%sblock%lmpot )

      call createBasisAtom(atom, 1, lpot_atom, nspin, entry(1)%sblock%IRT1P - entry(1)%sblock%IRNS, entry(1)%sblock%IRT1P)  ! create dummy basis atom

      inquire (iolength = reclen) atom%potential%VINS, &
                                  atom%potential%VISP, &
                                  atom%core%ECORE

      max_reclen = max(max_reclen, reclen)

      ! write file 'atoms'
      atom%atom_index = iatom
      atom%cell_index = NTCELL(iatom)
      atom%Z_nuclear = ZAT(iatom)

      atom%radius_muffin_tin = radius_muffin_tin(iatom)

      atom%core%NCORE = 0
      atom%core%LCORE = 0
      atom%core%ITITLE = 0

      do ispin = 1, nspin
        atom%core%NCORE(ispin) = entry(ispin)%csblock%NCORE
        atom%core%LCORE(:, ispin) = entry(ispin)%csblock%LCORE
        atom%core%ITITLE(:, ispin) = entry(ispin)%header%ITITLE
      enddo

      call writeBasisAtomDA(atom, 37, iatom)

      ! determine maximal record length for meshes.0 file
      ! this is a bit of a hack
      cell_index = NTCELL(iatom)
      CHECKASSERT (cell_index <= sfile%NCELL .and. cell_index > 0)
      call createRadialMeshData(mesh, entry(1)%sblock%IRT1P, sfile%mesh(cell_index)%npan + 1, sfile%mesh(cell_index)%meshn, sfile%shapes(iatom)%nfu)
      max_reclen_mesh = max(getMinReclenMesh(mesh), max_reclen_mesh)
      ! cleanup
      call destroyRadialMeshData(mesh)

      call destroyBasisAtom(atom)

      do ispin = 1, nspin
        call destroy_PotentialEntry(entry(ispin))
      enddo

    enddo ! end loop over atoms
  close(UNIT)

  call closeBasisAtomDAFile(37)

end subroutine

integer function lmpot_to_lpot(lmpot)
  implicit none
  integer, intent(in) :: lmpot

  lmpot_to_lpot = int( sqrt(dble(lmpot)) - 1.0d0 + 0.01d0)

end function


subroutine write_binary_potential(alat, NSPIN, &
                                  NTCELL, &
                                  NAEZD, sfile, max_reclen, max_reclen_mesh)

  use read_formatted_mod
  use read_formatted_shapefun_mod
  use BasisAtom_mod
  use RadialMeshData_mod
  implicit none

  ! Arguments
  double precision, intent(in) :: alat
  INTEGER, intent(in) :: NSPIN
  integer, intent(in) :: NAEZD
  INTEGER, intent(in), dimension(*) :: NTCELL
  type (Shapefunfile), intent(in) :: sfile
  integer, intent(in) :: max_reclen
  integer, intent(in) :: max_reclen_mesh

  integer :: iatom
  integer :: ispin

  type (PotentialEntry) :: entry(2)
  type (BasisAtom) :: atom
  type (RadialMeshData) :: meshdata
  integer, parameter :: UNIT = 13
  integer :: cell_index

  integer lpot, irmd, irmind, ipand, irns, irid
  integer, external :: lmpot_to_lpot

  open (UNIT, file='potential', status='old', form='formatted')

    do iatom = 1, naezd

      do ispin = 1, nspin
        call create_read_PotentialEntry(entry(ispin), UNIT)
      enddo

      cell_index = NTCELL(iatom)

      lpot = lmpot_to_lpot(entry(1)%sblock%lmpot)
      irmd = entry(1)%sblock%irt1p
      irns = entry(1)%sblock%irns
      irmind = irmd - irns
      ipand = sfile%mesh(cell_index)%npan + 1
      irid = sfile%mesh(cell_index)%meshn

      call createBasisAtom(atom, iatom, lpot, nspin, irmind, irmd)
      call createRadialMeshData(meshdata, irmd, ipand, irid, sfile%shapes(cell_index)%nfu)

      ! set potential
      do ispin = 1, nspin
        atom%potential%VINS(:,:,ispin) = entry(ispin)%nsblocks%VINS
        atom%potential%VISP(:, ispin) = entry(ispin)%sblock%VISP
        atom%core%ECORE(:,ispin) = entry(ispin)%csblock%ECORE
      enddo

      ! initialise radial mesh
      call initRadialMesh(meshdata, alat, sfile%mesh(cell_index)%xrn, sfile%mesh(cell_index)%drn, sfile%mesh(cell_index)%nm, irmd - irid, irns, &
                          sfile%shapes(cell_index)%nfu, sfile%shapes(cell_index)%llmsp, sfile%shapes(cell_index)%thetas)

      if (iatom == 1) then
#ifndef TASKLOCAL_FILES
        call openBasisAtomPotentialIndexDAFile(atom, 38, 'vpotnew.0.idx')
        call openRadialMeshDataIndexDAFile(meshdata, 93, "meshes.0.idx")
#endif
        call openBasisAtomPotentialDAFile(atom, 39, 'vpotnew.0', max_reclen)
        call openRadialMeshDataDAFile(meshdata, 94 , "meshes.0", max_reclen_mesh)
      endif

#ifndef TASKLOCAL_FILES
      call writeBasisAtomPotentialIndexDA(atom, 38, iatom, max_reclen)
      call writeRadialMeshDataIndexDA(meshdata, 93, iatom, max_reclen_mesh)
#endif

      call writeBasisAtomPotentialDA(atom, 39, iatom)
      call writeRadialMeshDataDA(meshdata, 94, iatom)

      ! cleanup
      call destroyBasisAtom(atom)
      call destroyRadialMeshData(meshdata)

      do ispin = 1, nspin
        call destroy_PotentialEntry(entry(ispin))
      enddo

    enddo ! end loop over atoms

  call closeBasisAtomPotentialDAFile(39)
  call closeRadialMeshDataDAFile(94)

  close(UNIT)

#ifndef TASKLOCAL_FILES
  call closeRadialMeshDataIndexDAFile(93)
  call closeBasisAtomPotentialIndexDAFile(38)
#endif

end subroutine
