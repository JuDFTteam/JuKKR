! move core to potential??

! cell_index, cluster_index --> probably unnecessary, only for additional safety?

! TODO: output of Atom:
! historical reasons: write to two files: vpotnew and atoms

! potential and atom information should be written separately
  ! ==> separation in changing and nonchanging part

! -============= rather not ===============
! Potential not as member of datastructure??? and core as member of potential?
! new potential file format including all core info  (also nonchanging?) ?
! ==========================================
! reference cluster: own file: ref_clusters - or move creation to kkr2


! eliminates: ZAT, VREF, REFPOT?
! TODO: establish BasisAtomGroup???
! references to t-matrix and G_ref ???

! build in checks for compatibility of shapes with potential etc...

module BasisAtom_mod
  use CellData_mod
  use RefClusterData_mod
  use PotentialData_mod
  use AtomicCoreData_mod
  implicit none

  type BasisAtom

    integer :: atom_index !todo (position in atominfo/rbasis file)
    integer :: cell_index
    ! TODO: do the same with cluster
    integer :: cluster_index
    double precision :: Z_nuclear !todo
    !double precision :: Vref      !todo  ! or move to reference cluster? yes
    !double precision :: RMTref    ! move to ref cluster
    ! double precision, dimension(3) :: position !todo
    ! position, etc...
    ! lmax ???

    type (PotentialData) :: potential
    type (AtomicCoreData) :: core

    type (CellData), pointer :: cell_ptr => null()
    type (RefClusterData), pointer :: cluster_ptr => null()

  end type

CONTAINS

  !----------------------------------------------------------------------------
  subroutine createBasisAtom(atom, atom_index, Z_nuclear, lpot, nspin, irmind, irmd)
    use PotentialData_mod
    use AtomicCoreData_mod
    implicit none

    type (BasisAtom), intent(inout) :: atom
    integer, intent(in) :: lpot
    integer, intent(in) :: nspin
    integer, intent(in) :: irmind
    integer, intent(in) :: irmd  !< number of mesh points
    integer, intent(in) :: atom_index
    double precision, intent(in) :: Z_nuclear !< nuclear charge

    atom%atom_index = atom_index
    atom%cell_index = -1
    atom%cluster_index = -1
    atom%Z_nuclear = Z_nuclear

    call createPotentialData(atom%potential, lpot, nspin, irmind, irmd)

    ! does irmd have to be the same?
    call createAtomicCoreData(atom%core, Z_nuclear, irmd)

  end subroutine

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its cell data.
  !> This is done that way because different atoms can share the same
  !> cell data.
  subroutine associateBasisAtomCell(atom, cell_ptr)
    use CellData_mod
    implicit none
    type (BasisAtom), intent(inout) :: atom
    type (CellData), pointer :: cell_ptr

    atom%cell_ptr => cell_ptr
    atom%cell_index = cell_ptr%cell_index

  end subroutine

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its reference cluster.
  !> This is done that way because different atoms can share the same
  !> reference cluster data.
  subroutine associateBasisAtomRefCluster(atom, cluster_ptr)
    use RefClusterData_mod
    implicit none
    type (BasisAtom), intent(inout) :: atom
    type (RefClusterData), pointer :: cluster_ptr

    atom%cluster_ptr => cluster_ptr
    atom%cluster_index = cluster_ptr%cluster_index  ! TODO

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyBasisAtom(atom)
    use PotentialData_mod
    use AtomicCoreData_mod
    implicit none

    type (BasisAtom), intent(inout) :: atom

    nullify(atom%cell_ptr)
    nullify(atom%cluster_ptr)

    call destroyPotentialData(atom%potential)
    call destroyAtomicCoreData(atom%core)
  end subroutine

  !----------------------------------------------------------------------------
  !> Write to unformatted direct-access (DA) file.
  !> The format is compatible to the 'vpotnew' file
  !> A record number has to be specified - argument 'recnr'
  !> File has to be opened and closed by user. No checks
  subroutine writeBasisAtomPotentialDA(atom, fileunit, recnr)
    implicit none
    type (BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    write (fileunit, rec=recnr) atom%potential%VINS, &
                                atom%potential%VISP, &
                                atom%core%ECORE

  end subroutine

  !----------------------------------------------------------------------------
  !> Opens unformatted direct-access (DA) file.
  !> The format is compatible to the 'vpotnew' file
  subroutine openBasisAtomPotentialDAFile(atom, fileunit, filename)
    implicit none
    type (BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

    integer :: reclen

    inquire (iolength = reclen)   atom%potential%VINS, &
                                  atom%potential%VISP, &
                                  atom%core%ECORE

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')

  end subroutine

  !----------------------------------------------------------------------------
  !> Reads from unformatted direct-access (DA) file.
  !> The format is compatible to the 'vpotnew' file
  !> A record number has to be specified - argument 'recnr'
  !> File has to be opened and closed by user. No checks
  subroutine readBasisAtomPotentialDA(atom, fileunit, recnr)
    implicit none
    type (BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    read (fileunit, rec=recnr) atom%potential%VINS, &
                                atom%potential%VISP, &
                                atom%core%ECORE

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes unformatted direct-access (DA) potential-file.
  !> The format is compatible to the 'vpotnew' file
  subroutine closeBasisAtomPotentialDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit

    close(fileunit)

  end subroutine

  !----------------------------------------------------------------------------
  !> Write basis atom data to direct access file 'fileunit' at record 'recnr'
  subroutine writeBasisAtomDA(atom, fileunit, recnr)

    implicit none
    type (BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = 385306

    write (fileunit, rec=recnr) MAGIC_NUMBER, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%cluster_index, &
                                atom%Z_nuclear, &
                                atom%core%NCORE_atom, &
                                atom%core%LCORE_atom, &
                                atom%core%LCOREMAX, &
                                atom%core%ITITLE, &
                                MAGIC_NUMBER
  end subroutine

  !> Read basis atom data from direct access file 'fileunit' at record 'recnr'
  subroutine readBasisAtomDA(atom, fileunit, recnr)
    implicit none

    type (BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = 385306
    integer :: magic, magic2

    read  (fileunit, rec=recnr) magic, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%cluster_index, &
                                atom%Z_nuclear, &
                                atom%core%NCORE_atom, &
                                atom%core%LCORE_atom, &
                                atom%core%LCOREMAX, &
                                atom%core%ITITLE, &
                                magic2

    if (magic /= MAGIC_NUMBER .or. magic2 /= MAGIC_NUMBER) then
      write (*,*) "ERROR: Invalid basis atom data read. ", __FILE__, __LINE__
      STOP
    end if
  end subroutine

  !----------------------------------------------------------------------------
  !> Opens BasisAtom direct access file.
  subroutine openBasisAtomDAFile(atom, fileunit, filename)
    implicit none

    type (BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

    !------
    integer :: reclen

    integer, parameter :: MAGIC_NUMBER = 385306

    inquire (iolength = reclen) MAGIC_NUMBER, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%cluster_index, &
                                atom%Z_nuclear, &
                                atom%core%NCORE_atom, &
                                atom%core%LCORE_atom, &
                                atom%core%LCOREMAX, &
                                atom%core%ITITLE, &
                                MAGIC_NUMBER

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes BasisAtom direct access file.
  subroutine closeBasisAtomDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit

    close(fileunit)

  end subroutine


end module BasisAtom_mod
