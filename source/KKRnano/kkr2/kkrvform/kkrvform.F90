!> Converts unformatted potential file to the JM formatted potential file format
!> NOTE: VBC entries are just dummy values!!! RMT, RMTNEW also?

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

program kkrvform
  use DimParams_mod
  use BasisAtom_mod
  use RadialMeshData_mod
  implicit none

  type (DimParams)      :: dims
  type (RadialMeshData) :: mesh
  type (BasisAtom)      :: atomdata
  integer :: iatom

  double precision :: alat = 0.0d0
  double precision :: efermi = 9999.0d0
  double precision :: getFermiEnergy
  double precision :: vbc(2) = 0.0d0
  integer :: kxc = 2

  efermi = getFermiEnergy()
  call getStuffFromInputCard(alat, kxc)

  call createDimParams(dims) ! read dim. parameters from 'inp0.unf'

  call createBasisAtom(atomdata, iatom, dims%lpot, &
                       dims%nspind, dims%irmind, dims%irmd)
  call createRadialMeshData(mesh, dims%irmd, dims%ipand)

  do iatom = 1, dims%naez

    write(*,*) "Writing potential ", iatom

    call openBasisAtomDAFile(atomdata, 37, "atoms")
    call readBasisAtomDA(atomdata, 37, iatom)
    call closeBasisAtomDAFile(37)
    CHECKASSERT( atomdata%atom_index == iatom )

    call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew")
    call readBasisAtomPotentialDA(atomdata, 37, iatom)
    call closeBasisAtomPotentialDAFile(37)

    call openRadialMeshDataDAFile(mesh, 37 , "meshes")
    call readRadialMeshDataDA(mesh, 37, iatom)
    call closeRadialMeshDataDAFile(37)

    call associateBasisAtomMesh(atomdata, mesh)

    call writeFormattedPotential(Efermi, ALAT, VBC, KXC, atomdata)

  end do



end program kkrvform


!------------------------------------------------------------------------------
!> Wraps writeFormattedPotentialImpl
subroutine writeFormattedPotential(Efermi, ALAT, VBC, KXC, atomdata)
  use BasisAtom_mod
  use RadialMeshData_mod
  implicit none

  double precision, intent(in) :: Efermi
  double precision, intent(in) :: VBC(2)
  integer, intent(in) :: KXC
  double precision, intent(in) :: ALAT
  type (BasisAtom), intent(in) :: atomdata

  !-------- locals
  integer :: nspind
  integer :: irnsd
  type (RadialMeshData), pointer :: mesh

  nspind = atomdata%nspin

  mesh => atomdata%mesh_ptr

  CHECKASSERT( associated(atomdata%mesh_ptr) )

  irnsd = atomdata%potential%irmd - atomdata%potential%irmind

  CHECKASSERT( atomdata%potential%irmd == mesh%irmd )

  call writeFormattedPotentialImpl(Efermi,VBC,NSPIND, &
            KXC,atomdata%potential%LPOT,mesh%A,mesh%B,mesh%IRC, &
            atomdata%potential%VINS,atomdata%potential%VISP,mesh%DRDI,mesh%IRNS,mesh%R,mesh%RWS,mesh%RMT,ALAT, &
            atomdata%core%ECORE,atomdata%core%LCORE(:,1:NSPIND),atomdata%core%NCORE(1:NSPIND),atomdata%Z_nuclear,atomdata%core%ITITLE(:,1:NSPIND), &
            atomdata%atom_index, mesh%irmd, irnsd)

end subroutine

!----------------------------------------------------------------------------
!> get Fermi energy from file 'energy_mesh'
  double precision function getFermiEnergy()
    implicit none
    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
    double complex, allocatable :: EZ(:)
    integer :: IELAST
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    double precision :: TK
    double complex, allocatable :: WEZ(:)

    open (67,file='energy_mesh',form='unformatted')
    read (67) IELAST
    close (67)

    allocate(WEZ(ielast))
    allocate(EZ(ielast))

    open (67,file='energy_mesh',form='unformatted')
    read (67) IELAST,EZ,WEZ,E1,E2
    read (67) NPOL,TK,NPNT1,NPNT2,NPNT3

    read (67) EFERMI
    close (67)

    getFermiEnergy = EFERMI
  end function

!---------------------------------------------------------------------------
  subroutine getStuffFromInputcard(alat, kxc)
    implicit none

    double precision, intent(out) :: alat
    integer, intent(out) :: kxc

    character(len=80) :: uio
    integer :: ier

    CALL IoInput('KEXCOR    ',UIO,1,7,IER)
    READ (UNIT=UIO,FMT=*) kxc
    CALL IoInput('ALATBASIS ',UIO,1,7,IER)
    READ (UNIT=UIO,FMT=*) ALAT

  end subroutine
