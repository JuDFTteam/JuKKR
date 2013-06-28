!-------------------------------------------------------------------------------
!> A wrapper for the subroutine STARTB1
subroutine STARTB1_wrapper(alat, LPOT,NSPIN, &
                           NTCELL, &
                           EFERMI,ZAT, RMTref, &
                           IPAND,IRID,NFUND,IRMD,NCELLD,NAEZD,IRNSD)

  use RadialMeshData_mod
  use CellData_mod
  implicit none

  ! Parameters

  ! Arguments
  double precision, intent(in) :: alat
  INTEGER, INTENT(IN) :: IPAND
  INTEGER, INTENT(IN) :: IRID
  INTEGER, INTENT(IN) :: NFUND
  INTEGER, INTENT(IN) :: IRMD
  INTEGER, INTENT(IN) :: NCELLD
  INTEGER, INTENT(IN) :: NAEZD
  INTEGER, INTENT(IN) :: IRNSD
  DOUBLE PRECISION :: EFERMI
  INTEGER :: LPOT
  INTEGER :: NSPIN
  DOUBLE PRECISION, dimension(*) :: ZAT
  double precision, dimension(naezd), intent(in) :: RMTref
  INTEGER, dimension(*) :: NTCELL

  DOUBLE PRECISION, parameter :: HFIELD = 0.0d0
  DOUBLE PRECISION, dimension(2), parameter :: VBC = (/ 0.0d0, 0.0d0 /)
  DOUBLE PRECISION, parameter :: VCONST = 0.0d0
  INTEGER, parameter :: IFILE = 13
  INTEGER, parameter :: IPE = 0
  INTEGER, parameter :: IPF = 6
  INTEGER, parameter :: IPFE = 9
  INTEGER, parameter :: KHFELD = 0
  integer :: max_reclen

  ! --- locals ---
  type (RadialMeshData) :: meshdata
  integer :: ii

  ! the following arrays serve as local dummies

  !     .. core states ..
  integer, dimension(:,:), allocatable :: ITITLE
  integer, dimension(:,:), allocatable :: LCORE
  integer, dimension(:),   allocatable :: NCORE

  double precision, dimension(:,:,:), allocatable :: THETAS
  integer, dimension(:,:), allocatable :: IFUNM
  integer, dimension(:),   allocatable :: IPAN
  integer, dimension(:,:), allocatable :: LLMSP
  integer, dimension(:,:), allocatable :: LMSP
  integer, dimension(:),   allocatable :: NFU

  double precision, dimension(:),   allocatable  :: A
  double precision, dimension(:),   allocatable  :: B
  double precision, dimension(:,:), allocatable  :: DRDI
  double precision, dimension(:,:), allocatable  :: R
  integer, dimension(:),  allocatable  :: IMT
  integer, dimension(:),  allocatable  :: IRC
  integer, dimension(:,:),allocatable  :: IRCUT
  integer, dimension(:),  allocatable  :: IRMIN
  integer, dimension(:),  allocatable  :: IRNS
  integer, dimension(:),  allocatable  :: IRWS
  double precision, dimension(:),   allocatable :: RWS
  double precision, dimension(:),   allocatable :: RMT

  double precision, dimension(:), allocatable :: RMTNEW
  integer, dimension(:), allocatable :: INIPOL

  allocate(ITITLE(20,naezd*nspin))
  allocate(LCORE(20,naezd*nspin))
  allocate(NCORE(naezd*nspin))

  allocate(THETAS(IRID,NFUND,NCELLD))
  allocate(IFUNM((2*LPOT+1)**2,NAEZD))
  allocate(IPAN(NAEZD))
  allocate(LLMSP(NFUND,NAEZD))
  allocate(LMSP((2*LPOT+1)**2,NAEZD))
  allocate(NFU(NAEZD))

  ! Radial mesh(es)
  allocate(A(NAEZD))
  allocate(B(NAEZD))
  allocate(DRDI(IRMD,NAEZD))
  allocate(R(IRMD,NAEZD))
  allocate(IMT(NAEZD))
  allocate(IRC(NAEZD))
  allocate(IRCUT(0:IPAND,NAEZD))
  allocate(IRMIN(NAEZD))
  allocate(IRNS(NAEZD))
  allocate(IRWS(NAEZD))
  allocate(RWS(NAEZD))
  allocate(RMT(NAEZD))

  allocate(RMTNEW(NAEZD))
  allocate(INIPOL(NAEZD))

  call createRadialMeshData(meshdata, irmd, ipand)

  NCORE = 0
  LCORE = 0
  ITITLE = 0

  open (19,file='shapefun',status='old',form='formatted')
  open (13,file='potential',status='old',form='formatted')

  call STARTB1(alat, IFILE,IPF,IPFE,IPE,KHFELD,1,naezd,RMTNEW,RMT, &
               ITITLE,HFIELD,IMT,IRC,VCONST,IRNS,LPOT,NSPIN,IRMIN, &
               NTCELL,IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP, &
               EFERMI,VBC,RWS,LCORE,NCORE,DRDI,R,ZAT,A,B,IRWS, &
               INIPOL,1,IPAND,IRID,NFUND,IRMD,NCELLD,NAEZD,IRNSD)

  close(13)
  close(19)

  call dochecks()

  call writeAtomData()

  call openRadialMeshDataDAFile(meshdata, 37 , "meshes")
  call openRadialMeshDataIndexDAFile(meshdata, 38, "meshes.idx")

  inquire(unit=37, recl = max_reclen)

  do ii = 1, naezd

    meshdata%R = R(:, ii)
    meshdata%DRDI = DRDI(:, ii)
    meshdata%IRCUT = IRCUT(:, ii)

    meshdata%A = A(ii)
    meshdata%B = B(ii)

    meshdata%IPAN = IPAN(ii)
    meshdata%IRC = IRC(ii)
    meshdata%IMT = IMT(ii)
    meshdata%IRNS = IRNS(ii)
    meshdata%IRWS = IRWS(ii)
    meshdata%IRMIN = IRMIN(ii)
    meshdata%RWS = RWS(ii)
    meshdata%RMT = RMT(ii)

    call writeRadialMeshDataDA(meshdata, 37, ii)
    call writeRadialMeshDataIndexDA(meshdata, 38, ii, max_reclen)

  end do


  call closeRadialMeshDataIndexDAFile(38)
  call closeRadialMeshDataDAFile(37)

  call destroyRadialMeshData(meshdata)


  CONTAINS

    subroutine writeAtomData()
      use BasisAtom_mod
      implicit none

      type (BasisAtom) :: atom
      integer :: ii, ispin, ipot
      integer :: max_reclen

      call createBasisAtom(atom, 1, lpot, nspin, (irmd-irnsd), irmd)  ! create dummy basis atom

      call openBasisAtomDAFile(atom, 37, 'atoms')
      call openBasisAtomPotentialIndexDAFile(atom, 38, 'vpotnew.idx')

      inquire (iolength = max_reclen) atom%potential%VINS, &
                                      atom%potential%VISP, &
                                      atom%core%ECORE

      do ii = 1, naezd
        atom%atom_index = ii
        atom%cell_index = NTCELL(ii)
        atom%cluster_index = -1 !TODO
        atom%Z_nuclear = ZAT(ii)
        atom%RMTref = RMTref(ii)

        atom%core%NCORE = 0
        atom%core%LCORE = 0
        atom%core%ITITLE = 0

        do ispin = 1, nspin
          ipot = NSPIN * (ii-1) + ispin
          atom%core%NCORE(ispin) = NCORE(ipot)
          atom%core%LCORE(:, ispin) = LCORE(:, ipot)
          atom%core%ITITLE(:, ispin) = ITITLE(:, ipot)
        enddo

        call writeBasisAtomDA(atom, 37, ii)
        call writeBasisAtomPotentialIndexDA(atom, 38, ii, max_reclen)

      enddo

      call closeBasisAtomPotentialIndexDAFile(38)
      call closeBasisAtomDAFile(37)

      call destroyBasisAtom(atom)

    end subroutine

!------------------------------------------------------------------------------
  subroutine dochecks()
    implicit none

    integer :: ierror
    integer :: I1, IE

!------------------- Do some tests on startb1 output -------------------
    ierror = 0
    do I1 = 1, NAEZD
      if ((IRMD - IRCUT(1, I1)) > IRID) then
        write(*,*) "Error: Not enough radial points for shape-function available."
        write(*,*) "Atom, points needed, points given ", I1, (IRMD - IRCUT(1, I1)), IRID
        ierror = 1
      end if
    end do

    do I1 = 1, NAEZD
      do IE = 0, IPAN(I1) - 1
        if (IE >= IPAND) then
          write(*,*) "Array IPAN contains bad number of panels for atom ", IE
          stop
        end if
        if ((IRCUT(IE+1,I1) -  IRCUT(IE,I1)) < 5) then
          write(*,*) "Error: Not enough points in panel. At least 5 points needed for REGSOL"
          write(*,*) "Atom, panel, points ", I1, IE, (IRCUT(IE+1,I1) -  IRCUT(IE,I1))
          ierror = 1
        end if
      end do
    end do
    if (ierror /= 0) stop
!-----------------------------------------------------------------------
  end subroutine

end subroutine

