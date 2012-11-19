!-------------------------------------------------------------------------------
!> A wrapper for the subroutine STARTB1
subroutine STARTB1_wrapper(IFILE,IPF,IPFE,IPE,KHFELD, &
                           ITITLE,HFIELD,VCONST,LPOT,NSPIN, &
                           NTCELL, &
                           EFERMI,VBC,LCORE,NCORE,ZAT, &
                           IPAND,IRID,NFUND,IRMD,NCELLD,NAEZD,IRNSD)

  use RadialMeshData_mod
  use CellData_mod
  implicit none

  ! Parameters

  ! Arguments
  INTEGER, INTENT(IN) :: IPAND
  INTEGER, INTENT(IN) :: IRID
  INTEGER, INTENT(IN) :: NFUND
  INTEGER, INTENT(IN) :: IRMD
  INTEGER, INTENT(IN) :: NCELLD
  INTEGER, INTENT(IN) :: NAEZD
  INTEGER, INTENT(IN) :: IRNSD
  DOUBLE PRECISION :: EFERMI
  DOUBLE PRECISION :: HFIELD
  DOUBLE PRECISION, dimension(*) :: VBC
  DOUBLE PRECISION :: VCONST
  INTEGER :: IFILE
  INTEGER :: IINFO
  INTEGER :: IPE
  INTEGER :: IPF
  INTEGER :: IPFE
  INTEGER :: KHFELD
  INTEGER :: LPOT
  INTEGER :: NBEG
  INTEGER :: NEND
  INTEGER :: NSPIN

  DOUBLE PRECISION, dimension(*) :: ZAT

  INTEGER, dimension(20,*) :: ITITLE
  INTEGER, dimension(20,*) :: LCORE
  INTEGER, dimension(*) :: NCORE

  INTEGER, dimension(*) :: NTCELL

  ! --- locals ---
  type (CellData) :: cell
  type (RadialMeshData) :: meshdata
  integer :: ii

  double precision, dimension(:,:,:), allocatable :: THETAS !DEL
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
  double precision, dimension(:),   allocatable :: RWS !DEL
  double precision, dimension(:),   allocatable :: RMT !DEL


  double precision, dimension(:), allocatable :: RMTNEW !DEL
  integer, dimension(:), allocatable :: INIPOL !DEL

  allocate(THETAS(IRID,NFUND,NCELLD)) !DEL
  allocate(IFUNM((2*LPOT+1)**2,NAEZD))
  allocate(IPAN(NAEZD))
  allocate(LLMSP(NFUND,NAEZD))
  allocate(LMSP((2*LPOT+1)**2,NAEZD))
  allocate(NFU(NAEZD))

  ! Radial mesh(es)
  allocate(A(NAEZD)) !DEL
  allocate(B(NAEZD)) !DEL
  allocate(DRDI(IRMD,NAEZD)) !DEL
  allocate(R(IRMD,NAEZD)) !DEL
  allocate(IMT(NAEZD)) !DEL
  allocate(IRC(NAEZD)) !DEL
  allocate(IRCUT(0:IPAND,NAEZD)) !DEL
  allocate(IRMIN(NAEZD)) !DEL
  allocate(IRNS(NAEZD)) !DEL
  allocate(IRWS(NAEZD)) !DEL
  allocate(RWS(NAEZD)) !DEL
  allocate(RMT(NAEZD))

  allocate(RMTNEW(NAEZD))
  allocate(INIPOL(NAEZD))

  call createCellData(cell, irmd, ipand, irid, (2*LPOT+1)**2, nfund)
  call createRadialMeshData(meshdata, irmd, ipand)

  call STARTB1(IFILE,IPF,IPFE,IPE,KHFELD,1,naezd,RMTNEW,RMT, &
               ITITLE,HFIELD,IMT,IRC,VCONST,IRNS,LPOT,NSPIN,IRMIN, &
               NTCELL,IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP, &
               EFERMI,VBC,RWS,LCORE,NCORE,DRDI,R,ZAT,A,B,IRWS, &
               INIPOL,1,IPAND,IRID,NFUND,IRMD,NCELLD,NAEZD,IRNSD)

  call openCellDataDAFile(cell, 37 , "cells")

  do ii = 1, ncelld  ! NCELL or NCELLD ???
    cell%cell_index = ii
    cell%shdata%THETA(:,:) = THETAS(:,:,ii)
    cell%shdata%LLMSP = LLMSP(:,ii)
    cell%shdata%IFUNM = IFUNM(:,ii)
    cell%shdata%LMSP = LMSP(:,ii)

    !!! NFU !!!!!!!!!!!! ????????????????????????????????????????????????????????????

    call writeCellDataDA(cell, 37, ii)

  end do

  call closeCellDataDAFile(37)

  call openRadialMeshDataDAFile(meshdata, 37 , "meshes")

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

    call writeRadialMeshDataDA(meshdata, 37, ii)

  end do

  call closeRadialMeshDataDAFile(37)

  call destroyCellData(cell)
  call destroyRadialMeshData(meshdata)

end subroutine

