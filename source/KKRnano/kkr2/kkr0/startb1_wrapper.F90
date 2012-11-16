!-------------------------------------------------------------------------------
!> A wrapper for the subroutine STARTB1
subroutine STARTB1_wrapper(IFILE,IPF,IPFE,IPE,KHFELD,NBEG,NEND,RMTNEW,RMT, &
                           ITITLE,HFIELD,IMT,IRC,VCONST,IRNS,LPOT,NSPIN,IRMIN, &
                           NTCELL,IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP, &
                           EFERMI,VBC,RWS,LCORE,NCORE,DRDI,R,ZAT,A,B, IRWS,INIPOL, &
                           IINFO,IPAND,IRID,NFUND,IRMD,NCELLD,NAEZD,IRNSD)

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
  DOUBLE PRECISION, dimension(*) :: A
  DOUBLE PRECISION, dimension(*) :: B
  DOUBLE PRECISION, dimension(IRMD,*) :: DRDI
  DOUBLE PRECISION, dimension(IRMD,*) :: R
  DOUBLE PRECISION, dimension(*) :: RMT
  DOUBLE PRECISION, dimension(*) :: RMTNEW
  DOUBLE PRECISION, dimension(*) :: RWS
  DOUBLE PRECISION, dimension(IRID,NFUND,NCELLD) :: THETAS
  DOUBLE PRECISION, dimension(*) :: ZAT
  INTEGER, dimension((2*LPOT+1)**2,NAEZD) :: IFUNM
  INTEGER, dimension(*) :: IMT
  INTEGER, dimension(*) :: INIPOL
  INTEGER, dimension(*) :: IPAN
  INTEGER, dimension(*) :: IRC
  INTEGER, dimension(0:IPAND,*) :: IRCUT
  INTEGER, dimension(*) :: IRMIN
  INTEGER, dimension(*) :: IRNS
  INTEGER, dimension(*) :: IRWS
  INTEGER, dimension(20,*) :: ITITLE
  INTEGER, dimension(20,*) :: LCORE
  INTEGER, dimension(NFUND,NAEZD) :: LLMSP
  INTEGER, dimension((2*LPOT+1)**2,NAEZD) :: LMSP
  INTEGER, dimension(*) :: NCORE
  INTEGER, dimension(NCELLD) :: NFU
  INTEGER, dimension(*) :: NTCELL

  ! --- locals ---
  type (CellData) :: cell
  type (RadialMeshData) :: meshdata
  integer :: ii

  call createCellData(cell, irmd, ipand, irid, (2*LPOT+1)**2, nfund)
  call createRadialMeshData(meshdata, irmd, ipand)

  call STARTB1(IFILE,IPF,IPFE,IPE,KHFELD,NBEG,NEND,RMTNEW,RMT, &
               ITITLE,HFIELD,IMT,IRC,VCONST,IRNS,LPOT,NSPIN,IRMIN, &
               NTCELL,IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP, &
               EFERMI,VBC,RWS,LCORE,NCORE,DRDI,R,ZAT,A,B,IRWS,INIPOL,IINFO,IPAND,IRID,NFUND,IRMD,NCELLD,NAEZD,IRNSD)

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

