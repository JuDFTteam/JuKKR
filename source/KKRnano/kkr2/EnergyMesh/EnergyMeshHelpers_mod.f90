module EnergyMeshHelpers_mod

contains

  ! VALENCE CONTOUR ONLY!
  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh.0'
  subroutine readEnergyMeshImpl(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ)
    implicit none
    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
    double complex :: EZ(:)
    integer :: IELAST
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    double precision :: TK
    double complex :: WEZ(:)

    open (67,file='energy_mesh.0',form='unformatted')
    read (67) IELAST,EZ,WEZ,E1,E2
    read (67) NPOL,TK,NPNT1,NPNT2,NPNT3

    !if ( NPOL==0 ) read(67) EFERMI
    read(67) EFERMI
    close (67)
  end subroutine

  ! VALENCE CONTOUR ONLY!
  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMeshImpl(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ)
    implicit none
    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
    double complex :: EZ(:)
    integer :: IELAST
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    double precision :: TK
    double complex :: WEZ(:)

    open (67,file='energy_mesh',form='unformatted')
    write (67) IELAST,EZ,WEZ,E1,E2
    write (67) NPOL,TK,NPNT1,NPNT2,NPNT3
    write (67) EFERMI
    close (67)
  end subroutine

  ! VALENCE CONTOUR ONLY!
  !------------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EMESHT
  subroutine updateEnergyMeshImpl(EZ,WEZ,IELAST,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3)
    implicit none
    double complex :: WEZ(:)
    double precision :: E1
    double precision :: E2
    double complex :: EZ(:)
    integer :: IE
    integer :: IELAST
    integer :: IEMXD
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    double precision :: PI
    double precision :: TK

    PI = 4.0D0*ATAN(1.0D0)
    IEMXD = IELAST

    ! --> update energy contour

    call EMESHT(EZ,WEZ,IELAST,E1,E2,E2,TK, &
    NPOL,NPNT1,NPNT2,NPNT3,IEMXD)

    do IE = 1,IELAST
      WEZ(IE) = -2.D0/PI*WEZ(IE)
    end do
  end subroutine


  ! VALENCE CONTOUR ONLY!
  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from rank 'BCRANK' to all other ranks
  subroutine broadcastEnergyMeshImpl_com(ACTVCOMM, BCRANK, E1, E2, EZ, IEMXD, WEZ)
    implicit none
    integer, intent(in) :: ACTVCOMM
    integer, intent(in) :: BCRANK
    double precision :: E1
    double precision :: E2
    double complex :: EZ(:)
    integer, intent(in) :: IEMXD
    double complex :: WEZ(:)

    !---------------
    integer :: IERR

    include 'mpif.h'

    call MPI_BCAST(EZ,IEMXD,MPI_DOUBLE_COMPLEX, &
    BCRANK,ACTVCOMM,IERR)

    call MPI_BCAST(WEZ,IEMXD,MPI_DOUBLE_COMPLEX, &
    BCRANK,ACTVCOMM,IERR)

    call MPI_BCAST(E1,1,MPI_DOUBLE_PRECISION, &
    BCRANK,ACTVCOMM,IERR)

    call MPI_BCAST(E2,1,MPI_DOUBLE_PRECISION, &
    BCRANK,ACTVCOMM,IERR)
  end subroutine


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh.0'
  subroutine readEnergyMeshImplSemi(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, &
                                    TK, WEZ, EBOTSEMI, EMUSEMI, FSEMICORE, IESEMICORE, N1SEMI, N2SEMI, &
                                    N3SEMI)
    implicit none

    ! valence contour parameters
    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
    double complex :: EZ(:)
    integer :: IELAST
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    double precision :: TK
    double complex :: WEZ(:)

    ! semicore contour parameters
    double precision :: EBOTSEMI
    double precision :: EMUSEMI
    double precision :: FSEMICORE
    integer :: IESEMICORE
    integer :: N1SEMI
    integer :: N2SEMI
    integer :: N3SEMI


    open (67,file='energy_mesh.0',form='unformatted')
    read (67) IELAST,EZ,WEZ,E1,E2
    read (67) NPOL,TK,NPNT1,NPNT2,NPNT3
    read (67) EFERMI
    read (67) IESEMICORE,FSEMICORE,EBOTSEMI
    read (67) EMUSEMI
    read (67) N1SEMI,N2SEMI,N3SEMI

    close (67)
  end subroutine

  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMeshImplSemi(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, &
                                    TK, WEZ, EBOTSEMI, EMUSEMI, FSEMICORE, IESEMICORE, N1SEMI, N2SEMI, &
                                    N3SEMI)
    implicit none

    ! valence contour parameters
    double precision :: E1
    double precision :: E2
    double precision :: EFERMI
    double complex :: EZ(:)
    integer :: IELAST
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    double precision :: TK
    double complex :: WEZ(:)

    ! semicore contour parameters
    double precision :: EBOTSEMI
    double precision :: EMUSEMI
    double precision :: FSEMICORE
    integer :: IESEMICORE
    integer :: N1SEMI
    integer :: N2SEMI
    integer :: N3SEMI

    open (67,file='energy_mesh',form='unformatted')
    write (67) IELAST,EZ,WEZ,E1,E2
    write (67) NPOL,TK,NPNT1,NPNT2,NPNT3
    write (67) EFERMI
    write (67) IESEMICORE,FSEMICORE,EBOTSEMI
    write (67) EMUSEMI
    write (67) N1SEMI,N2SEMI,N3SEMI

    close (67)
  end subroutine

  ! VALENCE AND SEMICORE CONTOUR!
  !------------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EMESHT
  subroutine updateEnergyMeshImplSemi(EZ,WEZ,IELAST,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3, &
                                      EBOTSEMI,EMUSEMI,IESEMICORE,FSEMICORE,N1SEMI,N2SEMI,N3SEMI)
    implicit none

    ! valence contour parameters
    double complex :: WEZ(:)
    double precision :: E1
    double precision :: E2
    double complex :: EZ(:)
    integer :: IE
    integer :: IELAST
    integer :: IEMXD
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    double precision :: PI
    double precision :: TK

    ! semicore contour parameters
    double precision :: EBOTSEMI
    double precision :: EMUSEMI
    double precision :: FSEMICORE
    integer :: IESEMICORE
    integer :: N1SEMI
    integer :: N2SEMI
    integer :: N3SEMI


    integer :: I

    PI = 4.0D0*ATAN(1.0D0)
    IEMXD = IELAST

    ! --> update energy contour

    call EPATHTB(EZ,WEZ,E2,IELAST,iesemicore,1, &
                 E1,E2,TK,npol,npnt1,npnt2,npnt3, &
                 ebotsemi,emusemi,tk,npol,n1semi,n2semi,n3semi, &
                 IEMXD)

    write(*,*) 'FSEMICORE=' , FSEMICORE
    do IE = 1,IELAST
      WEZ(IE) = -2.D0/PI*WEZ(IE)
      IF ( IE.LE.IESEMICORE ) WEZ(IE) = WEZ(IE)*FSEMICORE
    end do

  do I = 1, IELAST
  end do

  end subroutine


  ! VALENCE AND SEMICORE CONTOUR!
  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from rank 'BCRANK' to all other ranks
  subroutine broadcastEnergyMeshImplSemi_com(ACTVCOMM, BCRANK, E1, E2, EZ, IEMXD, WEZ, EBOTSEMI, EMUSEMI)
    implicit none

    ! valence contour parameters
    integer, intent(in) :: ACTVCOMM
    integer, intent(in) :: BCRANK
    double precision :: E1
    double precision :: E2
    double complex :: EZ(:)
    integer, intent(in) :: IEMXD
    double complex :: WEZ(:)

    ! valence contour parameters
    double precision :: EBOTSEMI
    double precision :: EMUSEMI


    !---------------
    integer :: IERR

    include 'mpif.h'

    call MPI_BCAST(EZ,IEMXD,MPI_DOUBLE_COMPLEX, &
    BCRANK,ACTVCOMM,IERR)

    call MPI_BCAST(WEZ,IEMXD,MPI_DOUBLE_COMPLEX, &
    BCRANK,ACTVCOMM,IERR)

    call MPI_BCAST(E1,1,MPI_DOUBLE_PRECISION, &
    BCRANK,ACTVCOMM,IERR)

    call MPI_BCAST(E2,1,MPI_DOUBLE_PRECISION, &
    BCRANK,ACTVCOMM,IERR)

    call MPI_BCAST(EBOTSEMI,1,MPI_DOUBLE_PRECISION, &
    BCRANK,ACTVCOMM,IERR)

    call MPI_BCAST(EMUSEMI,1,MPI_DOUBLE_PRECISION, &
    BCRANK,ACTVCOMM,IERR)
  end subroutine





end module
