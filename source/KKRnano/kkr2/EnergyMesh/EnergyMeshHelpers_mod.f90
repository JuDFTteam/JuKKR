module EnergyMeshHelpers_mod
  implicit none
  private
  public :: readEnergyMeshImpl, writeEnergyMeshImpl, updateEnergyMeshImpl, broadcastEnergyMeshImpl_com
  public :: readEnergyMeshImplSemi, writeEnergyMeshImplSemi, updateEnergyMeshImplSemi

  contains

  ! VALENCE CONTOUR ONLY!
  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh.0'
  subroutine readEnergyMeshImpl(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ)
    double precision, intent(out) :: E1
    double precision, intent(out) :: E2
    double precision, intent(out) :: EFERMI
    double complex, intent(out) :: EZ(:)
    integer, intent(out) :: IELAST
    integer, intent(out) :: NPNT1
    integer, intent(out) :: NPNT2
    integer, intent(out) :: NPNT3
    integer, intent(out) :: NPOL
    double precision, intent(out) :: TK
    double complex, intent(out) :: WEZ(:)

    open (67, file='energy_mesh.0', form='unformatted')
    read (67) IELAST,EZ,WEZ,E1,E2
    read (67) NPOL,TK,NPNT1,NPNT2,NPNT3

    !if ( NPOL==0 ) read(67) EFERMI
    read (67) EFERMI
    close(67)
  end subroutine

  ! VALENCE CONTOUR ONLY!
  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMeshImpl(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ)
    double precision, intent(in) :: E1
    double precision, intent(in) :: E2
    double precision, intent(in) :: EFERMI
    double complex, intent(in) :: EZ(:)
    integer, intent(in) :: IELAST
    integer, intent(in) :: NPNT1
    integer, intent(in) :: NPNT2
    integer, intent(in) :: NPNT3
    integer, intent(in) :: NPOL
    double precision, intent(in) :: TK
    double complex, intent(in) :: WEZ(:)

    open  (67, file='energy_mesh', form='unformatted')
    write (67) IELAST,EZ,WEZ,E1,E2
    write (67) NPOL,TK,NPNT1,NPNT2,NPNT3
    write (67) EFERMI
    close (67)
  end subroutine

  ! VALENCE CONTOUR ONLY!
  !------------------------------------------------------------------------------
  !> Update Energy mesh. Essentially a wrapper for EMESHT
  subroutine updateEnergyMeshImpl(EZ,WEZ,IELAST,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3)
  
    external :: EMESHT
    
    double complex, intent(out) :: EZ(:)
    double complex, intent(out) :: WEZ(:)
    integer, intent(inout) :: IELAST
    double precision, intent(in) :: E1
    double precision, intent(in) :: E2
    double precision, intent(in) :: TK
    integer, intent(in) :: NPOL
    integer, intent(in) :: NPNT1
    integer, intent(in) :: NPNT2
    integer, intent(in) :: NPNT3
    
    integer :: IE, IEMXD
    double precision :: PI
    
    PI = 4.0D0*ATAN(1.0D0)
    IEMXD = IELAST

    ! --> update energy contour
    call EMESHT(EZ,WEZ,IELAST,E1,E2,E2,TK, NPOL,NPNT1,NPNT2,NPNT3,IEMXD)
    ! IELAST will be overwritten

    do IE = 1,IELAST
      WEZ(IE) = -2.D0/PI*WEZ(IE)
    end do
    
  end subroutine


  ! VALENCE CONTOUR ONLY!
  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from rank 'BCRANK' to all other ranks
  subroutine broadcastEnergyMeshImpl_com(ACTVCOMM, BCRANK, E1, E2, EZ, IEMXD, WEZ)
    integer, intent(in) :: ACTVCOMM
    integer, intent(in) :: BCRANK
    double precision, intent(inout) :: E1
    double precision, intent(inout) :: E2
    double complex, intent(inout) :: EZ(:)
    integer, intent(in) :: IEMXD
    double complex, intent(inout) :: WEZ(:)

    include 'mpif.h'
    integer :: IERR

    call MPI_BCAST(EZ,IEMXD,MPI_DOUBLE_COMPLEX, BCRANK,ACTVCOMM,IERR)
    call MPI_BCAST(WEZ,IEMXD,MPI_DOUBLE_COMPLEX, BCRANK,ACTVCOMM,IERR)
    call MPI_BCAST(E1,1,MPI_DOUBLE_PRECISION, BCRANK,ACTVCOMM,IERR)
    call MPI_BCAST(E2,1,MPI_DOUBLE_PRECISION, BCRANK,ACTVCOMM,IERR)
  end subroutine


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> read energy mesh data from file 'energy_mesh.0'
  subroutine readEnergyMeshImplSemi(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, &
                   TK, WEZ, EBOTSEMI, EMUSEMI, FSEMICORE, IESEMICORE, N1SEMI, N2SEMI, N3SEMI)
    ! valence contour parameters
    double precision, intent(out) :: E1
    double precision, intent(out) :: E2
    double precision, intent(out) :: EFERMI
    double complex, intent(out) :: EZ(:)
    integer, intent(out) :: IELAST
    integer, intent(out) :: NPNT1
    integer, intent(out) :: NPNT2
    integer, intent(out) :: NPNT3
    integer, intent(out) :: NPOL
    double precision, intent(out) :: TK
    double complex, intent(out) :: WEZ(:)

    ! semicore contour parameters
    double precision, intent(out) :: EBOTSEMI
    double precision, intent(out) :: EMUSEMI
    double precision, intent(out) :: FSEMICORE
    integer, intent(out) :: IESEMICORE
    integer, intent(out) :: N1SEMI
    integer, intent(out) :: N2SEMI
    integer, intent(out) :: N3SEMI

    open (67,file='energy_mesh.0',form='unformatted')
    read (67) IELAST,EZ,WEZ,E1,E2
    read (67) NPOL,TK,NPNT1,NPNT2,NPNT3
    read (67) EFERMI
    read (67) IESEMICORE,FSEMICORE,EBOTSEMI
    read (67) EMUSEMI
    read (67) N1SEMI,N2SEMI,N3SEMI
    close(67)
  end subroutine

  ! VALENCE AND SEMICORE CONTOUR!
  !----------------------------------------------------------------------------
  !> write energy mesh data to file 'energy_mesh'
  subroutine writeEnergyMeshImplSemi(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, &
                    TK, WEZ, EBOTSEMI, EMUSEMI, FSEMICORE, IESEMICORE, N1SEMI, N2SEMI, N3SEMI)
    ! valence contour parameters
    double precision, intent(in) :: E1
    double precision, intent(in) :: E2
    double precision, intent(in) :: EFERMI
    double complex, intent(in) :: EZ(:)
    integer, intent(in) :: IELAST
    integer, intent(in) :: NPNT1
    integer, intent(in) :: NPNT2
    integer, intent(in) :: NPNT3
    integer, intent(in) :: NPOL
    double precision, intent(in) :: TK
    double complex, intent(in) :: WEZ(:)

    ! semicore contour parameters
    double precision, intent(in) :: EBOTSEMI
    double precision, intent(in) :: EMUSEMI
    double precision, intent(in) :: FSEMICORE
    integer, intent(in) :: IESEMICORE
    integer, intent(in) :: N1SEMI
    integer, intent(in) :: N2SEMI
    integer, intent(in) :: N3SEMI

    open  (67,file='energy_mesh',form='unformatted')
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
  !> Update Energy mesh. Essentially a wrapper for EPATHTB
  subroutine updateEnergyMeshImplSemi(EZ,WEZ,IELAST,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3, &
                                      EBOTSEMI,EMUSEMI,IESEMICORE,FSEMICORE,N1SEMI,N2SEMI,N3SEMI)
    ! valence contour parameters
    double complex, intent(out) :: EZ(:)
    double complex, intent(out) :: WEZ(:)
    integer, intent(inout) :: IELAST
    double precision, intent(in) :: E1
    double precision, intent(in) :: E2
    double precision, intent(in) :: TK
    integer, intent(in) :: NPOL
    integer, intent(in) :: NPNT1
    integer, intent(in) :: NPNT2
    integer, intent(in) :: NPNT3
    
    ! semicore contour parameters
    double precision, intent(in) :: EBOTSEMI
    double precision, intent(in) :: EMUSEMI
    double precision, intent(in) :: FSEMICORE
    integer, intent(in) :: IESEMICORE
    integer, intent(in) :: N1SEMI
    integer, intent(in) :: N2SEMI
    integer, intent(in) :: N3SEMI

    double precision :: PI
    integer :: IEMXD
    integer :: IE

    PI = 4.0D0*ATAN(1.0D0)
    IEMXD = IELAST

    ! --> update energy contour
    call EPATHTB(EZ,WEZ,E2,IELAST,iesemicore,1, &
                 E1,E2,TK,npol,npnt1,npnt2,npnt3, &
                 ebotsemi,emusemi,tk,npol,n1semi,n2semi,n3semi, &
                 IEMXD)

    do IE = 1,IELAST
      WEZ(IE) = -2.D0/PI*WEZ(IE)
      IF ( IE.LE.IESEMICORE ) WEZ(IE) = WEZ(IE)*FSEMICORE
    end do

  end subroutine


  ! VALENCE AND SEMICORE CONTOUR!
  !---------------------------------------------------------------------------------
  !> Distribute EnergyMesh from rank 'BCRANK' to all other ranks
  subroutine broadcastEnergyMeshImplSemi_com(ACTVCOMM, BCRANK, E1, E2, EZ, IEMXD, WEZ, EBOTSEMI, EMUSEMI)
    ! valence contour parameters
    integer, intent(in) :: ACTVCOMM
    integer, intent(in) :: BCRANK
    double precision, intent(inout) :: E1
    double precision, intent(inout) :: E2
    double complex, intent(inout) :: EZ(:)
    integer, intent(in) :: IEMXD
    double complex, intent(inout) :: WEZ(:)

    ! valence contour parameters
    double precision, intent(inout) :: EBOTSEMI
    double precision, intent(inout) :: EMUSEMI
    
    include 'mpif.h'
    integer :: IERR

    call MPI_BCAST(EZ,IEMXD,MPI_DOUBLE_COMPLEX, BCRANK,ACTVCOMM,IERR)
    call MPI_BCAST(WEZ,IEMXD,MPI_DOUBLE_COMPLEX, BCRANK,ACTVCOMM,IERR)
    call MPI_BCAST(E1,1,MPI_DOUBLE_PRECISION, BCRANK,ACTVCOMM,IERR)
    call MPI_BCAST(E2,1,MPI_DOUBLE_PRECISION, BCRANK,ACTVCOMM,IERR)
    call MPI_BCAST(EBOTSEMI,1,MPI_DOUBLE_PRECISION, BCRANK,ACTVCOMM,IERR)
    call MPI_BCAST(EMUSEMI,1,MPI_DOUBLE_PRECISION, BCRANK,ACTVCOMM,IERR)
  end subroutine

end module
