module ShapeIntegration_mod
  implicit none

  contains

!------------------------------------------------------------------------------
!> Performs the shape-function integration.
!> @brief
!> Needs information from module/common block replacement "tetrahedra_common"
!> Needs information from module/common block replacement "angles_common"
!> @param[in]  LMAX calculate shape function up to lmax
!> @param[in]  NFACE number of cell faces
!> @param[out] MESHN number of radial mesh points
!> @param[in]  XRN radial mesh points
!> @param[in]  DLT step size for Gauss-Legendre integration
!> @param[out] THETAS_S on output it contains the non-zero shape-functions
!!                      dimension MESHND x IBMAXD
!> @param[out] LMIFUN_S array with lm-indeices for which the shapefunction is non-zero
!> @param[out] NFUN number of non-zero shape-functions
!> @param[in]  MESHND
!> @param[in]  IBMAXD

SUBROUTINE shapeIntegration(LMAX, NFACE, MESHN, XRN, DLT, THETAS_S, LMIFUN_S, NFUN, MESHND, IBMAXD)

  use shape_constants_mod, only: PI, LMAXD1, ISUMD, ICD, ICED
  use tetrahedra_common
  use angles_common
  use ShapeIntegrationHelpers_mod
  implicit none

  integer :: LMAX
  integer :: NFACE
  integer :: MESHN
  REAL*8 ::  XRN(MESHND)
  real*8 :: DLT
  real*8 ::  THETAS_S(MESHND,IBMAXD)
  integer :: LMIFUN_S(IBMAXD)
  integer, intent(in) :: MESHND
  integer, intent(in) :: IBMAXD
  integer, intent(out) :: NFUN

  REAL*8 ::  DMATL(ISUMD)

  REAL*8 ::     CL(ICD)
  REAL*8 ::     C(ICED)

  REAL*8 ::  R
  REAL*8 ::  RAP
  REAL*8 ::  RDOWN

  REAL*8 ::  ARG1
  REAL*8 ::  ARG2

  integer :: IVTOT
  integer :: IFACE
  integer :: NTET
  integer :: ITET
  integer :: I
  integer :: M
  integer :: IC
  integer :: ICE
  integer :: IB
  integer :: ISU
  integer :: L
  integer :: MP
  integer :: K0
  integer :: K
  integer :: IS
  integer :: IMAX
  integer :: MO
  integer :: IP
  integer :: IPMAX
  integer :: icount

  integer :: IBM
  integer :: IBMAX
  integer :: N

  REAL*8, allocatable ::  RUPSQ(:) ! dim NVTOTD

  REAL*8 ::     S(-LMAXD1:LMAXD1,0:LMAXD1)
  REAL*8 ::    S1(-LMAXD1:LMAXD1,0:LMAXD1)
  REAL*8 ::    S2(-LMAXD1:LMAXD1,0:LMAXD1)
  REAL*8 ::    S3(-LMAXD1:LMAXD1,0:LMAXD1)
  REAL*8 ::    SUM(0:LMAXD1,2)
  REAL*8 ::    FK
  REAL*8 ::    FL
  REAL*8 ::    FPISQ

  ! local automatic arrays
  INTEGER ::   LOFM(IBMAXD)
  integer ::   MOFM(IBMAXD)
  REAL*8 ::    B(IBMAXD)
  integer ::   ISW(IBMAXD)
  ! ISW index array, value is 1 for lm for which shapefunction is non-zero
  ! otherwise the value is 0

  allocate(RUPSQ(size(RD))) ! dim: NVTOTD

  FPISQ=DSQRT(4.D0*PI)

  IBMAX=(LMAX+1)*(LMAX+1)

  THETAS_S = 0.0d0

  !.......................................................................
  !     E X P A N S I O N    C O E F F I C I E N T S
  !.......................................................................
  CALL CCOEF(LMAX,CL,C)
  IVTOT=0
  DO IFACE=1,NFACE
    NTET=NTT(IFACE)
    DO ITET=1,NTET
      IVTOT=IVTOT+1
      RUPSQ(IVTOT)=SQRT((RD(IVTOT)-R0(IFACE))*(RD(IVTOT)+R0(IFACE)))
    END DO
  END DO
  DO IBM=1,IBMAX
    ISW(IBM)=0
  END DO

  !===================== split ??? =======================================

  !.......................................................................
  !     L O O P    O V E R    R A D I A L    M E S H    P O I N T S
  !.......................................................................
  meshloop: DO N=1,MESHN
    R=XRN(N)
    DO IBM=1,IBMAX
      B(IBM)=0.D0
    END DO
    IVTOT=0
    !.......................................................................
    !     L O O P    O V E R    P Y R A M I D S
    !.......................................................................
py: DO IFACE=1,NFACE
      NTET=NTT(IFACE)

      IF(R > R0(IFACE))  GOTO 31
      IVTOT=IVTOT+NTET
      DO I=0,LMAX
        S(0,I) =0.D0
      END DO
      DO M=1,LMAX
        DO I=0,LMAX-M
          S(-M,I)=0.D0
          S( M,I)=0.D0
        END DO
      END DO
      GOTO 13
31    CONTINUE
      !IF(NEWSCH(IFACE) == 1) GOTO 35
      !IVTOT=IVTOT+NTET
      !GOTO 32
35    ARG1=R0(IFACE)/R
      RDOWN=SQRT((R-R0(IFACE))*(R+R0(IFACE)))
      DO I=0,LMAX
        S(0,I) =0.D0
      END DO
      DO M=1,LMAX
        DO I=0,LMAX-M
          S(-M,I)=0.D0
          S( M,I)=0.D0
        END DO
      END DO
      !.......................................................................
      !     L O O P     O V E R     T E T R A H E D R A
      !.......................................................................
      DO ITET=1,NTET
        IVTOT=IVTOT+1
        IF(R <= RD(IVTOT))      THEN
          CALL PINTG(FA(IVTOT),FB(IVTOT),DLT,S1,LMAX,ISIGNU(IVTOT), &
          ARG1,FD(IVTOT),0)
          DO I=0,LMAX
            S(0,I)=S(0,I)+S1(0,I)
          END DO
          DO M=1,LMAX
            DO I=0,LMAX-M
              S(-M,I)=S(-M,I)+S1(-M,I)
              S( M,I)=S( M,I)+S1( M,I)
            END DO
          END DO
        ELSE
          RAP =RUPSQ(IVTOT)/RDOWN
          ARG2=RUPSQ(IVTOT)/R0(IFACE)
          FK=FD(IVTOT)-ACOS(RAP)
          FL=FD(IVTOT)+ACOS(RAP)

          FK=DMAX1(FA(IVTOT),FK)
          FL=DMAX1(FA(IVTOT),FL)
          FK=DMIN1(FB(IVTOT),FK)
          FL=DMIN1(FB(IVTOT),FL)
          CALL PINTG(FA(IVTOT),FK,DLT,S1,LMAX,ISIGNU(IVTOT), &
          ARG1,FD(IVTOT),0)
          CALL PINTG(FK       ,FL,DLT,S2,LMAX,ISIGNU(IVTOT), &
          ARG2,FD(IVTOT),1)
          CALL PINTG(FL,FB(IVTOT),DLT,S3,LMAX,ISIGNU(IVTOT), &
          ARG1,FD(IVTOT),0)
          DO I=0,LMAX
            S(0,I)=S(0,I)+S1(0,I)+S2(0,I)+S3(0,I)
          END DO
          DO M=1,LMAX
            DO I=0,LMAX-M
              S(-M,I)=S(-M,I)+S1(-M,I)+S2(-M,I)+S3(-M,I)
              S( M,I)=S( M,I)+S1( M,I)+S2( M,I)+S3( M,I)
            END DO
          END DO
        END IF

      END DO  ! tetraeder loop

32    CONTINUE

      !.......................................................................
      !     I N T E G R A L   E X P A N S I O N        B A C K - R O T A T I O
      !.......................................................................

      IB=0
      IC=0
      ICE=0
      ! calculate transformation matrices for spherical harmonics
      CALL DREAL(LMAX,ALPHA(IFACE),BETA(IFACE),GAMMA(IFACE),DMATL,ISUMD,LMAXD1)

      ISU=0
      DO L=0,LMAX
        IB=IB+L+1
        DO MP=L,1,-1
          SUM(MP,1)=0.D0
          SUM(MP,2)=0.D0
          ICE=ICE+1
          K0=(L+MP+1)/2
          DO K=L,K0,-1
            IS=2*K-L-MP
            IC=IC+1
            SUM(MP,2)=SUM(MP,2)+CL(IC)*S(-MP,IS)
            SUM(MP,1)=SUM(MP,1)+CL(IC)*S( MP,IS)
          END DO
          SUM(MP,2)=SUM(MP,2)*C(ICE)
          SUM(MP,1)=SUM(MP,1)*C(ICE)
        END DO
        SUM(0,1)=0.D0
        ICE=ICE+1
        K0=(L+1)/2
        DO K=L,K0,-1
          IS=2*K-L
          IC=IC+1
          SUM(0,1)=SUM(0,1)+CL(IC)*S(0,IS)
        END DO
        SUM(0,1)=SUM(0,1)*C(ICE)
        IMAX=1
        M=0
  8     CONTINUE
        DO I=1,IMAX
          MO=(3-2*I)*M
          IBM=IB+MO
          LOFM(IBM)=L
          MOFM(IBM)=MO
          IPMAX=1
          MP=0
    16    CONTINUE
          DO IP=1,IPMAX
            ISU=ISU+1
            B(IBM)=B(IBM)+SUM(MP,IP)*DMATL(ISU)
          END DO
          IPMAX=2
          MP=MP+1
          IF(MP <= L) GOTO 16
        END DO
        IMAX=2
        M=M+1
        IF(M <= L) GOTO 8
        IB=IB+L
      END DO ! loop over L

  13 CONTINUE
    END DO py
    !.......................................................................
    !     D E F I N E S   A N D    S A V E S   S H A P E    F U N C T I O N ???
    !.......................................................................
    B(1)=FPISQ-B(1)/FPISQ
    DO IBM=2,IBMAX
      B(IBM)=-B(IBM)/FPISQ
    END DO
    DO IBM=1,IBMAX
      !     write(6,*) ibm,b(ibm)
      IF(ABS(B(IBM)) > 1.D-6) ISW(IBM)=1
      !IREC=(IBM-1)*MESHN+N
      !WRITE(11,REC=IREC) B(IBM)
      THETAS_S(N, IBM) = B(IBM)
    END DO
  END DO meshloop

  !Now rearrange Thetas_s array that it contains only non-zero shapefunctions
  !This is done "in-place"

  LMIFUN_S = 0

  icount = 1
  DO IBM=1,IBMAX
    IF (ISW(IBM) == 1) THEN

      LMIFUN_S(ICOUNT) = LOFM(IBM)*LOFM(IBM)+LOFM(IBM)+MOFM(IBM)+1

      IF (icount .ne. IBM) THEN
        DO N = 1, MESHN
          THETAS_S(N, icount) = THETAS_S(N, IBM)
        END DO
      END IF
      icount = icount + 1

    ELSE

      DO N = 1, MESHN
        THETAS_S(N, IBM) = 0.0d0
      END DO
    END IF
  END DO

  ! count non-zero shape functions
  NFUN=0
  DO IBM=1,IBMAX
    IF(ISW(IBM) == 1)  NFUN=NFUN+1
  END DO

END subroutine

end module ShapeIntegration_mod
