      MODULE MOD_BASISTRANSFORM
c **************************************************************************
c
c  This module transforms the wave function from the relativistic kappa-mu
c  representation in spin spherical harmonics representation to the real
c  l,m,s basis in real spherical harmonics representation
c
c  all subroutines have been taken from the Jülich-München code except for
c  "wavefunc_transpose"
c
c created July 2011 by Pascal Kordt
c
c calling:
c CALL SINGLE_TRANSFORM(LCUT,(LCUT+1)**2),RLL) 
c
c **************************************************************************

      CONTAINS



!       SUBROUTINE WAVEFUNC_TRANSFORM(LCUT,LMCUT,NRMAX,RLLVEC)
! 
!       IMPLICIT NONE
! 
!       INTEGER, INTENT(IN) :: LCUT, LMCUT, NRMAX
!       INTEGER             :: RNUM
!       COMPLEX*16          :: RLLVEC(LMCUT,LMCUT,NRMAX) 
! 
! !      DO RNUM=1,NRMAX
!       CALL SINGLE_TRANSFORM(LCUT,LMCUT,RLLVEC(:,:,50))
!       write(*,*) "alles in Budda"
!       stop
! c      WRITE(*,*) RLLVEC(1,1,1)
! !      END DO
! 
!       END SUBROUTINE

      SUBROUTINE RLL_TRANSFORM(RLL,lmaxatom,mode)
      IMPLICIT NONE
      DOUBLE COMPLEX :: RLL(:,:,:)
      INTEGER        :: lmaxatom
      CHARACTER(len=*) :: mode
      INTEGER        :: lmsize,ir,nrmax

      IF (MODE/='RLM>REL' .and. MODE/='REL>RLM') THEN
        STOP '[RLL_TRANSFORM] mode not known'
      END IF

      LMSIZE=2*(lmaxatom+1)**2

      IF ( 2*LMSIZE/=UBOUND(RLL,1) ) THEN 
        STOP '[BasisTransform] LMSIZE in RLL not equal'
      END IF
      IF (   LMSIZE/=UBOUND(RLL,2) ) THEN
        STOP '[BasisTransform] LMSIZE in RLL not equal'
      END IF
 
      NRMAX=UBOUND(RLL,3)

      DO IR=1,NRMAX
        call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,
     +                        RLL(1:lmsize,1:lmsize,ir),mode)
        call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,
     +               RLL(lmsize+1:2*lmsize,1:lmsize,ir),mode)
      END DO

      END SUBROUTINE RLL_TRANSFORM !(LCUT,LMCUT,RLL)

      SUBROUTINE VLL_TRANSFORM(VLL,lmaxatom)
      IMPLICIT NONE
      DOUBLE COMPLEX :: VLL(:,:,:)
      INTEGER        :: lmaxatom

      INTEGER        :: lmsize,ir,nrmax

      LMSIZE=2*(lmaxatom+1)**2

      IF ( 2*LMSIZE/=UBOUND(VLL,1) ) THEN 
        STOP '[BasisTransform] LMSIZE in VLL not equal'
      END IF
      IF (  2*LMSIZE/=UBOUND(VLL,2) ) THEN
        STOP '[BasisTransform] LMSIZE in VLL not equal'
      END IF
 
      NRMAX=UBOUND(VLL,3)

!       write(6264,'(50000E)') VLL(:,:,1)


      DO IR=1,NRMAX
        call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,
     +                        VLL(1:lmsize,1:lmsize,ir),'RLM>REL')
        call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,
     +               VLL(lmsize+1:2*lmsize,1:lmsize,ir),'RLM>REL')
        call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,
     +               VLL(1:lmsize,lmsize+1:2*lmsize,ir),'RLM>REL')
        call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2,
     +      VLL(lmsize+1:2*lmsize,lmsize+1:2*lmsize,ir),'RLM>REL')
      END DO

!       write(6265,'(50000E)') VLL(:,:,1)
!       stop

      END SUBROUTINE VLL_TRANSFORM !(LCUT,LMCUT,RLL)


      SUBROUTINE JLK_TRANSFORM(RLL,lmaxatom,loflm)
      IMPLICIT NONE
      DOUBLE COMPLEX,allocatable :: RLL(:,:)
      INTEGER        :: lmaxatom
      INTEGER        :: loflm(:)

      INTEGER        :: lmsize,ir,nrmax,LM
      DOUBLE COMPLEX,allocatable :: RLL2(:,:),RLL3(:,:)
      DOUBLE COMPLEX,allocatable :: MAT(:,:),MAT2(:,:)

!       LSIZE = 2*(lmaxatom+1)
      LMSIZE= 2*(lmaxatom+1)**2
!       print *,LMSIZE,ubound(LOFLM,1)
!       print *,'rll',ubound(rll)
      IF (2*LMSIZE/=ubound(LOFLM,1)) STOP '[WAVEFN2] loflm bound'
      NRMAX=ubound(RLL,2)
!       print *,'rll2',nrmax
      ALLOCATE ( RLL2(LMSIZE,NRMAX),  RLL3(LMSIZE,NRMAX)  )
      ALLOCATE ( MAT (LMSIZE,LMSIZE), MAT2(LMSIZE,LMSIZE) )
      MAT=(0.0D0,0.0D0)
      MAT2=(0.0D0,0.0D0)

      DO IR=1,NRMAX
        DO LM=1,LMSIZE
          MAT(LM,LM)=RLL(LOFLM(LM),IR)
        END DO !LM
        DO LM=LMSIZE+1,2*LMSIZE
          MAT2(LM-LMSIZE,LM-LMSIZE)=RLL(LOFLM(LM),IR)
        END DO !LM

!         write(888,'(50000E)') MAT
        call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2, MAT,'RLM>REL')
!         write(889,'(50000E)') MAT


        call SINGLE_TRANSFORM(lmaxatom,(lmaxatom+1)**2, MAT2,'RLM>REL')

!         DO LM=1,LMSIZE
!           write(6162,'(50000F)') MAT(:,LM)
!           write(6163,'(50000F)') MAT2(:,LM)
!         END DO !LM
!          stop

        DO LM=1,LMSIZE
          RLL2(LM,IR)=MAT(LM,LM)
        END DO !LM
        DO LM=1,LMSIZE
          RLL3(LM,IR)=MAT2(LM,LM)
        END DO !LM

        MAT=(0.0D0,0.0D0)
        MAT2=(0.0D0,0.0D0)


      END DO

      DEALLOCATE(RLL)
      ALLOCATE(RLL(2*LMSIZE,nrmax))

      DO IR=1,NRMAX
        DO LM=1,LMSIZE
          RLL(LM,IR)=RLL2(LM,IR)
        END DO !LM
        DO LM=LMSIZE+1,2*LMSIZE
          RLL(LM,IR)=RLL3(LM-LMSIZE,IR)
        END DO !LM
      END DO !IR


!       print *, ubound(rll)
      END SUBROUTINE JLK_TRANSFORM!(RLL,lmaxatom,loflm)




      SUBROUTINE SINGLE_TRANSFORM(LCUT,LMCUT,RLL,MODE)

      IMPLICIT NONE

      INTEGER,INTENT(IN) ::  LCUT,LMCUT

C local variables
      INTEGER            :: NLMAX,NKMMAX,NMUEMAX,
     +                      NKMPMAX,NKMAX,LINMAX 

      COMPLEX*16         :: RC(2*LMCUT,2*LMCUT),
     +                      CREL(2*LMCUT,2*LMCUT),
     +                      RREL(2*LMCUT,2*LMCUT),
     +                      RLL_REL(2*LMCUT,2*LMCUT),
     +                      RLL_NR(2*LMCUT,2*LMCUT),
     +                      RLL(2*LMCUT,2*LMCUT)
      CHARACTER(LEN=7)  :: MODE

!       IF (.not. present(mode) ) mode='REL>RLM'
C  prepare the transformation
      NLMAX     = LCUT+1 
      NKMMAX    = 2*NLMAX**2
      NMUEMAX   = 2*NLMAX
      NKMPMAX   = (NKMMAX+2*NLMAX)
      NKMAX     = 2*NLMAX-1
      LINMAX    = (2*NLMAX*(2*NLMAX-1))

      CALL DRVBASTRANS(RC,CREL,RREL,
     &                       NLMAX,NKMMAX,NMUEMAX,NKMPMAX,NKMAX,LINMAX)


      RLL_REL=0.0d0
      RLL_NR=0.0d0



C Relativistic kappa-mu (REL, spin spherical harmonics)
C to real l-m basis (RLM, real spherical hamonics)

      CALL CHANGEREP(RLL(1,1),mode,RLL_NR(1,1),
     +           NKMMAX,NKMMAX,RC,CREL,RREL,'T_Real',0)

      RLL(:,:) = RLL_NR(:,:)

      END SUBROUTINE












      SUBROUTINE CHANGEREP(A,MODE,B,N,M,RC,CREL,RREL,TEXT,LTEXT)
C   ********************************************************************
C   *                                                                  *
C   *   change the representation of matrix A and store in B           *
C   *   according to MODE:                                             *
C   *                                                                  *
C   *   RLM>REL   non-relat. REAL spher. harm.  >   (kappa,mue)        *
C   *   REL>RLM   (kappa,mue)  > non-relat. REAL spher. harm.          *
C   *   CLM>REL   non-relat. CMPLX. spher. harm.  >   (kappa,mue)      *
C   *   REL>CLM   (kappa,mue)  > non-relat. CMPLX. spher. harm.        *
C   *   RLM>CLM   non-relat. REAL spher. harm.  >  CMPLX. spher. harm. *
C   *   CLM>RLM   non-relat. CMPLX. spher. harm.  >  REAL spher. harm. *
C   *                                                                  *
C   *   the non-relat. representations include the  spin index         *
C   *                                                                  *
C   *   for LTEXT > 0 the new matrix  B  is printed                    *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 C1,C0
      PARAMETER (C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER LTEXT,M,N
      CHARACTER*7 MODE
      CHARACTER*(*) TEXT
      COMPLEX*16 A(M,M),B(M,M),CREL(M,M),RC(M,M),RREL(M,M)
C
C Local variables
C
      INTEGER KEY
      COMPLEX*16 W1(M,M)
C
C*** End of declarations rewritten by SPAG
C


C---------------------- transform MAT from (kappa,mue) to REAL (l,ml,ms)
      IF ( MODE.EQ.'REL>RLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,RREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,RREL,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'RLM>REL' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,RREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,RREL,M,C0,B,M)
         KEY = 3
      ELSE IF ( MODE.EQ.'REL>CLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,CREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,CREL,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'CLM>REL' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,CREL,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,CREL,M,C0,B,M)
         KEY = 3
      ELSE IF ( MODE.EQ.'CLM>RLM' ) THEN
         CALL ZGEMM('N','N',N,N,N,C1,RC,M,A,M,C0,W1,M)
         CALL ZGEMM('N','C',N,N,N,C1,W1,M,RC,M,C0,B,M)
         KEY = 2
      ELSE IF ( MODE.EQ.'RLM>CLM' ) THEN
         CALL ZGEMM('C','N',N,N,N,C1,RC,M,A,M,C0,W1,M)
         CALL ZGEMM('N','N',N,N,N,C1,W1,M,RC,M,C0,B,M)
         KEY = 2
      ELSE
         WRITE (*,*) ' MODE = ',MODE
         STOP 'in <ROTATE>  MODE not allowed'
      END IF
C
C show the input matrix
C      IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,A,N,M,KEY,KEY,0,1D-8,6)
C show the output matrix
      IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-8,6)
C     IF ( LTEXT.GT.0 ) CALL CMATSTR(TEXT,LTEXT,B,N,M,KEY,KEY,0,1D-12,6)
      END SUBROUTINE










      SUBROUTINE BASTRMAT(LMAX,CGC,RC,CREL,RREL,NKMMAX,NKMPMAX)
C   ********************************************************************
C   *                                                                  *
C   *    INITIALIZE TRANSFORMATION MATRIX THAT TAKES MATRICES FROM     *
C   *    RELATIVISTIC  TO  REAL SPERICAL HARM.  REPRESENTATION         *
C   *                                                                  *
C   *    this is a special version of <STRSMAT> passing the            *
C   *    full BASis TRansformation MATrices  RC, CREL and RREL         *
C   *                                                                  *
C   * 13/01/98  HE                                                     *
C   ********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CI,C1,C0
      PARAMETER (CI=(0.0D0,1.0D0),C1=(1.0D0,0.0D0),C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER LMAX,NKMMAX,NKMPMAX
      REAL*8 CGC(NKMPMAX,2)
      COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX),
     &           RREL(NKMMAX,NKMMAX)
C
C Local variables
C
      INTEGER I,IKM,J,JP05,K,L,LM,LNR,M,MUEM05,MUEP05,NK,NKM,NLM
      REAL*8 W
C
C*** End of declarations rewritten by SPAG
C
      NK = 2*(LMAX+1) + 1
      NLM = (LMAX+1)**2
      NKM = 2*NLM
C     ===================================================
C     INDEXING:
C     IKM  = L*2*(J+1/2) + J + MUE + 1
C     LM   = L*(L+1)     +     M   + 1
C     ===================================================
C
C ----------------------------------------------------------------------
C CREL  transforms from  COMPLEX (L,M,S)  to  (KAP,MUE) - representation
C                 |LAM> = sum[LC] |LC> * CREL(LC,LAM)
C ----------------------------------------------------------------------
      CALL CINIT(NKMMAX*NKMMAX,CREL)
C
      LM = 0
      DO LNR = 0,LMAX
         DO M = -LNR,LNR
            LM = LM + 1
C
            IKM = 0
            DO K = 1,NK
               L = K/2
               IF ( 2*L.EQ.K ) THEN
                  JP05 = L
               ELSE
                  JP05 = L + 1
               END IF
C
               DO MUEM05 = -JP05,(JP05-1)
                  MUEP05 = MUEM05 + 1
                  IKM = IKM + 1
C
                  IF ( L.EQ.LNR ) THEN
                     IF ( MUEP05.EQ.M ) CREL(LM,IKM) = CGC(IKM,1)
                     IF ( MUEM05.EQ.M ) CREL(LM+NLM,IKM) = CGC(IKM,2)
                  END IF
C
               END DO
            END DO
C
         END DO
      END DO
C
C ----------------------------------------------------------------------
C    RC  transforms from  REAL to  COMPLEX (L,M,S) - representation
C                 |LC> = sum[LR] |LR> * RC(LR,LC)
C ----------------------------------------------------------------------
      CALL CINIT(NKMMAX*NKMMAX,RC)
C
      W = 1.0D0/SQRT(2.0D0)
C
      DO L = 0,LMAX
         DO M = -L,L
            I = L*(L+1) + M + 1
            J = L*(L+1) - M + 1
C
            IF ( M.LT.0 ) THEN
               RC(I,I) = -CI*W
               RC(J,I) = W
               RC(I+NLM,I+NLM) = -CI*W
               RC(J+NLM,I+NLM) = W
            END IF
            IF ( M.EQ.0 ) THEN
               RC(I,I) = C1
               RC(I+NLM,I+NLM) = C1
            END IF
            IF ( M.GT.0 ) THEN
               RC(I,I) = W*(-1.0D0)**M
               RC(J,I) = CI*W*(-1.0D0)**M
               RC(I+NLM,I+NLM) = W*(-1.0D0)**M
               RC(J+NLM,I+NLM) = CI*W*(-1.0D0)**M
            END IF
         END DO
      END DO
C
C ----------------------------------------------------------------------
C RREL  transforms from   REAL (L,M,S)  to  (KAP,MUE) - representation
C                 |LAM> = sum[LR] |LR> * RREL(LR,LAM)
C ----------------------------------------------------------------------
C
      CALL ZGEMM('N','N',NKM,NKM,NKM,C1,RC,NKMMAX,CREL,NKMMAX,C0,RREL,
     &           NKMMAX)
C
      END SUBROUTINE










      SUBROUTINE CALCCGC(LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX)
C   ********************************************************************
C   *                                                                  *
C   *   CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)                    *
C   *                                                                  *
C   *   IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE                  *
C   *   IKM  = L*2*(J+1/2) + J + MUE + 1                               *
C   *   IS= 1/2  SPIN DOWN/UP                                          *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER NKMAX,NKMPMAX,NMUEMAX
      REAL*8 CGC(NKMPMAX,2)
      INTEGER KAPTAB(NMUEMAX),LTAB(NMUEMAX),NMUETAB(NMUEMAX)
C
C Local variables
C
      INTEGER IKM,K,KAPPA,M
      REAL*8 J,L,MUE,TWOLP1
C
      IKM = 0
      DO K = 1,(NKMAX+1)
         L = LTAB(K)
         KAPPA = KAPTAB(K)
         J = ABS(KAPPA) - 0.5D0
         MUE = -J - 1.0D0
         TWOLP1 = 2.0D0*L + 1.0D0
C
         IF ( KAPPA.LT.0 ) THEN
C
C     J = L + 1/2
            DO M = 1,NMUETAB(K)
C
               MUE = MUE + 1.0D0
               IKM = IKM + 1
               CGC(IKM,1) = DSQRT((L-MUE+0.5D0)/TWOLP1)
               CGC(IKM,2) = DSQRT((L+MUE+0.5D0)/TWOLP1)
            END DO
         ELSE
C     J = L - 1/2
            DO M = 1,NMUETAB(K)
C
               MUE = MUE + 1.0D0
               IKM = IKM + 1
               CGC(IKM,1) = DSQRT((L+MUE+0.5D0)/TWOLP1)
               CGC(IKM,2) = -DSQRT((L-MUE+0.5D0)/TWOLP1)
C
            END DO
         END IF
C
C
      END DO
C
      END SUBROUTINE











      SUBROUTINE DRVBASTRANS(RC,CREL,RREL,
     &                       NLMAX,NKMMAX,NMUEMAX,NKMPMAX,NKMAX,LINMAX)
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LINMAX,NKMAX,NKMMAX,NKMPMAX,NLMAX,NMUEMAX
      COMPLEX*16 CREL(NKMMAX,NKMMAX),RC(NKMMAX,NKMMAX),
     &           RREL(NKMMAX,NKMMAX),SRREL(2,2,NKMMAX)
      INTEGER IRREL(2,2,NKMMAX),NRREL(2,NKMMAX)
C
C Local variables
C
      REAL*8 CGC(NKMPMAX,2)
      INTEGER I,IKM1LIN(LINMAX),IKM2LIN(LINMAX),IL,IMUE,IPRINT,
     &        KAPTAB(NMUEMAX),LTAB(NMUEMAX),MMAX,NMUETAB(NMUEMAX),
     &        NSOLLM(NLMAX,NMUEMAX)
C
C*** End of declarations rewritten by SPAG
C
      IF (NKMMAX.NE.2*NLMAX**2) 
     &     STOP ' Check NLMAX,NKMMAX in < DRVBASTRANS > '
      IF (NMUEMAX.NE.2*NLMAX) 
     &     STOP ' Check NLMAX,NMUEMAX in < DRVBASTRANS > '
      IF (NKMPMAX.NE.(NKMMAX+2*NLMAX)) 
     &     STOP ' Check NLMAX,NKMMAX,NKMPMAX in < DRVBASTRANS > '
      IF (NKMAX.NE.2*NLMAX-1) 
     &     STOP ' Check NLMAX,NKMAX in < DRVBASTRANS > '
      IF (LINMAX.NE.(2*NLMAX*(2*NLMAX-1)))
     &     STOP ' Check NLMAX,LINMAX in < DRVBASTRANS > '
C
      IPRINT = 0
C
      DO I = 1,NMUEMAX
         LTAB(I) = I/2
         IF ( 2*LTAB(I).EQ.I ) THEN
            KAPTAB(I) = LTAB(I)
         ELSE
            KAPTAB(I) = -LTAB(I) - 1
         END IF
         NMUETAB(I) = 2*ABS(KAPTAB(I))
      END DO
C
      DO IL = 1,NLMAX
         MMAX = 2*IL
         DO IMUE = 1,MMAX
            IF ( (IMUE.EQ.1) .OR. (IMUE.EQ.MMAX) ) THEN
               NSOLLM(IL,IMUE) = 1
            ELSE
               NSOLLM(IL,IMUE) = 2
            END IF
         END DO
      END DO
C
c      CALL IKMLIN(IPRINT,NSOLLM,IKM1LIN,IKM2LIN,NLMAX,NMUEMAX,LINMAX,
c     &            NLMAX)
C
      CALL CALCCGC(LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX)
C
C ---------------------------- now calculate the transformation matrices
C
C      CALL STRSMAT(NLMAX-1,CGC,SRREL,NRREL,IRREL,NKMMAX,NKMPMAX)
C
      CALL BASTRMAT(NLMAX-1,CGC,RC,CREL,RREL,NKMMAX,NKMPMAX)
C
      RETURN
      END SUBROUTINE
C







C ************************************************************************
      SUBROUTINE CINIT(N,A)
C ************************************************************************
c-----------------------------------------------------------------------
c     initialize the first n values of a complex array a with zero
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX A(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
      DOUBLE COMPLEX CZERO
      parameter(CZERO=(0.0d0,0.0d0))
C     ..
      DO 10 I = 1,N
        A(I) = CZERO
   10 CONTINUE
      END SUBROUTINE









      SUBROUTINE CMATSTR(STR,LSTR,A,N,M,MLIN,MCOL,IJQ,TOLP,K_FMT_FIL)
C   ********************************************************************
C   *                                                                  *
C   *   writes structure of COMPLEX   NxN   matrix   A                 *
C   *                                                                  *
C   *   M           is the actual array - size used for   A            *
C   *   MLIN/COL    MODE for line and column indexing                  *
C   *               0: plain, 1: (l,ml), 2: (l,ml,ms), 3: (kap,mue)    *
C   *   TOL         tolerance for difference                           *
C   *   IJQ         if IJQ > 1000    pick  IQ-JQ-block matrix          *
C   *               assuming  IJQ = IQ*1000 + JQ                       *
C   *               else: no IQ-JQ-indexing                            *
C   *   K_FMT_FIL   output channel                                     *
C   *               a negative sign suppresses table at the end        *
C   *                                                                  *
C   *   any changes should be done in RMATSTR as well !!!!!!!!!!!!!!!  *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT COMPLEX*16(A-H,O-Z)
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      COMPLEX*16 CI
      PARAMETER (CI=(0.0D0,1.0D0))
C
C Dummy arguments
C
      INTEGER IJQ,K_FMT_FIL,LSTR,M,MCOL,MLIN,N
      CHARACTER*(*) STR
      REAL*8 TOLP
      COMPLEX*16 A(M,M)
C
C Local variables
C
      COMPLEX*16 B(N,N),CA,CB,ARG,DTAB(0:N*N)
      CHARACTER CHAR
      LOGICAL SAME,SMALL
      CHARACTER*1 CTAB(0:N*N),VZ(-1:+1)
      DOUBLE PRECISION DBLE
      CHARACTER*150 FMT1,FMT2,FMT3,FMT4
      INTEGER I,I1,IC0,ID,IL,ILSEP(20),IPT(218),IQ,ISL,IW(M),J,
     &        J0,JP,JQ,K,L3,LF,MM,N1,N2,N3,NC,ND,NFIL,NK,NM,NM1,NM2,NM3,
     &        NNON0,NSL
      INTEGER ICHAR,ISIGN,NINT
      REAL*8 TOL
C
C*** End of declarations rewritten by SPAG
C
      DATA VZ/'-',' ',' '/
C
      SMALL(ARG) = ABS(ARG*TOL).LT.1.0D0
C
      SAME(CA,CB) = SMALL(1.0D0-CA/CB)
C
      NFIL = ABS(K_FMT_FIL)
C
      TOL = 1.0D0/TOLP
C
C----------------------------------------------- set block indices IQ JQ
C
      IF ( IJQ.GT.1000 ) THEN
         IQ = IJQ/1000
         JQ = IJQ - IQ*1000
         IF ( IQ*N.GT.M .OR. IQ*N.GT.M ) THEN
            WRITE (6,99002) IJQ,IQ,JQ,IQ*N,JQ*N,N,M
            RETURN
         END IF
      ELSE
         IQ = 1
         JQ = 1
      END IF
C
C----------------------------------------------------- copy matrix block
C
      J0 = N*(JQ-1)
      DO J = 1,N
         I1 = N*(IQ-1)+1
         JP = J0 + J
         CALL ZCOPY(N,A(I1,JP),1,B(1,J),1)
      END DO
C
C------------------------------------------------ set up character table
C
      NC = 0
      DO I = 1,26
         NC = NC + 1
         IPT(NC) = 62 + I
      END DO
      DO I = 1,8
         NC = NC + 1
         IPT(NC) = 96 + I
      END DO
      DO I = 10,26
         NC = NC + 1
         IPT(NC) = 96 + I
      END DO
      DO I = 191,218
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 35,38
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 40,42
         NC = NC + 1
         IPT(NC) = I
      END DO
      DO I = 91,93
         NC = NC + 1
         IPT(NC) = I
      END DO
C
C---------------------------------------------------------------- header
      IC0 = ICHAR('0')
      N3 = N/100
      N2 = N/10 - N3*10
      N1 = N - N2*10 - N3*100
C
      IF ( N.LE.18 ) THEN
         FMT1 = '(8X,I3,''|'','
         FMT2 = '( 9X,''--|'','
         FMT3 = '( 9X,'' #|'','
         FMT4 = '( 9X,''  |'','
      ELSE
         FMT1 = '(   I4,''|'','
         FMT2 = '( 2X,''--|'','
         FMT3 = '( 2X,'' #|'','
         FMT4 = '( 2X,''  |'','
      END IF
C
      LF = 11
      L3 = 11
      IF ( MCOL.EQ.0 ) THEN
         FMT1 = FMT1(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)//CHAR(IC0+N1)
     &          //'( 2A1),''|'',I3)'
         FMT2 = FMT2(1:LF)//CHAR(IC0+N3)//CHAR(IC0+N2)//CHAR(IC0+N1)
     &          //'(''--''),''|'',I3)'
         FMT3 = FMT3(1:LF)//'60(2X,I2))'
         FMT4 = FMT4(1:LF)//'60(I2,2X))'
         LF = 21
      ELSE
         IF ( MCOL.EQ.1 ) THEN
            NK = NINT(SQRT(DBLE(N)))
         ELSE IF ( MCOL.EQ.2 ) THEN
            NK = NINT(SQRT(DBLE(N/2)))
         ELSE IF ( MCOL.EQ.3 ) THEN
            NK = 2*NINT(SQRT(DBLE(N/2))) - 1
         END IF
         DO K = 1,NK
            IF ( MCOL.LE.2 ) THEN
               NM = 2*K - 1
            ELSE
               NM = 2*((K+1)/2)
            END IF
            NM2 = NM/10
            NM1 = NM - NM2*10
            NM3 = NM/2
            FMT1 = FMT1(1:LF)//CHAR(IC0+NM2)//CHAR(IC0+NM1)
     &             //'( 2A1),''|'','
            FMT2 = FMT2(1:LF)//CHAR(IC0+NM2)//CHAR(IC0+NM1)
     &             //'(''--''),''|'','
C
            IF ( MCOL.LE.2 ) THEN
               DO MM = 1,NM
                  IF ( MOD(MM,2).EQ.MOD(K,2) ) THEN
                     FMT3 = FMT3(1:L3)//'2X,'
                     FMT4 = FMT4(1:L3)//'I2,'
                  ELSE
                     FMT3 = FMT3(1:L3)//'I2,'
                     FMT4 = FMT4(1:L3)//'2X,'
                  END IF
                  L3 = L3 + 3
               END DO
               FMT3 = FMT3(1:L3)//'''|'','
               FMT4 = FMT4(1:L3)//'''|'','
               L3 = L3 + 4
            ELSE
               FMT3 = FMT3(1:LF)//CHAR(IC0+NM3)//'(2X,I2),''|'','
               FMT4 = FMT4(1:LF)//CHAR(IC0+NM3)//'(I2,2X),''|'','
               L3 = L3 + 13
            END IF
            LF = LF + 13
         END DO
         IF ( MCOL.EQ.2 ) THEN
            FMT1 = FMT1(1:LF)//FMT1(12:LF)
            FMT2 = FMT2(1:LF)//FMT2(12:LF)
            FMT3 = FMT3(1:L3)//FMT3(12:L3)
            FMT4 = FMT4(1:L3)//FMT4(12:L3)
            LF = 2*LF - 11
            L3 = 2*L3 - 11
         END IF
         FMT1 = FMT1(1:LF)//'I3)'
         FMT2 = FMT2(1:LF)//'I3)'
         FMT3 = FMT3(1:L3)//'I3)'
         FMT4 = FMT4(1:L3)//'I3)'
      END IF
      IF ( MLIN.EQ.0 ) THEN
         NSL = 1
         ILSEP(1) = N
      ELSE IF ( MLIN.EQ.1 ) THEN
         NSL = NINT(SQRT(DBLE(N)))
         DO IL = 1,NSL
            ILSEP(IL) = IL**2
         END DO
      ELSE IF ( MLIN.EQ.2 ) THEN
         NSL = NINT(SQRT(DBLE(N/2)))
         DO IL = 1,NSL
            ILSEP(IL) = IL**2
         END DO
         DO IL = 1,NSL
            ILSEP(NSL+IL) = ILSEP(NSL) + IL**2
         END DO
         NSL = 2*NSL
      ELSE IF ( MLIN.EQ.3 ) THEN
         NSL = 2*NINT(SQRT(DBLE(N/2))) - 1
         ILSEP(1) = 2
         DO K = 2,NSL
            ILSEP(K) = ILSEP(K-1) + 2*((K+1)/2)
         END DO
      END IF
C
C
      WRITE (NFIL,99001) STR(1:LSTR)
      IF ( IJQ.GT.1000 ) WRITE (NFIL,99003) IQ,JQ
      WRITE (NFIL,FMT3) (I,I=2,N,2)
      WRITE (NFIL,FMT4) (I,I=1,N,2)
      WRITE (NFIL,FMT=FMT2)
C------------------------------------------------------------ header end
      NNON0 = 0
      ND = 0
      CTAB(0) = ' '
      DTAB(0) = 9999D0
C
      DO I = 1,N
         DO J = 1,N
            IF ( .NOT.SMALL(B(I,J)) ) THEN
               NNON0 = NNON0 + 1
               DO ID = 1,ND
                  IF ( SAME(B(I,J),+DTAB(ID)) ) THEN
                     IW(J) = +ID
                     GOTO 50
                  END IF
                  IF ( SAME(B(I,J),-DTAB(ID)) ) THEN
                     IW(J) = -ID
                     GOTO 50
                  END IF
               END DO
C----------------------------------------------------------- new element
               ND = ND + 1
               IW(J) = ND
               DTAB(ND) = B(I,J)
               IF ( ABS(DTAB(ND)-1.0D0)*TOL.LT.1.0D0 ) THEN
                  CTAB(ND) = '1'
               ELSE IF ( ABS(DTAB(ND)+1.0D0)*TOL.LT.1.0D0 ) THEN
                  DTAB(ND) = +1.0D0
                  CTAB(ND) = '1'
                  IW(J) = -ND
               ELSE IF ( ABS(DTAB(ND)-CI)*TOL.LT.1.0D0 ) THEN
                  CTAB(ND) = 'i'
               ELSE IF ( ABS(DTAB(ND)+CI)*TOL.LT.1.0D0 ) THEN
                  DTAB(ND) = +CI
                  CTAB(ND) = 'i'
                  IW(J) = -ND
               ELSE
                  CTAB(ND) = CHAR(IPT(1+MOD((ND+1),NC)))
               END IF
            ELSE
               IW(J) = 0
            END IF
 50      END DO
C------------------------------------------------------------ write line
         WRITE (NFIL,FMT=FMT1) I,
     &                         (VZ(ISIGN(1,IW(J))),CTAB(ABS(IW(J))),J=1,
     &                         N),I
C
         DO ISL = 1,NSL
            IF ( I.EQ.ILSEP(ISL) ) WRITE (NFIL,FMT=FMT2)
         END DO
      END DO
C
C------------------------------------------------------------------ foot
C
      WRITE (NFIL,FMT4) (I,I=1,N,2)
      WRITE (NFIL,FMT3) (I,I=2,N,2)
C
      IF ( K_FMT_FIL.GT.0 ) THEN
         WRITE (NFIL,99004) (ID,CTAB(ID),DTAB(ID),ID=1,ND)
         WRITE (NFIL,99005) NNON0,TOLP,N*N - NNON0,TOLP
      ELSE
         WRITE (NFIL,*) ' '
      END IF
C
99001 FORMAT (/,8X,A,/)
99002 FORMAT (/,1X,79('*'),/,10X,'inconsistent call of <CMATSTR>',/,10X,
     &        'argument IJQ =',I8,'  implies IQ=',I3,'   JQ=',I3,/,10X,
     &        'IQ*N=',I6,' > M   or   JQ*N=',I6,' > M   for N =',I4,
     &        ' M=',I4,/,1X,79('*'),/)
99003 FORMAT (8X,'IQ-JQ-block  for  IQ = ',I3,'   JQ = ',I3,/)
99004 FORMAT (/,8X,'symbols used:',/,(8X,I3,3X,A1,2X,2F20.12))
99005 FORMAT (/,8X,I5,' elements   >',1PE9.1,/,
     &          8X,I5,' elements   <',1PE9.1,/)
      END SUBROUTINE



      END MODULE