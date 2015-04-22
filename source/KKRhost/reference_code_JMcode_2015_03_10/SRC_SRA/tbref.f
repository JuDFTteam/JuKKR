      SUBROUTINE TBREF(EZ,IELAST,ALATC,VREF,IEND,LMAX,NCLS,NINEQ,NREF,
     +                 CLEB,RCLS,ATOM,CLS,ICLEB,LOFLM,NACLS,
     +                 REFPOT,RMTREF,TOLRDIF,TMPDIR,ITMPDIR,ILTMP,
     &                 NAEZ,LLY) ! LLY Lloyd
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
C *          function, set up in the spin-independent non-relativstic *
C *          (l,m_l)-representation                                   *
C *                                                                   *
C *********************************************************************
C
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAXD+1)**2)
      INTEGER LM2D
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      INTEGER LRECGRF,LRECGRF1
!       PARAMETER (LRECGRF=4*NACLSD*LMGF0D*LMGF0D*NCLSD)
      PARAMETER (LRECGRF=WLENGTH*4*NACLSD*LMGF0D*LMGF0D*NCLSD) 
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO=(0.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALATC
      DOUBLE PRECISION TOLRDIF ! Set free GF to zero for r < tolrdif
      INTEGER IELAST,IEND,LMAX,NCLS,NINEQ,NREF,NAEZ
      INTEGER ITMPDIR,ILTMP
      INTEGER LLY ! LLY <> 0 : alpply Lloyds formula
      CHARACTER*80 TMPDIR
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX EZ(IEMXD)
      DOUBLE PRECISION CLEB(NCLEB,2),RCLS(3,NACLSD,NCLSD),
     +                 RMTREF(NREFD),VREF(NREFD)
      INTEGER ATOM(NACLSD,NAEZD+NEMBD),CLS(NAEZD+NEMBD),ICLEB(NCLEB,4),
     +        LOFLM(LM2D),NACLS(NCLSD),REFPOT(NAEZD+NEMBD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ERYD
      INTEGER I1,IC,ICLS,IE,LM1,LM2,LRECGF1,NACLSMAX
      DOUBLE COMPLEX   LLY_G0TR_IE ! LLY
CMPI  INTEGER MAPBLOCK
CMPI  INTEGER MYRANK,NROFNODES
CMPI  COMMON /MPI/MYRANK,NROFNODES
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX TREFLL(LMGF0D,LMGF0D,NREFD),
     &               DTREFLL(LMGF0D,LMGF0D,NREFD)                     ! LLY 
      DOUBLE COMPLEX,ALLOCATABLE:: GINP(:,:,:),
     &                             DGINP(:,:,:)                       ! LLY
      DOUBLE COMPLEX ALPHAREF(0:LMAXD,NREFD),DALPHAREF(0:LMAXD,NREFD) ! LLY Lloyd Alpha matrix and deriv.
      DOUBLE COMPLEX   LLY_G0TR(IEMXD,NCLSD) ! LLY
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
C     .. External Subroutines ..
      EXTERNAL CALCTREF13,GLL13
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      
      NACLSMAX = 1
      DO IC = 1,NCLS
         IF (NACLS(IC).GT.NACLSMAX) NACLSMAX = NACLS(IC)
      ENDDO
      LRECGRF1 = WLENGTH*4*NACLSMAX*LMGF0D*LMGF0D*NCLS 

      ALLOCATE (  GINP(NACLSMAX*LMGF0D,LMGF0D,NCLS) )
      ALLOCATE ( DGINP(NACLSMAX*LMGF0D,LMGF0D,NCLS) )


      CALL OPENDAFILE(68,'gref',4,LRECGRF1,TMPDIR,ITMPDIR,ILTMP)
      IF (LLY.NE.0) THEN
         CALL OPENDAFILE(681,'dgrefde',7,LRECGRF1,TMPDIR,ITMPDIR,ILTMP)
         OPEN(682,FILE='lly_g0tr_ie.ascii',FORM='FORMATTED')
      ENDIF

C
C ======================================================================

      DO IE = 1,IELAST
         WRITE(*,FMT='(A23,I5,2F14.8)') 
     &                         'TBREF: GREF for energy:',IE,EZ(IE)
C
CMPI      IF(MYRANK.EQ.MAPBLOCK(IE,1,IELAST,1,0,NROFNODES-1)) THEN
C
         ERYD = EZ(IE)
         DO I1 = 1,NREF
            CALL CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAX,LM1,              ! LLY Lloyd
     &                    TREFLL(1,1,I1),DTREFLL(1,1,I1),                   ! LLY Lloyd
     &                    ALPHAREF(0,I1),DALPHAREF(0,I1),LMAXD+1,LMGF0D)    ! LLY Lloyd

            IF (TEST('flow    ')) WRITE (6,FMT=*) 'tll(ref),  i1 = ',I1
         END DO
C
         IF (TEST('flow    ')) WRITE (6,FMT=*) 't-mat(Ref) o.k.', IE
C ----------------------------------------------------------------------
         DO ICLS = 1,NCLS
            I1 = 1
            IC = 0
            DO WHILE (IC.EQ.0 .AND. I1.LE.NINEQ)
               IF (CLS(I1).EQ.ICLS) IC = I1
               I1 = I1 + 1
            END DO
C     
            IF (IC.EQ.0) STOP 'Error in CLS(*) array in tbref'
            IF (TEST('flow    ')) WRITE (6,FMT=*) 'CLUSTER ',ICLS,
     +           ' at ATOM ',IC
C     
            CALL GLL13(ERYD,CLEB(1,2),ICLEB,LOFLM,IEND,TREFLL,DTREFLL,
     +                 ATOM(1,IC),REFPOT,RCLS(1,1,ICLS),NACLS(ICLS),
     &                 TOLRDIF,
     +                 ALATC,0,GINP(1,1,ICLS),
     &                 DGINP(1,1,ICLS),NACLSMAX,LLY_G0TR(IE,ICLS),LLY)
         END DO

         IF (LLY.NE.0) THEN                                    ! LLY Lloyd
            LLY_G0TR_IE = CZERO                                ! LLY Lloyd
            DO I1 = 1,NAEZ                                     ! LLY Lloyd
               ICLS = CLS(I1)                                  ! LLY Lloyd
               LLY_G0TR_IE = LLY_G0TR_IE + LLY_G0TR(IE,ICLS)   ! LLY Lloyd
            ENDDO                                              ! LLY Lloyd
         ENDIF                                                 ! LLY Lloyd
C ----------------------------------------------------------------------
         WRITE (68,REC=IE) GINP                 ! Gref
         IF (LLY.NE.0) WRITE (681,REC=IE) DGINP ! dGref/dE     ! LLY Lloyd
         IF (LLY.NE.0) WRITE (682,FMT='(2E24.16)') LLY_G0TR_IE ! LLY Lloyd
C
         IF (TEST('flow    ')) WRITE (6,FMT=*) 'G(n,lm,n,lm) (Ref) o.k.'

CMPI      END IF
      END DO
C ======================================================================
      CLOSE (68)
      IF (LLY.NE.0) CLOSE (681)
      DEALLOCATE (GINP,DGINP)

      END
