      SUBROUTINE SYMETRMAT(NSYM,CPREF,DSYMLL,SYMUNITARY,MATQ,
     &                     IQS,MATSYM,LMMAXD)
C **********************************************************************
C *                                                                    *
C *  Symmetrising the t/G matrix (or their inverses):                  *
C *                                                                    *
C *                       nsym                                         *
C *     MATSYM = CPREF *  SUM  [ DLL(i) * MATQ(iqs(i)) * DLL(i)^T ]    *
C *                      i = 1                                         *
C *                                                                    *
C *  IQS - set outside the routine - has either the same value         *
C *        regardless of i (e.g. in case of the single-site matrices)  *
C *        or is taking on the value of i => MATSYM is a cummulative   *
C *        sum over MATQ[1...nsym] (e.g. in case of the BZ-integration *
C *        of G)                                                       *
C *                                                                    *
C * CPREF = 1/NSYM  or 1/VBZ                                           *
C *                                                                    *
C *                                    v.popescu, munich nov. 2004     *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Parameters ..
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER ( CZERO = (0D0,0D0), CONE = (1D0,0D0) )
C     ..
C     .. Arguments ..
      INTEGER LMMAXD,NSYM
      INTEGER IQS(*)
      DOUBLE COMPLEX CPREF
      LOGICAL SYMUNITARY(*)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,*),MATSYM(LMMAXD,LMMAXD),
     &               MATQ(LMMAXD,LMMAXD,*)
C     ..
C     .. Locals ..
      INTEGER ISYM,L1,L2
      CHARACTER*1 CNT
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD),TS(LMMAXD,LMMAXD)
C     ..
      EXTERNAL ZCOPY,ZGEMM
C     ..
      CALL ZCOPY(LMMAXD*LMMAXD,MATQ(1,1,IQS(1)),1,TS,1)
C
C ----------------------------------------------------------------------
      DO ISYM = 2,NSYM
C
C --> look if the symmetry operation is unitary / anti-unitary 
C
         CNT = 'N'
         IF ( .NOT.SYMUNITARY(ISYM) ) CNT = 'T'
C
         CALL ZGEMM('N',CNT,LMMAXD,LMMAXD,LMMAXD,CONE,DSYMLL(1,1,ISYM),
     &              LMMAXD,MATQ(1,1,IQS(ISYM)),LMMAXD,CZERO,W1,LMMAXD)
C
         CALL ZGEMM('N','C',LMMAXD,LMMAXD,LMMAXD,CONE,W1,
     &              LMMAXD,DSYMLL(1,1,ISYM),LMMAXD,CONE,TS,LMMAXD)
      END DO
C ----------------------------------------------------------------------
C
      DO L2 = 1,LMMAXD
         DO L1 = 1,LMMAXD
            MATSYM(L1,L2) = CPREF * TS(L1,L2)
         END DO
      END DO
C
      END
