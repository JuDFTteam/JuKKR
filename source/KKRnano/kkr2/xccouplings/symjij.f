      SUBROUTINE SYMJIJ(
     >                  ALAT,TAUVBZ,
     >                  NSYMAT,DSYMLL,
     >                  NXIJ,IXCP,
     >                  TMATLL,MSSQ,
     >                  GSXIJ,
     <                  GMATXIJ,
C                       new input parameters after inc.p removal
     &                  naez, lmmaxd, nxijd)
C =======================================================================
C
C  a) symmetrize GSXIJ
C
C  b) GMATXIJ = - Delta_t^(-1) * GLL * Delta_t^(-1)
C                                                     A.Thiess Sep'09
C =======================================================================
C
      IMPLICIT NONE

      INTEGER naez
      INTEGER lmmaxd
      INTEGER nxijd

      INTEGER          NSYMAXD
      PARAMETER        (NSYMAXD=48)

      DOUBLE COMPLEX   CONE,CZERO
      PARAMETER        (CONE  = ( 1.0D0,0.0D0))
      PARAMETER        (CZERO  = ( 0.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
C     ..
      DOUBLE COMPLEX   TAUVBZ
      DOUBLE PRECISION ALAT
      INTEGER          NSYMAT,XIJ,NXIJ
C     ..
C     .. Array Arguments ..
C     ..

      DOUBLE COMPLEX   GMATXIJ(LMMAXD,LMMAXD,NXIJD)
      DOUBLE COMPLEX   GSXIJ  (LMMAXD,LMMAXD,NSYMAXD,NXIJD)
      DOUBLE COMPLEX   TMATLL (LMMAXD,LMMAXD,NAEZ)
      DOUBLE COMPLEX   MSSQ   (LMMAXD,LMMAXD)
      DOUBLE COMPLEX   DSYMLL (LMMAXD,LMMAXD,NSYMAXD)

      INTEGER          IXCP(NXIJD)
C     ..
C     .. Local Scalars ..
C     ..
      DOUBLE PRECISION RFCTOR
      INTEGER          IU,LM1,LM2,INFO
C     ..
C     .. Local Arrays ..
C     ..

      DOUBLE COMPLEX       GLL(LMMAXD,LMMAXD)
      DOUBLE COMPLEX   MSSXCPL(LMMAXD,LMMAXD)
      DOUBLE COMPLEX       TPG(LMMAXD,LMMAXD)
      DOUBLE COMPLEX        XC(LMMAXD,LMMAXD)
      DOUBLE COMPLEX        W1(LMMAXD,LMMAXD)
      INTEGER             IPVT(LMMAXD)
C     ..
      EXTERNAL ZCOPY,ZAXPY,ZGETRF,ZGETRS,ZGETRI,ZGEMM,ZSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE

C     ..
C     ! = ALAT/(2*PI)
      RFCTOR = ALAT/(8.D0*ATAN(1.0D0))



C================================
       DO XIJ = 2, NXIJ
C================================
         DO LM2 = 1,LMMAXD
           DO LM1 = 1,LMMAXD
             MSSXCPL(LM1,LM2) = TMATLL(LM1,LM2,IXCP(XIJ))
           END DO
         END DO
C
C ---> inversion 
C
         CALL ZGETRF(LMMAXD,LMMAXD,MSSXCPL,LMMAXD,IPVT,INFO)
         CALL ZGETRI(LMMAXD,MSSXCPL,LMMAXD,IPVT,W1,
     &               LMMAXD*LMMAXD,INFO)
C
C
C-------------------------------------------------------- SYMMETRISE GLL
C
         DO IU = 1,NSYMAT
C
C --->    GLL = sum(i=1,iumax)(tauvbz * DLL(i) * GS * DLL(i)^T)
C
            IF ( IU.EQ.1 ) THEN
C
C --->    ull(1) is equal to unity matrix
C
               CALL ZCOPY(LMMAXD*LMMAXD,GSXIJ(1,1,1,XIJ),1,GLL,1)
               CALL ZSCAL(LMMAXD*LMMAXD,TAUVBZ,GLL,1)
C
            ELSE
C
C --->      tpg = tauvbz * DLL * GS
C
               CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,TAUVBZ,
     &         DSYMLL(1,1,IU),LMMAXD,GSXIJ(1,1,IU,XIJ),LMMAXD,
     &                    CZERO,TPG,LMMAXD)
C
C --->    GLL = GLL + TPG * DLL(i)^T
C
               CALL ZGEMM('N','C',LMMAXD,LMMAXD,LMMAXD,CONE,TPG,LMMAXD,
     &                    DSYMLL(1,1,IU),LMMAXD,CONE,GLL,LMMAXD)
            END IF
C
         END DO

C-------------------------------------------------------- IU = 1,NSYMAT
C
C In case of more than one atom per unit cell different MSSQ's have to
C be prepared
C
C --->  XC = Delta_t^(-1) * GLL
C
        CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,MSSQ,
     &             LMMAXD,GLL,LMMAXD,CZERO,XC,LMMAXD)
C
C --->  GLL = - Delta_t^(-1) * GLL * Delta_t^(-1)
C
        CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD,
     &             MSSXCPL,LMMAXD,-CZERO,GLL,LMMAXD)
C
C --->  GMATXIJ = GMATLL = GLL/RFCTOR
C
        DO LM1 = 1,LMMAXD
          DO LM2 = 1,LMMAXD
            GMATXIJ(LM2,LM1,XIJ) = GLL(LM2,LM1)/RFCTOR
          END DO
        END DO

C================================
      ENDDO
C================================
Cxccpl
C
      RETURN

      END
