      SUBROUTINE GIJDMAT(TAUQ,TSST,MSSQ,DMAT,DTIL,CFCTORINV,IPRINT,
     &                   IE,IT,KREL,LMMAXD)
C **********************************************************************
C * Subroutine to get the projection matrices                          *
C *                                                                    *
C *                ii     i       i    (-1)                            *
C *    D  = [ 1 + G   * (t    -  t  ) ]                                *
C *     a          CPA    CPA     a                                    *
C *                                                                    *
C *    _            i       i      ii   (-1)                           *
C *    D  = [ 1 + (t    -  t  ) * G    ]                               *
C *     a           CPA     a      CPA                                 *
C *                                                                    *
C * Used is made of the same routine GETDMAT as in the case of TAU     *
C * projection. Note, however, that there the matrices are defined for *
C * example as                                                         *
C *                                                                    *
C *                ij     i       i      (-1)                          *
C *    D  = [ 1 + tau * (m    -  m    ) ]                              *
C *     a          CPA    a       CPA                                  *
C *                                                                    *
C * with m = (t)**(-1)                                                 *
C *                                                                    *
C *                                      v.popescu Oct. 2004           *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Arguments ..
      INTEGER IPRINT,LMMAXD
      INTEGER IE,IT,KREL
      DOUBLE COMPLEX CFCTORINV
      DOUBLE COMPLEX TAUQ(LMMAXD,LMMAXD),TSST(LMMAXD,LMMAXD),
     &               MSSQ(LMMAXD,LMMAXD)
      DOUBLE COMPLEX DMAT(LMMAXD,LMMAXD),DTIL(LMMAXD,LMMAXD)
C     ..
C     .. Locals ..
      INTEGER IK,INFO,I1,I2
      INTEGER IPVT(LMMAXD)
      DOUBLE COMPLEX CONE,CZERO
      DOUBLE COMPLEX GLL(LMMAXD,LMMAXD),TSSQ(LMMAXD,LMMAXD),
     &               TPG(LMMAXD,LMMAXD),XC(LMMAXD,LMMAXD)
      CHARACTER*18 BANNER
C     ..
C     .. Externals
      EXTERNAL CMATSTR,GETDMAT,ZCOPY,ZGEMM,ZGETRF,ZGETRI,ZSCAL
C     ..
C     .. Data
      DATA CONE / (1D0,0D0) /
      DATA CZERO / (0D0,0D0) /
C 
C --> get G(CPA) using the same procedure as for GMATLL in < KLOOPZ >
C           G(CPA) = -MSSQ - MSSQ * TAUQ * MSSQ 
C
      DO I2 = 1,LMMAXD
         DO I1 = 1,LMMAXD
            TPG(I1,I2) = -CONE*TAUQ(I1,I2)
         END DO
      END DO
C
      CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,MSSQ,
     &           LMMAXD,TPG,LMMAXD,CZERO,XC,LMMAXD)
C
      CALL ZCOPY(LMMAXD*LMMAXD,MSSQ,1,GLL,1)
C
      CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD,
     &                 MSSQ,LMMAXD,-CONE,GLL,LMMAXD)
C
      CALL ZSCAL(LMMAXD*LMMAXD,CFCTORINV,GLL,1)
C
C --> invert MSSQ to get TSSQ 
C
      DO I2 = 1,LMMAXD
         DO I1 = 1,LMMAXD
            TSSQ(I1,I2) = MSSQ(I1,I2) * CFCTORINV
         END DO
      END DO
      CALL ZGETRF(LMMAXD,LMMAXD,TSSQ,LMMAXD,IPVT,INFO)
      CALL ZGETRI(LMMAXD,TSSQ,LMMAXD,IPVT,XC,LMMAXD*LMMAXD,INFO)
C            
C --> get projection matrices from G(CPA) 
C     call getdmat with t,c for the local arg. c,t
C
      CALL GETDMAT(GLL,DMAT,DTIL,XC,LMMAXD,TSST,TSSQ,LMMAXD)
C
      IF ( IPRINT.LE.1 ) RETURN
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      IK = 2*KREL + 1
      WRITE(BANNER,'("DMAT IE=",I3," IT=",I3)') IE,IT
      CALL CMATSTR(BANNER,18,DMAT,LMMAXD,LMMAXD,IK,IK,0,1D-10,6)
C
      WRITE(BANNER,'("DTIL IE=",I3," IT=",I3)') IE,IT
      CALL CMATSTR(BANNER,18,DTIL,LMMAXD,LMMAXD,IK,IK,0,1D-10,6)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C
      END
