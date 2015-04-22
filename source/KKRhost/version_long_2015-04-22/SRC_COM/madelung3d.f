C*==madelung3d.f    processed by SPAG 6.05Rc at 11:55 on 19 May 2004
      SUBROUTINE MADELUNG3D(LPOT,YRG,WG,NAEZ,ALAT,VOLUME0,
     &                      BRAVAIS,RECBV,RBASIS,RMAX,GMAX,
     &                      NAEZD,LMXSPD,LASSLD,LPOTD,LMPOTD,
     &                      NMAXD,ISHLD,NEMBD,WLENGTH)
C **********************************************************************
C *                                                                    *
C * This subroutine calculates the Madelung potential coefficients     *
C * in the 3D case and stores them in the DA-file abvmad.unformatted   *
C * The record index is simply (IQ1-1)*NAEZ + IQ2 for record (IQ1,IQ2) *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER LPOT,NAEZ,WLENGTH
      INTEGER NAEZD,LMXSPD,LASSLD,LPOTD,LMPOTD,NMAXD,ISHLD,NEMBD
      DOUBLE PRECISION ALAT,VOLUME0,RMAX,GMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION YRG(LASSLD,0:LASSLD,0:LASSLD),WG(LASSLD)
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION RBASIS(3,NAEZD+NEMBD)
C     ..
C     .. Local Scalars ..
      INTEGER IEND,IPRINT,IQ1,IQ2,LFMT,LM1,LM2,NCLEBD
      INTEGER NGMAX,NRMAX,NSHLG,NSHLR
      INTEGER LRECABMAD,IREC
      CHARACTER*80 FMT
C     ..
C     .. Local Arrays ..
C     .. Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
      DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD),BVMAD(LMPOTD)
      DOUBLE PRECISION CLEB(LMXSPD*LMPOTD)
      DOUBLE PRECISION MADELSMAT(LMXSPD,NAEZD,NAEZD)
      DOUBLE PRECISION SMAT1(6,6),SMAT2(6,6)
      DOUBLE PRECISION GN(3,NMAXD),RM(3,NMAXD)
      INTEGER ICLEB(LMXSPD*LMPOTD,3)
      INTEGER NSG(ISHLD),NSR(ISHLD)
C     ..
C     .. External subroutines
      EXTERNAL LATTICE3D,STRMAT,MADELGAUNT,MADELCOEF
C ......................................................................
      IPRINT = 0
      NCLEBD = LMXSPD*LMPOTD
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,'(79(1H=))')
      WRITE (6,'(18X,A)') 
     &                  'MADELUNG3D: setting bulk Madelung coefficients'
      WRITE (6,'(79(1H=))')
      WRITE (6,*)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ======================================================================
      CALL LATTICE3D(ALAT,BRAVAIS,RECBV,NGMAX,NRMAX,NSHLG,NSHLR,NSG,NSR,
     &               GN,RM,RMAX,GMAX,IPRINT,NMAXD,ISHLD)
C
      CALL STRMAT(ALAT,LPOT,NAEZ,NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM,
     &            RBASIS,MADELSMAT,VOLUME0,IPRINT,LASSLD,LMXSPD,NAEZD)
C ======================================================================
C
      LRECABMAD = WLENGTH*2*LMPOTD*LMPOTD + WLENGTH*2*LMPOTD
      OPEN (69,ACCESS='direct',RECL=LRECABMAD,FILE='abvmad.unformatted',
     &      FORM='unformatted')
C
C --> calculate the gaunt coefficients
C
      CALL MADELGAUNT(LPOT,YRG,WG,CLEB,ICLEB,IEND,LASSLD,NCLEBD)
C
C --> calculate the madelung coefficients to be used for VMAD
C     call MADELCOEF with first arg. .FALSE. = 3D case
C
      DO IQ1 = 1,NAEZ
         DO IQ2 = 1,NAEZ
            CALL MADELCOEF(.FALSE.,LPOT,AVMAD,BVMAD,MADELSMAT(1,IQ1,IQ2)
     &                     ,CLEB,ICLEB,IEND,LPOTD,LMPOTD,LMXSPD,NCLEBD)
C
            IREC = IQ2 + NAEZ*(IQ1-1)
            WRITE (69,REC=IREC) AVMAD,BVMAD
C-----------------------------------------------------------------------
            IF ( (IQ1.LE.6) .AND. (IQ2.LE.6) ) THEN
               SMAT1(IQ1,IQ2) = AVMAD(1,1)
               SMAT2(IQ1,IQ2) = BVMAD(1)
            END IF
C-----------------------------------------------------------------------
         END DO
      END DO
      CLOSE (69)
C
      IF ( IPRINT.LT.1 ) RETURN
C ======================================================================
C
      CALL MADEL3OUT(IPRINT,NAEZ,LRECABMAD,SMAT1,SMAT2,LMPOTD)
C
      END
