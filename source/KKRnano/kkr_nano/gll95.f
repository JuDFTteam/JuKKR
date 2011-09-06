c**********************************************************************
      SUBROUTINE GLL95(E,CLEB,ICLEB,LOFLM,IEND,TREFLL,DTREFLL,ATOM,
     +                 REFPOT,RATOM,NATOM,ALAT,GREF0,DGDEOUT,
     +                 LLY_G0TR,ICLS,CLS,I3 )
c **********************************************************************
c
c     solution of the DYSON equation for a cluster of potentials
c     (TREFLL) centered at positions RATOM in free space,
c
c ----------------------------------------------------------------------
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'
c
      INTEGER LMGF0D,NGD
      PARAMETER (LMGF0D= (LMAXD+1)**2,NGD=LMGF0D*NACLSD)
      INTEGER LLYNGD
      PARAMETER (LLYNGD=LLY*(LMGF0D*NACLSD-1)+1)
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX E,LLY_G0TR
      DOUBLE PRECISION ALAT
      INTEGER I3,IEND,NATOM,ICLS,IPVT(NGD),INFO
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DTREFLL(LMGF0D,LMGF0D,NREFD)
      DOUBLE COMPLEX GREF0(NGD,LMGF0D),TREFLL(LMGF0D,LMGF0D,NREFD),
     +               DGTDE(LLYNGD,LMGF0D),
     +               DGTDE0(LLYNGD,LLYNGD),
     +               DGDE(LLYNGD,LLYNGD),DGDEOUT(LLYNGD,LMGF0D)
      DOUBLE PRECISION CLEB(*),RATOM(3,*)
      INTEGER ATOM(*),ICLEB(NCLEB,3),LOFLM(*),REFPOT(*),CLS(NAEZD)
C     ..
C     .. Local Scalars ..
      INTEGER I,LM1,LM2,N1,N2,NDIM,NLM1,NLM2
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX DGLLDE(LMGF0D,LMGF0D)
      DOUBLE COMPLEX GLL(LMGF0D,LMGF0D),GREF(NGD,NGD),GTREF(NGD,LMGF0D)
      DOUBLE PRECISION RDIFF(3)
C     ..
C     .. External Subroutines ..
      EXTERNAL GFREE,GREFSY,ZCOPY,ZGEMM
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE
C     ..
      IF (TEST('flow    ')) WRITE (6,FMT=*) '>>> GLL95'

      NDIM = LMGF0D*NATOM

c
c ---> construct free Green's function
c
      DO 70 N1 = 1,NATOM
        DO 60 N2 = 1,NATOM
          DO 10 I = 1,3
c            RDIFF(I) = (RATOM(I,N1) - RATOM(I,N2))*ALAT
c           changed P.Z. 4.7.97
            RDIFF(I) = - (RATOM(I,N1)-RATOM(I,N2))*ALAT
   10     CONTINUE
          IF (N1.NE.N2) THEN

            CALL GFREE(RDIFF,E,GLL,CLEB,ICLEB,LOFLM,IEND)

            DO 30 LM2 = 1,LMGF0D
              NLM2 = (N2-1)*LMGF0D + LM2
              DO 20 LM1 = 1,LMGF0D
                NLM1 = (N1-1)*LMGF0D + LM1
                GREF(NLM1,NLM2) = GLL(LM1,LM2)
   20         CONTINUE
   30       CONTINUE
          ELSE
            DO 50 LM2 = 1,LMGF0D
              NLM2 = (N2-1)*LMGF0D + LM2
              DO 40 LM1 = 1,LMGF0D
                NLM1 = (N1-1)*LMGF0D + LM1
                GREF(NLM1,NLM2) = CZERO
   40         CONTINUE
   50       CONTINUE
          END IF
   60   CONTINUE
   70 CONTINUE

      IF (TEST('flow    ')) WRITE (6,FMT=*) 'GFREE o.k.'
c ----------------------------------------------------------------------
      
      IF (LLY.EQ.1) THEN
c
c ---> construct derivative of free Green's function
c
      DO N1 = 1,NATOM
        DO N2 = 1,NATOM
          DO I = 1,3
            RDIFF(I) = - (RATOM(I,N1)-RATOM(I,N2))*ALAT
          END DO
          IF (N1.NE.N2) THEN

            CALL DGFREE(RDIFF,E,DGLLDE,CLEB,ICLEB,LOFLM,IEND)

            DO LM2 = 1,LMGF0D
              NLM2 = (N2-1)*LMGF0D + LM2
              DO LM1 = 1,LMGF0D
                NLM1 = (N1-1)*LMGF0D + LM1
                DGDE(NLM1,NLM2) = DGLLDE(LM1,LM2)
              END DO
            END DO
          ELSE
            DO LM2 = 1,LMGF0D
              NLM2 = (N2-1)*LMGF0D + LM2
              DO LM1 = 1,LMGF0D
                NLM1 = (N1-1)*LMGF0D + LM1
                DGDE(NLM1,NLM2) = CZERO
              END DO
            END DO
          END IF
        END DO
      END DO

      ENDIF


      CALL ZCOPY(NGD*LMGF0D,GREF,1,GREF0,1)

      IF (LLY.EQ.1) THEN

      DO N2 = 1,NATOM
        NLM2 = (N2-1)*LMGF0D + 1
        CALL ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,DGDE(1,NLM2),NGD,
     +             TREFLL(1,1,REFPOT(ABS(ATOM(N2)))),LMGF0D,
     +             CZERO,GTREF,NGD)
        CALL ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,GREF(1,NLM2),NGD,
     +             DTREFLL(1,1,REFPOT(ABS(ATOM(N2)))),LMGF0D,
     +             CONE,GTREF,NGD)
        CALL ZCOPY(NGD*LMGF0D,GTREF,1,DGTDE0(1,NLM2),1)
      END DO
      END IF

      DO 80 N2 = 1,NATOM
        NLM2 = (N2-1)*LMGF0D + 1
        CALL ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,GREF(1,NLM2),NGD,
     +             TREFLL(1,1,REFPOT(ABS(ATOM(N2)))),LMGF0D,
     +             CZERO,GTREF,NGD)
        CALL ZCOPY(NGD*LMGF0D,GTREF,1,GREF(1,NLM2),1)
   80 CONTINUE

      IF (LLY.EQ.1) THEN
        DO N2 = 1,LMGF0D
          DO N1 = 1,NGD 
            DGTDE(N1,N2) = DGTDE0(N1,N2)
          ENDDO
        ENDDO
      ENDIF


      CALL GREFSY(GREF,GREF0,IPVT,NDIM,ICLS,DGTDE,
     +            I3,CLS,LLY_G0TR)

      IF (LLY.EQ.1) THEN

      CALL ZGEMM('N','N',NDIM,LMGF0D,NDIM,-CONE,DGTDE0,NGD,
     +           GREF0,NGD,CONE,DGDE,NGD)

      CALL ZGETRS('N',NDIM,LMGF0D,GREF,NGD,IPVT,DGDE,NGD,INFO)

        DO N2 = 1,LMGF0D
          DO N1 = 1,NGD 
            DGDEOUT(N1,N2) = DGDE(N1,N2)
          ENDDO
        ENDDO

      ENDIF

      IF (TEST('flow    ')) WRITE (6,FMT=*) 'GREFSY o.k.'

c ----------------------------------------------------------------------
      RETURN
 
      END
