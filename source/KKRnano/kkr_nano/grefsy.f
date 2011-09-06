C ************************************************************************
      SUBROUTINE GREFSY(GTMAT,GMAT,IPVT,NDIM,ICLS,DGTDE,
     +                  I3,CLS,LLY_G0TR)
C ************************************************************************
C
C---> SOLVE THE DYSON EQUATION TO GET REFERENCE GREEN FUNCTION
C
C     .. PARAMETERS ..
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'
C
      INTEGER LMGF0D,NGD
      PARAMETER (LMGF0D= (LMAXD+1)**2,NGD=LMGF0D*NACLSD)
      INTEGER LLYNGD
      PARAMETER (LLYNGD=LLY*(LMGF0D*NACLSD-1)+1)
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO= (0.D0,0.D0))
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,INFO
C     ..
C     .. LOCAL ARRAYS ..
      INTEGER IPVT(NGD)
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL ZGETRF,ZGETRS
C     ..
C     .. SAVE STATEMENT ..
      SAVE
C     ..
C     .. SCALAR ARGUMENTS ..
      INTEGER NDIM,ICLS,I3
      DOUBLE COMPLEX LLY_G0TR
C     ..
C     .. ARRAY ARGUMENTS ..
      DOUBLE COMPLEX GMAT(NGD,LMGF0D),GTMAT(NGD,NGD),
     +               DGTDE(LLYNGD,LMGF0D)
      INTEGER CLS(NAEZD)
C     ..
C
 
      DO 10 I = 1,NDIM
        GTMAT(I,I) = CONE + GTMAT(I,I) ! GTMAT= 1 - G * T
   10 CONTINUE
C
C---> SOLVE THE SYSTEM OF LINEAR EQUATIONS
C
      CALL ZGETRF(NDIM,NDIM,GTMAT,NGD,IPVT,INFO)

      CALL ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,GMAT,NGD,INFO)


C ..  Lloyd
C      IF (LLY.EQ.1.AND.ICLS.EQ.CLS(I3)) THEN
      IF (LLY.EQ.1) THEN

        CALL ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,DGTDE,NGD,INFO)

        LLY_G0TR = CZERO
        TRACEG0TR = CZERO

        DO I = 1,LMGF0D
          LLY_G0TR = LLY_G0TR - DGTDE(I,I)
        ENDDO

      ELSE 

        LLY_G0TR = CZERO

      ENDIF
C .. Lloyd



      RETURN

      END
