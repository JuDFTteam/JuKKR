C ************************************************************************
      SUBROUTINE GREFSY(GTMAT,GMAT,IPVT,NDIM,DGTDE,
     +                  LLY_G0TR,
C                       new input parameters after inc.p removal
     &                  LMMAXD, LLY)
C ************************************************************************
      IMPLICIT NONE


C     Solves (1 - g0 \Delta t) G_ref = g0 for G_ref. (Full inversion)
C
C     GTMAT ... on input it has to contain (-1) * g0 * \Delta t
C               on output it contains G_ref
C               dimension NDIM x NDIM
C     GMAT  ... input: g0 free-space Green's function for the central
C                      reference cluster atom: g0^{(1)N'}_{LL'}
C     IPVT  ... integer work array of dimension NDIM
C     NDIM  ... GTMAT is a matrix of dimension NDIM x NDIM
C           ... NDIM = NACLSD * LMMAXD (#cluster atoms * maximal LM)
C     commented by E. Rabel, Nov 2011

C
C---> SOLVE THE DYSON EQUATION TO GET REFERENCE GREEN FUNCTION
C
C
      INTEGER LMMAXD
C     Lloyd's formula switch 0 (inactive)/ 1 (active)
      INTEGER LLY

      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO= (0.D0,0.D0))

C     ..
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL ZGETRF,ZGETRS
C     ..
C     .. SAVE STATEMENT ..
C     SAVE
C     ..
C     .. SCALAR ARGUMENTS ..
      INTEGER NDIM
      DOUBLE COMPLEX LLY_G0TR
C     ..
C     .. ARRAY ARGUMENTS ..
C     DOUBLE COMPLEX GMAT(NGD,LMGF0D),GTMAT(NGD,NGD),
C    +               DGTDE(LLYNGD,LMGF0D)
      DOUBLE COMPLEX GMAT (NDIM,LMMAXD)
      DOUBLE COMPLEX GTMAT(NDIM,NDIM)
      DOUBLE COMPLEX DGTDE(LLY*(NDIM-1)+1,LMMAXD)

C     ..
C
C     .. LOCAL ARRAYS ..
      INTEGER IPVT(NDIM)

C     .. LOCAL SCALARS ..
      INTEGER I,INFO

      INTEGER LMGF0D,NGD

      LMGF0D = LMMAXD
      NGD = NDIM
C     NGD=NACLSD*(LMAX+1)**2
C     LLYNGD=LLY*(NGD-1)+1
C     LLYNGD=LLY*(NACLSD*((LMAX+1)**2)-1)+1

 
      DO 10 I = 1,NDIM
        GTMAT(I,I) = CONE + GTMAT(I,I) ! GTMAT= 1 - G0 * \Delta T
   10 CONTINUE
C
C---> SOLVE THE SYSTEM OF LINEAR EQUATIONS
C
      CALL ZGETRF(NDIM,NDIM,GTMAT,NGD,IPVT,INFO)

      CALL ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,GMAT,NGD,INFO)


C ..  Lloyd
      IF (LLY.EQ.1) THEN

        CALL ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,DGTDE,NGD,INFO)

        LLY_G0TR = CZERO

        DO I = 1,LMGF0D
          LLY_G0TR = LLY_G0TR - DGTDE(I,I)
        ENDDO

      ELSE 

        LLY_G0TR = CZERO

      ENDIF
C .. Lloyd

      END
