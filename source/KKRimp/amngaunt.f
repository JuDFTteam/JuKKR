!------------------------------------------------------------------------------------
!> Summary: Module handling the Gaunt coefficients for the structure constants used in the intersite potential
!> Author:
!> For details see vinters2010
!------------------------------------------------------------------------------------
      MODULE MOD_AMNGAUNT

      CONTAINS

!-------------------------------------------------------------------------------
!> Summary: Gaunt coefficients for the structure constants used in the intersite potential
!> Author:
!> Category: electrostatics, potential, special-functions, KKRimp
!> Deprecated: False 
!> For details see vinters2010
!> This sub calculates the gaunt coefs in the ordering neaded
!> for calculating the electrostatic AMAT^(nn')_(LL') 
!>                                              17.09.2001
!-------------------------------------------------------------------------------
      SUBROUTINE AMNGAUNT(LMAX,CLEB,ICLEB,IEND,W,YR,n,lassld,ncleb,lm3d)
      use mod_ymy, only: ymy        
      use mod_types, only: t_inc
      implicit none
!       include 'gaunt.param'
!       include 'parameters.file'
c      INTEGER LMAXD,LMX,LPOTD
c      PARAMETER (LMAXD=3,LPOTD=6)
!       INTEGER LMMAXD,L3D,LM3D
!       PARAMETER (LMMAXD= (LPOTD+1)**2,L3D=2*LPOTD,LM3D= (L3D+1)**2)
      INTEGER N,LASSLD,LM3D
!       PARAMETER (N=4* LMAXD,LASSLD=N)
!       INTEGER LMPOTD
!       PARAMETER (LMPOTD=(LPOTD+1)**2)
      INTEGER NCLEB
!       PARAMETER (NCLEB=LM3D*LMMAXD)
C     ..
      INTEGER LMAX
C     .. Array Arguments ..
      REAL*8 W(N),YR(N,0:LASSLD,0:LASSLD)

C     .. Local Scalars ..
      REAL*8 CLECG,EPI,FACTOR,FPI,PI,S
      INTEGER I,IEND,J,L1,L2,L3,L3MAX,
     +        LMMAX,M1,M1A,M1S,M2,M2A,M2S,M3,
     +        M3A,M3S
C     ..
C     .. Local Arrays ..
      REAL*8  CLEB(NCLEB,2)
      INTEGER ICLEB(NCLEB,4)
C     ..
C     .. Statement Functions ..
C      INTEGER MOFLM
C     ..
C     .. Save statement ..
      SAVE PI
C     ..
C     .. External Subroutines ..
      EXTERNAL ROTCOEF
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,REAL,SIGN
C     ..
C     .. Data statements ..
      DATA PI/3.14159265359D0/
C     ..
C     .. Statement Function definitions ..
c
c---> calculate the m-value for given lm and l
c
!       MOFLM(ILM,IL) = ILM - IL**2 - IL - 1
C     ..
c
      FPI = 4.0D0*PI
      EPI = 8.0D0*PI
      L3MAX = 2*LMAX
      LMMAX = (LMAX+1)**2
!       write(*,*) n,lassld,ncleb,lm3d
c
c---> set up of the gaunt coefficients with an index field
c     recognize that they are needed here only for l3=l1+l2
c
      I = 1
      DO 100 L1 = 0,LMAX
        DO 90 L2 = 0,LMAX
          L3 = L1 + L2
          DO 80 M1 = -L1,L1
            DO 70 M2 = -L2,L2
              DO 60 M3 = -L3,L3
                M1S = SIGN(1,M1)
                M2S = SIGN(1,M2)
                M3S = SIGN(1,M3)
c
                IF (M1S*M2S*M3S.GE.0) THEN
c
                  M1A = ABS(M1)
                  M2A = ABS(M2)
                  M3A = ABS(M3)
c
                  FACTOR = 0.0D0
c
                  IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR +
     +                REAL(3*M3S+SIGN(1,-M3))/8.0D0
                  IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR + REAL(M1S)/4.0D0
                  IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR + REAL(M2S)/4.0D0
c
                  IF (FACTOR.NE.0.0D0) THEN
c
                    IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR.
     +                  M1S*M3S.NE.1) FACTOR = -FACTOR
c
                    S = 0.0D0
                    DO 50 J = 1,N
                      S = S + W(J)*YR(J,L1,M1A)*YR(J,L2,M2A)*
     +                    YR(J,L3,M3A)
   50               CONTINUE
                    CLECG = S*FACTOR
                    IF (ABS(CLECG).GT.1.D-10) THEN
                      CLEB(I,1) = CLECG
                      ICLEB(I,1) = L1* (L1+1) + M1 + 1
                      ICLEB(I,2) = L2* (L2+1) + M2 + 1
                      ICLEB(I,3) = L3* (L3+1) + M3 + 1
                      I = I + 1
                    END IF

                  END IF

                END IF

   60         CONTINUE
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
      IEND = I - 1
      IF (NCLEB.LT.IEND) THEN
        STOP 13

      ELSE
        if (t_inc%i_write>0) WRITE (1337,FMT=900) IEND
      END IF
 900  FORMAT ('Gaunt coefs for VINTERS :',I6)
      END SUBROUTINE      
      END MODULE MOD_AMNGAUNT
