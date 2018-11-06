      MODULE MOD_INTERPOLPOT

      CONTAINS

!-------------------------------------------------------------------------------
!> Summary: This subroutine does the interpolation
!> Author: Who wrote this subroutine
!> Category: interpolation, potential, kkrimp
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------


      SUBROUTINE INTERPOLPOT(R0,R1,V0,V1,N0,N1,NLOC)
c  ******************************************************************
c  *
C  * This subroutine does the interpolation
C  ******************************************************************
C
c take every NNLSQ points and interpolate with lsq.
c
      IMPLICIT NONE
C     .. Parameters ..
C     .. Scalar Arguments ..
      INTEGER N0,N1, NLOC
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION R0(N0),R1(N1),V0(N0),V1(N1)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION FACTOR
      INTEGER I,ICON,IINT,ILAST,IP,IR,J,J1,NN,NNLSQ,NORDER
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION C(NLOC),R(NLOC),V(NLOC)
C     ..
C     .. External Subroutines ..
!       EXTERNAL INTERPOLPOT_LEASQR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      NNLSQ = 8
      NORDER = 3
      DO I = 1,N1
        V1(I) = 0.0
      END DO
c
      IINT = (N0-1)/NNLSQ
!       ICON = N1ST Bauer
      ICON = 0
      DO 40 J = 1,IINT
        FACTOR = V0(NNLSQ*J+1-NNLSQ/2)
        IF (ABS(FACTOR).LT.1.D-10) FACTOR = 1.D0
        NN = NNLSQ
        J1 = 1
        DO 10 I = 1,NNLSQ
          IR = (J-1)*NNLSQ + I
          IF ((R0(IR+1)-R0(IR)).GT.1.D-10) THEN
            R(J1) = R0(IR+1)
            V(J1) = V0(IR+1)/FACTOR
            J1 = J1 + 1
          ELSE
            NN = NN - 1
          END IF
   10   CONTINUE
        CALL INTERPOLPOT_LEASQR(R,V,NN,NORDER,C)
        ILAST = J*NNLSQ
        DO 30 I = ICON + 1,N1
          IF (R1(I).LT.R0(ILAST+1)) THEN
            DO 20 IP = 1,NORDER + 1
              V1(I) = V1(I) + C(IP)* (R1(I))** (IP-1)
   20       CONTINUE
            V1(I) = V1(I)*FACTOR
            ICON = I
c               write(6,*) i,r1(i),v1(i)
          END IF
   30   CONTINUE
   40 CONTINUE
      DO 50 I = ICON + 1,N1
        DO IP = 1,NORDER + 1
          V1(I) = V1(I) + C(IP)* (R1(I))** (IP-1)
        END DO
        V1(I) = V1(I)*FACTOR
   50 CONTINUE
      END SUBROUTINE
!-------------------------------------------------------------------------------
!> Summary: This subroutine is used to fit a polynomial to a set of data 
!> Author: Who wrote this subroutine
!> Category: interpolation, potential
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------


      SUBROUTINE INTERPOLPOT_LEASQR(R,V,N,NORDER,C)
      IMPLICIT NONE
C  ******************************************************************
C  *
C  *  This program is used in fitting a polynomial to a set of data
C  *  The program reads in N pairs of X and Y values and computes the
C  *  coefficients of the normal equations for the least squares
C  *  method.
C  *
C  *------------------------------------------------------------------
C  *  Parameters are:
C  *
C  *        X, Y   - Array of X and Y values
C  *        N      - Number of data pairs
C  *        MS,MF  - The range of the degree of polynomials
C  *                 to be computed. The MAXIMUM degree is 9
C  *        A      - Augmented array of the coefficients
C  *                 of the normal equations
C  *        C      - array of the coeficients of the least squares
C  *                 polynomials
C  *
C  *
C  *------------------------------------------------------------------
C  *
C  *
C  *
C
C
C
C
C
C
C      Read in N, then the X and Y values
C
C     READ(*,*) N,(X(I),Y(I),I=1,N)
C
c     DATA N/11/
c     DATA X/0.05,0.11,0.15,0.31,0.46,0.52,0.7,0.74,0.82,
c    &       0.98,1.17,89*0.0/
c     DATA Y/0.956,0.89,0.832,0.717,0.571,0.539,0.378,
c    &       0.37,0.306,0.242,0.104,89*0.0/
C
C
C   Read in MS, MF The program will find coefficients for each
C   degree of polynomial from degree MS to degree MF
C
c     DATA MS,MF/1,7/
C
C   Compute matrix of coefficients and R.H.S. for MF degree.
C   However, first check to see if max degree requested is too
C   large. Ti cannot exceed N-1. If it does reduce to equal to
C   N-1 and print message.
C
C ------------------------------
C     .. Parameters ..
      INTEGER IRMD
      PARAMETER (IRMD=1000)
C     ..
C     .. Scalar Arguments ..
      INTEGER N,NORDER
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION C(IRMD),R(IRMD),V(IRMD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BETA,SUM
      INTEGER I,ICOEF,IM1,IPT,J,JCOEF,MF,MFP1,MFP2,MS,MSP1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(10,11),X(IRMD),XN(IRMD),Y(IRMD)
C     ..
C     .. External Subroutines ..
!       EXTERNAL INTERPOLPOT_LUDCMQ,INTERPOLPOT_SOLNQ
C     ..
      DO I = 1,N
        X(I) = R(I)
        Y(I) = V(I)
      END DO
C    Initialize A
      DO I = 1,10
        DO J = 1,11
          A(I,J) = 0.0D0
        END DO
      END DO
      DO I = N + 1,IRMD
        X(I) = 0.0D0
        Y(I) = 0.0D0
      END DO
      MS = NORDER
      MF = NORDER
      IF (NORDER.GT.9) THEN
        STOP ' *** Order of polynomial greater than 9 ***'
      END IF
C -----------------------------
      IF (MF.GT. (N-1)) THEN
        MF = N - 1
c        WRITE(6,200) MF
      END IF
   10 MFP1 = MF + 1
      MFP2 = MF + 2
C
C
C   Put ones into a new array. This will hold the powers of the
C   X values as we procced.
C
      DO 20 I = 1,N
        XN(I) = 1.0D0
   20 CONTINUE
C
C
C   ------------------------------------------------------------------
C
C   compute first column and N+1 st column of A. I moves down the
C   rows, J sums over the N values.
C
      DO 40 I = 1,MFP1
        A(I,1) = 0.0D0
        A(I,MFP2) = 0.0D0
        DO 30 J = 1,N
          A(I,1) = A(I,1) + XN(J)
          A(I,MFP2) = A(I,MFP2) + Y(J)*XN(J)
          XN(J) = XN(J)*X(J)
   30   CONTINUE
   40 CONTINUE
C
C   Compute the last row of A. I moves across the columns, J
C   sums over the N values
C
      DO 60 I = 2,MFP1
        A(MFP1,I) = 0.0D0
        DO 50 J = 1,N
          A(MFP1,I) = A(MFP1,I) + XN(J)
          XN(J) = XN(J)*X(J)
   50   CONTINUE
   60 CONTINUE
C
C
C   Now fill in the rest of the A matrix. I moves down the rows
C   J moves across the columns.
C
      DO 80 J = 2,MFP1
        DO 70 I = 1,MF
          A(I,J) = A(I+1,J-1)
   70   CONTINUE
   80 CONTINUE
C
C
C   Write out the matrix of normal equations
C
c     WRITE(6,'(///)')
c     WRITE(6,*)'         THE NORMAL EQUATIONS ARE'
c     WRITE(6,*)
c     WRITE(6,201) ((A(I,J),J=1,MFP2), I=1,MFP1)
c     WRITE(6,*)
C
C   Now call a subroutine to solve the system. Do this for each
C   degree from MS to MF. Get the LU decomposition of A.
C
      CALL INTERPOLPOT_LUDCMQ(A,MFP1,10)
C
C   Reset the R.H.S. into C. We need to do this for each degree.
C
      MSP1 = MS + 1
      DO 120 I = MSP1,MFP1
        DO 90 J = 1,I
          C(J) = A(J,MFP2)
   90   CONTINUE
        CALL INTERPOLPOT_SOLNQ(A,C,I,10)
        IM1 = I - 1
C
C   Now write the coefficients of the least square polynomial.
C
c     WRITE(6,202) IM1, (C(J),J=1,I)
C
C   Compute and print the value of beta = Sum of Dev squared /
C   ( N - M - 1).
C
        BETA = 0.0D0
        DO 110 IPT = 1,N
          SUM = 0.0D0
          DO 100 ICOEF = 2,I
            JCOEF = I - ICOEF + 2
            SUM = (SUM+C(JCOEF))*X(IPT)
  100     CONTINUE
          SUM = SUM + C(1)
          BETA = BETA + (Y(IPT)-SUM)**2
  110   CONTINUE
        BETA = BETA/ (N-I)
c      WRITE(6,203) BETA
  120 CONTINUE
C
C
C
 9000 FORMAT (/,/,'DEGREE OF POLYNOMIAL CANNOT EXCEED N - 1.',/,
     +       ' REQUESTED MAXIMUM DEGREE TOO LARGE - ','REDUCED TO ',I3)
 9010 FORMAT (1X,9D12.4)
 9020 FORMAT (/,' FOR DEGREE OF ',I2,' COEFFICIENTS ARE',/,/,' ',5X,
     +       11D16.9)
 9030 FORMAT (9X,' BETA IS ',E15.7,/,/)
      END SUBROUTINE


      SUBROUTINE INTERPOLPOT_LUDCMQ(A,N,NDIM)
C
C  *******************************************************************
C  *  SUBROUTINE LUDCMQ :
C  *                     This subroutine forms the LU equivalent of
C  *  the qsuare coefficient matrix A. The LU, in compact form, is
C  *  returned in the A matrix space. The Upper triangular matrix U
C  *  has ones on its diagonal - these values are not included in
C  *  the result.
C  *
C  *******************************************************************
C
C
C     .. Scalar Arguments ..
      INTEGER N,NDIM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NDIM,NDIM)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SUM
      INTEGER I,IM1,J,JM1,K
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
      DO 50 I = 1,N
        DO 40 J = 2,N
          SUM = 0.0D0
          IF (J.LE.I) THEN
            JM1 = J - 1
            DO 10 K = 1,JM1
              SUM = SUM + A(I,K)*A(K,J)
   10       CONTINUE
            A(I,J) = A(I,J) - SUM
          ELSE
            IM1 = I - 1
            IF (IM1.NE.0) THEN
              DO 20 K = 1,IM1
                SUM = SUM + A(I,K)*A(K,J)
   20         CONTINUE
            END IF
C
C   Test for small value on the diagonal
C
   30       IF (DABS(A(I,I)).LT.1.D-20) THEN
              A(I,I) = 1.D-20
              WRITE (6,FMT=9000) I,A(I,I)
              RETURN
            ELSE
              A(I,J) = (A(I,J)-SUM)/A(I,I)
            END IF
          END IF
   40   CONTINUE
   50 CONTINUE
      RETURN
 9000 FORMAT (' REDUCTION NOT COMPLITED BECAUSE SMALL VALUE',
     +       ' FOUND FOR DIVISOR IN ROW ',I3,D11.4)
      END SUBROUTINE


      SUBROUTINE INTERPOLPOT_SOLNQ(A,B,N,NDIM)
C
C  *******************************************************************
C  *
C  *  SUBROUTINE SOLNQ :
C  *                    This subroutine finds the solution to a set
C  *  of N linear equations that corresponds to the right hand side
C  *  vector B. The A matrix is the LU decomposition equivalent to
C  *  the coefficient matrix of the original equations, as produced by
C  *  LUDCMQ. The solution vector is returned in the vector B.
C  *
C*******************************************************************
C
C   Do the reduction step
C
C     .. Scalar Arguments ..
      INTEGER N,NDIM
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(NDIM,NDIM),B(NDIM)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION SUM
      INTEGER I,IM1,J,K,NMJP1,NMJP2
C     ..
      B(1) = B(1)/A(1,1)
      DO 20 I = 2,N
        IM1 = I - 1
        SUM = 0.0D0
        DO 10 K = 1,IM1
          SUM = SUM + A(I,K)*B(K)
   10   CONTINUE
        B(I) = (B(I)-SUM)/A(I,I)
   20 CONTINUE
C
C   Now we are ready for the back substitution. remember that the
C   elements of U on the diagonal are all ones.
C
C
      DO 40 J = 2,N
        NMJP2 = N - J + 2
        NMJP1 = N - J + 1
        SUM = 0.0D0
        DO 30 K = NMJP2,N
          SUM = SUM + A(NMJP1,K)*B(K)
   30   CONTINUE
        B(NMJP1) = B(NMJP1) - SUM
   40 CONTINUE
      RETURN
      END SUBROUTINE

      END MODULE
