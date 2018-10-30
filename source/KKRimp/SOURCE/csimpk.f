!-------------------------------------------------------------------------------
!> Summary: Integration of a complex function with extended 3-point-Simpson
!> Author: 
!> Date: 
!-------------------------------------------------------------------------------
!> This subroutine does an integration up to rcut of an
!> complex function cf with an extended 3-point-simpson :
!>
!> \begin{equation}
!> cf_{int} = \int_{0}^{r_{cut}} cf\left(r'\right) dr'
!> \end{equation}
!>
!> Modified for functions with kinks: at each kink the integration is restarted.
!-------------------------------------------------------------------------------
!> @warning Input cf is destroyed!
!> @endwarning
!-------------------------------------------------------------------------------
      MODULE MOD_CSIMPK
      CONTAINS
!-------------------------------------------------------------------------------
!> Summary: Integration of a complex function with extended 3-point-Simpson
!> Author: 
!> Date: 
!> Category: KKRimp, numerical-tools, radial-grid
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!-------------------------------------------------------------------------------
!> This subroutine does an integration up to rcut of an
!> complex function cf with an extended 3-point-simpson :
!>
!> \begin{equation}
!> cf_{int} = \int_{0}^{r_{cut}} cf\left(r'\right) dr'
!> \end{equation}
!>
!> Modified for functions with kinks: at each kink the integration is restarted.
!-------------------------------------------------------------------------------
      SUBROUTINE CSIMPK(CF,CFINT,IPAN,IRCUT,DRDI)
C     .. Scalar Arguments ..
      DOUBLE COMPLEX CFINT
      INTEGER IPAN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CF(*)
      DOUBLE PRECISION DRDI(*)
      INTEGER IRCUT(0:IPAN)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,A2
      INTEGER I,IEN,IP,IST,N
C     ..
C     .. External Functions ..
      DOUBLE COMPLEX CSUM
      EXTERNAL CSUM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      A1 = 4.0D0/3.0D0
      A2 = 2.0D0/3.0D0
      CFINT = 0.0D0
c
      DO 20 IP = 1,IPAN
c
c---> loop over kinks
c
        IST = IRCUT(IP-1) + 1
        IEN = IRCUT(IP)
c
        DO 10 I = IST,IEN
          CF(I) = CF(I)*DRDI(I)
   10   CONTINUE
c
        IF (MOD(IEN-IST,2).EQ.0) THEN
          CFINT = CFINT + (CF(IST)-CF(IEN))/3.0D0
          IST = IST + 1
          N = (IEN-IST+1)/2

        ELSE
c---> four point lagrange integration for the first step
          CFINT = CFINT + (9.0D0*CF(IST)+19.0D0*CF(IST+1)-
     +            5.0D0*CF(IST+2)+CF(IST+3))/24.0D0 +
     +            (CF(IST+1)-CF(IEN))/3.0D0
          IST = IST + 2
          N = (IEN-IST+1)/2
        END IF
c
c---> calculate with an extended 3-point-simpson
c
        CFINT = CFINT + A1*THIS_CSUM(N,CF(IST),2)
     +                + A2*THIS_CSUM(N,CF(IST+1),2)
   20 CONTINUE
c
      END SUBROUTINE CSIMPK

      COMPLEX*16 FUNCTION THIS_CSUM(N,V,IV)
      IMPLICIT NONE
c **********************************************************************
c        sum up the first N elements of the double complex
c        array V(*) with a stepwidth of IV
c ----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IV,N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX V(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX VSUM
      INTEGER I,IBOT,ITOP
C     ..
      IF (IV.GE.0) THEN
        IBOT = 1
        ITOP = 1 + (N-1)*IV

      ELSE
        IBOT = 1 - (N-1)*IV
        ITOP = 1
      END IF

      VSUM = (0D0,0D0)
      DO 10 I = IBOT,ITOP,IV
        VSUM = VSUM + V(I)
   10 CONTINUE
      THIS_CSUM = VSUM
      RETURN
      END FUNCTION THIS_CSUM

      END MODULE MOD_CSIMPK
