MODULE MOD_SIMPK

CONTAINS

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: KKRimp, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  SUBROUTINE SIMPK(F,FINT,IPAN,IRCUT,DRDI,IPAND)
    ! **********************************************************************
    !     this subroutine does an integration up to rcut of an
    !     real function f with an extended 3-point-simpson :
    !
    !                            rcut
    !                      fint = { f(r') dr'
    !                             0
    !
    !     modified for functions with kinks - at each kink the
    !     integration is restarted .
    !
    !     attention : input f is destroyed !
    !
    !-----------------------------------------------------------------------
    IMPLICIT NONE

    !       INCLUDE 'inc.p'
    !     .. Scalar Arguments ..
    DOUBLE PRECISION FINT
    INTEGER IPAN, IPAND
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION DRDI(*),F(*)
    INTEGER IRCUT(0:IPAND)
    !     ..
    !     .. Local Scalars ..
    DOUBLE PRECISION A1,A2
    INTEGER I,IEN,IP,IST,N
    !     ..
    !     .. External Functions ..
    !       DOUBLE PRECISION SIMPK_SSUM
    !       EXTERNAL SSUM
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC MOD
    !     ..
    A1 = 4.0D0/3.0D0
    A2 = 2.0D0/3.0D0
    FINT = 0.0D0
    
    DO IP = 1,IPAN
    
      !---> loop over kinks
      
      IST = IRCUT(IP-1) + 1
      IEN = IRCUT(IP)
      
      DO  I = IST,IEN
        F(I) = F(I)*DRDI(I)
      END DO
      IF (MOD(IEN-IST,2).EQ.0) THEN
        FINT = FINT + (F(IST)-F(IEN))/3.0D0
        IST = IST + 1
        N = (IEN-IST+1)/2

      ELSE
        !---> four point lagrange integration for the first step
        FINT = FINT + (9.0D0*F(IST)+19.0D0*F(IST+1)-5.0D0*F(IST+2)+ F(IST+3))/24.0D0 + (F(IST+1)-F(IEN))/3.0D0
        IST = IST + 2
        N = (IEN-IST+1)/2
      END IF
      
      !---> calculate with an extended 3-point-simpson
      
      FINT = FINT + A1*SIMPK_SSUM(N,F(IST),2) + A2*SIMPK_SSUM(N,F(IST+1),2)
    END DO
  
  END SUBROUTINE SIMPK

  DOUBLE PRECISION FUNCTION SIMPK_SSUM(N,V,IV)
    ! **********************************************************************
    !        sum up the first N elements of the double precision
    !        array V(*) with a stepwidth of IV
    ! ----------------------------------------------------------------------
    !     .. Scalar Arguments ..
    INTEGER IV,N
    !     ..
    !     .. Array Arguments ..
    DOUBLE PRECISION V(*)
    !     ..
    !     .. Local Scalars ..
    DOUBLE PRECISION VSUM
    INTEGER I,IBOT,ITOP
    !     ..
    IF (IV.GE.0) THEN
      IBOT = 1
      ITOP = 1 + (N-1)*IV

    ELSE
      IBOT = 1 - (N-1)*IV
      ITOP = 1
    END IF

    VSUM = 0.0D0
    DO I = IBOT,ITOP,IV
      VSUM = VSUM + V(I)
    END DO
    SIMPK_SSUM = VSUM
  END FUNCTION SIMPK_SSUM


END MODULE MOD_SIMPK


