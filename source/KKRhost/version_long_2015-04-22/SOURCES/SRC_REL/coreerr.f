C*==coreerr.f    processed by SPAG 6.05Rc at 12:37 on 29 Apr 2001
      SUBROUTINE COREERR(ERR,VAR,S,NSOL,POW,QOW,PIW,QIW)
C   ********************************************************************
C   *                                                                  *
C   *   CALCULATE THE MISMATCH OF THE RADIAL WAVE FUNCTIONS AT THE     *
C   *   POINT  NMATCH  FOR OUT- AND INWARD INTEGRATION                 *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER NSOL,S
      REAL*8 ERR(4),PIW(2,2),POW(2,2),QIW(2,2),QOW(2,2),VAR(4)
C
C Local variables
C
      INTEGER T
C
C*** End of declarations rewritten by SPAG
C
      ERR(1) = POW(S,S) - PIW(S,S)*VAR(2)
      ERR(2) = QOW(S,S) - QIW(S,S)*VAR(2)
C
      IF ( NSOL.EQ.1 ) RETURN
C
      T = 3 - S
C
      ERR(1) = ERR(1) + POW(S,T)*VAR(3) - PIW(S,T)*VAR(2)*VAR(4)
      ERR(2) = ERR(2) + QOW(S,T)*VAR(3) - QIW(S,T)*VAR(2)*VAR(4)
      ERR(3) = POW(T,S) - PIW(T,S)*VAR(2) + POW(T,T)*VAR(3) - PIW(T,T)
     &         *VAR(2)*VAR(4)
      ERR(4) = QOW(T,S) - QIW(T,S)*VAR(2) + QOW(T,T)*VAR(3) - QIW(T,T)
     &         *VAR(2)*VAR(4)
C
      END
