! ************************************************************************
SUBROUTINE reg2(m,x,f,w,c,xm,fm)
!*****************************************************************
!                                                                *
!  Das Unterprogramm REG2 berechnet die Koeffizienten eines      *
!  Ausgleichspolynoms 2-ten Grades nach der diskreten Fehler-    *
!  quadratmethode von Gauss.                                     *
!                                                                *
!                                                                *
!  EINGABEPARAMETER:                                             *
!  =================                                             *
!  M        : Nummer der letzten Stuetzstelle.                    *
!  X        : 1-dim. Feld X(M) mit den Stuetzstellen.           *
!  F        : 1-dim. Feld F(M) mit den Funktionswerten zu den  *
!             Stuetzstellen.                                     *
!  W        : 1-dim. Feld W(M) mit den Gewichten.              *
!                                                                *
!                                                                *
!  AUSGABEPARAMETER:                                             *
!  =================                                             *
!  C        : 1-dim. Feld C(1:N) mit den Koeffizienten des Aus-  *
!             gleichspolynoms.                                   *
!                                                                *
!                                                                *
!  HILFSPARAMETER:                                               *
!  ===============                                               *
!  A        : 2-dim. Feld A(1:LDA,1:N+1).                        *
!  B        : 1-dim. Feld B(1:N+1).                              *
!                                                                *
!----------------------------------------------------------------*
!                                                                *
!  Autor     : Guido Dubois                                      *
!  Datum     : 30.05.87                                          *
!  Quellcode : FORTRAN 77                                        *
!                                                                *
!*****************************************************************

! Deklarationen.

IMPLICIT NONE
INTEGER N,LDA
! Polynom 2. Grades
PARAMETER(N=2,LDA=N+1)
DOUBLE PRECISION XM,FM
DOUBLE PRECISION X(*),F(*),W(*),C(*)

DOUBLE PRECISION DUMMY
DOUBLE PRECISION A(LDA,N+1),B(N+1)
INTEGER IPIV(N+1),INFO
INTEGER I,J,J1,K,K1,L,M

!  Berechnung der ersten Spalte der Matrix und der rechten Seite
!  des Gleichungssystems.

a(1,1)=0.0D0
b(1)=0.0D0
DO  i=1,m
  a(1,1)=a(1,1)+w(i)
  b(1)=b(1)+w(i)*f(i)
END DO
DO  j=1,n
  j1=j+1
  a(j1,1)=0.0D0
  b(j1)=0.0D0
  DO  i=1,m
    dummy=w(i)*x(i)**j
    a(j1,1)=a(j1,1)+dummy
    b(j1)=b(j1)+dummy*f(i)
  END DO
END DO

!  Berechnung der letzten Zeile der Matrix.

DO  k=1,n
  k1=k+1
  l=k+n
  a(j1,k1)=0.0D0
  DO  i=0,m
    a(j1,k1)=a(j1,k1)+w(i)*x(i)**l
  END DO
END DO

!  Vervollstaendigung der Matrix.

DO  k=1,n
  DO  i=1,n
    a(i,k+1)=a(i+1,k)
  END DO
END DO

!*****************************************************************
!     Loesung eines linearen Gleichungssystems                   *
!                    A * X = B                                   *
! ****************************************************************

!      CALL DGEF(A,LDA,N+1,IPIV)
CALL dgetrf(n+1,n+1,a,lda,ipiv,info)
!      CALL DGESM('N',A,LDA,N+1,IPIV,B,N+1,1)
CALL dgetrs('N',n+1,1,a,lda,ipiv,b,n+1,info)

!      write(6,FMT='(f12.4  )') (b(i),i=1,3)

DO  j=1,n+1
  c(j)=b(j)
END DO

xm = -c(2)/(2.d0*c(3))
fm = c(3)*xm*xm + c(2)*xm + c(1)

RETURN
END SUBROUTINE reg2
