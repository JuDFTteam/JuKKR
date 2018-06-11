! ************************************************************************
    Subroutine reg2(m, x, f, w, c, xm, fm)
      Use mod_datatypes, Only: dp
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

      Implicit None
      Integer :: n, lda
! Polynom 2. Grades
      Parameter (n=2, lda=n+1)
      Real (Kind=dp) :: xm, fm
      Real (Kind=dp) :: x(*), f(*), w(*), c(*)

      Real (Kind=dp) :: dummy
      Real (Kind=dp) :: a(lda, n+1), b(n+1)
      Integer :: ipiv(n+1), info
      Integer :: i, j, j1, k, k1, l, m

!  Berechnung der ersten Spalte der Matrix und der rechten Seite
!  des Gleichungssystems.

      a(1, 1) = 0.0E0_dp
      b(1) = 0.0E0_dp
      Do i = 1, m
        a(1, 1) = a(1, 1) + w(i)
        b(1) = b(1) + w(i)*f(i)
      End Do
      Do j = 1, n
        j1 = j + 1
        a(j1, 1) = 0.0E0_dp
        b(j1) = 0.0E0_dp
        Do i = 1, m
          dummy = w(i)*x(i)**j
          a(j1, 1) = a(j1, 1) + dummy
          b(j1) = b(j1) + dummy*f(i)
        End Do
      End Do

!  Berechnung der letzten Zeile der Matrix.

      Do k = 1, n
        k1 = k + 1
        l = k + n
        a(j1, k1) = 0.0E0_dp
        Do i = 0, m
          a(j1, k1) = a(j1, k1) + w(i)*x(i)**l
        End Do
      End Do

!  Vervollstaendigung der Matrix.

      Do k = 1, n
        Do i = 1, n
          a(i, k+1) = a(i+1, k)
        End Do
      End Do

!*****************************************************************
!     Loesung eines linearen Gleichungssystems                   *
!                    A * X = B                                   *
! ****************************************************************

!      CALL DGEF(A,LDA,N+1,IPIV)
      Call dgetrf(n+1, n+1, a, lda, ipiv, info)
!      CALL DGESM('N',A,LDA,N+1,IPIV,B,N+1,1)
      Call dgetrs('N', n+1, 1, a, lda, ipiv, b, n+1, info)

!      write(6,FMT='(f12.4  )') (b(i),i=1,3)

      Do j = 1, n + 1
        c(j) = b(j)
      End Do

      xm = -c(2)/(2.E0_dp*c(3))
      fm = c(3)*xm*xm + c(2)*xm + c(1)

      Return
    End Subroutine
