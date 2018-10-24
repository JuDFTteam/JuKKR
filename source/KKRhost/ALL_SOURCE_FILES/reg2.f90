!------------------------------------------------------------------------------------
!> Summary: The subroutine `REG2` calculates the coefficients of a equalization polynomial of the 2nd degree after the discrete error square method of Gauss
!> Author:  Guido Dubois
!> The subroutine `REG2` calculates the coefficients of a equalization polynomial of 
!> the 2nd degree after the discrete error square method of Gauss
!------------------------------------------------------------------------------------
!> @note This module seems to be used nowhere
!> @endnote
!------------------------------------------------------------------------------------
module mod_reg2
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: The subroutine `REG2` calculates the coefficients of a equalization polynomial of the 2nd degree after the discrete error square method of Gauss
  !> Author: Guido Dubois
  !> Category: numerical-tools, deprecated, KKRhost 
  !> Deprecated: True
  !> The subroutine `REG2` calculates the coefficients of a equalization polynomial of 
  !> the 2nd degree after the discrete error square method of Gauss
  !-------------------------------------------------------------------------------
  !> @note This module seems to be used nowhere
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine reg2(m, x, f, w, c, xm, fm)
    ! *****************************************************************
    ! *
    ! Das Unterprogramm REG2 berechnet die Koeffizienten eines      *
    ! Ausgleichspolynoms 2-ten Grades nach der diskreten Fehler-    *
    ! quadratmethode von Gauss.                                     *
    ! *
    ! *
    ! EINGABEPARAMETER:                                             *
    ! =================                                             *
    ! M        : Nummer der letzten Stuetzstelle.                    *
    ! X        : 1-dim. Feld X(M) mit den Stuetzstellen.           *
    ! F        : 1-dim. Feld F(M) mit den Funktionswerten zu den  *
    ! Stuetzstellen.                                     *
    ! W        : 1-dim. Feld W(M) mit den Gewichten.              *
    ! *
    ! *
    ! AUSGABEPARAMETER:                                             *
    ! =================                                             *
    ! C        : 1-dim. Feld C(1:N) mit den Koeffizienten des Aus-  *
    ! gleichspolynoms.                                   *
    ! *
    ! *
    ! HILFSPARAMETER:                                               *
    ! ===============                                               *
    ! A        : 2-dim. Feld A(1:LDA,1:N+1).                        *
    ! B        : 1-dim. Feld B(1:N+1).                              *
    ! *
    ! ----------------------------------------------------------------*
    ! *
    ! Autor     : Guido Dubois                                      *
    ! Datum     : 30.05.87                                          *
    ! Quellcode : FORTRAN 77                                        *
    ! *
    ! *****************************************************************

    ! Deklarationen.

    implicit none
    integer :: n, lda
    ! Polynom 2. Grades
    parameter (n=2, lda=n+1)
    real (kind=dp) :: xm, fm
    real (kind=dp), dimension(*) :: x, f, w, c

    real (kind=dp) :: dummy
    real (kind=dp), dimension(n+1) :: b
    real (kind=dp), dimension(lda, n+1) :: a
    integer :: info
    integer, dimension(n+1) :: ipiv
    integer :: i, j, j1, k, k1, l, m

    ! Berechnung der ersten Spalte der Matrix und der rechten Seite
    ! des Gleichungssystems.

    a(1, 1) = 0.0e0_dp
    b(1) = 0.0e0_dp
    do i = 1, m
      a(1, 1) = a(1, 1) + w(i)
      b(1) = b(1) + w(i)*f(i)
    end do
    do j = 1, n
      j1 = j + 1
      a(j1, 1) = 0.0e0_dp
      b(j1) = 0.0e0_dp
      do i = 1, m
        dummy = w(i)*x(i)**j
        a(j1, 1) = a(j1, 1) + dummy
        b(j1) = b(j1) + dummy*f(i)
      end do
    end do

    ! Berechnung der letzten Zeile der Matrix.

    do k = 1, n
      k1 = k + 1
      l = k + n
      a(j1, k1) = 0.0e0_dp
      do i = 0, m
        a(j1, k1) = a(j1, k1) + w(i)*x(i)**l
      end do
    end do

    ! Vervollstaendigung der Matrix.

    do k = 1, n
      do i = 1, n
        a(i, k+1) = a(i+1, k)
      end do
    end do

    ! *****************************************************************
    ! Loesung eines linearen Gleichungssystems                   *
    ! A * X = B                                   *
    ! ****************************************************************

    ! CALL DGEF(A,LDA,N+1,IPIV)
    call dgetrf(n+1, n+1, a, lda, ipiv, info)
    ! CALL DGESM('N',A,LDA,N+1,IPIV,B,N+1,1)
    call dgetrs('N', n+1, 1, a, lda, ipiv, b, n+1, info)

    ! write(6,FMT='(f12.4  )') (b(i),i=1,3)

    do j = 1, n + 1
      c(j) = b(j)
    end do

    xm = -c(2)/(2.e0_dp*c(3))
    fm = c(3)*xm*xm + c(2)*xm + c(1)

    return
  end subroutine reg2

end module mod_reg2
