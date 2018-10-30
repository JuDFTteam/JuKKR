c ************************************************************************
      SUBROUTINE REG2(M,X,F,W,C,XM,FM)
C*****************************************************************
C                                                                *
C  Das Unterprogramm REG2 berechnet die Koeffizienten eines      *
C  Ausgleichspolynoms 2-ten Grades nach der diskreten Fehler-    *
C  quadratmethode von Gauss.                                     *
C                                                                *
C                                                                *
C  EINGABEPARAMETER:                                             *
C  =================                                             *
C  M        : Nummer der letzten Stuetzstelle.                    *
C  X        : 1-dim. Feld X(M) mit den Stuetzstellen.           *
C  F        : 1-dim. Feld F(M) mit den Funktionswerten zu den  *
C             Stuetzstellen.                                     *
C  W        : 1-dim. Feld W(M) mit den Gewichten.              *
C                                                                *
C                                                                *
C  AUSGABEPARAMETER:                                             *
C  =================                                             *
C  C        : 1-dim. Feld C(1:N) mit den Koeffizienten des Aus-  *
C             gleichspolynoms.                                   *
C                                                                *
C                                                                *
C  HILFSPARAMETER:                                               *
C  ===============                                               *
C  A        : 2-dim. Feld A(1:LDA,1:N+1).                        *
C  B        : 1-dim. Feld B(1:N+1).                              *
C                                                                *
c----------------------------------------------------------------*
C                                                                *
C  Autor     : Guido Dubois                                      *
C  Datum     : 30.05.87                                          *
C  Quellcode : FORTRAN 77                                        *
C                                                                *
c*****************************************************************
C
C  Deklarationen.
      IMPLICIT NONE
      INTEGER N,LDA
c     Polynom 2. Grades
      PARAMETER(N=2,LDA=N+1)
      DOUBLE PRECISION XM,FM
      DOUBLE PRECISION X(*),F(*),W(*),C(*)

      DOUBLE PRECISION DUMMY
      DOUBLE PRECISION A(LDA,N+1),B(N+1)
      INTEGER IPIV(N+1),INFO
      INTEGER I,J,J1,K,K1,L,M
C
C
C  Berechnung der ersten Spalte der Matrix und der rechten Seite
C  des Gleichungssystems.
C
      A(1,1)=0.0D0
      B(1)=0.0D0
      DO 10 I=1,M
         A(1,1)=A(1,1)+W(I)
         B(1)=B(1)+W(I)*F(I)
   10 CONTINUE
      DO 20 J=1,N
         J1=J+1
         A(J1,1)=0.0D0
         B(J1)=0.0D0
         DO 30 I=1,M
            DUMMY=W(I)*X(I)**J
            A(J1,1)=A(J1,1)+DUMMY
            B(J1)=B(J1)+DUMMY*F(I)
   30    CONTINUE
   20 CONTINUE
C
C  Berechnung der letzten Zeile der Matrix.
C
      DO 40 K=1,N
         K1=K+1
         L=K+N
         A(J1,K1)=0.0D0
         DO 50 I=0,M
            A(J1,K1)=A(J1,K1)+W(I)*X(I)**L
   50    CONTINUE
   40 CONTINUE
C
C  Vervollstaendigung der Matrix.
C
      DO 60 K=1,N
         DO 70 I=1,N
            A(I,K+1)=A(I+1,K)
   70    CONTINUE
   60 CONTINUE
C
C*****************************************************************
C     Loesung eines linearen Gleichungssystems                   *
C                    A * X = B                                   *
c ****************************************************************
c
c      CALL DGEF(A,LDA,N+1,IPIV)
      CALL DGETRF(N+1,N+1,A,LDA,IPIV,INFO)
c      CALL DGESM('N',A,LDA,N+1,IPIV,B,N+1,1)
      CALL DGETRS('N',N+1,1,A,LDA,IPIV,B,N+1,INFO)

c      write(6,FMT='(f12.4  )') (b(i),i=1,3)

      DO 80 J=1,N+1
        C(J)=B(J)
 80   CONTINUE

      XM = -C(2)/(2.D0*C(3))
      FM = C(3)*XM*XM + C(2)*XM + C(1)

      RETURN
      END
