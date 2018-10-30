subroutine calcJij(Gmat,tmatspinup,tmatspindown)

!
! ==>  get DTJXIJ = T(UP) - T(DOWN) for all atoms
!


        DO XIJ = 2, NXIJ ! loop XIJ = 1, NXIJ(I1)
!
! ==>  get the off-diagonal Green function matrix Gij(UP) and Gji(DOWN)
!
         DO ISPIN = 1,2
           IF ( ISPIN.EQ.1 ) THEN
             CALL ZCOPY(LMMAXD*LMMAXD,GMATXIJ(1,1,XIJ,ISPIN),
     +                  1,GMIJ,1)
           ELSE
             DO LM2 = 1,LMMAXD
               DO LM1 = 1,LMMAXD
!
! -> use Gji = Gij^T
!
                  GMJI(LM1,LM2) = GMATXIJ(LM2,LM1,XIJ,ISPIN)

               ENDDO
             ENDDO
           ENDIF
         ENDDO

! ----------------------------------------------------------------------
!
! ==> calculate the exchange coupling constant J_ij via Eq. (19)
!     modified for G instead of tau:
!          J_ij ~ Trace [ (t_i(D)-t_i(U)) * Gij(U)
!                       * (t_j(D)-t_j(U)) * Gji(D)]
!
! -------------------------------------------------- loop over occupants
! --> Delta_j * Gjt,it
!
         CALL CMATMUL(LMMAXD,LMMAXD,DTNXIJ(1,1,IXCP(XIJ)),GMJI,W2)
!
! --> Delta_i * Git,jt
!
         CALL CMATMUL(LMMAXD,LMMAXD,DTIXIJ,GMIJ,W3)
!
! --> Delta_i * Git,jt * Delta_j * Gjt,it
!
         CALL CMATMUL(LMMAXD,LMMAXD,W3,W2,W1)
!
         CSUM = CZERO
         DO LM1 = 1,LMMAXD
           CSUM = CSUM + W1(LM1,LM1)
         ENDDO
!
         JOUT = -DIMAG(WGTE*CSUM*JSCAL)
!


END SUBROUTINE CALCJIJ



      SUBROUTINE CMATMUL(N,M,A,B,C)
C   ********************************************************************
C   *                                                                  *
C   *   perform  the matrix-matrix operation           C = A * B       *
C   *                                                                  *
C   *   A,B,C   complex  SQUARE  N x N - matrices                      *
C   *   N       dimension of A, B and C                                *
C   *   M       array size of A, B, C with M >= N                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT DOUBLE COMPLEX(A-H,O-Z)
C
C PARAMETER definitions
C
      DOUBLE COMPLEX C0
      PARAMETER (C0=(0.0D0,0.0D0))
C
C Dummy arguments
C
      INTEGER M,N
      DOUBLE COMPLEX A(M,M),B(M,M),C(M,M)
C
C Local variables
C
      DOUBLE COMPLEX BLJ
      INTEGER I,J,L
C
      DO J = 1,N
         DO I = 1,N
            C(I,J) = C0
         END DO
      END DO
C
      DO J = 1,N
         DO L = 1,N
            BLJ = B(L,J)
            IF ( BLJ.NE.C0 ) THEN
               DO I = 1,N
                  C(I,J) = C(I,J) + A(I,L)*BLJ
               END DO
            END IF
         END DO
      END DO
C
      END SUBROUTINE CMATMUL















