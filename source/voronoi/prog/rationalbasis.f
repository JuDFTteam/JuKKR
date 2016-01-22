      SUBROUTINE RATIONALBASIS(
     >     BRAVAIS,RECBV,NBASIS,
     X     RBASIS)
      implicit none
c     This subroutine rationalises the basis vectors
c     by shifting them by appropriate lattice vectors
c     so that they are closest to the origin.

c     Input:
      INTEGER NBASIS                 ! Number of basis vectors
      REAL*8 BRAVAIS(3,3),RECBV(3,3) ! Bravais vectors of real and rec. lattice
                                     ! 1st index: xyz, second index: abc

c     Input and output:
      REAL*8 RBASIS(3,*)             ! Basis vectors, to be changed on output.

c     Local:
      INTEGER IBASIS,IX,NB1,NB2,NB3
      REAL*8 RLATT(3)
      


      DO IBASIS = 1,NBASIS
c Find closest lattice point to basis site. 
c Use reciprocal lattice vectors (remember they are given in units 2 pi/a).
c Remember, NINT(x) is the closest integer to x.
         NB1 = NINT( RBASIS(1,IBASIS)*RECBV(1,1) + 
     &               RBASIS(2,IBASIS)*RECBV(2,1) + 
     &               RBASIS(3,IBASIS)*RECBV(3,1)   )
         NB2 = NINT( RBASIS(1,IBASIS)*RECBV(1,2) + 
     &               RBASIS(2,IBASIS)*RECBV(2,2) + 
     &               RBASIS(3,IBASIS)*RECBV(3,2)   )
         NB3 = NINT( RBASIS(1,IBASIS)*RECBV(1,3) + 
     &               RBASIS(2,IBASIS)*RECBV(2,3) + 
     &               RBASIS(3,IBASIS)*RECBV(3,3)   )
c The closest lattice point is 
c RLATTT = NB1 * Bravais1 + NB2 * Bravais2 + NB3 * Bravais3.
c Shift the basis site by this lattice vector.
         DO IX = 1,3
            RLATT(IX) = NB1 * BRAVAIS(IX,1) + 
     &                  NB2 * BRAVAIS(IX,2) + NB3 * BRAVAIS(IX,3)
            RBASIS(IX,IBASIS) = RBASIS(IX,IBASIS) - RLATT(IX)
         ENDDO
      ENDDO

      RETURN
      END
