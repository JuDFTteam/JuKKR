C>    Sums up core state energies taking into account degeneracies.
C>
C>     Note: LCORE can have values only from 0 to 3
C>    @param[out] ESPC core contribution to single particle energies
C>                l and spin resolved
c 13.10.95 ***************************************************************
      SUBROUTINE ESPCB_NEW(ESPC,NSPIN,ECORE,LCORE,LCOREMAX,NCORE)
c ************************************************************************
c
c     attention : energy zero ---> electro static zero
c
c                 since input potential and single particle energies
c                 are using muffin tin zero as zero the energy shift
c                 is cancelled in the kinetic energy contribution !
c
c     calculate the core contribution of the single particle energies
c     l and spin dependent .
c     attention : here are the results of the subroutine corel (stored
c                 in the common block core) used .
c                                        (see notes by b.drittler)
c
c                 modified for bandstructure code
c                               b.drittler   jan 1990
c-----------------------------------------------------------------------
      IMPLICIT NONE

C     .. Scalar Arguments ..
      INTEGER NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ECORE(20,2),ESPC(0:3,NSPIN)
      INTEGER LCORE(20,2),NCORE(2)  ! changed dimensions!!!
      INTEGER LCOREMAX
C     ..
C     .. Local Scalars ..
      INTEGER ISPIN,L,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
c
c---> loop over reference atoms
c
        LCOREMAX = 0

c
c---> initialize espc
c
        ESPC = 0.0D0

        DO 30 ISPIN = 1,NSPIN
c
c---> loop over all core states
c
          DO 20 N = 1,NCORE(ISPIN)
            L = LCORE(N,ISPIN)
            LCOREMAX = MAX(LCOREMAX,L)
            ESPC(L,ISPIN) = ESPC(L,ISPIN) +
     +                      ECORE(N,ISPIN)*DBLE(2*L+1)*DBLE(3-NSPIN)
   20     CONTINUE
   30   CONTINUE
      END
