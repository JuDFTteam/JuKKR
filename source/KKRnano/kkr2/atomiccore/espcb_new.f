C             TODO: MODIFY RHOCORE and
C             COREL first that it gives output only
C             for one atom
C             new call for main2
C             call ESPCB_NEW(ESPC,NSPIN,I1,ECORE,LCORE,LCOREMAX,NCORE)

C>    @param[out] ESPC core contribution to single particle energies
C>                l and spin resolved
c 13.10.95 ***************************************************************
      SUBROUTINE ESPCB_NEW(ESPC,NSPIN,I1,ECORE,LCORE,LCOREMAX,NCORE)
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
      INTEGER LCORE(20,*),NCORE(*),LCOREMAX
C     ..
C     .. Local Scalars ..
      INTEGER I1,IPOT,ISPIN,L,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
c
c---> loop over reference atoms
c
        LCOREMAX = 0
        DO 30 ISPIN = 1,NSPIN

c
c---> determine correct potential indices
c
          IPOT = NSPIN* (I1-1) + ISPIN
c
c---> initialize espc
c
          DO 10 L = 0,3
            ESPC(L,ISPIN) = 0.0D0
   10     CONTINUE
c
c---> loop over all core states
c
          DO 20 N = 1,NCORE(IPOT)
            L = LCORE(N,IPOT)
            LCOREMAX = MAX(LCOREMAX,L)
            ESPC(L,ISPIN) = ESPC(L,ISPIN) +
     +                      ECORE(N,ISPIN)*DBLE(2*L+1)*DBLE(3-NSPIN)
   20     CONTINUE
   30   CONTINUE
      END
