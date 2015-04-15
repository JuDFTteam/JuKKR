c 13.10.95 ***************************************************************
      SUBROUTINE ESPCB(ESPC,NSPIN,NATYP,ECORE,LCORE,LCOREMAX,NCORE)
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
C     .. Parameters ..
      include 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
C *          function, set up in the spin-independent non-relativstic *
C *          (l,m_l)-representation                                   *
C *                                                                   *
C *********************************************************************
C
      INTEGER NPOTD
      PARAMETER (NPOTD=(2*KREL + (1-KREL)*NSPIND)*NATYPD)
      INTEGER LMAXD1
      PARAMETER (LMAXD1= LMAXD+1)
C     .. Scalar Arguments ..
      INTEGER NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ECORE(20,*),ESPC(0:3,NPOTD)
      INTEGER LCORE(20,*),NCORE(*),LCOREMAX(*)
C     ..
C     .. Local Scalars ..
      INTEGER I1,IPOT,IS,L,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DBLE
C     ..
c
c---> loop over reference atoms
c
      DO 40 I1 = 1,NATYP
        LCOREMAX(I1) = 0
        DO 30 IS = 1,NSPIN

c
c---> determine correct potential indices
c
          IPOT = NSPIN* (I1-1) + IS
c
c---> initialize espc
c
          DO 10 L = 0,3
            ESPC(L,IPOT) = 0.0D0
   10     CONTINUE
c
c---> loop over all core states
c
          DO 20 N = 1,NCORE(IPOT)
            L = LCORE(N,IPOT)
            LCOREMAX(I1) = MAX(LCOREMAX(I1),L)
            ESPC(L,IPOT) = ESPC(L,IPOT) +
     +                      ECORE(N,IPOT)*DBLE(2*L+1)*DBLE(3-NSPIN)
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
      END
