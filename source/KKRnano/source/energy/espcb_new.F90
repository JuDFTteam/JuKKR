!>    Sums up core state energies taking into account degeneracies.
!>
!>     Note: LCORE can have values only from 0 to 3
!>    @param[out] ESPC core contribution to single particle energies
!>                l and spin resolved
! 13.10.95 ***************************************************************
      SUBROUTINE ESPCB_NEW(ESPC,NSPIN,ECORE,LCORE,LCOREMAX,NCORE)
! ************************************************************************
!
!     attention : energy zero ---> electro static zero
!
!                 since input potential and single particle energies
!                 are using muffin tin zero as zero the energy shift
!                 is cancelled in the kinetic energy contribution !
!
!     calculate the core contribution of the single particle energies
!     l and spin dependent .
!     attention : here are the results of the subroutine corel (stored
!                 in the common block core) used .
!                                        (see notes by b.drittler)
!
!                 modified for bandstructure code
!                               b.drittler   jan 1990
!-----------------------------------------------------------------------
      IMPLICIT NONE

!     .. Scalar Arguments ..
      INTEGER NSPIN
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION ECORE(20,2),ESPC(0:3,NSPIN)
      INTEGER LCORE(20,2),NCORE(2)  ! changed dimensions!!!
      INTEGER LCOREMAX
!     ..
!     .. Local Scalars ..
      INTEGER ISPIN,L,N
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
!     ..
!
!---> loop over reference atoms
!
        LCOREMAX = 0

!
!---> initialize espc
!
        ESPC = 0.0D0

        DO 30 ISPIN = 1,NSPIN
!
!---> loop over all core states
!
          DO 20 N = 1,NCORE(ISPIN)
            L = LCORE(N,ISPIN)
            LCOREMAX = MAX(LCOREMAX,L)
            ESPC(L,ISPIN) = ESPC(L,ISPIN) + ECORE(N,ISPIN)*DBLE(2*L+1)*DBLE(3-NSPIN)
   20     CONTINUE
   30   CONTINUE
      ENDsubroutine
