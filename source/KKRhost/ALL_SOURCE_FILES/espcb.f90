! 13.10.95 ***************************************************************
SUBROUTINE espcb(espc,nspin,natyp,ecore,lcore,lcoremax,ncore)
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

implicit none

!     .. Parameters ..
INCLUDE 'inc.p'
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************

INTEGER NPOTD
PARAMETER (NPOTD=(2*KREL + (1-KREL)*NSPIND)*NATYPD)
INTEGER LMAXD1
PARAMETER (LMAXD1= LMAXD+1)

!.. Scalar Arguments ..
INTEGER NATYP,NSPIN

!.. Array Arguments ..
DOUBLE PRECISION ECORE(20,*),ESPC(0:3,NPOTD)
INTEGER LCORE(20,*),NCORE(*),LCOREMAX(*)

!.. Local Scalars ..
INTEGER I1,IPOT,IS,L,N

!.. Intrinsic Functions ..
INTRINSIC DBLE

!loop over reference atoms
DO I1 = 1,NATYP
  LCOREMAX(I1) = 0
  DO IS = 1,NSPIN


    ! determine correct potential indices
    IPOT = NSPIN* (I1-1) + IS

    ! initialize espc
    DO L = 0,3
      ESPC(L,IPOT) = 0.0D0
    END DO

    !loop over all core states
    DO N = 1,NCORE(IPOT)
      L = LCORE(N,IPOT)
      LCOREMAX(I1) = MAX(LCOREMAX(I1),L)
      write(1337,*) 'in espcb:',N,L,IPOT,I1
      ESPC(L,IPOT) = ESPC(L,IPOT) + &
                      ECORE(N,IPOT)*DBLE(2*L+1)*DBLE(3-NSPIN)
    END DO
  END DO
END DO
END
