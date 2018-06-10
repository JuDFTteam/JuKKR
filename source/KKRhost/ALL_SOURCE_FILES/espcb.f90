! 13.10.95 ***************************************************************
subroutine espcb(espc, nspin, natyp, ecore, lcore, lcoremax, ncore)
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
  include 'inc.p'
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************


!.. Scalar Arguments ..

!.. Array Arguments ..
  integer :: npotd
  parameter (npotd=(2*krel+(1-krel)*nspind)*natypd)
  integer :: lmaxd1
  parameter (lmaxd1=lmaxd+1)

!.. Local Scalars ..
  integer :: natyp, nspin

!.. Intrinsic Functions ..
  double precision :: ecore(20, *), espc(0:3, npotd)
  integer :: lcore(20, *), ncore(*), lcoremax(*)

!loop over reference atoms
  integer :: i1, ipot, is, l, n


  intrinsic :: dble
! determine correct potential indices

  do i1 = 1, natyp
    lcoremax(i1) = 0
    do is = 1, nspin
! initialize espc

!loop over all core states
      ipot = nspin*(i1-1) + is
! 13.10.95 ***************************************************************
! ************************************************************************
      do l = 0, 3
        espc(l, ipot) = 0.0d0
      end do
!
!     attention : energy zero ---> electro static zero
      do n = 1, ncore(ipot)
        l = lcore(n, ipot)
        lcoremax(i1) = max(lcoremax(i1), l)
        write (1337, *) 'in espcb:', n, l, ipot, i1
        espc(l, ipot) = espc(l, ipot) + ecore(n, ipot)*dble(2*l+1)*dble(3- &
          nspin)
      end do
    end do
  end do
end subroutine
