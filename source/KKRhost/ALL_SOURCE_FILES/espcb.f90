! 13.10.95 ***************************************************************
subroutine espcb(espc, nspin, natyp, ecore, lcore, lcoremax, ncore)
  ! ************************************************************************

  ! attention : energy zero ---> electro static zero

  ! since input potential and single particle energies
  ! are using muffin tin zero as zero the energy shift
  ! is cancelled in the kinetic energy contribution !

  ! calculate the core contribution of the single particle energies
  ! l and spin dependent .
  ! attention : here are the results of the subroutine corel (stored
  ! in the common block core) used .
  ! (see notes by b.drittler)

  ! modified for bandstructure code
  ! b.drittler   jan 1990
  ! -----------------------------------------------------------------------

  use :: mod_datatypes, only: dp
  use global_variables
  implicit none

  ! .. Local Scalars ..
  integer :: natyp, nspin

  ! .. Intrinsic Functions ..
  real (kind=dp) :: ecore(20, *), espc(0:3, npotd)
  integer :: lcore(20, *), ncore(*), lcoremax(*)

  ! loop over reference atoms
  integer :: i1, ipot, is, l, n


  intrinsic :: real
  ! determine correct potential indices

  do i1 = 1, natyp
    lcoremax(i1) = 0
    do is = 1, nspin
      ! initialize espc

      ! loop over all core states
      ipot = nspin*(i1-1) + is
      ! 13.10.95
      ! ***************************************************************
      ! ************************************************************************
      do l = 0, 3
        espc(l, ipot) = 0.0e0_dp
      end do

      ! attention : energy zero ---> electro static zero
      do n = 1, ncore(ipot)
        l = lcore(n, ipot)
        lcoremax(i1) = max(lcoremax(i1), l)
        write (1337, *) 'in espcb:', n, l, ipot, i1
        espc(l, ipot) = espc(l, ipot) + ecore(n, ipot)*real(2*l+1, kind=dp)* &
          real(3-nspin, kind=dp)
      end do
    end do
  end do
end subroutine espcb
