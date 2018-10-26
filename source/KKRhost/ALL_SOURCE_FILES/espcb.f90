!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_espcb

contains

  !-------------------------------------------------------------------------------
  !> Summary: Collects single-particle core energy
  !> Author: B. Drittler
  !> Category: KKRhost, total-energy
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Attention : energy zero ---> electro static zero
  !>
  !> Since input potential and single particle energies
  !> are using muffin tin zero as zero the energy shift
  !> is cancelled in the kinetic energy contribution!
  !>
  !> Calculate the core contribution of the single particle energies
  !> l and spin dependent.
  !> Attention : here are the results of the subroutine corel (stored
  !> in the common block core) used.
  !> (see notes by B. Drittler)
  !>
  !> modified for bandstructure code
  !> B. Drittler   Jan 1990
  !-------------------------------------------------------------------------------
  subroutine espcb(espc, nspin, natyp, ecore, lcore, lcoremax, ncore)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: npotd
    implicit none

    ! .. Local Scalars ..
    integer :: natyp, nspin

    ! .. Intrinsic Functions ..
    real (kind=dp) :: ecore(20, *), espc(0:3, npotd)
    integer :: lcore(20, *), ncore(*), lcoremax(*)

    integer :: i1, ipot, is, l, n


    ! loop over reference atoms
    do i1 = 1, natyp
      lcoremax(i1) = 0
      do is = 1, nspin
        ! determine correct potential indices
        ipot = nspin*(i1-1) + is

        ! initialize espc
        do l = 0, 3
          espc(l, ipot) = 0.0e0_dp
        end do

        ! loop over all core states
        do n = 1, ncore(ipot)
          l = lcore(n, ipot)
          lcoremax(i1) = max(lcoremax(i1), l)
          write (1337, *) 'in espcb:', n, l, ipot, i1
          espc(l, ipot) = espc(l, ipot) + ecore(n, ipot)*real(2*l+1, kind=dp)*real(3-nspin, kind=dp)
        end do

      end do
    end do

  end subroutine espcb

end module mod_espcb
