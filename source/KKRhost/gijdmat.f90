!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_gijdmat
  
  private
  public :: gijdmat

contains

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to get the projection matrices to get G_ab from G_CPA
  !> Author: V. Popescu
  !> Date: Oct. 2004
  !> Category: KKRhost, structural-greensfunction, coherent-potential-approx
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Subroutine to get the projection matrices
  !>
  !>               ii     i       i    (-1)
  !>   D  = [ 1 + G   * (t    -  t  ) ]
  !>    a          CPA    CPA     a
  !>
  !>   _            i       i      ii   (-1)
  !>   D  = [ 1 + (t    -  t  ) * G    ]
  !>    a           CPA     a      CPA
  !>
  !> Used is made of the same routine GETDMAT as in the case of TAU    
  !> projection. Note, however, that there the matrices are defined for
  !> example as
  !>
  !>               ij     i       i      (-1)
  !>   D  = [ 1 + tau * (m    -  m    ) ]
  !>    a          CPA    a       CPA
  !>
  !> with m = (t)**(-1)
  !-------------------------------------------------------------------------------
  subroutine gijdmat(tauq, tsst, mssq, dmat, dtil, cfctorinv, iprint, ie, it, krel, lmmaxd)
    use :: mod_datatypes, only: dp
    use :: mod_getdmat, only: getdmat
    use :: mod_cmatstr, only: cmatstr
    use :: mod_constants, only: cone, czero
    implicit none

    integer :: iprint, lmmaxd
    integer :: ie, it, krel
    complex (kind=dp) :: cfctorinv
    complex (kind=dp) :: tauq(lmmaxd, lmmaxd), tsst(lmmaxd, lmmaxd), mssq(lmmaxd, lmmaxd)
    complex (kind=dp) :: dmat(lmmaxd, lmmaxd), dtil(lmmaxd, lmmaxd)

    integer :: ik, info, i1, i2
    integer :: ipvt(lmmaxd)
    complex (kind=dp) :: gll(lmmaxd, lmmaxd), tssq(lmmaxd, lmmaxd), tpg(lmmaxd, lmmaxd), xc(lmmaxd, lmmaxd)
    character (len=18) :: banner

    ! --> get G(CPA) using the same procedure as for GMATLL in < KLOOPZ >
    ! G(CPA) = -MSSQ - MSSQ * TAUQ * MSSQ

    do i2 = 1, lmmaxd
      do i1 = 1, lmmaxd
        tpg(i1, i2) = -cone*tauq(i1, i2)
      end do
    end do

    call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, mssq, lmmaxd, tpg, lmmaxd, czero, xc, lmmaxd)

    call zcopy(lmmaxd*lmmaxd, mssq, 1, gll, 1)

    call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, -cone, xc, lmmaxd, mssq, lmmaxd, -cone, gll, lmmaxd)

    call zscal(lmmaxd*lmmaxd, cfctorinv, gll, 1)

    ! --> invert MSSQ to get TSSQ

    do i2 = 1, lmmaxd
      do i1 = 1, lmmaxd
        tssq(i1, i2) = mssq(i1, i2)*cfctorinv
      end do
    end do
    call zgetrf(lmmaxd, lmmaxd, tssq, lmmaxd, ipvt, info)
    call zgetri(lmmaxd, tssq, lmmaxd, ipvt, xc, lmmaxd*lmmaxd, info)

    ! --> get projection matrices from G(CPA)
    ! call getdmat with t,c for the local arg. c,t

    call getdmat(gll, dmat, dtil, xc, lmmaxd, tsst, tssq, lmmaxd)

    if (iprint<=1) return

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    ik = 2*krel + 1
    write (banner, '("DMAT IE=",I3," IT=",I3)') ie, it
    call cmatstr(banner, 18, dmat, lmmaxd, lmmaxd, ik, ik, 0, 1e-10_dp, 6)

    write (banner, '("DTIL IE=",I3," IT=",I3)') ie, it
    call cmatstr(banner, 18, dtil, lmmaxd, lmmaxd, ik, ik, 0, 1e-10_dp, 6)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

  end subroutine gijdmat

end module mod_gijdmat
