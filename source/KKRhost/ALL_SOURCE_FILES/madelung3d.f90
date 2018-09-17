module mod_madelung3d

contains

  subroutine madelung3d(lpot, yrg, wg, naez, alat, volume0, bravais, recbv, rbasis, rmax, gmax, naezd, lmxspd, lassld, lpotd, lmpotd, nmaxd, ishld, nembd, wlength)
    ! **********************************************************************
    ! *                                                                    *
    ! * This subroutine calculates the Madelung potential coefficients     *
    ! * in the 3D case and stores them in the DA-file abvmad.unformatted   *
    ! * The record index is simply (IQ1-1)*NAEZ + IQ2 for record (IQ1,IQ2) *
    ! *                                                                    *
    ! **********************************************************************
    use :: mod_datatypes, only: dp
    use :: mod_madelgaunt
    use :: mod_madelcoef
    use :: mod_madelout, only: madel3out
    use :: mod_lattice3d
    use :: mod_strmat

    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer :: lpot, naez, wlength
    integer :: naezd, lmxspd, lassld, lpotd, lmpotd, nmaxd, ishld, nembd
    real (kind=dp) :: alat, volume0, rmax, gmax
    ! ..
    ! .. Array Arguments ..
    real (kind=dp) :: yrg(lassld, 0:lassld, 0:lassld), wg(lassld)
    real (kind=dp) :: bravais(3, 3), recbv(3, 3)
    real (kind=dp) :: rbasis(3, naezd+nembd)
    ! ..
    ! .. Local Scalars ..
    integer :: iend, iprint, iq1, iq2, nclebd
    integer :: ngmax, nrmax, nshlg, nshlr
    integer :: lrecabmad, irec
    ! ..
    ! .. Local Arrays ..
    ! .. Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
    real (kind=dp) :: avmad(lmpotd, lmpotd), bvmad(lmpotd)
    real (kind=dp) :: cleb(lmxspd*lmpotd)
    real (kind=dp) :: madelsmat(lmxspd, naezd, naezd)
    real (kind=dp) :: smat1(6, 6), smat2(6, 6)
    real (kind=dp) :: gn(3, nmaxd), rm(3, nmaxd)
    integer :: icleb(lmxspd*lmpotd, 3)
    integer :: nsg(ishld), nsr(ishld)
    ! ......................................................................
    iprint = 0
    nclebd = lmxspd*lmpotd

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    write (1337, '(79("="))')
    write (1337, '(18X,A)') 'MADELUNG3D: setting bulk Madelung coefficients'
    write (1337, '(79("="))')
    write (1337, *)
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    ! ======================================================================
    call lattice3d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsg, nsr, gn, rm, rmax, gmax, iprint, nmaxd, ishld)

    call strmat(alat, lpot, naez, ngmax, nrmax, nsg, nsr, nshlg, nshlr, gn, rm, rbasis, madelsmat, volume0, iprint, lassld, lmxspd, naezd)
    ! ======================================================================

    lrecabmad = wlength*2*lmpotd*lmpotd + wlength*2*lmpotd
    ! lrecabmad = wlength*kind(0.0_dp)*lmpotd*lmpotd + wlength*kind(0.0_dp)*lmpotd
    open (69, access='direct', recl=lrecabmad, file='abvmad.unformatted', form='unformatted')

    ! --> calculate the gaunt coefficients

    call madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)

    ! --> calculate the madelung coefficients to be used for VMAD
    ! call MADELCOEF with first arg. .FALSE. = 3D case

    do iq1 = 1, naez
      do iq2 = 1, naez
        call madelcoef(.false., lpot, avmad, bvmad, madelsmat(1,iq1,iq2), cleb, icleb, iend, lpotd, lmpotd, lmxspd, nclebd)

        irec = iq2 + naez*(iq1-1)
        write (69, rec=irec) avmad, bvmad
        ! -----------------------------------------------------------------------
        if ((iq1<=6) .and. (iq2<=6)) then
          smat1(iq1, iq2) = avmad(1, 1)
          smat2(iq1, iq2) = bvmad(1)
        end if
        ! -----------------------------------------------------------------------
      end do
    end do
    close (69)

    if (iprint<1) return
    ! ======================================================================

    call madel3out(iprint, naez, lrecabmad, smat1, smat2, lmpotd)

  end subroutine madelung3d

end module mod_madelung3d
