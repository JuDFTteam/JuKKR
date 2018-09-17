module mod_vxcdrv

contains

  subroutine vxcdrv(exc, kte, kxc, lpot, nspin, nstart, nend, rho2ns, vons, r, drdi, a, irws, ircut, ipan, ntcell, kshape, gsh, ilm, imaxsh, ifunm, thetas, lmsp)
    use :: global_variables

    use :: mod_datatypes, only: dp
    use :: mod_sphere_nogga
    use :: mod_sphere_gga
    use :: mod_vxcgga
    use :: mod_vxclm
    implicit none
    ! INCLUDE 'inc.p'
    ! Parameters ..
    integer :: ijd                 ! ,LMPOTD,LMXSPD
    ! parameter (LMPOTD= (LPOTD+1)**2,LMXSPD= (2*LPOTD+1)**2)
    parameter (ijd=434)

    ! Scalar Arguments ..
    integer :: kshape, kte, kxc, lpot, nend, nspin, nstart

    ! Array Arguments ..
    real (kind=dp) :: a(natypd), drdi(irmd, *), exc(0:lpotd, *), gsh(*), r(irmd, *), rho2ns(irmd, lmpotd, natypd, *), thetas(irid, nfund, *), vons(irmd, lmpotd, *)
    integer :: ifunm(natypd, *), ilm(ngshd, 3), imaxsh(0:lmpotd), ipan(*), ircut(0:ipand, *), irws(*), lmsp(natypd, *), ntcell(*)

    ! Local Arrays ..
    real (kind=dp) :: dylmf1(ijd, lmpotd), dylmf2(ijd, lmpotd), dylmt1(ijd, lmpotd), dylmt2(ijd, lmpotd), dylmtf(ijd, lmpotd), rho2iat(irmd, lmpotd, 2), rij(ijd, 3), thet(ijd), &
      wtyr(ijd, lmpotd), ylm(ijd, lmpotd), yr(ijd, lmpotd)
    integer :: ifunmiat(lmxspd), lmspiat(lmxspd)

    ! Local Scalars ..
    integer :: iatyp, icell, ipot, lmx1

    if (kxc<3) then
      call sphere_nogga(lpot, yr, wtyr, rij, ijd)
    else
      call sphere_gga(lpot, yr, wtyr, rij, ijd, lmpotd, thet, ylm, dylmt1, dylmt2, dylmf1, dylmf2, dylmtf)
    end if
    do iatyp = nstart, nend
      icell = ntcell(iatyp)
      ipot = nspin*(iatyp-1) + 1
      do lmx1 = 1, lmxspd
        ifunmiat(lmx1) = ifunm(icell, lmx1)
        lmspiat(lmx1) = lmsp(icell, lmx1)
      end do
      call dcopy(irmd*lmpotd, rho2ns(1,1,iatyp,1), 1, rho2iat(1,1,1), 1)
      if (nspin==2 .or. krel==1) then
        call dcopy(irmd*lmpotd, rho2ns(1,1,iatyp,2), 1, rho2iat(1,1,2), 1)
      end if
      if (kxc<3) then
        call vxclm(exc, kte, kxc, lpot, nspin, iatyp, rho2iat, vons(1,1,ipot), r(1,iatyp), drdi(1,iatyp), irws(iatyp), ircut(0,iatyp), ipan(iatyp), kshape, gsh, ilm, imaxsh, &
          ifunmiat, thetas(1,1,icell), yr, wtyr, ijd, lmspiat)
      else

        ! GGA EX-COR POTENTIAL

        call vxcgga(exc, kte, kxc, lpot, nspin, iatyp, rho2iat, vons(1,1,ipot), r(1,iatyp), drdi(1,iatyp), a(iatyp), irws(iatyp), ircut(0,iatyp), ipan(iatyp), kshape, gsh, ilm, &
          imaxsh, ifunmiat, thetas(1,1,icell), wtyr, ijd, lmspiat, thet, ylm, dylmt1, dylmt2, dylmf1, dylmf2, dylmtf)
      end if
    end do
  end subroutine vxcdrv

end module mod_vxcdrv
