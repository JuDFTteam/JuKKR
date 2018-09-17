module mod_vxcgga

contains

  subroutine vxcgga(exc, kte, kxc, lmax, nspin, iatyp, rho2ns, v, r, drdi, a, irws, ircut, ipan, kshape, gsh, ilm, imaxsh, ifunm, thetas, wtyr, ijend, lmsp, thet, ylm, dylmt1, &
    dylmt2, dylmf1, dylmf2, dylmtf)
    ! -----------------------------------------------------------------------
    ! add the exchange-correlation-potential to the given potential
    ! and if total energies should be calculated (kte=1) the exchange-
    ! correlation-energies are calculated .
    ! use as input the charge density times r**2 (rho2ns(...,1)) and
    ! in the spin-polarized case (nspin=2) the spin density times r**2
    ! (rho2ns(...,2)) .
    ! the density times 4 pi is generated at an angular mesh .
    ! the exchange-correlation potential and the exchange-correlation
    ! energy are calculated at those mesh points with a subroutine .
    ! in the non-spin-polarized case the "spin-density" is
    ! set equal zero .
    ! after that the exchange-correlation potential and in the case of
    ! total energies (kte=1) the exchange-correlation energy are
    ! expanded into spherical harmonics .
    ! the ex.-cor. potential is added to the given potential .
    ! the expansion into spherical harmonics uses the orthogonality
    ! of these harmonics . - therefore a gauss-legendre integration
    ! for "theta" and a gauss-tschebyscheff integration for "phi"
    ! is used .
    ! all needed values for the angular mesh and angular integration
    ! are generate in the subroutine sphere .

    ! the ex.-cor. potential is extrapolated to the origin only
    ! for the lm=1 value .

    ! b.drittler   june 1987

    ! modified for shape functions
    ! b. drittler oct. 1989
    ! simplified and modified for Paragon X/PS
    ! R. Zeller Nov. 1993
    ! -----------------------------------------------------------------------
    ! INCLUDE 'inc.p'
    ! Parameters ..
    ! integer LMPOTD
    ! parameter (LMPOTD= (LPOTD+1)**2)
    ! integer LMXSPD
    ! parameter (LMXSPD= (2*LPOTD+1)**2)
    use :: mod_datatypes, only: dp
    use :: global_variables
    use :: mod_mkxcpe2
    use :: mod_mkxcpe
    use :: mod_gradrl
    use :: mod_simpk
    use :: mod_simp3
    implicit none

    ! Scalar Arguments ..
    real (kind=dp) :: a
    integer :: iatyp, ijend, ipan, irws, kshape, kte, kxc, lmax, nspin

    ! Array Arguments ..
    real (kind=dp) :: drdi(irmd), r(irmd), dylmf1(ijend, lmpotd), dylmf2(ijend, lmpotd), dylmt1(ijend, lmpotd), dylmt2(ijend, lmpotd), dylmtf(ijend, lmpotd), exc(0:lpotd, *), &
      gsh(*), rho2ns(irmd, lmpotd, 2), thet(ijend), wtyr(ijend, *), thetas(irid, nfund), v(irmd, lmpotd, 2), ylm(ijend, lmpotd)
    integer :: ifunm(lmxspd)
    integer :: ilm(ngshd, 3), imaxsh(0:lmpotd), ircut(0:ipand), lmsp(lmxspd)

    ! Local Scalars ..
    real (kind=dp) :: chgden, dx, elmxc, fpi, r1, r2, rpoint, spiden, vlmxc, vxc1, vxc2, vxc3, zero, zero1
    integer :: ifun, ipan1, ipot, ir, irc0, irc1, irh, irs1, ispin, j, l, l1max, lm, lm2, lmmax, m, mesh, nspin2

    ! Local Arrays ..
    real (kind=dp) :: ddrrl(irmd, lmpotd), ddrrul(irmd, lmpotd), drrl(irmd, lmpotd), drrul(irmd, lmpotd), er(irmd, 0:lpotd), estor(irmd, lmpotd), excij(ijend), &
      rhol(irmd, 2, lmpotd), rholm(lmpotd, 2), vxc(ijend, 2), vxcr(2:3, 2)

    ! External Functions ..
    real (kind=dp) :: ddot
    external :: ddot

    ! Data statements ..
    data zero, zero1/0.d0, 1.d-12/

    write (1337, fmt=*) ' GGA CALCULATION '
    fpi = 16.0d0*atan(1.0d0)
    lmmax = (lmax+1)*(lmax+1)

    ! loop over given representive atoms

    if (kshape/=0) then
      ipan1 = ipan
      irc1 = ircut(ipan)
      irs1 = ircut(1)
      irc0 = 2
      if (krel==1) stop ' REL + FULL POTENTIAL N/A '
    else

      irc1 = irws
      irs1 = irc1
      ipan1 = 1
      irc0 = 2
      if (krel==1) irc0 = 2 + mod(ircut(1), 2)
    end if

    do ispin = 1, nspin
      vxcr(2, ispin) = 0.0d0
      vxcr(3, ispin) = 0.0d0
    end do

    ! initialize for ex.-cor. energy

    if (kte==1) then
      do l = 0, lmax
        exc(l, iatyp) = 0.0d0
        do ir = 1, irc1
          er(ir, l) = 0.0d0
        end do
      end do

      do lm = 1, lmmax
        do ir = 1, irc1
          estor(ir, lm) = 0.0d0
        end do
      end do
    end if

    l1max = lmax + 1
    mesh = irws
    dx = a

    if (nspin==2) then
      do lm = 1, lmmax
        do ir = 2, mesh
          r1 = r(ir)
          r2 = r1*r1
          chgden = rho2ns(ir, lm, 1)/r2
          spiden = rho2ns(ir, lm, 2)/r2
          if (abs(chgden)<=zero1) chgden = zero
          if (abs(spiden)<=zero1) spiden = zero
          rhol(ir, 2, lm) = (chgden+spiden)/2.d0
          rhol(ir, 1, lm) = (chgden-spiden)/2.d0
        end do

        ! extrapolate

        rhol(1, 1, lm) = rhol(2, 1, lm)
        rhol(1, 2, lm) = rhol(2, 2, lm)
      end do

    else

      do lm = 1, lmmax
        do ir = 2, mesh
          r1 = r(ir)
          r2 = r1*r1

          chgden = rho2ns(ir, lm, 1)/r2
          if (abs(chgden)<=zero1) chgden = zero
          rhol(ir, 1, lm) = chgden/2.d0
          rhol(ir, 2, lm) = chgden/2.d0
        end do

        ! extrapolate
        rhol(1, 1, lm) = rhol(2, 1, lm)
        rhol(1, 2, lm) = rhol(2, 2, lm)
      end do
    end if


    call gradrl(nspin, mesh, l1max, dx, rhol, r, drdi, ipan1, ipand, ircut, drrl, ddrrl, drrul, ddrrul, irmd, lmpotd)


    ! loop over radial mesh


    do ir = irc0, irc1
      rpoint = r(ir)

      ! calculate the ex.-cor. potential

      nspin2 = 2

      do ispin = 1, nspin2
        do lm = 1, lmmax
          rholm(lm, ispin) = rhol(ir, ispin, lm)
        end do
      end do

      ! only for spin-polarized

      ! PW91 functional
      if (kxc==3) then
        call mkxcpe(nspin2, ir, ijend, l1max, rpoint, rholm, vxc, excij, thet, ylm, dylmt1, dylmt2, dylmf1, dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irmd, lmpotd)
        ! PBE functional
      else if (kxc==4) then
        call mkxcpe2(ir, ijend, rpoint, rholm, vxc, excij, ylm, dylmt1, dylmf1, dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irmd, lmpotd, lmmax, .false.)
        ! PBEsol functional
      else if (kxc==5) then
        call mkxcpe2(ir, ijend, rpoint, rholm, vxc, excij, ylm, dylmt1, dylmf1, dylmf2, dylmtf, drrl, ddrrl, drrul, ddrrul, irmd, lmpotd, lmmax, .true.)
      else
        write (1337, *) ' KXC ???'
        stop
      end if




      ! expand the ex.-cor. potential into spherical harmonics ,
      ! using the orthogonality

      do ispin = 1, nspin

        ! determine the corresponding potential number

        ipot = ispin
        do lm = 1, lmmax
          vlmxc = ddot(ijend, vxc(1,ispin), 1, wtyr(1,lm), 1)
          v(ir, lm, ipot) = v(ir, lm, ipot) + vlmxc

          ! store the ex.-c. potential of ir=2 and =3 for the extrapolation

          if (lm==1 .and. (ir==2 .or. ir==3)) vxcr(ir, ispin) = vlmxc
        end do
      end do

      ! file er in case of total energies

      if (kte==1) then

        ! expand ex.-cor. energy into spherical harmonics
        ! using the orthogonality

        do l = 0, lmax
          do m = -l, l
            lm = l*l + l + m + 1
            elmxc = ddot(ijend, excij, 1, wtyr(1,lm), 1)

            ! multiply the lm-component of the ex.-cor. energy with the same
            ! lm-component of the charge density times r**2 and sum over lm
            ! this corresponds to a integration over the angular .

            if ((kshape/=0) .and. (ir>irs1)) then
              estor(ir, lm) = elmxc

            else

              er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*elmxc
            end if

          end do

        end do

      end if

    end do

    ! integrate er in case of total energies to get exc

    if (kte==1) then
      if (kshape==0) then
        do l = 0, lmax
          call simp3(er(1,l), exc(l,iatyp), 1, irs1, drdi)
        end do

      else

        do l = 0, lmax
          do m = -l, l
            lm = l*l + l + m + 1

            ! convolute with shape function

            do j = imaxsh(lm-1) + 1, imaxsh(lm)
              lm2 = ilm(j, 2)
              if (lmsp(ilm(j,3))>0) then
                ifun = ifunm(ilm(j,3))
                do ir = irs1 + 1, irc1
                  irh = ir - irs1
                  er(ir, l) = er(ir, l) + rho2ns(ir, lm, 1)*gsh(j)*thetas(irh, ifun)*estor(ir, lm2)
                end do
              end if
            end do
          end do
          call simpk(er(1,l), exc(l,iatyp), ipan1, ircut, drdi)
        end do
      end if

    end if

    ! extrapolate ex.-cor potential to the origin only for lm=1

    do ispin = 1, nspin
      ipot = ispin

      vxc2 = vxcr(2, ispin)
      vxc3 = vxcr(3, ispin)
      vxc1 = vxc2 - r(2)*(vxc3-vxc2)/(r(3)-r(2))

      v(1, 1, ipot) = v(1, 1, ipot) + vxc1
    end do

  end subroutine vxcgga

end module mod_vxcgga
