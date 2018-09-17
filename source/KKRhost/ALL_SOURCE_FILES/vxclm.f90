module mod_vxclm

contains

  subroutine vxclm(exc, kte, kxc, lmax, nspin, iatyp, rho2ns, v, r, drdi, irws, ircut, ipan, kshape, gsh, ilm, imaxsh, ifunm, thetas, yr, wtyr, ijend, lmsp)
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
    ! in the paramagnetic case the "spin-density" is set equal zero .
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
    ! cor error 23/6/1996
    ! -----------------------------------------------------------------------
    ! INCLUDE 'inc.p'
    ! .. Parameters ..
    ! integer LMPOTD
    ! parameter (LMPOTD= (LPOTD+1)**2)
    ! integer LMXSPD
    ! parameter (LMXSPD= (2*LPOTD+1)**2)
    use :: mod_datatypes, only: dp
    use :: global_variables
    use :: mod_vosko
    use :: mod_vxcspo
    use :: mod_simpk
    use :: mod_simp3
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    integer :: iatyp, ijend, ipan, irws, kshape, kte, kxc, lmax, nspin
    ! ..
    ! .. Array Arguments ..
    real (kind=dp) :: drdi(irmd), exc(0:lpotd, *), gsh(*), r(irmd), rho2ns(irmd, lmpotd, 2), thetas(irid, nfund), v(irmd, lmpotd, 2), wtyr(ijend, *), yr(ijend, *)
    integer :: ifunm(lmxspd)
    integer :: ilm(ngshd, 3), imaxsh(0:lmpotd), ircut(0:ipand), lmsp(lmxspd)
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: elmxc, fpi, fpipr2, vlmxc, vxc1, vxc2, vxc3, factor
    integer :: ifun, ij, ipot, ir, irc1, irh, irs1, is, ispin, j, l, lm, lm2, lmmax, m
    ! ..
    ! .. Local Arrays ..
    real (kind=dp) :: er(irmd, 0:lpotd), estor(irmd, lmpotd), excij(ijend), fprho(ijend, 2), vxc(ijend, 2), vxcr(2:3, 2)
    ! ..
    ! .. External Functions ..
    real (kind=dp) :: ddot
    external :: ddot

    write (1337, *) 'Including cutoff of vxc for small density'
    fpi = 16.0d0*atan(1.0d0)
    lmmax = (lmax+1)*(lmax+1)


    ! loop over given representive atoms

    if (kshape/=0) then
      irc1 = ircut(ipan)
      irs1 = ircut(1)

    else

      irc1 = irws
      irs1 = irc1
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

    ! loop over radial mesh


    do ir = 2, irc1


      ! generate the densities on an angular mesh

      do is = 1, 2
        do ij = 1, ijend
          fprho(ij, is) = 0.d0
        end do
      end do

      fpipr2 = fpi/r(ir)**2
      do ispin = 1, nspin
        do lm = 1, lmmax
          call daxpy(ijend, rho2ns(ir,lm,ispin)*fpipr2, yr(1,lm), 1, fprho(1,ispin), 1)
        end do
      end do

      ! calculate the ex.-cor. potential

      if (kxc<=1) then
        call vxcspo(excij, fprho, vxc, kxc, ijend, ijend)
      else
        call vosko(excij, fprho, vxc, ijend, ijend)
      end if

      do ij = 1, ijend
        factor = (1.d0-exp(-abs(fprho(ij,1))*1000.d0))
        do ispin = 1, nspin
          vxc(ij, ispin) = vxc(ij, ispin)*factor ! cutoff
        end do
      end do


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
          call simpk(er(1,l), exc(l,iatyp), ipan, ircut, drdi)
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

  end subroutine vxclm

end module mod_vxclm
