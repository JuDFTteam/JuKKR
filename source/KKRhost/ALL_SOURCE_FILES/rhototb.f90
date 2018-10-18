!------------------------------------------------------------------------------------
!> Summary: Add core and valence density expanded in spherical harmonics (convention see subroutine rholm )
!> Author: B. Drittler
!> In the paramagnetic case (nspin=1) the core valence charge times
!> \(r^2\) is added to the valence charge density times \(r^2\)
!> then only rho2ns(irmd,lmxtsq,natypd,1) is used .
!> In the spin-polarized case (nspin=2) the spin-splitted core
!> charge density times \(r^2\) is converted into core charge
!> density times \(r^2\) and core spin density times \(r^2\).
!> then these parts are added to corresponding parts of
!> the valence densities times \(r^2\), that are `rho2ns(...,1)`
!> which contains the charge density and `rho2ns(...,2)` which
!> contains in that case the spin density .
!> (see notes by b.drittler)
!------------------------------------------------------------------------------------
!> @note -V. Popescu March 2002: Total orbital moment within the WS sphere is
!> also calculated in the relativistic case; orbital density is normalised in the
!> same way as the charge density.
!> @endnote
!> @warning The core density is spherically averaged and multiplied by \(4\pi\)
!> therefore the core density is only added to l=0 part.
!> @endwarning
!------------------------------------------------------------------------------------
module mod_rhototb

contains

  !-------------------------------------------------------------------------------
  !> Summary: Add core and valence density expanded in spherical harmonics (convention see subroutine rholm )
  !> Author: B. Drittler
  !> Category: core-electrons, physical-observables, KKRhost
  !> Deprecated: False 
  !> In the paramagnetic case (nspin=1) the core valence charge times
  !> \(r^2\) is added to the valence charge density times \(r^2\)
  !> then only `rho2ns(irmd,lmxtsq,natypd,1)` is used .
  !> In the spin-polarized case (nspin=2) the spin-splitted core
  !> charge density times \(r^2\) is converted into core charge
  !> density times \(r^2\) and core spin density times \(r^2\).
  !> then these parts are added to corresponding parts of
  !> the valence densities times \(r^2\), that are `rho2ns(...,1)`
  !> which contains the charge density and `rho2ns(...,2)` which
  !> contains in that case the spin density .
  !> (see notes by b.drittler)
  !-------------------------------------------------------------------------------
  !> @note -V. Popescu March 2002: Total orbital moment within the WS sphere is
  !> also calculated in the relativistic case; orbital density is normalised in the
  !> same way as the charge density.
  !> @endnote
  !> @warning The core density is spherically averaged and multiplied by \(4\pi\)
  !> therefore the core density is only added to l=0 part.
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine rhototb(ipf,natyp,naez,nspin,rho2ns,rhoc,rhoorb,z,drdi,irws,ircut,nfu, &
    llmsp,thetas,ntcell,kshape,ipan,chrgnt,itc,nshell,noq,conc,kaoez,catom,irm,nemb,&
    lmpot)

    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_simpk
    use :: mod_simp3

    implicit none

    ! .. Parameters ..
    ! .. Input variables
    integer, intent (in) :: itc
    integer, intent (in) :: ipf
    integer, intent (in) :: irm    !! Maximum number of radial points
    integer, intent (in) :: nemb   !! Number of 'embedding' positions
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: kshape !! Exact treatment of WS cell
    ! .. Array Arguments ..
    integer, dimension (naez), intent (in) :: noq !! Number of diff. atom types located
    integer, dimension (*), intent (in) :: nfu !! number of shape function components in cell 'icell'
    integer, dimension (*), intent (in) :: ipan !! Number of panels in non-MT-region
    integer, dimension (*), intent (in) :: irws !! R point at WS radius
    integer, dimension (0:nsheld), intent (in) :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0)= number of shells
    integer, dimension (*), intent (in) :: ntcell !! Index for WS cell

    integer, dimension (0:ipand, *), intent (in) :: ircut !! R points of panel borders
    integer, dimension (natyp, *), intent (in) :: llmsp !! lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
    integer, dimension (natyp, naez+nemb), intent (in) :: kaoez !! Kind of atom at site in elem. cell
    real (kind=dp), dimension (*), intent (in) :: z
    real (kind=dp), dimension (natyp), intent (in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension (irm, *), intent (in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension (irm, *), intent (in) :: rhoc !! core charge density
    real (kind=dp), dimension (irm*krel+(1-krel), natyp), intent (in) :: rhoorb !! Orbital density
    real (kind=dp), dimension (irid, nfund, *), intent (in) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion

    ! .. In/Out variables
    real (kind=dp), dimension (irm, lmpot, natyp, *), intent (inout) :: rho2ns
    ! .. Output variables
    real (kind=dp), intent (out) :: chrgnt
    real (kind=dp), dimension (natyp, 2*krel+(1-krel)*nspin), intent (out) :: catom
    ! .. Local variables
    integer :: i, i1, iatyp, icell, ifun, ipan1, ipotd, ipotu, irc1, irs1, ispin, lm, iqez, ioez
    real (kind=dp) :: diff, factor, rfpi, sum, totsmom, totomom, sumo
    real (kind=dp), dimension (natyp) :: omom !! Orbital moment
    real (kind=dp), dimension (irm) :: rho
    real (kind=dp), dimension (naez, 2*krel+(1-krel)*nspin) :: csite
    real (kind=dp), dimension (krel*naez+(1-krel)) :: muosite

    rfpi = sqrt(16.0e0_dp*atan(1.0e0_dp))

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Loop over atomic sites
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do iqez = 1, naez
      ! -------------------------------------------------------------------------
      ! Loop over atoms located on IQEZ
      ! -------------------------------------------------------------------------
      do ispin = 1, nspin
        csite(iqez, ispin) = 0.0e0_dp
      end do

      if (krel==1) muosite(iqez) = 0.0e0_dp

      do ioez = 1, noq(iqez)

        iatyp = kaoez(ioez, iqez)
        ! ---------------------------------------------------------------------
        ! Determine the right potential numbers for rhoc
        ! ---------------------------------------------------------------------
        if (nspin==2) then
          ipotd = 2*iatyp - 1
          ipotu = 2*iatyp
          factor = 1.0e0_dp
        else
          ipotd = iatyp
          ipotu = iatyp
          factor = 0.5e0_dp
        end if

        if (kshape/=0) then
          ipan1 = ipan(iatyp)
          irs1 = ircut(1, iatyp)
          irc1 = ircut(ipan1, iatyp)
        else
          irs1 = irws(iatyp)
          irc1 = irs1
        end if
        ! -----------------------------------------------------------------------
        do i = 2, irs1
          ! -------------------------------------------------------------------
          ! Convert core density
          ! -------------------------------------------------------------------
          sum = (rhoc(i,ipotd)+rhoc(i,ipotu))*factor/rfpi
          diff = (rhoc(i,ipotu)-rhoc(i,ipotd))/rfpi
          ! -------------------------------------------------------------------
          ! Add this to the lm=1 component of rho2ns
          ! -------------------------------------------------------------------
          rho2ns(i, 1, iatyp, 1) = rho2ns(i, 1, iatyp, 1) + sum
          rho2ns(i, 1, iatyp, nspin) = rho2ns(i, 1, iatyp, nspin) + diff
        end do
        ! ----------------------------------------------------------------------
        ! Calculate  charge and moment of the atom
        ! ----------------------------------------------------------------------
        do ispin = 1, nspin

          if (kshape==0) then
            ! ----------------------------------------------------------------
            ! Integrate over wigner seitz sphere - no shape correction
            ! ----------------------------------------------------------------
            call simp3(rho2ns(1,1,iatyp,ispin), sum, 1, irs1, drdi(1,iatyp))
            ! ----------------------------------------------------------------
            ! The result has to be multiplied by sqrt(4 pi)
            ! (4 pi for integration over angle and 1/sqrt(4 pi) for
            ! the spherical harmonic y(l=0))
            ! ----------------------------------------------------------------
            sum = sum*rfpi
          else                     ! (KSHAPE.EQ.0)
            ! ----------------------------------------------------------------
            ! convolute charge density with shape function to get the
            ! charge in the exact cell - if kshape .gt. 0
            ! ----------------------------------------------------------------
            icell = ntcell(iatyp)

            do i = 1, irs1
              rho(i) = rho2ns(i, 1, iatyp, ispin)*rfpi
            end do

            do i = irs1 + 1, irc1
              rho(i) = 0.0e0_dp
            end do

            do ifun = 1, nfu(icell)
              lm = llmsp(icell, ifun)
              if (lm<=lmpot .and. lm>0) then
                do i = irs1 + 1, irc1
                  rho(i) = rho(i) + rho2ns(i, lm, iatyp, ispin)*thetas(i-irs1, ifun, icell)
                end do
              end if
            end do
            ! ----------------------------------------------------------------
            ! Integrate over circumscribed sphere
            ! ----------------------------------------------------------------
            call simpk(rho, sum, ipan1, ircut(0,iatyp), drdi(1,iatyp))
          end if                   ! (KSHAPE.EQ.0)

          catom(iatyp, ispin) = sum
          csite(iqez, ispin) = csite(iqez, ispin) + catom(iatyp, ispin)*conc(iatyp)

          if (ispin/=1) then
            ! ----------------------------------------------------------------
            ! Calculate orbital moment (ASA) and add it to the total
            ! ----------------------------------------------------------------
            if ((krel==1) .and. (kshape==0)) then
              call simp3(rhoorb(1,iatyp), sumo, 1, irs1, drdi(1,iatyp))
              sumo = sumo*rfpi
              omom(iatyp) = sumo
              muosite(iqez) = muosite(iqez) + omom(iatyp)*conc(iatyp)
            end if

            if (kshape/=0) then
              write (ipf, fmt=110) sum
            else
              write (ipf, fmt=130) sum
              if (krel==1) then
                write (ipf, fmt=140) omom(iatyp)
                write (ipf, fmt=150) sum + omom(iatyp)
              end if
            end if
          else                     ! (ISPIN.NE.1)
            if (kshape/=0) then
              write (ipf, fmt=100) iatyp, sum
            else
              write (ipf, fmt=120) iatyp, sum
            end if
          end if                   ! (ISPIN.NE.1)
        end do                     ! ISPIN = 1,NSPIN
        ! ----------------------------------------------------------------------
        if (ioez/=noq(iqez)) write (ipf, '(2X,77("-"))')
      end do
      ! -------------------------------------------------------------------------
      ! IOEZ = 1, NOQ(IQEZ)
      ! -------------------------------------------------------------------------
      if (noq(iqez)>1) then
        write (ipf, '(2X,77("="))')
        write (ipf, fmt=200) iqez, csite(iqez, 1)
        if (nspin==2) then
          write (ipf, fmt=210) csite(iqez, nspin)
          if (krel==1) then
            write (ipf, fmt=220) muosite(iqez)
            write (ipf, fmt=230) csite(iqez, nspin) + muosite(iqez)
          end if
        end if
        if (iqez/=naez) write (ipf, '(2X,77("="))')
      else
        if (iqez/=naez) write (ipf, '(2X,77("="))')
      end if
    end do
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! IQEZ = 1, NAEZ
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write (ipf, *)

    chrgnt = 0.0e0_dp
    do i1 = 1, natyp
      chrgnt = chrgnt + real(nshell(i1), kind=dp)*(catom(i1,1)-z(i1))*conc(i1)
    end do

    write (ipf, '(79("+"))')
    write (ipf, fmt=160) itc, chrgnt
    write (6, fmt=160) itc, chrgnt

    if (nspin==2) then
      totsmom = 0.0e0_dp
      if (krel==1) totomom = 0.0e0_dp
      do i1 = 1, natyp
        totsmom = totsmom + real(nshell(i1), kind=dp)*catom(i1, nspin)*conc(i1)
        if (krel==1) totomom = totomom + real(nshell(i1), kind=dp)*omom(i1)*conc(i1)
      end do

      if (krel==0) then
        write (ipf, fmt=170) totsmom
        write (6, fmt=170) totsmom
      else
        write (ipf, fmt=170) totsmom + totomom
        write (ipf, fmt=180) totsmom
        write (ipf, fmt=190) totomom
        write (6, fmt=170) totsmom + totomom
        write (6, fmt=180) totsmom
        write (6, fmt=190) totomom
      end if
    end if
    write (ipf, *)

    return

100 format ('  Atom ', i4, ' charge in wigner seitz cell =', f10.6)
110 format (7x, 'spin moment in wigner seitz cell =', f10.6)
120 format ('  Atom ', i4, ' charge in wigner seitz sphere =', f10.6)
130 format (7x, 'spin moment in wigner seitz sphere =', f10.6)
140 format (7x, 'orb. moment in wigner seitz sphere =', f10.6)
150 format (7x, 'total magnetic moment in WS sphere =', f10.6)
160 format ('      ITERATION', i4, ' charge neutrality in unit cell = ', f12.6)
170 format ('                   ', ' TOTAL mag. moment in unit cell = ', f12.6)
180 format ('                   ', '           spin magnetic moment = ', f12.6)
190 format ('                   ', '        orbital magnetic moment = ', f12.6)
200 format ('      Site ', i3, ' total charge =', f10.6)
210 format ('         ', ' total spin moment =', f10.6)
220 format ('         ', ' total orb. moment =', f10.6)
230 format ('      total magnetic moment =', f10.6)

  end subroutine rhototb

end module mod_rhototb
