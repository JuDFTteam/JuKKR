!------------------------------------------------------------------------------------
!> Summary: Calculate the electron-intracell-potentials and the charge-moments 
!> of given charge densities. ( For each spin-direction the potential is the
!> same in the polarized case.)
!> Author: B. Drittler
!> Initialize the potential $$ V$$ with the electron-intracell-potentials
!> the intracell-potential is expanded into spherical harmonics .
!> the lm-term of the intracell-potential of the representive atom i is given by
!>
!> $$V\left(r,lm,i\right)=\frac{8\pi}{2l+1}\left(\int_{0}^{r}dr'\frac{r'^{l}}{r^{l+1}} rho2ns(r',lm,i,1) +\int_{r}^{r_{cut}}dr' \frac{r^l}{r'^{l+1}}rho2ns(r',lm,i,1) ) \right)$$
!>
!> the lm contribution of the charge moment of the representive atom i is given by
!>
!> $$ cmom\left(lm,i\right)=\int_{0}^{r_{cut}} dr'!r'^{l}rho2ns(r',lm,i,1)$$
!>
!> (see notes by b.drittler and u.klemradt) $$ r_{cut}$$ is muffin tin or
!> Wigner-Seitz sphere radius, depending on kshape turned on or off
!> @note Attention : $$ rho2ns(...,1)$$ is the real charge density times $$r^2$$
!> developed into spherical harmonics . (see deck rholm)
!------------------------------------------------------------------------------------
module mod_vintras

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the electron-intracell-potentials and the charge-moments 
  !> of given charge densities. ( For each spin-direction the potential is the
  !> same in the polarized case.)
  !> Author: B. Drittler
  !> Category: potential, KKRhost
  !> Deprecated: False 
  !> Initialize the potential $$ V$$ with the electron-intracell-potentials
  !> the intracell-potential is expanded into spherical harmonics .
  !> the lm-term of the intracell-potential of the representive atom i is given by
  !>
  !> $$V\left(r,lm,i\right)=\frac{8\pi}{2l+1}\left(\int_{0}^{r}dr'\frac{r'^{l}}{r^{l+1}} rho2ns(r',lm,i,1) +\int_{r}^{r_{cut}}dr' \frac{r^l}{r'^{l+1}}rho2ns(r',lm,i,1) ) \right)$$
  !>
  !> the lm contribution of the charge moment of the representive atom i is given by
  !>
  !> $$ cmom\left(lm,i\right)=\int_{0}^{r_{cut}} dr'!r'^{l}rho2ns(r',lm,i,1)$$
  !>
  !> (see notes by b.drittler and u.klemradt) $$ r_{cut}$$ is muffin tin or
  !> Wigner-Seitz sphere radius, depending on kshape turned on or off
  !> @note Attention : $$ rho2ns(...,1)$$ is the real charge density times $$r^2$$
  !> developed into spherical harmonics . (see deck rholm)
  !-------------------------------------------------------------------------------
  subroutine vintras(cmom,cminst,lmax,nspin,nstart,nend,rho2ns,v,r,drdi,irws,ircut, &
    ipan,kshape,ntcell,ilm_map,ifunm,imaxsh,gsh,thetas,lmsp,lmpot,natyp)

    use :: constants
    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_sinwk
    use :: mod_soutk

    implicit none

    ! .. Input Variables
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: nend
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nstart
    integer, intent (in) :: kshape !! Exact treatment of WS cell
    integer, dimension (natyp), intent (in) :: irws !! R point at WS radius
    integer, dimension (natyp), intent (in) :: ipan !! Number of panels in non-MT-region
    integer, dimension (natyp), intent (in) :: ntcell !! Index for WS cell
    integer, dimension (0:lmpot), intent (in) :: imaxsh
    integer, dimension (ngshd, 3), intent (in) :: ilm_map
    integer, dimension (natyp, lmxspd), intent (in) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (natyp, lmxspd), intent (in) :: ifunm
    integer, dimension (0:ipand, natyp), intent (in) :: ircut !! R points of panel borders
    real (kind=dp), dimension (ngshd), intent (in) :: gsh
    real (kind=dp), dimension (irmd, natyp), intent (in) :: r !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (irmd, natyp), intent (in) :: drdi !! Derivative dr/di
    real (kind=dp), dimension (irid, nfund, ncelld), intent (in) :: thetas !! shape function THETA=0 outer space THETA=1 inside WS cell in spherical harmonics expansion
    real (kind=dp), dimension (irmd, lmpot, natyp, 2), intent (in) :: rho2ns !! radial density
    ! .. Output variables
    real (kind=dp), dimension (lmpot, natyp), intent (out) :: cmom !! LM moment of total charge
    real (kind=dp), dimension (lmpot, natyp), intent (out) :: cminst
    real (kind=dp), dimension (irmd, lmpot, npotd), intent (out) :: v
    ! .. Local Variables
    real (kind=dp) :: fac, rl
    integer :: i, iatyp, icell, iend, ifun, ipot, irc1, irs1, istart, j, l, lm, lm2, lm3, m
    ! .. Local Arrays
    integer, dimension (0:ipand) :: ircutm
    real (kind=dp), dimension (irmd) :: v1
    real (kind=dp), dimension (irmd) :: v2
    real (kind=dp), dimension (irmd) :: vint1
    real (kind=dp), dimension (irmd) :: vint2
    ! ..
    do iatyp = nstart, nend
      if (kshape/=0) then
        irs1 = ircut(1, iatyp)
        irc1 = ircut(ipan(iatyp), iatyp)
        icell = ntcell(iatyp)
        do i = 0, ipan(iatyp)
          ircutm(i) = ircut(i, iatyp)
        end do                     ! I
      else
        irs1 = irws(iatyp)
        irc1 = irs1
        ircutm(0) = ircut(0, iatyp)
        ircutm(1) = irc1
      end if
      ! -------------------------------------------------------------------------
      ! Determine the right potential numbers
      ! -------------------------------------------------------------------------
      ipot = nspin*iatyp

      do l = 0, lmax
        fac = 8.0e0_dp*pi/real(2*l+1, kind=dp)
        do m = -l, l
          lm = l*l + l + m + 1
          ! -------------------------------------------------------------------
          ! Set up of the integrands v1 and v2
          ! -------------------------------------------------------------------
          v1(1) = 0.0e0_dp
          v2(1) = 0.0e0_dp
          do i = 2, irs1
            rl = r(i, iatyp)**l
            v1(i) = rho2ns(i, lm, iatyp, 1)*rl*drdi(i, iatyp)
            v2(i) = rho2ns(i, lm, iatyp, 1)/r(i, iatyp)/rl*drdi(i, iatyp)
          end do                   ! I
          ! -------------------------------------------------------------------
          ! Convolute charge density of interstial with shape function if
          ! kshape.gt.0
          ! -------------------------------------------------------------------
          if (kshape/=0) then
            do i = irs1 + 1, irc1
              v1(i) = 0.0e0_dp
            end do                 ! I
            istart = imaxsh(lm-1) + 1
            iend = imaxsh(lm)
            do j = istart, iend
              lm2 = ilm_map(j, 2)
              lm3 = ilm_map(j, 3)
              if (lmsp(icell,lm3)>0) then
                ifun = ifunm(icell, lm3)
                do i = irs1 + 1, irc1
                  v1(i) = v1(i) + gsh(j)*rho2ns(i, lm2, iatyp, 1)*thetas(i-irs1, ifun, icell)
                end do             ! I
              end if
            end do                 ! J

            do i = irs1 + 1, irc1
              rl = r(i, iatyp)**l
              v2(i) = v1(i)/r(i, iatyp)/rl*drdi(i, iatyp)
              v1(i) = v1(i)*rl*drdi(i, iatyp)
            end do                 ! I
          end if
          ! -------------------------------------------------------------------
          ! Now integrate v1 and v2
          ! -------------------------------------------------------------------
          call soutk(v1, vint1, ipan(iatyp), ircutm)
          call sinwk(v2, vint2, ipan(iatyp), ircutm)
          ! -------------------------------------------------------------------
          ! Gather all parts
          ! -------------------------------------------------------------------
          if (lm==1) then
            v(1, lm, ipot) = fac*vint2(1)
          else
            v(1, lm, ipot) = 0.0e0_dp
          end if

          do i = 2, irc1
            rl = r(i, iatyp)**l
            v(i, lm, ipot) = fac*(vint1(i)/r(i,iatyp)/rl+vint2(i)*rl)
          end do                   ! I
          ! -------------------------------------------------------------------
          ! Store charge moment - in case of kshape.gt.0 this is the moment
          ! of the charge in the muffin tin sphere
          ! -------------------------------------------------------------------
          cmom(lm, iatyp) = vint1(irs1)
          ! -------------------------------------------------------------------
          ! Store charge moment of interstial in case of kshape.gt.0
          ! -------------------------------------------------------------------
          if (kshape/=0) cminst(lm, iatyp) = vint1(irc1) - vint1(irs1)

          if (nspin==2) then
            do i = 1, irc1
              v(i, lm, ipot-1) = v(i, lm, ipot)
            end do                 ! I
          end if
        end do                     ! M
      end do                       ! L
    end do                         ! IATYP
    return

  end subroutine vintras

end module mod_vintras
