module mod_vintras_kkrimp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the electron-intracell-potentials and the charge-moments of given charge densities
  !> Author: B. Drittler, U. Klemradt
  !> Date: May 1987
  !> Category: KKRimp, potential
  !> Deprecated: False ! this needs to be set to true for deprecated subroutines
  !>
  !> Calculate the electron-intracell-potentials and the charge-
  !> moments of given charge densities. ( for each spin-direc-
  !> tion the potential is the same in the polarized case. )
  !> initialize the potential v with the electron-intracell-potentials
  !> the intracell-potential is expanded into spherical harmonics.
  !> the lm-term of the intracell-potential of the representive atom i
  !> is given by
  !>                8pi        r      r'** l
  !>  v(r,lm,i) =  ----- *  (  s dr' --------   rho2ns(r',lm,i,1)
  !>               2*l+1       0     r **(l+1)
  !
  !>                             rcut    r ** l
  !>                           +  s dr' ---------   rho2ns(r',lm,i,1) )
  !>                              r     r' **(l+1)
  !
  !> the lm contribution of the charge moment of the representive
  !> atom i is given by
  !
  !>                         rcut
  !>          cmom(lm,i) =    s dr' r'** l rho2ns(r',lm,i,1)
  !>                          0
  !
  !>         (see notes by B. Drittler and U. Klemradt)
  !
  !>          rcut is muffin tin or wigner seitz sphere radius,
  !>          depending on kshape turned on or off
  !
  !> @warning rho2ns(...,1) is the real charge density times r**2
  !> developed into spherical harmonics. (see deck rholm) @endwarning
  !-------------------------------------------------------------------------------
  subroutine vintras(natom, nspin, nrmaxd,lmaxd, lmaxatom,cell, vpot_out, shapefun, gauntshape, density,cmom,cmom_interst,ins)

  use mod_sinwk, only: sinwk
  use mod_soutk, only: soutk
  use type_cell, only: cell_type
  use type_shapefun, only: shapefun_type
  use type_gauntshape, only: gauntshape_type
  use type_density, only: density_type
  implicit none

  integer :: natom
  integer :: nspin
  integer :: nrmaxd
  integer :: lmaxd
  integer :: lmaxatom(natom)
  real*8,allocatable  :: cmom(:,:)
  real*8,allocatable  :: cmom_interst(:,:)
  real*8                ::  vpot_out(nrmaxd,(2*lmaxd+1)**2,nspin,natom) !thetas(iri,nfund,*),

  type(cell_type)                     :: cell(natom)
  type(shapefun_type)                 :: shapefun(natom)
  type(gauntshape_type)               :: gauntshape(lmaxd)
  type(density_type)                  :: density(natom)

  real*8 fac,pi,rl


  integer i,iatom,iend,ifun,irc1,irs1,istart,j,lval,lm,lm2,lm3,mval,lmax
  real*8 v1(nrmaxd),v2(nrmaxd),vint1(nrmaxd),vint2(nrmaxd)
  integer ins
  integer ipand
  integer,allocatable :: ircutm(:)


  ! ######################################################
  ! calculate the maximum number of panels
  ! ######################################################
  ipand=0
  do iatom=1,natom
    if (cell(iatom)%npan> ipand) ipand=cell(iatom)%npan
  end do !natom

  allocate(ircutm(0:ipand))

  pi = 4.d0*datan(1.d0)

  do iatom = 1,natom

    if (ins.ne.0) then
    irs1 = cell(iatom)%nrcut(      1           ) !ircut(1,iatom)

    irc1 = cell(iatom)%nrcut( cell(iatom)%npan ) !ircut(ipan(iatom),iatom)

    ! error in here
    !          do 10 i = 0,ipan(icell)
    ! next line is correct one (t.korhonen, nov 1997)
    ircutm=0
    do i = 0, cell(iatom)%npan  !ipan(iatom)
    ircutm(i) = cell(iatom)%nrcut(i) !ircut(i,iatom)
    end do

    else
    irs1 = cell(iatom)%nrmax !irws(iatom) ????
    irc1 = irs1
    ircutm(0) = cell(iatom)%nrcut(0) !ircut(0,iatom)
    ircutm(1) = irc1
    end if
    !---> determine the right potential numbers
    !        ipot = nspin*iatom

    lmax=lmaxatom(iatom)
    do lval = 0,2*lmaxatom(iatom)
      v1=0.0d0
      v2=0.0d0
      vint1=0.0d0
      vint2=0.0d0

      fac = 8.0d0*pi/real(2*lval+1)
      do mval = -lval,lval
        lm = lval*lval + lval + mval + 1
        
        !---> set up of the integrands v1 and v2
        v1(1) = 0.0d0
        v2(1) = 0.0d0
        do i = 2,irs1
          rl = cell(iatom)%rmesh(i)**lval
          v1(i) = density(iatom)%rho2ns(i,lm,1)*rl*cell(iatom)%drmeshdi(i)
          v2(i) = density(iatom)%rho2ns(i,lm,1)/cell(iatom)%rmesh(i)/rl*cell(iatom)%drmeshdi(i)
        end do ! i

        !---> convolute charge density of interstial with shape function
        !        if ins.gt.0
        if (ins.ne.0) then
          do i = irs1 + 1,irc1
            v1(i) = 0.0d0
          end do !i
          istart = gauntshape(lmax)%imaxsh(lm-1) + 1
          iend   = gauntshape(lmax)%imaxsh(lm)
          do j = istart,iend
            lm2 = gauntshape(lmax)%ilm(j,2)
            lm3 = gauntshape(lmax)%ilm(j,3)
            if (shapefun(iatom)%lmused(lm3).gt.0) then !lmsp(icell,lm3).gt.0
              ifun = shapefun(iatom)%lm2index(lm3) !ifunm(icell,lm3)
              do i = irs1 + 1,irc1
                v1(i) = v1(i) + gauntshape(lmax)%gsh(j)*density(iatom)%rho2ns(i,lm2,1)* shapefun(iatom)%thetas(i-irs1,ifun)
              end do !i
            end if
          end do !j

          do i = irs1 + 1,irc1
            rl = cell(iatom)%rmesh(i)**lval
            v2(i) = v1(i)/cell(iatom)%rmesh(i)/rl*cell(iatom)%drmeshdi(i)
            v1(i) = v1(i)*rl*cell(iatom)%drmeshdi(i)
          end do !i
        end if

        !---> now integrate v1 and v2
        call soutk(v1,vint1,cell(iatom)%npan,ircutm)

        call sinwk(v2,vint2,cell(iatom)%npan,ircutm)

        !---> gather all parts
        if (lm.eq.1) then
          vpot_out(1,lm,1,iatom) = fac*vint2(1)
        else
          vpot_out(1,lm,1,iatom) = 0.0d0
        end if

        do i = 2,irc1
          rl = cell(iatom)%rmesh(i)**lval
          vpot_out(i,lm,1,iatom) = fac* (vint1(i)/cell(iatom)%rmesh(i)/rl+vint2(i)*rl)
        end do

        !---> store charge moment - in case of kshape.gt.0 this is the moment
        !      of the charge in the muffin tin sphere
        !            cmom(lm,iatyp) = vint1(irs1)


        !---> store charge moment of interstial in case of kshape.gt.0

        ! if (kshape.ne.0) cminst(lm,iatom) = vint1(irc1) - vint1(irs1)
        if (ins.ne.0) then
          cmom_interst(lm,iatom) = vint1(irc1) - vint1(irs1)
        end if
        ! else
        cmom(lm,iatom) = vint1(irc1)
        ! end if

        if (nspin.eq.2) then
          do i = 1,irc1
            vpot_out(i,lm,2,iatom) = vpot_out(i,lm,1,iatom)
          end do  
        end if

      end do !m

    end do !l

  end do !iatom

  end subroutine vintras

end module mod_vintras_kkrimp
