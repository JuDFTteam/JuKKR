
!------------------------------------------------------------------------------
!> Module with routines for alternative calculation of the
!> total energy.
!
! Weinert: -(E_coulomb + E_Z_interstitial) + E_madelung - 2*E_VXC + E_single_particle + E_core + E_XC + E_LDAU...
!
! E_VXC = 1/2 \int V_XC \rho
module total_energy_mod
#include "macros.h"
  use Constants_mod, only : pi
  implicit none
  private
  public :: madelung_energy, madelung_ref_radius_correction, energy_electrostatic_wrapper, energy_electrostatic_L_resolved_wrapper
  
  contains

  !----------------------------------------------------------------------------
  !> Calculates electrostatic energy 'ecou' for arbitrary potential 'v_potential'
  !> and charge density 'rho2ns'.
  !>
  !> 1/2 \int \rho V
  double precision function energy_electrostatic_wrapper(v_potential, Z_nuclear, RHO2NS, shgaunts, atomdata) result(ecou)
    use BasisAtom_mod, only: BasisAtom
    use ShapeGauntCoefficients_mod, only: ShapeGauntCoefficients

    double precision, intent(in) :: v_potential(:,:,:)
    double precision, intent(in) :: Z_nuclear
    double precision, intent(in) :: RHO2NS(:,:,:)
    type(BasisAtom), intent(in) :: atomdata
    type(ShapeGauntCoefficients), intent(in) :: shgaunts

    integer :: lmax_potential
    double precision, allocatable :: energy(:)

    lmax_potential = int(sqrt(size(v_potential, 2) - 1 + 0.1d0))

    allocate(energy(0:lmax_potential))
    call energy_electrostatic_L_resolved_wrapper(energy, v_potential, Z_nuclear, RHO2NS, shgaunts, atomdata)

    ecou = sum(energy)

    deallocate(energy)
  endfunction ! energy

  !----------------------------------------------------------------------------
  !> Calculates electrostatic energy 'ecou' for arbitrary potential 'v_potential'
  !> and charge density 'rho2ns'.
  !>
  !> 1/2 \int \rho V
  subroutine energy_electrostatic_L_resolved_wrapper(energy, v_potential, Z_nuclear, RHO2NS, shgaunts, atomdata)
    use BasisAtom_mod, only: BasisAtom
    use ShapeGauntCoefficients_mod, only: ShapeGauntCoefficients

    double precision, intent(out) :: energy(0:)
    double precision, intent(in) :: v_potential(:,:,:)
    double precision, intent(in) :: Z_nuclear
    double precision, intent(in) :: RHO2NS(:,:,:)
    type(BasisAtom), intent(in) :: atomdata
    type(ShapeGauntCoefficients), intent(in) :: shgaunts

    integer :: nspind, lpot, lmax_potential

    nspind = atomdata%nspin

#define mesh atomdata%mesh_ptr
#define cell atomdata%cell_ptr

    CHECKASSERT( associated(atomdata%mesh_ptr) )
    CHECKASSERT( associated(atomdata%cell_ptr) )

    lpot = int(sqrt(size(rho2ns, 2) - 1 + 0.1d0))
    lmax_potential = int(sqrt(size(v_potential, 2) - 1 + 0.1d0))

    CHECKASSERT( (lpot+1)**2 == size(rho2ns, 2) )
    CHECKASSERT( (lmax_potential + 1)**2 == size(v_potential, 2) )
    CHECKASSERT( size(energy) == lmax_potential + 1)

    call energy_electrostatic_L_resolved(energy, lpot, lmax_potential, &
                  nspind, rho2ns, v_potential, Z_nuclear, mesh%r, mesh%drdi, &
                  mesh%ircut, mesh%ipan, shgaunts%imaxsh, cell%ifunm, shgaunts%ilm, &
                  shgaunts%gsh, cell%theta, cell%lmsp, &
                  mesh%irmd, cell%irid, cell%nfund, mesh%ipand, shgaunts%ngshd)
#undef cell
#undef mesh
  endsubroutine ! energy

  !------------------------------------------------------------------------------
  !>    Modified ecoub to calculate electrostatic energy.
  !>
  !>     Calculates 1/2 * \int \rho * V
  !>     for arbitrary potentials E.R.
  !>
  !>     potential 'vons' is given up to lmax_potential
  !>     charge density up to lpot
  !>
  !>     the effect of a Coulomb singularity on the interstitial potential can be taken into account
  !>     by setting Z_nuclear to a value not equal to 0.0 - this is needed for the Coulomb energy
  !>     Note: effect only calculated in interstital!
  double precision function energy_electrostatic(lpot, lmax_potential, &
                  nspin,rho2ns,vons, Z_nuclear, r, drdi, &
                  ircut,ipan,imaxsh,ifunm,ilm, &
                  gsh,thetas,lmsp, &
                  irmd, irid, nfund, ipand, ngshd) result(ecou)

    ! Arguments
    integer, intent(in) :: lpot
    integer, intent(in) :: lmax_potential
    integer :: nspin
    double precision :: rho2ns(irmd,(lpot+1)**2,2)
    double precision :: vons(irmd,(lmax_potential+1)**2,2)
    double precision, intent(in) :: Z_nuclear
    double precision :: r(irmd)
    double precision :: drdi(irmd)
    integer :: ircut(0:ipand)
    integer :: ipan
    integer :: imaxsh(0:(lmax_potential+1)**2)
    integer :: ifunm(*)
    integer :: ilm(ngshd,3)
    double precision :: gsh(ngshd)
    double precision :: thetas(irid,nfund)
    integer :: lmsp(*)
    integer :: irmd
    integer :: irid
    integer :: nfund
    integer :: ipand
    integer :: ngshd

    ! local
    double precision :: energy(0:lmax_potential)

    call energy_electrostatic_L_resolved(energy, lpot, lmax_potential, &
                  nspin,rho2ns,vons, Z_nuclear, r, drdi, &
                  ircut,ipan,imaxsh,ifunm,ilm, &
                  gsh,thetas,lmsp, &
                  irmd, irid, nfund, ipand, ngshd)

    ecou = sum(energy)
  endfunction ! energy

  !------------------------------------------------------------------------------
  !>    Modified ecoub to calculate electrostatic energies - returns l-resolved
  !>    energy in 'energy'
  !>
  !>     Calculates 1/2 * \int \rho * V
  !>     for arbitrary potentials E.R.
  !>
  !>     potential 'vons' is given up to lmax_potential
  !>     charge density up to lpot
  !>
  !>     the effect of a Coulomb singularity on the interstitial potential can be taken into account
  !>     by setting Z_nuclear to a value not equal to 0.0 - this is needed for the Coulomb energy
  !>     Note: effect only calculated in interstital!
  subroutine energy_electrostatic_L_resolved(energy, lpot, lmax_potential, &
                  nspin,rho2ns,vons, Z_nuclear, r, drdi, &
                  ircut,ipan,imaxsh,ifunm,ilm, &
                  gsh,thetas,lmsp, &
                  irmd, irid, nfund, ipand, ngshd)
                  
  ! ************************************************************************
  !
  !
  !     calculate the electrostatic potential-energies
  !
  !                          rc
  !      ecou(i) =  1/2 (  {  s dr' vons(r',lm,i)*rho2ns(r',lm,i,1) }
  !                           0
  !
  !
  !
  !                                         ( {..} = summed over lm )
  !             (see notes by b.drittler)
  !     vons is the coulomb potential of the atom without the nuclear
  !             potential of the atom
  !     rho2ns(...,1) is the real charge density times r**2
  !
  !      both developed into spherical harmonics . (see deck rholm)
  !
  !
  !                 modified for band structure code
  !                               b.drittler   jan. 1990
  !-----------------------------------------------------------------------
    use Quadrature_mod, only: simpson

    double precision, intent(out) :: energy(0:lmax_potential)

    integer, intent(in) :: lpot
    integer, intent(in) :: lmax_potential
    integer :: nspin
    double precision :: rho2ns(irmd,(lpot+1)**2,2)
    double precision :: vons(irmd,(lmax_potential+1)**2,2)
    double precision, intent(in) :: Z_nuclear
    double precision :: r(irmd)
    double precision :: drdi(irmd)
    integer :: ircut(0:ipand)
    integer :: ipan
    integer :: imaxsh(0:(lmax_potential+1)**2)
    integer :: ifunm(*)
    integer :: ilm(ngshd,3)
    double precision :: gsh(ngshd)
    double precision :: thetas(irid,nfund)
    integer :: lmsp(*)
    integer :: irmd
    integer :: irid
    integer :: nfund
    integer :: ipand
    integer :: ngshd

    ! Local variables
    double precision :: er(irmd), rfpi, rhosp, sgn, factor
    integer :: ifun, ipan1, ipot, ir, irc1, irh, irs1, ispin, j, l, m, lm, lm2

    rfpi = 2.d0*sqrt(pi)

    ipan1 = ipan
    irs1 = ircut(1)
    irc1 = ircut(ipan1)

    ! determine the right potential numbers - the Coulomb potential is not spin dependend
    ipot = nspin

    energy = 0.d0

    do l = 0, lmax_potential
      er(:) = 0.d0
      do ispin = 1, nspin

        sgn = -1.d0; if (ispin == nspin) sgn = 1.d0

        do m = -l, l
          lm = l*l + l + m + 1

          if (l <= lpot) then
            do ir = 1, irs1
              rhosp = (rho2ns(ir,lm,1) + sgn*rho2ns(ir,lm,nspin))*0.25d0
              er(ir) = er(ir) + rhosp*vons(ir,lm,ipot)
            enddo ! ir
          endif

          ! convolute with shape function
          do j = imaxsh(lm-1) + 1,imaxsh(lm)

            if (lmsp(ilm(j,3)) > 0) then
              ifun = ifunm(ilm(j,3))
              lm2 = ilm(j,2)

              factor = 0.d0; if (lm == 1) factor = 1.d0

              do ir = irs1+1, irc1
                irh = ir - irs1
                rhosp = (rho2ns(ir,lm2,1) + sgn*rho2ns(ir,lm2,nspin))*0.5d0
                er(ir) = er(ir) + rhosp*gsh(j) * thetas(irh,ifun)*(vons(ir,lm,ipot)*0.5d0 - factor*RFPI*Z_nuclear/r(ir))
              enddo ! ir
            endif ! lmsp(ilm(j,3)) > 0
          enddo ! j
          
        enddo ! m
        
      enddo ! ispin

      energy(l) = simpson(er, ipan1, ircut, drdi) ! now integrate

    enddo ! l = 0, lmax_potential
  endsubroutine ! energy

  !------------------------------------------------------------------------------
  !> Calculates generalised Madelung energy except a contribution that cancels with
  !> electrostatic and double-counting terms
  !>
  !> Since the nucleus is point shaped, only the l=0 components of the potential and
  !> density are required.
  !>
  !> An electrostatic Dirichlet problem is solved by specifying the potential on a
  !> surface of a sphere with a radius determined by 'ind_reference' <= ind_muffin_tin
  double precision function madelung_energy(vons_spherical, rho2ns_spherical, &
                          r, drdi, irmd, Z_nuclear, ind_reference) result(e_madelung)
    use Quadrature_mod, only: simpson
      double precision, intent(in) :: vons_spherical(irmd)
      double precision, intent(in) :: rho2ns_spherical(irmd)
      double precision, intent(in) :: r(irmd)
      double precision, intent(in) :: drdi(irmd)
      integer, intent(in) :: irmd
      double precision, intent(in) :: Z_nuclear
      integer, intent(in) :: ind_reference

      double precision :: er(irmd), rfpi, vmad, charge_in_sphere

      rfpi = 2.d0*sqrt(pi)

      er = 0.d0
      er(1:ind_reference) = rho2ns_spherical(1:ind_reference)*rfpi

      charge_in_sphere = simpson(er, 1, ind_reference, drdi)

      vmad = vons_spherical(ind_reference)/rfpi - 2.d0*charge_in_sphere/r(ind_reference)

      e_madelung = -0.5d0*z_nuclear*vmad
  endfunction ! madelung

  !------------------------------------------------------------------------------
  !> Correction of the madelung energy that occurs when the reference radius is not equal to
  !> the muffin-tin radius and this term does not cancel anymore.
  double precision function madelung_ref_radius_correction(rho2ns_spherical, &
                          r, drdi, irmd, Z_nuclear, ind_reference, ind_muffin_tin) result(e_correction)
    use Quadrature_mod, only: simpson
      double precision, intent(in)  :: rho2ns_spherical(irmd)
      double precision, intent(in)  :: r(irmd)
      double precision, intent(in)  :: drdi(irmd)
      integer, intent(in) :: irmd
      double precision, intent(in)  :: Z_nuclear
      integer, intent(in) :: ind_reference
      integer, intent(in) :: ind_muffin_tin

      double precision :: rfpi, vm, er(irmd)
      integer :: ii

      rfpi = 2.d0*sqrt(pi)

      er = 0.d0
      vm = 0.d0

      if (ind_reference /= ind_muffin_tin) then
        do ii = ind_reference, ind_muffin_tin
          er(ii) = rho2ns_spherical(ii)/r(ii)
        enddo ! ii
        vm = simpson(er, ind_reference, ind_muffin_tin, drdi)
      endif

      vm = -2.d0*rfpi*vm
      e_correction = -0.5d0*z_nuclear*vm
  endfunction ! madelung

endmodule ! total_energy_mod
