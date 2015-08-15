#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

!------------------------------------------------------------------------------
!> Module with routines for alternative calculation of the
!> total energy.
!
! Weinert: -(E_coulomb + E_Z_interstitial) + E_madelung - 2*E_VXC + E_single_particle + E_core + E_XC + E_LDAU...
!
! E_VXC = 1/2 \int V_XC \rho
module total_energy_mod
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

  double precision, intent(in), dimension(:,:,:) :: v_potential
  double precision, intent(in) :: Z_nuclear
  double precision, intent(in) :: RHO2NS(:,:,:)
  type (BasisAtom), intent(in) :: atomdata
  type (ShapeGauntCoefficients), intent(in) :: shgaunts

  !-------- locals
  integer :: lmax_potential
  double precision, allocatable :: energy(:)

  lmax_potential = int(sqrt(size(v_potential, 2) - 1 + 0.1d0))

  allocate(energy(0:lmax_potential))
  call energy_electrostatic_L_resolved_wrapper(energy, v_potential, Z_nuclear, RHO2NS, shgaunts, atomdata)

  ecou = sum(energy)

  deallocate(energy)

end function

!----------------------------------------------------------------------------
!> Calculates electrostatic energy 'ecou' for arbitrary potential 'v_potential'
!> and charge density 'rho2ns'.
!>
!> 1/2 \int \rho V
subroutine energy_electrostatic_L_resolved_wrapper(energy, v_potential, Z_nuclear, RHO2NS, shgaunts, atomdata)
  use BasisAtom_mod, only: BasisAtom
  use ShapeGauntCoefficients_mod, only: ShapeGauntCoefficients
  
  use RadialMeshData_mod, only: RadialMeshData
  use CellData_mod, only: CellData

  double precision, intent(out) :: energy(0:)
  double precision, intent(in), dimension(:,:,:) :: v_potential
  double precision, intent(in) :: Z_nuclear
  double precision, intent(in) :: RHO2NS(:,:,:)
  type (BasisAtom), intent(in) :: atomdata
  type (ShapeGauntCoefficients), intent(in) :: shgaunts

  !-------- locals
  integer :: nspind
  type (RadialMeshData), pointer :: mesh
  type (CellData), pointer       :: cell
  integer :: lpot, lmax_potential

  nspind = atomdata%nspin

  mesh => atomdata%mesh_ptr
  cell => atomdata%cell_ptr

  CHECKASSERT( associated(atomdata%mesh_ptr) )
  CHECKASSERT( associated(atomdata%cell_ptr) )

  lpot = int(sqrt(size(rho2ns, 2) - 1 + 0.1d0))
  lmax_potential = int(sqrt(size(v_potential, 2) - 1 + 0.1d0))

  CHECKASSERT( (lpot+1)**2 == size(rho2ns, 2) )
  CHECKASSERT( (lmax_potential + 1)**2 == size(v_potential, 2) )
  CHECKASSERT( size(energy) == lmax_potential + 1)

  call energy_electrostatic_L_resolved(energy, lpot, lmax_potential, &
                 nspind, rho2ns, v_potential, Z_nuclear, mesh%r, mesh%drdi, &
                 mesh%ircut, mesh%ipan, shgaunts%imaxsh, cell%shdata%ifunm, shgaunts%ilm, &
                 shgaunts%gsh, cell%shdata%theta, cell%shdata%lmsp, &
                 mesh%irmd, cell%shdata%irid, cell%shdata%nfund, mesh%ipand, shgaunts%ngshd)

end subroutine

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
  double precision, dimension(irmd,(lpot + 1)**2,2) :: rho2ns
  double precision, dimension(irmd,(lmax_potential + 1)**2,2) :: vons
  double precision, intent(in) :: Z_nuclear
  double precision, dimension(irmd) :: r
  double precision, dimension(irmd) :: drdi
  integer, dimension(0:ipand) :: ircut
  integer :: ipan
  integer, dimension(0:(lmax_potential + 1)**2) :: imaxsh
  integer, dimension(*) :: ifunm
  integer, dimension(ngshd,3) :: ilm
  double precision, dimension(ngshd) :: gsh
  double precision, dimension(irid,nfund) :: thetas
  integer, dimension(*) :: lmsp
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

end function

!------------------------------------------------------------------------------
!>    Modified ecoub to calculate electrostatic energies - returns L-resolved
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

  double precision, intent(out) :: energy(0:lmax_potential)

  integer, intent(in) :: lpot
  integer, intent(in) :: lmax_potential
  integer :: nspin
  double precision, dimension(irmd,(lpot + 1)**2,2) :: rho2ns
  double precision, dimension(irmd,(lmax_potential + 1)**2,2) :: vons
  double precision, intent(in) :: Z_nuclear
  double precision, dimension(irmd) :: r
  double precision, dimension(irmd) :: drdi
  integer, dimension(0:ipand) :: ircut
  integer :: ipan
  integer, dimension(0:(lmax_potential + 1)**2) :: imaxsh
  integer, dimension(*) :: ifunm
  integer, dimension(ngshd,3) :: ilm
  double precision, dimension(ngshd) :: gsh
  double precision, dimension(irid,nfund) :: thetas
  integer, dimension(*) :: lmsp
  integer :: irmd
  integer :: irid
  integer :: nfund
  integer :: ipand
  integer :: ngshd

  ! Local variables
  double precision :: rfpi
  double precision :: rhosp
  double precision :: sign
  integer :: i
  integer :: ifun
  integer :: ipan1
  integer :: ipot
  integer :: ir
  integer :: irc1
  integer :: irh
  integer :: irs1
  integer :: ispin
  integer :: j
  integer :: L
  integer :: lm
  integer :: lm2
  integer :: m
  double precision, dimension(irmd) :: er
  double precision :: factor

  external :: simpk

  intrinsic :: atan,sqrt
  !     ..
  rfpi = sqrt(16.d0*atan(1.0d0))
  !
  ipan1 = ipan
  irs1 = ircut(1)
  irc1 = ircut(ipan1)
  !
  !--->   determine the right potential numbers - the coulomb potential
  !       is not spin dependend
  !
  ipot = nspin

  energy = 0.0d0

  do 80 L = 0, lmax_potential

    er = 0.0d0

    do 70 ispin = 1,nspin

      if (ispin.eq.nspin) then
        sign = 1.0d0
      else
        sign = -1.0d0
      end if

      do 60 m = -L,L
        lm = L*L + L + m + 1


        if (L <= lpot) then
        do 20 i = 1,irs1
          rhosp = (rho2ns(i,lm,1)+ &
                  sign*rho2ns(i,lm,nspin))/4.0d0
          er(i) = er(i) + rhosp*vons(i,lm,ipot)
  20       continue
        endif
  !
  !--->           convolute with shape function
  !
        do 50 j = imaxsh(lm-1) + 1,imaxsh(lm)
          lm2 = ilm(j,2)

          if (lmsp(ilm(j,3)).gt.0) then
            ifun = ifunm(ilm(j,3))

              if (lm == 1) then
                factor = 1.0d0
              else
                factor = 0.0d0
              endif

              do 40 ir = irs1 + 1,irc1
                irh = ir - irs1
                rhosp = (rho2ns(ir,lm2,1)+ &
                        sign*rho2ns(ir,lm2,nspin))/2.0d0
                er(ir) = er(ir) + rhosp*gsh(j)* &
                         thetas(irh,ifun)* (vons(ir,lm,ipot) / 2.0d0 - factor * Z_nuclear / r(ir) * RFPI)
  40          continue

          end if

  50       continue

  60     continue

  70   continue
  !
  !--->     now integrate
  !
    call simpk(er, energy(L),ipan1,ircut, drdi)

  80   continue                    ! L = 0, lmax_potential

end subroutine

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

    double precision, intent(in)  :: vons_spherical(irmd)
    double precision, intent(in)  :: rho2ns_spherical(irmd)
    double precision, intent(in)  :: r(irmd)
    double precision, intent(in)  :: drdi(irmd)
    integer, intent(in) :: irmd
    double precision, intent(in)  :: Z_nuclear
    integer, intent(in) :: ind_reference

    double precision :: ER(irmd)
    double precision :: rfpi
    double precision :: VMAD, charge_in_sphere

    rfpi = sqrt(16.0d0 * atan(1.0d0))

    ER = 0.0d0
    ER(1:ind_reference) = RHO2NS_spherical(1:ind_reference) * RFPI

    charge_in_sphere = 0.0d0
    CALL SIMP3(ER,charge_in_sphere,1,ind_reference,DRDI)

    VMAD = VONS_spherical(ind_reference)/RFPI - 2.0D0 * charge_in_sphere/R(ind_reference)

    e_madelung = -Z_nuclear/2.0d0 * VMAD

end function

!------------------------------------------------------------------------------
!> Correction of the madelung energy that occurs when the reference radius is not equal to
!> the muffin-tin radius and this term does not cancel anymore.
double precision function madelung_ref_radius_correction(rho2ns_spherical, &
                         r, drdi, irmd, Z_nuclear, ind_reference, ind_muffin_tin) result(e_correction)

    double precision, intent(in)  :: rho2ns_spherical(irmd)
    double precision, intent(in)  :: r(irmd)
    double precision, intent(in)  :: drdi(irmd)
    integer, intent(in) :: irmd
    double precision, intent(in)  :: Z_nuclear
    integer, intent(in) :: ind_reference
    integer, intent(in) :: ind_muffin_tin

    double precision :: ER(irmd)
    double precision :: rfpi
    double precision :: VM
    integer :: ii

    rfpi = sqrt(16.0d0 * atan(1.0d0))

    ER = 0.0D0
    VM = 0.0d0

    if (ind_reference /= ind_muffin_tin) then
      do ii = ind_reference, ind_muffin_tin
        ER(ii) = RHO2NS_spherical(ii)/R(ii)
      enddo

      CALL SIMP3(ER,VM,ind_reference,ind_muffin_tin,DRDI)
    endif

    VM = -2.0D0*RFPI*VM

    e_correction = -Z_nuclear/2.0d0 * VM

end function

end module total_energy_mod
