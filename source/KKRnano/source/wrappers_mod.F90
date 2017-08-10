
!------------------------------------------------------------------------------
!> This module provides wrappers for routines with enormous argument lists.
!>
!> Most of them are called in main2.
!> @todo check intents
module wrappers_mod
#include "macros.h"
  implicit none
  private
  
  public :: RHOVAL_wrapper, calcTmat_wrapper, calcdTmat_wrapper, MIXSTR_wrapper, VINTRAS_wrapper
  public :: RHOTOTB_wrapper, ESPCB_wrapper, EPOTINB_wrapper, VXCDRV_wrapper, MTZERO_wrapper
  public :: CONVOL_wrapper, RHOMOM_NEW_wrapper

  contains

#define mesh atomdata%mesh_ptr
#define chebmesh atomdata%chebmesh_ptr
#define cell atomdata%cell_ptr
  
  !----------------------------------------------------------------------------
  !> Calculate valence electron density.
  subroutine RHOVAL_wrapper(atomdata, ldorhoef, icst, nsra, rho2ns, r2nef, den, &
                            espv, gmatn, gaunts, emesh, ldau_data, method, &
                            korbit, theta_noco, phi_noco, theta_noco_old, &                          
                            phi_noco_old, angle_fixed, &
                            moment_x, moment_y, moment_z, &
                            muorb, iemxd, params) ! NOCO/SOC
    use BasisAtom_mod, only: BasisAtom
    use GauntCoefficients_mod, only: GauntCoefficients
    use EnergyMesh_mod, only: EnergyMesh
    use LDAUData_mod, only: LDAUData

    use ValenceDensity_mod, only: rhoval
    
    use NonCollinearMagnetism_mod, only: rhovalnew ! NOCO  
    use InputParams_mod, only: InputParams ! NOCO  
    
    logical, intent(in) :: ldorhoef
    integer, intent(in) :: icst !< num. born iterations
    integer, intent(in) :: nsra !< flag scalar relativistic
    double precision, intent(inout) :: rho2ns(:,:,:) ! inout?
    double precision, intent(inout) :: r2nef(:,:,:) ! inout?
    double complex, intent(inout) :: den(:,:,:)
    double precision, intent(inout) :: espv(:,:)
    double complex, intent(in) :: gmatn(:,:,:,:) !in
    type(BasisAtom), intent(inout) :: atomdata !in or inout?
    type(GauntCoefficients), intent(in) :: gaunts
    type(EnergyMesh), intent(inout) :: emesh !inout or in?
    type(LDAUData), intent(inout) :: ldau_data ! inout?
    integer, intent(in) :: method !< method for solving the single site problem, Volterra or Fredholm

    integer, intent(in)             :: korbit          ! NOCO
    double precision, intent(out)   :: theta_noco      ! NOCO
    double precision, intent(out)   :: phi_noco        ! NOCO
    integer (kind=1), intent(in)    :: angle_fixed     ! NOCO
    double precision, intent(out)   :: theta_noco_old  ! NOCO
    double precision, intent(out)   :: phi_noco_old    ! NOCO
    double precision, intent(out)   :: moment_x        ! NOCO
    double precision, intent(out)   :: moment_y        ! NOCO
    double precision, intent(out)   :: moment_z        ! NOCO
!    logical, intent(in)             :: soc            ! NOCO
!    double precision, intent(in)    :: socscale       ! NOCO
    double precision, intent(out)   :: muorb(0:,:)     ! NOCO
    integer, intent(in)             :: iemxd           ! NOCO
    type(InputParams), intent(in)   :: params          ! NOCO
    
    integer :: ispin, nspind, irmind, irnsd, lmaxd, l

    nspind = atomdata%nspin
    irmind = atomdata%potential%irmind
    irnsd = atomdata%potential%irmd - atomdata%potential%irmind
    lmaxd = atomdata%potential%lpot/2

    CHECKASSERT( 2*lmaxd == atomdata%potential%lpot )

    !-------------------------------------- NOCO ---------------------
    if (korbit == 1) then

        ! store angles from iteration before to compare
        theta_noco_old = theta_noco
        phi_noco_old   = phi_noco
      
        call rhovalnew(ldorhoef,emesh%ielast,nsra,nspind,lmaxd,emesh%ez,emesh%wez,atomdata%z_nuclear,  &
                      params%socscale,gaunts%cleb(:,1),gaunts%icleb,gaunts%iend, &
                      cell%ifunm,cell%lmsp,params%ncheb,  &
                      chebmesh%npan_tot,params%npan_log,params%npan_eq,mesh%r,mesh%irws,  &
                      chebmesh%rpan_intervall,chebmesh%ipan_intervall,  &
                      chebmesh%rnew,atomdata%potential%vinscheb,chebmesh%thetasnew, &
                      theta_noco,phi_noco,angle_fixed,moment_x,moment_y,moment_z,&
                      1,  &  ! ipot=1
                      den,espv,rho2ns,r2nef, gmatn(:,:,:,1), muorb,  & ! just one spin component of gmatn needed
                      atomdata%potential%lpot,lmaxd,mesh%irmd,chebmesh%irmd_new,iemxd, params%soc)
 
       ! calculate correct orbital moment
       do ispin=1,nspind
         do l=0,lmaxd+1
           muorb(l,3)=muorb(l,3)+muorb(l,ispin)
         enddo
       enddo
       do ispin=1,3
         do l=0,lmaxd+1
           muorb(lmaxd+2,ispin)=muorb(lmaxd+2,ispin)+muorb(l,ispin)
         enddo
       enddo
  !------------------------------------------- ---------------------

   else !KORBIT
    do ispin = 1, nspind

      ! has to be done after lloyd
      ! output: rho2ns, r2nef, den, espv
      call rhoval(ldorhoef, icst, emesh%ielast, &
                  nsra, ispin, nspind, emesh%ez, emesh%wezrn(1,ispin), &   ! unfortunately spin-dependent
                  mesh%drdi, mesh%r, mesh%irmin, &
                  atomdata%potential%vins(irmind,1,ispin), atomdata%potential%visp(1,ispin), &
                  atomdata%z_nuclear, mesh%ipan, mesh%ircut, &
                  cell%theta, cell%ifunm, cell%lmsp, &
                  rho2ns, r2nef, &
                  !den(0,1,ispin), espv(0,ispin), &
                  den(:,:,ispin), espv(:,ispin), &
                  gaunts%cleb, gaunts%loflm, gaunts%icleb, gaunts%iend, gaunts%jend, &
                  gmatn, &
                  ldau_data%ldau, ldau_data%nldau, ldau_data%lldau, ldau_data%phildau, ldau_data%wmldau, &
                  ldau_data%dmatldau, &
                  emesh%ielast, &
                  lmaxd, mesh%irmd, irnsd, cell%irid, mesh%ipand, cell%nfund, gaunts%ncleb, method)

    enddo ! ispin
   endif !KORBIT

  endsubroutine ! rhoval

  !------------------------------------------------------------------------------
  subroutine calcTmat_wrapper(atomdata, emesh, ie, ispin, icst, nsra, gaunts, tmatn, tr_alph, ldau_data, method)
    use BasisAtom_mod, only: BasisAtom
    use EnergyMesh_mod, only: EnergyMesh
    use LDAUData_mod, only: LDAUData
    use GauntCoefficients_mod, only: GauntCoefficients
    
    use SingleSite_mod, only: calcTmat

    type(BasisAtom), intent(in) :: atomdata
    type(GauntCoefficients), intent(in) :: gaunts
    type(EnergyMesh), intent(in) :: emesh
    type(LDAUData), intent(inout) :: ldau_data
    integer, intent(in) :: ie
    integer, intent(in) :: ispin
    integer, intent(in) :: icst
    integer, intent(in) :: nsra
    double complex, intent(inout) :: tmatn(:,:,:)
    double complex, intent(inout) :: tr_alph(:)
    integer, intent(in) :: method

    integer :: nspind, irmind, irnsd, lmaxd, lmmaxd

    nspind = atomdata%nspin
    irmind = atomdata%potential%irmind
    irnsd = atomdata%potential%irmd - atomdata%potential%irmind

    CHECKASSERT( atomdata%potential%irmd == mesh%irmd )

    lmaxd = atomdata%potential%lpot / 2
    lmmaxd = (lmaxd + 1)**2

    CHECKASSERT( lmaxd*2 == atomdata%potential%lpot )
    CHECKASSERT( size(tr_alph) == nspind )
    CHECKASSERT( size(tmatn, 3) == nspind )
    if( size(tmatn, 1) /= lmmaxd ) warn(6,"t-matrix does not have the same number of lm-components as given in the potential")
    if( size(tmatn, 2) /= lmmaxd ) warn(6,"t-matrix does not have the same number of lm-components as given in the potential"

    call calcTmat(ldau_data%ldau, ldau_data%nldau, icst, &
                  nsra, emesh%ez(ie), &
                  mesh%drdi, mesh%r, atomdata%potential%vins(irmind,1,ispin), &
                  atomdata%potential%visp(:,ispin), atomdata%z_nuclear, mesh%ipan, &
                  mesh%ircut, gaunts%cleb, gaunts%loflm, gaunts%icleb, gaunts%iend, &
                  tmatn(:,:,ispin), tr_alph(ispin), lmaxd, &
                  ldau_data%lldau, ldau_data%wmldau(:,:,:,ispin), &
                  gaunts%ncleb, mesh%ipand, mesh%irmd, irnsd, method)
  endsubroutine ! calcTmat

  !------------------------------------------------------------------------------
  subroutine calcdTmat_wrapper(atomdata, emesh, ie, ispin, icst, nsra, gaunts, dtde, tr_alph, ldau_data, method)
    use BasisAtom_mod, only: BasisAtom
    use RadialMeshData_mod, only: RadialMeshData
    use EnergyMesh_mod, only: EnergyMesh
    use LDAUData_mod, only: LDAUData
    use GauntCoefficients_mod, only: GauntCoefficients
    use SingleSite_mod, only: calcdTmat, calcdTmat_DeltaEz
    
    type(BasisAtom), intent(in) :: atomdata
    type(GauntCoefficients), intent(in) :: gaunts
    type(EnergyMesh), intent(in) :: emesh
    type(LDAUData), intent(inout) :: ldau_data
    integer, intent(in) :: ie
    integer, intent(in) :: ispin
    integer, intent(in) :: icst
    integer, intent(in) :: nsra
    double complex, intent(inout) :: dtde(:,:,:)
    double complex, intent(inout) :: tr_alph(:)
    integer, intent(in) :: method

    double complex :: delta_e_z
    integer :: nspind, irmind, irnsd, lmaxd, lmmaxd

    nspind = atomdata%nspin
    irmind = atomdata%potential%irmind
    irnsd = atomdata%potential%irmd - atomdata%potential%irmind

    CHECKASSERT( atomdata%potential%irmd == mesh%irmd )

    lmaxd = atomdata%potential%lpot / 2
    lmmaxd = (lmaxd + 1)**2

    CHECKASSERT( lmaxd*2 == atomdata%potential%lpot )
    CHECKASSERT( size(tr_alph) == nspind )
    CHECKASSERT( size(dtde, 3) == nspind )
    CHECKASSERT( size(dtde, 1) == lmmaxd )
    CHECKASSERT( size(dtde, 2) == lmmaxd )

    call calcdTmat_DeltaEz(delta_e_z, ie, emesh%npnt123(1), emesh%npnt123(2), emesh%npnt123(3), emesh%tk)

    call calcdTmat(ldau_data%ldau, ldau_data%nldau, icst, &
                  nsra, emesh%ez(ie), delta_e_z, &
                  mesh%drdi, mesh%r, atomdata%potential%vins(irmind,1,ispin), &
                  atomdata%potential%visp(:,ispin), atomdata%z_nuclear, mesh%ipan, &
                  mesh%ircut, gaunts%cleb, gaunts%loflm, gaunts%icleb, gaunts%iend, &
                  dtde(:,:,ispin), tr_alph(ispin), lmaxd, &
                  ldau_data%lldau, ldau_data%wmldau(:,:,:,ispin), &
                  gaunts%ncleb, mesh%ipand, mesh%irmd, irnsd, method)
  endsubroutine ! calcdTmat

  !------------------------------------------------------------------------------
  !> Wraps mixstr_new.
  subroutine MIXSTR_wrapper(atomdata, rmsavq, rmsavm, mixing, fcm)
    use BasisAtom_mod, only: BasisAtom

    double precision, intent(inout) :: rmsavq, rmsavm
    double precision, intent(in)    :: mixing, fcm
    type(BasisAtom), intent(inout) :: atomdata

    integer :: nspind, irnsd

    nspind = atomdata%nspin
    irnsd = atomdata%potential%irmd - atomdata%potential%irmind

    CHECKASSERT( atomdata%potential%irmd == mesh%irmd )

    call mixstr_new(rmsavq, rmsavm, atomdata%potential%lmpot, nspind, mixing, fcm, &
                    mesh%irc, mesh%irmin, mesh%r, mesh%drdi, atomdata%potential%vons, atomdata%potential%visp, atomdata%potential%vins, &
                    mesh%irmd, irnsd)

  endsubroutine ! mixstr

  !----------------------------------------------------------------------------
  !> Adds intracell potential. in and output: atomdata (changed)
  subroutine VINTRAS_wrapper(rho2ns, shgaunts, atomdata)
    use BasisAtom_mod, only: BasisAtom
    use ShapeGauntCoefficients_mod, only: ShapeGauntCoefficients

    double precision, intent(inout) :: rho2ns(:,:) ! inout?
    type(BasisAtom), intent(inout) :: atomdata
    type(ShapeGauntCoefficients), intent(in) :: shgaunts

    integer :: nspind

    nspind = atomdata%nspin


    !output: VONS
    call vintras_new(atomdata%potential%lpot, nspind, rho2ns, atomdata%potential%vons, &
    mesh%r, mesh%drdi, mesh%ircut, mesh%ipan, shgaunts%ilm, &
    cell%ifunm, shgaunts%imaxsh, shgaunts%gsh, &
    cell%theta, cell%lmsp, &
    mesh%irmd, cell%irid, cell%nfund, shgaunts%ngshd, mesh%ipand)

  endsubroutine ! vintras

  !----------------------------------------------------------------------------
  !> Calculates total charge and spin of cell.
  !>
  !> Core charge density must have been already treated.
  subroutine RHOTOTB_wrapper(catom, rho2ns, atomdata)
    use BasisAtom_mod, only: BasisAtom

    double precision, intent(inout) :: catom(:)
    double precision, intent(inout) :: rho2ns(:,:,:)
    type(BasisAtom), intent(inout) :: atomdata

    integer :: nspind

    nspind = atomdata%nspin

    CHECKASSERT( size(catom) == nspind )

    call rhototb_new(nspind, rho2ns, atomdata%core%rhocat, &
                    mesh%drdi, mesh%ircut, &
                    atomdata%potential%lpot, cell%nfu, cell%llmsp, cell%theta, mesh%ipan, &
                    catom, &
                    mesh%irmd, cell%irid, mesh%ipand, cell%nfund)

  endsubroutine ! rhototb


  !------------------------------------------------------------------------------
  !> Contribution of core electrons to total energy.
  !> @param[out] ESPC  1st component l-resolved energies s, p, d, f
  subroutine ESPCB_wrapper(espc, lcoremax, atomdata)
    use BasisAtom_mod, only: BasisAtom

    double precision, intent(out) :: espc(:,:) ! espc(0:3,nspin)
    integer, intent(out) :: lcoremax
    type(BasisAtom), intent(in) :: atomdata

    integer :: nspind

    nspind = atomdata%nspin

    call espcb_new(espc, nspind, atomdata%core%ecore, &
        atomdata%core%lcore(:,1:nspind), lcoremax, atomdata%core%ncore(1:nspind))

  endsubroutine ! espcb

  !------------------------------------------------------------------------------
  subroutine EPOTINB_wrapper(epotin, rho2ns, atomdata)
    use BasisAtom_mod, only: BasisAtom

    double precision, intent(out)   :: epotin
    double precision, intent(inout) :: rho2ns(:,:,:)
    type(BasisAtom), intent(in) :: atomdata

    integer :: nspind, irnsd

    nspind = atomdata%nspin

    irnsd = atomdata%potential%irmd - atomdata%potential%irmind

    CHECKASSERT( atomdata%potential%irmd == mesh%irmd )

    call epotinb_new(epotin, nspind, rho2ns, atomdata%potential%visp, mesh%r, mesh%drdi, &
                    mesh%irmin, mesh%irws, atomdata%potential%lpot, atomdata%potential%vins, &
                    mesh%ircut, mesh%ipan, atomdata%z_nuclear, &
                    mesh%irmd, irnsd, mesh%ipand)

  endsubroutine ! epotinb

  !------------------------------------------------------------------------------
  !> Add exchange-correlation to 'vons_potential' and calculate exchange correlation energy.
  subroutine VXCDRV_wrapper(vons_potential, exc, kxc, rho2ns, shgaunts, atomdata)
    use BasisAtom_mod, only: BasisAtom
    use ShapeGauntCoefficients_mod, only: ShapeGauntCoefficients

    double precision, intent(inout) :: vons_potential(:,:,:)
    double precision, intent(inout) :: exc(:)
    integer, intent(in)             :: kxc
    double precision, intent(in) :: rho2ns(:,:,:)
    type(BasisAtom), intent(in) :: atomdata
    type(ShapeGauntCoefficients), intent(in) :: shgaunts

    integer :: nspind, lpot
    integer, parameter :: KTE = 1 ! always calculate exchange energy

    nspind = atomdata%nspin
    lpot   = atomdata%potential%lpot

    CHECKASSERT( size(exc) == lpot+1 )
    CHECKASSERT( size(vons_potential, 2) == (lpot + 1)**2 )

    call vxcdrv_new(exc, KTE, kxc, lpot, nspind, rho2ns, &
              vons_potential, mesh%r, mesh%drdi, mesh%a, &
              mesh%irws, mesh%ircut, mesh%ipan, shgaunts%gsh, shgaunts%ilm, shgaunts%imaxsh, cell%ifunm, &
              cell%theta, cell%lmsp, &
              mesh%irmd, cell%irid, cell%nfund, shgaunts%ngshd, mesh%ipand)

  endsubroutine ! vxcdrv

  !----------------------------------------------------------------------------
  subroutine MTZERO_wrapper(vav0, vol0, atomdata)
    use BasisAtom_mod, only: BasisAtom

    type(BasisAtom), intent(inout) :: atomdata
    double precision, intent(inout) :: vav0, vol0

    integer :: nspind

    nspind = atomdata%nspin

    vav0 = 0.d0
    vol0 = 0.d0

    !output: vav0, vol0
    call mtzero_new(atomdata%potential%lmpot, nspind, atomdata%potential%vons, &
                    atomdata%z_nuclear, mesh%r, mesh%drdi, mesh%imt, mesh%ircut, &
                    mesh%ipan, cell%lmsp, cell%ifunm, &
                    cell%theta, mesh%irws, vav0, vol0, &
                    mesh%irmd, cell%irid, cell%nfund, mesh%ipand)
  endsubroutine ! MTzero


  !----------------------------------------------------------------------------
  subroutine CONVOL_wrapper(vbc, shgaunts, atomdata)
    use BasisAtom_mod, only: BasisAtom
    use ShapeGauntCoefficients_mod, only: ShapeGauntCoefficients

    type(BasisAtom), intent(inout) :: atomdata
    type(ShapeGauntCoefficients), intent(in) :: shgaunts
    double precision, intent(in)    :: vbc(2)

    integer :: ispin, nspind

    nspind = atomdata%nspin

    do ispin = 1, nspind

      call shiftpotential(atomdata%potential%vons(:,:,ispin), mesh%ircut(mesh%ipan), vbc(ispin))

      !output: vons (changed)
      call convol_new(mesh%ircut(1), mesh%irc, &
                      shgaunts%imaxsh(shgaunts%lmpotd), shgaunts%ilm, cell%ifunm, shgaunts%lmpotd, shgaunts%gsh, &
                      cell%theta, atomdata%z_nuclear, &
                      mesh%r, atomdata%potential%vons(1,1,ispin), cell%lmsp, &
                      cell%irid, cell%nfund, mesh%irmd, shgaunts%ngshd)
    enddo ! ispin
  endsubroutine ! convol


!   !------------------------------------------------------------------------------
!   !> Do core relaxation for all spin directions.
!   !>
!   !> @param[in]      E1       bottom energy of valence energy contour
!   !> @param[in]     NSRA      flag for scalar relativistic calculation
!   !> @param[in, out] atomdata  basis atom - changed on output
!   subroutine RHOCORE_wrapper(e1, nsra, atomdata)
!     use BasisAtom_mod, only: BasisAtom
!     use AtomicCore_mod, only: rhocore
! 
!     double precision, intent(in)    :: e1
!     integer, intent(in)             :: nsra
!     type(BasisAtom), intent(inout)  :: atomdata
! 
!     atomdata%core%qc_corecharge = rhocore(e1, nsra, atomdata%nspin, atomdata%atom_index, &  ! atom_index is used only for debugging output
!                   mesh%drdi, mesh%r, atomdata%potential%visp(:,:), &
!                   mesh%a, mesh%b, atomdata%z_nuclear, &
!                   mesh%ircut, atomdata%core%rhocat, &
!                   atomdata%core%ecore(:,:), atomdata%core%ncore(:), atomdata%core%lcore(:,:), &
!                   mesh%irmd)
!   endsubroutine

#undef  cell
#undef  mesh

  !-------------------------------------------------------------------------------
  !> A wrapper for the subroutine RHOMOM_NEW.
  subroutine RHOMOM_NEW_wrapper(cmom, cminst, rho2ns, cell, mesh, shgaunts)
    use ShapefunData_mod, only: ShapefunData
    use RadialMeshData_mod, only: RadialMeshData
    use ShapeGauntCoefficients_mod, only: ShapeGauntCoefficients

    type(ShapefunData), intent(in) :: cell
    type(RadialMeshData), intent(in) :: mesh
    type(ShapeGauntCoefficients), intent(in) :: shgaunts
    double precision, intent(out) :: cminst(:)
    double precision, intent(out) :: cmom(:)
    double precision, intent(in)  :: rho2ns(:,:)

    integer :: lpot

    lpot = shgaunts%lmax * 2

    CHECKASSERT( size(CMOM) == (LPOT+1)**2)
    CHECKASSERT( size(CMINST) == (LPOT+1)**2)
    CHECKASSERT( size(RHO2NS, 1) == mesh%irmd)
    CHECKASSERT( size(RHO2NS, 2) == (LPOT+1)**2)
    CHECKASSERT( cell%lmmax_shape == (2*LPOT + 1) ** 2)

    !call rhomom_new(cmom, cminst, lpot, rho2ns, r, drdi, ircut, ipan, ilm, ifunm, imaxsh, gsh, thetas, lmsp, irmd, irid, nfund, ipand, ngshd)

    !output: cmom, cminst
    call rhomom_new(cmom, cminst, lpot, rho2ns, &
              mesh%r, mesh%drdi, mesh%ircut, mesh%ipan, shgaunts%ilm, &
              cell%ifunm, shgaunts%imaxsh, shgaunts%gsh, &
              cell%theta, cell%lmsp, &
              mesh%irmd, cell%irid, cell%nfund, mesh%ipand, shgaunts%ngshd)

  endsubroutine ! rhomom

!==============================================================================
! Some helper routines for the wrappers
!==============================================================================

  !----------------------------------------------------------------------------
  !> Shift lm-decomposed potential by a constant sqrt(4 * pi) * VBC.
  subroutine shiftPotential(vons_ispin, index_rmax, vbc)
    double precision, intent(inout) :: vons_ispin(:,:)
    integer, intent(in) :: index_rmax
    double precision, intent(in) :: vbc

    integer :: ir
    double precision :: rfpi

    rfpi = sqrt(16.d0*atan(1.d0))

    do ir = 1, index_rmax
      ! a constant potential shift affects only lm=1 component
      vons_ispin(ir,1) = vons_ispin(ir,1) + rfpi*vbc
    enddo ! ir

  endsubroutine ! shift


endmodule ! wrappers_mod
