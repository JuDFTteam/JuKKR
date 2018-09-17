#define OUTPUT_CHARGE_AND_POTENTIAL_MOMENTS

!> Adds the Madelung Potential to all atoms.
module MadelungPotential_mod
#include "macros.h"
  implicit none
  private
  public :: addMadelungPotentialnew_com

  contains

  !----------------------------------------------------------------------------
  !> Add Madelung potential to VONS.
  !> Needs SMAT (Lattice sums from calculateMadelungLatticeSum)
  !> principal input: CMOM, CMINST, SMAT, VONS --> VONS (changed)
  !> Wrapper for VMADELBLK
  subroutine addMadelungPotentialnew_com(calc_data, zat, rbasis, communicator)
    use CalculationData_mod, only: CalculationData
    
    type(CalculationData), intent(inout) :: calc_data
    double precision, intent(in) :: zat(:), rbasis(:,:)
    integer, intent(in) :: communicator

#define mc calc_data%madelung_calc
    call VMADELBLK_new2_com(calc_data, mc%LPOT, size(zat), zat, rbasis, mc%lmpotd, &
         mc%clebsch%CLEB, mc%clebsch%ICLEB, mc%clebsch%IEND, &
         mc%LMXSPD,mc%clebsch%NCLEBD, mc%clebsch%LOFLM, mc%DFAC, communicator)
#undef mc
  endsubroutine ! add

  ! **********************************************************************
  !
  !>    calculate the madelung potentials and add these to the poten-
  !>    tial v  (in the spin-polarized case for each spin-direction
  !>    this is the same).
  !
  !>    it uses the structure dependent matrices AVMAD and BVMAD which
  !>    are calculated once in the subroutine MADELUNG3D
  !>    ( may 2004)
  !>    the charge-moments are calculated in the subroutine vintras,
  !>    therefore vintras has to be called first.
  !>    the madelung-potential is expanded into spherical harmonics.
  !>    the lm-term of the potential v of the atom i is given by
  !>
  !>     v(r,lm,i) =  (-r)**ell * {avmad(i,i2,lm,ell'emm')*cmom(i2,ell'emm')
  !>                                              +bvmad(i,i2,lm)*z(i2)}
  !>
  !>    summed over i2 (all atoms) and ell'emm'
  !>    (see notes by b.drittler)
  !>
  !>                              b.drittler   nov. 1989
  !>
  !>    adopted for the case of more atoms on the same site, summation is
  !>    done over the occupants of that site, the charge is weighted with
  !>    the appropriate concentration of the occupant  v.popescu feb. 2002
  !>
  !>    impurity-program adopted feb. 2004 (according to n.papanikalou)
  !>
  ! **********************************************************************
  subroutine VMADELBLK_new2_com(calc_data, lpot, naez, zat, rbasis, &
    lmpot, cleb, icleb, iend, lmxspd, nclebd, loflm, dfac, communicator)

    use CalculationData_mod, only: CalculationData, getEnergies, getDensities, getEnergies, getAtomData
    use EnergyResults_mod, only: EnergyResults
    use DensityResults_mod, only: DensityResults
    use BasisAtom_mod, only: BasisAtom
    use RadialMeshData_mod, only: RadialMeshData
    use ChunkIndex_mod, only: getRankAndLocalIndex
    include 'mpif.h'

    type(CalculationData), intent(in) :: calc_data
    integer, intent(in) :: iend, lpot, lmxspd, nclebd, lmpot, naez
    double precision, intent(in) :: zat(:)
    double precision, intent(in) :: rbasis(:,:)
    double precision, intent(in) :: cleb(nclebd)
    double precision, intent(in) :: dfac(0:lpot,0:lpot)
    integer, intent(in) :: icleb(nclebd,3)
    integer, intent(in) :: loflm(lmxspd)
    integer, intent(in) :: communicator

    type(BasisAtom), pointer :: atomdata
    type(DensityResults), pointer :: densities
    type(EnergyResults), pointer :: energies
    type(RadialMeshData), pointer :: mesh
    integer :: i2, ilm, ierr, nranks, myrank
    double precision :: cmom_save((lpot+1)**2), cminst_save((lpot+1)**2)
    integer :: num_local_atoms, ila, jla
    integer(kind=4) :: chunk_ind(2)
    integer :: rank, lmpotd
    
    lmpotd = (lpot+1)**2

    num_local_atoms = calc_data%num_local_atoms

    do ila = 1, num_local_atoms
      energies => getEnergies(calc_data, ila)
      energies%AC_madelung = 0.d0
    enddo
    
    call MPI_Comm_size(communicator, nranks, ierr)
    call MPI_Comm_rank(communicator, myrank, ierr)

    !===== begin loop over all atoms =========================================
    ! O(N**2) in calculation and communication !!!
#ifdef  OUTPUT_CHARGE_AND_POTENTIAL_MOMENTS
    if (0 == myrank) write(300, '(/,a,i0)') "# new iteration: charge moments"
#endif      

    do i2 = 1, naez

      ! use omp single for MPI part?
      ! if MPI comm. dominates over calculation, OpenMP won't do much

      chunk_ind(:) = getRankAndLocalIndex(i2, naez, nranks)
      rank = chunk_ind(1)
      jla  = chunk_ind(2)

      !-------- bcast information of cmom and cminst begin --------------------

      if (myrank == rank) then

        CHECKASSERT( i2 == calc_data%atom_ids(jla) )

        densities => getDensities(calc_data, jla)

        do ilm = 1, lmpot
          cmom_save(ilm)   = densities%CMOM(ilm)
          cminst_save(ilm) = densities%CMINST(ilm)
        enddo
      else
        cmom_save   = 0.d0
        cminst_save = 0.d0
      endif

      call MPI_Bcast(cmom_save,   lmpotd, MPI_DOUBLE_PRECISION, rank, communicator, ierr)
      call MPI_Bcast(cminst_save, lmpotd, MPI_DOUBLE_PRECISION, rank, communicator, ierr)

      !-------- bcast information of cmom and cminst end----------------------

#ifdef  OUTPUT_CHARGE_AND_POTENTIAL_MOMENTS
      if (0 == myrank) write(300, '(f8.3,3f16.9,999es24.16)') zat(i2), rbasis(1:3,i2), cmom_save(:) + cminst_save(:)
#endif      

      ! --> calculate avmad(lm1,lm2)

      do ila = 1, num_local_atoms  ! no OpenMP ?
        energies => getEnergies(calc_data, ila)

        call sumAC(energies%AC_madelung, cmom_save, cminst_save, zat(i2), &
                   calc_data%madelung_sum_a(ila)%smat(:,i2), &
                   lpot, cleb, icleb, iend, loflm, dfac)

      enddo ! ila

    enddo ! i2
    !===== endloop over all atoms =========================================

    !     contributions are accumulated in AC !!!

#ifdef  OUTPUT_CHARGE_AND_POTENTIAL_MOMENTS
    call MPI_Barrier(communicator, ierr)
    if (0 == myrank) write(400, '(/,a,i0)') "# new iteration: potential moments"
#endif      

    do ila = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ila)
      energies => getEnergies(calc_data, ila)
      mesh => atomdata%mesh_ptr
      
#ifdef  OUTPUT_CHARGE_AND_POTENTIAL_MOMENTS
      i2 = calc_data%atom_ids(ila)
      if (nranks > 1) then
        write(4000+myrank, '(i8.8,f8.3,3f16.9,999es24.16)') i2, zat(i2), rbasis(1:3,i2), energies%AC_madelung(:)
      else
        write(400, '(f8.3,3f16.9,999es24.16)') zat(i2), rbasis(1:3,i2), energies%AC_madelung(:)
      endif
#endif      

      call addPot(atomdata%potential%VONS, energies%VMAD, &
                  energies%AC_madelung, atomdata%potential%LPOT, &
                  mesh%R, mesh%IRCUT, mesh%IPAN, atomdata%potential%NSPIN)

    enddo ! ila

  endsubroutine ! VMADELBLK_new2_com

  !----------------------------------------------------------------------------
  subroutine addPot(vons, vmad, ac, lpot, r, ircut, ipan, nspin)
    use Constants_mod, only: pi
    double precision, intent(inout) :: vons(:,:,:)
    double precision, intent(out)   :: vmad
    double precision, intent(in)    :: ac(:)
    integer, intent(in) :: lpot
    double precision, intent(in) :: r(:)
    integer, intent(in) :: ircut(0:)
    integer, intent(in) :: ipan
    integer, intent(in) :: nspin

    integer :: irs1, ell, emm, lm, ir, ispin
    
    vmad = 0.d0

    irs1 = ircut(ipan)

    do ell = 0, lpot
      do emm = -ell, ell
        lm = ell*ell + ell + emm + 1

        if (lm == 1) vmad = ac(1)/sqrt(4.d0*pi)

        !---> add to v the intercell-potential

        do ispin = 1, nspin

          !---> in the case of ell=0 : r(1)**ell is not defined

          if (ell == 0) vons(1,1,ispin) = vons(1,1,ispin) + ac(lm)
          do ir = 2, irs1
            vons(ir,lm,ispin) = vons(ir,lm,ispin) + (-r(ir))**ell*ac(lm)
          enddo ! ir
        enddo ! ispin

      enddo ! emm
    enddo ! ell

  endsubroutine ! add


  !------------------------------------------------------------------------------
  subroutine sumAC(ac, cmom_save, cminst_save, zat_i2, smat_i2, lpot, cleb, icleb, iend, loflm, dfac)
  use Constants_mod, only: pi
    double precision, intent(inout) :: ac(:) ! multipole moments of the potential in the local atom i1
    double precision, intent(in) :: cmom_save(:) ! multipole moments inside the MT sphere of the remote atom i2
    double precision, intent(in) :: cminst_save(:) ! multipole moments in the interstitial of the remote atom i2
    double precision, intent(in) :: zat_i2 ! nuclear charge of the remote atom i2
    double precision, intent(in) :: smat_i2(:) ! Madelung matrix of the remote atom i2 
    double precision, intent(in) :: cleb(:) ! Clebsh-Gordon coefficients
    integer, intent(in) :: iend ! number of Clebsh-Gordon coefficients
    double precision, intent(in) :: dfac(0:,0:) ! precomputed factorials
! --> calculate:                                (2*(l+l')-1)!!
!                 dfac(l,l') = (4*pi)**2 *  ----------------------
!                                           (2*l+1)!! * (2*l'+1)!!
    integer, intent(in) :: icleb(:,:) ! Clebsh-Gordon indices
    integer, intent(in) :: loflm(:) ! retrieve ell(lm)
    integer, intent(in) :: lpot ! angular momentum truncation

    integer :: lm1, lm2, lm3, i, lmmax, lmpot, ell
    double precision :: avmad((lpot+1)**2,(lpot+1)**2)
    double precision :: bvmad((lpot+1)**2)

    lmpot = (lpot+1)**2
    lmmax = lmpot

    avmad(:,:) = 0.d0

    do i = 1, iend
      lm1 = icleb(i,1)
      lm2 = icleb(i,2)
      lm3 = icleb(i,3)

      ! --> this loop has to be calculated only for l1+l2=l3

      avmad(lm1,lm2) = avmad(lm1,lm2) + 2.d0*dfac(loflm(lm1),loflm(lm2))*cleb(i)*smat_i2(lm3)
      !             here, the prefactor 2.d0 is due to Rydberg units
    enddo ! i

    ! --> calculate bvmad(lm1)

    bvmad(:) = 0.d0

    do lm1 = 1, lmpot
      ell = loflm(lm1)
      ! bvmat takes care of the correct interaction of nuclei:
      ! a prefactor of 4pi/(2l+1) comes from the Poisson equation
      ! another prefactor of -2 is because of Rydberg units: V(r) = -2Z/r
      ! the interaction is diagonal in lm
      bvmad(lm1) = bvmad(lm1) - (8.d0*pi/(2*ell + 1.d0))*smat_i2(lm1)
    enddo ! lm1

    do lm1 = 1, lmpot
      ac(lm1) = ac(lm1) + bvmad(lm1)*zat_i2

      !---> take moments of sphere

      do lm2 = 1, lmmax
        ac(lm1) = ac(lm1) + avmad(lm1,lm2)*(cmom_save(lm2) + cminst_save(lm2)) ! matrix multiply with the total multipole moments
      enddo ! lm2
    enddo ! lm1

  endsubroutine ! sum

endmodule ! MadelungPotential_mod
