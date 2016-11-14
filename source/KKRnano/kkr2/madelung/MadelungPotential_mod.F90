
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
  subroutine addMadelungPotentialnew_com(calc_data, ZAT, rank, atoms_per_proc, communicator)
    use CalculationData_mod, only: CalculationData
    
    type(CalculationData), intent(inout) :: calc_data
    double precision, intent(in) :: ZAT(:)
    integer, intent(in) :: rank, atoms_per_proc
    integer, intent(in) :: communicator

    integer :: naez

    naez = size(ZAT)
#define mc calc_data%madelung_calc
    call VMADELBLK_new2_com(calc_data, mc%LPOT, naez, ZAT, mc%LMPOTD, &
         mc%clebsch%CLEB, mc%clebsch%ICLEB, mc%clebsch%IEND, &
         mc%LMXSPD,mc%clebsch%NCLEBD, mc%clebsch%LOFLM, mc%DFAC, &
         rank, atoms_per_proc, &
         communicator)
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
  !>     v(r,lm,i) =  (-r)**l * {avmad(i,i2,lm,l'm')*cmom(i2,l'm')
  !>                                              +bvmad(i,i2,lm)*z(i2)}
  !>
  !>    summed over i2 (all atoms) and l'm'
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
  subroutine VMADELBLK_new2_com(calc_data, lpot, naez, zat, &
    lmpot, cleb, icleb, iend, lmxspd, nclebd, loflm, dfac, mylrank, atoms_per_proc, communicator)

    use CalculationData_mod, only: CalculationData, getEnergies, getDensities, getEnergies, getAtomData
    use EnergyResults_mod, only: EnergyResults
    use DensityResults_mod, only: DensityResults
    use BasisAtom_mod, only: BasisAtom
    use RadialMeshData_mod, only: RadialMeshData

    include 'mpif.h'

    type(CalculationData), intent(inout) :: calc_data

    integer, intent(in) :: iend, lpot, lmxspd, nclebd, lmpot, naez

    type(BasisAtom), pointer :: atomdata
    type(DensityResults), pointer :: densities
    type(EnergyResults), pointer :: energies
    type(RadialMeshData), pointer :: mesh

    double precision, intent(in) :: zat(:)
    double precision, intent(in) :: cleb(nclebd)
    double precision, intent(in) :: dfac(0:lpot,0:lpot)
    integer, intent(in) :: icleb(nclebd,3)
    integer, intent(in) :: loflm(lmxspd)
    integer, intent(in) :: atoms_per_proc

    integer :: i2, ilm
    integer ierr, mylrank, communicator
    double precision :: cmom_save((lpot+1)**2)
    double precision :: cminst_save((lpot+1)**2)
    integer :: num_local_atoms
    integer :: ilocal
    integer :: ilocal2
    integer :: root
    integer :: lmpotd
    
    lmpotd = (lpot+1)**2

    num_local_atoms = calc_data%num_local_atoms

    do ilocal = 1, num_local_atoms
      energies => getEnergies(calc_data, ilocal)
      energies%AC_madelung = 0.d0
    enddo

    !===== begin loop over all atoms =========================================
    ! O(N**2) in calculation and communication !!!

    do I2 = 1, NAEZ

      ! use omp single for MPI part?
      ! if MPI comm. dominates over calculation, OpenMP won't do much

      root = (I2 - 1) / atoms_per_proc
      ilocal2 = mod( (I2 - 1), atoms_per_proc ) + 1

      !-------- bcast information of cmom and cminst begin --------------------

      if (MYLRANK == root) then

        CHECKASSERT(I2 == calc_data%atom_ids(ilocal2))

        densities => getDensities(calc_data, ilocal2)

        do ILM = 1, LMPOT
          CMOM_SAVE(ILM) = densities%CMOM(ILM)
          CMINST_SAVE(ILM) = densities%CMINST(ILM)
        enddo
      else
        CMOM_SAVE = 0.d0
        CMINST_SAVE = 0.d0
      endif

      call MPI_Bcast(CMOM_SAVE,   LMPOTD, MPI_DOUBLE_PRECISION, root, communicator, IERR)
      call MPI_Bcast(CMINST_SAVE, LMPOTD, MPI_DOUBLE_PRECISION, root, communicator, IERR)

      !-------- bcast information of cmom and cminst end----------------------


      ! --> calculate avmad(lm1,lm2)

      do ilocal = 1, num_local_atoms  ! no OpenMP ?
        energies => getEnergies(calc_data, ilocal)

        call sumAC(energies%AC_madelung, CMOM_SAVE, CMINST_SAVE, ZAT(I2), &
                   calc_data%madelung_sum_a(ilocal)%SMAT(:,I2), &
                   LPOT, CLEB, ICLEB, IEND, LOFLM, DFAC)

      enddo ! ilocal

    enddo ! I2
    !===== endloop over all atoms =========================================

    !     contributions are accumulated in AC !!!


    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)
      energies => getEnergies(calc_data, ilocal)
      mesh => atomdata%mesh_ptr

      call addPot(atomdata%potential%VONS, energies%VMAD, &
                  energies%AC_madelung, atomdata%potential%LPOT, &
                  mesh%R, mesh%IRCUT, mesh%IPAN, atomdata%potential%NSPIN)

    enddo ! ilocal

  endsubroutine ! VMADELBLK_new2_com

  !----------------------------------------------------------------------------
  subroutine addPot(vons, vmad, ac, lpot, r, ircut, ipan, nspin)
    double precision, intent(inout) :: vons(:,:,:)
    double precision, intent(inout) :: vmad
    double precision, intent(in)    :: ac(:)
    integer, intent(in) :: lpot
    double precision, intent(in) :: r(:)
    integer, intent(in) :: ircut(0:)
    integer, intent(in) :: ipan
    integer, intent(in) :: nspin

    integer :: irs1
    integer :: l, m, lm, i
    integer :: ispin
    double precision :: pi

    pi = 4.d0*atan(1.d0)
    vmad = 0.d0

    irs1 = ircut(ipan)

    do l = 0, lpot
      do m = -l, l
        lm = l*l + l + m + 1

        if (lm == 1) vmad = ac(1)/sqrt(4.d0*pi)

        !---> add to v the intercell-potential

        do ispin = 1,nspin

          !---> in the case of l=0 : r(1)**l is not defined

          if (l == 0) vons(1,1,ispin) = vons(1,1,ispin) + ac(lm)
          do i = 2, irs1
            vons(i,lm,ispin) = vons(i,lm,ispin) + (-r(i))**l*ac(lm)
          enddo ! i
        enddo ! ispin

      enddo ! m
    enddo ! l

  endsubroutine ! add


  !------------------------------------------------------------------------------
  subroutine sumAC(ac, cmom_save, cminst_save, zat_i2, smat_i2, lpot, cleb, icleb, iend, loflm, dfac)
    double precision, intent(inout) :: ac(:)
    double precision, intent(in) :: cmom_save(:)
    double precision, intent(in) :: cminst_save(:)
    double precision, intent(in) :: zat_i2
    double precision, intent(in) :: smat_i2(:)
    double precision, intent(in) :: cleb(:)
    integer, intent(in) :: iend
    double precision, intent(in) :: dfac(0:,0:)
    integer, intent(in) :: icleb(:,:)
    integer, intent(in) :: loflm(:)
    integer, intent(in) :: lpot

    integer :: lm, lm1, lm2, lm3, l1, l2, i, lmmax, lmpot
    double precision :: pi, fpi
    double precision :: avmad((lpot+1)**2,(lpot+1)**2)
    double precision :: bvmad((lpot+1)**2)

    pi = 4.d0*atan(1.d0)
    fpi = 4.d0*pi

    lmpot = (lpot+1)**2
    lmmax = lmpot

    avmad(:,:) = 0.d0

    do i = 1, iend
      lm1 = icleb(i,1)
      lm2 = icleb(i,2)
      lm3 = icleb(i,3)
      l1 = loflm(lm1)
      l2 = loflm(lm2)

      ! --> this loop has to be calculated only for l1+l2=l3

      avmad(lm1,lm2) = avmad(lm1,lm2) + 2.d0*dfac(l1,l2)*smat_i2(lm3)*cleb(i)
    enddo ! i

    ! --> calculate bvmad(lm1)

    bvmad(:) = 0.d0

    do lm1 = 1, lmpot
      l1 = loflm(lm1)
      bvmad(lm1) = bvmad(lm1) - 2.d0*fpi/dble(2*l1+1)*smat_i2(lm1)
    enddo ! lm1

    do lm = 1, lmpot
      ac(lm) = ac(lm) + bvmad(lm)*zat_i2

      !---> take moments of sphere

      do lm2 = 1, lmmax
        ac(lm) = ac(lm) + avmad(lm,lm2)*(cmom_save(lm2) + cminst_save(lm2))
      enddo ! lm2
    enddo ! lm

  endsubroutine ! sum

endmodule ! MadelungPotential_mod
