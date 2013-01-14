module MadelungPotential_mod

  private :: sumAC

CONTAINS

  !----------------------------------------------------------------------------
  !> Add Madelung potential to VONS.
  !> Needs SMAT (Lattice sums from calculateMadelungLatticeSum)
  !> principal input: CMOM, CMINST, SMAT, VONS --> VONS (changed)
  !> Wrapper for VMADELBLK
  subroutine addMadelungPotential_com(calc_data, ZAT, rank, atoms_per_proc, &
                                      communicator, comm_size)

    use CalculationData_mod
    use MadelungCalculator_mod
    implicit none

    type (CalculationData), intent(inout) :: calc_data

    double precision, intent(in) ::  ZAT(:)

    integer, intent(in) :: rank
    integer, intent(in) :: atoms_per_proc
    integer, intent(in) :: communicator
    integer, intent(in) :: comm_size
    !------------------------------------

    type (MadelungCalculator), pointer :: madelung_calc
    integer :: naez

    madelung_calc => getMadelungCalculator(calc_data)
    naez = size(ZAT)

    call VMADELBLK_new2_com(calc_data,madelung_calc%LPOT,naez,ZAT, &
         madelung_calc%LMPOTD,madelung_calc%clebsch%CLEB,madelung_calc%clebsch%ICLEB, &
         madelung_calc%clebsch%IEND, &
         madelung_calc%LMXSPD,madelung_calc%clebsch%NCLEBD,madelung_calc%clebsch%LOFLM, &
         madelung_calc%DFAC, &
         rank, atoms_per_proc, &
         communicator,comm_size)

  !      SUBROUTINE VMADELBLK_new_com(CMOM,CMINST,LPOT,NSPIN,NAEZ,
  !     &                     VONS,ZAT,R,
  !     &                     IRCUT,IPAN,VMAD,
  !     &                     LMPOT,SMAT,CLEB,ICLEB,IEND,
  !     &                     LMXSPD,NCLEBD,LOFLM,DFAC,
  !     >                     MYLRANK,
  !     >                     communicator,comm_size,
  !     &                     irmd, ipand)

  end subroutine

  subroutine VMADELBLK_new2_com(calc_data,LPOT,NAEZ, &
  ZAT, &
  LMPOT,CLEB,ICLEB,IEND, &
  LMXSPD,NCLEBD,LOFLM,DFAC, &
  MYLRANK, atoms_per_proc, &
  communicator,comm_size)
    ! **********************************************************************
    !
    !     calculate the madelung potentials and add these to the poten-
    !     tial v  (in the spin-polarized case for each spin-direction
    !     this is the same)
    !     it uses the structure dependent matrices AVMAD and BVMAD which
    !     are calculated once in the subroutine MADELUNG3D
    !     ( may 2004)
    !     the charge-moments are calculated in the subroutine vintras,
    !     therefore vintras has to be called first.
    !     the madelung-potential is expanded into spherical harmonics.
    !     the lm-term of the potential v of the atom i is given by
    !
    !      v(r,lm,i) =  (-r)**l * {avmad(i,i2,lm,l'm')*cmom(i2,l'm')
    !                                               +bvmad(i,i2,lm)*z(i2)}
    !
    !     summed over i2 (all atoms) and l'm'
    !     (see notes by b.drittler)
    !
    !                               b.drittler   nov. 1989
    !
    !     adopted for the case of more atoms on the same site, summation is
    !     done over the occupants of that site, the charge is weighted with
    !     the appropriate concentration of the occupant  v.popescu feb. 2002
    !
    !     impurity-program adopted feb. 2004 (according to n.papanikalou)
    !
    ! **********************************************************************
    use CalculationData_mod
    use MadelungCalculator_mod
    use EnergyResults_mod
    use DensityResults_mod
    use BasisAtom_mod
    use RadialMeshData_mod

    implicit none
    include 'mpif.h'

    type (CalculationData), intent(inout) :: calc_data

    integer IEND,LPOT,LMXSPD,NCLEBD,LMPOT,NAEZ

    type (BasisAtom), pointer :: atomdata
    type (DensityResults), pointer :: densities
    type (EnergyResults), pointer :: energies
    type (MadelungLatticeSum), pointer :: madelung_sum
    type (RadialMeshData), pointer :: mesh

    double precision ZAT(*)
    double precision CLEB(NCLEBD)
    double precision DFAC(0:LPOT,0:LPOT)
    integer ICLEB(NCLEBD,3)
    integer LOFLM(LMXSPD)
    !     ..
    !     .. Local Scalars ..
    integer I2,ILM

    !     .. Local Arrays ..
    !     Fortran 90 automatic arrays

    double precision CMOM_SAVE((LPOT+1)**2)
    double precision CMINST_SAVE((LPOT+1)**2)

    !----- MPI ---------------------------------------------------------------
    integer IERR
    !     .. L-MPI
    integer MYLRANK, communicator, comm_size
    integer, intent(in) :: atoms_per_proc

    integer :: num_local_atoms
    integer :: ilocal
    integer :: ilocal2
    integer :: root

    external MPI_BCAST
    !     .. Intrinsic Functions ..
    intrinsic ATAN,SQRT

    integer LMPOTD
    LMPOTD=(LPOT+1)**2
    !     ..................................................................

    num_local_atoms = getNumLocalAtoms(calc_data)

    do ilocal = 1, num_local_atoms
      energies => getEnergies(calc_data, ilocal)
      energies%AC_madelung = 0.0d0
    end do

    !===== begin loop over all atoms =========================================
    ! O(N**2) in calculation and communication !!!

    do I2 = 1,NAEZ

      root = (I2 - 1) / atoms_per_proc
      ilocal2 = mod( (I2 - 1), atoms_per_proc ) + 1

      !-------- bcast information of cmom and cminst begin --------------------

      if (MYLRANK == root) then

        ! TODO: check indices
        densities => getDensities(calc_data, ilocal2)

        do ILM = 1, LMPOT
          CMOM_SAVE(ILM) = densities%CMOM(ILM)
          CMINST_SAVE(ILM) = densities%CMINST(ILM)
        enddo
      else
        CMOM_SAVE = 0.0d0
        CMINST_SAVE = 0.0d0
      endif

      call MPI_BCAST(CMOM_SAVE,LMPOTD,MPI_DOUBLE_PRECISION, &
      root, communicator,IERR)
      call MPI_BCAST(CMINST_SAVE,LMPOTD,MPI_DOUBLE_PRECISION, &
      root, communicator,IERR)

      !-------- bcast information of cmom and cminst end ----------------------


      ! --> calculate avmad(lm1,lm2)

      do ilocal = 1, num_local_atoms  ! no OpenMP ?
        madelung_sum => getMadelungSum(calc_data, ilocal)
        energies => getEnergies(calc_data, ilocal)

        call sumAC(energies%AC_madelung, CMOM_SAVE, CMINST_SAVE, ZAT(I2), &
                   madelung_sum%SMAT(:,I2), &
                   LPOT, CLEB, ICLEB, IEND, LOFLM, DFAC)

      end do ! ilocal

    end do
    !===== end loop over all atoms =========================================

    !     contributions are accumulated in AC !!!


    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)
      energies => getEnergies(calc_data, ilocal)
      mesh => atomdata%mesh_ptr

      call addPot(atomdata%potential%VONS, energies%VMAD, &
                  energies%AC_madelung, atomdata%potential%LPOT, &
                  mesh%R, mesh%IRCUT, mesh%IPAN, atomdata%potential%NSPIN)

    end do   !ilocal

  end subroutine VMADELBLK_new2_com

  !----------------------------------------------------------------------------
  subroutine addPot(VONS, VMAD, AC, LPOT, R, IRCUT, IPAN, NSPIN)
    implicit none

    double precision, intent(inout) :: VONS(:,:,:)
    double precision, intent(inout) :: VMAD
    double precision, intent(in)    :: AC(:)
    integer, intent(in) :: LPOT
    double precision, intent(in) :: R(:)
    integer, intent(in) :: IRCUT(:)
    integer, intent(in) :: IPAN
    integer, intent(in) :: NSPIN

    integer :: IRS1
    integer :: L, M, LM, I
    integer :: ISPIN
    double precision :: PI

    PI = 4.0D0*atan(1.0D0)
    VMAD = 0.0d0

    IRS1 = IRCUT(IPAN)

    ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
    do L = 0,LPOT
      ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      do M = -L,L
        LM = L*L + L + M + 1

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if ( LM.eq.1 ) VMAD = AC(1)/sqrt(4.D0*PI)

        !---> add to v the intercell-potential

        ! ================================================================= SPIN
        do ISPIN = 1,NSPIN

          !---> in the case of l=0 : r(1)**l is not defined

          if ( L.eq.0 ) VONS(1,1,ISPIN) = VONS(1,1,ISPIN) + AC(LM)
          do I = 2,IRS1
            VONS(I,LM,ISPIN) = VONS(I,LM,ISPIN) + (-R(I))**L*AC(LM)
          end do
        end do
         ! ================================================================= SPIN

      end do
       ! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    end do
    ! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

  end subroutine


  !------------------------------------------------------------------------------
  subroutine sumAC(AC, CMOM_SAVE, CMINST_SAVE, ZAT_I2, SMAT_I2, &
  LPOT, CLEB, ICLEB, IEND, LOFLM, DFAC)

    implicit none

    double precision, intent(inout) :: AC(:)
    double precision, intent(in) :: CMOM_SAVE(:)
    double precision, intent(in) :: CMINST_SAVE(:)
    double precision, intent(in) :: ZAT_I2
    double precision, intent(in) :: SMAT_I2(:)

    double precision, intent(in) :: CLEB(:)
    integer, intent(in) :: IEND
    double precision, intent(in) :: DFAC(:,:)
    integer, intent(in) :: ICLEB(:,:)
    integer, intent(in) :: LOFLM(:)

    integer, intent(in) :: LPOT

    !------------------
    integer :: LM, LM1, LM2, LM3, L1, L2, I
    integer :: LMMAX, LMPOT

    double precision :: PI, FPI

    !     .. Local Arrays ..
    !     Fortran 90 automatic arrays
    double precision AVMAD((LPOT+1)**2,(LPOT+1)**2)
    double precision BVMAD((LPOT+1)**2)

    PI = 4.0D0*atan(1.0D0)
    FPI = 4.0D0*PI

    LMPOT = (LPOT+1)**2
    LMMAX = LMPOT

    do LM1 = 1,LMPOT
      do LM2 = 1,LMPOT
        AVMAD(LM1,LM2) = 0.0D0
      end do
    end do

    do I = 1,IEND
      LM1 = ICLEB(I,1)
      LM2 = ICLEB(I,2)
      LM3 = ICLEB(I,3)
      L1 = LOFLM(LM1)
      L2 = LOFLM(LM2)

      ! --> this loop has to be calculated only for l1+l2=l3

      AVMAD(LM1,LM2) = AVMAD(LM1,LM2) + &
      2.0D0*DFAC(L1,L2)*SMAT_I2(LM3)*CLEB(I)
    end do


    ! --> calculate bvmad(lm1)

    do LM1 = 1,LMPOT
      BVMAD(LM1) = 0.0D0
    end do
    !     IF(NAEZ.GT.1) THEN
    !---> lm = 1 component disappears if there is only one host atom
    !     WRONG E.R.

    do LM1 = 1,LMPOT
      L1 = LOFLM(LM1)
      BVMAD(LM1) = BVMAD(LM1) - &
      2.0D0*FPI/dble(2*L1+1)*SMAT_I2(LM1)
    end do

    !     END IF

    do LM = 1,LMPOT
      !if(I2.eq.1) AC(LM) = 0.0D0
      AC(LM) = AC(LM) + BVMAD(LM)*ZAT_I2

      !---> take moments of sphere

      do LM2 = 1,LMMAX
        AC(LM) = AC(LM) + AVMAD(LM,LM2)*CMOM_SAVE(LM2)
      end do

      do LM2 = 1,LMMAX
        AC(LM) = AC(LM) + AVMAD(LM,LM2)*CMINST_SAVE(LM2)
      end do
    end do

  end subroutine

end module
