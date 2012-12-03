module main2_aux_mod
  implicit none

  contains

  ! from dimension parameters, calculate some derived parameters
  !----------------------------------------------------------------------------
  subroutine getDerivedParameters(IGUESSD, IRMD, IRMIND, IRNSD, &
                                  LMAXD, LMAXD1, LMMAXD, LMPOTD, LMXSPD, &
                                  LRECRES2, MMAXD, NAEZD, &
                                  NGUESSD, NSPIND, NTIRD)
    implicit none
    integer :: IGUESSD
    integer :: IRMD
    integer :: IRMIND
    integer :: IRNSD

    integer :: LMAXD
    integer :: LMAXD1
    integer :: LMMAXD
    integer :: LMPOTD
    integer :: LMXSPD
    integer :: LRECRES2
    integer :: MMAXD
    integer :: NAEZD
    integer :: NGUESSD
    integer :: NSPIND
    integer :: NTIRD

    !--------- local
    integer :: LPOTD

    ! derived dimension parameters
    LPOTD = 2*LMAXD
    LMMAXD= (LMAXD+1)**2

    LMAXD1=LMAXD+1
    MMAXD  = 2*LMAXD + 1
    LMXSPD= (2*LPOTD+1)**2

    LMPOTD= (LPOTD+1)**2
    NTIRD=(IRMD+(IRNSD+1)*(LMPOTD-1))*NSPIND
    IRMIND=IRMD-IRNSD

    NGUESSD = 1 + IGUESSD * ( NAEZD * (LMAXD+1)**2 - 1 )

    ! Record lengths
    LRECRES2=4+8*(NSPIND*(LMAXD+7)+2*LPOTD+4+2)
  end subroutine

  !----------------------------------------------------------------------------
  subroutine consistencyCheck03(ATOM, CLS, EZOA, INDN0, NACLS, NACLSD, NAEZ, NCLSD, NR, NUMN0)
    implicit none
    integer :: ATOM(:,:)
    integer :: CLS(:)
    integer :: EZOA(:,:)
    integer :: INDN0(:,:)
    integer :: NACLS(:)
    integer :: NACLSD
    integer :: NAEZ
    integer :: NCLSD
    integer :: NR
    integer :: NUMN0(:)

    integer :: I1
    integer :: IE

    do I1 = 1, NAEZ
      if (CLS(I1) < 1 .or. CLS(I1) > NCLSD) then
        write (*,*) "main2: CLS defect, site ", I1, " value ", CLS(I1)
        stop
      end if
    end do

    do I1 = 1, NCLSD
      if (NACLS(I1) < 1 .or. NACLS(I1) > NACLSD) then
        write (*,*) "main2: NACLS defect, cluster ", I1, " value ", NACLS(I1)
        stop
      end if
    end do

    do I1 = 1, NAEZ
      do IE = 1, NACLS(CLS(I1))
        if (ATOM(IE, I1) < 1 .or. ATOM(IE, I1) > NAEZ) then
          write (*,*) "main2: ATOM defect value ", IE, I1, ATOM(IE, I1)
          stop
        end if

        ! EZOA - index array to point to lattice vector (0 is allowed value)
        if (EZOA(IE, I1) < 0 .or. EZOA(IE, I1) > NR) then
          write (*,*) "main2: EZOA defect value ", IE, I1, EZOA(IE, I1)
          stop
        end if
      end do
    end do

    do I1 = 1, NAEZ
      if (NUMN0(I1) < 1 .or. NUMN0(I1) > NAEZ) then
        write (*,*) "main2: NUMN0 inconsistent site ", I1, " value ", NUMN0(I1)
        stop
      end if

      do IE = 1, NUMN0(I1)
        if (INDN0(I1, IE) < 1 .or. INDN0(I1, IE) > NAEZ) then
          write (*,*) "main2: INDN0 inconsistent site ", I1, " value ", INDN0(I1, IE)
          stop
        end if
      end do
    end do
  end subroutine

  !----------------------------------------------------------------------------
  ! Read k-mesh file
  subroutine readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)
    implicit none
    double precision :: BZKP(:,:,:)
    integer :: MAXMESH
    integer :: NOFKS(:)
    double precision :: VOLBZ(:)
    double precision :: VOLCUB(:,:)

    ! -----------------------------
    integer :: I
    integer :: ID
    integer :: L

    open (52,file='kpoints',form='formatted')
    rewind (52)

    do L = 1,MAXMESH
      read (52,fmt='(I8,f15.10)') NOFKS(L),VOLBZ(L)
      read (52,fmt=*) (BZKP(ID,1,L),ID=1,3),VOLCUB(1,L)
      do I=2,NOFKS(L)
        read (52,fmt=*) (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L)
      end do
    end do

    close (52)
  end subroutine

  !----------------------------------------------------------------------------
  !> Print info about Energy-Point currently treated.
  !>
  subroutine printEnergyPoint(EZ_point, IE, ISPIN, NMESH)
    implicit none
    double complex :: EZ_point
    integer :: IE
    integer :: ISPIN
    integer :: NMESH
    write (6,'(A,I3,A,2(1X,F10.6),A,I3,A,I3)')  &
    ' ** IE = ',IE,' ENERGY =',EZ_point, &
    ' KMESH = ', NMESH,' ISPIN = ',ISPIN
  end subroutine

  !----------------------------------------------------------------------------
  !> Print Fermi-Energy information to screen.
  subroutine printFermiEnergy(DENEF, E2, E2SHIFT, EFOLD, NAEZ)
    implicit none
    double precision :: DENEF
    double precision :: E2
    double precision :: E2SHIFT
    double precision :: EFOLD
    integer :: NAEZ

    write (6,fmt=9020) EFOLD,E2SHIFT

    ! --> divided by NAEZ because the weight of each atom has been already
    !     taken into account in 1c

    write (6,fmt=9030) E2,DENEF/DBLE(NAEZ)
    write(6,'(79(1H+),/)')
9020 format ('                old', &
    ' E FERMI ',F12.6,' Delta E_F = ',f12.6)
9030 format ('                new', &
    ' E FERMI ',F12.6,'  DOS(E_F) = ',f12.6)
  end subroutine

  !----------------------------------------------------------------------------
  !> Print total number of Iterative solver iterations
  subroutine printSolverIterationNumber(ITER, NOITER_ALL)
    implicit none
    integer :: ITER
    integer :: NOITER_ALL

    write(6,'(79(1H=))')
    write(6,'(19X,A,I3,A,I10)') '       ITERATION : ', &
    ITER,' SUM of QMR ',NOITER_ALL
    write(6,'(79(1H=),/)')
  end subroutine

  !----------------------------------------------------------------------------
  !> Print timing information for SCF iteration.
  subroutine writeIterationTimings(ITER, TIME_I, TIME_S)
    implicit none
    integer :: ITER
    double precision :: TIME_I
    double precision :: TIME_S

    call OUTTIME(.true. ,'end .................', &
    TIME_I,ITER)
    write(6,'(79(1H=))')
    write(2,'(79(1H=))')
    call OUTTIME(.true. ,'finished in .........', &
    TIME_S,ITER)
    write(2,'(79(1H=))')
    write(6,'(79(1H=),/)')
  end subroutine

  !----------------------------------------------------------------------------
  !> Open Results2 File
  subroutine openResults2File(LRECRES2)
    implicit none
    integer :: LRECRES2
    open (72,access='direct',recl=LRECRES2,file='results2', &
    form='unformatted')
  end subroutine

  !----------------------------------------------------------------------------
  !> Write calculated stuff into historical 'results2' file
  subroutine writeResults2File(CATOM, ECOU, EDCLDAU, EPOTIN, ESPC, ESPV, EULDAU, EXC, I1, LCOREMAX, VMAD)
    implicit none
    double precision :: CATOM(:)
    double precision :: ECOU(:)
    double precision :: EDCLDAU
    double precision :: EPOTIN
    double precision :: ESPC(:,:)
    double precision :: ESPV(:,:)
    double precision :: EULDAU
    double precision :: EXC(:)
    integer :: I1
    integer :: LCOREMAX
    double precision :: VMAD

    write(72,rec=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX, &
                     EULDAU,EDCLDAU
  end subroutine

  !----------------------------------------------------------------------------
  !> Close the file 'results2'
  subroutine closeResults2File()
    implicit none
    close(72)
  end subroutine

  !----------------------------------------------------------------------------
  !> Open the file 'results1'
  subroutine openResults1File(IEMXD, LMAXD, NPOL)
    implicit none
    integer :: IEMXD
    integer :: LMAXD
    integer :: NPOL
    logical :: TEST

    integer :: LRECRES1

    LRECRES1 = 8*43 + 16*(LMAXD+2)
    if (NPOL==0 .or. TEST('DOS     ')) then
      LRECRES1 = LRECRES1 + 32*(LMAXD+2)*IEMXD
    end if

    open (71,access='direct',recl=LRECRES1,file='results1', &
    form='unformatted')
  end subroutine

  !----------------------------------------------------------------------------
  !> Write some stuff to the 'results1' file
  subroutine writeResults1File(CATOM, CHARGE, DEN, ECORE, I1, NPOL, QC)
    implicit none
    double precision :: CATOM(:)
    double precision :: CHARGE(:,:)
    double complex :: DEN(:,:,:)
    double precision :: ECORE(20,2)
    integer :: I1
    integer :: NPOL
    double precision :: QC

    logical :: TEST

    if (NPOL==0 .or. TEST('DOS     ')) then
      write(71,rec=I1) QC,CATOM,CHARGE,ECORE,DEN  ! write density of states (DEN) only when certain options set
    else
      write(71,rec=I1) QC,CATOM,CHARGE,ECORE
    end if
  end subroutine

  !---------------------------------------------------------------------------
  !> Closes the file 'results1'
  subroutine closeResults1File()
    implicit none
    close(71)
  end subroutine

  !----------------------------------------------------------------------------
  !> Calculate \Delta T_up - T_down for exchange couplings calculation.
  !> The result is stored in DTIXIJ(:,:,1)
  subroutine calcDeltaTupTdown(DTIXIJ)
    implicit none
    double complex, intent(inout) :: DTIXIJ(:,:,:)
    integer :: LMMAXD

    integer :: LM1
    integer :: LM2

    lmmaxd = size(DTIXIJ,1)

    do LM2 = 1,LMMAXD
      do LM1 = 1,LMMAXD
        DTIXIJ(LM1,LM2,1) = DTIXIJ(LM1,LM2,2) - DTIXIJ(LM1,LM2,1)
      enddo
    enddo

  end subroutine

  !----------------------------------------------------------------------------
  !> Substract diagonal reference T matrix of certain spin channel
  !> from real system's T matrix
  subroutine substractReferenceTmatrix(TMATN, TREFLL, LMMAXD)
    implicit none
    integer :: LM1
    integer :: LMMAXD
    double complex :: TMATN(:,:)
    double complex :: TREFLL(:,:)

    ! Note: TREFLL is diagonal! - spherical reference potential
    do LM1 = 1,LMMAXD
      TMATN(LM1,LM1) =  TMATN(LM1,LM1) - TREFLL(LM1,LM1)
    end do

  end subroutine

  !> Check if file 'STOP' exists. If yes, tell all ranks
  !> in communicator to abort by returning .true.
  !> The idea is to let only one rank to inquire if
  !> 'STOP' exists.
  logical function isManualAbort_com(rank, communicator)
    implicit none

    integer, intent(in) :: rank
    integer, intent(in) :: communicator

    integer :: ierr
    integer :: stop_integer
    integer :: master_rank
    logical :: STOPIT

    include 'mpif.h'

    isManualAbort_com = .false.
    stop_integer = 0
    master_rank = 0

    if (rank == master_rank) then
      inquire(file='STOP',exist=STOPIT)
      if (STOPIT) stop_integer = 1
    end if

    call MPI_BCAST(stop_integer,1,MPI_INTEGER, &
    master_rank,communicator,ierr)

    if (stop_integer == 1) then
      isManualAbort_com = .true.
    end if

  end function

  !----------------------------------------------------------------------------
  !> Prints a double line separator. (=========)
  !> @param unit_number  optional: write to file 'unit_number' otherwise to
  !>                               stdout
  subroutine printDoubleLineSep(unit_number)
    implicit none
    integer, intent(in), optional :: unit_number

    if (.not. present(unit_number)) then
      write(*,'(79("="))')
    else
      write(unit_number,'(79("="))')
    end if
  end subroutine

  !----------------------------------------------------------------------------
  !> Communicate and sum up contributions for charge neutrality and
  !> density of states at Fermi level.
  subroutine sumNeutralityDOSFermi_com(CHRGNT, DENEF, communicator)
    implicit none
    include 'mpif.h'
    double precision, intent(inout) :: CHRGNT
    double precision, intent(inout) :: DENEF
    integer, intent(in) :: communicator
    !----------------------------------

    double precision :: WORK1(2)
    double precision :: WORK2(2)
    integer :: ierr

    WORK1(1) = CHRGNT
    WORK1(2) = DENEF

    call MPI_ALLREDUCE(WORK1,WORK2,2,MPI_DOUBLE_PRECISION,MPI_SUM, &
    communicator,IERR)

    CHRGNT = WORK2(1)
    DENEF  = WORK2(2)

  end subroutine

  !----------------------------------------------------------------------------
  !> Extracts the DOS at the Fermi level and stores result in DENEF
  !> @param IELAST energy index of Fermi level
  !> @param LMAXD1 lmax+1
  double precision function calcDOSatFermi(DEN, IELAST, IEMXD, LMAXD1, NSPIN)
    implicit none

    double complex, intent(in) :: DEN(0:LMAXD1,IEMXD,NSPIN)
    integer, intent(in) :: IELAST
    integer, intent(in) :: IEMXD
    integer, intent(in) :: LMAXD1
    integer, intent(in) :: NSPIN

    integer :: ISPIN
    integer :: L
    double precision :: PI
    double precision :: DENEF

    PI = 4.0D0*ATAN(1.0D0)

    ! get density of states at Fermi-level
    DENEF = 0.0d0
    do ISPIN = 1,NSPIN
      do L = 0,LMAXD1
        DENEF = DENEF - 2.0D0 * &
        DIMAG(DEN(L,IELAST,ISPIN))/PI/DBLE(NSPIN)
      end do
    end do

    calcDOSatFermi = DENEF

  end function calcDOSatFermi

  !----------------------------------------------------------------------------
  !> Copy output potential to input potential for new iteration.
  subroutine resetPotentials(IRC1, IRMD, IRMIN1, IRMIND, LMPOTD, NSPIN, VINS, VISP, VONS)
    implicit none

    integer :: IRC1
    integer :: IRMIN1
    integer :: ISPIN
    integer :: J
    integer :: LM

    integer :: IRMD

    integer :: IRMIND
    integer :: LMPOTD
    integer :: NSPIN

    double precision, intent(out) :: VINS(IRMIND:IRMD,LMPOTD,2)
    double precision, intent(out) :: VISP(IRMD,2)
    double precision, intent(inout) :: VONS(IRMD,LMPOTD,2)

    !DEBUG
    if (IRMIN1 < IRMIND) then
      write(*,*) "resetPotentials: IRMIN1 < IRMIND"
      stop
    end if

    !DEBUG
    if (IRC1 > IRMD) then
      write(*,*) "resetPotentials: IRC1 > IRMD"
      stop
    end if

    !initialise VINS
    do ISPIN = 1,2
      do LM = 1, LMPOTD
        do J = IRMIND, IRMD
          VINS(J,LM,ISPIN) = 0.0D0
        enddo
      enddo
    enddo

    !initialise VISP
    do ISPIN = 1,2
      do J = 1, IRMD
        VISP(J,ISPIN) = 0.0D0
      enddo
    enddo

    ! copy output potential to input potential for new iteration
    do ISPIN = 1,NSPIN

      call DCOPY(IRC1,VONS(1,1,ISPIN),1,VISP(1,ISPIN),1)

      if (LMPOTD>1) then

        do LM = 2,LMPOTD
          do J = IRMIN1,IRC1
            VINS(J,LM,ISPIN) = VONS(J,LM,ISPIN)
          end do
        end do
      end if
    enddo
  end subroutine resetPotentials

  !----------------------------------------------------------------------------
  !> Calculates the l-resolved charges.
  !> This is done by energy integration in the complex plane over the imaginary
  !> part of the diagonal of the structural Green's function (DEN)

  subroutine calcChargesLres(CHARGE, DEN, IELAST, LMAXD1, NSPIN, WEZ, IEMXD)
    implicit none

    double precision, intent(out) :: CHARGE(0:LMAXD1,2)
    double complex, intent(in) :: DEN(0:LMAXD1,IEMXD,NSPIN)
    doublecomplex, intent(in) :: WEZ(IEMXD)

    integer, intent(in) :: IEMXD
    integer, intent(in) :: NSPIN
    integer, intent(in) :: LMAXD1
    integer, intent(in) :: IELAST

    integer :: ISPIN
    integer :: L
    integer :: IE

    ! ---> l/m_s/atom-resolved charges

    do ISPIN = 1,NSPIN
      do L = 0,LMAXD1
        CHARGE(L,ISPIN) = 0.0D0

        do IE = 1,IELAST
          CHARGE(L,ISPIN) = CHARGE(L,ISPIN) + &
          DIMAG(WEZ(IE)*DEN(L,IE,ISPIN))/ &
          DBLE(NSPIN)
        end do

      end do
    end do
  end subroutine calcChargesLres


  !--------------------------------------------------------------------------------------
  ! read input generated by kkr0 from file inpn.unf
  subroutine readKKR0InputNew(NSYMAXD, ALAT, ATOM, BCP, BRAVAIS, &
                           CLS, DSYMLL, EZOA, FCM, GMAX, ICST, &
                           IGUESS, IMIX, INDN0, &
                           ISYMINDEX, &
                           JIJ, KFORCE, KMESH, KPRE, KTE, KXC, &
                           LDAU, MAXMESH, &
                           MIXING, NACLS, NCLS, NR, NREF, &
                           NSRA, NSYMAT, NUMN0, OPTC, QMRBOUND, &
                           RBASIS, RCLS, RCUTJIJ, REFPOT, RMAX, RMTREF, &
                           RR, SCFSTEPS, TESTC, VREF, ZAT)

    implicit none
    integer :: NSYMAXD
    double precision :: ALAT
    integer, allocatable :: ATOM(:,:)
    integer :: BCP
    double precision :: BRAVAIS(3,3)
    integer, allocatable :: CLS(:)
    doublecomplex, allocatable :: DSYMLL(:,:,:)
    integer, allocatable :: EZOA(:,:)
    double precision :: FCM
    double precision :: GMAX
    integer :: ICST
    integer :: IGUESS
    integer :: IMIX
    integer, allocatable :: INDN0(:,:)
    integer :: ISHIFT
    integer :: ISYMINDEX(NSYMAXD)
    logical :: JIJ
    integer :: KFORCE
    integer, allocatable :: KMESH(:)
    integer :: KPRE
    integer :: KTE
    integer :: KVMADdummy
    integer :: KXC
    logical :: LDAU
    integer :: MAXMESH
    double precision :: MIXING
    integer, allocatable :: NACLS(:)
    integer :: NCLS
    integer :: NR
    integer :: NREF
    integer :: NSRA
    integer :: NSYMAT

    integer, allocatable :: NUMN0(:)
    character(len=8) :: OPTC(8)
    double precision :: QMRBOUND
    double precision, allocatable :: RBASIS(:,:)
    double precision, allocatable :: RCLS(:,:,:)
    double precision :: RCUTJIJ
    integer, allocatable :: REFPOT(:)
    double precision :: RMAX
    double precision, allocatable :: RMTREF(:)
    double precision, allocatable :: RR(:,:)
    integer :: SCFSTEPS
    character(len=8) :: TESTC(16)

    double precision, allocatable :: VREF(:)
    double precision, allocatable :: ZAT(:)

    ! ======================================================================
    ! =             read in variables from unformatted files               =
    ! ======================================================================
    open (67,file='inpn.unf',form='unformatted')
    read (67) KMESH
    read (67) MAXMESH
    read (67) RR
    read (67) EZOA
    read (67) NUMN0
    read (67) INDN0
    read (67) NSYMAT
    read (67) DSYMLL
    read (67) IMIX
    read (67) MIXING
    read (67) FCM
    read (67) KPRE
    read (67) KTE
    read (67) KVMADdummy
    read (67) KXC
    read (67) ISHIFT
    read (67) KFORCE
    read (67) IGUESS
    read (67) BCP
    read (67) QMRBOUND
    read (67) ZAT
    read (67) ALAT
    read (67) TESTC
    read (67) OPTC
    read (67) BRAVAIS
    read (67) RBASIS
    read (67) RMAX
    read (67) GMAX
    read (67) NSRA
    read (67) ICST
    read (67) NCLS
    read (67) NREF
    read (67) RMTREF
    read (67) VREF
    read (67) RCLS
    read (67) ATOM
    read (67) CLS
    read (67) NACLS
    read (67) REFPOT
    read (67) NR
    read (67) RCUTJIJ
    read (67) JIJ
    read (67) LDAU
    read (67) ISYMINDEX
    read (67) SCFSTEPS
    close (67)
  end subroutine

  subroutine doFermiEnergyCorrection(atomdata, output, naez, max_shift, CHRGNT, DENEF, R2NEF, &
                                     ESPV, RHO2NS, E2)
    use BasisAtom_mod
    use RadialMeshData_mod
    implicit none

    type (BasisAtom), intent(in) :: atomdata
    logical, intent(in) :: output !< output to stdout - yes/no
    integer, intent(in) :: naez
    double precision, intent(in) :: max_shift !< maximally allowed Fermi-Energy shift (good: 0.03d0)
    double precision, intent(in) :: CHRGNT
    double precision, intent(in) :: DENEF
    double precision, dimension(:,:,:), intent(in) ::  R2NEF

    ! in,out - arguments
    double precision, dimension(0:,:), intent(inout) :: ESPV
    double precision, dimension(:,:,:), intent(inout) ::  RHO2NS
    double precision, intent(inout) :: E2

    !-------- locals
    integer :: nspind
    type (RadialMeshData), pointer :: mesh

    double precision :: E2SHIFT
    double precision :: EFold
    double precision :: DF
    double precision :: PI
    integer :: ispin, lm, lmpotd

    PI = 4.0D0*ATAN(1.0D0)

    nspind = atomdata%nspin
    lmpotd = atomdata%potential%lmpot

    mesh => atomdata%mesh_ptr
! --> determine new Fermi level due to valence charge up to
!     old Fermi level E2 and density of states DENEF

        E2SHIFT = CHRGNT/DENEF
        E2SHIFT = DMIN1(DABS(E2SHIFT),max_shift)*DSIGN(1.0D0,E2SHIFT)
        EFOLD = E2

        E2 = E2 - E2SHIFT

        if( output ) then
          call printFermiEnergy(DENEF, E2, E2SHIFT, EFOLD, NAEZ)
        end if

! ----------------------------------------------------------------------
        DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIND)
! ----------------------------------------------------------------------

        do ISPIN = 1,NSPIND

! -->     get correct density and valence band energies

          ESPV(0,ISPIN) = ESPV(0,ISPIN) - &
          EFOLD*CHRGNT/DBLE(NSPIND*NAEZ)

          do LM = 1,LMPOTD
            call DAXPY(mesh%IRC,DF,R2NEF(1,LM,ISPIN),1, &
            RHO2NS(1,LM,ISPIN),1)
          end do
        end do
! ----------------------------------------------------------------------
  end subroutine

end module main2_aux_mod
