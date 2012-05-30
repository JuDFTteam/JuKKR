module main2_aux_mod
  implicit none

  contains

  ! from dimension parameters, calculate some derived parameters
  !----------------------------------------------------------------------------
  subroutine getDerivedParameters(IGUESSD, IRMD, IRMIND, IRNSD, LASSLD, LM2D, &
                                  LMAXD, LMAXD1, LMMAXD, LMPOTD, LMXSPD, &
                                  LRECRES2, MMAXD, NAEZD, NCLEB, &
                                  NGUESSD, NPOTD, NSPIND, NTIRD)
    implicit none
    integer :: IGUESSD
    integer :: IRMD
    integer :: IRMIND
    integer :: IRNSD
    integer :: LASSLD
    integer :: LM2D
    integer :: LMAXD
    integer :: LMAXD1
    integer :: LMMAXD
    integer :: LMPOTD
    integer :: LMXSPD
    integer :: LRECRES2
    integer :: MMAXD
    integer :: NAEZD
    integer :: NCLEB
    integer :: NGUESSD
    integer :: NPOTD
    integer :: NSPIND
    integer :: NTIRD

    !--------- local
    integer :: LPOTD

    ! derived dimension parameters
    LPOTD = 2*LMAXD
    LMMAXD= (LMAXD+1)**2
    NPOTD=NSPIND*NAEZD
    LMAXD1=LMAXD+1
    MMAXD  = 2*LMAXD + 1
    LM2D= (2*LMAXD+1)**2
    LMXSPD= (2*LPOTD+1)**2
    LASSLD=4*LMAXD
    LMPOTD= (LPOTD+1)**2
    NTIRD=(IRMD+(IRNSD+1)*(LMPOTD-1))*NSPIND
    IRMIND=IRMD-IRNSD

    NGUESSD = 1 + IGUESSD * ( NAEZD * (LMAXD+1)**2 - 1 )
    NCLEB = (LMAXD*2+1)**2 * (LMAXD+1)**2

    ! Record lengths
    LRECRES2=4+8*(NSPIND*(LMAXD+7)+2*LPOTD+4+2)
  end subroutine

  ! Consistency checks
  !----------------------------------------------------------------------------
  subroutine consistencyCheck01(IEMXD, LMAXD, NSPIND, SMPID)
    implicit none
    integer :: IEMXD
    integer :: LMAXD
    integer :: NSPIND
    integer :: SMPID

    ! -------------------------------------------------------------------------
    ! consistency checks
    if (IEMXD < 1) then
      write (*,*) "main2: IEMXD must be >= 1"
      stop
    end if

    if (LMAXD < 0) then
      write (*,*) "main2: LMAXD must be >= 0"
      stop
    end if

    if (SMPID /= 1 .and. SMPID /=2) then
      write (*,*) "main2: SMPID must be 1 or 2"
      stop
    end if

    if (NSPIND /= 1 .and. NSPIND /=2) then
      write (*,*) "main2: NSPIND must be 1 or 2"
      stop
    end if
  end subroutine

  !----------------------------------------------------------------------------
  subroutine consistencyCheck02(IELAST, IEMXD, IGUESS, IGUESSD, LMAX, LMAXD, &
                                NAEZ, NAEZD, NPNT1, NPNT2, NPNT3, NPOL, &
                                NR, NRD, NSPIN, NSPIND)
    implicit none
    integer :: IELAST
    integer :: IEMXD
    integer :: IGUESS
    integer :: IGUESSD
    integer :: LMAX
    integer :: LMAXD
    integer :: NAEZ
    integer :: NAEZD
    integer :: NPNT1
    integer :: NPNT2
    integer :: NPNT3
    integer :: NPOL
    integer :: NR
    integer :: NRD
    integer :: NSPIN
    integer :: NSPIND

    ! -------------- Consistency checks -------------------------------
    if (IGUESS /= IGUESSD .or. IGUESS < 0 .or. IGUESS > 1) then
      write (*,*) "main2: IGUESSD IGUESS inconsistent ", IGUESSD, IGUESS
      stop
    end if

    if (NSPIN /= NSPIND) then
      write (*,*) "main2: NSPIN /= NSPIND"
      stop
    end if

    if (NAEZ /= NAEZD) then
      write (*,*) "main2: NAEZ /= NAEZD"
      stop
    end if

    if (IEMXD /= IELAST) then
      write (*,*) "main2: IEMXD /= IELAST"
      stop
    end if

    if (NPOL /= 0) then
      if (NPNT1 + NPNT2 + NPNT3 + NPOL /= IELAST) then
        write(*,*) "main2: Energy point numbers inconsistent."
      end if
    end if

    if (LMAX /= LMAXD) then
      write (*,*) "main2: LMAX /= LMAXD"
      stop
    end if

    if (NSPIN /= NSPIND) then
      write (*,*) "main2: NSPIN /= NSPIND"
      stop
    end if

    if (NR > NRD .or. NR < 1) then
      write (*,*) "main2: NR inconsistent NR NRD ", NR, NRD
      stop
    end if
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


  !--------------------------------------------------------------------------------------
  ! read input generated by kkr0 from file inp.unf
  subroutine readKKR0Input(NSYMAXD, A, ALAT, ATOM, B, BCP, BRAVAIS, CLEB1C, &
                           CLS, DRDI, DSYMLL, EZOA, FCM, GMAX, GSH, ICLEB1C, ICST, &
                           IEND1, IFUNM, IGUESS, ILM, IMAXSH, IMIX, IMT, INDN0, IPAN, &
                           IRC, IRCUT, IRMIN, IRNS, IRWS, ISHIFT, ISYMINDEX, ITITLE, &
                           JEND, JIJ, KFORCE, KMESH, KPRE, KTE, KVMAD, KXC, LCORE, &
                           LDAU, LLMSP, LMSP, LOFLM1C, MAXMESH, &
                           MIXING, NACLS, NCLS, NCORE, NFU, NR, NREF, &
                           NSRA, NSYMAT, NTCELL, NUMN0, OPTC, QBOUND, QMRBOUND, R, &
                           RBASIS, RCLS, RCUTJIJ, RECBV, REFPOT, RMAX, RMT, RMTREF, &
                           RR, RWS, SCFSTEPS, TESTC, THETAS, VOLUME0, VREF, ZAT)

    implicit none
    integer :: NSYMAXD
    double precision, allocatable :: A(:)
    double precision :: ALAT
    integer, allocatable :: ATOM(:,:)
    double precision, allocatable :: B(:)
    integer :: BCP
    double precision :: BRAVAIS(3,3)
    double precision, allocatable :: CLEB1C(:,:)
    integer, allocatable :: CLS(:)
    double precision, allocatable :: DRDI(:,:)
    doublecomplex, allocatable :: DSYMLL(:,:,:)
    integer, allocatable :: EZOA(:,:)
    double precision :: FCM
    double precision :: GMAX
    double precision, allocatable :: GSH(:)
    integer, allocatable :: ICLEB1C(:,:)
    integer :: ICST
    integer :: IEND1
    integer, allocatable :: IFUNM(:,:)
    integer :: IGUESS
    integer, allocatable :: ILM(:,:)
    integer, allocatable :: IMAXSH(:)
    integer :: IMIX
    integer, allocatable :: IMT(:)
    integer, allocatable :: INDN0(:,:)
    integer, allocatable :: IPAN(:)
    integer, allocatable :: IRC(:)
    integer, allocatable :: IRCUT(:,:)
    integer, allocatable :: IRMIN(:)
    integer, allocatable :: IRNS(:)
    integer, allocatable :: IRWS(:)
    integer :: ISHIFT
    integer :: ISYMINDEX(NSYMAXD)
    integer, allocatable :: ITITLE(:,:)
    integer, allocatable :: JEND(:,:,:)
    logical :: JIJ
    integer :: KFORCE
    integer, allocatable :: KMESH(:)
    integer :: KPRE
    integer :: KTE
    integer :: KVMAD
    integer :: KXC
    integer, allocatable :: LCORE(:,:)
    logical :: LDAU
    integer, allocatable :: LLMSP(:,:)
    integer, allocatable :: LMSP(:,:)
    integer, allocatable :: LOFLM1C(:)
    integer :: MAXMESH
    double precision :: MIXING
    integer, allocatable :: NACLS(:)
    integer :: NCLS
    integer, allocatable :: NCORE(:)
    integer, allocatable :: NFU(:)
    integer :: NR
    integer :: NREF
    integer :: NSRA
    integer :: NSYMAT
    integer, allocatable :: NTCELL(:)
    integer, allocatable :: NUMN0(:)
    character(len=8) :: OPTC(8)
    double precision :: QBOUND
    double precision :: QMRBOUND
    double precision, allocatable :: R(:,:)
    double precision, allocatable :: RBASIS(:,:)
    double precision, allocatable :: RCLS(:,:,:)
    double precision :: RCUTJIJ
    double precision :: RECBV(3,3)
    integer, allocatable :: REFPOT(:)
    double precision :: RMAX
    double precision, allocatable :: RMT(:)
    double precision, allocatable :: RMTREF(:)
    double precision, allocatable :: RR(:,:)
    double precision, allocatable :: RWS(:)
    integer :: SCFSTEPS
    character(len=8) :: TESTC(16)
    double precision, allocatable :: THETAS(:,:,:)
    double precision :: VOLUME0
    double precision, allocatable :: VREF(:)
    double precision, allocatable :: ZAT(:)

    !----- unused dummy variables
    integer :: NSPIN
    integer :: LMAX
    integer :: LMPOT
    integer :: NAEZ_dummy
    integer :: LPOT

    ! ======================================================================
    ! =             read in variables from unformatted files               =
    ! ======================================================================
    open (67,file='inp.unf',form='unformatted')
    read (67) KMESH,MAXMESH,RR,EZOA,NUMN0,INDN0,NSYMAT,DSYMLL
    read (67) NAEZ_dummy,NSPIN,IPAN,IRNS,IRCUT,LCORE,NCORE,NTCELL, &
    LMAX,LPOT,LMPOT
    read (67) IMIX,MIXING,QBOUND,FCM,KPRE,KTE, &
    KVMAD,KXC,ISHIFT,KFORCE,IGUESS,BCP,QMRBOUND
    read (67) A,B,DRDI,R,THETAS,ZAT,IMT,IRC, &
    IRMIN,IRWS,RWS,RMT,ITITLE,LLMSP,NFU
    read (67) ALAT,GSH,ILM,IMAXSH,TESTC,OPTC
    read (67) BRAVAIS,RBASIS,RECBV,VOLUME0,RMAX,GMAX
    read (67) IEND1,CLEB1C,ICLEB1C,LOFLM1C,JEND,IFUNM,LMSP,NSRA,ICST
    read (67) NCLS,NREF,RMTREF,VREF,RCLS, &
    ATOM,CLS,NACLS,REFPOT
    read (67) NR,RCUTJIJ,JIJ,LDAU
    read (67) ISYMINDEX,SCFSTEPS
    close (67)
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
  !> Print timing information for SCF iteration
  subroutine writeIterationTimings(ITER, TIME_I, TIME_S)
    implicit none
    integer :: ITER
    real :: TIME_I
    real :: TIME_S

    call OUTTIME(0 ,'end .................', &
    TIME_I,ITER)
    write(6,'(79(1H=))')
    write(2,'(79(1H=))')
    call OUTTIME(0 ,'finished in .........', &
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

  !---------------------------------------------------------------------------
  !> Open file 'vpotnew'
  subroutine openPotentialFile(LMPOTD, IRNSD, IRMD)
    implicit none
    integer, intent(in) :: LMPOTD
    integer, intent(in) :: IRNSD
    integer, intent(in) :: IRMD

    integer :: LRECPOT

    LRECPOT=8*(LMPOTD*(IRNSD+1)+IRMD+20)

    open (66,access='direct',recl=LRECPOT*2,file='vpotnew', &
    form='unformatted')
  end subroutine

  !--------------------------------------------------------------------------
  !> Read potential for atom I1 from 'vpotnew' file.
  subroutine readPotential(I1, VISP, VINS, ECORE)
    implicit none
    double precision, intent(inout) :: ECORE(20,2)
    integer, intent(in) :: I1
    double precision, intent(inout) :: VINS(:,:,:)
    double precision, intent(inout) :: VISP(:,:)

    read(66,rec=I1) VINS,VISP,ECORE

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes file 'vpotnew'.
  subroutine closePotentialFile()
    implicit none

    close(66)
  end subroutine

  !----------------------------------------------------------------------------
  !> Write potential for atom I1 to 'vpotnew' file.
  subroutine writePotential(I1, VISP, VINS, ECORE)
    implicit none
    double precision, intent(inout) :: ECORE(20,2)
    integer, intent(in) :: I1
    double precision, intent(inout) :: VINS(:,:,:)
    double precision, intent(inout) :: VISP(:,:)

    write(66,rec=I1) VINS,VISP,ECORE

  end subroutine

  !----------------------------------------------------------------------------
  !> Calculate \Delta T_up - T_down for exchange couplings calculation.
  !> Call first for ISPIN=1 then for ISPIN=2 with same parameters.
  subroutine calcDeltaTupTdown(DTIXIJ, ISPIN, LMMAXD, TMATN)
    implicit none
    double complex, intent(inout) :: DTIXIJ(:,:)
    integer :: ISPIN
    integer :: LMMAXD
    double complex, intent(in) :: TMATN(:,:,:)

    integer :: LM1
    integer :: LM2

    if (ISPIN==1) then
      do LM1 = 1,LMMAXD
        do LM2 = 1,LMMAXD
          DTIXIJ(LM1,LM2) = TMATN(LM1,LM2,ISPIN)
        enddo
      enddo
    else
      do LM1 = 1,LMMAXD
        do LM2 = 1,LMMAXD
          DTIXIJ(LM1,LM2) = DTIXIJ(LM1,LM2)-TMATN(LM1,LM2,ISPIN)
        enddo
      enddo
    endif
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

end module main2_aux_mod
