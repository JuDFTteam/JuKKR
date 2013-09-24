module main2_aux_mod
  implicit none

  contains

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

    integer :: LRECRES1

    LRECRES1 = 8*43 + 16*(LMAXD+2)
    if (NPOL==0) then
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

    if (NPOL==0) then
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

  !-------------------------------------------------------------------------
  !> Checks for abort flag on rank 0.
  !>
  !> If rank 0 passes flag=1 this function returns .true. on all ranks,
  !> otherwise .false.
  !> The value of flag passed by ranks other than rank 0 is ignored.
  logical function is_abort_by_rank0(flag, communicator)
    implicit none

    integer, intent(in) :: flag
    integer, intent(in) :: communicator

    include 'mpif.h'
    integer :: ierr
    integer :: master_rank
    integer :: stop_integer

    master_rank = 0
    stop_integer = flag

    call MPI_BCAST(stop_integer,1,MPI_INTEGER, master_rank,communicator,ierr)

    is_abort_by_rank0 = (stop_integer == 1)

  end function

  !---------------------------------------------------------------------------
  !> Checks for file 'STOP' in current working directory, returns 1 if it
  !> exists, otherwise 0.
  !>
  !> Should be called only by master rank!!!
  integer function stopfile_flag()
    implicit none
    logical :: stopfile_exists
    stopfile_exists = .false.
    inquire(file='STOP',exist=stopfile_exists)
    stopfile_flag = 0
    if (stopfile_exists) then
      stopfile_flag = 1
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

  !----------------------------------------------------------------------------
  !> correct Fermi-energy (charge neutrality).
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
