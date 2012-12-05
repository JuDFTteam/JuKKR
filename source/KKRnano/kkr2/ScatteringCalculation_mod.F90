#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"

module ScatteringCalculation_mod

public  :: energyLoop
private :: printEnergyPoint
private :: calcDeltaTupTdown
private :: substractReferenceTmatrix

CONTAINS

!> Input:  *) ebalance_handler (properly initialised)
!>         *) jij_data, ldau_data (properly initialised !!!)
!>
!> Output: ebalance_handler (changed, updated timings and process distribution)
!>         arrays - KKRResults (todo)
!>         jij_data results of jij-calculation
!>         ldau_data LDA+U results
!>
!>         FILES WRITTEN:
!>         *) Logfiles (if requested)
!>         *) JIJ-Files (if requested)
!>         *) matrix dump (if requested)
subroutine energyLoop(iter, atomdata, emesh, params, dims, gaunts, &
                      ebalance_handler, my_mpi, arrays, kkr, jij_data, ldau_data)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use TimerMpi_mod
  use EBalanceHandler_mod

  use TEST_lcutoff_mod !TODO: remove

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod
  use KKRresults_mod

  use GauntCoefficients_mod
  use BasisAtom_mod
  use JijData_mod
  use LDAUData_mod
  use EnergyMesh_mod
  use KKRnanoParallel_mod
  use EBalanceHandler_mod

  use KKRnano_Comm_mod

  use wrappers_mod,     only: calctmat_wrapper, calcdtmat_wrapper

  implicit none

  integer, intent(in) :: iter
  type (KKRnanoParallel), intent(in)    :: my_mpi
  type (EBalanceHandler), intent(inout) :: ebalance_handler
  type (BasisAtom), intent(inout )      :: atomdata  ! in only?
  type (EnergyMesh), intent(in)         :: emesh
  type (LDAUData), intent(inout)        :: ldau_data
  type (JijData), intent(inout)         :: jij_data
  type (Main2Arrays), intent(inout)     :: arrays
  type (KKRresults), intent(inout)      :: kkr  ! out only?
  type (DimParams), intent(in)          :: dims
  type (InputParams), intent(in)        :: params
  type (GauntCoefficients), intent(in)  :: gaunts

  !---------- locals ----------------------------
  type (TimerMpi) :: mult_scattering_timer
  type (TimerMpi) :: single_site_timer
  integer :: ie
  integer :: rf
  integer :: ispin
  integer :: prspin
  integer :: nmesh
  integer :: ekm
  integer :: noiter
  logical :: xccpl
  double complex :: JSCAL ! scaling factor for Jij calculation
  integer :: I1

  ! TODO: FIXME to reenable jij-Calculation !!!!!!!!!!
  xccpl = .false. ! TODO !!!!!

  call resetTimer(mult_scattering_timer)
  call stopTimer(mult_scattering_timer)

  call resetTimer(single_site_timer)

  I1 = atomdata%atom_index

  EKM = 0
  prspin = 1
  kkr%noiter = 0

! IE ====================================================================
!     BEGIN do loop over energies (EMPID-parallel)
! IE ====================================================================
  do IE = 1, params%IELAST
! IE ====================================================================
    if (getMyEnergyId(my_mpi)==ebalance_handler%EPROC(IE)) then
! IE ====================================================================
      call startEBalanceTiming(ebalance_handler, IE)

      WRITELOG(2, *) "Working on energy point ", IE

      do RF = 1,params%NREF
        call TREF(emesh%EZ(IE),arrays%VREF(RF),arrays%LMAXD,arrays%RMTREF(RF), &
                  kkr%TREFLL(1,1,RF),kkr%DTREFLL(1,1,RF), dims%LLY)
      end do

      TESTARRAYLOG(3, kkr%TREFLL)
      TESTARRAYLOG(3, kkr%DTREFLL)

      call GREF_com(emesh%EZ(IE),params%ALAT,gaunts%IEND,params%NCLS,arrays%NAEZ, &
                    gaunts%CLEB,arrays%RCLS,arrays%ATOM,arrays%CLS,gaunts%ICLEB, &
                    gaunts%LOFLM,arrays%NACLS, &
                    arrays%REFPOT, &
                    kkr%TREFLL,kkr%DTREFLL,kkr%GREFN,kkr%DGREFN, &
                    kkr%LLY_G0TR(:,IE), &
                    getMyAtomRank(my_mpi),getMySEcommunicator(my_mpi),&
                    getNumAtomRanks(my_mpi), &
                    arrays%lmaxd, arrays%naclsd, gaunts%ncleb, kkr%nrefd, kkr%nclsd, &
                    dims%LLY)

      TESTARRAYLOG(3, kkr%GREFN)
      TESTARRAYLOG(3, kkr%DGREFN)

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN===================================================================
!------------------------------------------------------------------------------
!     beginning of SMPID-parallel section
!------------------------------------------------------------------------------
      spinloop: do ISPIN = 1,dims%NSPIND
        if(isWorkingSpinRank(my_mpi, ispin)) then

          if (dims%SMPID==1) then
            PRSPIN   = ISPIN
          else
            PRSPIN   = 1
          endif

          call CALCTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, gaunts, arrays%TMATN, kkr%TR_ALPH, ldau_data)

          jij_data%DTIXIJ(:,:,ISPIN) = arrays%TMATN(:,:,ISPIN)  ! save t-matrix for Jij-calc.

          if(dims%LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula
            call CALCDTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, gaunts, arrays%DTDE, kkr%TR_ALPH, ldau_data)
          end if

          RF = arrays%REFPOT(I1)
          call substractReferenceTmatrix(arrays%TMATN(:,:,ISPIN), kkr%TREFLL(:,:,RF), kkr%LMMAXD)
          ! do the same for derivative of T-matrix
          call substractReferenceTmatrix(arrays%DTDE(:,:,ISPIN), kkr%DTREFLL(:,:,RF), kkr%LMMAXD)

          ! TMATN now contains Delta t = t - t_ref !!!
          ! DTDE now contains Delta dt !!!

          ! renormalize TR_ALPH
          kkr%TR_ALPH(ISPIN) = kkr%TR_ALPH(ISPIN) - kkr%LLY_G0TR(arrays%CLS(I1), IE)

          NMESH = arrays%KMESH(IE)

          if( getMyAtomRank(my_mpi)==0 ) then
            if (params%KTE >= 0) call printEnergyPoint(emesh%EZ(IE), IE, ISPIN, NMESH)
          end if

          call stopTimer(single_site_timer)
          call resumeTimer(mult_scattering_timer)

! <<>> Multiple scattering part

!                  if (I1 == 1 .and. IE == IELAST .and. ISPIN == 1) then
!                     DEBUG_dump_matrix = .true.
!                  else
!                     DEBUG_dump_matrix = .false.
!                  endif

          TESTARRAYLOG(3, arrays%TMATN(:,:,ISPIN))

          call KLOOPZ1( &
          arrays%GMATN(1,1,1,ISPIN), &
          params%ALAT,IE,ITER,arrays%NAEZ, &
          arrays%NOFKS(NMESH),arrays%VOLBZ(NMESH), &
          arrays%BZKP(1,1,NMESH),arrays%VOLCUB(1,NMESH), &
          arrays%CLS,arrays%NACLS,arrays%RR, &
          arrays%EZOA,arrays%ATOM,kkr%GREFN,kkr%DGREFN, &
          params%NSYMAT,arrays%DSYMLL, &
          arrays%TMATN(:,:,ISPIN),arrays%DTDE(:,:,ISPIN), &
          arrays%NUMN0,arrays%INDN0,I1, &
          arrays%PRSC(1,1,PRSPIN), &
          EKM,NOITER, &
          params%QMRBOUND,dims%IGUESSD,dims%BCPD, &
          jij_data%NXIJ,XCCPL,jij_data%IXCP,jij_data%ZKRXIJ, &
          kkr%LLY_GRDT(IE,ISPIN),kkr%TR_ALPH(ISPIN), &
          jij_data%GMATXIJ(1,1,1,ISPIN), &
          getMySEcommunicator(my_mpi),getNumAtomRanks(my_mpi), &
          arrays%iemxd, &
          arrays%lmmaxd, arrays%naclsd, arrays%nclsd, dims%xdim, dims%ydim, dims%zdim, dims%natbld, dims%LLY, &
          jij_data%nxijd, arrays%nguessd, arrays%kpoibz, arrays%nrd, arrays%ekmd)

          TESTARRAYLOG(3, arrays%GMATN(:,:,IE,ISPIN))

          call stopTimer(mult_scattering_timer)
          call resumeTimer(single_site_timer)

        endif
      end do spinloop                          ! ISPIN = 1,NSPIN
!------------------------------------------------------------------------------
!        End of SMPID-parallel section
!------------------------------------------------------------------------------
! SPIN ==================================================================
!     END do loop over spins
! SPIN===================================================================

! =====================================================================
! Calculate Jij for the in CLSJIJ predefined atom pairs i,j
! xccpl

      if (XCCPL) then
        call jijSpinCommunication_com(my_mpi, jij_data%GMATXIJ, jij_data%DTIXIJ)

        ! calculate DTIXIJ = T_down - T_up
        call calcDeltaTupTdown(jij_data%DTIXIJ)

        JSCAL = emesh%WEZ(IE)/DBLE(jij_data%NSPIND)

        call jijLocalEnergyIntegration(my_mpi, JSCAL, jij_data%GMATXIJ, &
                                        jij_data%DTIXIJ(:,:,1), jij_data%RXIJ, jij_data%NXIJ, jij_data%IXCP, &
                                        jij_data%RXCCLS, jij_data%JXCIJINT)
      end if

! xccpl
! End of Jij calculation
! =====================================================================

      call stopEBalanceTiming(ebalance_handler, ie)

! IE ====================================================================
    endif
! IE ====================================================================

! for initial guess calculate sparse indices combining IE.KPT
    EKM = EKM + arrays%NOFKS(arrays%KMESH(IE))

  end do                   ! IE = 1,IELAST

! IE ====================================================================
!     END do loop over energies (EMPID-parallel)
! IE ====================================================================

  call stopTimer(single_site_timer)

!=======================================================================
!communicate information of 1..EMPID and 1..SMPID processors to MASTERGROUP
  call collectMSResults_com(my_mpi, arrays%GMATN, kkr%LLY_GRDT, ebalance_handler%EPROC)
!=======================================================================

! TIME
  call OUTTIME(isMasterRank(my_mpi),'Single Site took.....',getElapsedTime(single_site_timer),ITER)
  call OUTTIME(isMasterRank(my_mpi),'Mult. Scat. took.....',getElapsedTime(mult_scattering_timer),ITER)

!=======================================================================
!     output of Jij's
!=======================================================================
  if (XCCPL) then

    call jijReduceIntResults_com(my_mpi, jij_data%JXCIJINT)

    if (isInMasterGroup(my_mpi)) then
      call writeJiJs(I1,jij_data%RXIJ,jij_data%NXIJ,jij_data%IXCP, &
                     jij_data%RXCCLS,jij_data%JXCIJINT, jij_data%nxijd)
    end if
  endif

!=======================================================================
!     on the basis of new timings determine now new distribution of
!     work to 1 .. EMPID processors - all processes SYNCED
!=======================================================================
  call updateEBalance_com(ebalance_handler, my_mpi)

!=======================================================================
!     in case of IGUESS and EMPID > 1 initial guess arrays might
!     have to be adjusted to new distributions
!=======================================================================
  if ((dims%IGUESSD==1).and.(dims%EMPID>1)) then

    do ISPIN = 1,dims%NSPIND
      if(isWorkingSpinRank(my_mpi, ispin)) then

        if (dims%SMPID==1) then
          PRSPIN   = ISPIN
        else
          PRSPIN   = 1
        endif

        WRITELOG(3, *) "EPROC:     ", ebalance_handler%EPROC
        WRITELOG(3, *) "EPROC_old: ", ebalance_handler%EPROC_old

        call redistributeInitialGuess_com(my_mpi, arrays%PRSC(:,:,PRSPIN), &
             ebalance_handler%EPROC, ebalance_handler%EPROC_old, arrays%KMESH, arrays%NofKs)

      endif
    enddo

  endif  ! IGUESS == 1 .and. EMPID > 1
!=======================================================================

  TESTARRAYLOG(3, arrays%GMATN)
  TESTARRAYLOG(3, kkr%LLY_GRDT)

end subroutine

! =============================================================================
! Helper routines
! =============================================================================

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


end module
