#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"

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
                      ebalance_handler, my_mpi, arrays, jij_data, ldau_data)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use TimerMpi_mod
  use EBalanceHandler_mod

  use TEST_lcutoff_mod !TODO: remove

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod

  use GauntCoefficients_mod
  use BasisAtom_mod
  use JijData_mod
  use LDAUData_mod
  use EnergyMesh_mod
  use KKRnanoParallel_mod
  use EBalanceHandler_mod

  use main2_aux_mod,    only: substractReferenceTmatrix, printEnergyPoint, calcDeltatuptdown
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
                      arrays%TREFLL(1,1,RF),arrays%DTREFLL(1,1,RF), dims%LLY)
          end do

          TESTARRAYLOG(3, arrays%TREFLL)
          TESTARRAYLOG(3, arrays%DTREFLL)

          call GREF_com(emesh%EZ(IE),params%ALAT,gaunts%IEND,params%NCLS,arrays%NAEZ, &
                        gaunts%CLEB,arrays%RCLS,arrays%ATOM,arrays%CLS,gaunts%ICLEB, &
                        gaunts%LOFLM,arrays%NACLS, &
                        arrays%REFPOT, &
                        arrays%TREFLL,arrays%DTREFLL,arrays%GREFN,arrays%DGREFN, &
                        arrays%LLY_G0TR(:,IE), &
                        getMyAtomRank(my_mpi),getMySEcommunicator(my_mpi),&
                        getNumAtomRanks(my_mpi), &
                        arrays%lmaxd, arrays%naclsd, gaunts%ncleb, arrays%nrefd, arrays%nclsd, &
                        dims%LLY)

          TESTARRAYLOG(3, arrays%GREFN)
          TESTARRAYLOG(3, arrays%DGREFN)

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

              call CALCTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, gaunts, arrays%TMATN, arrays%TR_ALPH, ldau_data)

              jij_data%DTIXIJ(:,:,ISPIN) = arrays%TMATN(:,:,ISPIN)  ! save t-matrix for Jij-calc.

              if(dims%LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula
                call CALCDTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, params%NSRA, gaunts, arrays%DTDE, arrays%TR_ALPH, ldau_data)
              end if

              RF = arrays%REFPOT(I1)
              call substractReferenceTmatrix(arrays%TMATN(:,:,ISPIN), arrays%TREFLL(:,:,RF), arrays%LMMAXD)
              ! do the same for derivative of T-matrix
              call substractReferenceTmatrix(arrays%DTDE(:,:,ISPIN), arrays%DTREFLL(:,:,RF), arrays%LMMAXD)

              ! TMATN now contains Delta t = t - t_ref !!!
              ! DTDE now contains Delta dt !!!

              ! renormalize TR_ALPH
              arrays%TR_ALPH(ISPIN) = arrays%TR_ALPH(ISPIN) - arrays%LLY_G0TR(arrays%CLS(I1), IE)

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
              arrays%EZOA,arrays%ATOM,arrays%GREFN,arrays%DGREFN, &
              params%NSYMAT,arrays%DSYMLL, &
              arrays%TMATN(:,:,ISPIN),arrays%DTDE(:,:,ISPIN), &
              arrays%NUMN0,arrays%INDN0,I1, &
              arrays%PRSC(1,1,PRSPIN), &
              EKM,NOITER, &
              params%QMRBOUND,dims%IGUESSD,dims%BCPD, &
              jij_data%NXIJ,XCCPL,jij_data%IXCP,jij_data%ZKRXIJ, &
              arrays%LLY_GRDT(IE,ISPIN),arrays%TR_ALPH(ISPIN), &
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
      call collectMSResults_com(my_mpi, arrays%GMATN, arrays%LLY_GRDT, ebalance_handler%EPROC)
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
          call writeJiJs(I1,jij_data%RXIJ,jij_data%NXIJ,jij_data%IXCP,jij_data%RXCCLS,jij_data%JXCIJINT, jij_data%nxijd)
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
      TESTARRAYLOG(3, arrays%LLY_GRDT)

end subroutine
