! 12.05.2016 *************************************************************
SUBROUTINE normcoeff_so_torq(natom, ircut, lmmax,pns,ntcell,  &
    ifunm,ipan,lmsp,ksra,cleb,icleb,iend,drdi, irws,visp,nspin,vins,irmin,mode)
! ************************************************************************
!     Calculates the KKR matrix elements for the torque operator, i.e.,

!           INT dr [R^{mu}_{Ls}]^dagger T^{mu}(r) R^{mu}_{L's'}.

!     Details are in http://arxiv.org/pdf/1602.03417v1.pdf

!     This subroutine was adapted from NORMCOEFF_SO.

!                                            Guillaume Géranton, 2016
!-----------------------------------------------------------------------
#ifdef CPP_MPI
use mpi
#endif
use mod_mympi, only: myrank, master
#ifdef CPP_MPI
use mod_types, only: t_mpi_c_grid, t_inc, t_imp
#else
use mod_types, only: t_inc, t_imp
#endif

      Use mod_datatypes, Only: dp
use global_variables
      
      IMPLICIT NONE
      real(kind=dp), parameter :: eps=1.0D-12
!..
!.. Scalar Arguments ..
      INTEGER          NATOM, mode, IEND, KSRA, &
                       IRWS(*),NSPIN,IRMIN(*), I2
!< lmsize without spin degree of freedom (in contrast to lmmaxd) 
integer, intent(in) :: lmmax
!..
!.. Array Arguments ..
      complex (kind=dp)   PNS(NSPIND*LMMAX,NSPIND*LMMAX,IRMD,2,NATYPD)
      real (kind=dp) CLEB(*), &
                       DRDI(IRMD,NATYPD), &
                       VISP(IRMD,*),         &      ! spherical part of the potential
                       VINS(IRMIND:IRMD,LMPOTD,*), & ! non-spher. part of the potential
                       THETA,PHI,THETA_TMP,PHI_TMP, &
                       SQA(3)
      INTEGER          ICLEB(NCLEB,4),IFUNM(NATYPD,LMPOTD), &
                       LMSP(NATYPD,*),IRCUT(0:IPAND,NATYPD), &
                       IPAN(NATYPD),NTCELL(*)
!..
!.. Local Scalars ..
      complex (kind=dp)   CZERO,NORM
      INTEGER          LM1,LM2,LM1P,LM2P, &
                       IR,I1,I1SP1,I1SP2, &
                       LMSP1,LMSP2,ISIGMA,I2SP1,I2SP2,INSRA,NSRA
      LOGICAL          lexist
!     MPI stuff
      INTEGER :: ierr, ihelp, i1_start, i1_end
!..
!.. External Subroutines ..
      EXTERNAL         ZGEMM
!..
!.. Intrinsic Functions ..
      INTRINSIC  DATAN,aimag,DSQRT
!..
!.. Save statement ..
      SAVE             CZERO
!..
!..Local Arrays..
      complex (kind=dp), ALLOCATABLE  ::   TORQ(:,:,:,:), &
                                        RLL_12(:,:,:), &
                                        DENS(:,:,:,:,:,:,:), &
                                        RLL(:,:,:,:,:,:,:)
!     MPI stuff
      complex (kind=dp), ALLOCATABLE  :: work(:,:,:,:,:,:,:)
!..
!.. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
!..

IF(t_inc%i_write>0) WRITE(1337,*) "KSRA",ksra
IF (ksra >= 1) THEN    ! previously this was .GT. which is wrong for kvrel=1
  nsra = 2
ELSE
  nsra = 1
endif

IF(t_inc%i_write>0) THEN
  WRITE(1337,*) "NSRA",nsra
  WRITE(1337,*) "LMMAX",lmmax
  WRITE(1337,*) "LMMAXSO",lmmaxso
endif

allocate(rll(irmd,lmmax,lmmax,2,2,2,natom))
allocate(rll_12(irmd,lmmax,lmmax))
allocate(dens(lmmax,lmmax,2,2,2,2,natom))
allocate(torq(lmmaxso,lmmaxso,natom,3))

rll =czero
dens=czero


! determine MPI work division for loop over atoms
#ifdef CPP_MPI
i1_start = t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie) + 1
i1_end   = t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie) +  &
    t_mpi_c_grid%ntot_pt1(t_mpi_c_grid%myrank_ie)
#else
i1_start = 1
i1_end   = natom
#endif

! rewrite the wavefunctions in RLL arrays of 1,2*LMMAX
DO i1=i1_start, i1_end
  IF(t_inc%i_write>0) WRITE(1337,*) 'ATOM',i1, i1_start, i1_end
! use I2 as index to map for mode==1 each impurity position to the corresponding layer index of the host
  IF (mode==1) THEN
    i2 = t_imp%atomimp(i1)
  ELSE ! for mode==0 I2 and I1 are the same
    i2 = i1
  endif
  
  DO insra=1,nsra
    DO ir=1,irmd
      
      DO i1sp1=1,2
        DO i1sp2=1,2
          DO lm1 =1,lmmax
            lmsp1=(i1sp1-1)*lmmax+lm1
            DO lm2 =1,lmmax
              lmsp2=(i1sp2-1)*lmmax+lm2
              rll(ir,lm2,lm1,i1sp2,i1sp1,insra,i1)=  &
                  pns(lmsp2,lmsp1,ir,insra,i1)
            END DO      !LM1=1,LMMAX
          END DO      !LM1=1,LMMAX
        END DO      !ISP1=1,2
      END DO      !ISP1=1,2
      
    END DO      !IR
  END DO      !INSRA
  
  
! set up the array R*_L1L2 R_L3L4
  
  DO i1sp1=1,2
    DO i1sp2=1,2
      DO i2sp1=1,2
        DO i2sp2=1,2
          
          DO lm1 = 1,lmmax
            DO lm2 = 1,lmmax
              
              DO insra=1,nsra
                
                rll_12=czero
                
                DO lm1p = 1,lmmax
                  DO lm2p = 1,lmmax
                    
                    DO ir=1,irmd
                      rll_12(ir,lm1p,lm2p)=  &
                          CONJG(rll(ir,lm1p,lm1,i1sp1,i1sp2,insra,i1))*  &
                          rll(ir,lm2p,lm2,i2sp1,i2sp2,insra,i1)
                    END DO       !IR
                    
                  END DO         !LM2P
                END DO           !LM1P
                
                CALL calc_torq_ll_ss(lmmax,rll_12,  &
                    ircut(0:ipand,i2),ipan(i2),ntcell(i2),  &
                    cleb,icleb,iend,ifunm,lmsp,irws(i2),  &
                    drdi(:,i2),norm,visp,nspin,i1,vins,irmin(i2))
                
                dens(lm1,lm2,i1sp1,i1sp2,i2sp1,i2sp2,i1) =  &
                    dens(lm1,lm2,i1sp1,i1sp2,i2sp1,i2sp2,i1) + norm
              END DO   !NSRA
              
            END DO         !LM2
          END DO           !LM1
          
        END DO            !I2SP2
      END DO              !I2SP1
    END DO         !I1SP2
  END DO           !I1SP1
  
END DO             !I1

#ifdef CPP_MPI
! finally gather DENS on master in case of MPI run
allocate(work(lmmax,lmmax,2,2,2,2,natom), stat=ierr)
IF(ierr /= 0) STOP  &
    'Error allocating work for MPI comm of DENS in normcoeff_torq'
ihelp = lmmax*lmmax*2*2*2*2*natom
CALL mpi_reduce(dens, work, ihelp, mpi_double_complex,  &
    mpi_sum, master,t_mpi_c_grid%mympi_comm_ie,ierr)
IF(ierr /= mpi_success) STOP 'Error in MPI comm of DENS in normcoeff_torq'
dens(:,:,:,:,:,:,:) = work(:,:,:,:,:,:,:)
deallocate(work, stat=ierr)
IF(ierr /= 0) STOP  &
    'Error deallocating work for MPI comm of DENS in normcoeff_torq'
#endif


IF(myrank==master) THEN ! do last part and writeout only on master
  
  WRITE(*,*) 'collect terms and writeout'
  
  
! reads spin quantization axis from file
  INQUIRE(FILE='nonco_angle.dat', EXIST=lexist)
  IF(lexist) THEN
    OPEN(UNIT=11,FILE='nonco_angle.dat',FORM='FORMATTED')
    
    DO i1=1,natom
      READ(11,*) theta_tmp,phi_tmp
      IF ( i1 > 1) THEN
        IF ((abs(theta_tmp-theta)>eps ) .OR. (abs(phi_tmp-phi)>eps))  &
            STOP "It seems you want to compute the torque  &
            in a non-colinear system. This is not implemented yet."
      endif
      theta = theta_tmp
      phi = phi_tmp
    END DO
    CLOSE(11)
  else
    theta = 0.0_dp
    phi = 0.0_dp
  endif
  
  sqa(1) = SIN(theta)*COS(phi)
  sqa(2) = SIN(theta)*SIN(phi)
  sqa(3) = COS(theta)
  
  torq=czero
  
  DO isigma=1,3  !ISIGMA == 1 --> T_x
!ISIGMA == 2 --> T_y
!ISIGMA == 3 --> T_z
    
    WRITE(6,*) "ISIGMA",isigma
    DO i1=1,natom
      
! TEMPORARY IMPLEMENTATION OF THE TORQUE OPERATOR FOR B // z
      IF (isigma==1) THEN  !T_x
        
        DO i1sp1=1,2
          DO i1sp2=1,2
            DO lm1=1,lmmax
              DO lm2=1,lmmax
                torq((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                    (0D0,1D0)*(dens(lm2,lm1,1,i1sp2,2,i1sp1,i1)-  &
                    dens(lm2,lm1,2,i1sp2,1,i1sp1,i1))*sqa(3)  &
                    -(-dens(lm2,lm1,1,i1sp2,1,i1sp1,i1)+  &
                    dens(lm2,lm1,2,i1sp2,2,i1sp1,i1))*sqa(2)
              END DO !LM2
            END DO !LM1
          END DO !I1SP2
        END DO !I1SP1
        
      ELSE IF (isigma==2) THEN !T_y
        
        DO i1sp1=1,2
          DO i1sp2=1,2
            DO lm1=1,lmmax
              DO lm2=1,lmmax
                torq((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                    (-dens(lm2,lm1,1,i1sp2,1,i1sp1,i1)+  &
                    dens(lm2,lm1,2,i1sp2,2,i1sp1,i1))*sqa(1)  &
                    -(dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+  &
                    dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))*sqa(3)
              END DO !LM2
            END DO !LM1
          END DO !I1SP2
        END DO !I1SP1
        
      ELSE IF (isigma==3) THEN !T_z
        
        DO i1sp1=1,2
          DO i1sp2=1,2
            DO lm1=1,lmmax
              DO lm2=1,lmmax
                torq((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                    (dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+  &
                    dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))*sqa(2)  &
                    -(0D0,1D0)*(-dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+  &
                    dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))*sqa(1)
              END DO !LM2
            END DO !LM1
          END DO !I1SP2
        END DO !I1SP1
        
      endif ! (ISIGMA=1,2,3)
      
    END DO  !I1
  END DO    !ISIGMA
  
! writeout
  IF(mode==0) THEN
    OPEN(UNIT=12, FILE='TBkkr_torq.txt', FORM='formatted',action='write')
  ELSE ! mode==1
    OPEN(UNIT=12, FILE='TBkkr_torq_imp.txt', FORM='formatted',action='write')
  endif
  DO isigma=1,3
    DO i1=1,natom
      DO lm2=1,lmmaxso
        DO lm1=1,lmmaxso
          WRITE(12,'(2ES25.16)') torq(lm1,lm2,i1,isigma)
        END DO
      END DO
    END DO
  END DO
  CLOSE(12)
  
endif !(myrank==master)

deallocate(rll)
deallocate(dens)
deallocate(rll_12)
deallocate(torq)

END SUBROUTINE normcoeff_so_torq

