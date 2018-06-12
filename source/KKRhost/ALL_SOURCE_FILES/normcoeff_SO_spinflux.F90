! 05.10.10 ***************************************************************
SUBROUTINE normcoeff_so_spinflux(natom, ircut, lmmax,pns,  &
    ksra,drdi, mode)
! ************************************************************************
!     Calculates the KKR matrix elements for the spin flux operator, i.e.,

!           INT dr [R^{mu}_{Ls}]^dagger Q_s R^{mu}_{L's'}.

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

!..
!.. Scalar Arguments ..
      INTEGER          NATOM, mode, LMMAX,KSRA
!..
!.. Array Arguments ..
      complex (kind=dp)   PNS(NSPIND*LMMAXD,NSPIND*LMMAXD,IRMD,2,NATOM)
      real (kind=dp) DRDI(IRMD,NATYPD)
      INTEGER          IRCUT(0:IPAND,NATYPD)
!..
!.. Local Scalars ..
      complex (kind=dp)   CZERO
      INTEGER          LM1,LM2,LM1P, &
                       IR,I1,I1SP1,I1SP2, &
                       LMSP1,LMSP2,ISIGMA,I2SP1,I2SP2,INSRA,NSRA
      INTEGER          I2
      complex (kind=dp) :: DELTA1, DELTA2
!     MPI stuff
      INTEGER :: ierr, ihelp, i1_start, i1_end
!..
!.. External Subroutines ..
      EXTERNAL         ZGEMM
!..
!.. Intrinsic Functions ..
      INTRINSIC DATAN,aimag,DSQRT
!..
!.. Save statement ..
      SAVE             CZERO
!..
!..Local Arrays..
      complex (kind=dp), ALLOCATABLE  ::   SPINFLUX(:,:,:,:), &
                                        RLL_12(:), &
                                        DENS(:,:,:,:,:,:,:), &
                                        RLL(:,:,:,:,:,:,:)
!     MPI stuff
      complex (kind=dp), ALLOCATABLE  :: work(:,:,:,:,:,:,:)
!..
!.. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/

IF(t_inc%i_write>0) WRITE(1337,*) "KSRA",ksra
IF (ksra >= 1) THEN    ! previously this was .GT. which is wrong for kvrel=1
  nsra = 2
ELSE
  nsra = 1
endif

IF (t_inc%i_write>0) THEN
  WRITE(1337,*) "NSRA",nsra
  WRITE(1337,*) "LMMAX",lmmax
  WRITE(1337,*) "LMMAXD",lmmaxd
  WRITE(1337,*) "LMMAXSO",lmmaxso
endif


allocate(rll(irmd,lmmax,lmmax,2,2,2,natom))
allocate(rll_12(lmmax))
allocate(dens(lmmaxd,lmmaxd,2,2,2,2,natom))
allocate(spinflux(lmmaxso,lmmaxso,natom,3))

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

! rewrite the wavefunctions in RLL arrays of 1,2*LMMAXD
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
          DO lm1 =1,lmmaxd
            lmsp1=(i1sp1-1)*lmmaxd+lm1
            DO lm2 =1,lmmaxd
              lmsp2=(i1sp2-1)*lmmaxd+lm2
              rll(ir,lm2,lm1,i1sp2,i1sp1,insra,i1)=  &
                  pns(lmsp2,lmsp1,ir,insra,i1)
            END DO      !LM1=1,LMMAXD
          END DO      !LM1=1,LMMAXD
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
                  
                  delta1=(rll(ircut(1,i2),lm1p,lm2,i2sp1,  &
                      i2sp2,insra,i1)-rll(ircut(1,i2)-1,  &
                      lm1p,lm2,i2sp1,i2sp2,insra,i1)) /drdi(ircut(1,i2),i2)
                  delta2=(rll(ircut(1,i2),lm1p,lm1,i1sp1,  &
                      i1sp2,insra,i1)-rll(ircut(1,i2)-1,  &
                      lm1p,lm1,i1sp1,i1sp2,insra,i1)) /drdi(ircut(1,i2),i2)
                  
                  rll_12(lm1p)=  &
                      CONJG(rll(ircut(1,i2)-1,lm1p,lm1,i1sp1,i1sp2,  &
                      insra,i1))*delta1-  &
                      rll(ircut(1,i2)-1,lm1p,lm2,i2sp1,i2sp2,insra,i1)  &
                      *CONJG(delta2)
                  
                END DO!LM1P
                
                DO lm1p = 1,lmmax
                  dens(lm1,lm2,i1sp1,i1sp2,i2sp1,i2sp2,i1) =  &
                      dens(lm1,lm2,i1sp1,i1sp2,i2sp1,i2sp2,i1) + rll_12(lm1p)
                END DO!LM1P
                
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
allocate(work(lmmaxd,lmmaxd,2,2,2,2,natom), stat=ierr)
IF(ierr /= 0) STOP  &
    'Error allocating work for MPI comm of DENS in normcoeff_spinf'
ihelp = lmmaxd*lmmaxd*2*2*2*2*natom
CALL mpi_reduce(dens, work, ihelp, mpi_double_complex,  &
    mpi_sum, master,t_mpi_c_grid%mympi_comm_ie,ierr)
IF(ierr /= mpi_success) STOP  &
    'Error in MPI comm of DENS in normcoeff_spinflux'
dens(:,:,:,:,:,:,:) = work(:,:,:,:,:,:,:)
deallocate(work, stat=ierr)
IF(ierr /= 0) STOP  &
    'Error deallocating work for MPI comm of DENS in normcoeff_spin'
#endif


IF(myrank==master) THEN ! do last part and writeout only on master
  
  WRITE(*,*) 'collect terms and writeout'
  
  spinflux=czero
  
  DO isigma=1,3  !ISIGMA == 1 --> Q_x
!ISIGMA == 2 --> Q_y
!ISIGMA == 3 --> Q_z
    
    WRITE(6,*) "ISIGMA",isigma
    DO i1=1,natom
      
      IF (isigma==1) THEN  !Q_x
        
        DO i1sp1=1,2
          DO i1sp2=1,2
            DO lm1=1,lmmax
              DO lm2=1,lmmax
                spinflux((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                    -(0D0,1D0)*(dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+  &
                    dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))/2
              END DO !LM2
            END DO !LM1
          END DO !I1SP2
        END DO !I1SP1
        
      ELSE IF (isigma==2) THEN !Q_y
        
        DO i1sp1=1,2
          DO i1sp2=1,2
            DO lm1=1,lmmax
              DO lm2=1,lmmax
                spinflux((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                    -(0D0,1D0)*(-1)*(0D0,1D0)*(dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)-  &
                    dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))/2
              END DO !LM2
            END DO !LM1
          END DO !I1SP2
        END DO !I1SP1
        
      ELSE IF (isigma==3) THEN !Q_z
        
        DO i1sp1=1,2
          DO i1sp2=1,2
            DO lm1=1,lmmax
              DO lm2=1,lmmax
                spinflux((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                    (0D0,1D0)*(dens(lm2,lm1,1,i1sp2,1,i1sp1,i1)-  &
                    dens(lm2,lm1,2,i1sp2,2,i1sp1,i1))/2
              END DO !LM2
            END DO !LM1
          END DO !I1SP2
        END DO !I1SP1
        
      endif
      
    END DO  !I1
  END DO    !ISIGMA
  
! write to file
  IF(mode==0) THEN
    OPEN(UNIT=12, FILE='TBkkr_spinflux.txt', FORM='formatted',action='write')
  ELSE ! mode==1
    OPEN(UNIT=12, FILE='TBkkr_spinflux_imp.txt',  &
        FORM='formatted',action='write')
  endif
  DO isigma=1,3
    DO i1=1,natom
      DO lm2=1,lmmaxso
        DO lm1=1,lmmaxso
!minus sign to get the spin flux into the sphere :
          WRITE(12,'(2ES25.16)') -spinflux(lm1,lm2,i1,isigma)
        END DO
      END DO
    END DO
  END DO
  CLOSE(12)
  
endif !(myrank==master)

deallocate(rll)
deallocate(dens)
deallocate(rll_12)
deallocate(spinflux)

END SUBROUTINE normcoeff_so_spinflux




