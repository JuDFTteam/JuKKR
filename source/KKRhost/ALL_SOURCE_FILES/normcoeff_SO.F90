SUBROUTINE normcoeff_so(natom, ircut,  &
        lmmax,pns,thetas,ntcell,  &
        ifunm,ipan,lmsp,ksra,cleb,icleb,iend,drdi,  &
        irws,nspoh, mode)
! ************************************************************************
!
!     calculates the norm of the wavefunctions with full potential and
!     spin orbit coupling. this is needed for the normalization of the
!     coefficients c_Lks .
!
!     attention : the gaunt coeffients are stored in index array
!                   (see subroutine gaunt)
!
!     Added mode (can be 0/1) to determine if operator is written out using
!     host/impurity wavefunctions. In case of mode==1 the index array
!     t_imp%atomimp is used to determine which position in the host
!     corresponds to each impurity poisition (i.e. layer index)
!
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

IMPLICIT NONE

!.. Parameters ..
include 'inc.p'
INTEGER          LMMAXD
PARAMETER        (LMMAXD= (LMAXD+1)**2)
INTEGER          LMPOTD
PARAMETER        (LMPOTD= (LPOTD+1)**2)
INTEGER, PARAMETER :: NSPOD=1+KORBIT
INTEGER, PARAMETER :: NSPD=NSPIND
!..
!.. Scalar Arguments ..
INTEGER      ::   NATOM, mode, IEND,LMMAX,KSRA,IRWS(*),LMMAXSO
!..
!.. Array Arguments ..
complex (kind=dp)   PNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATOM)   ! non-sph. eigen states of single pot 
real (kind=dp) CLEB(*),THETAS(IRID,NFUND,*), &
                 DRDI(IRMD,NATYPD)                            ! derivative dr/di
INTEGER          ICLEB(NCLEB,4),IFUNM(NATYPD,LMPOTD), &
                 LMSP(NATYPD,*),IRCUT(0:IPAND,NATYPD), &
                 IPAN(NATYPD),NTCELL(*)
!..
!.. Local Scalars ..
complex (kind=dp)   CZERO,NORM
real (kind=dp) PI
INTEGER          IR,LM1,LM2,LM1P,LM2P, &
                 I1,I1SP1,I1SP2, &
                 LMSP1,LMSP2,I2SP1,I2SP2,INSRA,NSRA, &
                 NSPOH, &
                 ISIGMA, I2
! MPI stuff
INTEGER :: ierr, ihelp, i1_start, i1_end
!.. External Subroutines ..
EXTERNAL         ZGEMM
!..
!.. Intrinsic Functions ..
INTRINSIC        DATAN,aimag,SQRT
!..Local Arrays..
complex (kind=dp), ALLOCATABLE  :: DENS(:,:,:,:,:,:,:)
complex (kind=dp), ALLOCATABLE  ::   RLL_12(:,:,:), &
                                  RLL(:,:,:,:,:,:,:), &
                                  RHOD(:,:,:,:)
! MPI stuff
complex (kind=dp), ALLOCATABLE  :: work(:,:,:,:,:,:,:)
!..
!.. Data statements ..
DATA CZERO/ (0.0D0,0.0D0)/
      !..

pi=4.d0*DATAN(1.d0)
lmmaxso=2*lmmaxd
IF(t_inc%i_write>0) WRITE(1337,*) "KSRA",ksra
IF (ksra >= 1) THEN    ! previously this was .GT. which is wrong for kvrel=1
  nsra = 2
ELSE
  nsra = 1
endif
IF(t_inc%i_write>0) WRITE(1337,*) "NSRA",nsra

allocate(rll(irmd,lmmax,lmmax,nspoh,nspoh,nspoh,natom))
allocate(rll_12(irmd,lmmax,lmmax))
allocate(dens(lmmaxd,lmmaxd,nspd,nspd,nspd,nspd,natom))

rll=czero
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
      
      DO i1sp1=1,nspoh
        DO i1sp2=1,nspoh
          DO lm1 =1,lmmax
            lmsp1=(i1sp1-1)*lmmaxd+lm1
            DO lm2 =1,lmmax
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
  DO i1sp1=1,nspoh
    DO i1sp2=1,nspoh
      DO i2sp1=1,nspoh
        DO i2sp2=1,nspoh
          
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
                
                
                CALL calc_rho_ll_ss(lmmax,rll_12,  &
                    ircut(0:ipand,i2),ipan(i2),ntcell(i2),thetas,  &
                    cleb,icleb,iend,ifunm,lmsp,irws(i2), drdi(:,i2),norm)
                
                dens(lm1,lm2,i1sp1,i1sp2,i2sp1,i2sp2,i1)=  &
                    dens(lm1,lm2,i1sp1,i1sp2,i2sp1,i2sp2,i1) + norm
              END DO   !NSRA
              
            END DO         !LM2
          END DO           !LM1
          
        END DO            !I2SP2
      END DO              !I2SP1
      
    END DO         !I1SP2
  END DO           !I1SP1
  
END DO             !I1
deallocate(rll)
deallocate(rll_12)

#ifdef CPP_MPI
! finally gather DENS on master in case of MPI run
allocate(work(lmmaxd,lmmaxd,nspd,nspd,nspd,nspd,natom), stat=ierr)
IF(ierr /= 0) STOP  &
    'Error allocating work for MPI comm of DENS in normcoeff_SO'
ihelp = lmmaxd*lmmaxd*nspd*nspd*nspd*nspd*natom
CALL mpi_reduce(dens, work, ihelp, mpi_double_complex,  &
    mpi_sum, master,t_mpi_c_grid%mympi_comm_ie,ierr)
IF(ierr /= mpi_success) STOP 'Error in MPI comm of DENS in normcoeff_SO'
dens(:,:,:,:,:,:,:) = work(:,:,:,:,:,:,:)
deallocate(work, stat=ierr)
IF(ierr /= 0) STOP  &
    'Error deallocating work for MPI comm of DENS in normcoeff_SO'
#endif


IF(myrank==master) THEN ! do last part and writeout only on master
  
  WRITE(*,*) 'collect terms and writeout'
  
! calculate rho
  allocate(rhod(lmmaxso,lmmaxso,natom,4))
  IF (nspoh /= 1) THEN
    DO isigma=1,4
      DO i1=1,natom
        IF (isigma == 1) THEN
          DO i1sp1=1,nspod
            DO i1sp2=1,nspod
              DO lm1=1,lmmax
                DO lm2=1,lmmax
                  rhod((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                      dens(lm2,lm1,1,i1sp2,1,i1sp1,i1)+  &
                      dens(lm2,lm1,2,i1sp2,2,i1sp1,i1)
                END DO
              END DO
            END DO
          END DO
        ELSE IF (isigma == 2) THEN
          DO i1sp1=1,nspod
            DO i1sp2=1,nspod
              DO lm1=1,lmmax
                DO lm2=1,lmmax
                  rhod((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                      dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)+  &
                      dens(lm2,lm1,1,i1sp2,2,i1sp1,i1)
                END DO
              END DO
            END DO
          END DO
        ELSE IF (isigma == 3) THEN
          DO i1sp1=1,nspod
            DO i1sp2=1,nspod
              DO lm1=1,lmmax
                DO lm2=1,lmmax
                  rhod((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                      -(0D0,1D0)*(dens(lm2,lm1,2,i1sp2,1,i1sp1,i1)-  &
                      dens(lm2,lm1,1,i1sp2,2,i1sp1,i1))
                END DO
              END DO
            END DO
          END DO
        ELSE IF (isigma == 4) THEN
          DO i1sp1=1,nspod
            DO i1sp2=1,nspod
              DO lm1=1,lmmax
                DO lm2=1,lmmax
                  rhod((i1sp2-1)*lmmax+lm2,(i1sp1-1)*lmmax+lm1,i1,isigma)=  &
                      -dens(lm2,lm1,1,i1sp2,1,i1sp1,i1)+  &
                      dens(lm2,lm1,2,i1sp2,2,i1sp1,i1)
                END DO
              END DO
            END DO
          END DO
        endif
      END DO
    END DO
    
! write to the file
    IF (mode==0) THEN
      OPEN(UNIT=12,FILE='TBkkr_rhod.txt',FORM='formatted', action='write')
    ELSE ! mode == 1
      OPEN(UNIT=12,FILE='TBkkr_rhod_imp.txt',FORM='formatted', action='write')
    endif
    DO isigma=1,4
      DO i1=1,natom
        DO lm2=1,lmmaxso
          DO lm1=1,lmmaxso
            WRITE(12,'(2ES25.16)') rhod(lm1,lm2,i1,isigma)
          END DO
        END DO
      END DO
    END DO
    CLOSE(12)
    
  endif ! NSPOH!=1
  deallocate(dens)
  deallocate(rhod)
  
endif !(myrank==master)

END SUBROUTINE normcoeff_so

