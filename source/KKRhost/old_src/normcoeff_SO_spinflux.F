c 05.10.10 ***************************************************************
      SUBROUTINE NORMCOEFF_SO_SPINFLUX(NATOM, IRCUT,
     +                   LMMAX,PNS,
     +                   KSRA,DRDI, mode)
c ************************************************************************
c     Calculates the KKR matrix elements for the spin flux operator, i.e.,
c     
c           INT dr [R^{mu}_{Ls}]^dagger Q_s R^{mu}_{L's'}.
c
c     Details are in http://arxiv.org/pdf/1602.03417v1.pdf
c     
c     This subroutine was adapted from NORMCOEFF_SO.
c
c                                            Guillaume Géranton, 2016
c-----------------------------------------------------------------------
#ifdef CPP_MPI
      use mpi
#endif
      use mod_mympi, only: myrank, master
#ifdef CPP_MPI
      use mod_types, only: t_mpi_c_grid, t_inc, t_imp
#else
      use mod_types, only: t_inc, t_imp
#endif

      IMPLICIT NONE

C     .. Parameters ..
      include 'inc.p'
      INTEGER          LMMAXD
      PARAMETER        (LMMAXD= (LMAXD+1)**2)
      INTEGER          LMPOTD
      PARAMETER        (LMPOTD= (LPOTD+1)**2)
      INTEGER          IRMIND,IRLMD
      PARAMETER        (IRMIND=IRMD-IRNSD,IRLMD= (IRNSD+1)*LMMAXD)
      INTEGER, PARAMETER :: NSPD=NSPIND
C     ..
C     .. Scalar Arguments ..
      INTEGER          NATOM, mode, LMMAX,KSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   PNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATOM)
      DOUBLE PRECISION DRDI(IRMD,NATYPD)
      INTEGER          IRCUT(0:IPAND,NATYPD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   CZERO
      INTEGER          LM1,LM2,LM1P,
     +                 IR,I1,I1SP1,I1SP2,
     +                 LMSP1,LMSP2,ISIGMA,I2SP1,I2SP2,INSRA,NSRA
      INTEGER          LMMAXSO, I2
      DOUBLE COMPLEX :: DELTA1, DELTA2
!     MPI stuff
      INTEGER :: ierr, ihelp, i1_start, i1_end
C     ..
C     .. External Subroutines ..
      EXTERNAL         ZGEMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        DATAN,DIMAG,DSQRT
C     ..
C     .. Save statement ..
      SAVE             CZERO
C     ..
C     ..Local Arrays..
      DOUBLE COMPLEX, ALLOCATABLE  ::   SPINFLUX(:,:,:,:),
     +                                  RLL_12(:),
     +                                  DENS(:,:,:,:,:,:,:),
     +                                  RLL(:,:,:,:,:,:,:)
!     MPI stuff
      DOUBLE COMPLEX, ALLOCATABLE  :: work(:,:,:,:,:,:,:)
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..
c
      LMMAXSO=2*LMMAXD

      if(t_inc%i_write>0) WRITE(1337,*) "KSRA",KSRA
      IF (KSRA.GE.1) THEN    ! previously this was .GT. which is wrong for kvrel=1
         NSRA = 2
      ELSE
         NSRA = 1
      END IF

      if (t_inc%i_write>0) then
        WRITE(1337,*) "NSRA",NSRA
        WRITE(1337,*) "LMMAX",LMMAX
        WRITE(1337,*) "LMMAXD",LMMAXD
        WRITE(1337,*) "LMMAXSO",LMMAXSO
      end if


      ALLOCATE(RLL(IRMD,LMMAX,LMMAX,2,2,2,NATOM))
      ALLOCATE(RLL_12(LMMAX))
      ALLOCATE(DENS(LMMAXD,LMMAXD,2,2,2,2,NATOM))
      ALLOCATE(SPINFLUX(LMMAXSO,LMMAXSO,NATOM,3))

      RLL =CZERO
      DENS=CZERO


      ! determine MPI work division for loop over atoms
#ifdef CPP_MPI
      i1_start = t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie) + 1
      i1_end   = t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie) +
     &           t_mpi_c_grid%ntot_pT1(t_mpi_c_grid%myrank_ie)
#else
      i1_start = 1
      i1_end   = NATOM
#endif

      ! rewrite the wavefunctions in RLL arrays of 1,2*LMMAXD
      DO I1=i1_start, i1_end
        if(t_inc%i_write>0) WRITE(1337,*) 'ATOM',I1, i1_start, i1_end

        ! use I2 as index to map for mode==1 each impurity position to the corresponding layer index of the host
        if (mode==1) then
           I2 = t_imp%ATOMIMP(I1)
        else ! for mode==0 I2 and I1 are the same
           I2 = I1
        end if

        DO INSRA=1,NSRA
          DO IR=1,IRMD

            DO I1SP1=1,2
              DO I1SP2=1,2
                DO LM1 =1,LMMAXD
                  LMSP1=(I1SP1-1)*LMMAXD+LM1
                  DO LM2 =1,LMMAXD
                    LMSP2=(I1SP2-1)*LMMAXD+LM2
                    RLL(IR,LM2,LM1,I1SP2,I1SP1,INSRA,I1)=
     +                        PNS(LMSP2,LMSP1,IR,INSRA,I1)
                  END DO      !LM1=1,LMMAXD
                END DO      !LM1=1,LMMAXD
              END DO      !ISP1=1,2
            END DO      !ISP1=1,2

          END DO      !IR
        END DO      !INSRA


        ! set up the array R*_L1L2 R_L3L4 
        DO I1SP1=1,2
          DO I1SP2=1,2
            DO I2SP1=1,2
              DO I2SP2=1,2

                DO LM1 = 1,LMMAX
                  DO LM2 = 1,LMMAX

                    DO INSRA=1,NSRA

                      RLL_12=CZERO

                      DO LM1P = 1,LMMAX

                        DELTA1=(RLL(IRCUT(1,I2),LM1P,LM2,I2SP1,
     +                          I2SP2,INSRA,I1)-RLL(IRCUT(1,I2)-1,
     +                          LM1P,LM2,I2SP1,I2SP2,INSRA,I1))
     +                          /DRDI(IRCUT(1,I2),I2)
                        DELTA2=(RLL(IRCUT(1,I2),LM1P,LM1,I1SP1,
     +                          I1SP2,INSRA,I1)-RLL(IRCUT(1,I2)-1,
     +                          LM1P,LM1,I1SP1,I1SP2,INSRA,I1))
     +                          /DRDI(IRCUT(1,I2),I2)

                        RLL_12(LM1P)=
     +                  DCONJG(RLL(IRCUT(1,I2)-1,LM1P,LM1,I1SP1,I1SP2,
     +                  INSRA,I1))*DELTA1-
     +                  RLL(IRCUT(1,I2)-1,LM1P,LM2,I2SP1,I2SP2,INSRA,I1)
     +                  *DCONJG(DELTA2)

                      END DO!LM1P

                      DO LM1P = 1,LMMAX
                        DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1) =
     +                  DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1)
     +                  + RLL_12(LM1P)
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
      allocate(work(LMMAXD,LMMAXD,2,2,2,2,NATOM), 
     &         stat=ierr)
      if(ierr.ne.0) stop 
     &  'Error allocating work for MPI comm of DENS in normcoeff_spinf'
      ihelp = LMMAXD*LMMAXD*2*2*2*2*NATOM
      call MPI_REDUCE(DENS, work, ihelp, MPI_DOUBLE_COMPLEX,
     &                MPI_SUM, master,t_mpi_c_grid%myMPI_comm_ie,ierr)
      if(ierr.ne.MPI_SUCCESS) stop 
     &  'Error in MPI comm of DENS in normcoeff_spinflux'
      DENS(:,:,:,:,:,:,:) = work(:,:,:,:,:,:,:)
      deallocate(work, stat=ierr)
      if(ierr.ne.0) stop 
     &  'Error deallocating work for MPI comm of DENS in normcoeff_spin'
#endif


      if(myrank==master) then ! do last part and writeout only on master

        write(*,*) 'collect terms and writeout'

        SPINFLUX=CZERO
        
        DO ISIGMA=1,3  !ISIGMA == 1 --> Q_x
                       !ISIGMA == 2 --> Q_y
                       !ISIGMA == 3 --> Q_z
        
          WRITE(6,*) "ISIGMA",ISIGMA
          DO I1=1,NATOM
        
            IF (ISIGMA==1) THEN  !Q_x 
        
              DO I1SP1=1,2
                DO I1SP2=1,2
                  DO LM1=1,LMMAX
                    DO LM2=1,LMMAX
            SPINFLUX((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +               -(0d0,1d0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)+
     +               DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1))/2
                    END DO !LM2
                  END DO !LM1
                END DO !I1SP2
              END DO !I1SP1
        
            ELSE IF (ISIGMA==2) THEN !Q_y
        
              DO I1SP1=1,2
                DO I1SP2=1,2
                  DO LM1=1,LMMAX
                    DO LM2=1,LMMAX
            SPINFLUX((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +      -(0d0,1d0)*(-1)*(0d0,1d0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)-
     +      DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1))/2
                    END DO !LM2
                  END DO !LM1
                END DO !I1SP2
              END DO !I1SP1
        
            ELSE IF (ISIGMA==3) THEN !Q_z
        
              DO I1SP1=1,2
                DO I1SP2=1,2
                  DO LM1=1,LMMAX
                    DO LM2=1,LMMAX
            SPINFLUX((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +               (0d0,1d0)*(DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)-
     +               DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1))/2
                    END DO !LM2
                  END DO !LM1
                END DO !I1SP2
              END DO !I1SP1
        
            END IF
        
          END DO  !I1
        END DO    !ISIGMA
        
        ! write to file
        if(mode==0) then
          open(unit=12, file='TBkkr_spinflux.txt',
     +       form='formatted',action='write')
        else ! mode==1
          open(unit=12, file='TBkkr_spinflux_imp.txt',
     +       form='formatted',action='write')
        end if
        DO ISIGMA=1,3
          DO I1=1,NATOM
            DO LM2=1,LMMAXSO
              DO LM1=1,LMMAXSO
                   !minus sign to get the spin flux into the sphere :
                   WRITE(12,'(2ES25.16)') -SPINFLUX(LM1,LM2,I1,ISIGMA)
              END DO
            END DO
          END DO
        END DO
        close(12)

      end if !(myrank==master)

      DEALLOCATE(RLL)
      DEALLOCATE(DENS)
      DEALLOCATE(RLL_12)
      DEALLOCATE(SPINFLUX)

      END




