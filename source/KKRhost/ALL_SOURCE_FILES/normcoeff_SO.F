      SUBROUTINE NORMCOEFF_SO(IRCUT,
     +                   LMMAX,PNS,THETAS,NTCELL,
     +                   IFUNM,IPAN,LMSP,KSRA,CLEB,ICLEB,IEND,DRDI,
     +                   IRWS,NSPOH)
c ************************************************************************
c
c     calculates the norm of the wavefunctions with full potential and
c     spin orbit coupling. this is needed for the normalization of the
c     coefficients c_Lks .
c
c     attention : the gaunt coeffients are stored in index array
c                   (see subroutine gaunt)
c
c-----------------------------------------------------------------------
C     .. Parameters ..
#ifdef CPP_MPI
      use mpi
#endif
      use mod_mympi, only: myrank, master
#ifdef CPP_MPI
      use mod_types, only: t_mpi_c_grid
#endif
      
      IMPLICIT NONE

      include 'inc.p'
      INTEGER          LMMAXD
      PARAMETER        (LMMAXD= (LMAXD+1)**2)
      INTEGER          LMPOTD
      PARAMETER        (LMPOTD= (LPOTD+1)**2)
      INTEGER, PARAMETER :: NSPOD=1+KORBIT
      INTEGER, PARAMETER :: NSPD=NSPIND
C     ..
C     .. Scalar Arguments ..
      INTEGER      ::   IEND,LMMAX,KSRA,IRWS(*),LMMAXSO
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   PNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2,NATYPD)   ! non-sph. eigen states of single pot 
      DOUBLE PRECISION CLEB(*),THETAS(IRID,NFUND,*),
     +                 DRDI(IRMD,NATYPD)                            ! derivative dr/di
      INTEGER          ICLEB(NCLEB,4),IFUNM(NATYPD,LMPOTD),
     +                 LMSP(NATYPD,*),IRCUT(0:IPAND,NATYPD),
     +                 IPAN(NATYPD),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   CZERO,NORM
      DOUBLE PRECISION PI
      INTEGER          IR,LM1,LM2,LM1P,LM2P,
     +                 I1,I1SP1,I1SP2,I1SP,
     +                 LMSP1,LMSP2,I2SP1,I2SP2,INSRA,NSRA,
     +                 NSPOH,
     +                 LMDOS,ISIGMA
      LOGICAL          TEST
!     MPI stuff
      INTEGER :: ierr, ihelp, i1_start, i1_end
C     .. External Subroutines ..
      EXTERNAL         ZGEMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC        DATAN,DIMAG,DSQRT
C     ..Local Arrays..
      DOUBLE COMPLEX, ALLOCATABLE  :: DENS(:,:,:,:,:,:,:)
      DOUBLE COMPLEX, ALLOCATABLE  ::   RLL_12(:,:,:),
     +                                  RLL(:,:,:,:,:,:,:),
     +                                  RHOD(:,:,:,:)
!     MPI stuff
      DOUBLE COMPLEX, ALLOCATABLE  :: work(:,:,:,:,:,:,:)
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..
c
      PI=4.d0*DATAN(1.d0)
      LMMAXSO=2*LMMAXD
      if(myrank==master) WRITE(6,*) "KSRA",KSRA
      IF (KSRA.GE.1) THEN    ! previously this was .GT. which is wrong for kvrel=1
         NSRA = 2
      ELSE
         NSRA = 1
      END IF
      if(myrank==master) WRITE(6,*) "NSRA",NSRA

      ALLOCATE(RLL(IRMD,LMMAX,LMMAX,NSPOH,NSPOH,NSPOH,NATYPD))
      ALLOCATE(RLL_12(IRMD,LMMAX,LMMAX))
      ALLOCATE(DENS(LMMAXD,LMMAXD,NSPD,NSPD,NSPD,NSPD,NATYPD))

      RLL=CZERO
      DENS=CZERO


      ! determine MPI work division for loop over atoms
#ifdef CPP_MPI
      i1_start = t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie) + 1
      i1_end   = t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie) +
     &           t_mpi_c_grid%ntot_pT1(t_mpi_c_grid%myrank_ie)
#else
      i1_start = 1
      i1_end   = NATYPD
#endif

      ! rewrite the wavefunctions in RLL arrays of 1,2*LMMAXD
      DO I1=i1_start, i1_end
       WRITE(6,*) 'ATOM',I1, i1_start, i1_end, myrank

        DO INSRA=1,NSRA
          DO IR=1,IRMD

            DO I1SP1=1,NSPOH
              DO I1SP2=1,NSPOH
                DO LM1 =1,LMMAX
                  LMSP1=(I1SP1-1)*LMMAXD+LM1
                  DO LM2 =1,LMMAX
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
        DO I1SP1=1,NSPOH
          DO I1SP2=1,NSPOH
            DO I2SP1=1,NSPOH
              DO I2SP2=1,NSPOH

                DO LM1 = 1,LMMAX
                  DO LM2 = 1,LMMAX
         
                    DO INSRA=1,NSRA
 
                    RLL_12=CZERO
         
                      DO LM1P = 1,LMMAX
                        DO LM2P = 1,LMMAX
           
                          DO IR=1,IRMD
           
                            RLL_12(IR,LM1P,LM2P)=
     +                    DCONJG(RLL(IR,LM1P,LM1,I1SP1,I1SP2,INSRA,I1))*
     +                            RLL(IR,LM2P,LM2,I2SP1,I2SP2,INSRA,I1)
                          END DO       !IR
           
                        END DO         !LM2P
                      END DO           !LM1P


                      CALL CALC_RHO_LL_SS(LMMAX,RLL_12,
     +                     IRCUT(0:IPAND,I1),IPAN(I1),NTCELL(I1),THETAS,
     +                     CLEB,ICLEB,IEND,IFUNM,LMSP,IRWS(I1),
     +                     DRDI(:,I1),NORM)

                      DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1)=
     +                DENS(LM1,LM2,I1SP1,I1SP2,I2SP1,I2SP2,I1) + NORM
                    END DO   !NSRA

                  END DO         !LM2
                END DO           !LM1

              END DO            !I2SP2
            END DO              !I2SP1

          END DO         !I1SP2
        END DO           !I1SP1 

      END DO             !I1
      DEALLOCATE(RLL)
      DEALLOCATE(RLL_12) 

#ifdef CPP_MPI
      ! finally gather DENS on master in case of MPI run
      allocate(work(LMMAXD,LMMAXD,NSPD,NSPD,NSPD,NSPD,NATYPD), 
     &         stat=ierr)
      if(ierr.ne.0) stop 
     &  'Error allocating work for MPI comm of PNS_SO_ALL in main1a'
      ihelp = LMMAXD*LMMAXD*NSPD*NSPD*NSPD*NSPD*NATYPD
      write(*,*) 'comm dens', myrank, ihelp, shape(dens),shape(work)
      call MPI_REDUCE(DENS, work, ihelp, MPI_DOUBLE_COMPLEX,
!     &                MPI_SUM, master, MPI_COMM_WORLD,ierr)
     &                MPI_SUM, master,t_mpi_c_grid%myMPI_comm_ie,ierr)
      if(ierr.ne.MPI_SUCCESS) stop 
     &  'Error in MPI comm of DENS in normcoeff_SO'
      DENS(:,:,:,:,:,:,:) = work(:,:,:,:,:,:,:)
      deallocate(work, stat=ierr)
      if(ierr.ne.0) stop 
     &  'Error deallocating work for MPI comm of DENS in normcoeff_SO'
#endif


      if(myrank==master) then ! do last part and writeout only on master

        write(*,*) 'collect terms and writeout'

        ! calculate rho
        ALLOCATE(RHOD(LMMAXSO,LMMAXSO,NATYPD,4))
        IF (NSPOH.NE.1) THEN
          DO ISIGMA=1,4
           DO I1=1,NATYPD
            IF (ISIGMA.EQ.1) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMMAX
                DO LM2=1,LMMAX
          RHOD((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +           DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)+ 
     +           DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1) 
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ELSEIF (ISIGMA.EQ.2) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMMAX
                DO LM2=1,LMMAX
          RHOD((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +             DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)+ 
     +             DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)
                ENDDO
               ENDDO
              ENDDO
             ENDDO
            ELSEIF (ISIGMA.EQ.3) THEN
             DO I1SP1=1,NSPOD
              DO I1SP2=1,NSPOD
               DO LM1=1,LMMAX
                DO LM2=1,LMMAX
          RHOD((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +          -(0D0,1D0)*(DENS(LM2,LM1,2,I1SP2,1,I1SP1,I1)-
     +                     DENS(LM2,LM1,1,I1SP2,2,I1SP1,I1)) 
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ELSEIF (ISIGMA.EQ.4) THEN
              DO I1SP1=1,NSPOD
               DO I1SP2=1,NSPOD
                DO LM1=1,LMMAX
                 DO LM2=1,LMMAX
          RHOD((I1SP2-1)*LMMAX+LM2,(I1SP1-1)*LMMAX+LM1,I1,ISIGMA)=
     +           -DENS(LM2,LM1,1,I1SP2,1,I1SP1,I1)+
     +            DENS(LM2,LM1,2,I1SP2,2,I1SP1,I1) 
                 ENDDO
                ENDDO
               ENDDO
              ENDDO
             ENDIF
            ENDDO
           ENDDO

          ! write to the file
          OPEN(UNIT=12,FILE='TBkkr_rhod.txt',FORM='formatted',
     +       ACTION='write')
          DO ISIGMA=1,4
           DO I1=1,NATYPD
            DO LM2=1,LMMAXSO
             DO LM1=1,LMMAXSO
              WRITE(12,'(2ES25.16)') RHOD(LM1,LM2,I1,ISIGMA)
             ENDDO
            ENDDO
           ENDDO
          ENDDO        
          CLOSE(12)

        ENDIF ! NSPOH!=1
        DEALLOCATE(DENS)
        DEALLOCATE(RHOD)

      end if !(myrank==master)

      END SUBROUTINE NORMCOEFF_SO

