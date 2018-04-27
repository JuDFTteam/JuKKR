       SUBROUTINE TMATIMP_NEWSOLVER(NSPIN,LMAX,R,Z,IELAST,E,KSRA,
     &                     IRWS,IPAN,IRCUT,IRMIN,C,CLEB,ICLEB,IEND,
     &                     NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,VINS,VM2Z,
     &                     NATOMIMP,RCLSIMP,ATOMIMP,IHOST,HOSTIMP,RIMP,
     &                     ZIMP,IRWSIMP,IPANIMP,IRCUTIMP,IRMINIMP,
     &                     VINSIMP,VM2ZIMP,DTMTRX)
!-----------------------------------------------------------------
! calculate and write down impurity tmatrix and delta matrix
! first calculate t-matrix for the host corresponding to imp. cluster
! N. H. Long, Juelich, 05.2013
! adapted to new routines (mainly changed interfaces) to work in KKRcode
! also added MPI parallelization
! P.R., Juelich, 09.2017
!-----------------------------------------------------------------
#ifdef CPP_MPI
       use mpi
       use mod_mympi, only: myrank, master, nranks,
     &                      distribute_linear_on_tasks
#else
       use mod_mympi, only: myrank, master, nranks
#endif
       use mod_types, only: t_inc, t_imp
       use mod_create_newmesh
       use mod_version_info
       use mod_wunfiles, only: t_params
       IMPLICIT NONE
       include 'inc.p'
       INTEGER NSPIN,IELAST,IRMIN(NATYPD),IPAN(NATYPD),IRWS(NATYPD),
     &         KSRA,LMAX,IEND
       INTEGER IRMINIMP(NATOMIMP),IPANIMP(NATOMIMP),IRWSIMP(NATOMIMP),
     &         ATOMIMP(NATOMIMP)
       INTEGER NATOMIMP,IHOST,HOSTIMP(NATYPD)
       INTEGER LMMAXD
       PARAMETER (LMMAXD= (LMAXD+1)**2)
       INTEGER LMMAXSO
       PARAMETER (LMMAXSO=2*LMMAXD)
       INTEGER LMPOTD
       PARAMETER (LMPOTD= (LPOTD+1)**2)
       INTEGER IRMIND
       PARAMETER (IRMIND= IRMD-IRNSD)
       DOUBLE COMPLEX CZERO,CONE
       PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0))
       INTEGER IRCUT(0:IPAND,NATYPD),IRCUTIMP(0:IPAND,NATOMIMP)
       DOUBLE PRECISION C,Z(NATYPD),ZIMP(NATOMIMP),RCLSIMP(3,NATOMIMP)
       DOUBLE COMPLEX E
       DOUBLE PRECISION R(IRMD,NATYPD),CLEB(NCLEB,2),RIMP(IRMD,NATOMIMP)
       INTEGER ICLEB(NCLEB,4)
       DOUBLE PRECISION 
     &   VINS(IRMIND:IRMD,LMPOTD,NSPOTD),
     &   VM2Z(IRMD,NSPOTD),
     &   VINSIMP(IRMIND:IRMD,LMPOTD,NSPIN*NATOMIMP),
     &   VM2ZIMP(IRMD,NSPIN*NATOMIMP)
       DOUBLE COMPLEX  
     &   TMATLL((KORBIT+1)*LMMAXD,(KORBIT+1)*LMMAXD,IHOST),
     &   DTMTRX((KORBIT+1)*LMMAXD*NATOMIMP,(KORBIT+1)*LMMAXD*NATOMIMP),
     &   dummy_alphaget((KORBIT+1)*LMMAXD,(KORBIT+1)*LMMAXD)
       INTEGER I1,IR,NSRA,USE_SRATRICK,NVEC,LM1,LM2,ISPIN,I2,
     &         IL1,IL2
       INTEGER NPAN_LOG,NPAN_EQ,NPAN_INST(NATYPD),
     &         NCHEB,NPAN_TOT(NATYPD),IRMDNEWD
       integer, allocatable :: irmdnew(:)
       DOUBLE PRECISION R_LOG,THETA,PHI,THETAimp(NATOMIMP),
     &                  PHIimp(NATOMIMP),THETAhost(NATYPD),
     &                  PHIhost(NATYPD),PI
       DOUBLE COMPLEX GMATPREFACTOR
       DOUBLE PRECISION, ALLOCATABLE :: RNEW(:,:),RPAN_INTERVALL(:,:)
       INTEGER, ALLOCATABLE :: IPAN_INTERVALL(:,:)
       DOUBLE PRECISION, ALLOCATABLE :: VINSNEW(:,:,:)
       DOUBLE COMPLEX,ALLOCATABLE :: VNSPLL0(:,:,:),VNSPLL(:,:,:),
     &                               VNSPLL1(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: HLK(:,:),JLK(:,:),
     &                                HLK2(:,:),JLK2(:,:),
     &                                RLL(:,:,:),SLL(:,:,:),
     &                                RLLHOST(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: VNSIMP(:,:,:),VNSHOST(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: DELTABG(:,:,:),DELTASM(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: RADIALHOST(:,:),RADIALIMP(:,:),
     &                               VLLIMP(:,:),DELTAV(:,:),
     &                               TMATLLIMP(:,:,:),DELTAMTR(:,:,:),
     &                               DELTATMP(:),DELTAIMP(:,:)
       DOUBLE COMPLEX TMATSPH(2*(LMAX+1)), dummy_alpha(2*(LMAX+1))
       INTEGER JLK_INDEX(2*LMMAXSO),ipot,
     &         NPAN_LOGNEW(NATYPD),NPAN_EQNEW(NATYPD)
       !parallelization:
       INTEGER :: ierr, i1_start, i1_end, i1_start_imp, i1_end_imp
#ifdef CPP_MPI
       INTEGER :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
       double complex, allocatable :: temp(:,:,:),temp2(:,:,:,:) ! needed for MPI communication
#endif
       LOGICAL, EXTERNAL :: OPT

       PI=4.0D0*DATAN(1.0D0)
       if(myrank==master) WRITE(6,*) 'in tmatimp'
       IF (KSRA.GE.1) THEN
        NSRA=2
       ELSE
        NSRA=1
       ENDIF
       ALLOCATE(RLLHOST(NSRA*LMMAXSO,LMMAXSO,
     &                 IHOST,NTOTD*(NCHEB+1)))
       ALLOCATE(VNSHOST(NSRA*LMMAXSO,NSRA*LMMAXSO,
     &                  IHOST,NTOTD*(NCHEB+1)))
       ALLOCATE(DELTAV(LMMAXSO,LMMAXSO))
       ALLOCATE(TMATLLIMP(LMMAXSO,LMMAXSO,NATOMIMP))
       ALLOCATE(DELTAMTR(LMMAXSO,LMMAXSO,NATOMIMP))
       TMATLL    = CZERO
       RLLHOST   = CZERO
       VNSHOST   = CZERO
       TMATLLIMP = CZERO
       DELTAMTR  = CZERO
       DELTAV    = CZERO
       TMATSPH   = CZERO
       
       
       if(myrank==master) then
         ! read angles from nonco_ange files
         OPEN(UNIT=12,FILE='nonco_angle.dat',FORM='FORMATTED')
         DO I1=1,NATYPD
            READ(12,*) THETAhost(i1),PHIhost(i1)
            THETAhost(i1)=THETAhost(i1)/360.0D0*2.0D0*PI
            PHIhost(i1)=PHIhost(i1)/360.0D0*2.0D0*PI
         ENDDO
         CLOSE(12)
         OPEN(UNIT=13,FILE='nonco_angle_imp.dat',FORM='FORMATTED')
         DO I1=1,NATOMIMP
            READ(13,*) THETAimp(i1),PHIimp(i1)
            THETAimp(i1)=THETAimp(i1)/360.0D0*2.0D0*PI
            PHIimp(i1)=PHIimp(i1)/360.0D0*2.0D0*PI
         ENDDO
         CLOSE(13)
       endif
       
       
#ifdef CPP_MPI
       ! broadcast read-in values to all ranks (MPI_COMM_WORLD since
       ! atom dimension is solely used without energy parallelization)
       call MPI_Bcast(THETAhost, NATYPD, MPI_DOUBLE_PRECISION, master,
     &                MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error MPI_Bcast THETAhost in tmatimp'
       call MPI_Bcast(PHIhost, NATYPD, MPI_DOUBLE_PRECISION, master,
     &                MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error MPI_Bcast PHIhost in tmatimp'
       call MPI_Bcast(THETAimp, NATOMIMP, MPI_DOUBLE_PRECISION,master,
     &                MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error MPI_Bcast THETAimp in tmatimp'
       call MPI_Bcast(PHIimp, NATOMIMP, MPI_DOUBLE_PRECISION, master,
     &                MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error MPI_Bcast PHIimp in tmatimp'
       
       ! set start/end of parallel atom loops
       if(t_inc%i_write>0) write(1337, *) 'Parallelization host atoms'
       call distribute_linear_on_tasks(nranks, myrank, master, IHOST, 
     &                                  ntot_pT,ioff_pT,.true.)
       i1_start = ioff_pT(myrank) + 1
       i1_end   = ioff_pT(myrank) + ntot_pT(myrank)
       if(t_inc%i_write>0) write(1337, *) 'Parallelization imp. atoms'
       call distribute_linear_on_tasks(nranks,myrank,master,NATOMIMP,
     &                                  ntot_pT,ioff_pT,.true.)
       i1_start_imp = ioff_pT(myrank) + 1
       i1_end_imp   = ioff_pT(myrank) + ntot_pT(myrank)
#else  
       i1_start = 1
       i1_end = IHOST
       i1_start_imp = 1
       i1_end_imp = NATOMIMP
#endif

       if(OPT('GREENIMP')) then

         !11111111111111111111111111111111111111111111111111111111111111111111
         ! calculate tmat and radial wavefunctions of host atoms     
         !11111111111111111111111111111111111111111111111111111111111111111111
         
         ! create new mesh before loop starts
         ! data for the new mesh
         allocate(irmdnew(natypd), stat=ierr)
         if(ierr/=0) stop 'Error allocating irmdnew'
         IRMDNEWD = 0
         DO I1=1,NATYPD
           NPAN_INST(I1)= IPAN(I1)-1
           NPAN_TOT(I1)= NPAN_LOG+NPAN_EQ+NPAN_INST(I1)
           IF(NPAN_TOT(I1)*(NCHEB+1)>IRMDNEWD) THEN
              IRMDNEWD = NPAN_TOT(I1)*(NCHEB+1)
           ENDIF
           IRMDNEW(I1) = NPAN_TOT(I1)*(NCHEB+1)
         END DO
         ! new mesh
         ALLOCATE(RNEW(IRMDNEWD,NATYPD))
         ALLOCATE(RPAN_INTERVALL(0:NTOTD,NATYPD))
         ALLOCATE(IPAN_INTERVALL(0:NTOTD,NATYPD))
         ALLOCATE(VINSNEW(IRMDNEWD,LMPOTD,NSPOTD)) !NSPIND*max(NATYPD,NATOMIMP)))
         
         CALL CREATE_NEWMESH(NATYPD,LMAXD,LPOTD,IRMD,IRNSD,IPAND,
     &                       IRID,NTOTD,NFUND,NCHEB,IRMDNEWD,
     &                       NSPIN,R(:,:),IRMIN(:),IPAN(:),
     &                       IRCUT(0:IPAND,:),
     &                       R_LOG,NPAN_LOG,NPAN_EQ,NCHEB,
     &                       NPAN_LOGNEW(:),NPAN_EQNEW(:),
     &                       NPAN_TOT(:),RNEW(:,:),
     &                       RPAN_INTERVALL(0:NTOTD,:),
     &                       IPAN_INTERVALL(0:NTOTD,:), 1)
         
         !in second step interpolate potential (gain atom by atom with NATYPD==1)
         CALL INTERPOLATE_POTEN(LPOTD,IRMD,IRNSD,NATYPD,IPAND,NSPOTD,
     &          NTOTD,NCHEBD,IRMDNEWD,NSPIN,R(:,:),IRMIN(:),IRWS(:),
     &          IRCUT(0:IPAND,:),VINS(IRMIND:IRMD,1:LMPOTD,:),
     &          VM2Z(1:IRMD,:),NPAN_LOGNEW(:),NPAN_EQNEW(:),
     &          NPAN_TOT(:),RNEW(:,:),IPAN_INTERVALL(0:NTOTD,:),VINSNEW)

 
         ! calculate tmat and radial wavefunctions of host atoms      
         ! parallelized with MPI over atoms
         DO I2=i1_start, i1_end
          I1=HOSTIMP(I2)
          
          THETA = THETAhost(i1)
          PHI = PHIhost(i1)
            ISPIN=1
            IPOT=NSPIN*(I1-1)+1
            WRITE(6,*) 'HOST',I2,I1
         
          ! set up the non-spherical ll' matrix for potential VLL'
          IF (NSRA.EQ.2) THEN
           USE_SRATRICK=1
          ELSEIF (NSRA.EQ.1) THEN
           USE_SRATRICK=0
          ENDIF
         
          ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW(I1)))
          VNSPLL0=CZERO
          ! output potential onto which SOC is added
          ALLOCATE(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW(I1)))
          VNSPLL1=CZERO
         
          CALL VLLMAT(1,IRMDNEW(I1),IRMDNEW(I1),LMMAXD,LMMAXSO,VNSPLL0,
     &                VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &                CLEB,ICLEB,IEND,NSPIN,
     &                Z(I1),RNEW(1:IRMDNEW(I1),I1),USE_SRATRICK)

          ! contruct the spin-orbit coupling hamiltonian and add to potential
          CALL SPINORBIT_HAM(LMAX,LMMAXD,
     &          VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &          RNEW(1:IRMDNEW(I1),I1),E,Z(I1),C,t_params%SOCSCALE(I1),
     &          NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(0:NTOTD,I1),
     &          RPAN_INTERVALL(0:NTOTD,I1),NPAN_TOT(I1),NCHEB,
     &          IRMDNEW(I1),IRMDNEW(I1),VNSPLL0,
     &          VNSPLL1,'1')
          ! extend matrix for the SRA treatment
          IF (NSRA.EQ.2) THEN
           ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW(I1)))
           IF (USE_SRATRICK.EQ.0) THEN
            CALL VLLMATSRA(VNSPLL1,VNSPLL,RNEW(1:IRMDNEW(I1),I1),
     &               LMMAXSO,IRMDNEW(I1),IRMDNEW(I1),E,C,LMAX,0,'Ref=0')
           ELSEIF (USE_SRATRICK.EQ.1) THEN
              CALL VLLMATSRA(VNSPLL1,VNSPLL,RNEW(1:IRMDNEW(I1),I1),
     &            LMMAXSO,IRMDNEW(I1),IRMDNEW(I1),E,C,LMAX,0,'Ref=Vsph')
           ENDIF
          ELSE
           ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW(I1)))
           VNSPLL(:,:,:)=VNSPLL1(:,:,:)
          ENDIF
          
          ! calculate the source terms in the Lippmann-Schwinger equation
          ! these are spherical hankel and bessel functions
          ALLOCATE(HLK(1:4*(LMAX+1),IRMDNEW(I1)))
          ALLOCATE(JLK(1:4*(LMAX+1),IRMDNEW(I1)))
          ALLOCATE(HLK2(1:4*(LMAX+1),IRMDNEW(I1)))
          ALLOCATE(JLK2(1:4*(LMAX+1),IRMDNEW(I1)))
          HLK=CZERO
          JLK=CZERO
          HLK2=CZERO
          JLK2=CZERO
          GMATPREFACTOR=CZERO
          CALL RLLSLLSOURCETERMS(NSRA,NVEC,E,RNEW(1:IRMDNEW(I1),I1),
     &                           IRMDNEW(I1),IRMDNEW(I1),
     &                           LMAX,LMMAXSO,1,JLK_INDEX,HLK,JLK,HLK2,
     &                           JLK2,GMATPREFACTOR)
          ! using spherical potential as reference
          IF (USE_SRATRICK.EQ.1) THEN
          CALL CALCSPH(NSRA,IRMDNEW(I1),IRMDNEW(I1),LMAX,NSPIN,Z(I1),C,
     &              E,LMPOTD,LMMAXSO,RNEW(1:IRMDNEW(I1),I1),
     &              VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &              NCHEB,NPAN_TOT(I1),
     &              RPAN_INTERVALL(0:NTOTD,I1),JLK_INDEX,HLK,JLK,HLK2,
     &              JLK2,GMATPREFACTOR,TMATSPH,dummy_alpha,USE_SRATRICK)
          ENDIF
          
          ! calculate the tmat and wavefunctions
          ALLOCATE(RLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1)))
          ALLOCATE(SLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1)))
          RLL=CZERO
          SLL=CZERO
          
          ! right solutions
          CALL RLLSLL(RPAN_INTERVALL(0:NTOTD,I1),RNEW(1:IRMDNEW(I1),I1),
     &      VNSPLL,RLL,SLL,TMATLL(:,:,I2),NCHEB,
     &      NPAN_TOT(I1),LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW(I1),
     &      IRMDNEW(I1),NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,
     &      '1','1','0',USE_SRATRICK,dummy_alphaget)
          IF (NSRA.EQ.2) THEN
           DO IR=1,IRMDNEW(I1)
            DO LM1=1,LMMAXSO
             DO LM2=1,LMMAXSO
              RLL(LM1+LMMAXSO,LM2,IR)=
     &                RLL(LM1+LMMAXSO,LM2,IR)/C
              SLL(LM1+LMMAXSO,LM2,IR)=
     &                SLL(LM1+LMMAXSO,LM2,IR)/C
             ENDDO
            ENDDO
           ENDDO
          ENDIF
          ! save radial wavefunction for a host
          DO IR=1,IRMDNEW(I1)
           DO LM1=1,NVEC*LMMAXSO
            DO LM2=1,LMMAXSO
              RLLHOST(LM1,LM2,I2,IR)=RLL(LM1,LM2,IR)
            ENDDO
           ENDDO
          ENDDO
          ! add spherical contribution of tmatrix
          
          IF (USE_SRATRICK.EQ.1) THEN
           DO LM1=1,(KORBIT+1)*LMMAXD
            TMATLL(LM1,LM1,I2)=TMATLL(LM1,LM1,I2)+
     &                            TMATSPH(JLK_INDEX(LM1))
           ENDDO
          ENDIF
          
          ! rotate tmatrix and radial wavefunction to global frame
          CALL ROTATEMATRIX(TMATLL(1,1,I2),THETA,PHI,LMMAXD,0)
          
          ! create SRA potential for host
          ! set up the non-spherical ll' matrix for potential VLL'
          VNSPLL0=CZERO
          VNSPLL1=CZERO
          CALL VLLMAT(1,IRMDNEW(I1),IRMDNEW(I1),LMMAXD,LMMAXSO,VNSPLL0,
     &                VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &             CLEB,ICLEB,IEND,NSPIN,Z(I1),RNEW(1:IRMDNEW(I1),I1),0)
          
          ! contruct the spin-orbit coupling hamiltonian and add to potential
          CALL SPINORBIT_HAM(LMAX,LMMAXD,
     &             VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &           RNEW(1:IRMDNEW(I1),I1),E,Z(I1),C,t_params%SOCSCALE(I1),
     &             NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(0:NTOTD,I1),
     &             RPAN_INTERVALL(0:NTOTD,I1),NPAN_TOT(I1),NCHEB,
     &             IRMDNEW(I1),IRMDNEW(I1),VNSPLL0,
     &             VNSPLL1,'1')
          
          ! save potential for a host
          DO IR=1,IRMDNEW(I1)
           DO LM1=1,LMMAXSO
            DO LM2=1,LMMAXSO
             VNSHOST(LM1,LM2,I2,IR)=VNSPLL1(LM1,LM2,IR)
             IF (NSRA.EQ.2) THEN
              VNSHOST(LM1+LMMAXSO,LM2+LMMAXSO,I2,IR)=VNSPLL1(LM1,LM2,IR)
             ENDIF
            ENDDO
           ENDDO
          ENDDO
          
          DEALLOCATE(VNSPLL0)
          DEALLOCATE(VNSPLL1)
          DEALLOCATE(VNSPLL)
          DEALLOCATE(HLK)
          DEALLOCATE(JLK)
          DEALLOCATE(HLK2)
          DEALLOCATE(JLK2)
          DEALLOCATE(SLL,RLL)
         ENDDO ! I2
          
         DEALLOCATE(RNEW,RPAN_INTERVALL,IPAN_INTERVALL,VINSNEW)
          
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef   CPP_MPI
         ! collect results and write out only on master
         !communicate VNSHOST, RLLHOST, TMATLL
         i1 = NSRA*LMMAXSO*LMMAXSO*IHOST*NTOTD*(NCHEB+1)
         allocate(temp2(NSRA*LMMAXSO,LMMAXSO,IHOST,NTOTD*(NCHEB+1)), 
     &            stat=ierr)
         if(ierr/=0) stop 'Error allocating temp2 for RLLHOST' 
         temp2 = CZERO
         call MPI_ALLREDUCE(RLLHOST,temp2,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     &                      MPI_COMM_WORLD,ierr)
         if(ierr/=0) stop 'Error in MPI_Allreduce for RLLHOST'
         RLLHOST = temp2
         deallocate(temp2)
         
         i1 = NSRA*LMMAXSO*NSRA*LMMAXSO*IHOST*NTOTD*(NCHEB+1)
         allocate(temp2(NSRA*LMMAXSO,NSRA*LMMAXSO,IHOST,NTOTD*(NCHEB+1))
     &            , stat=ierr)
         if(ierr/=0) stop 'Error allocating temp2 for VNSHOST' 
         temp2 = CZERO
         call MPI_ALLREDUCE(VNSHOST,temp2,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     &                      MPI_COMM_WORLD,ierr)
         if(ierr/=0) stop 'Error in MPI_Allreduce for VNSHOST'
         VNSHOST = temp2
         deallocate(temp2)
         
         i1 = (KORBIT+1)*LMMAXD*(KORBIT+1)*LMMAXD*IHOST
         allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=ierr)
         if(ierr/=0) stop 'Error allocating temp for TMATLL MPI comm.' 
         temp = CZERO
         call MPI_ALLREDUCE(TMATLL,temp,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     &                      MPI_COMM_WORLD,ierr)
         if(ierr/=0) stop 'Error in MPI_Allreduce for TMATLL in tmatimp'
         TMATLL = temp
         deallocate(temp)
#endif

         ! write out DTMTRX file containgin Delta_T and Delta-matrices
         if(myrank==master) then
            IF (IELAST.EQ.1) THEN
             OPEN(UNIT=20,FILE='DTMTRX',FORM='FORMATTED')
             WRITE(20,'(I5)') NATOMIMP
             DO I1=1,NATOMIMP
              WRITE(20,'(3e17.9)') (RCLSIMP(I2,I1),I2=1,3)
             ENDDO
            ENDIF
         end if
         
         deallocate(irmdnew)

         !111 end 11111 end 11111 end 11111 end 11111 end 11111 end 11111 end 
         ! calculate tmat and radial wavefunctions of host atoms     
         !111 end 11111 end 11111 end 11111 end 11111 end 11111 end 11111 end 

       elseif(myrank==master) then

         write(*,*) 'skipping host atom loop in tmatimp_newsolver'

       end if ! (OPT('GREENIMP'))


       !22222222222222222222222222222222222222222222222222222222222222222222
       ! calculate tmat and radial wavefunctions of impurity atoms      
       !22222222222222222222222222222222222222222222222222222222222222222222

       ! create new mesh before loop starts
       ! data for the new mesh
       allocate(irmdnew(natomimp), stat=ierr)
       if(ierr/=0) stop 'Error allocating irmdnew'
       IRMDNEWD = 0
       DO I1=1,NATOMIMP
         NPAN_INST(I1)= IPANIMP(I1)-1
         NPAN_TOT(I1)= NPAN_LOG+NPAN_EQ+NPAN_INST(I1)
         IF(NPAN_TOT(I1)*(NCHEB+1)>IRMDNEWD) THEN
            IRMDNEWD = NPAN_TOT(I1)*(NCHEB+1)
         ENDIF
         IRMDNEW(I1) = NPAN_TOT(I1)*(NCHEB+1)
       END DO
       ! new mesh
       ALLOCATE(RNEW(IRMDNEWD,NATOMIMP))
       ALLOCATE(RPAN_INTERVALL(0:NTOTD,NATOMIMP))
       ALLOCATE(IPAN_INTERVALL(0:NTOTD,NATOMIMP))
       ALLOCATE(VINSNEW(IRMDNEWD,LMPOTD,NSPOTD))

       ! initialize with zeros
       TMATLLIMP = CZERO
       TMATSPH = CZERO

       CALL CREATE_NEWMESH(NATOMIMP,LMAXD,LPOTD,IRMD,IRNSD,IPAND,
     &            IRID,NTOTD,NFUND,NCHEB,IRMDNEWD,
     &            NSPIN,RIMP(:,1:NATOMIMP),
     &            IRMINIMP(1:NATOMIMP),IPANIMP(1:NATOMIMP),
     &            IRCUTIMP(0:IPAND,1:NATOMIMP),
     &            R_LOG,NPAN_LOG,NPAN_EQ,NCHEB,
     &            NPAN_LOGNEW(1:NATOMIMP),
     &            NPAN_EQNEW(1:NATOMIMP),
     &            NPAN_TOT(1:NATOMIMP),RNEW(1:IRMDNEWD,1:NATOMIMP),
     &            RPAN_INTERVALL(0:NTOTD,1:NATOMIMP),
     &            IPAN_INTERVALL(0:NTOTD,1:NATOMIMP), 1)

       !in second step interpolate potential 
       CALL INTERPOLATE_POTEN(LPOTD,IRMD,IRNSD,NATOMIMP,IPAND,NSPOTD,
     &            NTOTD,NCHEBD,IRMDNEWD,
     &            NSPIN,RIMP(:,1:NATOMIMP),IRMINIMP(1:NATOMIMP),
     &            IRWSIMP(1:NATOMIMP),
     &            IRCUTIMP(0:IPAND,1:NATOMIMP),
     &            VINSIMP(IRMIND:IRMD,1:LMPOTD,1:NATOMIMP),
     &            VM2ZIMP(1:IRMD,1:NATOMIMP),NPAN_LOGNEW(1:NATOMIMP),
     &            NPAN_EQNEW(1:NATOMIMP),NPAN_TOT(1:NATOMIMP),
     &            RNEW(1:IRMDNEWD,1:NATOMIMP),
     &            IPAN_INTERVALL(0:NTOTD,1:NATOMIMP),VINSNEW)


       ! now start loop over atoms
       DO I1=i1_start_imp,i1_end_imp
        THETA = THETAimp(i1)
        PHI = PHIimp(i1)
        ISPIN=1
        IPOT=NSPIN*(I1-1)+1
        WRITE(6,*) 'IMP',I1

        ALLOCATE(VNSIMP(NSRA*LMMAXSO,NSRA*LMMAXSO,IRMDNEW(I1)))
        VNSIMP=CZERO
        ! set up the non-spherical ll' matrix for potential VLL'
        IF (NSRA.EQ.2) THEN
         USE_SRATRICK=1
        ELSEIF (NSRA.EQ.1) THEN
         USE_SRATRICK=0
        ENDIF
        ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW(I1)))
        VNSPLL0=CZERO
        ! output potential onto which SOC is added
        ALLOCATE(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW(I1)))
        VNSPLL1=CZERO
        
        CALL VLLMAT(1,IRMDNEW(I1),IRMDNEW(I1),LMMAXD,LMMAXSO,VNSPLL0,
     &              VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &              CLEB,ICLEB,IEND,NSPIN,
     &              ZIMP(I1),RNEW(1:IRMDNEW(I1),I1),USE_SRATRICK)
        
        ! contruct the spin-orbit coupling hamiltonian and add to potential
        CALL SPINORBIT_HAM(LMAX,LMMAXD,
     &       VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &       RNEW(1:IRMDNEW(I1),I1),E,ZIMP(I1),C,t_params%SOCSCALE(I1),
     &       NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(0:NTOTD,I1),
     &       RPAN_INTERVALL(0:NTOTD,I1),NPAN_TOT(I1),NCHEB,
     &       IRMDNEW(I1),IRMDNEW(I1),VNSPLL0,
     &       VNSPLL1,'1')
        
        ! extend matrix for the SRA treatment
        IF (NSRA.EQ.2) THEN
         ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW(I1)))
         IF (USE_SRATRICK.EQ.0) THEN
          CALL VLLMATSRA(VNSPLL1,VNSPLL,
     &                  RNEW(1:IRMDNEW(I1),I1),LMMAXSO,IRMDNEW(I1),
     &                  IRMDNEW(I1),E,C,LMAX,0,'Ref=0')
         ELSEIF (USE_SRATRICK.EQ.1) THEN                    
          CALL VLLMATSRA(VNSPLL1,VNSPLL,
     &                  RNEW(1:IRMDNEW(I1),I1),LMMAXSO,IRMDNEW(I1),
     &                  IRMDNEW(I1),E,C,LMAX,0,'Ref=Vsph')
         ENDIF
        ELSE
         ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW(I1)))
         VNSPLL(:,:,:)=VNSPLL1(:,:,:)
        ENDIF
        
        ! calculate the source terms in the Lippmann-Schwinger equation
        ! these are spherical hankel and bessel functions
        ALLOCATE(HLK(1:4*(LMAX+1),IRMDNEW(I1)))
        ALLOCATE(JLK(1:4*(LMAX+1),IRMDNEW(I1)))
        ALLOCATE(HLK2(1:4*(LMAX+1),IRMDNEW(I1)))
        ALLOCATE(JLK2(1:4*(LMAX+1),IRMDNEW(I1)))
        HLK=CZERO
        JLK=CZERO
        HLK2=CZERO
        JLK2=CZERO
        GMATPREFACTOR=CZERO
        CALL RLLSLLSOURCETERMS(NSRA,NVEC,E,RNEW(1:IRMDNEW(I1),I1),
     &                         IRMDNEW(I1),IRMDNEW(I1),
     &                         LMAX,LMMAXSO,1,JLK_INDEX,
     &                         HLK,JLK,HLK2,JLK2,GMATPREFACTOR)
        ! using spherical potential as reference
        IF (USE_SRATRICK.EQ.1) THEN
        CALL CALCSPH(NSRA,IRMDNEW(I1),IRMDNEW(I1),LMAX,NSPIN,ZIMP(I1),C,
     &               E,LMPOTD,LMMAXSO,RNEW(1:IRMDNEW(I1),I1),
     &               VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &               NCHEB,NPAN_TOT(I1),RPAN_INTERVALL(0:NTOTD,I1),
     &               JLK_INDEX,HLK,JLK,HLK2,JLK2,
     &               GMATPREFACTOR,TMATSPH,dummy_alpha,USE_SRATRICK)
        ENDIF
        
        ! calculate the tmat and wavefunctions
        ALLOCATE(RLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1)))
        ALLOCATE(SLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW(I1)))
        RLL=CZERO
        SLL=CZERO
        
        ! right solutions
        CALL RLLSLL(RPAN_INTERVALL(0:NTOTD,I1),RNEW(1:IRMDNEW(I1),I1),
     &       VNSPLL,RLL,SLL,TMATLLIMP(:,:,I1),NCHEB,
     &       NPAN_TOT(I1),LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW(I1),
     &       IRMDNEW(I1),NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,
     &       '1','1','0',USE_SRATRICK,dummy_alphaget)
        IF (NSRA.EQ.2) THEN
         DO IR=1,IRMDNEW(I1)
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            RLL(LM1+LMMAXSO,LM2,IR)=
     &             RLL(LM1+LMMAXSO,LM2,IR)/C
            SLL(LM1+LMMAXSO,LM2,IR)=
     &             SLL(LM1+LMMAXSO,LM2,IR)/C
           ENDDO
          ENDDO
         ENDDO
        ENDIF

        ! for OPERATOR option save impurity wavefuncitons
        if (OPT('OPERATOR')) then
          t_imp%RLLIMP(:,:,:,i1) = RLL(:,:,:)
        end if
             
        ! add spherical contribution of tmatrix
        
        IF (USE_SRATRICK.EQ.1) THEN
         DO LM1=1,(KORBIT+1)*LMMAXD
          TMATLLIMP(LM1,LM1,I1)=TMATLLIMP(LM1,LM1,I1)+
     &                          TMATSPH(JLK_INDEX(LM1))
         ENDDO
        ENDIF
        
        ! rotate tmatrix and radial wavefunction to global frame
        
        CALL ROTATEMATRIX(TMATLLIMP(:,:,I1),THETA,PHI,LMMAXD,0)
        
        ! create SRA potential for impurity
        ! set up the non-spherical ll' matrix for potential VLL'
        VNSPLL0=CZERO
        CALL VLLMAT(1,IRMDNEW(I1),IRMDNEW(I1),LMMAXD,LMMAXSO,VNSPLL0,
     &              VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &              CLEB,ICLEB,IEND,NSPIN,
     &              ZIMP(I1),RNEW(1:IRMDNEW(I1),I1),0)
        
        ! contruct the spin-orbit coupling hamiltonian and add to potential
        CALL SPINORBIT_HAM(LMAX,LMMAXD,
     &       VINSNEW(1:IRMDNEW(I1),1:LMPOTD,IPOT:IPOT+NSPIN-1),
     &       RNEW(1:IRMDNEW(I1),I1),E,ZIMP(I1),C,t_params%SOCSCALE(I1),
     &       NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(0:NTOTD,I1),
     &       RPAN_INTERVALL(0:NTOTD,I1),NPAN_TOT(I1),NCHEB,
     &       IRMDNEW(I1),IRMDNEW(I1),VNSPLL0,
     &       VNSPLL1,'1')
        DO IR=1,IRMDNEW(I1)
         DO LM1=1,LMMAXSO
          DO LM2=1,LMMAXSO
           VNSIMP(LM1,LM2,IR)=VNSPLL1(LM1,LM2,IR)
           IF (NSRA.EQ.2) THEN
            VNSIMP(LM1+LMMAXSO,LM2+LMMAXSO,IR)=VNSPLL1(LM1,LM2,IR)
           ENDIF
          ENDDO
         ENDDO
        ENDDO
        
        ! calculate delta_t_imp matrix written in TMATLLIMP
        
        DO I2=1,IHOST
         IF (ATOMIMP(I1).EQ.HOSTIMP(I2)) THEN
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            TMATLLIMP(LM1,LM2,I1)=TMATLLIMP(LM1,LM2,I1)-
     &                            TMATLL(LM1,LM2,I2)
           ENDDO
          ENDDO
          DO LM1=1,NSRA*LMMAXSO
           DO LM2=1,NSRA*LMMAXSO
            DO IR=1,IRMDNEW(I1)
             VNSIMP(LM1,LM2,IR)=VNSIMP(LM1,LM2,IR)-
     &                             VNSHOST(LM1,LM2,I2,IR)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO
        
        ! calculate delta matrix \delta=int{R_imp*\deltaV*R_host}
        
        IF (IELAST.EQ.1) THEN
         ALLOCATE(DELTABG(LMMAXSO,LMMAXSO,IRMDNEW(I1)))
         ALLOCATE(DELTASM(LMMAXSO,LMMAXSO,IRMDNEW(I1)))
         DELTABG=CZERO
         DELTASM=CZERO
         ALLOCATE(DELTATMP(IRMDNEW(I1)))
         ALLOCATE(RADIALHOST(LMMAXSO,LMMAXSO))
         ALLOCATE(RADIALIMP(LMMAXSO,LMMAXSO))
         ALLOCATE(VLLIMP(LMMAXSO,LMMAXSO))
         DELTATMP=CZERO
        
         ! big component for SRA stored in DELTABG
         DO IR=1,IRMDNEW(I1)
          RADIALHOST=CZERO
          RADIALIMP=CZERO
          VLLIMP=CZERO 
          DELTAV=CZERO
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            DO I2=1,IHOST
             IF (ATOMIMP(I1).EQ.HOSTIMP(I2)) THEN
              RADIALHOST(LM1,LM2)=RLLHOST(LM1,LM2,I2,IR)
             ENDIF
            ENDDO !I2
            RADIALIMP(LM1,LM2)=RLL(LM1,LM2,IR)
            VLLIMP(LM1,LM2)=VNSIMP(LM1,LM2,IR)
           ENDDO !LM2
          ENDDO !LM1
        
          CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VLLIMP,
     &               LMMAXSO,RADIALIMP,LMMAXSO,CZERO,DELTAV,LMMAXSO)
          CALL ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTABG(:,:,IR),LMMAXSO)
        
          ! small component for SRA stored in DELTASM
          IF (NSRA.EQ.2) THEN
          RADIALHOST=CZERO
          RADIALIMP=CZERO
          VLLIMP=CZERO
          DELTAV=CZERO
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            DO I2=1,IHOST
             IF (ATOMIMP(I1).EQ.HOSTIMP(I2)) THEN
              RADIALHOST(LM1,LM2)=RLLHOST(LM1+LMMAXSO,LM2,I2,IR)
             ENDIF
            ENDDO
            RADIALIMP(LM1,LM2)=RLL(LM1+LMMAXSO,LM2,IR)
            VLLIMP(LM1,LM2)=VNSIMP(LM1+LMMAXSO,LM2+LMMAXSO,IR)
           ENDDO
          ENDDO 
          CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VLLIMP,
     &               LMMAXSO,RADIALIMP,LMMAXSO,CZERO,DELTAV,LMMAXSO)
          CALL ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTASM(:,:,IR),LMMAXSO)
        
          ! sum up big and small component stored in DELTABG
          DO LM1=1,LMMAXSO
            DO LM2=1,LMMAXSO
             DELTABG(LM1,LM2,IR)=DELTABG(LM1,LM2,IR)+DELTASM(LM1,LM2,IR)
            ENDDO
          ENDDO
        
          ENDIF ! NSRA
         ENDDO ! IR
        
         ! integrate 
         DO LM1=1,LMMAXSO
          DO LM2=1,LMMAXSO
           DO IR=1,IRMDNEW(I1)
            DELTATMP(IR)=DELTABG(LM1,LM2,IR)
           ENDDO
           CALL INTCHEB_CELL(DELTATMP,DELTAMTR(LM1,LM2,I1),
     &                       RPAN_INTERVALL(0:NTOTD,I1),
     &                       IPAN_INTERVALL(0:NTOTD,I1),
     &                       NPAN_TOT(I1),NCHEB,IRMDNEW(I1))
          ENDDO
         ENDDO
        
         DEALLOCATE(DELTATMP)
         DEALLOCATE(RADIALHOST)
         DEALLOCATE(RADIALIMP)
         DEALLOCATE(VLLIMP)
         DEALLOCATE(DELTABG,DELTASM)
       
        ENDIF ! IELAST.EQ.1
           
        DEALLOCATE(VNSPLL0)
        DEALLOCATE(VNSPLL1)
        DEALLOCATE(VNSPLL)
        DEALLOCATE(HLK)
        DEALLOCATE(JLK)
        DEALLOCATE(HLK2)
        DEALLOCATE(JLK2)
        DEALLOCATE(RLL,SLL)
        DEALLOCATE(VNSIMP)

       ENDDO ! I1 impurity

       DEALLOCATE(RNEW,RPAN_INTERVALL,IPAN_INTERVALL,VINSNEW)


       !222 end 22222 end 22222 end 22222 end 22222 end 22222 end 22222 end 
       ! calculate tmat and radial wavefunctions of impurity atoms      
       !222 end 22222 end 22222 end 22222 end 22222 end 22222 end 22222 end 


#ifdef CPP_MPI
       !collect results and write out only on master
       !collect TMATLLIMP, DELTAMTR 
       !communicate TMATLLIMP, DELTAMTR 
       i1 = LMMAXSO*LMMAXSO*NATOMIMP
       allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=ierr)
       temp = CZERO
       if(ierr/=0) stop 'Error allocating temp for TMATLLIMP' 
       call MPI_ALLREDUCE(TMATLLIMP,temp,i1,MPI_DOUBLE_COMPLEX,
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error in MPI_Allreduce TMATLLIMP in tmatimp'
       TMATLLIMP = temp
       deallocate(temp)
       
       i1 = LMMAXSO*LMMAXSO*NATOMIMP
       allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=ierr)
       temp = CZERO
       if(ierr/=0) stop 'Error allocating temp for DELTAMTR' 
       call MPI_ALLREDUCE(DELTAMTR,temp,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     &                    MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error in MPI_Allreduce for DELTAMTR'
       DELTAMTR = temp
       deallocate(temp)
#endif

       ! collect results and writeout only for GREENIMP option
       if(OPT('GREENIMP') .and. myrank==master) then

         DO I1=1,NATOMIMP
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            IL1=LMMAXSO*(I1-1)+LM1
            IL2=LMMAXSO*(I1-1)+LM2
            DTMTRX(IL1,IL2)=TMATLLIMP(LM1,LM2,I1)
           ENDDO
          ENDDO
         ENDDO

         IF (IELAST.EQ.1) THEN
          ALLOCATE(DELTAIMP((KORBIT+1)*LMMAXD*NATOMIMP,
     &                      (KORBIT+1)*LMMAXD*NATOMIMP))
          DELTAIMP=CZERO
          DO I1=1,NATOMIMP
           DO LM1=1,LMMAXSO
            DO LM2=1,LMMAXSO
             IL1=LMMAXSO*(I1-1)+LM1
             IL2=LMMAXSO*(I1-1)+LM2
             DELTAIMP(IL1,IL2)=DELTAMTR(LM1,LM2,I1)
            ENDDO
           ENDDO
          ENDDO

         ! write down to the file DTMTRX
          DO LM1=1,LMMAXSO*NATOMIMP
           DO LM2=1,LMMAXSO*NATOMIMP
             WRITE(20,'((2I5),(4e17.9))') LM2,LM1,DTMTRX(LM2,LM1),
     &                                    DELTAIMP(LM2,LM1)
           ENDDO
          ENDDO
          DEALLOCATE(DELTAIMP)
          if(myrank==master) WRITE(6,*) 'end of delta t'
         ENDIF ! IELAST.EQ.1
         
         CLOSE(20) ! output file DTMTRX

       end if ! myrank==master  
         
         
       DEALLOCATE(VNSHOST)
       DEALLOCATE(RLLHOST)
       DEALLOCATE(TMATLLIMP)
       DEALLOCATE(DELTAMTR)
       DEALLOCATE(DELTAV)
       deallocate(irmdnew)

       END 

