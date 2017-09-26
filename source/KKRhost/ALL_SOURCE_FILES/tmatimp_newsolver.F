       SUBROUTINE TMATIMP_NEWSOLVER(NSPIN,LMAX,R,Z,IELAST,E,KSRA,
     +                     IRWS,IPAN,IRCUT,IRMIN,C,CLEB,ICLEB,IEND,
     +                     NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,VINS,VM2Z,
     +                     NATOMIMP,RCLSIMP,ATOMIMP,IHOST,HOSTIMP,RIMP,
     +                     ZIMP,IRWSIMP,IPANIMP,IRCUTIMP,IRMINIMP,
     +                     VINSIMP,VM2ZIMP,DTMTRX)
#ifdef CPP_MPI
       use mpi
       use mod_mympi, only: myrank, master, nranks,
     +                      distribute_linear_on_tasks
#else
       use mod_mympi, only: myrank, master, nranks
#endif
       use mod_types, only: t_inc
       use mod_create_newmesh
       use mod_version_info
       use mod_wunfiles, only: t_params
       IMPLICIT NONE
c-----------------------------------------------------------------
c calculate and write down impurity tmatrix and delta matrix
c first calculate t-matrix for the host corresponding to imp. cluster
c N. H. Long, Juelich, 05.2013
c adapted to new routines (mainly changed interfaces) to work in KKRcode
c also added MPI parallelization
c P.R., Juelich, 09.2017
c-----------------------------------------------------------------
       include 'inc.p'
       INTEGER NSPIN,IELAST,IRMIN(NATYPD),IPAN(NATYPD),IRWS(NATYPD),
     +         KSRA,LMAX,IEND
       INTEGER IRMINIMP(NATOMIMP),IPANIMP(NATOMIMP),IRWSIMP(NATOMIMP),
     +         ATOMIMP(NATOMIMP)
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
     +   VINS(IRMIND:IRMD,LMPOTD,NSPOTD),
     +   VM2Z(IRMD,NSPOTD),
     +   VINSIMP(IRMIND:IRMD,LMPOTD,NSPIN*NATOMIMP),
     +   VM2ZIMP(IRMD,NSPIN*NATOMIMP)
       DOUBLE COMPLEX  
     +   TMATLL((KORBIT+1)*LMMAXD,(KORBIT+1)*LMMAXD,IHOST),
     +   DTMTRX((KORBIT+1)*LMMAXD*NATOMIMP,(KORBIT+1)*LMMAXD*NATOMIMP),
     +   dummy_alphaget((KORBIT+1)*LMMAXD,(KORBIT+1)*LMMAXD)
       INTEGER I1,IR,NSRA,USE_SRATRICK,NVEC,LM1,LM2,ISPIN,I2,
     +         IL1,IL2
       INTEGER NPAN_LOG,NPAN_EQ,NPAN_INST(NATYPD),
     +         NCHEB,NPAN_TOT(NATYPD),IRMDNEW
       DOUBLE PRECISION R_LOG,THETA,PHI,THETAimp(NATOMIMP),
     +                  PHIimp(NATOMIMP),THETAhost(NATYPD),
     +                  PHIhost(NATYPD),PI
       DOUBLE COMPLEX GMATPREFACTOR
       DOUBLE PRECISION, ALLOCATABLE :: RNEW(:,:),RPAN_INTERVALL(:,:)
       INTEGER, ALLOCATABLE :: IPAN_INTERVALL(:,:)
       DOUBLE PRECISION, ALLOCATABLE :: VINSNEW(:,:,:)
       DOUBLE COMPLEX,ALLOCATABLE :: VNSPLL0(:,:,:),VNSPLL(:,:,:),
     &                               VNSPLL1(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: HLK(:,:),JLK(:,:),
     +                                HLK2(:,:),JLK2(:,:),
     +                                RLL(:,:,:),SLL(:,:,:),
     +                                RLLHOST(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: VNSIMP(:,:,:),VNSHOST(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: DELTABG(:,:,:),DELTASM(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: RADIALHOST(:,:),RADIALIMP(:,:),
     +                               VLLIMP(:,:),DELTAV(:,:),
     +                               TMATLLIMP(:,:,:),DELTAMTR(:,:,:),
     +                               DELTATMP(:),DELTAIMP(:,:)
       DOUBLE COMPLEX TMATSPH(2*(LMAX+1)), dummy_alpha(2*(LMAX+1))
       INTEGER JLK_INDEX(2*LMMAXSO),ipot,
     +         NPAN_LOGNEW(NATYPD),NPAN_EQNEW(NATYPD)
       !parallelization:
       INTEGER :: ierr, i1_start, i1_end, i1_start_imp, i1_end_imp
#ifdef CPP_MPI
       INTEGER :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
       double complex, allocatable :: temp(:,:,:),temp2(:,:,:,:) ! needed for MPI communication
#endif


       PI=4.0D0*DATAN(1.0D0)
       WRITE(6,*) 'in tmatimp'
       IF (KSRA.GE.1) THEN
        NSRA=2
       ELSE
        NSRA=1
       ENDIF
       ALLOCATE(RLLHOST(NSRA*LMMAXSO,LMMAXSO,
     +                 IHOST,NTOTD*(NCHEB+1)))
       ALLOCATE(VNSHOST(NSRA*LMMAXSO,NSRA*LMMAXSO,
     +                  IHOST,NTOTD*(NCHEB+1)))
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
     +                MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error MPI_Bcast THETAhost in tmatimp'
       call MPI_Bcast(PHIhost, NATYPD, MPI_DOUBLE_PRECISION, master,
     +                MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error MPI_Bcast PHIhost in tmatimp'
       call MPI_Bcast(THETAimp, NATOMIMP, MPI_DOUBLE_PRECISION, master,
     +                MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error MPI_Bcast THETAimp in tmatimp'
       call MPI_Bcast(PHIimp, NATOMIMP, MPI_DOUBLE_PRECISION, master,
     +                MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error MPI_Bcast PHIimp in tmatimp'

       ! set start/end of parallel atom loops
       if(t_inc%i_write>0) write(1337, *) 'Parallelization host atoms:'
       call distribute_linear_on_tasks(nranks, myrank, master, NATYPD, 
     &                                  ntot_pT,ioff_pT,.true.)
       i1_start = ioff_pT(myrank) + 1
       i1_end   = ioff_pT(myrank) + ntot_pT(myrank)
       if(t_inc%i_write>0) write(1337, *) 'Parallelization imp. atoms:'
       call distribute_linear_on_tasks(nranks, myrank, master, NATOMIMP,
     &                                  ntot_pT,ioff_pT,.true.)
       i1_start_imp = ioff_pT(myrank) + 1
       i1_end_imp   = ioff_pT(myrank) + ntot_pT(myrank)
#else
       i1_start = 1
       i1_end = NATYPD
       i1_start_imp = 1
       i1_end_imp = NATOMIMP
#endif

!11111111111111111111111111111111111111111111111111111111111111111111
! calculate tmat and radial wavefunctions of host atoms     
!11111111111111111111111111111111111111111111111111111111111111111111

c create new mesh before loop starts
c data for the new mesh
       IRMDNEW = 0
       DO I1=1,NATYPD
         NPAN_INST(I1)= IPAN(I1)-1
         NPAN_TOT(I1)= NPAN_LOG+NPAN_EQ+NPAN_INST(I1)
         IF(NPAN_TOT(I1)*(NCHEB+1)>IRMDNEW) THEN
            IRMDNEW = NPAN_TOT(I1)*(NCHEB+1)
         ENDIF
       END DO
c new mesh
       ALLOCATE(RNEW(IRMDNEW,NATYPD))
       ALLOCATE(RPAN_INTERVALL(0:NTOTD,NATYPD))
       ALLOCATE(IPAN_INTERVALL(0:NTOTD,NATYPD))
       ALLOCATE(VINSNEW(IRMDNEW,LMPOTD,NSPOTD)) !NSPIND*max(NATYPD,NATOMIMP)))
       ! call create_newmesh here with NATYPD==1 to use it atom by atom
!       open(756149,file='input_create_newmesh', form='formatted')
!       write(756149,*) 'params',NATYPD,LMAXD,LPOTD,IRMD,IRNSD,IPAND,
!     +                     IRID,NTOTD,NFUND,NCHEBD,NSPIN,IRMD
!       write(756149,*) 'atom loop'
!       do i2=1,IHOST
!         i1 = hostimp(i2)
!         write(756149,*) 'R',i1,IRMDNEW,NPAN_LOG(1),NPAN_EQ(1),R(:,i1)
!         write(756149,*) 'IRMIN',i1,IRMIN(i1)
!         write(756149,*) 'IPAN',i1,IPAN(i1)
!         write(756149,*) 'IRCUT',i1,IRCUT(0:,i1)
!         write(756149,*) 'NPAN_EQNEW',i1,NPAN_EQNEW(i1)
!         write(756149,*) 'NPAN_LOGNEW',i1,NPAN_LOGNEW(i1)
!         write(756149,*) 'NPAN_TOT',i1,NPAN_TOT(i1)
!       end do
!       close(756149)


!       open(757547,file='calling_list_create_newmesh',form='formatted')
!       write(757547,*)NATYPD,LMAXD,LPOTD,IRMD,IRNSD,IPAND
!       write(757547,*)'next...'
!       write(757547,*)IRID,NTOTD,NFUND,NCHEBD,IRMDNEW
!       write(757547,*)'next...'
!       write(757547,*)NSPIN,R,IRMIN,IPAN,IRCUT
!       write(757547,*)'next...'
!       write(757547,*)R_LOG,NPAN_LOG(1),NPAN_EQ(1),NCHEB
!       write(757547,*)'next...'
!       write(757547,*)NPAN_LOGNEW,NPAN_EQNEW
!       write(757547,*)'next...'
!       write(757547,*)NPAN_TOT,RNEW,RPAN_INTERVALL,IPAN_INTERVALL
!       write(757547,*)'next...'
!       close(757547)


       CALL CREATE_NEWMESH(NATYPD,LMAXD,LPOTD,IRMD,IRNSD,IPAND,
     +                     IRID,NTOTD,NFUND,NCHEB,IRMDNEW,
     +                     NSPIN,R(:,:),IRMIN(:),IPAN(:),
     +                     IRCUT(0:IPAND,:),
     +                     R_LOG,NPAN_LOG,NPAN_EQ,NCHEB,
     +                     NPAN_LOGNEW(:),NPAN_EQNEW(:),
     +                     NPAN_TOT(:),RNEW(:,:),
     +                     RPAN_INTERVALL(0:,:),
     +                     IPAN_INTERVALL(0:,:), 1)
!       open(756149,file='output_create_newmesh', form='formatted')
!       write(756149,*) 'atom loop'
!       do i2=1,IHOST
!         i1 = hostimp(i2)
!         write(756149,*) 'R',i1,NPAN_LOG(1),NPAN_EQ(1),R(:,i1)
!         write(756149,*)
!         write(756149,*) 'IRMIN',i1,IRMIN(i1)
!         write(756149,*)
!        write(756149,*) 'IPAN',i1,IPAN(i1)
!        write(756149,*)
!        write(756149,*) 'IRCUT',i1,IRCUT(0:,i1)
!        write(756149,*)
!        write(756149,*) 'NPAN_EQNEW',i1,NPAN_EQNEW(i1)
!        write(756149,*)
!        write(756149,*) 'NPAN_LOGNEW',i1,NPAN_LOGNEW(i1)
!        write(756149,*)
!        write(756149,*) 'NPAN_TOT',i1,NPAN_TOT(i1)
!        write(756149,*)
!        write(756149,*) 'RNEW',i1,
!    + RNEW(1:(NPAN_LOG(1)+NPAN_EQ(1)+IPAN(I1)-1)*(NCHEB+1),i1)
!        write(756149,*)
!        write(756149,*) 'RPAN_INTERVALL',i1,
!    + RPAN_INTERVALL(0:NPAN_LOG(1)+NPAN_EQ(1)+IPAN(I1)-1,i1)
!        write(756149,*)
!         write(756149,*) 'IPAN_INTERVALL',i1,
!     + IPAN_INTERVALL(0:NPAN_LOG(1)+NPAN_EQ(1)+IPAN(I1)-1,i1)
!       end do
!       close(756149)

       !in second step interpolate potential (gain atom by atom with NATYPD==1)
       CALL INTERPOLATE_POTEN(LPOTD,IRMD,IRNSD,NATYPD,IPAND,NSPOTD,
     +             NTOTD,NCHEBD,IRMDNEW,
     +             NSPIN,R(:,:),IRMIN(:),IRWS(:),
     +             IRCUT(0:IPAND,:),
     +             VINS(IRMIND:IRMD,1:LMPOTD,:),
     +             VM2Z(1:IRMD,:),NPAN_LOGNEW(:),NPAN_EQNEW(:),
     +             NPAN_TOT(:),RNEW(:,:),IPAN_INTERVALL,VINSNEW)
!     +             VM2Z(1:IRMD,:),NPAN_LOG(:),NPAN_EQ(:),
!     +             NPAN_TOT(:),RNEW(:,:),IPAN_INTERVALL,VINSNEW)

       
!      open(756149,file='output_interpolat_poten', form='formatted')
!      write(756149,*) 'atom loop'
!      do i2=1,IHOST
!        i1 = hostimp(i2)
!        IPOT=NSPIN*(I1-1)+1
!        write(756149,*) 'RNEW',i1,
!    + RNEW(1:(NPAN_LOG(1)+NPAN_EQ(1)+IPAN(I1)-1)*(NCHEB+1),i1)
!        write(756149,*) 'VINS',i1,IRMIND,IRMD,LMPOTD,IPOT,NSPIN
!        write(756149,*)
!    +             vins(IRMIND:IRMD,1:LMPOTD,IPOT:IPOT+NSPIN-1)
!        write(756149,*) 'next...'
!        write(756149,*)
!        write(756149,*) 'VM2Z',i1,IRMD,IPOT,NSPIN
!        write(756149,*)
!    +          vm2z(1:IRMD,IPOT:IPOT+NSPIN-1)
!        write(756149,*) 'next...'
!        write(756149,*)
!        write(756149,*) 'VINSNEW',i1,IRMDNEW,LMPOTD,NSPIND
!        write(756149,*)
!    +          vinsnew(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1)
!        write(756149,*)
!      end do
!      close(756149)
 
 
c calculate tmat and radial wavefunctions of host atoms      
! parallelized with MPI over atoms
      DO I1=i1_start, i1_end
       THETA = THETAhost(i1)
       PHI = PHIhost(i1)
       DO I2=1,IHOST
        IF (I1.EQ.HOSTIMP(I2)) THEN
         ISPIN=1
         IPOT=NSPIN*(I1-1)+1
         WRITE(6,*) 'HOST',I2,I1,IPOT,allocated(vinsnew)

c set up the non-spherical ll' matrix for potential VLL'
       IF (NSRA.EQ.2) THEN
        USE_SRATRICK=1
       ELSEIF (NSRA.EQ.1) THEN
        USE_SRATRICK=0
       ENDIF

       ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW))
       VNSPLL0=CZERO
       ! output potential onto which SOC is added
       ALLOCATE(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW))
       VNSPLL1=CZERO

!      if(i1==hostimp(1))
!    +     open(66666,file='vllmat_input1',form='formatted')
!      write(66666,*) 'params'
!      write(66666,*) IRMDNEW,IRMDNEW,LMMAXD,LMMAXSO
!      write(66666,*) 'next...'
!      write(66666,*)
!      write(66666,*) 'vnspll0'
!      write(66666,*) VNSPLL0
!      write(66666,*) 'next...'
!      write(66666,*)
!      write(66666,*) 'vinsnew',i1,IRMDNEW,LMPOTD,NSPIN
!      write(66666,*) 
!    +          vinsnew(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1)
!      write(66666,*) 'next...'
!      write(66666,*)
!      write(66666,*) 'icleb,iend.nspin'
!      write(66666,*) ICLEB,IEND,NSPIN
!      write(66666,*) 'next...'
!      write(66666,*)
!      write(66666,*) 'Z,rnew,use_sratrick'
!      write(66666,*) Z(I1),RNEW,USE_SRATRICK
!      write(66666,*) 'next...'
!       write(66666,*)

       CALL VLLMAT(1,IRMDNEW,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,
     +             VINSNEW(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1),
     +             CLEB,ICLEB,IEND,NSPIN,
     +             Z(I1),RNEW(:,I1),USE_SRATRICK)

!      if(i1==hostimp(1)) 
!    +     open(66667,file='spinorbitham_input1',form='formatted')
!      write(66667,*) 'params'
!      write(66667,*) LMAX,LMMAXD
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'vinsnew',i1,IRMDNEW,LMPOTD,NSPIN
!      write(66667,*) 
!    +          vinsnew(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1)
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'rnew'
!      write(66667,*)  RNEW(:,I1)
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'e,z,c,socscale'
!      write(66667,*) E,Z(I1),C,t_params%SOCSCALE(I1)
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'nspin,lmpotd,...'
!      write(66667,*) NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(:,I1)
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'rpan_intervall,...'
!      write(66667,*) RPAN_INTERVALL(:,I1),NPAN_TOT(I1),NCHEB
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'npan_tot...,irmdnew'
!      write(66667,*) NPAN_TOT(I1)*(NCHEB+1),IRMDNEW
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'vnspll0'
!      write(66667,*) VNSPLL0
!       write(66667,*) 'next...'
!       write(66667,*)
!       write(66667,*) 'vnspll1'
!       write(66667,*) VNSPLL1
!       write(66667,*) 'next...'
!       write(66667,*)
c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW(:,:,IPOT:IPOT+NSPIN-1),
     +                    RNEW(:,I1),E,Z(I1),C,t_params%SOCSCALE(I1),
     +                    NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(:,I1),
     +                    RPAN_INTERVALL(:,I1),NPAN_TOT(I1),NCHEB,
     +                    NPAN_TOT(I1)*(NCHEB+1),IRMDNEW,VNSPLL0,
     +                    VNSPLL1,'1')
c extend matrix for the SRA treatment
       IF (NSRA.EQ.2) THEN
!      if(i1==hostimp(1)) 
!    +     open(66668,file='vllmatsra_input1',form='formatted')
!      write(66668,*) 'params'
!      write(66668,*) LMMAXSO,IRMDNEW,E,C,LMAX
!      write(66668,*) 'next...'
!      write(66668,*)
!      write(66668,*) 'vnspll1'
!      write(66668,*) vnspll1
!      write(66668,*) 'next...'
!      write(66668,*)
        ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW))
        IF (USE_SRATRICK.EQ.0) THEN
         CALL VLLMATSRA(VNSPLL1,VNSPLL,RNEW(:,I1),LMMAXSO,IRMDNEW,
     +                  IRMDNEW,E,C,LMAX,0,'Ref=0')
        ELSEIF (USE_SRATRICK.EQ.1) THEN
         CALL VLLMATSRA(VNSPLL1,VNSPLL,RNEW(:,I1),LMMAXSO,IRMDNEW,
     +                  IRMDNEW,E,C,LMAX,0,'Ref=Vsph')
        ENDIF
       ELSE
        ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW))
        VNSPLL(:,:,:)=VNSPLL1(:,:,:)
       ENDIF

c calculate the source terms in the Lippmann-Schwinger equation
c these are spherical hankel and bessel functions
       ALLOCATE(HLK(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(JLK(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(HLK2(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(JLK2(1:4*(LMAX+1),IRMDNEW))
       HLK=CZERO
       JLK=CZERO
       HLK2=CZERO
       JLK2=CZERO
       GMATPREFACTOR=CZERO
       CALL RLLSLLSOURCETERMS(NSRA,NVEC,E,RNEW(:,I1),IRMDNEW,IRMDNEW,
     +                        LMAX,LMMAXSO,1,JLK_INDEX,HLK,JLK,HLK2,
     +                        JLK2,GMATPREFACTOR)
c using spherical potential as reference
       IF (USE_SRATRICK.EQ.1) THEN
!      if(i1==hostimp(1)) 
!    +     open(66669,file='calcsph_input1',form='formatted')
!      write(66669,*) 'params'
!      write(66669,*) NSRA,IRMDNEW,LMAX,NSPIN,Z(I1),C,E,LMPOTD,LMMAXSO
!      write(66669,*) 'next...'
!      write(66669,*)
!      write(66669,*) 'HLK'
!      write(66669,*) HLK
!       write(66669,*) 'next...'
!      write(66669,*)
!      write(66669,*) 'JLK'
!      write(66669,*) JLK
!      write(66669,*) 'next...'
!      write(66669,*)
!      write(66669,*) 'HLK2'
!      write(66669,*) HLK2
!      write(66669,*) 'next...'
!      write(66669,*)
!      write(66669,*) 'JLK2'
!      write(66669,*) JLK2
!      write(66669,*) 'next...'
!      write(66669,*)
!      write(66669,*) 'GMATPREFACTOR'
!      write(66669,*) GMATPREFACTOR
!      write(66669,*) 'next...'
!       write(66669,*)
       CALL CALCSPH(NSRA,IRMDNEW,IRMDNEW,LMAX,NSPIN,Z(I1),C,E,LMPOTD,
     +              LMMAXSO,RNEW(:,I1),VINSNEW(:,:,IPOT:IPOT+NSPIN-1),
     +              NCHEB,NPAN_TOT(I1),
     +              RPAN_INTERVALL(:,I1),JLK_INDEX,HLK,JLK,HLK2,JLK2,
     +              GMATPREFACTOR,TMATSPH,dummy_alpha,USE_SRATRICK)
       ENDIF

c calculate the tmat and wavefunctions
       ALLOCATE(RLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW))
       ALLOCATE(SLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW))
       RLL=CZERO
       SLL=CZERO

c right solutions
       CALL RLLSLL(RPAN_INTERVALL,RNEW(:,I1),VNSPLL,RLL,
     +          SLL,TMATLL(1,1,I2),NCHEB,
     +          NPAN_TOT(I1),LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +          IRMDNEW,NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,
     +          '1','1','0',USE_SRATRICK,dummy_alphaget)
       IF (NSRA.EQ.2) THEN
        DO IR=1,IRMDNEW
         DO LM1=1,LMMAXSO
          DO LM2=1,LMMAXSO
           RLL(LM1+LMMAXSO,LM2,IR)=
     +             RLL(LM1+LMMAXSO,LM2,IR)/C
           SLL(LM1+LMMAXSO,LM2,IR)=
     +             SLL(LM1+LMMAXSO,LM2,IR)/C
          ENDDO
         ENDDO
        ENDDO
       ENDIF
c save radial wavefunction for a host
       DO IR=1,IRMDNEW
        DO LM1=1,NVEC*LMMAXSO
         DO LM2=1,LMMAXSO
           RLLHOST(LM1,LM2,I2,IR)=RLL(LM1,LM2,IR)
         ENDDO
        ENDDO
       ENDDO
c add spherical contribution of tmatrix

       IF (USE_SRATRICK.EQ.1) THEN
        DO LM1=1,(KORBIT+1)*LMMAXD
         TMATLL(LM1,LM1,I2)=TMATLL(LM1,LM1,I2)+
     &                         TMATSPH(JLK_INDEX(LM1))
        ENDDO
       ENDIF

c rotate tmatrix and radial wavefunction to global frame
       CALL ROTATEMATRIX(TMATLL(1,1,I2),THETA,PHI,LMMAXD,0)

c create SRA potential for host
c set up the non-spherical ll' matrix for potential VLL'
       VNSPLL0=CZERO
       VNSPLL1=CZERO
       CALL VLLMAT(1,IRMDNEW,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,
!     +             VINSNEW(:,:,IPOT:IPOT+NSPIN-1),
     +             VINSNEW(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1),
     +             CLEB,ICLEB,IEND,NSPIN,Z(I1),RNEW(:,I1),0)
!     +             CLEB,ICLEB,IEND,NSPIN,Z(I1),RNEW(:,I1),USE_SRATRICK)

c contruct the spin-orbit coupling hamiltonian and add to potential
!     if(i1==hostimp(1)) 
!    +     open(76667,file='spinorbitham_input3',form='formatted')
!     write(76667,*) 'params'
!     write(76667,*) LMAX,LMMAXD
!     write(76667,*) 'next...'
!     write(76667,*)
!     write(76667,*) 'vinsnew',i1,IRMDNEW,LMPOTD,NSPIN
!     write(76667,*) 
!    +          vinsnew(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1)
!     write(76667,*) 'next...'
!     write(76667,*)
!     write(76667,*) 'rnew'
!     write(76667,*)  RNEW(:,I1)
!     write(76667,*) 'next...'
!     write(76667,*)
!     write(76667,*) 'e,z,c,socscale'
!     write(76667,*) E,Z(I1),C,t_params%SOCSCALE(I1)
!     write(76667,*) 'next...'
!     write(76667,*)
!     write(76667,*) 'nspin,lmpotd,...'
!     write(76667,*) NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(:,I1)
!     write(76667,*) 'next...'
!     write(76667,*)
!     write(76667,*) 'rpan_intervall,...'
!     write(76667,*) RPAN_INTERVALL(:,I1),NPAN_TOT(I1),NCHEB
!     write(76667,*) 'next...'
!     write(76667,*)
!     write(76667,*) 'npan_tot...,irmdnew'
!     write(76667,*) NPAN_TOT(I1)*(NCHEB+1),IRMDNEW
!     write(76667,*) 'next...'
!     write(76667,*)
!     write(76667,*) 'vnspll0'
!     write(76667,*) VNSPLL0
!      write(76667,*) 'next...'
!      write(76667,*)
       CALL SPINORBIT_HAM(LMAX,LMMAXD,
     +                    VINSNEW(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1),
!    +          vinsnew(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1)
     +                    RNEW(:,I1),E,Z(I1),C,t_params%SOCSCALE(I1),
     +                    NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(:,I1),
     +                    RPAN_INTERVALL(:,I1),NPAN_TOT(I1),NCHEB,
     +                    NPAN_TOT(I1)*(NCHEB+1),IRMDNEW,VNSPLL0,
     +                    VNSPLL1,'1')
!     if(i1==hostimp(1)) 
!    +     open(76668,file='vllmatsra_input3',form='formatted')
!     write(76668,*) 'params'
!     write(76668,*) LMMAXSO,IRMDNEW,E,C,LMAX
!     write(76668,*) 'next...'
!     write(76668,*)
!     write(76668,*) 'vnspll1'
!     write(76668,*) vnspll1
!     write(76668,*) 'next...'
!     write(76668,*)
c save potential for a host
       DO IR=1,IRMDNEW
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          VNSHOST(LM1,LM2,I2,IR)=VNSPLL1(LM1,LM2,IR)
          IF (NSRA.EQ.2) THEN
           VNSHOST(LM1+LMMAXSO,LM2+LMMAXSO,I2,IR)=VNSPLL1(LM1,LM2,IR)
!           if(ir<20) then
!            write(75613,*) i1,ipot,VNSHOST(LM1,LM2,I2,IR)
!            write(75613,*) VNSHOST(LM1+LMMAXSO,LM2+LMMAXSO,I2,IR)
!           end if
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
         ENDIF
        ENDDO ! I1
       ENDDO ! I2

       DEALLOCATE(RNEW,RPAN_INTERVALL,IPAN_INTERVALL,VINSNEW)

       !test
!       close(66666)
!       close(66667)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CPP_MPI
       ! collect results and write out only on master
       !communicate VNSHOST, RLLHOST, TMATLL
       i1 = NSRA*LMMAXSO*LMMAXSO*IHOST*NTOTD*(NCHEB+1)
       allocate(temp2(NSRA*LMMAXSO,LMMAXSO,IHOST,NTOTD*(NCHEB+1)), 
     +          stat=ierr)
       if(ierr/=0) stop 'Error allocating temp2 for RLLHOST MPI comm.' 
       temp2 = CZERO
       call MPI_ALLREDUCE(RLLHOST,temp2,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                    MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error in MPI_Allreduce for RLLHOST in tmatimp'
       RLLHOST = temp2
       deallocate(temp2)

       i1 = NSRA*LMMAXSO*NSRA*LMMAXSO*IHOST*NTOTD*(NCHEB+1)
       allocate(temp2(NSRA*LMMAXSO,NSRA*LMMAXSO,IHOST,NTOTD*(NCHEB+1)),
     +          stat=ierr)
       if(ierr/=0) stop 'Error allocating temp2 for VNSHOST MPI comm.' 
       temp2 = CZERO
       call MPI_ALLREDUCE(VNSHOST,temp2,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                    MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error in MPI_Allreduce for VNSHOST in tmatimp'
       VNSHOST = temp2
       deallocate(temp2)

       i1 = (KORBIT+1)*LMMAXD*(KORBIT+1)*LMMAXD*IHOST
       allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=ierr)
       if(ierr/=0) stop 'Error allocating temp for TMATLL MPI comm.' 
       temp = CZERO
       call MPI_ALLREDUCE(TMATLL,temp,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                    MPI_COMM_WORLD,ierr)
       if(ierr/=0) stop 'Error in MPI_Allreduce for TMATLL in tmatimp'
       TMATLL = temp
       deallocate(temp)

       Write(*,*) myrank, TMATLL(32,1,IHOST)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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

!111 end 11111 end 11111 end 11111 end 11111 end 11111 end 11111 end 
! calculate tmat and radial wavefunctions of host atoms     
!111 end 11111 end 11111 end 11111 end 11111 end 11111 end 11111 end 


!22222222222222222222222222222222222222222222222222222222222222222222
! calculate tmat and radial wavefunctions of impurity atoms      
!22222222222222222222222222222222222222222222222222222222222222222222

c create new mesh before loop starts
c data for the new mesh
       IRMDNEW = 0
       DO I1=1,NATOMIMP
         NPAN_INST(I1)= IPANIMP(I1)-1
         NPAN_TOT(I1)= NPAN_LOG+NPAN_EQ+NPAN_INST(I1)
         IF(NPAN_TOT(I1)*(NCHEB+1)>IRMDNEW) THEN
            IRMDNEW = NPAN_TOT(I1)*(NCHEB+1)
         ENDIF
       END DO
c new mesh
       ALLOCATE(RNEW(IRMDNEW,NATOMIMP))
       ALLOCATE(RPAN_INTERVALL(0:NTOTD,NATOMIMP))
       ALLOCATE(IPAN_INTERVALL(0:NTOTD,NATOMIMP))
       ALLOCATE(VINSNEW(IRMDNEW,LMPOTD,NSPOTD))

       ! initialize with zeros
       TMATLLIMP = CZERO
       TMATSPH = CZERO

!      open(756149,file='input_create_newmesh2', form='formatted')
!      write(756149,*) 'params',NATYPD,LMAXD,LPOTD,IRMD,IRNSD,IPAND,
!    +                     IRID,NTOTD,NFUND,NCHEBD,NSPIN,IRMD
!      write(756149,*) 'atom loop'
!      do i1=1,natomimp
!        write(756149,*) 'R',i1,IRMDNEW,NPAN_LOG(1),NPAN_EQ(1),
!    +                    RIMP(:,i1)
!        write(756149,*) 'IRMIN',i1,IRMINIMP(i1)
!        write(756149,*) 'IPAN',i1,IPANIMP(i1)
!        write(756149,*) 'IRCUT',i1,IRCUTIMP(0:,i1)
!        write(756149,*) 'NPAN_EQNEW',i1,NPAN_EQNEW(i1)
!         write(756149,*) 'NPAN_LOGNEW',i1,NPAN_LOGNEW(i1)
!         write(756149,*) 'NPAN_TOT',i1,NPAN_TOT(i1)
!       end do
!       close(756149)

       CALL CREATE_NEWMESH(NATOMIMP,LMAXD,LPOTD,IRMD,IRNSD,IPAND,
     +                     IRID,NTOTD,NFUND,NCHEB,IRMDNEW,
     +                     NSPIN,RIMP(:,1:NATOMIMP),
     +                     IRMINIMP(1:NATOMIMP),IPANIMP(1:NATOMIMP),
     +                     IRCUTIMP(0:IPAND,1:NATOMIMP),
     +                     R_LOG,NPAN_LOG,NPAN_EQ,NCHEB,
     +                     NPAN_LOGNEW(1:NATOMIMP),
     +                     NPAN_EQNEW(1:NATOMIMP),
     +                     NPAN_TOT(1:NATOMIMP),RNEW(:,1:NATOMIMP),
     +                     RPAN_INTERVALL(0:NTOTD,1:NATOMIMP),
     +                     IPAN_INTERVALL(0:NTOTD,1:NATOMIMP), 1)

!      open(756149,file='output_create_newmesh2', form='formatted')
!      write(756149,*) 'atom loop'
!      do i1=1,Natomimp
!        write(756149,*) 'R',i1,NPAN_LOG(1),NPAN_EQ(1),RIMP(:,i1)
!        write(756149,*)
!        write(756149,*) 'IRMIN',i1,IRMINIMP(i1)
!        write(756149,*)
!        write(756149,*) 'IPAN',i1,IPANIMP(i1)
!        write(756149,*)
!        write(756149,*) 'IRCUT',i1,IRCUTIMP(0:,i1)
!        write(756149,*)
!        write(756149,*) 'NPAN_EQNEW',i1,NPAN_EQNEW(i1)
!        write(756149,*)
!        write(756149,*) 'NPAN_LOGNEW',i1,NPAN_LOGNEW(i1)
!        write(756149,*)
!        write(756149,*) 'NPAN_TOT',i1,NPAN_TOT(i1)
!        write(756149,*)
!        write(756149,*) 'RNEW',i1,
!    + RNEW(1:(NPAN_LOG(1)+NPAN_EQ(1)+IPANIMP(I1)-1)*(NCHEB+1),i1)
!        write(756149,*)
!        write(756149,*) 'RPAN_INTERVALL',i1,
!    + RPAN_INTERVALL(0:NPAN_LOG(1)+NPAN_EQ(1)+IPANIMP(I1)-1,i1)
!        write(756149,*)
!        write(756149,*) 'IPAN_INTERVALL',i1,
!    + IPAN_INTERVALL(0:NPAN_LOG(1)+NPAN_EQ(1)+IPANIMP(I1)-1,i1)
!      end do
!      close(756149)
 
       !in second step interpolate potential 
       CALL INTERPOLATE_POTEN(LPOTD,IRMD,IRNSD,NATOMIMP,IPAND,NSPOTD,
     +            NTOTD,NCHEBD,IRMDNEW,
     +            NSPIN,RIMP(:,1:NATOMIMP),IRMINIMP(1:NATOMIMP),
     +            IRWSIMP(1:NATOMIMP),
     +            IRCUTIMP(0:IPAND,1:NATOMIMP),
     +            VINSIMP(IRMIND:IRMD,1:LMPOTD,1:NATOMIMP),
!     +            VM2ZIMP(1:IRMD,1:NATOMIMP),NPAN_LOG(1:NATOMIMP),
!     +            NPAN_EQ(1:NATOMIMP),NPAN_TOT(1:NATOMIMP),
     +            VM2ZIMP(1:IRMD,1:NATOMIMP),NPAN_LOGNEW(1:NATOMIMP),
     +            NPAN_EQNEW(1:NATOMIMP),NPAN_TOT(1:NATOMIMP),
     +            RNEW(:,1:NATOMIMP),IPAN_INTERVALL(0:NTOTD,1:NATOMIMP),
     +            VINSNEW)

!      open(756149,file='output_interpolat_poten2', form='formatted')
!      write(756149,*) 'atom loop'
!      do i1=1,Natomimp
!        IPOT=NSPIN*(I1-1)+1
!        write(756149,*) 'VINS',i1,IRMIND,IRMD,LMPOTD,IPOT,NSPIN
!        write(756149,*)
!    +             vinsimp(IRMIND:IRMD,1:LMPOTD,IPOT:IPOT+NSPIN-1)
!        write(756149,*) 'next...'
!        write(756149,*)
!        write(756149,*) 'VM2Z',i1,IRMD,IPOT,NSPIN
!        write(756149,*)
!    +          vm2zimp(1:IRMD,IPOT:IPOT+NSPIN-1)
!        write(756149,*) 'next...'
!        write(756149,*)
!        write(756149,*) 'VINSNEW',i1,IRMDNEW,LMPOTD,NSPIND
!        write(756149,*)
!    +          vinsnew(1:IRMDNEW,1:LMPOTD,IPOT:IPOT+NSPIN-1)
!        write(756149,*)
!      end do
!      close(756149)
 

       DO I1=i1_start_imp,i1_end_imp
        THETAimp = THETAhost(i1)
        PHIimp = PHIhost(i1)
        ISPIN=1
        IPOT=NSPIN*(I1-1)+1
        WRITE(6,*) 'IMP',I1,ipot

       ALLOCATE(VNSIMP(NSRA*LMMAXSO,NSRA*LMMAXSO,IRMDNEW))
       VNSIMP=CZERO
c set up the non-spherical ll' matrix for potential VLL'
       IF (NSRA.EQ.2) THEN
        USE_SRATRICK=1
       ELSEIF (NSRA.EQ.1) THEN
        USE_SRATRICK=0
       ENDIF
       ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW))
       VNSPLL0=CZERO
       ! output potential onto which SOC is added
       ALLOCATE(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW))
       VNSPLL1=CZERO

!      if(i1==1) open(66666,file='vllmat_input2',form='formatted')
!      write(66666,*) IRMDNEW,IRMDNEW,LMMAXD,LMMAXSO
!      write(66666,*) 'next...'
!      write(66666,*)
!      write(66666,*) VNSPLL0
!      write(66666,*) 'next...'
!      write(66666,*)
!      write(66666,*) VINSNEW(:,:,IPOT:IPOT+NSPIN-1)
!      write(66666,*) 'next...'
!      write(66666,*)
!      write(66666,*) ICLEB,IEND,NSPIN
!      write(66666,*) 'next...'
!      write(66666,*)
!      write(66666,*) ZIMP(I1),RNEW,USE_SRATRICK
!      write(66666,*) 'next...'
!      write(66666,*)

       CALL VLLMAT(1,IRMDNEW,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,
     +             VINSNEW(:,:,IPOT:IPOT+NSPIN-1),CLEB,ICLEB,IEND,NSPIN,
     +             ZIMP(I1),RNEW(:,I1),USE_SRATRICK)

!      if(i1==1) open(66667,file='spinorbitham_input2',form='formatted')
!      write(66667,*) LMAX,LMMAXD
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'VINSNEW',i1
!      write(66667,*) VINSNEW(:,:,IPOT:IPOT+NSPIN-1)
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*)  RNEW(:,I1)
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) E,ZIMP(I1),C,t_params%SOCSCALE(I1)
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(:,I1)
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) RPAN_INTERVALL(:,I1),NPAN_TOT(I1),NCHEB
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) NPAN_TOT(I1)*(NCHEB+1),IRMDNEW
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) 'vnspll0',i1
!      write(66667,*) VNSPLL0
!      write(66667,*) 'next...'
!      write(66667,*)
!      write(66667,*) VNSPLL1
!      write(66667,*) 'next...'
!      write(66667,*)
c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW(:,:,IPOT:IPOT+NSPIN-1),
     +                    RNEW(:,I1),E,ZIMP(I1),C,t_params%SOCSCALE(I1),
     +                    NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(:,I1),
     +                    RPAN_INTERVALL(:,I1),NPAN_TOT(I1),NCHEB,
     +                    NPAN_TOT(I1)*(NCHEB+1),IRMDNEW,VNSPLL0,
     +                    VNSPLL1,'1')
!      write(66667,*) 'vnspll1',i1
!      write(66667,*) VNSPLL1
!      write(66667,*) 'next...'
!      write(66667,*)
c extend matrix for the SRA treatment
       IF (NSRA.EQ.2) THEN
        ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW))
        IF (USE_SRATRICK.EQ.0) THEN
         CALL VLLMATSRA(VNSPLL1,VNSPLL,RNEW(:,I1),LMMAXSO,IRMDNEW,
     +                 IRMDNEW,E,C,LMAX,0,'Ref=0')                     
        ELSEIF (USE_SRATRICK.EQ.1) THEN                    
         CALL VLLMATSRA(VNSPLL1,VNSPLL,RNEW(:,I1),LMMAXSO,IRMDNEW,
     +                 IRMDNEW,E,C,LMAX,0,'Ref=Vsph')
        ENDIF
       ELSE
        ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW))
        VNSPLL(:,:,:)=VNSPLL1(:,:,:)
       ENDIF

c calculate the source terms in the Lippmann-Schwinger equation
c these are spherical hankel and bessel functions
       ALLOCATE(HLK(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(JLK(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(HLK2(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(JLK2(1:4*(LMAX+1),IRMDNEW))
       HLK=CZERO
       JLK=CZERO
       HLK2=CZERO
       JLK2=CZERO
       GMATPREFACTOR=CZERO
       CALL RLLSLLSOURCETERMS(NSRA,NVEC,E,RNEW(:,I1),IRMDNEW,IRMDNEW,
     +                        LMAX,LMMAXSO,1,JLK_INDEX,
     +                        HLK,JLK,HLK2,JLK2,GMATPREFACTOR)
c using spherical potential as reference
       IF (USE_SRATRICK.EQ.1) THEN
       CALL CALCSPH(NSRA,IRMDNEW,IRMDNEW,LMAX,NSPIN,ZIMP(I1),C,E,LMPOTD,
     +              LMMAXSO,RNEW(:,I1),VINSNEW(:,:,IPOT:IPOT+NSPIN-1)
     +              ,NCHEB,NPAN_TOT(I1),RPAN_INTERVALL(:,I1),
     +              JLK_INDEX,HLK,JLK,HLK2,JLK2,
     +              GMATPREFACTOR,TMATSPH,dummy_alpha,USE_SRATRICK)
       ENDIF

c calculate the tmat and wavefunctions
       ALLOCATE(RLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW))
       ALLOCATE(SLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW))
       RLL=CZERO
       SLL=CZERO

c right solutions
       CALL RLLSLL(RPAN_INTERVALL(:,I1),RNEW(:,I1),VNSPLL,RLL,
     +          SLL,TMATLLIMP(:,:,I1),NCHEB,
     +          NPAN_TOT(I1),LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +          IRMDNEW,NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,
     +          '1','1','0',USE_SRATRICK,dummy_alphaget)
       IF (NSRA.EQ.2) THEN
        DO IR=1,IRMDNEW
         DO LM1=1,LMMAXSO
          DO LM2=1,LMMAXSO
           RLL(LM1+LMMAXSO,LM2,IR)=
     +            RLL(LM1+LMMAXSO,LM2,IR)/C
           SLL(LM1+LMMAXSO,LM2,IR)=
     +            SLL(LM1+LMMAXSO,LM2,IR)/C
          ENDDO
         ENDDO
        ENDDO
       ENDIF
            
c add spherical contribution of tmatrix
       
       IF (USE_SRATRICK.EQ.1) THEN
        DO LM1=1,(KORBIT+1)*LMMAXD
         TMATLLIMP(LM1,LM1,I1)=TMATLLIMP(LM1,LM1,I1)+
     &                         TMATSPH(JLK_INDEX(LM1))
        ENDDO
       ENDIF

c rotate tmatrix and radial wavefunction to global frame

       CALL ROTATEMATRIX(TMATLLIMP(:,:,I1),THETA,PHI,LMMAXD,0)

c create SRA potential for impurity
c set up the non-spherical ll' matrix for potential VLL'
       VNSPLL0=CZERO
       CALL VLLMAT(1,IRMDNEW,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,
     +             VINSNEW(:,:,IPOT:IPOT+NSPIN-1),CLEB,ICLEB,IEND,NSPIN,
     +             ZIMP(I1),RNEW(:,I1),0)
!     +             ZIMP(I1),RNEW(:,I1),USE_SRATRICK)

c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW(:,:,IPOT:IPOT+NSPIN-1),
     +                    RNEW(:,I1),E,ZIMP(I1),C,t_params%SOCSCALE(I1),
     +                    NSPIN,LMPOTD,THETA,PHI,IPAN_INTERVALL(:,I1),
     +                    RPAN_INTERVALL(:,I1),NPAN_TOT(I1),NCHEB,
     +                    NPAN_TOT(I1)*(NCHEB+1),IRMDNEW,VNSPLL0,
     +                    VNSPLL1,'1')
       DO IR=1,IRMDNEW
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          VNSIMP(LM1,LM2,IR)=VNSPLL1(LM1,LM2,IR)
          IF (NSRA.EQ.2) THEN
           VNSIMP(LM1+LMMAXSO,LM2+LMMAXSO,IR)=VNSPLL1(LM1,LM2,IR)
          ENDIF
         ENDDO
        ENDDO
        !if(ir<20) write(789,*) vnsimp(:,:,ir)
       ENDDO

c calculate delta_t_imp matrix written in TMATLLIMP

        DO I2=1,IHOST
         IF (ATOMIMP(I1).EQ.HOSTIMP(I2)) THEN
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            TMATLLIMP(LM1,LM2,I1)=TMATLLIMP(LM1,LM2,I1)-
     +                            TMATLL(LM1,LM2,I2)
           ENDDO
          ENDDO
          DO LM1=1,NSRA*LMMAXSO
           DO LM2=1,NSRA*LMMAXSO
            DO IR=1,IRMDNEW
             VNSIMP(LM1,LM2,IR)=VNSIMP(LM1,LM2,IR)-
     +                             VNSHOST(LM1,LM2,I2,IR)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

c calculate delta matrix \delta=int{R_imp*\deltaV*R_host}
 
       IF (IELAST.EQ.1) THEN
        ALLOCATE(DELTABG(LMMAXSO,LMMAXSO,IRMDNEW))
        ALLOCATE(DELTASM(LMMAXSO,LMMAXSO,IRMDNEW))
        DELTABG=CZERO
        DELTASM=CZERO
        ALLOCATE(DELTATMP(IRMDNEW))
        ALLOCATE(RADIALHOST(LMMAXSO,LMMAXSO))
        ALLOCATE(RADIALIMP(LMMAXSO,LMMAXSO))
        ALLOCATE(VLLIMP(LMMAXSO,LMMAXSO))
        DELTATMP=CZERO

c big component for SRA stored in DELTABG
        DO IR=1,IRMDNEW
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
!        if(ir<20) then
!           if(ir==1) open(77776,file='deltabg0',form='formatted')
!           write(77776,*) 'radialhost'
!           write(77776,*) RADIALHOST
!           write(77776,*) 'radialIMP'
!           write(77776,*) RADIALIMP
!           write(77776,*) 'VLLIMP'
!           write(77776,*) VLLIMP
!        end if
         CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VLLIMP,
     &              LMMAXSO,RADIALIMP,LMMAXSO,CZERO,DELTAV,LMMAXSO)
c         CALL ZGEMM('T','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
c     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTABG(1,1,IR),LMMAXSO)
         CALL ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTABG(:,:,IR),LMMAXSO)

c small component for SRA stored in DELTASM
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
!       if(ir<20) then
!          if(ir==1) open(77775,file='deltabg1',form='formatted')
!          write(77775,*) 'radialhost'
!          write(77775,*) RADIALHOST
!          write(77775,*) 'radialIMP'
!          write(77775,*) RADIALIMP
!          write(77775,*) 'VLLIMP'
!          write(77775,*) VLLIMP
!       end if
         CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VLLIMP,
     &              LMMAXSO,RADIALIMP,LMMAXSO,CZERO,DELTAV,LMMAXSO)
c         CALL ZGEMM('T','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
c     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTASM(1,1,IR),LMMAXSO)
         CALL ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTASM(:,:,IR),LMMAXSO)

c sum up big and small component stored in DELTABG
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            DELTABG(LM1,LM2,IR)=DELTABG(LM1,LM2,IR)+DELTASM(LM1,LM2,IR)
           ENDDO
          ENDDO

         ENDIF ! NSRA
        ENDDO ! IR

!        open(77777,file='deltabg',form='formatted')
!        write(77777,*) DELTABG

c integrate 
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          DO IR=1,IRMDNEW
           DELTATMP(IR)=DELTABG(LM1,LM2,IR)
          ENDDO
          CALL INTCHEB_CELL(DELTATMP,DELTAMTR(LM1,LM2,I1),
     &                      RPAN_INTERVALL(:,I1),IPAN_INTERVALL(:,i1),
     +                      NPAN_TOT(I1),NCHEB,IRMDNEW)
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



! final writeout ony on master
#ifdef CPP_MPI
      !collect results and write out only on master
      !collect TMATLLIMP, DELTAMTR 
      !communicate TMATLLIMP, DELTAMTR 
      i1 = LMMAXSO*LMMAXSO*NATOMIMP
      allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=ierr)
      temp = CZERO
      if(ierr/=0) stop 'Error allocating temp for TMATLLIMP MPI comm.' 
      call MPI_ALLREDUCE(TMATLLIMP,temp,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                   MPI_COMM_WORLD,ierr)
      if(ierr/=0) stop 'Error in MPI_Allreduce for TMATLLIMP in tmatimp'
      TMATLLIMP = temp
      deallocate(temp)

      i1 = LMMAXSO*LMMAXSO*NATOMIMP
      allocate(temp(LMMAXSO,LMMAXSO,NATOMIMP),stat=ierr)
      temp = CZERO
      if(ierr/=0) stop 'Error allocating temp for DELTAMTR MPI comm.' 
      call MPI_ALLREDUCE(DELTAMTR,temp,i1,MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                   MPI_COMM_WORLD,ierr)
      if(ierr/=0) stop 'Error in MPI_Allreduce for DELTAMTR in tmatimp'
      DELTAMTR = temp
      deallocate(temp)
#endif

      if(myrank==master) then

         DO I1=1,NATOMIMP
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            IL1=LMMAXSO*(I1-1)+LM1
            IL2=LMMAXSO*(I1-1)+LM2
            DTMTRX(IL1,IL2)=TMATLLIMP(LM1,LM2,I1)
           ENDDO
          ENDDO
         ENDDO

c write down to the file DTMTRX
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
          DO LM1=1,LMMAXSO*NATOMIMP
           DO LM2=1,LMMAXSO*NATOMIMP
             WRITE(20,'((2I5),(4e17.9))') LM2,LM1,DTMTRX(LM2,LM1),
     &                                    DELTAIMP(LM2,LM1)
           ENDDO
          ENDDO
          DEALLOCATE(DELTAIMP)
          WRITE(6,*) 'end of delta t'
         ENDIF ! IELAST.EQ.1
         
         CLOSE(20) ! output file DTMTRX

       end if ! myrank==master  
         
         
       DEALLOCATE(VNSHOST)
       DEALLOCATE(RLLHOST)
       DEALLOCATE(TMATLLIMP)
       DEALLOCATE(DELTAMTR)
       DEALLOCATE(DELTAV)
       END 

