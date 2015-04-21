       SUBROUTINE RHOVALNEW(LDORHOEF,IELAST,NSRA,NSPIN,LMAX,EZ,WEZ,ZAT,
     +                      SOCSCALE,CLEB,ICLEB,IEND,IFUNM,LMSP,NCHEB,
     +                      NPAN_TOT,NPAN_LOG,NPAN_EQ,RMESH,IRWS,
     +                      RPAN_INTERVALL,IPAN_INTERVALL,
     +                      RNEW,VINSNEW,THETASNEW,THETA,PHI,I1,IPOT,
     +                      DEN_out,ESPV,RHO2NS,R2NEF,MUORB)
       use omp_lib
       use mod_types, only: t_tgmat
       IMPLICIT NONE
       include 'inc.p'
       INTEGER NSPIN,NSRA,LMAX,IEND,IPOT,IELAST,NPAN_TOT,NCHEB,
     +         NPAN_LOG,NPAN_EQ,IRWS
       INTEGER LMMAXD
       PARAMETER (LMMAXD= (LMAXD+1)**2)
       INTEGER LMAXD1
       PARAMETER (LMAXD1= LMAXD+1)
       INTEGER LMMAXSO
       PARAMETER (LMMAXSO=2*LMMAXD)
       INTEGER LMPOTD
       PARAMETER (LMPOTD= (LPOTD+1)**2)
       INTEGER LMXSPD
       PARAMETER (LMXSPD= (2*LPOTD+1)**2)
       DOUBLE PRECISION CVLIGHT
       PARAMETER (CVLIGHT=274.0720442D0)
       DOUBLE COMPLEX CZERO,CONE
       PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0))
       INTEGER NRMAXD
       PARAMETER (NRMAXD=NTOTD*(NCHEBD+1))
       DOUBLE PRECISION ZAT,RMESH(IRMD)
       DOUBLE COMPLEX EZ(IEMXD),ERYD,WEZ(IEMXD),EK,DF
       DOUBLE PRECISION SOCSCALE
       DOUBLE PRECISION CLEB(*)
       INTEGER ICLEB(NCLEB,4),IFUNM(LMXSPD),LMSP(LMXSPD)
       DOUBLE PRECISION  RNEW(NRMAXD),
     +                   RPAN_INTERVALL(0:NTOTD),
     +                   VINSNEW(NRMAXD,LMPOTD,NSPOTD),
     +                   THETASNEW(NRMAXD,NFUND)
       INTEGER           IPAN_INTERVALL(0:NTOTD)
       DOUBLE COMPLEX  
     +   TMATLL(LMMAXSO,LMMAXSO),
     +   TMATTEMP(LMMAXSO,LMMAXSO)
       DOUBLE COMPLEX GMATLL(LMMAXSO,LMMAXSO,IEMXD),
     +                GMAT0(LMMAXSO,LMMAXSO)
       INTEGER I1,IR,IREC,USE_SRATRICK,NVEC,LM1,LM2,IE,IRMDNEW,IMT1,
     +         JSPIN,IDIM,IORB,L1
       DOUBLE PRECISION THETA,PHI,PI,THETANEW,PHINEW
       DOUBLE COMPLEX GMATPREFACTOR
       DOUBLE PRECISION, ALLOCATABLE :: VINS(:,:,:)
       DOUBLE COMPLEX,ALLOCATABLE :: VNSPLL0(:,:,:),VNSPLL1(:,:,:,:),
     +                               VNSPLL(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: HLK(:,:,:),JLK(:,:,:),
     +                                HLK2(:,:,:),JLK2(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: RLL(:,:,:,:),SLL(:,:,:,:),
     +   RLLLEFT(:,:,:,:),SLLLEFT(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: TMATSPH(:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: CDEN(:,:,:,:),
     +   CDENLM(:,:,:,:),CDENNS(:,:,:),RHO2NSC(:,:,:),R2NEFC(:,:,:),
     +   RHO2NSNEW(:,:,:),R2NEFNEW(:,:,:),R2ORBC(:,:,:,:),
     +   GFLLE_PART(:,:,:),GFLLE(:,:,:,:),RHO2NSC_loop(:,:,:,:),
     +   R2NEFC_loop(:,:,:,:)
       DOUBLE PRECISION RHO2NS(IRMD,LMPOTD,4),R2NEF(IRMD,LMPOTD,4)
       DOUBLE COMPLEX, ALLOCATABLE:: DEN(:,:,:,:),DENLM(:,:,:,:)
       DOUBLE COMPLEX RHO2(4),RHO2INT(4),TEMP1,DEN_out(0:LMAXD1,IEMXD,2)
       DOUBLE PRECISION ESPV(0:LMAXD1,2)
       DOUBLE COMPLEX RHO2NS_TEMP(2,2),DENTEMP
       DOUBLE PRECISION MOMENT(3),TOTMOMENT,TOTXYMOMENT
       DOUBLE PRECISION MUORB(0:LMAXD1+1,3),
     +                  DENORBMOM(3),DENORBMOMSP(2,4),
     +                  DENORBMOMLM(0:LMAXD,3),DENORBMOMNS(3)
       DOUBLE COMPLEX, ALLOCATABLE :: CDENTEMP(:,:),
     +                                RHOTEMP(:,:),RHONEWTEMP(:,:)
       INTEGER JLK_INDEX(2*LMMAXSO)
       LOGICAL LDORHOEF
       LOGICAL TEST,OPT
       EXTERNAL TEST,OPT
       DOUBLE PRECISION QVEC(:,:)       ! qdos ruess: q-vectors for qdos
       ALLOCATABLE QVEC                 ! qdos ruess
       DOUBLE COMPLEX DENTOT(2),DENTMP(0:LMAXD1,2) ! qdos ruess
       INTEGER IQ,NQDOS ! qdos ruess: number of qdos points
       INTEGER IX,M1,LMSHIFT1(4),LMSHIFT2(4)       ! qdos ruess
       INTEGER LRECGFLLE,IERR                           ! lmlm-dos
       ! OMP - number of threads, thread id
       integer nth,ith
! determine if omp is used
!$omp parallel shared(nth,ith) 
!$omp single 
      ith = 0
      nth = 1
#ifdef CPP_OMP
      nth = omp_get_num_threads()
#endif
!$omp end single
!$omp end parallel
       write(*,*) 'nth =',nth

       PI=4d0*DATAN(1d0)
       IRMDNEW= NPAN_TOT*(NCHEB+1)
       IMT1=IPAN_INTERVALL(NPAN_LOG+NPAN_EQ)+1
       ALLOCATE(VINS(IRMDNEW,LMPOTD,NSPIN))
       VINS=0d0
       DO LM1=1,LMPOTD
        DO IR=1,IRMDNEW
         VINS(IR,LM1,1)=VINSNEW(IR,LM1,IPOT)
         VINS(IR,LM1,NSPIN)=VINSNEW(IR,LM1,IPOT+NSPIN-1)
        ENDDO
       ENDDO

cc set up the non-spherical ll' matrix for potential VLL'
       USE_SRATRICK=1
       ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW))
       ALLOCATE(VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1))
       VNSPLL0=CZERO
       CALL VLLMAT(1,NRMAXD,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINS,
     +                    CLEB,ICLEB,IEND,NSPIN,ZAT,RNEW,USE_SRATRICK)

c initial allocate
       IF (NSRA.EQ.2) THEN
        ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW,0:nth-1))
       ELSE
        ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1))
       ENDIF

       ALLOCATE(HLK(4*(LMAX+1),IRMDNEW,0:nth-1))
       ALLOCATE(JLK(4*(LMAX+1),IRMDNEW,0:nth-1))
       ALLOCATE(HLK2(4*(LMAX+1),IRMDNEW,0:nth-1))
       ALLOCATE(JLK2(4*(LMAX+1),IRMDNEW,0:nth-1))
       ALLOCATE(TMATSPH(2*(LMAX+1),0:nth-1))
       ALLOCATE(RLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1))
       ALLOCATE(SLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1))
       ALLOCATE(RLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1))
       ALLOCATE(SLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1))
       ALLOCATE(CDEN(IRMDNEW,0:LMAXD,4,0:nth-1))
       ALLOCATE(CDENLM(IRMDNEW,LMMAXD,4,0:nth-1))
       ALLOCATE(CDENNS(IRMDNEW,4,0:nth-1))
       ALLOCATE(RHO2NSC(IRMDNEW,LMPOTD,4))
       ALLOCATE(RHO2NSC_loop(IRMDNEW,LMPOTD,4,ielast))
       ALLOCATE(RHO2NSNEW(IRMD,LMPOTD,4))
       ALLOCATE(R2NEFC(IRMDNEW,LMPOTD,4))
       ALLOCATE(R2NEFC_loop(IRMDNEW,LMPOTD,4,0:nth-1))
       ALLOCATE(R2NEFNEW(IRMD,LMPOTD,4))
       ALLOCATE(R2ORBC(IRMDNEW,LMPOTD,4,0:nth-1))
       ALLOCATE(CDENTEMP(IRMDNEW,0:nth-1))
       ALLOCATE(GFLLE_PART(LMMAXSO,LMMAXSO,0:nth-1))
       ALLOCATE(GFLLE(LMMAXSO,LMMAXSO,IELAST,1))
       ALLOCATE(DEN(0:LMAXD1,IEMXD,2,1),DENLM(LMMAXD,IEMXD,2,1))
       RHO2NSC=CZERO
       RHO2NSC_loop=CZERO
       R2NEFC=CZERO
       R2NEFC_loop=CZERO
       R2ORBC=CZERO
       RHO2NS=0.D0  ! fivos 19.7.2014, this was CZERO
       R2NEF=0.D0   ! fivos 19.7.2014, this was CZERO
       RHO2NSNEW=CZERO
       R2NEFNEW=CZERO
       DEN=CZERO
       ESPV=0d0
       RHO2INT(:)=CZERO
       DENORBMOM=0d0
       DENORBMOMSP=0d0
       DENORBMOMLM=0d0
       DENORBMOMNS=0d0
       THETANEW=0d0
       PHINEW=0d0
       GFLLE_PART=CZERO
       GFLLE=CZERO
c LM shifts for correct density summation
       LMSHIFT1(1)=0                                                   ! qdos ruess
       LMSHIFT1(2)=LMMAXD                                              ! qdos ruess
       LMSHIFT1(3)=0                                                   ! qdos ruess
       LMSHIFT1(4)=LMMAXD                                              ! qdos ruess
       LMSHIFT2(1)=0                                                   ! qdos ruess
       LMSHIFT2(2)=LMMAXD                                              ! qdos ruess
       LMSHIFT2(3)=LMMAXD                                              ! qdos ruess
       LMSHIFT2(4)=0                                                   ! qdos ruess

       DO IR=1,3
        DO LM1=0,LMAXD1+1
         MUORB(LM1,IR)=0d0
        ENDDO
       ENDDO 

      NQDOS = 1                                                         ! qdos ruess 
      IF (OPT('qdos    ')) THEN                                         ! qdos ruess
C        Read BZ path for qdos calculation:                             ! qdos ruess
         OPEN(67,FILE='qvec.dat',STATUS='old',IOSTAT=IERR,ERR=3000)     ! qdos ruess
         READ(67,*) NQDOS                                               ! qdos ruess
         ALLOCATE(QVEC(3,NQDOS))                                        ! qdos ruess
         DO IQ = 1,NQDOS                                                ! qdos ruess
            READ(67,*) (QVEC(IX,IQ),IX=1,3)                             ! qdos ruess
         ENDDO                                                          ! qdos ruess
         CLOSE(67)                                                      ! qdos ruess
C        Change allocation for GFLLE to be suitabel for qdos run        ! qdos ruess
         DEALLOCATE(GFLLE,DEN,DENLM)                                    ! qdos ruess
         ALLOCATE(GFLLE(LMMAXSO,LMMAXSO,IELAST,NQDOS))                  ! qdos ruess
         ALLOCATE(DEN(0:LMAXD1,IEMXD,2,NQDOS),
     +            DENLM(LMMAXD,IEMXD,2,NQDOS))
3000  IF (IERR.NE.0) STOP 'ERROR READING ''qvec.dat'''                  ! qdos ruess
      END IF  ! OPT('qdos    ')                                         ! qdos ruess

      IF ((OPT('lmlm-dos')).AND.(I1.EQ.1)) THEN                         ! lmlm-dos ruess
         LRECGFLLE = 4*LMMAXSO*LMMAXSO*IELAST*NQDOS                     ! lmlm-dos ruess
         OPEN(91,ACCESS='direct',RECL=LRECGFLLE,FILE='gflle',           ! lmlm-dos ruess
     &         FORM='unformatted',STATUS='replace',ERR=3001,IOSTAT=IERR)! lmlm-dos ruess
 3001 IF (IERR.NE.0) STOP 'ERROR CREATING ''gflle'''                    ! lmlm-dos ruess
      ENDIF                                                             ! lmlm-dos ruess

c energy loop
       WRITE(6,*) 'atom: ',I1
! omp: start parallel region here
!$omp parallel do default(none)
!$omp& private(eryd,ie,ir,irec,lm1,lm2,gmatprefactor,nvec)
!$omp& private(jlk_index,tmatll,ith)
!$omp& shared(nspin,nsra,iend,ipot,ielast,npan_tot,ncheb,lmax)
!$omp& shared(zat,socscale,ez,rmesh,cleb,rnew,nth,icleb,thetasnew,i1)
!$omp& shared(rpan_intervall,vinsnew,ipan_intervall,r2nefc_loop)
!$omp& shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)
!$omp& shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,cdentemp)
!$omp& shared(tmatsph,den,denlm,gflle,gflle_part,rllleft,sllleft)
!$omp& private(iq,df,ek,tmattemp,gmatll,gmat0,iorb,dentemp)
!$omp& private(rho2ns_temp,dentot,dentmp,rho2,temp1,jspin)
!$omp& shared(ldorhoef,nqdos,lmshift1,lmshift2,wez,lmsp,imt1,ifunm)
!$omp& shared(r2orbc,r2nefc,cden,cdenlm,cdenns,rho2nsc_loop)
!$omp& reduction(+:rho2int,espv) reduction(-:muorb) 
!$omp& reduction(-:denorbmom,denorbmomsp,denorbmomlm,denorbmomns)
!$omp& shared(t_tgmat)
       DO IE=1,IELAST
         if (nth>=1) then
#ifdef CPP_OMP
            ith = omp_get_thread_num()
#endif
         else
            ith = 0
         endif

        ERYD=EZ(IE)
        EK=SQRT(ERYD)
        DF=WEZ(IE)/DBLE(NSPIN)
        IF (NSRA.EQ.2) EK = SQRT( ERYD + ERYD*ERYD/(CVLIGHT*CVLIGHT) ) *
     &                          ( 1d0 + ERYD/(CVLIGHT*CVLIGHT) )
!$omp critical
        WRITE(6,*) 'energy:',IE,'',ERYD
!$omp end critical
! 
!         IREC=IE+IELAST*(I1-1)
!         READ(69,REC=IREC) GMAT0
! 
! c rotate gmat from global frame to local frame
!         CALL ROTATEMATRIX(GMAT0,THETA,PHI,LMMAXD,1)
! 
!         DO LM1=1,LMMAXSO
!          DO LM2=1,LMMAXSO
!           GMATLL(LM1,LM2,IE)=GMAT0(LM1,LM2)
!          ENDDO
!         ENDDO

c recalculate wavefuntions, also include left solution
c contruct the spin-orbit coupling hamiltonian and add to potential
        CALL SPINORBIT_HAM(LMAX,LMMAXD,VINS,RNEW,
     +                     ERYD,ZAT,CVLIGHT,SOCSCALE,NSRA,NSPIN,LMPOTD,
     +                     THETA,PHI,IPAN_INTERVALL,RPAN_INTERVALL,
     +                     NPAN_TOT,NCHEB,IRMDNEW,NRMAXD,
     &                     VNSPLL0,VNSPLL1(:,:,:,ith),'1')

cc extend matrix for the SRA treatment
        VNSPLL(:,:,:,ith)=CZERO
        IF (NSRA.EQ.2) THEN
         IF (USE_SRATRICK.EQ.0) THEN
          CALL VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,
     +               LMMAXSO,IRMDNEW,NRMAXD,ERYD,CVLIGHT,LMAX,0,'Ref=0')
         ELSEIF (USE_SRATRICK.EQ.1) THEN
          CALL VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,
     +            LMMAXSO,IRMDNEW,NRMAXD,ERYD,CVLIGHT,LMAX,0,'Ref=Vsph')
         ENDIF
        ELSE
         VNSPLL(:,:,:,ith)=VNSPLL1(:,:,:,ith)
        ENDIF

cc calculate the source terms in the Lippmann-Schwinger equation
cc these are spherical hankel and bessel functions
        HLK(:,:,ith)=CZERO
        JLK(:,:,ith)=CZERO
        HLK2(:,:,ith)=CZERO
        JLK2(:,:,ith)=CZERO
        GMATPREFACTOR=CZERO
        JLK_INDEX=0
        CALL RLLSLLSOURCETERMS(NSRA,NVEC,ERYD,RNEW,IRMDNEW,NRMAXD,LMAX,
     +                         LMMAXSO,1,JLK_INDEX,HLK(:,:,ith),
     +                         JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),
     +                         GMATPREFACTOR)

c using spherical potential as reference
        IF (USE_SRATRICK.EQ.1) THEN
         CALL CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,CVLIGHT,ERYD,
     +           LMPOTD,LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +                JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),HLK2(:,:,ith),
     +                JLK2(:,:,ith),GMATPREFACTOR,TMATSPH(:,ith),
     +                USE_SRATRICK)
        ENDIF

cc calculate the tmat and wavefunctions
        RLLLEFT(:,:,:,ith)=CZERO
        SLLLEFT(:,:,:,ith)=CZERO

cc right solutions
        TMATLL=CZERO
        CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),
     +              RLL(:,:,:,ith),SLL(:,:,:,ith),TMATLL,
     +              NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),
     +          IRMDNEW,NRMAXD,NSRA,JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),
     +              HLK2(:,:,ith),JLK2(:,:,ith),
     +              GMATPREFACTOR,'1','1','0',USE_SRATRICK)
        IF (NSRA.EQ.2) THEN
         RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
     +            RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
         SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
     +            SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
        ENDIF

c left solutions
c contruct the TRANSPOSE spin-orbit coupling hamiltonian and add to potential
        CALL SPINORBIT_HAM(LMAX,LMMAXD,VINS,RNEW,ERYD,ZAT,
     +                     CVLIGHT,SOCSCALE,NSRA,NSPIN,LMPOTD,THETA,PHI,
     +                     IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,
     +                     IRMDNEW,NRMAXD,VNSPLL0,VNSPLL1(:,:,:,ith),
     +                     'transpose')

cc extend matrix for the SRA treatment
        VNSPLL(:,:,:,ith)=CZERO
        IF (NSRA.EQ.2) THEN
         IF (USE_SRATRICK.EQ.0) THEN
          CALL VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,
     +               LMMAXSO,IRMDNEW,NRMAXD,ERYD,CVLIGHT,LMAX,0,'Ref=0')
         ELSEIF (USE_SRATRICK.EQ.1) THEN
          CALL VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,
     +            LMMAXSO,IRMDNEW,NRMAXD,ERYD,CVLIGHT,LMAX,0,'Ref=Vsph')
         ENDIF
        ELSE
         VNSPLL(:,:,:,ith)=VNSPLL1(:,:,:,ith)
        ENDIF

cc calculate the source terms in the Lippmann-Schwinger equation
cc these are spherical hankel and bessel functions
        HLK(:,:,ith)=CZERO
        JLK(:,:,ith)=CZERO
        HLK2(:,:,ith)=CZERO
        JLK2(:,:,ith)=CZERO
        GMATPREFACTOR=CZERO
        JLK_INDEX=0
        CALL RLLSLLSOURCETERMS(NSRA,NVEC,ERYD,RNEW,IRMDNEW,NRMAXD,LMAX,
     +                         LMMAXSO,1,JLK_INDEX,HLK(:,:,ith),
     +                         JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),
     +                         GMATPREFACTOR)

cc using spherical potential as reference
c notice that exchange the order of left and right hankel/bessel functions
        IF (USE_SRATRICK.EQ.1) THEN
         CALL CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,CVLIGHT,ERYD,
     +           LMPOTD,LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +                JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),
     +                HLK(:,:,ith),JLK(:,:,ith),GMATPREFACTOR,
     +                TMATSPH(:,ith),USE_SRATRICK)
        ENDIF

cc calculate the tmat and wavefunctions
        RLLLEFT(:,:,:,ith)=CZERO
        SLLLEFT(:,:,:,ith)=CZERO

cc left solutions
c notice that exchange the order of left and right hankel/bessel functions
        TMATTEMP=CZERO
        CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),
     +              RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),TMATTEMP,
     +              NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),
     +        IRMDNEW,NRMAXD,NSRA,JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),
     +              HLK(:,:,ith),JLK(:,:,ith),
     +              GMATPREFACTOR,'1','1','0',USE_SRATRICK)
        IF (NSRA.EQ.2) THEN
         RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
     +            RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
         SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
     +            SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
        ENDIF

        DO 200 IQ = 1,NQDOS                                       ! qdos
        DEN(:,ie,:,IQ)=CZERO
c read in gf
        IREC = IQ + NQDOS * (IE-1) +  NQDOS * IELAST * (I1-1)     ! qdos
        !$omp critical
        if (t_tgmat%gmat_to_file) then
         READ(69,REC=IREC) GMAT0 
        else
         GMAT0(:,:) = t_tgmat%gmat(:,:,irec)
        end if
        !$omp end critical

c rotate gmat from global frame to local frame
        CALL ROTATEMATRIX(GMAT0,THETA,PHI,LMMAXD,1)

        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          GMATLL(LM1,LM2,IE)=GMAT0(LM1,LM2)
         ENDDO
        ENDDO
c calculate density 

! !$omp critical
!        open(555555,file="rhoout_input_ith"//char(48+ith)//'.ie'
!      +//char(48+ie))
!        write(555555,*) 'line1',NSRA,LMMAXD,LMMAXSO,LMAX,GMATLL(:,:,IE)
!      +     ,'line2',EK,LMPOTD,DF,NPAN_TOT,NCHEB,ICLEB,IEND,
!      +      'line3',IRMDNEW,NRMAXD,THETASNEW,IFUNM,RNEW,IMT1,LMSP,
!      +      'line4',RLL(:,:,:,ith),SLL(:,:,:,ith),
!      +      'line5',RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),
!      +      'line6',CDEN(:,:,:,ith),CDENLM(:,:,:,ith),
!      +      'line7',CDENNS(:,:,ith),RHO2NSC_loop(:,:,:,ie),0,
!      +      'line8',GFLLE(:,:,IE,IQ),RPAN_INTERVALL,IPAN_INTERVALL
!        close(555555)
! !$omp end critical

       CALL RHOOUTNEW(NSRA,LMMAXD,LMMAXSO,LMAX,GMATLL(1,1,IE),EK,
     +                LMPOTD,DF,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,
     +                IRMDNEW,NRMAXD,THETASNEW,IFUNM,RNEW,IMT1,LMSP,
     +                RLL(:,:,:,ith),SLL(:,:,:,ith),
     +                RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),
     +                CDEN(:,:,:,ith),CDENLM(:,:,:,ith),
     +                CDENNS(:,:,ith),RHO2NSC_loop(:,:,:,ie),0,
     +                GFLLE(:,:,IE,IQ),RPAN_INTERVALL,IPAN_INTERVALL)

       DO JSPIN=1,4

        DO LM1 = 0,LMAX
         CDENTEMP(:,ith)=CZERO
         DENTEMP=CZERO
         DO IR=1,IRMDNEW
          CDENTEMP(IR,ith)=CDEN(IR,LM1,JSPIN,ith)
         ENDDO
         CALL INTCHEB_CELL(CDENTEMP(:,ith),DENTEMP,
     +             RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
         RHO2(JSPIN)=DENTEMP
         RHO2INT(JSPIN)=RHO2INT(JSPIN)+RHO2(JSPIN)*DF
         IF (JSPIN.LE.2) THEN
          DEN(LM1,IE,JSPIN,IQ)=RHO2(JSPIN)
         ENDIF
        ENDDO
                 
        IF (JSPIN.LE.2) THEN
         DO LM1 = 1,LMMAXD
          CDENTEMP(:,ith)=CZERO
          DENTEMP=CZERO
          DO IR=1,IRMDNEW
           CDENTEMP(IR,ith)=CDENLM(IR,LM1,JSPIN,ith)
          ENDDO
          CALL INTCHEB_CELL(CDENTEMP(:,ith),DENTEMP,
     +            RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
          DENLM(LM1,IE,JSPIN,IQ)=DENTEMP
         ENDDO
          CDENTEMP(:,ith)=CZERO
          DENTEMP=CZERO
          DO IR=1,IRMDNEW
           CDENTEMP(IR,ith)=CDENNS(IR,JSPIN,ith)
          ENDDO
          CALL INTCHEB_CELL(CDENTEMP(:,ith),DENTEMP,
     +           RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
          DEN(LMAXD1,IE,JSPIN,IQ)=DENTEMP
          RHO2INT(JSPIN)=RHO2INT(JSPIN)+DEN(LMAXD1,IE,JSPIN,IQ)*DF
         ENDIF
        ENDDO ! JSPIN

        DO JSPIN=1,4
         IF (JSPIN.LE.2) THEN
          DO LM1=0,LMAXD1
           ESPV(LM1,JSPIN)=ESPV(LM1,JSPIN)+
     +                       DIMAG( ERYD * DEN(LM1,IE,JSPIN,IQ) * DF )
          ENDDO
         ENDIF 
        ENDDO
 200   END DO   ! IQ = 1,NQDOS

c get charge at the Fermi energy (IELAST)

       IF (IE.EQ.IELAST.AND.LDORHOEF) THEN
       CALL RHOOUTNEW(NSRA,LMMAXD,LMMAXSO,LMAX,GMATLL(1,1,IE),EK,
     +                LMPOTD,CONE,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,
     +                IRMDNEW,NRMAXD,THETASNEW,IFUNM,RNEW,IMT1,LMSP,
     +                RLL(:,:,:,ith),SLL(:,:,:,ith),
     +                RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),
     +                CDEN(:,:,:,ith),CDENLM(:,:,:,ith),
     +                CDENNS(:,:,ith),R2NEFC_loop(:,:,:,ith),0,
     +                GFLLE_PART(:,:,ith),RPAN_INTERVALL,IPAN_INTERVALL)
       ENDIF


c get orbital moment
       DO IORB=1,3
        CALL RHOOUTNEW(NSRA,LMMAXD,LMMAXSO,LMAX,GMATLL(1,1,IE),EK,
     +                 LMPOTD,CONE,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,
     +                 IRMDNEW,NRMAXD,THETASNEW,IFUNM,RNEW,IMT1,LMSP,
     +                 RLL(:,:,:,ith),SLL(:,:,:,ith),
     +                 RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),
     +                 CDEN(:,:,:,ith),CDENLM(:,:,:,ith),
     +                 CDENNS(:,:,ith),R2ORBC(:,:,:,ith),IORB,
     +                GFLLE_PART(:,:,ith),RPAN_INTERVALL,IPAN_INTERVALL)
         DO JSPIN=1,4
          IF (JSPIN.LE.2) THEN
           DO LM1=0,LMAX
            CDENTEMP(:,ith)=CZERO
            DENTEMP=CZERO
            DO IR=1,IRMDNEW
             CDENTEMP(IR,ith)=CDEN(IR,LM1,JSPIN,ith)
            ENDDO 
            CALL INTCHEB_CELL(CDENTEMP(:,ith),DENTEMP,RPAN_INTERVALL,
     +                        IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
            RHO2(JSPIN)=DENTEMP
            MUORB(LM1,JSPIN)=MUORB(LM1,JSPIN)-DIMAG(RHO2(JSPIN)*DF)
            DENORBMOM(IORB)=DENORBMOM(IORB)-DIMAG(RHO2(JSPIN)*DF)
            DENORBMOMSP(JSPIN,IORB)=DENORBMOMSP(JSPIN,IORB)-
     +                   DIMAG(RHO2(JSPIN)*DF)
            DENORBMOMLM(LM1,IORB)=DENORBMOMLM(LM1,IORB)-
     +                   DIMAG(RHO2(JSPIN)*DF)
            CDENTEMP(:,ith)=CZERO
            DO IR=1,IRMDNEW
             CDENTEMP(IR,ith)=CDENNS(IR,JSPIN,ith)
            ENDDO
            CALL INTCHEB_CELL(CDENTEMP(:,ith),TEMP1,
     +            RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
            DENORBMOMNS(IORB)=DENORBMOMNS(IORB)-DIMAG(TEMP1*DF)
           ENDDO
          ENDIF
         ENDDO
        ENDDO ! IORB
       ENDDO ! IE loop 
!$omp end parallel do
! omp: move sum from rhooutnew here after parallel calculation
      DO IR=1,IRMDNEW
        DO LM1=1,LMPOTD
          DO JSPIN=1,4
            DO IE=1,IELAST
             RHO2NSC(IR,LM1,JSPIN) = RHO2NSC(IR,LM1,JSPIN) + 
     +                                    RHO2NSC_loop(IR,LM1,JSPIN,IE)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
! omp: don't forget to do the same with density at fermi energy:
      DO ith=0,nth-1
         R2NEFC(:,:,:) = R2NEFC(:,:,:) + R2NEFC_loop(:,:,:,ith)
      ENDDO

! omp: moved write-out of dos files out of parallel energy loop
c Write out qdos and lm-dos:                                            ! lm-dos
       DO IE=1,IELAST                                                   ! lm-dos
       DO IQ=1,NQDOS                                                    ! lm-dos
       IF ((IQ.EQ.1).AND.(IE.EQ.1)) THEN                                ! lm-dos
           OPEN(29,                                                     ! lm-dos
     &     FILE="lmdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//    ! lm-dos
     &     char(48+1)//".dat")                                          ! lm-dos
           WRITE (29,*) ' '                                             ! lm-dos
           WRITE (29,8600) '# ISPIN=',1,' I1=',I1                       ! lm-dos
           OPEN(30,                                                     ! lm-dos
     &     FILE="lmdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//    ! lm-dos
     &     char(48+2)//".dat")                                          ! lm-dos
           WRITE (30,*) ' '                                             ! lm-dos
           WRITE (30,8600) '# ISPIN=',2,' I1=',I1                       ! lm-dos
       ENDIF                                                            ! lm-dos

       IF (OPT('qdos    ')) THEN                                        ! qdos ruess
         IF ((IQ.EQ.1).AND.(IE.EQ.1)) THEN                              ! qdos ruess
            OPEN(31,                                                    ! qdos ruess
     +        FILE="qdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//  ! qdos ruess
     +        char(48+1)//".dat")                                       ! qdos ruess
            WRITE (31,*) ' '                                            ! qdos ruess
            WRITE (31,8600) '# ISPIN=',1,' I1=',I1                      ! qdos ruess
            WRITE(31,'(7(A,3X))') '#   Re(E)','Im(E)','k_x','k_y','k_z',! qdos
     &                      'DEN_tot','DEN_s,p,...'                     ! qdos
            OPEN(32,                                                    ! qdos ruess
     +        FILE="qdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//  ! qdos ruess
     +        char(48+2)//".dat")                                       ! qdos ruess
            WRITE (32,*) ' '                                            ! qdos ruess
            WRITE (32,8600) '# ISPIN=',2,' I1=',I1                      ! qdos ruess
            WRITE(32,'(7A)') '#   Re(E)','Im(E)','k_x','k_y','k_z',     ! qdos
     &                      'DEN_tot','DEN_s,p,...'                     ! qdos

 8600 FORMAT (a8,I3,a4,I5)                                              ! qdos ruess
         ENDIF   ! IQ.EQ.1                                              ! qdos ruess
            DO JSPIN =1,2                                               ! qdos ruess
              DENTOT(JSPIN) = DCMPLX(0.D0,0.D0)                         ! qdos ruess
              DO L1 = 0,LMAXD1                                          ! qdos ruess
                DENTOT(JSPIN) = DENTOT(JSPIN) + DEN(L1,IE,JSPIN,IQ)     ! qdos ruess
              ENDDO                                                     ! qdos ruess
            ENDDO                                                       ! qdos ruess
c    write qdos.nn.s.dat                                                ! qdos ruess
c    and lmdos.nn.s.dat                                                 ! qdos ruess
            WRITE(29,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),     ! qdos ruess
     &               (-DIMAG(DENLM(L1,IE,1,IQ))/PI,L1=1,LMMAXD)         ! qdos ruess
            WRITE(30,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),     ! qdos ruess
     &               (-DIMAG(DENLM(L1,IE,2,IQ))/PI,L1=1,LMMAXD)         ! qdos ruess
            WRITE(31,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),     ! qdos ruess
     &     -DIMAG(DENTOT(1))/PI,(-DIMAG(DEN(L1,IE,1,IQ))/PI,L1=0,LMAXD1)! qdos ruess
            WRITE(32,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),     ! qdos ruess
     &     -DIMAG(DENTOT(2))/PI,(-DIMAG(DEN(L1,IE,2,IQ))/PI,L1=0,LMAXD1)! qdos ruess
       ELSE                                                             ! lm-dos
            WRITE(29,9001) EZ(IE),                                      ! lm-dos
     &               (-DIMAG(DENLM(L1,IE,1,IQ))/PI,L1=1,LMMAXD)         ! lm-dos
            WRITE(30,9001) EZ(IE),                                      ! lm-dos
     &               (-DIMAG(DENLM(L1,IE,2,IQ))/PI,L1=1,LMMAXD)         ! lm-dos
 9001       FORMAT(30E12.4)                                             ! lm-dos
       ENDIF      ! OPT('qdos    ')                                     ! qdos ruess
 9000      FORMAT(5F10.6,40E16.8)                                       ! qdos ruess
       ENDDO !IQ
       ENDDO !IE

c write
       IF (OPT('lmlm-dos')) THEN                                         ! lmlm-dos
!          DO JSPIN = 1,2                                                  ! lmlm-dos
!           OPEN(90,                                                       ! lmlm-dos
!      &    FILE="lmlmdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//    ! lmlm-dos
!      &                                          char(48+JSPIN)//".dat")  ! lmlm-dos
!           DO IE = 1,IELAST                                               ! lmlm-dos
!              DO LM1 = 1,LMMAXD                                           ! lmlm-dos
!                 IF (.NOT.(OPT('qdos    '))) THEN                                   ! qdos
!                    WRITE(90,1000) EZ(IE),                                ! lmlm-dos
!      &                            (-DIMAG(GFLLE(LM1+LMSHIFT1(JSPIN),     ! lmlm-dos
!      &                     LM2+LMSHIFT2(JSPIN),IE,1))/PI,LM2 = 1,LMMAXD) ! lmlm-dos
!                 ELSE                                                               ! qdos
!                   DO IQ=1,NQDOS                                                    ! qdos
!                    WRITE(90,1000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),                    ! qdos
!      &                            QVEC(3,IQ),(-DIMAG(GFLLE(LM1+                    ! qdos
!      &                            LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),             ! qdos
!      &                            IE,IQ))/PI,LM2 = 1,LMMAXD)                       ! qdos
!                   ENDDO ! IQ=1,NQDOS                                               ! qdos
!                 ENDIF                                                              ! qdos
!              ENDDO                                                       ! lmlm-dos
!           ENDDO !IE                                                      ! lmlm-dos
!           CLOSE(90)                                                      ! lmlm-dos
!          ENDDO !JSPIN                                                    ! lmlm-dos
!  1000  FORMAT(5F10.6,I3,40E16.8)                                         ! lmlm-dos
c write gflle to file                                                    ! lmlm-dos
       WRITE(91,REC=I1) GFLLE                                            ! lmlm-dos
       ENDIF                                                             ! lmlm-dos

       ALLOCATE(RHOTEMP(IRMDNEW,LMPOTD))
       ALLOCATE(RHONEWTEMP(IRWS,LMPOTD))
       DO JSPIN=1,4
        RHOTEMP=CZERO
        RHONEWTEMP=CZERO
        DO LM1=1,LMPOTD
         DO IR=1,IRMDNEW
          RHOTEMP(IR,LM1)=RHO2NSC(IR,LM1,JSPIN)  
         ENDDO
        ENDDO
        CALL CHEB2OLDGRID(IRWS,IRMDNEW,LMPOTD,RMESH,NCHEB,NPAN_TOT,
     +                    RPAN_INTERVALL,IPAN_INTERVALL,
     +                    RHOTEMP,RHONEWTEMP,IRMD)
        DO LM1=1,LMPOTD
         DO IR=1,IRWS
          RHO2NSNEW(IR,LM1,JSPIN)=RHONEWTEMP(IR,LM1)        
         ENDDO
        ENDDO 
         
        RHOTEMP=CZERO
        RHONEWTEMP=CZERO
        DO LM1=1,LMPOTD
         DO IR=1,IRMDNEW
          RHOTEMP(IR,LM1)=R2NEFC(IR,LM1,JSPIN)        
         ENDDO
        ENDDO
        CALL CHEB2OLDGRID(IRWS,IRMDNEW,LMPOTD,RMESH,NCHEB,NPAN_TOT,
     +                     RPAN_INTERVALL,IPAN_INTERVALL,
     +                     RHOTEMP,RHONEWTEMP,IRMD)
        DO LM1=1,LMPOTD
         DO IR=1,IRWS
          R2NEFNEW(IR,LM1,JSPIN)=RHONEWTEMP(IR,LM1)        
         ENDDO
        ENDDO          
       ENDDO
       DEALLOCATE(RHOTEMP)
       DEALLOCATE(RHONEWTEMP)
c calculate new THETA and PHI for non-colinear 
       IF (.NOT.TEST('FIXMOM  ')) THEN
        RHO2NS_TEMP(1,1)=RHO2INT(1)
        RHO2NS_TEMP(2,2)=RHO2INT(2)
        RHO2NS_TEMP(1,2)=RHO2INT(3)
        RHO2NS_TEMP(2,1)=RHO2INT(4)

        CALL ROTATEMATRIX(RHO2NS_TEMP,THETA,PHI,1,0)

        RHO2INT(1)=RHO2NS_TEMP(1,1)
        RHO2INT(2)=RHO2NS_TEMP(2,2)
        RHO2INT(3)=RHO2NS_TEMP(1,2)
        RHO2INT(4)=RHO2NS_TEMP(2,1)
        
        
        MOMENT(1)=DIMAG(RHO2INT(3)+RHO2INT(4))
        MOMENT(2)=-REAL(RHO2INT(3)-RHO2INT(4))
        MOMENT(3)=DIMAG(-RHO2INT(1)+RHO2INT(2))
      
        TOTMOMENT=SQRT(MOMENT(1)**2+MOMENT(2)**2+MOMENT(3)**2)
        TOTXYMOMENT=SQRT(MOMENT(1)**2+MOMENT(2)**2)

        IF (ABS(TOTXYMOMENT).GT.1d-05) THEN
         IF (ABS(MOMENT(3).LT.1d-05)) THEN
          THETANEW=PI/2d0
         ELSE
          THETANEW=ACOS(MOMENT(3)/TOTMOMENT)
         ENDIF
         IF (TOTXYMOMENT.LT.1d-05) THEN
          PHINEW=0d0
         ELSE
          PHINEW=DATAN2(MOMENT(2),MOMENT(1))
         ENDIF
        ENDIF
c          THETANEW=ACOS(MOMENT(3)/TOTMOMENT)
c          PHINEW=DATAN2(MOMENT(2),MOMENT(1))
        WRITE(6,*) 'moment',MOMENT(1),MOMENT(2),MOMENT(3)
c        WRITE(6,*) 'total moment',TOTMOMENT,TOTXYMOMENT
        WRITE(6,*) THETANEW,PHINEW
        WRITE(11,*) THETANEW,PHINEW
        WRITE(12,*) THETANEW,PHINEW
        CALL ROTATEVECTOR(RHO2NSNEW,RHO2NS,IRWS,LMPOTD,THETANEW,PHINEW,
     +                    THETA,PHI,IRMD)
        CALL ROTATEVECTOR(R2NEFNEW,R2NEF,IRWS,LMPOTD,THETANEW,PHINEW,
     +                    THETA,PHI,IRMD)
       ELSE
        RHO2NS(:,:,:)=DIMAG(RHO2NSNEW(:,:,:))       
        R2NEF(:,:,:)=DIMAG(R2NEFNEW(:,:,:))       
       ENDIF

       IDIM = IRMD*LMPOTD
       CALL DSCAL(IDIM,2.D0,RHO2NS(1,1,1),1)
       CALL DAXPY(IDIM,-0.5D0,RHO2NS(1,1,1),1,RHO2NS(1,1,2),1)
       CALL DAXPY(IDIM,1.0D0,RHO2NS(1,1,2),1,RHO2NS(1,1,1),1)
C
C --> do the same at the Fermi energy
C
       CALL DSCAL(IDIM,2.D0,R2NEF(1,1,1),1)
       CALL DAXPY(IDIM,-0.5D0,R2NEF(1,1,1),1,R2NEF(1,1,2),1)
       CALL DAXPY(IDIM,1.0D0,R2NEF(1,1,2),1,R2NEF(1,1,1),1)

       DO LM1=0,LMAXD1
       DO IE=1,IEMXD
       DO JSPIN=1,NSPIN
       DEN_out(LM1,IE,JSPIN) =  DEN(LM1,IE,JSPIN,1)
       ENDDO
       ENDDO
       ENDDO

       DEALLOCATE(VINS)
       DEALLOCATE(VNSPLL0)
       DEALLOCATE(VNSPLL1)
       DEALLOCATE(VNSPLL)
       DEALLOCATE(HLK)
       DEALLOCATE(JLK)
       DEALLOCATE(HLK2)
       DEALLOCATE(JLK2)
       DEALLOCATE(TMATSPH)
       DEALLOCATE(RLL)
       DEALLOCATE(SLL)
       DEALLOCATE(RLLLEFT)
       DEALLOCATE(SLLLEFT)
       DEALLOCATE(CDEN)
       DEALLOCATE(CDENLM)
       DEALLOCATE(CDENNS)
       DEALLOCATE(RHO2NSC,RHO2NSC_loop)
       DEALLOCATE(RHO2NSNEW)
       DEALLOCATE(R2NEFC,R2NEFC_loop)
       DEALLOCATE(R2NEFNEW)
       DEALLOCATE(R2ORBC)
       DEALLOCATE(CDENTEMP)
       DEALLOCATE(GFLLE_PART)
       DEALLOCATE(GFLLE)
       DEALLOCATE(DEN,DENLM)
       END 

