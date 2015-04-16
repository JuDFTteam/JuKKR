       SUBROUTINE TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,RMESH,ZAT,SOCSCALE,
     +                        EZ,NSRA,CLEB,ICLEB,IEND,NCHEB,NPAN_TOT,
     +                        RPAN_INTERVALL,IPAN_INTERVALL,
     +                        RNEW,VINSNEW,THETA,PHI,I1,IPOT,
     &                        LLY,DELTAE)
       use omp_lib        ! necessary for omp functions
       IMPLICIT NONE
       include 'inc.p'
       INTEGER NSPIN,NSRA,LMAX,IEND,IPOT,IELAST,NPAN_TOT,NCHEB
       INTEGER LMMAXD
       PARAMETER (LMMAXD= (LMAXD+1)**2)
       INTEGER LMMAXSO
       PARAMETER (LMMAXSO=2*LMMAXD)
       INTEGER LMPOTD
       PARAMETER (LMPOTD= (LPOTD+1)**2)
       DOUBLE PRECISION CVLIGHT
       PARAMETER (CVLIGHT=274.0720442D0)
       DOUBLE COMPLEX CZERO
       PARAMETER (CZERO=(0d0,0d0))
       INTEGER NRMAXD
       PARAMETER (NRMAXD=NTOTD*(NCHEBD+1))
       DOUBLE PRECISION ZAT
       DOUBLE PRECISION SOCSCALE
       DOUBLE COMPLEX EZ(IEMXD),ERYD
       DOUBLE PRECISION RMESH(IRMD),CLEB(*)
       INTEGER ICLEB(NCLEB,4)
       DOUBLE PRECISION  RNEW(NRMAXD),
     +                   RPAN_INTERVALL(0:NTOTD),
     +                   VINSNEW(NRMAXD,LMPOTD,NSPOTD)
       INTEGER           IPAN_INTERVALL(0:NTOTD)
       DOUBLE COMPLEX  
     +   TMATLL(LMMAXSO,LMMAXSO,ielast) 
       DOUBLE COMPLEX DTMATLL(LMMAXSO,LMMAXSO),TMATLLAV(LMMAXSO,LMMAXSO) ! LLY
   
       INTEGER I1,IR,IREC,USE_SRATRICK,NVEC,LM1,LM2,IE,IRMDNEW
       DOUBLE PRECISION THETA,PHI
       DOUBLE COMPLEX GMATPREFACTOR
       DOUBLE PRECISION, ALLOCATABLE :: VINS(:,:,:)
       DOUBLE COMPLEX,ALLOCATABLE :: VNSPLL0(:,:,:),VNSPLL1(:,:,:,:),
     +                               VNSPLL(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: HLK(:,:,:),JLK(:,:,:),
     +                                HLK2(:,:,:),JLK2(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: RLL(:,:,:,:),SLL(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: TMATSPH(:,:), tmat_out(:,:,:)   !tmat_out necessary for parallel ie loop
       INTEGER JLK_INDEX(2*LMMAXSO)
      ! LLoyd:
      INTEGER LLY,IDERIV,SIGNDE        ! LLY
      DOUBLE COMPLEX DELTAE            ! LLY
!     .. OMP ..
      integer nth,ith         ! total number of threads and thread id
! determine if omp parallelisation is used (compiled with -openmp flag and OMP_NUM_THREADS>1)
!$omp parallel shared(nth,ith) 
!$omp single 
      nth = 1
      ith = 0
      nth = omp_get_num_threads()
!$omp end single
!$omp end parallel
      ! write(*,*) 'nth =',nth

       IRMDNEW= NPAN_TOT*(NCHEB+1)
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

       ALLOCATE(HLK(1:4*(LMAX+1),IRMDNEW,0:nth-1))
       ALLOCATE(JLK(1:4*(LMAX+1),IRMDNEW,0:nth-1))
       ALLOCATE(HLK2(1:4*(LMAX+1),IRMDNEW,0:nth-1))
       ALLOCATE(JLK2(1:4*(LMAX+1),IRMDNEW,0:nth-1))
       ALLOCATE(TMATSPH(2*(LMAX+1),0:nth-1))
       ALLOCATE(RLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1))
       ALLOCATE(SLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1))
       ALLOCATE(tmat_out(lmmaxso,lmmaxso,ielast))
          
c energy loop
       WRITE(6,*) 'atom: ',I1,' NSRA:',NSRA

!$omp parallel do default(none) 
!$omp& private(eryd,ie,i1,ir,irec,nvec,lm1,lm2,gmatprefactor) 
!$omp& private(jlk_index,tmatll,ith) 
!$omp& shared(nspin,nsra,lmax,iend,ipot,ielast,npan_tot,ncheb) 
!$omp& shared(zat,socscale,ez,rmesh,cleb,rnew,nth) 
!$omp& shared(rpan_intervall,vinsnew,ipan_intervall) 
!$omp& shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0) 
!$omp& shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,tmat_out)
!$omp& shared(tmatsph)
       DO IE=1,IELAST
         ! get current thread
         if (nth>=1) then
            ith = omp_get_thread_num()
         else
            ith = 0
         endif
        ERYD = EZ(IE)
!$omp critical
        WRITE(6,*) 'energy:',IE,'',ERYD
        !write(*,*) 'nested omp?',omp_get_nested()
!$omp end critical

c contruct the spin-orbit coupling hamiltonian and add to potential
        CALL SPINORBIT_HAM(LMAX,LMMAXD,VINS,RNEW,
     +                     ERYD,ZAT,CVLIGHT,SOCSCALE,NSRA,NSPIN,LMPOTD,
     +                     THETA,PHI,IPAN_INTERVALL,RPAN_INTERVALL,
     +                     NPAN_TOT,NCHEB,IRMDNEW,NRMAXD,
     +                     VNSPLL0(:,:,:),VNSPLL1(:,:,:,ith),'1')
cc extend matrix for the SRA treatment
        VNSPLL(:,:,:,ith)=CZERO
        IF (NSRA.EQ.2) THEN
         IF (USE_SRATRICK.EQ.0) THEN
          CALL VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,
     +        LMMAXSO,IRMDNEW,NRMAXD,ERYD,CVLIGHT,LMAX,0,'Ref=0')
         ELSEIF (USE_SRATRICK.EQ.1) THEN
          CALL VLLMATSRA(VNSPLL1(:,:,:,ith),VNSPLL(:,:,:,ith),RNEW,
     +        LMMAXSO,IRMDNEW,NRMAXD,ERYD,CVLIGHT,LMAX,0,'Ref=Vsph')
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
        CALL RLLSLLSOURCETERMS(NSRA,NVEC,ERYD,RNEW,IRMDNEW,NRMAXD,LMAX,
     +                         LMMAXSO,1,JLK_INDEX,HLK(:,:,ith),
     +                         JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),
     +                         GMATPREFACTOR)
cc using spherical potential as reference
        IF (USE_SRATRICK.EQ.1) THEN
         CALL CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,CVLIGHT,ERYD,
     +           LMPOTD,LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +                JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),HLK2(:,:,ith),
     +                JLK2(:,:,ith),GMATPREFACTOR,TMATSPH(:,ith),
     +                USE_SRATRICK)
        ENDIF

cc calculate the tmat and wavefunctions
        RLL(:,:,:,ith)=CZERO
        SLL(:,:,:,ith)=CZERO

cc right solutions
        TMATLL=CZERO
        CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),
     +              RLL(:,:,:,ith),SLL(:,:,:,ith),TMATLL(:,:,ie),NCHEB,
     +              NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +              NRMAXD,NSRA,JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),
     +              HLK2(:,:,ith),JLK2(:,:,ith),GMATPREFACTOR,
     +              '1','1','0',USE_SRATRICK)
!     &              ,ith) ! test fivos
        IF (NSRA.EQ.2) THEN
c         RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
c     +            RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/C
c         SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
c     +            SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/C
        ENDIF


c add spherical contribution of tmatrix
         IF (USE_SRATRICK.EQ.1) THEN
          DO LM1=1,LMMAXSO
       TMATLL(LM1,LM1,ie)=TMATLL(LM1,LM1,ie)+TMATSPH(JLK_INDEX(LM1),ith)
          ENDDO
         ENDIF
         tmat_out(:,:,ie) = tmatll(:,:,ie)
       ENDDO ! IE loop 
!$omp end parallel do

! serial write out after parallel calculation of tmat
       DO IE=1,IELAST 
         IREC = IE + IELAST*(I1-1)
         WRITE(69,REC=IREC) TMAT_out(:,:,ie)
c         write(696969,*) TMAT_out(:,:,ie)
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
       END 
