       SUBROUTINE TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT,SOCSCALE,
     +                        EZ,NSRA,CLEB,ICLEB,IEND,NCHEB,NPAN_TOT,
     +                        RPAN_INTERVALL,IPAN_INTERVALL,
     +                        RNEW,VINSNEW,THETA,PHI,I1,IPOT,
     &                        LLY,DELTAE,IDOLDAU,LOPT,WLDAU,
     &                        t_dtmatJij_at)
#ifdef CPP_OMP
       use omp_lib        ! necessary for omp functions
#endif
#ifdef CPP_MPI
       use mpi
#endif
       use mod_mympi, only: myrank, nranks, master
#ifdef CPP_MPI
     &                      ,distribute_linear_on_tasks, MPIadapt
       use mod_timing
#endif
       use mod_types, only: t_tgmat,t_inc,t_mpi_c_grid,init_tgmat,
     &                     t_lloyd,init_tlloyd, type_dtmatJijDij,
     &                     init_t_dtmatJij_at
       use mod_wunfiles, only: t_params
       use mod_jijhelp, only: calc_dtmatJij
       use mod_save_wavefun, only: t_wavefunctions, find_isave_wavefun,
     &                             save_wavefunc

       IMPLICIT NONE
       include 'inc.p'
       INTEGER NSPIN,NSRA,LMAX,IEND,IPOT,IELAST,NPAN_TOT,NCHEB
       INTEGER LMMAXD
       PARAMETER (LMMAXD= (LMAXD+1)**2)
       INTEGER LMMAXSO
       PARAMETER (LMMAXSO=2*LMMAXD)
       INTEGER LMPOTD
       PARAMETER (LMPOTD= (LPOTD+1)**2)
       INTEGER MMAXD
       PARAMETER (MMAXD=2*LMAXD+1)
       DOUBLE PRECISION CVLIGHT
       PARAMETER (CVLIGHT=274.0720442D0)
       DOUBLE COMPLEX CZERO,CONE,CI
       PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0),CI=(0d0,1d0))
       INTEGER NRMAXD
       PARAMETER (NRMAXD=NTOTD*(NCHEBD+1))
       DOUBLE PRECISION ZAT
       DOUBLE PRECISION SOCSCALE
       DOUBLE COMPLEX EZ(IEMXD),ERYD
       DOUBLE PRECISION CLEB(*)
       INTEGER ICLEB(NCLEB,4)
       DOUBLE PRECISION  RNEW(NRMAXD),
     +                   RPAN_INTERVALL(0:NTOTD),
     +                   VINSNEW(NRMAXD,LMPOTD,NSPOTD)
       INTEGER           IPAN_INTERVALL(0:NTOTD)
       DOUBLE COMPLEX  TMATLL(LMMAXSO,LMMAXSO) 
       type(type_dtmatJijDij) t_dtmatJij_at

   
       INTEGER I1,IR,IREC,USE_SRATRICK,NVEC,LM1,LM2,IE,IRMDNEW
       DOUBLE PRECISION THETA,PHI
       DOUBLE COMPLEX GMATPREFACTOR
       DOUBLE PRECISION, ALLOCATABLE :: VINS(:,:,:)
       DOUBLE COMPLEX,ALLOCATABLE :: VNSPLL0(:,:,:),VNSPLL1(:,:,:,:),
     +                               VNSPLL(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: HLK(:,:,:),JLK(:,:,:),
     +                                HLK2(:,:,:),JLK2(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: RLL(:,:,:,:),SLL(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: RLLLEFT(:,:,:,:),SLLLEFT(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: TMATSPH(:,:)
       DOUBLE COMPLEX DTMATLL(LMMAXSO,LMMAXSO),TMAT0(LMMAXSO,LMMAXSO) ! LLY
       DOUBLE COMPLEX ALPHALL(LMMAXSO,LMMAXSO), ! LLY
     +                DALPHALL(LMMAXSO,LMMAXSO), ! LLY
     +                ALPHA0(LMMAXSO,LMMAXSO), ! LLY
     +                AUX(LMMAXSO,LMMAXSO)         ! LLY
     +               ,ALPHASPH(2*(LMAX+1)) ! LLY
         
       INTEGER JLK_INDEX(2*LMMAXSO)
      ! LLoyd:
      INTEGER LLY,IDERIV,SIGNDE,IPIV(LMMAXSO)    ! LLY
      DOUBLE COMPLEX DELTAE,TRALPHA              ! LLY
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND)             ! LDAU
      INTEGER IDOLDAU,LOPT,LMLO,LMHI,MMAX                    ! LDAU
!     .. OMP ..
      integer nth,ith         ! total number of threads and thread id
#ifdef CPP_MPI
!       integer :: mytot
      integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
#endif
      integer :: ie_end, ie_num, ie_start, ierr
      
      !rhoqtest
      logical, external :: test, opt
      integer :: mu0, nscoef, ifinished!, lmgf0d
!       double complex, allocatable :: rllleft(:,:,:,:),sllleft(:,:,:,:)
      
      ! initialize inifished
      ifinished = 0

#ifdef CPP_OMP
! determine if omp parallelisation is used (compiled with -openmp flag and OMP_NUM_THREADS>1)
!$omp parallel shared(nth,ith) 
!$omp single 
      nth = omp_get_num_threads()
      if(t_inc%i_write>0) write(1337,*) 'nth =',nth
!$omp end single
!$omp end parallel
#else
      nth = 1
      ith = 0
#endif

       IRMDNEW= NPAN_TOT*(NCHEB+1)
       ALLOCATE(VINS(IRMDNEW,LMPOTD,NSPIN), stat=ierr)
       if(ierr/=0) stop '[tmat_newsolver] Error allocating VINS'
       VINS=0d0
       DO LM1=1,LMPOTD
        DO IR=1,IRMDNEW
         VINS(IR,LM1,1)=VINSNEW(IR,LM1,IPOT)
         VINS(IR,LM1,NSPIN)=VINSNEW(IR,LM1,IPOT+NSPIN-1)
        ENDDO
       ENDDO
cc set up the non-spherical ll' matrix for potential VLL'
       IF (NSRA.EQ.2) THEN
        USE_SRATRICK=1
       ELSEIF (NSRA.EQ.1) THEN
        USE_SRATRICK=0
       ENDIF
       ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW),
     &          VNSPLL1(LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1), stat=ierr)
       if(ierr/=0) stop '[tmat_newsolver] Error allocating vnspll0/1'
       VNSPLL0=CZERO
       CALL VLLMAT(1,NRMAXD,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINS,
     +                    CLEB,ICLEB,IEND,NSPIN,ZAT,RNEW,USE_SRATRICK)
! LDAU
       IF (IDOLDAU.EQ.1) THEN
        LMLO=LOPT**2+1
        LMHI=(LOPT+1)**2
         DO IR=1,IRMDNEW
          VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)=
     +                           VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)+
     +                           WLDAU(1:MMAX,1:MMAX,1)
         ENDDO
        LMLO=LMLO+LMMAXD
        LMHI=LMHI+LMMAXD
         DO IR=1,IRMDNEW
          VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)=
     +                           VNSPLL0(LMLO:LMHI,LMLO:LMHI,IR)+
     +                           WLDAU(1:MMAX,1:MMAX,2)
         ENDDO
       ENDIF
! LDAU

c initial allocate
       ! potential array
       IF (NSRA.EQ.2) THEN
        ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW,0:nth-1), stat=ierr)
       ELSE
        ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1), stat=ierr)
       ENDIF
       if(ierr/=0) stop '[tmat_newsolver] Error allocating vnspll'
       
       ! source terms (bessel and hankel functions)
       ALLOCATE(HLK(1:4*(LMAX+1),IRMDNEW,0:nth-1),
     &          JLK(1:4*(LMAX+1),IRMDNEW,0:nth-1),
     &          HLK2(1:4*(LMAX+1),IRMDNEW,0:nth-1),
     &          JLK2(1:4*(LMAX+1),IRMDNEW,0:nth-1), stat=ierr)
       if(ierr/=0) stop '[tmat_newsolver] Error allocating hlk etc.'
       
       ! spherical part of tmatrix (used with SRATRICK)
       ALLOCATE(TMATSPH(2*(LMAX+1),0:nth-1), stat=ierr)
       if(ierr/=0) stop '[tmat_newsolver] Error allocating tmatsph'
       
       !regular and irregular wavefunctions
       ALLOCATE(RLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),
     &          SLL(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1), stat=ierr)
       if(ierr/=0) stop '[tmat_newsolver] Error allocating rll/sll'
       
       !left regular and irregular wavefunctions (used here only in case of XCPL or saving of left wavefunctions)
       if( opt('XCPL    ') .or. (t_wavefunctions%save_rllleft
     &                         .or.t_wavefunctions%save_sllleft)
     &                      .or.test('rhoqtest')) then
          ALLOCATE(RLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),
     &             SLLLEFT(NSRA*LMMAXSO,LMMAXSO,IRMDNEW,0:nth-1),
     &             stat=ierr)
          if(ierr/=0) stop '[tmat_newsolver] Error allocating r/sllleft'
       else
          ALLOCATE(RLLLEFT(1,1,1,0:nth-1), SLLLEFT(1,1,1,0:nth-1),
     &             stat=ierr)
          if(ierr/=0) stop 
     &            '[tmat_newsolver] Error allocating r/sllleft dummies'
       end if ! ( opt('XCPL    ') .or. ... )
          
c energy loop
       if(myrank==master.and.(t_inc%i_write>0)) 
     &        WRITE(1337,*) 'atom: ',I1,' NSRA:',NSRA

#ifdef CPP_MPI
       call distribute_linear_on_tasks(t_mpi_c_grid%nranks_at, 
     &           t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at+(i1-1), ! print this info only for first atom at master
     &                  master,IELAST,ntot_pT,ioff_pT,.true.,.true.)
     
       ie_start = ioff_pT(t_mpi_c_grid%myrank_at)
       ie_end   = ntot_pT(t_mpi_c_grid%myrank_at)

       t_mpi_c_grid%ntot2=ie_end   !t_mpi_c_grid%dims(1)
       if(.not. (allocated(t_mpi_c_grid%ntot_pT2) .or. 
     &           allocated(t_mpi_c_grid%ioff_pT2))) 
     &    allocate(t_mpi_c_grid%ntot_pT2(0:t_mpi_c_grid%nranks_at-1),
     &             t_mpi_c_grid%ioff_pT2(0:t_mpi_c_grid%nranks_at-1))
       t_mpi_c_grid%ntot_pT2 = ntot_pT
       t_mpi_c_grid%ioff_pT2 = ioff_pT
       ! now initialize arrays for tmat, gmat, and gref
       call init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
       if(lly.ne.0) call init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)
       
       
       if(test('rhoqtest')) then
          if(ielast/=3) stop 'Error: wrong energy contour for rhoqtest'
          ie_start=1
          ie_end=1
       end if
       
#else
       if(.not. (allocated(t_mpi_c_grid%ntot_pT2) .or. 
     &           allocated(t_mpi_c_grid%ioff_pT2))) 
     & allocate(t_mpi_c_grid%ntot_pT2(1),t_mpi_c_grid%ioff_pT2(1))
       t_mpi_c_grid%ntot2=IELAST
       t_mpi_c_grid%ntot_pT2 = IELAST 
       t_mpi_c_grid%ioff_pT2 = 0      
       ! now initialize arrays for tmat, gmat, and gref
       call init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
       if(lly.ne.0) call init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)

       ie_start = 0
       ie_end = IELAST
#endif

       !for Jij-tensor calculation: allocate array to hold additional t-matrices
       call init_t_dtmatJij_at(t_inc, t_mpi_c_grid, t_dtmatJij_at)
       
                
       !initialize wfsave
       if(t_inc%i_iteration==0) then
            call find_isave_wavefun(t_wavefunctions)
            ! reset Nwfsavemax to 0 if test option 'STOP1B  ' is found
            ! to prevent unnessesary storing of wavefunctions
            if(test('STOP1B  ')) t_wavefunctions%Nwfsavemax = 0
       endif
         

#ifdef CPP_OMP
!$omp parallel do default(none) 
!$omp& private(eryd,ie,ir,nvec,lm1,lm2,gmatprefactor) 
!$omp& private(jlk_index,tmatll,ith,irec, ie_num)
!$omp& private(tralpha, aux, ideriv, ipiv) 
!$omp& private(alpha0)
!$omp& private(alphall)
!$omp& private(tmat0) 
!$omp& private(alphasph) 
!$omp& private(dtmatll) 
!$omp& private(dalphall) 
!$omp& firstprivate(t_inc)
!$omp& shared(nspin,nsra,lmax,iend,ipot,ielast,npan_tot,ncheb) 
!$omp& shared(zat,socscale,ez,cleb,rnew,nth) 
!$omp& shared(rpan_intervall,vinsnew,ipan_intervall) 
!$omp& shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0) 
!$omp& shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,rllleft,sllleft)
!$omp& shared(tmatsph, ie_end,t_tgmat,t_lloyd, ie_start, t_dtmatjij_at)
!$omp& shared(lly,deltae,i1,t_mpi_c_grid, t_wavefunctions, icleb)
#endif
       DO ie_num=1,ie_end
       
         IE = ie_start+ie_num
         
#ifdef CPP_MPI
         !start timing measurement for this pair of ie and i1, needed for MPIadapt
         call timing_start('time_1a_ieiatom')
#endif

         ! get current thread
         if (nth>=1) then
#ifdef CPP_OMP
             ith = omp_get_thread_num()
#endif
         else
            ith = 0
         endif

! In case of Lloyds formula the derivative of t is needed. 
! Then calculate t at E+dE, E-dE and average for t, subtract for dt/dE
        TMATLL=CZERO
        ALPHALL=CZERO ! LLY
        DTMATLL=CZERO ! LLY
        DALPHALL=CZERO ! LLY
        IDERIV=0
        IF (LLY.NE.0) IDERIV=1
        DO SIGNDE=-IDERIV,IDERIV,2
         ERYD = EZ(IE)+DFLOAT(SIGNDE)*DELTAE/2.D0 ! LLY
#ifdef CPP_OMP
!$omp critical
#endif
        if(t_inc%i_write>0) WRITE(1337,*) 'energy:',IE,'',ERYD
#ifdef CPP_OMP        
        if(ie==1.and.(t_inc%i_write>0)) 
     &      write(1337,*) 'nested omp?',omp_get_nested()
!$omp end critical
#endif

c contruct the spin-orbit coupling hamiltonian and add to potential
        CALL SPINORBIT_HAM(LMAX,LMMAXD,VINS,RNEW,
     +                     ERYD,ZAT,CVLIGHT,SOCSCALE,NSPIN,LMPOTD,
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
        TMATSPH(:,ith)=CZERO
         CALL CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,CVLIGHT,ERYD,
     +           LMPOTD,LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +                JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),HLK2(:,:,ith),
     +                JLK2(:,:,ith),GMATPREFACTOR,TMATSPH(:,ith),
     +                ALPHASPH,USE_SRATRICK)
        ENDIF
cc calculate the tmat and wavefunctions
        RLL(:,:,:,ith)=CZERO
        SLL(:,:,:,ith)=CZERO

cc right solutions
        TMAT0=CZERO
        ALPHA0=CZERO ! LLY
C faster calculation of RLL.
C no irregular solutions are needed in self-consistent iterations
C because the t-matrix depends only on RLL
        IF(OPT('RLL-SLL ') .and. .not.OPT('XCPL    ')) THEN
        CALL rll_global_solutions(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),
     +              RLL(:,:,:,ith),TMAT0(:,:),NCHEB,
     +              NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +              NRMAXD,NSRA,JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),
     +              HLK2(:,:,ith),JLK2(:,:,ith),GMATPREFACTOR,
     +              '1',USE_SRATRICK,ALPHA0(:,:))
        ELSE
        CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),
     +              RLL(:,:,:,ith),SLL(:,:,:,ith),TMAT0(:,:),NCHEB,
     +              NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +              NRMAXD,NSRA,JLK_INDEX,HLK(:,:,ith),JLK(:,:,ith),
     +              HLK2(:,:,ith),JLK2(:,:,ith),GMATPREFACTOR,
     +              '1','1','0',USE_SRATRICK,ALPHA0(:,:))
        END IF

        ! devide small component by cvlight so that it actually is small
        IF (NSRA.EQ.2) THEN
          RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
     +            RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
          SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
     +            SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
        ENDIF

c add spherical contribution of tmatrix
         IF (USE_SRATRICK.EQ.1) THEN
          DO LM1=1,LMMAXSO
           TMAT0(LM1,LM1)=TMAT0(LM1,LM1)+
     +                               TMATSPH(JLK_INDEX(LM1),ith)
          ENDDO
          DO LM2=1,LMMAXSO
           DO LM1=1,LMMAXSO
            ALPHA0(LM1,LM2)=ALPHASPH(JLK_INDEX(LM1))*
     +                                ALPHA0(LM1,LM2)        ! LLY
           ENDDO
          ENDDO
         ENDIF
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            TMATLL(LM1,LM2)=TMATLL(LM1,LM2)+TMAT0(LM1,LM2)
           ENDDO
          ENDDO          
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            ALPHALL(LM1,LM2)=ALPHALL(LM1,LM2)+ALPHA0(LM1,LM2)
           ENDDO
          ENDDO          
         IF (LLY.NE.0) THEN
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            DTMATLL(LM1,LM2)=DTMATLL(LM1,LM2)+
     +                  DFLOAT(SIGNDE)*TMAT0(LM1,LM2) ! LLY
            DALPHALL(LM1,LM2)=DALPHALL(LM1,LM2)+
     +                  DFLOAT(SIGNDE)*ALPHA0(LM1,LM2) ! LLY
           ENDDO
          ENDDO
         ENDIF
        ENDDO ! SIGNDE=-IDERIV,IDERIV,2 ! LLY

! average values of t-matrix and alpha at e+de and e-de
        DO LM1=1,LMMAXSO        
         DO LM2=1,LMMAXSO        
          TMATLL(LM1,LM2)=TMATLL(LM1,LM2)/DFLOAT(1+IDERIV) ! LLY
          ALPHALL(LM1,LM2)=ALPHALL(LM1,LM2)/DFLOAT(1+IDERIV) ! LLY
         IF (LLY.NE.0) THEN
! Contruct derivative of t-matrix and alpha
          DTMATLL(LM1,LM2)=DTMATLL(LM1,LM2)/DELTAE ! LLY
          DALPHALL(LM1,LM2)=DALPHALL(LM1,LM2)/DELTAE ! LLY
         ENDIF
         ENDDO
        ENDDO
        IF (LLY.NE.0) THEN
! calculate Tr[alpha^-1*dalpha/de] for LLoyd's formula
        ALPHA0=CZERO ! LLY
        AUX=CZERO ! LLY
        CALL ZGEINV1(ALPHALL,ALPHA0,AUX,IPIV,LMMAXSO)
        CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,ALPHA0,LMMAXSO,
     +             DALPHALL,LMMAXSO,CZERO,AUX,LMMAXSO) ! LLY
! Trace of AUX 
        TRALPHA=CZERO ! LLY
        DO LM1=1,LMMAXSO 
         TRALPHA=TRALPHA+AUX(LM1,LM1) ! LLY
        ENDDO
        ENDIF ! LLY


        if(test('rhoqtest') .and. ie==2) then
           ! read in mu0 atom index
           open(9999,file='mu0')
           read(9999,*) mu0, nscoef
           close(9999)
        end if

! Calculate additional t-matrices for Jij-tensor calculation
        if(t_dtmatJij_at%calculate 
     &    .or. ( t_wavefunctions%isave_wavefun(i1, ie)>0 .and. 
     &          (t_wavefunctions%save_rllleft .or.
     &           t_wavefunctions%save_sllleft) )
     &    .or. ((test('rhoqtest') .and. ie==2).and.(i1==mu0)) !rhoqtest
     &     ) then
c        LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
c        LLLLLLLLLLL Calculate the left-hand side solution LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
c        this needs to be done for the calculation of t-matrices for Jij tensor or if wavefunctions should be saved
c
c        contruct the spin-orbit coupling hamiltonian and add to potential
         CALL SPINORBIT_HAM(LMAX,LMMAXD,VINS,RNEW,
     +                     ERYD,ZAT,CVLIGHT,SOCSCALE,NSPIN,LMPOTD,
     +                     THETA,PHI,IPAN_INTERVALL,RPAN_INTERVALL,
     +                     NPAN_TOT,NCHEB,IRMDNEW,NRMAXD,
     +                     VNSPLL0(:,:,:),VNSPLL1(:,:,:,ith),
     +                     'transpose')

cc       extend matrix for the SRA treatment
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

cc       calculate the source terms in the Lippmann-Schwinger equation
cc       these are spherical hankel and bessel functions
         HLK(:,:,ith)=CZERO
         JLK(:,:,ith)=CZERO
         HLK2(:,:,ith)=CZERO
         JLK2(:,:,ith)=CZERO
         GMATPREFACTOR=CZERO
         jlk_index = 0
         CALL RLLSLLSOURCETERMS(NSRA,NVEC,ERYD,RNEW,IRMDNEW,NRMAXD,LMAX,
     +                         LMMAXSO,1,JLK_INDEX,HLK(:,:,ith),
     +                         JLK(:,:,ith),HLK2(:,:,ith),JLK2(:,:,ith),
     +                         GMATPREFACTOR)

cc       using spherical potential as reference
c        notice that exchange the order of left and right hankel/bessel functions
         IF (USE_SRATRICK.EQ.1) THEN
          TMATSPH(:,ith)=CZERO
          CALL CALCSPH(NSRA,IRMDNEW,NRMAXD,LMAX,NSPIN,ZAT,CVLIGHT,ERYD,
     +           LMPOTD,LMMAXSO,RNEW,VINS,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +               JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),HLK(:,:,ith),
     +               JLK(:,:,ith),GMATPREFACTOR,ALPHASPH,TMATSPH(:,ith),
     +               USE_SRATRICK)
         ENDIF
         
cc       calculate the tmat and wavefunctions
         RLLLEFT(:,:,:,ith)=CZERO
         SLLLEFT(:,:,:,ith)=CZERO

cc       left solutions
c        notice that exchange the order of left and right hankel/bessel functions
         TMAT0=CZERO
         ALPHA0=CZERO ! LLY
C faster calculation of RLL.
C no left solutions are needed in self-consistent iterations
C because the t-matrix depends only on RLL
         IF(OPT('RLL-SLL ') .and. .not.OPT('XCPL    ')) THEN
         ! do nothing
         ELSE
         CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL(:,:,:,ith),
     +             RLLLEFT(:,:,:,ith),SLLLEFT(:,:,:,ith),TMAT0,NCHEB,
     +             NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +             NRMAXD,NSRA,JLK_INDEX,HLK2(:,:,ith),JLK2(:,:,ith),
     +             HLK(:,:,ith),JLK(:,:,ith),GMATPREFACTOR,
     +             '1','1','0',USE_SRATRICK,ALPHA0)
         END IF
         IF (NSRA.EQ.2) THEN
          RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
     +            RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
          SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)=
     +            SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:,ith)/CVLIGHT
         ENDIF
         
       
         if(test('rhoqtest')) then
           !
           open(9999, file='params.txt')
           write(9999,*) lmmaxso, t_params%natyp
           write(9999,*) t_params%naez, t_params%nclsd, t_params%nrd, 
     &                   t_params%nembd1-1, t_params%lmax
           write(9999,*) t_params%alat
           close(9999)
           !
           open(9999, file='host.txt')
           write(9999,*) t_params%rbasis(1:3,1:t_params%natyp)
           write(9999,*) t_params%rcls(1:3,1:t_params%nclsd,
     &                                 1:t_params%nclsd), 
     &                   t_params%rr(1:3,0:t_params%nrd), 
     &                   t_params%atom(1:t_params%nclsd,
     &                              1:t_params%naez+t_params%nembd1-1)
!      
!      &t_params%ATOM(t_params%NACLSD,t_params%NEMBD2),
!      
         write(9999,*) t_params%cls(1:t_params%naez+t_params%nembd1-1), 
     &                   t_params%ezoa(1:t_params%nclsd,
     &                              1:t_params%naez+t_params%nembd1-1), 
     &                   t_params%nacls(1:t_params%nclsd)
           close(9999)
           !
           open(9999, file='wavefunctions.txt')
           write(9999,'(100I9)') ntotd, npan_tot, ncheb, nsra, irmdnew
           write(9999,'(1000E26.17)') rnew(1:irmdnew)
!            write(*,*) ntotd, npan_tot, ncheb, nsra, irmdnew, lmmaxso, ie
           do ir=1,irmdnew
             do lm1=1,nsra*lmmaxso
               do lm2=1,lmmaxso
                  write(9999,'(20000E16.7)') Rll(lm1, lm2, ir, ith), 
     &                                       Rllleft(lm1, lm2, ir, ith)
               end do
             end do
           enddo
           do lm1=0,npan_tot
             write(9999,'(E16.7,I9)') rpan_intervall(lm1),
     &                                ipan_intervall(lm1)
           enddo
           close(9999)
           !
           ifinished = 1
         end if ! test('rhoqtest')
        
c        LLLLLLLLLLL Calculate the left-hand side solution LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
c        LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
        end if!t_dtmatJij_at%calculate .or. t_wavefunctions%Nwfsavemax>0
        
        
        ! save_wavefuncions         
        if(t_wavefunctions%Nwfsavemax>0) then
          ! here all four (left, right, regular and irregular) are stored, the memory demand cound be reduced by a factor 2 if only the right solution would be computed here and saved and the left solution would be calculated later in main1c
          call save_wavefunc(t_wavefunctions, rll, rllleft, 
     &          sll, sllleft, i1, ie, NSRA, LMMAXSO, IRMDNEW, ith)
        end if
        
        if(t_dtmatJij_at%calculate) then
         call calc_dtmatJij(LMMAXD,LMMAXSO,LMPOTD,NTOTD,NRMAXD,
     +        NSRA,IRMDNEW,NSPIN,VINS,RLLLEFT(:,:,:,ith),RLL(:,:,:,ith),
     +        RPAN_INTERVALL,
     +        IPAN_INTERVALL,NPAN_TOT,NCHEB,CLEB,ICLEB,IEND,NCLEB,RNEW,
     +        t_dtmatJij_at%dtmat_xyz(:,:,:,ie_num))

        end if!t_dtmatJij_at%calculate
         
         ! writeout
#ifdef CPP_OMP        
!$omp critical
#endif
         if (t_tgmat%tmat_to_file) then
            IREC = IE + IELAST*(I1-1)
                if(test('rhoqtest')) then
!                 write(*,*) 'writing tmat', irec, lmmaxso

                  if(ie_num==1.and.i1==1) then
                    write(*,*)                      ! status bar
                    write(*,*) 'rhoq: write-out t-mat', 
     &                               ie_end, t_params%natyp
                    write(*,
     &   '("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")')
     &                          0, 20, 40, 60, 80, 100
                    write(*,FMT=190) !beginning of statusbar
                  endif
                  
                  if(t_params%NATYP*ie_end>=50) then
                    if(mod( I1+t_params%natyp*(ie_num-1),
     &                 (t_params%NATYP*ie_end/50))==0 ) write(6,FMT=200)
                  else
                    write(6,FMT=200)
                  end if
                           
!                 write(*,*) 'rotating with', theta, phi
!             LMGF0D= (LMAXD+1)**2
!             CALL ROTATEMATRIX(TMATLL,THETA,PHI,LMGF0D,0)
                end if     
                           
            WRITE(69,REC=IREC) TMATLL(:,:)
            ! human readable writeout if test option is hit
            if(test('fileverb')) then
               WRITE(696969,'(i9,20000F15.7)') irec, TMATLL(:,:)
            end if
         else
#ifdef CPP_MPI
            irec = ie_num +
     +       ie_end*(I1-t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie)-1)
#else
            irec = ie_num + ie_end * (i1-1)
#endif
            t_tgmat%tmat(:,:,irec) = TMATLL(:,:)
         end if
         IF (LLY.NE.0) THEN
           if(t_lloyd%dtmat_to_file) then
             IREC = IE + IELAST*(I1-1)
             WRITE(691,REC=IREC) DTMATLL(:,:)    ! LLY      
           else
             irec = ie_num + ie_end*(i1-1)
             t_lloyd%dtmat(:,:,irec) = DTMATLL(:,:)
           end if
           if(t_lloyd%tralpha_to_file) then
             IREC = IE + IELAST*(I1-1)
             WRITE(692,REC=IREC) TRALPHA                              ! LLY
           else
             irec = ie_num + ie_end*(i1-1)
             t_lloyd%tralpha(irec) = TRALPHA
           end if
         ENDIF
#ifdef CPP_OMP
!$omp end critical         
#endif

#ifdef CPP_MPI
         !stop timing measurement for this pair of ie and i1, needed for MPIadapt
         if(MPIadapt)
     & call timing_stop('time_1a_ieiatom',save_out=timings_1a(ie, i1))
#endif

       ENDDO ! IE loop 
#ifdef CPP_OMP
!$omp end parallel do
#endif


190     FORMAT('                 |'$)   ! status bar
200     FORMAT('|'$)                    ! status bar
        if(test('rhoqtest').and.i1==t_params%natyp) write(6,*) ! status bar
        ! finished kpts status bar

       DEALLOCATE(VINS, VNSPLL0, VNSPLL1, VNSPLL, stat=ierr)
       if(ierr/=0) stop '[tmat_newsolver] Error dellocating vins etc.'
       DEALLOCATE(HLK, JLK, HLK2, JLK2, stat=ierr)
       if(ierr/=0) stop '[tmat_newsolver] Error dellocating hlk etc.'
       DEALLOCATE(TMATSPH, RLL, SLL, RLLLEFT, SLLLEFT, stat=ierr)
       if(ierr/=0) stop '[tmat_newsolver] Error dellocating rll etc.'

       END 

