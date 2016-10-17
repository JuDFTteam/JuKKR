      SUBROUTINE KKRMAT01(BZKP,NOFKS,GS,VOLCUB,TINVLL,RROT,
     &                    NSHELL,NSDIA,ALAT,NSYMAT,
     &                    NAEZ,CLS,NACLS,NACLSMAX,RR,EZOA,ATOM,
     &                    NSH1,NSH2,GINP,RBASIS,RCLS,
     &                    TINVBUP,TINVBDOWN,VACFLAG,NLBASIS,NRBASIS,
     &                    FACTL,ICHECK,INVMOD,IDECI,SRREL,IRREL,NRREL,
     &     DTREFLL,DTMATLL,DGINP,REFPOT,LLY_GRTR,TRACET,CFCTOR,LLY, ie)   ! LLY
! **********************************************************************
! * Performs k-space integration, determines scattering path operator  *
! *                tau = (g(k,e)-t**-1)**-1                            *
! * and Greens function of the real system -> GS(*,*,*,*)              *
! *                                                                    *
! * new version 10.99: up -> left , down -> right, for decimation      *
! *                                                                    *
! * For KREL = 1 (relativistic mode)                                   *
! *  NPOTD = 2 * NATYPD                                                *
! *  LMMAXD = 2 * (LMAXD+1)^2                                          *
! *  NSPIND = 1                                                        *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green      *
! *          function, set up in the spin-independent non-relativstic  *
! *          (l,m_l)-representation                                    *
! *                                                                    *
! *  Modifications according to H. Hoehler ( July 2002)                *
! *     define Fourier transformation as                               *
! *                                                                    *
! *             /         n   0          n                             *
! *   G mu mu'= | sum_n G mu mu' exp(-iKR ) +                          *
! *     L  L'   \         L   L'                                       *
! *                                                                    *
! *                                 n   0            n   \    1        *
! *                         sum_n G mu mu' exp(-iK(-R )) | * ---       *
! *                                 L   L'               /    2        *
! *                                                                    *
! *     this operation has to be done to satisfy the point symmetry;   *
! *     the application of the fourier transformation is just an       *
! *     approximation for the tb system, since the transl. invariance  *
! *     is not satisfied --> force it by R, -R                         *
! *                                                                    *
! **********************************************************************
#ifdef CPP_MPI
      use mpi
      
      use mod_types, only: t_mpi_c_grid
      use mod_mympi, only: myrank, nranks, master,
     &                     distribute_linear_on_tasks
#else
      use mod_mympi, only: myrank, nranks, master 
#endif
      use mod_types, only: t_inc
#ifdef CPP_TIMING
      use mod_timing
#endif
      IMPLICIT NONE
!     ..
!     .. Parameters ..
      include 'inc.p'
      INTEGER LMAX,NSYMAXD
      PARAMETER (LMAX=LMAXD,NSYMAXD=48)
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAXD+1)**2)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+KORBIT+1) * (LMAX+1)**2)
      INTEGER ALM
      PARAMETER (ALM = NAEZD*LMMAXD)
      INTEGER ALMGF0
      PARAMETER (ALMGF0 = NAEZD*LMGF0D)
      DOUBLE COMPLEX CZERO,CMI,CONE
      PARAMETER ( CZERO=(0D0,0D0), CMI=(0D0,-1D0), CONE=(1D0,0D0) )
!     ..
!     .. Scalar arguments ..
      DOUBLE PRECISION ALAT
      INTEGER NAEZ,NOFKS,NSHELL,NSYMAT,NSDIA,NACLSMAX
      INTEGER IDECI,INVMOD,NLBASIS,NRBASIS
      INTEGER LLY ! LLY <> 0 --> use Lloyds formula
!     ..
!     .. Array arguments ..
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD),REFPOT(*) ! REFPOT(NAEZD+NEMBD) 
      DOUBLE COMPLEX GINP(LMGF0D*NACLSMAX,LMGF0D,*),  ! Gref
     &               DGINP(LMGF0D*NACLSMAX,LMGF0D,*), ! LLY dGref/dE
     +               GS(LMMAXD,LMMAXD,NSYMAXD,*),
     +               TINVLL(LMMAXD,LMMAXD,NAEZ),FACTL(LMMAXD,LMMAXD)
      DOUBLE COMPLEX TINVBUP(LMMAXD,LMMAXD,*),TINVBDOWN(LMMAXD,LMMAXD,*)
     &              ,DTREFLL(LMMAXD,LMMAXD,NREFD), ! LLY dtref/dE
     &               DTMATLL(LMMAXD,LMMAXD,NAEZD), ! LLY  dt/dE (should be av.-tmatrix in CPA)
     &               T_AUX(LMMAXD,LMMAXD,NAEZD),   ! LLY auxiliary array for t-matrix manipulation
     &               GAUX1(LMMAXD,LMMAXD),GAUX2(LMMAXD,LMMAXD),! LLY
     &               GAUX3(LMMAXD,LMMAXD) ! LLY
      DOUBLE PRECISION BZKP(3,*),RROT(48,3,*),VOLCUB(*),RBASIS(3,*), 
     +                 RR(3,0:NRD),RRM(3,0:NRD),RCLS(3,NACLSD,*)
      DOUBLE COMPLEX SRREL(2,2,LMMAXD)
      INTEGER IRREL(2,2,LMMAXD),NRREL(2,LMMAXD)
      INTEGER IQ1,IQ2,IOFF1,IOFF2,JOFF1,JOFF2
      INTEGER IKM1,IKM2,IS,N1,N2,J1,J2,I2
      DOUBLE COMPLEX CSUM1,CSUM2,TRACE,TRACET  ! LLY Lloyd
      DOUBLE COMPLEX LLY_GRTR_K,LLY_GRTR ! Trace Eq.5.38 PhD Thiess  (k-dependent and integrated) ! LLY Lloyd
      INTEGER ATOM(NACLSD,*),CLS(*),EZOA(NACLSD,*),
     &        NACLS(*),NSH1(*),NSH2(*) 
      LOGICAL VACFLAG(2)
!     ..
!     .. Local scalars ..
      DOUBLE COMPLEX CARG,CITPI,CFCTOR
      DOUBLE PRECISION ZKTR
      INTEGER I,I1,ILM,ISYM,IU,J,JLM,IL1,KPT,LM,LM1,LM2,
     +        NS,IL2,
     &        JL1,JL2
!     ..
!     .. Local arrays ..
!----------------------------------------------------------------
      DOUBLE COMPLEX GLLKE(:,:),GLLKEM(:,:),GLLKEN(:,:)
      DOUBLE COMPLEX DGLLKE(:,:),DGLLKEM(:,:),DGLLKEN(:,:),GREFLLKE(:,:) ! LLY
      ALLOCATABLE DGLLKE,DGLLKEM,DGLLKEN,GREFLLKE ! LLY

      DOUBLE COMPLEX GLLKE0V(:,:),GLLKE0V2(:,:),GLLKETV(:,:) ! for VIRTUAL ATOMS
      ALLOCATABLE GLLKE0V,GLLKE0V2,GLLKETV
      DOUBLE COMPLEX GLLKETV_new(:,:) ! for VIRTUAL ATOMS
      ALLOCATABLE GLLKETV_new

      DOUBLE COMPLEX GLLKE0(:,:),GLLKE0M(:,:)
      ALLOCATABLE GLLKE,GLLKEM,GLLKEN,GLLKE0,GLLKE0M
 
      DOUBLE COMPLEX ETAIKR(NSYMAXD,NSHELD),G(LMMAXD,LMMAXD)
      DOUBLE PRECISION BZKPK(6),KP(3)
      INTEGER NDIM
#ifdef CPP_MPI
      integer :: myMPI_comm_grid, myMPI_comm_at, myMPI_comm_ie,
     &           myrank_at, myrank_ie, nranks_at,nranks_ie
      integer :: ntot1, mytot, ii
      integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1),
     &           ntot_all(0:nranks-1), ioff_all(0:nranks-1)
      integer :: k_start, k_end
#endif


      LOGICAL TEST,OPT
!     ..
!     .. External subroutines ..
      EXTERNAL CINIT,DLKE0,OPT,TEST,GTDYSON
!     ..
!     .. Intrinsic functions ..
      INTRINSIC ATAN,EXP
!     ..
#ifdef CPP_MPI
      DOUBLE COMPLEX WORK(LMMAXD,LMMAXD,NSYMAXD)
      INTEGER IERR,IWORK
      integer :: ix,jx, mu, irec, nscoef, imin, ie
      integer, allocatable :: iatomimp(:)
      
      double precision, allocatable :: rhoq_kmask(:,:) ! only in reduced number of kpts
      DOUBLE PRECISION RECBV(3,3)
      logical, allocatable :: kmask(:) ! logical array over all kpts (determine if kpt=1,nofks is in reduced set)
#endif
!     ..

!      NDIM=LMGF0D*NAEZ
      NDIM=LMMAXD*NAEZ

      IF ( TEST('flow     ') .and. (t_inc%i_write>0))
     +     WRITE(1337,*) '>>> kkrmat1: loop over k-points'
!
      CITPI = CMI*8.D0*ATAN(1.D0)    ! = -i*2*PI 
!
      DO NS = 1,NSHELL
         DO IU = 1,NSYMAXD
            CALL CINIT(LMMAXD*LMMAXD,GS(1,1,IU,NS))
         END DO
      END DO



      LLY_GRTR = CZERO                           ! LLY Lloyd
!-----------------------------------------------------------------------
      ALLOCATE(GLLKE(ALM,ALM),STAT=IU)
      IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate GLLKE'
      IF (LLY.NE.0) ALLOCATE(DGLLKE(ALM,ALM),STAT=IU) ! LLY Lloyd
      IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate DGLLKE'
      IF (LLY.NE.0) ALLOCATE(GREFLLKE(ALM,ALM),STAT=IU) ! LLY Lloyd
      IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate GREFLLKE'

      IF ( OPT('VIRATOMS') ) THEN
         ALLOCATE (GLLKE0V(ALM,ALM),GLLKE0V2(ALM,ALM),
     +        GLLKETV(ALM,LMMAXD),GLLKETV_new(LMMAXD,ALM))
      END IF !( OPT('VIRATOMS') ) THEN


!-----------------------------------------------------------------------


      if(test('rhoqtest')) then

      if(myrank==master) then
        open(8888,file='mu0',form='formatted')
        read(8888,*) mu, nscoef
        allocate(iatomimp(nscoef))
        do i1=1,nscoef
          read(8888,*) iatomimp(i1)
        end do
        close(8888)
      end if
      
      call MPI_Bcast(mu, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'Error Bcast mu0'
      call MPI_Bcast(nscoef,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      if(ierr/=MPI_SUCCESS) stop 'Error Bcast nscoef'
      if(.not.allocated(iatomimp)) allocate(iatomimp(nscoef))
      call MPI_Bcast(iatomimp, nscoef, MPI_INTEGER, master, 
     &               MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'Error Bcast iatomimp'
      
      
      !find imin
      imin = 1000
      do i1=1,nscoef
        if(iatomimp(i1)<imin) imin = iatomimp(i1)
      end do
      nscoef = nscoef-1

        
        ! done in main1b already:
!       open(9889, access='direct', file='tau0_k', &
!      &form='unformatted', recl=(LMMAXD*LMMAXD+2)*4) ! lm blocks

      end if

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ K-points loop
#ifdef CPP_MPI
! MPI:     
      ntot1 = NOFKS

      if(myrank==master) write(1337,*) 'kkrmat k loop:', NOFKS, 
     &                          t_mpi_c_grid%nranks_ie
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, 
     &                 t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,
     &                           master,ntot1,ntot_pT,ioff_pT,.true.)

      k_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
      k_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+
     &           ntot_pT(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)

      t_mpi_c_grid%ntot_pT1 = ntot_pT
      t_mpi_c_grid%ioff_pT1 = ioff_pT
      
!       if(test('rhoqtest')) then
!          k_start = 1
!          k_end = NOFKS
!       end if
      
#else
      k_start = 1
      k_end = NOFKS 
#endif

      !print header of statusbar
      if(test('rhoqtest')) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !read recbv  
        open(8888, file='kpts.txt', form='formatted')
        read(8888,*)
        read(8888,*)
        do kpt=1,nofks
          read(8888,*) 
        end do
        read(8888,'(100E16.7)') RECBV(1:3,1:3)
        close(8888)
        
        allocate(kmask(nofks))
        
        ! use these as counters
        k_start = 1
        k_end = 0
        
        ! find kmask and number points in box
        do kpt=1,nofks
          !findig kmask
          kmask(kpt) = .false.
          do i=-1,1,1
            do j=-1,1,1
              kp(1:3) = bzkp(1:3,kpt)+i*recbv(1:3,1)+j*recbv(1:3,2)
              if(dsqrt(kp(1)**2+kp(2)**2)<0.4d0) kmask(kpt)=.true.
              if(dsqrt(kp(1)**2+kp(2)**2)<0.05d0) kmask(kpt)=.false.
            end do
          end do
          ! end finding kmask
          ! count number of kpts in reduced part
          if(kmask(kpt)) k_end = k_end+1
        end do
        
        ! fill rhoq_kmask (on reduced set of kpts)
        allocate(rhoq_kmask(1:4,k_end))
        do kpt=1, nofks
          if(kmask(kpt)) then
            rhoq_kmask(1:3,k_start) = bzkp(1:3,kpt)
            rhoq_kmask(4,k_start) = dfloat(kpt)
            k_start = k_start+1
          end if
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        ! reset k_start to 1 (was used as a counter in the loop above)
        k_start = 1
        
        
        write(*,*) 'rhow: k-loop', k_start, k_end 
        write(6,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")')
     &      0, 20, 40, 60, 80, 100
        write(6,FMT=190) !beginning of statusbar
        
      end if
      
      ! kpts loop
      DO KPT = k_start,k_end
         GLLKE(:,:) = CZERO
         IF (LLY.NE.0) DGLLKE(:,:) = CZERO
         KP(1:3) = BZKP(1:3,KPT)
         
         ! overwrite kpt in case of rhoqtest (take only reduced set of kpts)
         if(test('rhoqtest')) kp(1:3) = rhoq_kmask(1:3,kpt)
 
         ETAIKR(1:NSYMAT,1:NSHELL) = VOLCUB(KPT)
 
         ! --> first NAEZ/NATYP elements of GS() are site-diagonal
 
         DO NS = NSDIA+1,NSHELL
            I = NSH1(NS)
            J = NSH2(NS)
            DO ISYM  = 1,NSYMAT
               CARG = CZERO
               DO I1 = 1,3
                  ZKTR = RROT(ISYM,I1,NS) - RBASIS(I1,J) + RBASIS(I1,I)
                  ZKTR = KP(I1)*ZKTR
                  CARG =  CARG + ZKTR
               END DO
               ETAIKR(ISYM,NS) = ETAIKR(ISYM,NS) * EXP(CARG*CITPI)
            END DO
         END DO
 
         BZKPK(1:3) = KP(1:3)
         BZKPK(4:6) = 0.D0

         ! -> Fourier transformation
#ifdef CPP_TIMING
            if(t_inc%i_time>0) call timing_start('main1b_fourier')
#endif
 
         RRM(1:3,1:NRD) = -RR(1:3,1:NRD)
 
         ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::: KREL .EQ. 0/1
         IF (KREL.EQ.0) THEN

            !--------------------------------------------------------------------
            ALLOCATE(GLLKEN(ALMGF0,ALMGF0),STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate GLLKEN'
            GLLKEN(:,:) = CZERO
            !--------------------------------------------------------------------
            CALL DLKE0(GLLKEN,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RR,
     &                 EZOA,ATOM,BZKPK,RCLS,GINP)
            !--------------------------------------------------------------------
            ALLOCATE(GLLKEM(ALMGF0,ALMGF0),STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate GLLKEM'
            GLLKEM(:,:) = CZERO
            !--------------------------------------------------------------------
            CALL DLKE0(GLLKEM,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RRM,
     &                 EZOA,ATOM,BZKPK,RCLS,GINP)

            !--------------------------------------------------------------------
            ! LLY Lloyd
            ! Fourier for dGref/dE for Lloyds formula, repeat the above allocation 
            ! and Fourier transform for the derivatives.
            IF (LLY.NE.0) THEN
               !--------------------------------------------------------------------
               ALLOCATE(DGLLKEN(ALMGF0,ALMGF0),STAT=IU)
               IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate DGLLKEM'
               DGLLKEN(:,:) = CZERO
               !--------------------------------------------------------------------
               CALL DLKE0(DGLLKEN,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RR,
     &                    EZOA,ATOM,BZKPK,RCLS,DGINP)
               !--------------------------------------------------------------------
               ALLOCATE(DGLLKEM(ALMGF0,ALMGF0),STAT=IU)
               IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate DGLLKEM'
               DGLLKEM(:,:) = CZERO
               !--------------------------------------------------------------------
               CALL DLKE0(DGLLKEM,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RRM,
     &                    EZOA,ATOM,BZKPK,RCLS,DGINP)
            ENDIF
            ! LLY Lloyd
            !--------------------------------------------------------------------
            IF (.NOT.OPT('NEWSOSOL')) THEN

               DO I2=1,ALM
                  DO I1=1,ALM
                     GLLKE(I1,I2)= (GLLKEN(I1,I2) + GLLKEM(I2,I1))*0.5D0
                     IF (LLY.NE.0)  DGLLKE(I1,I2) =                      ! LLY Lloyd
     &                         ( DGLLKEN(I1,I2) + DGLLKEM(I2,I1) )*0.5D0 ! LLY Lloyd
                  ENDDO           
               ENDDO

            ELSE                ! (.NOT.OPT('NEWSOSOL')) 

               DO I2=1,ALMGF0
                  DO I1=1,ALMGF0
                     GLLKEN(I1,I2)=(GLLKEN(I1,I2) + GLLKEM(I2,I1))*0.5D0  
                     IF (LLY.NE.0)  DGLLKEN(I1,I2) =                      ! LLY Lloyd
     &                     ( DGLLKEN(I1,I2) + DGLLKEM(I2,I1) )*0.5D0      ! LLY Lloyd
                  ENDDO           
               ENDDO

               ! bigger GLLKE matrix and rearrange with atom block
               DO IQ1=1,NAEZ
                  DO IQ2=1,NAEZ

                     IOFF1 = LMMAXD*(IQ1-1)
                     JOFF1 = LMGF0D*(IQ1-1)
                     IOFF2 = LMMAXD*(IQ2-1)
                     JOFF2 = LMGF0D*(IQ2-1)

                     DO LM1=1,LMGF0D
                        DO LM2=1,LMGF0D
                           GLLKE(IOFF1+LM1,IOFF2+LM2) = 
     &                          GLLKEN(JOFF1+LM1,JOFF2+LM2)
                           GLLKE(IOFF1+LM1+LMGF0D,IOFF2+LM2+LMGF0D) =
     &                          GLLKEN(JOFF1+LM1,JOFF2+LM2)
                           IF (LLY.NE.0) THEN                             ! LLY Lloyd
                              DGLLKE(IOFF1+LM1,IOFF2+LM2) =               ! LLY Lloyd
     &                             DGLLKEN(JOFF1+LM1,JOFF2+LM2)           ! LLY Lloyd
                              DGLLKE(IOFF1+LM1+LMGF0D,IOFF2+LM2+LMGF0D)=  ! LLY Lloyd
     &                             DGLLKEN(JOFF1+LM1,JOFF2+LM2)           ! LLY Lloyd
                           ENDIF ! (LLY.NE.0)                             ! LLY Lloyd
                        ENDDO
                     ENDDO

                  ENDDO
               ENDDO
 
            ENDIF               ! (.NOT.OPT('NEWSOSOL'))
            !-----------------------------------------------------------------------
            DEALLOCATE(GLLKEM,STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > deallocate GLLKEM'
            DEALLOCATE(GLLKEN,STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > deallocate GLLKEN'
            IF (LLY.NE.0) DEALLOCATE(DGLLKEM,STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > deallocate DGLLKEM'
            IF (LLY.NE.0) DEALLOCATE(DGLLKEN,STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > deallocate DGLLKEN'
            !-----------------------------------------------------------------------

            ! LLY Lloyd At this point DGLLKE contains the Fourier transform of the dGref/dE

         ELSE                   !  (KREL.EQ.0) 
            ! LLY Lloyd Not implementing Lloyds formula for KREL=1 (Dirac ASA)
            !-----------------------------------------------------------------------
            ALLOCATE(GLLKE0(ALMGF0,ALMGF0),STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate GLLKE0'
            ALLOCATE(GLLKE0M(ALMGF0,ALMGF0),STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > allocate GLLKE0M'
            !-----------------------------------------------------------------------
            CALL DLKE0(GLLKE0,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RR,
     &                 EZOA,ATOM,BZKPK,RCLS,GINP)
            CALL DLKE0(GLLKE0M,ALAT,NAEZ,CLS,NACLS,NACLSMAX,RRM,
     &                 EZOA,ATOM,BZKPK,RCLS,GINP)

            DO I2=1,ALMGF0
               DO I1=1,ALMGF0
                  GLLKE0(I1,I2)=(GLLKE0(I1,I2) + GLLKE0M(I2,I1))*0.5D0
               ENDDO           
            ENDDO

            ! -> double the GLLKE0 matrix and transform to the REL representation
            !    ==> GLLKE

            ! ======================================================================
            DO IQ1=1,NAEZ
               DO IQ2=1,NAEZ   
                  IOFF1 = LMMAXD*(IQ1-1)
                  JOFF1 = LMGF0D*(IQ1-1)

                  IOFF2 = LMMAXD*(IQ2-1)
                  JOFF2 = LMGF0D*(IQ2-1)
                  ! ----------------------------------------------------------------------
                  DO IKM2 = 1,LMMAXD
                     DO IKM1 = 1,LMMAXD
     
                        CSUM1 = CZERO
                        DO IS = 1,2
                           N1 = NRREL(IS,IKM1)
                           N2 = NRREL(IS,IKM2)
                           DO I1 = 1,N1
                              J1 = IRREL(I1,IS,IKM1) + JOFF1
    
                              CSUM2 = CZERO
                              DO I2 = 1,N2
                                 J2 = IRREL(I2,IS,IKM2) + JOFF2
                                 CSUM2 = CSUM2 + 
     &                                GLLKE0(J1,J2)*SRREL(I2,IS,IKM2)
                              END DO
     
                              CSUM1 = CSUM1 + 
     &                             DCONJG(SRREL(I1,IS,IKM1))*CSUM2
                           END DO 
                        END DO
                        GLLKE(IOFF1+IKM1,IOFF2+IKM2) = CSUM1
                     END DO
                  END DO         
                  ! ----------------------------------------------------------------------
               END DO
            END DO
            !-----------------------------------------------------------------------
            DEALLOCATE(GLLKE0,GLLKE0M,STAT=IU)
            IF ( IU.NE.0 ) STOP '      <kkrmat01 > deallocate GLLKE0'
            !-----------------------------------------------------------------------
            ! ======================================================================
         END IF !  (KREL.EQ.0) 
#ifdef CPP_TIMING
            if(t_inc%i_time>0) call timing_pause('main1b_fourier')
#endif
         ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

         IF (LLY.NE.0) GREFLLKE(1:ALM,1:ALM) = GLLKE(1:ALM,1:ALM) ! LLY Save k-dependent Gref

         IF ( IDECI.EQ.1 ) CALL DECIMATE(GLLKE,NAEZ,TINVBUP,TINVBDOWN,
     &                                   VACFLAG,FACTL,NLBASIS,NRBASIS)

         ! -> Construct the matrix M=[-(t)^-1 + G^r] and store it 
         !    in the same matrix GLLKE where G^r was stored.

         IF ( .not. OPT('VIRATOMS') ) THEN

            DO I1=1,NAEZ
               DO LM1 = 1,LMMAXD
                  DO LM2 = 1,LMMAXD
                     IL1 = LMMAXD*(I1-1)+LM1
                     IL2 = LMMAXD*(I1-1)+LM2
                     GLLKE(IL1,IL2)= GLLKE(IL1,IL2) - TINVLL(LM1,LM2,I1)
                  ENDDO
               ENDDO
            ENDDO
      
#ifdef CPP_TIMING
            if(t_inc%i_time>0) call timing_start('main1b_inversion')
#endif
            !     --> Perform the inversion of matrix M
            !     the output is the scattering path operator TAU stored in GLLKE
            !     Actually -TAU, because TAU = (Deltat^-1 - Gref)^-1
            IF (LLY.NE.0) THEN ! If LLY, full inversion is needed
            CALL INVERSION(GLLKE,0,ICHECK) ! LLY
            ELSE
            CALL INVERSION(GLLKE,INVMOD,ICHECK)
            ENDIF
#ifdef CPP_TIMING
            if(t_inc%i_time>0) call timing_pause('main1b_inversion')
#endif
            ! ----------------------------------------------------------
            ! LLY Lloyd ----------------------------------------------------------
            IF (LLY.NE.0) THEN 

               ! LLY  Prepare quantities for Lloyds formula.
               ! LLY  Needed is Trace[ (1-Gref * Deltat)^-1 * d(1-Gref * Deltat)/dE ] (PhD Thiess Eq.5.38)
               ! LLY  where Deltat = t-tref. This is re-written as:
               ! LLY  -Trace[ Tau * ( dGref/dE + Gref * (dt/dE - dtref/dE) Deltat^-1 ) ]
               ! LLY  where Tau is the scattering path operator Tau = (Deltat^-1 - Gref)^-1 
               ! LLY  (negative of array GLLKE) and (t-tref)^-1 is in array TINVLL.
               ! LLY  The quantities Gref, dGref/dE, dt/dE have been prepared by main1a.
               ! LLY  Quantity dtref/dE is in array DTREFLL

               ! First set up (dt/dE - dtref/dE) Deltat^-1, store in array t_aux
               DO I1 = 1,NAEZ
                  ! GAUX1 = dt/dE-dtref/dE
                  GAUX1(1:LMMAXD,1:LMMAXD) = (1.D0/CFCTOR) *
     &              (  DTMATLL(1:LMMAXD,1:LMMAXD,I1) - 
     &                 DTREFLL(1:LMMAXD,1:LMMAXD,REFPOT(I1)) )
                  GAUX2(1:LMMAXD,1:LMMAXD) =TINVLL(1:LMMAXD,1:LMMAXD,I1)
                  ! T_AUX = (dt/dE-dtref/dE)* Deltat^-1
                  CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,
     &                    GAUX1,LMMAXD,GAUX2,LMMAXD,CZERO,GAUX3,LMMAXD)
                  T_AUX(1:LMMAXD,1:LMMAXD,I1) = GAUX3(1:LMMAXD,1:LMMAXD)
               ENDDO

               ! Now perform dGref/dE + Gref * t_aux 
               ! (Gref is ALM*ALM ; t_aux site-diagonal LMMAXD*LMMAXD)
               DO J1 = 1,NAEZ               ! Loop over columns of Gref
                  JL1 = LMMAXD*(J1-1) + 1
                  JL2 = LMMAXD*(J1-1) + LMMAXD
                  GAUX3(1:LMMAXD,1:LMMAXD) = T_AUX(1:LMMAXD,1:LMMAXD,J1)
                  DO I1 = 1,NAEZ            ! Loop over rows of Gref
                     IL1 = LMMAXD*(I1-1) + 1
                     IL2 = LMMAXD*(I1-1) + LMMAXD
                     ! Copy to small matrices
                     GAUX1(1:LMMAXD,1:LMMAXD) =GREFLLKE(IL1:IL2,JL1:JL2)
                     GAUX2(1:LMMAXD,1:LMMAXD) = DGLLKE(IL1:IL2,JL1:JL2)
                     ! GAUX2 = GAUX2 + GAUX1 * T_AUX
                     CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,GAUX1,
     &                    LMMAXD,GAUX3,LMMAXD,CONE,GAUX2,LMMAXD)
                     ! Copy back to large matrix, use again array DGLLKE 
                     ! (the I1-J1 block is not needed any more)
                     DGLLKE(IL1:IL2,JL1:JL2) = GAUX2(1:LMMAXD,1:LMMAXD)
                  ENDDO
               ENDDO

c full matrix multiple
c            ALLOCATE(GLLKE0(ALM,ALM))
c            GLLKE0=CZERO
c            DO I1=1,NAEZ
c               DO LM1 = 1,LMMAXD
c                  DO LM2 = 1,LMMAXD
c                     IL1 = LMMAXD*(I1-1)+LM1
c                     IL2 = LMMAXD*(I1-1)+LM2
c                     GLLKE0(IL1,IL2)= T_AUX(LM1,LM2,I1)
c                  ENDDO
c               ENDDO
c            ENDDO
c            CALL ZGEMM('N','N',ALM,ALM,ALM,CONE,GREFLLKE,
c     &                 ALM,GLLKE0,ALM,CONE,DGLLKE,ALM)
c            DEALLOCATE(GLLKE0)

               ! Now array DGLLKE contains 
               ! ( dGref/dE + Gref * (dt/dE - dtref/dE) Deltat^-1 )
               ! Build trace of tau * DGLLKE, -tau is conained in GLLKE.
               TRACE = CZERO
               DO IL1 = 1,ALM
                  DO IL2 = 1,ALM
                     TRACE = TRACE + GLLKE(IL1,IL2) * DGLLKE(IL2,IL1)
                  ENDDO
               ENDDO
               LLY_GRTR_K = TRACE

            ENDIF ! (LLY.NE.0)
            ! LLY Lloyd ----------------------------------------------------------
            ! ----------------------------------------------------------

         ELSE                   !  .not. OPT('VIRATOMS') 

            ! LLY Lloyd formula not built in yet for viratoms
            GLLKE0V(1:ALM,1:ALM) = GLLKE(1:ALM,1:ALM) 

            DO I1 = 1,NAEZ
               IL1 = (I1-1)*LMMAXD + 1
               
               ! GLLKETV = -GLLKE0V * TINVLL, 
               ! where TINVLL contains (t-tref) and not 1/(t-tref) in case of opt VIRATOMS
               ! tref=0 for each vir. atom.
               CALL ZGEMM('N','N',NDIM,LMMAXD,LMMAXD,-CONE,
     +             GLLKE0V(1,IL1),ALM,TINVLL(1,1,I1),LMGF0D,
     +             CZERO,GLLKETV(1,1),ALM)
               CALL ZCOPY(ALM*LMMAXD,GLLKETV(1,1),1,GLLKE0V2(1,IL1),1)
            END DO

            ! Solve (1-gt)G=g instead of [Gref - t^-1]^-1 for viratoms 
            ! because for a virtual atom t=0, t^-1 undefined.
            CALL GTDYSON(GLLKE0V2,GLLKE,NDIM,ALM,ALM)
            
         END IF                 !  .not. OPT('VIRATOMS') 
        

         ! --> global sum on array gs

         ! ======================================================================
         DO NS = 1,NSHELL
            I = NSH1(NS)
            J = NSH2(NS)
            ILM = LMMAXD*(I-1) + 1
            JLM = LMMAXD*(J-1)
     
            DO LM = 1,LMMAXD
               CALL ZCOPY(LMMAXD,GLLKE(ILM,JLM+LM),1,G(1,LM),1)
            END DO
            ! ----------------------------------------------------------------------
            DO ISYM = 1,NSYMAT
               DO LM2=1,LMMAXD
                  DO LM1=1,LMMAXD
                     GS(LM1,LM2,ISYM,NS) = GS(LM1,LM2,ISYM,NS)
     &                                   + ETAIKR(ISYM,NS) * G(LM1,LM2)
                  END DO
               END DO

               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               !>>>>>>>>>>>>>>>>>>>>>> rhoqtest >>>>>>>>>>>>>>>>>>>>>>>>
               if(test('rhoqtest')) then
!                
!                  write(*,'(9I9)') myrank, ie, nscoef, nofks, kpt, mu,  
!      &                            i,j,lmmaxd

                 irec = (nscoef*2)*(int(rhoq_kmask(4,kpt))-1) + 
     &                  (nscoef*2)*nofks*(IE-1-1)
!                  irec = (nscoef*2)*(kpt-1) + (nscoef*2)*nofks*(IE-1-1)
                 if( ((i==mu) .and. any(j==iatomimp(1:nscoef))) ) then
                   ix=0
                   jx=0
                   lm1 = 1
                   do while (ix==0 .and. lm1<=nscoef)
                     if(iatomimp(lm1)==j) ix = j - imin + 1
                     lm1 = lm1 + 1
                   end do
                   irec = irec + nscoef + ix
!                    write(*,'(A,2I9,1000000ES15.7)') 'wrt',myrank, 
!       &       irec, kp(1:2), GLLKE(ILM+1:ILM+LMMAXD,JLM+LM+1:JLM+LM+LMMAXD)!g(1:LMMAXD,1:LMMAXD)
!                    write(9889,rec=irec) G(1:LMMAXD,1:LMMAXD)!*ETAIKR(ISYM,NS)
                   write(9889,rec=irec) KP(1:2), G(1:LMMAXD,1:LMMAXD)!*ETAIKR(ISYM,NS)
!                    write(998899,'(10000ES15.7)') KP(1:2), G(1:LMMAXD,1:LMMAXD)!*ETAIKR(ISYM,NS)
                 end if
               
                 irec = (nscoef*2)*(int(rhoq_kmask(4,kpt))-1) + 
     &                  (nscoef*2)*nofks*(IE-1-1)
!                  irec = (nscoef*2)*(kpt-1) + (nscoef*2)*nofks*(IE-1-1)
                 if( ((j==mu) .and. any(i==iatomimp(1:nscoef))) ) then
                   ix=0
                   jx=0
                   lm1 = 1
                   do while (jx==0 .and. lm1<=nscoef+1)
                     if(iatomimp(lm1)==i) jx = i - imin + 1
                     lm1 = lm1 + 1
                   end do
                   irec = irec + jx
!                    write(*,'(A,2I9,1000ES15.7)') 'wrt',myrank, 
!       &                                  irec, kp(1:2), g(1:LMMAXD,1:LMMAXD)
!                    write(9889,rec=irec) g(1:LMMAXD,1:LMMAXD)!*ETAIKR(ISYM,NS)
                   write(9889,rec=irec) kp(1:2), g(1:LMMAXD,1:LMMAXD)!*ETAIKR(ISYM,NS)
!                    write(998888,'(10000ES15.7)') KP(1:2), G(1:LMMAXD,1:LMMAXD)!*ETAIKR(ISYM,NS)
                 end if ! iii==mu ...

               end if ! rhoqtest
               !<<<<<<<<<<<<<<<<<<<<<< rhoqtest <<<<<<<<<<<<<<<<<<<<<<<<
               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
            END DO ! ISYM = 1,NSYMAT
            
            ! ----------------------------------------------------------------------
            
         END DO !NS = 1,NSHELL
         ! ======================================================================
         IF (LLY.NE.0) LLY_GRTR =                                   ! LLY Lloyd Integration
     &                 LLY_GRTR + LLY_GRTR_K * VOLCUB(KPT) * NSYMAT ! LLY Lloyd Integration
         

        if(test('rhoqtest')) then
         if(k_end>=50) then
           if(mod(KPT-k_start,k_end/50)==0) write(6,FMT=200)
         else
           write(6,FMT=200)
         end if
        end if

      END DO ! KPT = 1,NOFKS   end K-points loop
         
      if(test('rhoqtest')) then
190     FORMAT('                 |'$)   ! status bar
200     FORMAT('|'$)                    ! status bar
        write(6,*)                      ! status bar
        ! finished kpts status bar
        
        write(*,*)                      ! status bar
        write(*,*) 'rhoq: write-out loop'
        write(*,'("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")')
     &      0, 20, 40, 60, 80, 100
        write(*,FMT=190) !beginning of statusbar
      
        ! write out fort.998899, fort.998888
        ! ======================================================================
        do kpt=1, nofks
          DO NS = 1,NSHELL
            I = NSH1(NS)
            J = NSH2(NS)
            DO ISYM = 1,NSYMAT
              irec = (nscoef*2)*(kpt-1) + (nscoef*2)*nofks*(IE-1-1)
              if( ((i==mu) .and. any(j==iatomimp(1:nscoef))) ) then
                ix=0
                jx=0
                lm1 = 1
                do while (ix==0 .and. lm1<=nscoef)
                  if(iatomimp(lm1)==j) ix = j - imin + 1
                  lm1 = lm1 + 1
                end do
                irec = irec + nscoef + ix 
                if(kmask(kpt)) then
                   read(9889,rec=irec) KP(1:2), G(1:LMMAXD,1:LMMAXD)
                else
                   KP(1:3) = bzkp(1:3,kpt)
                   G(1:LMMAXD,1:LMMAXD) = CZERO
                end if
                write(998899,'(10000ES15.7)') KP(1:2), 
     &             G(1:LMMAXD,1:LMMAXD)!*ETAIKR(ISYM,NS)
              end if
              irec = (nscoef*2)*(kpt-1) + (nscoef*2)*nofks*(IE-1-1)
              if( ((j==mu) .and. any(i==iatomimp(1:nscoef))) ) then
                ix=0
                jx=0
                lm1 = 1
                do while (jx==0 .and. lm1<=nscoef+1)
                  if(iatomimp(lm1)==i) jx = i - imin + 1
                  lm1 = lm1 + 1
                end do
                irec = irec + jx
                if(kmask(kpt)) then
                   read(9889,rec=irec) KP(1:2), G(1:LMMAXD,1:LMMAXD)
                else
                   KP(1:3) = bzkp(1:3,kpt)
                   G(1:LMMAXD,1:LMMAXD) = CZERO
                end if
                write(998888,'(10000ES15.7)') KP(1:2), 
     &             G(1:LMMAXD,1:LMMAXD)!*ETAIKR(ISYM,NS)
              end if ! iii==mu ...
            END DO ! ISYM = 1,NSYMAT
          END DO !NS = 1,NSHELL
          
          if(nofks>=50) then
            if(mod(KPT-0,nofks/50)==0) write(6,FMT=200)
          else
            write(6,FMT=200)
          end if
          
        end do !kpt=1,nofks
        ! ======================================================================
      
        write(6,*)                      ! status bar
      end if !test('rhoqtest')

      TRACET = CZERO
      IF (LLY.EQ.2) THEN 
      ! Add trace of (t-tref)^-1 * d(t-tref)/dE. Remember that in this case 
      ! Tr(alpha^-1 d alpha/dE) should be subtracted and 
      ! Tr(alpha_ref^-1 d alpha_ref/dE) should be added.
         DO I1 = 1,NAEZ
            GAUX1(1:LMMAXD,1:LMMAXD) = CFCTOR*!(1.D0/CFCTOR) *
     &           (  DTMATLL(1:LMMAXD,1:LMMAXD,I1) - 
     &           DTREFLL(1:LMMAXD,1:LMMAXD,REFPOT(I1)) )
            DO LM1 = 1,LMMAXD
               DO LM2 = 1,LMMAXD
                  TRACET = TRACET + GAUX1(LM1,LM2) * TINVLL(LM2,LM1,I1)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      DEALLOCATE(GLLKE,STAT=IU)
      IF ( IU.NE.0 ) STOP '      <kkrmat01 > deallocate GLLKE'
      IF (LLY.NE.0) DEALLOCATE(DGLLKE,STAT=IU)                      ! LLY Lloyd
      IF ( IU.NE.0 ) STOP '      <kkrmat01 > deallocate DGLLKE'     ! LLY Lloyd
      IF (LLY.NE.0) DEALLOCATE(GREFLLKE,STAT=IU)                    ! LLY Lloyd
      IF ( IU.NE.0 ) STOP '      <kkrmat01 > deallocate GREFLLKE'   ! LLY Lloyd
!-----------------------------------------------------------------------
! 
#ifdef CPP_MPI
        DO NS = 1,NSHELL
          IWORK = LMMAXD*LMMAXD*NSYMAXD
          CALL MPI_ALLREDUCE(GS(1,1,1,NS),WORK,IWORK,
     +                       MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                       t_mpi_c_grid%mympi_comm_ie,IERR)
          CALL ZCOPY(IWORK,WORK,1,GS(1,1,1,NS),1)
        END DO
          
        IF (LLY.NE.0) then
          IWORK = 1
          CALL MPI_ALLREDUCE(LLY_GRTR,WORK(1,1,1),IWORK,
     +                       MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                       t_mpi_c_grid%mympi_comm_ie,IERR)
          CALL ZCOPY(IWORK,WORK(1,1,1),1,LLY_GRTR,1)
        endif
          
        if(lly.eq.2) then
          IWORK = 1
          CALL MPI_ALLREDUCE(TRACET,WORK(1,1,1),IWORK,
     +                       MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                       t_mpi_c_grid%mympi_comm_ie,IERR)
          CALL ZCOPY(IWORK,WORK(1,1,1),1,TRACET,1)
        endif
#endif

      IF ( TEST('flow    ') .and. (t_inc%i_write>0)) 
     &       WRITE(1337,*) '<<< KKRMAT1'

      END
!-----------------------------------------------------------------------
      SUBROUTINE GTDYSON(GTMAT,GMAT,NDIM,LMGF0D,NGD)
! **********************************************************************
! * Solve the Dyson equation (1-g*t) * G = g                           *
! **********************************************************************

      IMPLICIT NONE
!     .. PARAMETERS ..
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
!     ..
!     .. SCALAR ARGUMENTS ..
      INTEGER NDIM,NGD,LMGF0D
!     ..
!     .. ARRAY ARGUMENTS ..
      DOUBLE COMPLEX GMAT(NGD,LMGF0D),GTMAT(NGD,NGD)
!     ..
!     .. LOCAL SCALARS ..
      INTEGER I,INFO
!     ..
!     .. LOCAL ARRAYS ..
      INTEGER IPVT(NGD)
!     ..
!     .. EXTERNAL SUBROUTINES ..
      EXTERNAL ZGETRF,ZGETRS
!     ..
 
      DO 10 I = 1,NDIM
        GTMAT(I,I) = CONE + GTMAT(I,I) ! GTMAT= 1 - G * T
   10 CONTINUE
!
!---> SOLVE THE SYSTEM OF LINEAR EQUATIONS
!
      CALL ZGETRF(NDIM,NDIM,GTMAT,NGD,IPVT,INFO)
      CALL ZGETRS('N',NDIM,LMGF0D,GTMAT,NGD,IPVT,GMAT,NGD,INFO)
      RETURN

      END
