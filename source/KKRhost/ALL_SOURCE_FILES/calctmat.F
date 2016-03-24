C*==calctmat.f    processed by SPAG 6.05Rc at 11:17 on 10 May 2004
      SUBROUTINE CALCTMAT(ICST,INS,IELAST,
     +                   NSRA,ISPIN,NSPIN,I1,EZ,
     +                   DRDI,RMESH,VINS,VISP,ZAT,IRMIN,IPAN,                ! Added IRMIN 1.7.2014
     +                   IRCUT,CLEB,LOFLM,ICLEB,IEND,SOLVER,SOCTL,CTL,
     +                   VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,
     +                   ZREL,JWSREL,IDOLDAU,LOPT,WLDAU,
     &                   LLY,DELTAE) ! LLY 
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *  LDA+U implementation     Mar. 2002-Dec.2004                      *
C *                           ph.mavropoulos, h. ebert, v. popescu    *
C * Notes:                                                            *
C *  average WLDAU for spherical wavefunctions:                       *
C *  The spherical part of the d or f wavefunction is found by adding *
C *  the average interaction potential WLDAUAV to the spherical       *
C *  potential. Then the non-spherical parts are found by using only  *
C *  the deviation of WLDAU from the average. This speeds up the      *
C *  convergence of the Born series. See also subroutines             *
C *  regsol, pnstmat and pnsqns                                       *
C *                                                                   *
C *********************************************************************
#ifdef CPP_MPI
      use mpi
#endif
      use mod_mympi, only: myrank, nranks, master
#ifdef CPP_MPI
     &                      ,distribute_linear_on_tasks
#endif
      use mod_types, only: t_tgmat,t_inc,t_mpi_c_grid,init_tgmat,
     &                     t_lloyd,init_tlloyd 
      IMPLICIT NONE
C
C     .. Parameters ..
      INCLUDE 'inc.p'
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER MMAXD
      PARAMETER ( MMAXD = 2*LMAXD+1 )
      INTEGER LM2D
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      DOUBLE PRECISION CVLIGHT
      PARAMETER (CVLIGHT=274.0720442D0)
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ZAT
      INTEGER I1,ICST,IELAST,IEND,INS,IPAN,ISPIN,NSPIN,NSRA
      INTEGER IDOLDAU,JWSREL,LOPT,ZREL
      INTEGER LLY ! LLY <> 0 for applying Lloyds formula
      INTEGER SIGNDE,IDERIV
      DOUBLE COMPLEX DELTAE  ! Energy difference for numerical derivative
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX EZ(IEMXD),
     +               TMAT0(LMMAXD,LMMAXD),
     +               TMATLL(LMMAXD,LMMAXD),
     &               DTMATLL(LMMAXD,LMMAXD) ! LLY t-matrix derivative
      DOUBLE PRECISION CLEB(NCLEB,2),DRDI(IRMD),RMESH(IRMD),
     +                 VINS(IRMIND:IRMD,LMPOTD),
     +                 VISP(IRMD)
      CHARACTER*10 SOLVER
      DOUBLE PRECISION SOCTL(KREL*LMAXD+1)
      DOUBLE PRECISION CTL(KREL*LMAXD+1)
      INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAND),
     +        LOFLM(LM2D)
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL)),
     &                 R2DRDIREL(IRMD*KREL+(1-KREL)),
     &                 RMREL(IRMD*KREL+(1-KREL))
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX ERYD,EK
      INTEGER IE,IREC,LM1,LM2,LMHI,LMLO,M1,MMAX,IRMIN
      DOUBLE PRECISION WLDAUAV
C     .. this routine does not need irregular wavefunctions
      LOGICAL LIRRSOL
      PARAMETER ( LIRRSOL = .TRUE. )
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),
     +               FZ(IRMD,0:LMAXD),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
     +               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),TMAT(0:LMAXD),
     &               ALPHALL(LMMAXD,LMMAXD),DALPHALL(LMMAXD,LMMAXD),  ! LLY alpha matrix and derivative
     &               ALPHA0(LMMAXD,LMMAXD),AUX(LMMAXD,LMMAXD),         ! LLY alpha-matrix
     &               TRALPHA,TRALPHA1  ! LLY Trace of alpha^-1 * d alpha/dE for Lloyds formula
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND)
      DOUBLE PRECISION RS(IRMD,0:LMAXD),SL(0:LMAXD)
      DOUBLE PRECISION CUTOFF(IRMD)
      INTEGER IPIV(LMMAXD) ! LLY
      CHARACTER*9 TXTS(2)
#ifdef CPP_MPI
      integer :: ie_start, mytot
      integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
#endif
      integer :: ie_end, ie_num
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
C     .. External Subroutines ..
      EXTERNAL CRADWF,PNSTMAT,WFMESH,CMATSTR,DRVRELTMAT,ZGEINV1,ZGEMM
C     ..
C     .. Data Statements
      DATA TXTS /'spin   UP','spin DOWN'/
C     ..................................................................
C ==LDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAU
C 
      IF ( IDOLDAU.EQ.1 ) THEN
         WLDAUAV = 0.D0                                        
         LMLO = LOPT*LOPT + 1
         LMHI = (LOPT+1)*(LOPT+1)
         MMAX = LMHI - LMLO + 1
         DO M1 = 1,MMAX                                        
            WLDAUAV = WLDAUAV + WLDAU(M1,M1,ISPIN)         
         ENDDO                                                 
         WLDAUAV = WLDAUAV/DBLE(MMAX)                        

C
C -> Note: Application if WLDAU makes the potential discontinuous.
C    A cutoff can be used if necessary to make the potential continuous
C    for example (array bounds should be adjusted):
C
Cccc            CUTOFF(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) *
Cccc     &                   ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
Cccc            CUTOFF(IR) = 1D0/CUTOFF(IR)
C
         DO M1 = 1,IRMD
            CUTOFF(M1) = 1.D0
         END DO
      END IF                                                    
C
C ==LDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAU
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      if(myrank==master .and. t_inc%i_write>0) WRITE(1337,*) 'atom: ',I1
#ifdef CPP_MPI
        call distribute_linear_on_tasks(t_mpi_c_grid%nranks_at, 
     &                   t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,
     &                             master,IELAST,ntot_pT,ioff_pT,.true.)
     
        ie_start = ioff_pT(t_mpi_c_grid%myrank_at)
        ie_end   = ntot_pT(t_mpi_c_grid%myrank_at)

        t_mpi_c_grid%ntot2=ie_end !t_mpi_c_grid%dims(1)
        if(.not. (allocated(t_mpi_c_grid%ntot_pT2) .or. 
     &            allocated(t_mpi_c_grid%ioff_pT2))) 
     &     allocate(t_mpi_c_grid%ntot_pT2(0:t_mpi_c_grid%nranks_at-1),
     &              t_mpi_c_grid%ioff_pT2(0:t_mpi_c_grid%nranks_at-1))
        t_mpi_c_grid%ntot_pT2 = ntot_pT
        t_mpi_c_grid%ioff_pT2 = ioff_pT
        ! now initialize arrays for tmat, gmat, and gref
        call init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
        if(lly.ne.0) call init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)
        
!         write(*,*) 'calctmatloop',myrank,i1,ie,
!      &      t_mpi_c_grid%myrank_ie,t_mpi_c_grid%myrank_at


       DO ie_num=1,ie_end
       
         IE = ie_start+ie_num
#else
      ! now initialize arrays for tmat, gmat, and gref
      call init_tgmat(t_inc,t_tgmat,t_mpi_c_grid)
      if(lly.ne.0) call init_tlloyd(t_inc,t_lloyd,t_mpi_c_grid)

      DO IE = 1,IELAST
         ie_num = IE
         ie_end = IELAST
#endif
         if(t_inc%i_write>0) WRITE(1337,*)'CALCTMAT: IE=',IE,' ATOM:',I1

         TMATLL(1:LMMAXD,1:LMMAXD) = CZERO
         DTMATLL(1:LMMAXD,1:LMMAXD) = CZERO

         ALPHALL(1:LMMAXD,1:LMMAXD) = CZERO
         DALPHALL(1:LMMAXD,1:LMMAXD) = CZERO

! In case of Lloyds formula the derivative of t is needed. 
! Then calculate t at E+dE, E-dE and average for t, subtract for dt/dE
         IDERIV = 0                         ! LLY
         IF (LLY.NE.0) IDERIV = 1           ! LLY
         DO SIGNDE = -IDERIV,IDERIV,2       ! LLY
            PZ(:,:) = CZERO
            QZ(:,:) = CZERO
            FZ(:,:) = CZERO
            SZ(:,:) = CZERO
            PNS(:,:,IRMIND:IRMD,:) = CZERO
            TMAT0(1:LMMAXD,1:LMMAXD) = CZERO
            ALPHA0(1:LMMAXD,1:LMMAXD) = CZERO

            ERYD = EZ(IE) + SIGNDE * DELTAE / 2.D0 ! LLY
            if(t_inc%i_write>0) WRITE(1337,*) 'energy:',IE,'',ERYD
C
C=======================================================================
C non/scalar-relativistic OR relativistic
C
            IF ( KREL.EQ.0 ) THEN
               CALL WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,RMESH,SL,RS,
     &                    IRCUT(IPAN),IRMD,LMAXD)
               CALL CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,SL,
     &                     PZ,FZ,QZ,SZ,TMAT,VISP,DRDI,RMESH,ZAT,LIRRSOL,
     &                     IDOLDAU,LOPT,WLDAUAV,CUTOFF)

               DO LM1 = 1,LMMAXD                      ! LLY
                  ALPHA0(LM1,LM1) = ALPHA(LOFLM(LM1)) ! LLY spherical (diag.) part of alpha
               END DO                                 ! LLY
C-----------------------------------------------------------------------
C spherical/non-spherical
C
               IF ( INS.EQ.0 ) THEN
                  DO LM1 = 1,LMMAXD
                     TMAT0(LM1,LM1) = TMAT(LOFLM(LM1))
                  END DO
               ELSE
                  CALL PNSTMAT(DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS,TMAT0,
     &                 VINS,IRMIN,IPAN,IRCUT,NSRA,CLEB,ICLEB,IEND,LOFLM,     ! Added IRMIN 1.7.2014
     &                      TMAT,LMAXD,IDOLDAU,LOPT,LMLO,LMHI,
     &                      WLDAU(1,1,ISPIN),WLDAUAV,CUTOFF,
     &                      ALPHA0)                 ! LLY In goes diag. alpha, out comes full alpha
               END IF
C
C-----------------------------------------------------------------------
C=======================================================================
!commented this out because compiling with debug options did not work
!due to an erro in the interface to drvreltmat > talk to long,
!philipp 20_08_2015
             ELSE
! C=======================================================================
               CALL DRVRELTMAT(ERYD,TMAT0,VTREL,BTREL,RMREL,
     &                      DRDIREL,R2DRDIREL,ZREL,JWSREL,SOLVER,SOCTL,
     &                      CTL,LMMAXD,LMAXD,IRMD)
            END IF
C
C non/scalar-relativistic OR relativistic
C=======================================================================
C
            ! In case of derivative calculation, average the t-matrix
            ! of the two energies and calculate the difference
            TMATLL(1:LMMAXD,1:LMMAXD) = TMATLL(1:LMMAXD,1:LMMAXD)   ! LLY
     &                                   + TMAT0(1:LMMAXD,1:LMMAXD)       ! LLY
            ALPHALL(1:LMMAXD,1:LMMAXD) =  ALPHALL(1:LMMAXD,1:LMMAXD)      ! LLY
     &                                   + ALPHA0(1:LMMAXD,1:LMMAXD)      ! LLY
            IF (LLY.NE.0) THEN
            DTMATLL(1:LMMAXD,1:LMMAXD) =DTMATLL(1:LMMAXD,1:LMMAXD)  ! LLY
     &                              + SIGNDE * TMAT0(1:LMMAXD,1:LMMAXD)   ! LLY
            DALPHALL(1:LMMAXD,1:LMMAXD) = DALPHALL(1:LMMAXD,1:LMMAXD)     ! LLY
     &                              + SIGNDE * ALPHA0(1:LMMAXD,1:LMMAXD)  ! LLY
            ENDIF

         ENDDO  ! SIGNDE = -IDERIV,IDERIV,2      ! LLY

         ! Average values of t-matrix and alpha at E+DE and E-DE
         TMATLL(1:LMMAXD,1:LMMAXD) = 
     &        TMATLL(1:LMMAXD,1:LMMAXD) / DFLOAT(1+IDERIV)  !/ 1 or 2 (IDERIV=1 or 2)! LLY
         ALPHALL(1:LMMAXD,1:LMMAXD) = 
     &        ALPHALL(1:LMMAXD,1:LMMAXD) / DFLOAT(1+IDERIV)    !/ 1 or 2 ! LLY

         IF (LLY.NE.0) THEN
            ! Construct derivative of t-matrix and aplha
            DTMATLL(1:LMMAXD,1:LMMAXD) = 
     &           DTMATLL(1:LMMAXD,1:LMMAXD) / DELTAE
            DALPHALL(1:LMMAXD,1:LMMAXD) = 
     &           DALPHALL(1:LMMAXD,1:LMMAXD) / DELTAE
            ! LLY Calculate Trace [alpha^-1 * d alpha/dE] for Lloyd's formula
            ALPHA0(1:LMMAXD,1:LMMAXD) = CZERO
            AUX(1:LMMAXD,1:LMMAXD) = CZERO
            ! ALPHA0 = ALPHALL^-1
            CALL ZGEINV1(ALPHALL,ALPHA0,AUX,IPIV,LMMAXD)
            ! AUX = ALPHALL^-1 * DALPHALL
            CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,ALPHA0,LMMAXD,
     &              DALPHALL,LMMAXD,CZERO,AUX,LMMAXD)
            ! Trace of AUX is Trace of [alpha^-1 * d alpha/dE]
            TRALPHA = CZERO
            DO LM1 = 1,LMMAXD
               TRALPHA = TRALPHA + AUX(LM1,LM1)
            ENDDO

            IF (LLY.EQ.-1) THEN ! Calculate Tr(t^-1 dt/dE) instead of Tr(alpha^-1 * d alpha/dE)
               ALPHA0(1:LMMAXD,1:LMMAXD) = TMATLL(1:LMMAXD,1:LMMAXD)
               TMAT0(1:LMMAXD,1:LMMAXD) = CZERO
               AUX(1:LMMAXD,1:LMMAXD) = CZERO
               ! ALPHA0 = TMATLL^-1
               CALL ZGEINV1(ALPHA0,TMAT0,AUX,IPIV,LMMAXD)
               ! Build trace
               TRALPHA1 = CZERO
               DO LM1 = 1,LMMAXD
                  DO LM2 = 1,LMMAXD
                     TRALPHA1 = TRALPHA1 + 
     &                          TMAT0(LM1,LM2) * DTMATLL(LM2,LM1)
                  ENDDO
               ENDDO
               TRALPHA1 = CZERO
               do lm1 = 1,lmmaxd
                  TRALPHA1 = TRALPHA1 + 
     &                       dtmatll(lm1,lm1)/tmatll(lm1,lm1)
               enddo

            ENDIF
            
         ENDIF ! (LLY.NE.0) 


         TMAT0(1:LMMAXD,1:LMMAXD) = TMATLL(1:LMMAXD,1:LMMAXD)
         IREC = IE + IELAST*(ISPIN-1) + IELAST*NSPIN* (I1-1)
         if (t_tgmat%tmat_to_file) then
            WRITE(69,REC=IREC) TMAT0
         else
#ifdef CPP_MPI
            IREC = ie_num + ie_end*(ISPIN-1) + ie_end*NSPIN*
     *             (I1-t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie)-1)
#else
            IREC = ie_num + ie_end*(ISPIN-1) + ie_end*NSPIN*I1
#endif
            t_tgmat%tmat(:,:,irec) = tmat0
         end if
         IF (LLY.NE.0) THEN                                          ! LLY
            if(t_lloyd%dtmat_to_file) then
               IREC = IE + IELAST*(ISPIN-1) + IELAST*NSPIN* (I1-1)
               WRITE(691,REC=IREC) DTMATLL                                ! LLY
            else
               irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(i1-1)
               t_lloyd%dtmat(:,:,irec) = DTMATLL(:,:)
            end if ! t_lloyd%dtmat_to_file
            if(t_lloyd%tralpha_to_file) then
               IREC = IE + IELAST*(ISPIN-1) + IELAST*NSPIN* (I1-1)
               WRITE(692,REC=IREC) TRALPHA                              ! LLY
            else
               irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(i1-1)
               t_lloyd%tralpha(irec) = TRALPHA
            end if
         ENDIF                                                       ! LLY
         
C ----------------------------------------------------------------------
         IF ( TEST('tmat    ') .and. (t_inc%i_write>0)) THEN
            WRITE (1337,*)
            WRITE (1337,99001) '-----> t matrix for atom: ',I1
            IF ( KREL.EQ.0 ) WRITE (1337,99002) TXTS(ISPIN)
            WRITE (1337,99003) ', energy: ',ERYD
            CALL CMATSTR(' ',1,TMAT0,LMMAXD,LMMAXD,
     &                   2*KREL+1,2*KREL+1,0,1D-8,6)
            WRITE (1337,*)
         END IF
C ----------------------------------------------------------------------
      END DO ! IE = 1,IELAST
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
99001 FORMAT (A,I3,$)
99002 FORMAT (', ',A,$)
99003 FORMAT (A,2F10.6)
C
      RETURN
      END
