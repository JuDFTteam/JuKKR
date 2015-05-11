C*==rhoval.f    processed by SPAG 6.05Rc at 11:22 on 10 May 2004
      SUBROUTINE RHOVAL(IHOST,LDORHOEF,ICST,INS,IELAST,NSRA,
     &                  ISPIN,NSPIN,NSPINPOT,
     &                  I1,EZ,WEZ,DRDI,R,VINS,VISP,ZAT,IPAN,IRCUT,IRMIN,
     &                  THETAS,IFUNM,LMSP,RHO2NS,R2NEF,RHOORB,DEN,DENLM,
     &                  MUORB,ESPV,CLEB,LOFLM,ICLEB,IEND,JEND,SOLVER,
     &                  SOCTL,CTL,VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,
     &                  ZREL,JWSREL,IRSHIFT,ITERMVDIR,QMTET,QMPHI,
     &                  MVEVIL,MVEVILEF,NMVECMAX,
     &                  IDOLDAU,LOPT,PHILDAU,WLDAU,DENMATC,
     &                  LLY,NATYP,NQDOS)          ! LLY Lloyd
C
C **********************************************************************
C * For KREL = 1 (relativistic mode)                                   *
C *                                                                    *
C *  NPOTD = 2 * NATYPD                                                *
C *  LMMAXD = 2 * (LMAXD+1)^2                                          *
C *  NSPIND = 1                                                        *
C *                                                                    *
C *  LDA+U implementation     Mar. 2002-Dec.2004                       *
C *                           ph.mavropoulos, h. ebert, v. popescu     *
C * Notes:                                                             *
C *  average WLDAU for spherical wavefunctions:                        *
C *  The spherical part of the d or f wavefunction is found by adding  *
C *  the average interaction potential WLDAUAV to the spherical        *
C *  potential. Then the non-spherical parts are found by using only   *
C *  the deviation of WLDAU from the average. This speeds up the       *
C *  convergence of the Born series. See also subroutines              *
C *  regsol, pnstmat and pnsqns                                        *
C *                                                                    *
C **********************************************************************
#ifdef CPP_MPI
      use mpi
#endif
      use mod_mympi, only: myrank, nranks, master,
     &                      distribute_linear_on_tasks
      use mod_types, only: t_tgmat,t_inc,t_mpi_c_grid, init_tgmat
      
      IMPLICIT NONE
C
C     .. Parameters ..
      INCLUDE 'inc.p'
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
      INTEGER LMAXD1
      PARAMETER (LMAXD1= LMAXD+1)
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
      DOUBLE COMPLEX CONE 
      PARAMETER ( CONE=(1D0,0D0) )
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ZAT
      INTEGER I1,ICST,IELAST,IEND,INS,IPAN,ISPIN,NSPIN,NSPINPOT,NSRA
      INTEGER JWSREL,ZREL,IRSHIFT,IDOLDAU,LOPT,NATYP
      INTEGER LLY,IRMAX1,IELAST1,NSPIN1,NATYP1 ! LLY Lloyd
      INTEGER IHOST,NMVECMAX,IRMIN,NQDOS
      LOGICAL LDORHOEF
C
C     IHOST = 1   < -- this routine is called by the HOST tbkkr-program
C     IHOST <> 1  < --                 called by the IMPURITY program
C
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DEN(0:LMAXD1,IEMXD*(1+KREL),NQDOS),EZ(IEMXD),
     +               WEZ(IEMXD),DENMATC(MMAXD,MMAXD),
     +               DENLM(LMMAXD,IEMXD*(1+KREL),NQDOS),
     &               DUM_DENLM(LMMAXD)
C     .. first 2 indices in dmuorb are the spin-resolved contributions,
C     .. the 3rd one should be the sum of them 
      DOUBLE COMPLEX DMUORB(0:KREL*LMAXD+(1-KREL),3)
      DOUBLE PRECISION MUORB(0:LMAXD1+1,3)
      DOUBLE COMPLEX PHILDAU(IRMD)
      DOUBLE PRECISION CLEB(NCLEB,2),DRDI(IRMD),
     +                 ESPV(0:LMAXD1,2),    ! changed for REL case
     +                 R(IRMD),RHO2NS(IRMD,LMPOTD,2),
     +                 R2NEF(IRMD,LMPOTD,2),   ! at fermi energy
     +                 THETAS(IRID,NFUND),VINS(IRMIND:IRMD,LMPOTD),
     +                 VISP(IRMD)
      DOUBLE PRECISION RHOORB(IRMD*KREL + (1-KREL))
      DOUBLE PRECISION SOCTL(KREL*LMAXD+1),CTL(KREL*LMAXD+1)
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL)),
     &                 R2DRDIREL(IRMD*KREL+(1-KREL)),
     &                 RMREL(IRMD*KREL+(1-KREL))
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND)
      CHARACTER*10 SOLVER
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C      ITERMDIR variables
C
      LOGICAL ITERMVDIR
      DOUBLE PRECISION QMTET,QMPHI
      COMPLEX*16 MVEVIL(0:LMAXD,3,NMVECMAX) ! OUTPUT
      COMPLEX*16 MVEVILEF(0:LMAXD,3,NMVECMAX) ! OUTPUT
C
C      ITERMDIR variables
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      INTEGER ICLEB(NCLEB,4),IFUNM(LMXSPD),IRCUT(0:IPAND),
     +        JEND(LMPOTD,0:LMAXD,0:LMAXD),
     +        LMSP(LMXSPD),LOFLM(LM2D)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION WLDAUAV
      DOUBLE COMPLEX DF,ERYD,EK
      DOUBLE PRECISION PI
      DOUBLE COMPLEX DENTOT ! qdos
      INTEGER IDIM,IE,IR,L,LM1,LM2,LMHI,LMLO,IREC,ISPINPOT,
     &        LASTEZ,M1,MMAX
      INTEGER IQ!,NQDOS ! qdos number of qdos points
      INTEGER IX       ! qdos
      INTEGER LRECGFLLE,LMSHIFT1(4),LMSHIFT2(4),JSPIN,L1 ! lmlm-dos

C     .. this routine needs irregular wavefunctions
      LOGICAL LIRRSOL
      PARAMETER ( LIRRSOL = .TRUE. )
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),
     +               DR(LMMAXD,LMMAXD),
     +               EKL(0:LMAXD),FZ(IRMD,0:LMAXD),
     +               GMATLL(LMMAXD,LMMAXD,IEMXD),
     +               GMAT0(LMMAXD,LMMAXD),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD),
     +               SZ(IRMD,0:LMAXD),TMAT(0:LMAXD)
      DOUBLE PRECISION RS(IRMD,0:LMAXD),S(0:LMAXD)
      DOUBLE PRECISION CUTOFF(IRMD)
      DOUBLE COMPLEX DENDUM(0:LMAXD1),GFLLE(:,:,:,:),DENTOT2(2),
     +               DENTMP2(LMAXD,2),DUM_GFLLE(:,:)
      DOUBLE PRECISION QVEC(:,:)       ! qdos, q-vectors for qdos
      ALLOCATABLE QVEC,GFLLE,DUM_GFLLE      ! qdos, lmlm-dos
#ifdef CPP_MPI
      integer :: ie_start, ie_end, mytot, ie_num
      integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
#endif
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DSCAL,CRADWF,DRVRHO,PNSQNS,RHOLM,
     +         RHONS,WFMESH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG,SQRT
C     ..
C     .. Arrays in Common ..
      CHARACTER*8 OPTC(32),TESTC(32)             ! qdos
C     ..
C     .. Common blocks ..
      COMMON /OPTC/OPTC                         ! qdos
      COMMON /TESTC/TESTC                       ! qdos
C     .. External Functions ..
      LOGICAL OPT,TEST                          ! qdos
      EXTERNAL OPT,TEST                         ! qdos

CMPI  INTEGER MYRANK,NROFNODES
CMPI  INTEGER MAPBLOCK
CMPI  EXTERNAL MAPBLOCK
CMPI  COMMON /MPI/MYRANK,NROFNODES

      PI = 4.D0 * DATAN(1.D0)                   ! qdos

C     ..
C     ..................................................................
C
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
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C initialise variables 
C
C ======================================================================
      IF ( KREL.EQ.0 ) THEN
         DO LM1 = 1,LMPOTD
            DO IR = 1,IRMD
               RHO2NS(IR,LM1,ISPIN) = 0.0D0
            END DO
         END DO
C     
         DO L = 0,LMAXD1
            ESPV(L,ISPIN) = 0.0D0
         END DO
C ======================================================================
      ELSE
C ======================================================================
         DO ISPINPOT = 1,2
            DO LM1 = 1,LMPOTD
               DO IR = 1,IRMD
                  RHO2NS(IR,LM1,ISPINPOT) = 0.0D0
                  R2NEF(IR,LM1,ISPINPOT)  = 0.0D0
               END DO
            END DO
C     
            DO L = 0,LMAXD1
               ESPV(L,ISPINPOT) = 0.0D0
            END DO
C     
         END DO
C     
         DO IR = 1,IRMD
            RHOORB(IR) = 0.0D0
         END DO
C     
         DO IR = 1,3 
            DO L = 0,LMAXD
               DMUORB(L,IR) = (0.0D0,0.0D0)
            END DO
C
            DO L = 0,LMAXD1 + 1
               MUORB(L,IR) = 0.0D0
            END DO
         END DO
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR
C
         IF (ITERMVDIR) THEN
            DO LM1 = 1,3
               DO LM2 = 1, NMVECMAX
                  DO L = 0, LMAXD
                     MVEVIL(L,LM1,LM2) = (0.0D0,0.0D0)
                     MVEVILEF(L,LM1,LM2) = (0.0D0,0.0D0)
                  END DO
               END DO
            END DO
         END IF
C     
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END IF                    ! KREL = 0/1
C ======================================================================
      LASTEZ = IELAST
C
C end initialise variables 
C IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
C
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
#ifndef CPP_MPI
      IF (OPT('qdos    ')) THEN
      OPEN(31,   ! lm-dos
     +FILE="qdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//   ! lm-dos
     +char(48+ISPIN)//".dat")   ! lm-dos
      ENDIF
      OPEN(30,   ! lm-dos
     +FILE="lmdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//   ! lm-dos
     +char(48+ISPIN)//".dat")   ! lm-dos
      WRITE (30,*) ' '   ! lm-dos
      WRITE (30,8600) '# ISPIN=',ISPIN,' I1=',I1   ! lm-dos
 8600 FORMAT (a8,I3,a4,I5)  ! lm-dos
#endif

c set LMSHIFT value which is need to construct dentmp
       LMSHIFT1(1)=0
       LMSHIFT1(2)=LMMAXD
       LMSHIFT1(3)=0
       LMSHIFT1(4)=LMMAXD
       LMSHIFT2(1)=0
       LMSHIFT2(2)=LMMAXD
       LMSHIFT2(3)=LMMAXD
       LMSHIFT2(4)=0

      NQDOS = 1                                         ! qdos 
      IF (OPT('qdos    ')) THEN                         ! qdos
C        Read BZ path for qdos calculation:
         OPEN(67,FILE='qvec.dat')                       ! qdos
         READ(67,*) NQDOS                               ! qdos
         ALLOCATE(QVEC(3,NQDOS))                        ! qdos
         DO IQ = 1,NQDOS                                ! qdos
            READ(67,*) (QVEC(IX,IQ),IX=1,3)             ! qdos
         ENDDO                                          ! qdos
         CLOSE(67)                                      ! qdos
      END IF

      ALLOCATE(GFLLE(LMMAXD,LMMAXD,IELAST,NQDOS))
      ALLOCATE(DUM_GFLLE(LMMAXD,LMMAXD))

      if(myrank==master) WRITE(*,*) 'atom',I1
#ifdef CPP_MPI
      ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_col)
      ie_end = t_mpi_c_grid%ntot2
!       write(*,*) myrank,ie_start,ie_end

      DO ie_num=1,ie_end
       
         IE = ie_start+ie_num
#else
      DO IE = 1,IELAST
#endif
         WRITE(*,*) 'energy',IE,EZ(IE)
CMPI     IF(MYRANK.EQ.MAPBLOCK(IE,1,IELAST,1,0,NROFNODES-1)) THEN
C
!          IREC = IE + IELAST*(ISPIN-1) + IELAST*NSPIN* (I1-1)
!          READ(69,REC=IREC) GMAT0
! C
!          DO LM2 = 1,LMMAXD
!             DO LM1 = 1,LMMAXD
!                GMATLL(LM1,LM2,IE) = GMAT0(LM1,LM2)
!             END DO
!          END DO
C
         ERYD = EZ(IE)
         DF = WEZ(IE)/DBLE(NSPINPOT)
C
C=======================================================================
C non/scalar-relativistic OR relativistic
C
         IF ( KREL.EQ.0 ) THEN
            CALL WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,R,S,RS,IRCUT(IPAN),
     &                  IRMD,LMAXD)
            CALL CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,S,
     &                  PZ,FZ,QZ,SZ,TMAT,VISP,DRDI,R,ZAT,LIRRSOL,
     &                  IDOLDAU,LOPT,WLDAUAV,CUTOFF)
C-----------------------------------------------------------------------
C non-spherical
C
            IF (INS.GT.0) CALL PNSQNS(AR,CR,DR,DRDI,EK,ICST,PZ,QZ,FZ,SZ,
     &                               PNS,QNS,NSRA,VINS,IPAN,IRMIN,IRCUT, ! Added IRMIN 1.7.2014
     &                                CLEB,ICLEB,IEND,LOFLM,LMAXD,
     &                                IDOLDAU,LOPT,LMLO,LMHI,
     &                                WLDAU(1,1,ISPIN),WLDAUAV,CUTOFF)
C
            DO L = 0,LMAXD
               EKL(L) = EK*DBLE(2*L+1)
            END DO

            DO 200 IQ = 1,NQDOS                                       ! qdos
C-----------------------------------------------------------------------
C Read in Green function
            IREC = IQ + NQDOS * (IE-1) + NQDOS * IELAST * (ISPIN-1) + ! qdos (without qdos, IQ=NQDOS=1)
     &                                NQDOS * IELAST * NSPIN * (I1-1) ! qdos 
            if (t_tgmat%gmat_to_file) then
               READ(69,REC=IREC) GMAT0
            else
               IREC = IQ + NQDOS * (ie_num-1) + NQDOS * 
     &                t_mpi_c_grid%ntot2 * (ISPIN-1) + 
     &                NQDOS * t_mpi_c_grid%ntot2 * NSPIN * (I1-1)
               GMAT0(:,:) = t_tgmat%gmat(:,:,irec)
            end if
            IF (TEST('GMAT=0  ')) THEN
               WRITE(*,*) 'TEST GMAT=0, setting GMAT to zero'
               GMAT0 = (0.D0,0.D0) 
            ENDIF
C-----------------------------------------------------------------------
C spherical/non-spherical input potential
C
            IF ( INS.EQ.0 ) THEN
               CALL RHOLM(DEN(0,IE,IQ),DF,GMAT0,NSRA,
     +              RHO2NS(1,1,ISPIN),DRDI,IPAN,IRCUT,PZ,FZ,QZ,SZ,  
     +              CLEB(1,1),ICLEB,IEND,JEND,EKL,DENLM(1,IE,IQ),EK)
            ELSE
               CALL RHONS(DEN(0,IE,IQ),DF,DRDI,GMAT0,EK,
     +             RHO2NS(1,1,ISPIN),IPAN,IRCUT,IRMIN,THETAS,IFUNM,LMSP,  ! Added IRMIN 1.7.2014
     +              NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,
     +              JEND,IEND,EKL,DENLM(1,IE,IQ),GFLLE(:,:,IE,IQ))
            END IF


#ifndef CPP_MPI
c Write out qdos:
            IF (OPT('qdos    ')) THEN
               DENTOT = DCMPLX(0.D0,0.D0)
               DO L = 0,LMAXD1
                  DENTOT = DENTOT + DEN(L,IE,IQ)
               ENDDO
               WRITE(30,9000) ERYD,QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),      ! lmdos
     &         -DIMAG(DENTOT)/PI,(-DIMAG(DENLM(L,IE,IQ))/PI,L=1,LMMAXD)      ! lmdos
 9000          FORMAT(5F10.6,40E16.8)
               WRITE(31,9000) ERYD,QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),      ! lmdos
     &         -DIMAG(DENTOT)/PI,(-DIMAG(DENLM(L,IE,IQ))/PI,L=1,LMMAXD)      ! lmdos
            ENDIF
#endif

 200        END DO                                                   ! qdos
C
C-----------------------------------------------------------------------
            DO L = 0,LMAXD1
               ESPV(L,ISPIN) = ESPV(L,ISPIN)+DIMAG(ERYD*DEN(L,IE,1)*DF)
            END DO
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     get the charge at the Fermi energy (IELAST)
C     call RHOLM/RHONS with the energy weight CONE --> not overwrite DF
C                      with the dummy DENDUM       --> not overwrite DEN
C
            IF ( (IE.EQ.IELAST) .AND. LDORHOEF ) THEN
               IF (INS.EQ.0) THEN
                  CALL RHOLM(DENDUM,CONE,GMAT0,NSRA,
     +                 R2NEF(1,1,ISPIN),DRDI,IPAN,IRCUT,PZ,FZ,QZ,SZ,
     +                 CLEB(1,1),ICLEB,IEND,JEND,EKL)
               ELSE
                 CALL RHONS(DENDUM,CONE,DRDI,GMAT0,EK,
     +              R2NEF(1,1,ISPIN),IPAN,IRCUT,IRMIN,THETAS,IFUNM,LMSP,  ! Added IRMIN 1.7.2014 
     +                 NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,
     +                 JEND,IEND,EKL,DUM_DENLM,DUM_GFLLE)
               END IF

            END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C=======================================================================
         ELSE ! ( KREL.EQ.0 )
C=======================================================================
            CALL DRVRHO_QDOS(LDORHOEF,RHO2NS,R2NEF,DEN,DMUORB,RHOORB,
     &                  IE,ERYD,DF,LASTEZ,
     &                  GMATLL,VTREL,BTREL,RMREL,DRDIREL,
     &                  R2DRDIREL,ZREL,JWSREL,IRSHIFT,SOLVER,SOCTL,CTL,
     &                  QMTET,QMPHI,ITERMVDIR,MVEVIL,MVEVILEF,LMMAXD,
     &                  LMAXD,IRMD,LMPOTD,IEMXD,NMVECMAX,
     &                  I1,QVEC,NQDOS)            ! qdos
C
            DO L = 0,LMAXD1
               ESPV(L,1) = ESPV(L,1) + DIMAG(ERYD*DEN(L,IE,IQ)*DF)
               ESPV(L,2) = ESPV(L,2) + DIMAG(ERYD*DEN(L,IE+IEMXD,IQ)*DF)
            END DO
C
            DO IR = 1,3
               DO L = 0,LMAXD
                  MUORB(L,IR) = MUORB(L,IR) + DIMAG(DMUORB(L,IR)*DF)
               END DO
            END DO
         END IF
C
C non/scalar-relativistic OR relativistic
C=======================================================================
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C LDA+U calculation 
         IF ( ( IDOLDAU.EQ.1 ).AND.( LOPT.GE.0 ) ) 
     &        CALL DENSITYMAT(DF,PZ,QZ,PNS,QNS,AR,CR,DR,GMATLL(1,1,IE),
     &                        IPAN,IRCUT,DRDI,EK,
     &                        IRMIN,LOPT,MMAX,LMLO,LMHI,PHILDAU,DENMATC
     &        ,den,ie) ! test fivos 19.9.08
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
CMPI     END IF




      END DO                    ! IE = 1,IELAST
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
c Write out gflle
            IF (OPT('lmlm-dos')) THEN                                     ! lmlm-dos
              IF (ISPIN.EQ.1) THEN                                        ! lmlm-dos
              LRECGFLLE = WLENGTH*4*LMMAXD*LMMAXD*IELAST*NQDOS  !4 words = 16 bytes / complex number
              open(96,ACCESS='direct',RECL=LRECGFLLE,FILE="gflle",        ! lmlm-dos
     &                FORM='unformatted')                                 ! lmlm-dos
              ENDIF                                                       ! lmlm-dos
              IREC = I1 + NATYPD * (ISPIN-1)                              ! lmlm-dos
              WRITE(96,REC=IREC) GFLLE(:,:,:,:)                           ! lmlm-dos
            ENDIF



      DEALLOCATE( GFLLE, DUM_GFLLE )

C        
      IF ( IHOST.NE.1) RETURN
C
! Transformation of ISPIN=1,2 from (spin-down,spin-up) to (charge-density,spin-density)
      IF (ISPIN.EQ.2) THEN
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


      END IF


C
C Test qdos/gflle
!       IF (OPT('qdos    ')) THEN
!        DO IQ = 1,NQDOS                                                    ! qdos
!         DO IE=1,IELAST                                                    ! qdos
!          IF ((IQ.EQ.1).AND.(IE.EQ.1)) THEN                                ! qdos
!             OPEN(33,                                                      ! qdos
!      +        FILE="qdos."//char(48+I1/10)//char(48+mod(I1,10))//         ! qdos
!      +        "."//char(48+ISPIN)//".dat")                                ! qdos
!             WRITE (33,*) ' '                                              ! qdos
!             WRITE (33,8600) '# ISPIN=',ISPIN,' I1=',I1                    ! qdos
!          ENDIF   ! IQ.EQ.1                                                ! qdos
!             DENTOT2(ISPIN) = DCMPLX(0.D0,0.D0)                            ! qdos
!             DENTMP2(:,ISPIN) = DCMPLX(0.D0,0.D0)                          ! qdos
!               DO L1 = 0,LMAXD                                             ! qdos
!                 DO M1 = -L1,L1                                            ! qdos
!                   LM1 = L1*(L1+1)+M1+1                                    ! qdos
!                     DENTMP2(L1,ISPIN) = DENTMP2(L1,ISPIN) +               ! qdos
!      &              GFLLE(LM1+LMSHIFT1(ISPIN),LM1+LMSHIFT2(ISPIN),IE,IQ)  ! qdos
!                 ENDDO                                                     ! qdos
!                 DENTOT2(ISPIN) = DENTOT2(ISPIN) + DENTMP2(L1,ISPIN)       ! qdos
!               ENDDO                                                       ! qdos
! c    write qdos.nn.s.dat                                                  ! qdos
!             WRITE(33,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),       ! qdos
!      &       -DIMAG(DENTOT2(ISPIN))/PI,(-DIMAG(GFLLE(LM1+LMSHIFT1(ISPIN)
!      &                     ,LM1+LMSHIFT2(ISPIN),IE,IQ))/PI,LM1=1,LMMAXD)  ! qdos
!         ENDDO ! IE
!        ENDDO ! IQ
!       ENDIF ! qdos
      END
