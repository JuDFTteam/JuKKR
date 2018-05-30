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
     &                  NATYP,NQDOS)
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
      use mod_mympi, only: myrank, master
#ifdef CPP_MPI
     &                     ,distribute_linear_on_tasks
#endif
      use mod_types, only: t_tgmat,t_inc,t_mpi_c_grid, init_tgmat
      use mod_version_info
      
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
      DOUBLE COMPLEX CONE,CZERO,CI 
      PARAMETER ( CONE=(1D0,0D0),CZERO=(0d0,0d0),CI=(0d0,1d0) )
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ZAT
      INTEGER I1,ICST,IELAST,IEND,INS,IPAN,ISPIN,NSPIN,NSPINPOT,NSRA
      INTEGER JWSREL,ZREL,IRSHIFT,IDOLDAU,LOPT,NATYP
      INTEGER IHOST,NMVECMAX,IRMIN,NQDOS
      LOGICAL LDORHOEF
C
C     IHOST = 1   < -- this routine is called by the HOST tbkkr-program
C     IHOST <> 1  < --                 called by the IMPURITY program
C
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DEN(0:LMAXD1,IELAST*(1+KREL),NQDOS),EZ(IEMXD),
     +               WEZ(IEMXD),DENMATC(MMAXD,MMAXD),
     +               DENLM(LMMAXD,IELAST*(1+KREL),NQDOS),
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
      
!             DOUBLE PRECISION VT(NRMAX,NTMAX),BT(NRMAX,NTMAX)
!       DOUBLE PRECISION R(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX)
!       DOUBLE PRECISION DRDI(NRMAX,NMMAX),SOCTL(NTMAX,NLMAX)
!       DOUBLE PRECISION CTL(NTMAX,NLMAX)
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
#ifndef CPP_MPI
      DOUBLE COMPLEX DENTOT ! qdos
#endif
      INTEGER IDIM,IE,IR,L,LM1,LM2,LMHI,LMLO,IREC,ISPINPOT,
     &        LASTEZ,M1,MMAX
      INTEGER IQ!,NQDOS ! qdos number of qdos points
      INTEGER IX       ! qdos
      INTEGER LRECGFLLE,LMSHIFT1(4),LMSHIFT2(4) ! lmlm-dos

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
      DOUBLE COMPLEX DENDUM(0:LMAXD1),GFLLE(:,:,:,:),DUM_GFLLE(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: GLDAU(:,:)
      DOUBLE PRECISION QVEC(:,:)       ! qdos, q-vectors for qdos
      ALLOCATABLE QVEC,GFLLE,DUM_GFLLE      ! qdos, lmlm-dos
#ifdef CPP_MPI
      integer :: ie_start
#endif
      integer :: ie_end, ie_num
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DSCAL,CRADWF,DRVRHO,PNSQNS,RHOLM,
     +         RHONS,WFMESH
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG,SQRT
C     ..
C     .. External Functions ..
      LOGICAL OPT,TEST                          ! qdos
      EXTERNAL OPT,TEST                         ! qdos


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
C           DO L = 0,LMAXD1 + 1       !zimmer: initialization shifted to main1c
C              MUORB(L,IR) = 0.0D0
C           END DO
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
      IF (OPT('qdos    ')) THEN                                       ! qdos
         if(NATYP.ge.100) then                                        ! qdos
           OPEN(31,                                                   ! qdos
     +     FILE="qdos."//char(48+I1/100)//char(48+mod(I1/10,10))//    ! qdos
     +     char(48+mod(I1,10))//"."//char(48+ISPIN)//".dat")          ! qdos
         else                                                         ! qdos
           OPEN(31,                                                   ! qdos
     +     FILE="qdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//   ! qdos
     +     char(48+ISPIN)//".dat")                                    ! qdos
         end if                                                       ! qdos
         call version_print_header(31)                                ! qdos
         WRITE (31,'(7(A,3X))') '#   Re(E)','Im(E)',                  ! qdos
     &                  'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'     ! qdos
      endif                                                           ! qdos

      OPEN(30,                                                        ! lmdos
     +FILE="lmdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//       ! lmdos
     +char(48+ISPIN)//".dat")                                         ! lmdos
      call version_print_header(31)                                   ! lmdos
      WRITE (30,*) ' '                                                ! lmdos
      WRITE (30,8600) '# ISPIN=',ISPIN,' I1=',I1                      ! lmdos
 8600 FORMAT (a8,I3,a4,I5)                                            ! lmdos
      
      ! write out complex qdos for interpolation to the real axis     ! complex qdos
      if(test('compqdos')) then                                       ! complex qdos
         if(NATYP.ge.100) then                                        ! complex qdos
           OPEN(3031,                                                 ! complex qdos
     +     FILE="cqdos."//char(48+I1/100)//char(48+mod(I1/10,10))//   ! complex qdos
     +     char(48+mod(I1,10))//"."//char(48+ISPIN)//".dat")          ! complex qdos
         else                                                         ! complex qdos
           OPEN(3031,                                                 ! complex qdos
     +     FILE="cqdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//  ! complex qdos
     +     char(48+ISPIN)//".dat")                                    ! complex qdos
         end if                                                       ! complex qdos
         call version_print_header(3031)                              ! complex qdos
         WRITE (3031,'(A)') '# lmax, natyp, nspin, nqdos, ielast:'    ! complex qdos
         WRITE (3031,'(5I9)') lmaxd, natyp, nspin, nqdos, ielast      ! complex qdos
         WRITE (3031,'(7(A,3X))') '#   Re(E)','Im(E)',                ! complex qdos
     &                  'k_x','k_y','k_z','DEN_tot','DEN_s,p,...'     ! complex qdos
      ENDIF                                                           ! qdos
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
C        Read BZ path for qdos calculation:             ! qdos
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
      IF (IDOLDAU.EQ.1) THEN
       ALLOCATE(GLDAU(LMMAXD,LMMAXD))
       GLDAU=CZERO
      ENDIF

      if(myrank==master) WRITE(1337,*) 'atom',I1
#ifdef CPP_MPI
      ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
      ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)

      DO ie_num=1,ie_end
         IE = ie_start+ie_num
#else
! ! omp: start parallel region here
! !$omp parallel do default(none)
! !$omp& private(eryd,ie,ir,irec,lm1,lm2)
! !$omp& private(jlk_index,tmatll,ith)
! !$omp& shared(nspin,nsra,iend,ipot,ielast,npan_tot,ncheb,lmaxd)
! !$omp& shared(zat,socscale,ez,rmesh,cleb,rnew,nth,icleb,thetasnew,i1)
! !$omp& shared(rpan_intervall,vinsnew,ipan_intervall,r2nefc_loop)
! !$omp& shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)
! !$omp& shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,cdentemp)
! !$omp& shared(tmatsph,den,denlm,gflle,gflle_part,rllleft,sllleft)
! !$omp& private(iq,df,ek,tmattemp,gmatll,gmat0,iorb,dentemp)
! !$omp& private(rho2ns_temp,dentot,dentmp,rho2,temp1)
! !$omp& shared(ldorhoef,nqdos,lmshift1,lmshift2,wez,lmsp,imt1,ifunm)
! !$omp& shared(r2orbc,r2nefc,cden,cdenlm,cdenns,rho2nsc_loop)
! !$omp& reduction(+:rho2int,espv) reduction(-:muorb) 
! !$omp& reduction(-:denorbmom,denorbmomsp,denorbmomlm,denorbmomns)
! !$omp& shared(t_tgmat)
! !$omp& private(alphasph,alphall)
      DO IE=1,IELAST
         ie_num = ie
         ie_end = ielast
#endif
         if(t_inc%i_write>0) WRITE(1337,*) 'energy',IE,EZ(IE)

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
     +              CLEB(1,1),ICLEB,IEND,JEND,EKL)
            ELSE
               CALL RHONS(DEN(0,IE,IQ),DF,DRDI,GMAT0,EK,
     +             RHO2NS(1,1,ISPIN),IPAN,IRCUT,IRMIN,THETAS,IFUNM,LMSP,  ! Added IRMIN 1.7.2014
     +              NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,
     +              JEND,IEND,EKL,DENLM(1,IE,IQ),GFLLE(:,:,IE,IQ))
            END IF
c LDA+U
            IF (IDOLDAU.EQ.1) THEN
             DO LM1=1,LMMAXD
              DO LM2=1,LMMAXD
               GLDAU(LM1,LM2)=GLDAU(LM1,LM2)+
     +                             DF*GFLLE(LM1,LM2,IE,IQ)
              ENDDO
             ENDDO
            ENDIF
           

#ifndef CPP_MPI
c Write out qdos:
            IF (OPT('qdos    ')) THEN                                     ! qdos
               DENTOT = DCMPLX(0.D0,0.D0)                                 ! qdos
               DO L = 0,LMAXD1                                            ! qdos
                  DENTOT = DENTOT + DEN(L,IE,IQ)                          ! qdos
               ENDDO                                                      ! qdos
               WRITE(30,9000) ERYD,QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),      ! lmdos
     &                       -DIMAG(DENTOT)/PI,                           ! lmdos
     &                      (-DIMAG(DENLM(L,IE,IQ))/PI,L=1,LMMAXD)        ! lmdos
               WRITE(31,9000) ERYD,QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),      ! qdos
     &                       -DIMAG(DENTOT)/PI,                           ! qdos
     &                      (-DIMAG(DENLM(L,IE,IQ))/PI,L=1,LMMAXD)        ! qdos
 9000          FORMAT(5F10.6,40E16.8)                                     ! qdos
               ! writeout complex qdos for interpolation                  ! complex qdos
               if(test('compqdos')) then                                  ! complex qdos
                 WRITE(3031,9001) ERYD,QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),  ! complex qdos
     &                           DENTOT,(DENLM(L,IE,IQ),L=1,LMMAXD)       ! complex qdos
               end if                                                     ! complex qdos
 9001          FORMAT(6F10.6,80E16.8)                                     ! qdos
            ENDIF                                                         ! qdos
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

            IQ = 1 ! reset IQ to zero, problem with qdos!!!

!              !!!! PROBLEM WITH ARRAY DIMENSIONS FOR VTREL ETC. !!!!
! #ifdef CPP_MPI
!              call MPI_FINALIZE(L)
! #endif
!              stop '[rhoval] ERROR array dimensions need to be checked!'
            CALL DRVRHO_QDOS(LDORHOEF,RHO2NS,R2NEF,DEN,DMUORB,RHOORB,
     &                  IE,ERYD,DF,LASTEZ,
     &                  GMATLL,VTREL,BTREL,RMREL,DRDIREL,
     &                  R2DRDIREL,ZREL,JWSREL,IRSHIFT,SOLVER,SOCTL,CTL,
     &                  QMTET,QMPHI,ITERMVDIR,MVEVIL,MVEVILEF,LMMAXD,
     &                  LMAXD,IRMD,LMPOTD,IEMXD,NMVECMAX,
     &                  I1,NQDOS)            ! qdos
     
C
            DO L = 0,LMAXD1
               ESPV(L,1) = ESPV(L,1) + DIMAG(ERYD*DEN(L,IE,IQ)*DF)
               ESPV(L,2) = ESPV(L,2) + DIMAG(ERYD*DEN(L,IE+IEMXD,IQ)*DF)
            END DO
C
            DO IR = 1,3
               DO L = 0,LMAXD
!                   write(1337,*) DF
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
C         IF ( ( IDOLDAU.EQ.1 ).AND.( LOPT.GE.0 ) ) 
C     &        CALL DENSITYMAT(DF,PZ,QZ,PNS,QNS,AR,CR,DR,GMATLL(1,1,IE),
C     &                        IPAN,IRCUT,DRDI,EK,
C     &                        IRMIN,LOPT,MMAX,LMLO,LMHI,PHILDAU,DENMATC
C     &        ,den,ie) ! test fivos 19.9.08
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C




      END DO                    ! IE = 1,IELAST

C LDA+U
      IF (IDOLDAU.EQ.1) THEN
        DO M1=1,MMAX
         LM1=LMLO-1+M1
         DENMATC(1:MMAX,M1)=(1.0/(2.0*CI))*
     +                       (GLDAU(LMLO:LMHI,LM1)-
     +                        DCONJG(GLDAU(LM1,LMLO:LMHI)))
        ENDDO
      ENDIF
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
      IF (IDOLDAU.EQ.1) DEALLOCATE(GLDAU)

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
      END
