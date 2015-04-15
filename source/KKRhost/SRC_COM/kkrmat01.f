      SUBROUTINE KKRMAT01(BZKP,NOFKS,GS,VOLCUB,TINVLL,RROT,
     &                    NSHELL,NSDIA,ALAT,NSYMAT,
     &                    NAEZ,CLS,NACLS,NACLSMAX,RR,EZOA,ATOM,
     &                    NSH1,NSH2,GINP,RBASIS,RCLS,
     &                    TINVBUP,TINVBDOWN,VACFLAG,NLBASIS,NRBASIS,
     &                    FACTL,ICHECK,INVMOD,IDECI,SRREL,IRREL,NRREL,
     &          DTREFLL,DTMATLL,DGINP,REFPOT,LLY_GRTR,TRACET,CFCTOR,LLY)   ! LLY
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
      IMPLICIT NONE
!     ..
!     .. Parameters ..
      include 'inc.p'
!OMPI include 'mpif.h'
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
     +        NS,COUNTER1,COUNTER2,IL2,IL1T,IL2T,IS1,IS2,
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
!----------------------------------------------------------------
!F77!F77-------------------------------------------------------------
!F77      DOUBLE COMPLEX GLLKE(ALM,ALM),GLLKEM(ALM,ALM)
!F77      DOUBLE COMPLEX GLLKE0(KREL*ALMGF0+(1-KREL),
!F77     &                      KREL*ALMGF0+(1-KREL))
!F77      DOUBLE COMPLEX GLLKE0M(KREL*ALMGF0+(1-KREL),
!F77     &                       KREL*ALMGF0+(1-KREL))
!F77!F77-------------------------------------------------------------
!---!----------------------------------------------------------------
!---      DOUBLE COMPLEX GLLKE(:,:),GLLKEM(:,:)
!---      DOUBLE COMPLEX GLLKE0(:,:),GLLKE0M(:,:)
!---      ALLOCATABLE GLLKE,GLLKEM,GLLKE0,GLLKE0M
!---!----------------------------------------------------------------
 
      DOUBLE COMPLEX ETAIKR(NSYMAXD,NSHELD),G(LMMAXD,LMMAXD)
      DOUBLE PRECISION BZKPK(6),KP(3)
      INTEGER NDIM


      LOGICAL TEST,OPT
!     ..
!     .. External subroutines ..
      EXTERNAL CINIT,DLKE0,OPT,TEST,GTDYSON
!     ..
!     .. Intrinsic functions ..
      INTRINSIC ATAN,EXP
!     ..
!OMPI INTEGER MYRANK,NROFNODES
!OMPI COMMON /MPI/MYRANK,NROFNODES
!OMPI DOUBLE COMPLEX WORK(LMMAXD,LMMAXD,NSYMAXD)
!OMPI INTEGER IERR,IWORK,MAPBLOCK
!     ..

!      NDIM=LMGF0D*NAEZ
      NDIM=LMMAXD*NAEZ

      IF ( TEST('flow     ') )
     +     WRITE(6,*) '>>> kkrmat1: loop over k-points'
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

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ K-points loop
      DO KPT = 1,NOFKS 
         GLLKE(:,:) = CZERO
         IF (LLY.NE.0) DGLLKE(:,:) = CZERO
!OMPI   IF(MYRANK.EQ.MAPBLOCK(K,1,NOFKS,1,0,NROFNODES-1)) THEN
         KP(1:3) = BZKP(1:3,KPT)
 
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
      
            !     --> Perform the inversion of matrix M
            !     the output is the scattering path operator TAU stored in GLLKE
            !     Actually -TAU, because TAU = (Deltat^-1 - Gref)^-1
 

            CALL INVERSION(GLLKE,INVMOD,ICHECK)

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
            END DO
            ! ----------------------------------------------------------------------
         END DO
         ! ======================================================================
         IF (LLY.NE.0) LLY_GRTR =                                   ! LLY Lloyd Integration
     &                 LLY_GRTR + LLY_GRTR_K * VOLCUB(KPT) * NSYMAT ! LLY Lloyd Integration



!OMPI END IF

      END DO ! KPT = 1,NOFKS   end K-points loop

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
!OMPI   DO NS = 1,NSHELL
!OMPI     IWORK = LMMAXD*LMMAXD*NSYMAXD
!OMPI     CALL MPI_ALLREDUCE(GS(1,1,1,NS),WORK,IWORK,
!OMPI+                       MPI_DOUBLE_COMPLEX,MPI_SUM,
!OMPI+                       MPI_COMM_WORLD,IERR)
!OMPI     CALL ZCOPY(IWORK,WORK,1,GS(1,1,1,NS),1)
!OMPI   END DO

      IF ( TEST('flow    ') ) WRITE(6,*) '<<< KKRMAT1'


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
