      SUBROUTINE KKRMAT01(BZKP,NOFKS,GS,VOLCUB,VOLBZ,
     +                 TMATLL,MSSQ,
     +                 IE,IELAST,ITER,
     +                 ALAT,NSYMAT,NAEZ,CLS,NACLS,RR,EZOA,ATOM,
     +                 GINP,DGINP,
     >                 NUMN0,INDN0,IAT,
     >                 NATRC,ATTRC,EZTRC,NUTRC,INTRC,
     +                 SPRS,PRSC,EKM,NOITER,
     +                 EZ,QMRBOUND,IGUESS,BCP,CNVFAC,
     +                 DTDE_LOCAL,
     <                 GSXIJ,
     >                 NXIJ,XCCPL,IXCP,ZKRXIJ,
     <                 LLY_GRDT,TR_ALPH,
     >                 LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE,
     >                 LSMYRANK,LSRANK,LSMPIB,LSMPIC)
C
C
      IMPLICIT NONE
c ************************************************************************
c   performs k-space integration,
c   determines scattering path operator (g(k,e)-t**-1)**-1 and
c   Greens function of the real system -> GS(*,*,*,*),
c
c   NEW VERSION 10.99
c   up -> left , down -> right, for decimation 
c ------------------------------------------------------------------------
      include 'mpif.h'
c     .. parameters ..
      include 'inc.p'
      include 'inc.cls'
C
      INTEGER        LMAX,NSYMAXD
      PARAMETER      (LMAX=LMAXD,NSYMAXD=48)
      INTEGER        LMGF0D
      PARAMETER      (LMGF0D= (LMAXD+1)**2)
      INTEGER        LMMAXD
      PARAMETER      (LMMAXD= (LMAX+1)**2)
      INTEGER        ALM
      PARAMETER      (ALM = NAEZD*LMMAXD)
      INTEGER        ALMGF0
      PARAMETER      (ALMGF0 = NAEZD*LMGF0D)
      INTEGER        NGTBD
      PARAMETER      (NGTBD = NACLSD*LMMAXD)
      INTEGER        NBLCKD
      PARAMETER     (NBLCKD = XDIM*YDIM*ZDIM)
      INTEGER        LLYALM
      PARAMETER      (LLYALM=LLY*(NAEZD*LMMAXD-1)+1)
      DOUBLE COMPLEX CI
      PARAMETER      (CI=(0.D0,1.D0))
      DOUBLE COMPLEX CONE
      PARAMETER      (CONE  = ( 1.0D0,0.0D0))
      DOUBLE COMPLEX CIONE
      PARAMETER      (CIONE  = ( 0.0D0,-1.0D0))
      DOUBLE COMPLEX CZERO
      PARAMETER      (CZERO=(0.0D0,0.0D0))
c     ..
c     .. GLOBAL SCALER ARGUMENTS ..
      DOUBLE PRECISION ALAT,VOLBZ
      INTEGER          NAEZ,NOFKS,NSYMAT,IGUESS,BCP,
     +                 IE,IELAST,ITER,NXIJ,IXCP(NXIJD),EKM,
     +                 NOITER,IAT
c     ..
c     .. GLOBAL ARRAY ARGUMENTS ..
      DOUBLE COMPLEX
     +                 DGINP(LMGF0D,LMGF0D,NACLSD,NCLSD),
     +                 GINP(LMGF0D,LMGF0D,NACLSD,NCLSD),
     +                 GS(LMMAXD,LMMAXD,NSYMAXD),
     +                 GSXIJ(LMMAXD,LMMAXD,NSYMAXD,NXIJD),
     +                 EIKRP(NACLSD),EIKRM(NACLSD)
C     ..                 .. Lloyd
      DOUBLE COMPLEX   DTDE_LOCAL(LMMAXD,LMMAXD),
     +                 DGDE(LLYALM,LMMAXD),
     +                 GLLKE_X(LLYALM,LMMAXD),
     +                 MSSQ(LMMAXD,LMMAXD),TR_ALPH,
     +                 DPDE_LOCAL(LLYALM,LMMAXD),LLY_GRDT
C     ..                 .. precond
      DOUBLE COMPLEX   EZ(IEMXD)
      COMPLEX          PRSC(NGUESSD*LMMAXD,EKMD)
      DOUBLE PRECISION BZKP(3,KPOIBZ),VOLCUB(KPOIBZ),RR(3,0:NRD),
     +                 ZKRXIJ(48,3,NXIJD),
     +                 CNVFAC(EKMD)
      INTEGER          NUMN0(NAEZD),
     +                 INDN0(NAEZD,NACLSD),
     +                 ATOM(NACLSD,*),
     +                 CLS(*),
     +                 EZOA(NACLSD,*),
     +                 NACLS(*),
     +                 SPRS(NGUESSD*LMMAXD+1,EKMD+1)
      INTEGER          NUTRC,              ! number of inequivalent atoms in the cluster
     +                 INTRC(NATRCD),      ! pointer to atoms in the unit cell
     +                 NATRC,              ! number of atoms in cluster
     +                 ATTRC(NATRCD),      ! index to atom in elem/cell at site in cluster
     +                 EZTRC(NATRCD)       ! index to bravais lattice  at site in cluster
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE COMPLEX   TRACE,TRACEK,GTDPDE,BZTR2,CFCTORINV
      DOUBLE PRECISION TPI,QMRBOUND,DZNRM2
      INTEGER          IC,I1,I2,ILM,ISYM,IU,IL1,IL1B,IL2,IL2B,
     +                 KPT,LM,LM1,LM2,LM3,XIJ,ITCOUNT
      LOGICAL          XCCPL
C     ..
C     .. LOCAL ARRAYS ..
      DOUBLE COMPLEX   G(LMMAXD,LMMAXD),TGH(LMMAXD),
     +                 GLLKE1(ALM,LMMAXD),
     +                 DUMMY(ALM,LMMAXD),
     +                 GLLH(LMMAXD,NGTBD,NAEZD),
     +                 TMATLL(LMMAXD,LMMAXD,NAEZD),
     +                 GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*NBLCKD),
     +                 X0(ALM,LMMAXD)
      DOUBLE PRECISION N2B(LMMAXD)
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL CINIT,DLKE0,OUTTIME,DZNRM2
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC ATAN,EXP
C     ..
C     .. MPI ..
      INTEGER, dimension(MPI_STATUS_SIZE) :: STATUS
C     .. L-MPI
      INTEGER      MYLRANK(LMPID*SMPID*EMPID),
     +             LCOMM(LMPID*SMPID*EMPID),
     +             LGROUP(LMPID*SMPID*EMPID),
     +             LSIZE(LMPID*SMPID*EMPID),
     +             LMPI,LMPIC
C     .. LS-MPI
      INTEGER      LSMYRANK(LMPID,NAEZD*SMPID*EMPID),
     +             LSRANK(LMPID,NAEZD*SMPID*EMPID),
     +             LSMPI,LSMPIB,LSMPIC
C     .. N-MPI
      INTEGER MYRANK,NROFNODES,MAPBLOCK,IERR
      COMMON /MPI/MYRANK,NROFNODES
      SAVE
C
c     ..
c-----------------------------------------------------------------------
c
      TPI = 8.D0*ATAN(1.D0)    ! = 2*PI
      CFCTORINV = (CONE*TPI)/ALAT
c
      BZTR2 = CZERO
c
      DO 20 IU = 1,NSYMAXD
        CALL CINIT(LMMAXD**2,GS(1,1,IU))
 20   CONTINUE
c
      IF (XCCPL) THEN
        DO XIJ = 1, NXIJ
          DO ISYM = 1,NSYMAT
            CALL CINIT(LMMAXD*LMMAXD,GSXIJ(1,1,ISYM,XIJ))
          ENDDO             ! ISYM = 1,NSYMAT
        ENDDO
      ENDIF
c
c
c ---> use sit
c      G(n,n',L,L')(-k) = G(n',n,L',L)(k)
C
C
C=======================================================================
      DO 300 KPT = 1, NOFKS                               ! K-POINT-LOOP
C=======================================================================
c
c ---> fourier transformation
c
c     added by h.hoehler 3.7.2002
c
c                                                     n   0          n
c     define fourier transform as g mu mu'= ( sum_n g mu mu' exp(-iKR )
c                                   L  L'             L   L'
c
c                                             n   0           n
c                                 +   sum_n g mu'mu exp(-iK(-R ))) *0.5
c                                             L'  L
c
c     this operation has to be done to satisfy e.g. the point symmetry! 
c     application of fourier transformation is just an approximation
c     for the tb system, since the transl. invariance is not satisfied.
c

      IF (LLY.EQ.1) THEN 
      CALL CINIT(ALM*NGTBD,GLLH)

      DO I1 = 1,NAEZ

        IC = CLS(I1)
        CALL DLKE1(ALAT,NACLS,RR,EZOA(1,I1),
     +             BZKP(1,KPT),IC,EIKRM,EIKRP)
        CALL DLKE0(I1,GLLH,EIKRP,EIKRM,
     +             IC,NACLS,ATOM(1,I1),NUMN0,INDN0,DGINP(1,1,1,IC))
      END DO

        DO I1=1,NAEZ
          DO I2=1,NUMN0(I1)
            DO LM2=1,LMMAXD
              IL2=LMMAXD*(I2-1)+LM2
              IF (INDN0(I1,I2).EQ.IAT) THEN
                DO LM1=1,LMMAXD
                  IL1=LMMAXD*(I1-1)+LM1
                  DGDE(IL1,LM2)= GLLH(LM1,IL2,I1)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
        CALL CINIT(ALM*NGTBD,GLLH)

        DO I1 = 1,NAEZ
          IC = CLS(I1)
          CALL DLKE1(ALAT,NACLS,RR,EZOA(1,I1),
     +               BZKP(1,KPT),IC,EIKRM,EIKRP)
          CALL DLKE0(I1,GLLH,EIKRP,EIKRM,
     +               IC,NACLS,ATOM(1,I1),NUMN0,INDN0,
     +               GINP(1,1,1,IC))
        END DO

      IF (LLY.EQ.1) THEN 
        DO I1=1,NAEZ
          DO I2=1,NUMN0(I1)
            DO LM2=1,LMMAXD
              IL2=LMMAXD*(I2-1)+LM2
              IF (INDN0(I1,I2).EQ.IAT) THEN
                DO LM1=1,LMMAXD
                  IL1=LMMAXD*(I1-1)+LM1
                  GLLKE_X(IL1,LM2)= GLLH(LM1,IL2,I1)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      END IF 

C=======================================================================

        CALL CINIT(ALM*NGTBD,GLLH)

        DO I1 = 1,NAEZ
          IC = CLS(I1)
          CALL DLKE1(ALAT,NACLS,RR,EZOA(1,I1),
     +               BZKP(1,KPT),IC,EIKRM,EIKRP)
          CALL DLKE0(I1,GLLH,EIKRM,EIKRP,
     +               IC,NACLS,ATOM(1,I1),NUMN0,INDN0,
     +               GINP(1,1,1,IC))
        END DO

        DO I1=1,NAEZ
           IL1B=LMMAXD*(I1-1)
          DO I2=1,NUMN0(I1)
             DO LM2=1,LMMAXD
                IL2=LMMAXD*(I2-1)+LM2
                IL2B=LMMAXD*(INDN0(I1,I2)-1)+LM2
                DO LM1=1,LMMAXD
                   TGH(LM1) = CZERO
                   DO LM3=1,LMMAXD
                   TGH(LM1)=TGH(LM1)+TMATLL(LM1,LM3,I1)*GLLH(LM3,IL2,I1)
                   ENDDO
                ENDDO
                DO LM1=1,LMMAXD
                   IL1=IL1B+LM1
                   GLLH(LM1,IL2,I1) = TGH(LM1)
                   IF (IL1.EQ.IL2B) THEN
                     GLLH(LM1,IL2,I1) = GLLH(LM1,IL2,I1) - CONE
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
        ENDDO


C ####################################################################
C Start Lloyd's formula here
C
C enter this part exclusively for the last of the three energy points
C ####################################################################
C
C dP(E,k)   dG(E,k)                   dT(E)
C ------- = ------- * T(E) + G(E,k) * -----
C   dE        dE                       dE

      IF (LLY.EQ.1) THEN 

       CALL CINIT(LLYALM*LMMAXD,DPDE_LOCAL)

       CALL ZGEMM('N','N',ALM,LMMAXD,LMMAXD,CONE,
     +             DGDE,ALM,
     +             TMATLL(1,1,IAT),LMMAXD,CZERO,
     +             DPDE_LOCAL,ALM)
       CALL ZGEMM('N','N',ALM,LMMAXD,LMMAXD,CFCTORINV,
     +             GLLKE_X,ALM,
     +             DTDE_LOCAL,LMMAXD,CONE,DPDE_LOCAL,ALM)

      ENDIF
C#######################################################################
C LLOYD
C#######################################################################


C===================================================================
C===================================================================
C
C solve linear matrix equation in 5 subsequent steps
C
C===================================================================
C
C 1) find true residual tolerance by calculation of |b| 
C    using GLLKE1 as DUMMY-array ..
C
      CALL CINIT(ALM*LMMAXD,DUMMY)
C
      DO LM1=1,LMMAXD
        IL1=LMMAXD*(IAT-1)+LM1
        DO LM2=1,LMMAXD
          DUMMY(IL1,LM2)=-TMATLL(LM1,LM2,IAT)
        ENDDO
      ENDDO
C
      DO LM2=1,LMMAXD
        N2B(LM2) = DZNRM2(NAEZD*LMMAXD,DUMMY(1,LM2),1)
      ENDDO
C
C ..
C===================================================================
C
C===================================================================
C 2) if IGUESS is activated perform intitial step 'I' of 
C    intial guess - store b in DUMMY and set up modified b' ..
C
      IF (IGUESS.EQ.1) THEN
C
        CALL CINIT(ALM*LMMAXD,X0)
        CALL CINIT(ALM*LMMAXD,DUMMY)
C
        DO I1=1,NAEZD
        DO LM1=1,LMMAXD
          IL1=LMMAXD*(I1-1)+LM1
          DO LM2=1,LMMAXD
            DUMMY(IL1,LM2)=TMATLL(LM1,LM2,I1)
          ENDDO
        ENDDO
        ENDDO
C
        IF (ITER.GT.1) THEN
          CALL PRINVIT(
     >                 'I',IAT,
     >                 NUMN0,INDN0,
     >                 TMATLL,GLLH,X0,
     >                 PRSC(1,EKM+KPT),SPRS(1,EKM+KPT),
     <                 GLLKE1)
        ENDIF
      ENDIF

C ..
C===================================================================
C
C===================================================================
C 3) if BCP is activated determine preconditioning matrix
C    GLLHBLCK ..
C
      CALL CINIT(LMMAXD*NATBLD*LMMAXD*NATBLD*NBLCKD,GLLHBLCK)
C
      IF (BCP.EQ.1)
     > CALL BCPWUPPER(GLLH,GLLHBLCK,NAEZ,NUMN0,INDN0)
C ..
C===================================================================

C===================================================================
C 4) solve linear set of equation by iteration using BLAS-3 ..
C
          CALL CINIT(ALM*LMMAXD,GLLKE1)
C
          CALL MMINVMOD(
     >                  GLLH,GLLKE1,TMATLL,NUMN0,INDN0,N2B,
     >                  IAT,ITER,ITCOUNT,
     >                  GLLHBLCK,BCP,IGUESS,CNVFAC(EKM+KPT),
     >                  QMRBOUND)
C
          NOITER = NOITER + ITCOUNT
C
C ..
C===================================================================

C===================================================================
C 5) if IGUESS is activated perform final step 'F' of
C    intial guess           ~
C                  X = X  + X
C                       0          ..
C
      IF (IGUESS.EQ.1) THEN
        CALL PRINVIT(
     >               'F',IAT,
     >               NUMN0,INDN0,
     >               TMATLL,GLLH,X0,
     >               PRSC(1,EKM+KPT),SPRS(1,EKM+KPT),
     <               GLLKE1)
C
        DO I1=1,NAEZD
        DO LM1=1,LMMAXD
          IL1=LMMAXD*(I1-1)+LM1
          DO LM2=1,LMMAXD
            TMATLL(LM1,LM2,I1)=DUMMY(IL1,LM2)
          ENDDO
        ENDDO
        ENDDO
C
      ENDIF
C ..
C===================================================================
C
C solved. Result in GLLKE1
C
C===================================================================
C===================================================================



C#######################################################################
C LLOYD calculate Trace of matrix ...
C#######################################################################
        IF (LLY.EQ.1) THEN
C===================================================================
C                /  -1    dM  \
C calculate  Tr  | M   * ---- | 
C                \        dE  /
C===================================================================
C

      TRACEK=CZERO

      DO LM1=1,LMMAXD
        DO LM2=1,LMMAXD
          GTDPDE = CZERO
          DO IL1 = 1,LLYALM
            GTDPDE = GTDPDE + GLLKE1(IL1,LM2)*DPDE_LOCAL(IL1,LM1)
          ENDDO
          TRACEK = TRACEK + MSSQ(LM1,LM2)*GTDPDE
        ENDDO
      ENDDO
C
          BZTR2 = BZTR2 + TRACEK*VOLCUB(KPT)
C
        ENDIF
C#######################################################################
C LLOYD .
C#######################################################################
C
C===================================================================
C
C
C
C
C
          ILM = LMMAXD*(IAT-1) + 1

          DO 140 LM = 1,LMMAXD
             CALL ZCOPY(LMMAXD,GLLKE1(ILM,LM),1,G(1,LM),1)
 140      CONTINUE


          DO 110 ISYM = 1,NSYMAT
             DO LM1=1,LMMAXD
                DO LM2=1,LMMAXD
                   GS(LM1,LM2,ISYM) = GS(LM1,LM2,ISYM) +
     &                                   VOLCUB(KPT) * G(LM1,LM2)
                END DO
             END DO
 110      CONTINUE             ! ISYM = 1,NSYMAT
C
C
C
C ================================================================
        IF (XCCPL) THEN
C ================================================================
C       XCCPL communicate off-diagonal elements and multiply with
C       exp-factor
C ================================================================
C
          CALL KKRJIJ(
     >                BZKP,VOLCUB,KPT,
     >                NSYMAT,NAEZ,IAT,
     >                NXIJ,IXCP,ZKRXIJ,
     >                GLLKE1,
     <                GSXIJ,
     >                LMPIC,MYLRANK,
     >                LGROUP,LCOMM,LSIZE)
C
C ================================================================
        ENDIF
C ================================================================
C
C
C
C===================================================================
C
C
C=======================================================================
 300  END DO                    ! KPT = 1,NOFKS
C=======================================================================
C
C
      IF(LLY.EQ.1)  THEN
      BZTR2 = BZTR2*NSYMAT/VOLBZ + TR_ALPH
      TRACE=CZERO
      CALL MPI_ALLREDUCE(BZTR2,TRACE,1,
     +                   MPI_DOUBLE_COMPLEX,MPI_SUM,
     +                   LCOMM(LMPIC),IERR)
      LLY_GRDT = TRACE
      ENDIF

      RETURN

 9100 FORMAT (2e24.5)

      END
