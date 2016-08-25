C ************************************************************************
       SUBROUTINE TMATRX01_SO(WFFLAG,ALPHA,ARSP,DRDI,E,ICST,INS,NSPO,
     &                        NSPOH,NSPIN,LMAX,PZ,QZ,FZ,
     +                        SZ,PNS_SO,TMATLL,R,VINS,VM2Z,Z,IRWS,IPAN,
     +                        IRCUT,IRMIN,IRMINSO,KSRA,C,DROR,RS,S,CLEB,
     +                        LOFLM,ICLEB,IEND,IE,IELAST,EF,REALFLAG,
     +                        PHASE_SHIFT,LSM,I1,IMPLAYER,TSTAR)
C ************************************************************************
c
c     this is the driver for the wavefunctions and t - matrices
c
c     in case of non spherical input potential the non spher.
c     wavefunctions are approximated inside a given sphere
c     with the nearly r - independent matrices :
c
c
c           the regular one (ir < irmin = irc - irns) :
c
c              pns(ir,lm1,lm2) = pz(ir,l1) * alfmat(lm1,lm2)
c
c          where pz is the regular wavefct of the spherically symmetric
c          part of the potential and alfmat the alpha matrix .
c
c
c           the irregular one (ir < irmin = irc - irns) :
c
c              qns(ir,lm1,lm2) = pz(ir,l1) * cr(lm1,lm2)
c
c                                    + qz(ir,l1) * dr(lm1,lm2)
c
c          where pz is the regular and qz is the irregular
c          wavefct of the spherically symmetric part of the
c          potential and cr , dr the matrices calculated
c          at the point irmin = irc - irns .
c
c     to save storage the non spherical wavefunctions are stored only
c        from the sphere boundary up to irc - irns
c
c             (see notes by b.drittler)
c
c             changed for band structure code
c
c                               b.drittler   nov. 1989
c                               valerio 14.6.99
C ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
c      INTEGER NCLEB
c      PARAMETER (NCLEB=LMPOTD*LMMAXD)
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER (CONE= (1.0D0,0.0D0),CZERO= (0d0,0d0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DET,E,EF
      DOUBLE PRECISION C,Z
      INTEGER ICST,IEND,INS,IPAN,IRMIN,IRWS,KSRA,LMAX,IELAST,IE,
     +        NSPO,IRMINSO,NSPIN,NSPOH,I1
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX    ALPHA(0:LMAXD,NSPIND),
     +                  FZ(IRMD,0:LMAXD,NSPIND),
     +                  PZ(IRMD,0:LMAXD,NSPIND),
     +                  QZ(IRMD,0:LMAXD,NSPIND),SZ(IRMD,0:LMAXD,NSPIND),
     +                  PNS_SO(NSPD*LMMAXD,NSPD*LMMAXD,IRMINSO:IRMD,2),
     +                  TMATLL(NSPD*LMMAXD,NSPD*LMMAXD),
     +                  ARSP(NSPD*LMMAXD,NSPD*LMMAXD),
     +                  PHASE_SHIFT(0:LMAXD),LSM(2*LMMAXD,2*LMMAXD)
      DOUBLE PRECISION  CLEB(*),DRDI(IRMD),DROR(IRMD),R(IRMD),
     +                  RS(IRMD,0:LMAXD),S(0:LMAXD),
     +                  VINS(IRMIND:IRMD,LMPOTD,NSPIND),
     +                  VM2Z(IRMD,NSPIND)
      INTEGER           ICLEB(NCLEB,4),IRCUT(0:IPAND),LOFLM(*)
      DOUBLE COMPLEX    TSTAR(NSPD*LMMAXD,NSPD*LMMAXD),
     +                  SUMR(NSPD*LMMAXD,NSPD*LMMAXD,IRMD)
C     ..
C     .. Local Scalars ..
      LOGICAL           TEST,NOSOC
      DOUBLE COMPLEX    EK
      INTEGER           I,INFO,IRC1,LM1,LMMAX,NSRA
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX,ALLOCATABLE ::   CR(:,:,:),
     +                                DR(:,:,:),
     +                                PNS(:,:,:,:),
     +                                QNS(:,:,:,:),
     +                                TMATLL1(:,:),
     +                                AR(:,:),
     +                                CMAT(:,:,:),
     +                                DMAT(:,:,:),
     +                                EFAC(:),PZEKDR(:,:,:,:),
     +                                PZLM(:,:,:,:),
     +                                QZEKDR(:,:,:,:),
     +                                QZLM(:,:,:,:),
     +                                TMAT(:,:)!,LSM(:,:)
      DOUBLE PRECISION,ALLOCATABLE :: VNSPLL(:,:,:,:),HSOFAC(:),RDVDR(:)
      DOUBLE PRECISION :: SQPZ(LMMAXD,IRMD),NORMPZ(LMMAXD)
      DOUBLE PRECISION :: SQPZ1(LMMAXD,IRMD),NORMPZ1(LMMAXD)
      DOUBLE COMPLEX :: NORMPZLM(LMMAXD,IRMD)
      INTEGER :: IRMT,SOS

      INTEGER IPVT(LMMAXD)
      LOGICAL WFFLAG,REALFLAG
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,CPLXWB01,IRWNS,REGNS,VLLNS,WFTSCA,ZGETRF,TEST
C     ..
      LMMAX = (LMAX+1)* (LMAX+1)

      ALLOCATE(CR(LMMAXD,LMMAXD,NSPIND))
      ALLOCATE(DR(LMMAXD,LMMAXD,NSPIND))
      ALLOCATE(PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2))
      ALLOCATE(QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2))
      ALLOCATE(TMATLL1(LMMAXD,LMMAXD))
      ALLOCATE(AR(LMMAXD,LMMAXD))
      ALLOCATE(CMAT(LMMAXD,LMMAXD,IRMIND:IRMD))
      ALLOCATE(DMAT(LMMAXD,LMMAXD,IRMIND:IRMD))
      ALLOCATE(EFAC(LMMAXD))
      ALLOCATE(PZEKDR(LMMAXD,IRMINSO:IRMD,2,NSPIND))
      ALLOCATE(PZLM(LMMAXD,IRMINSO:IRMD,2,NSPIND))
      ALLOCATE(QZEKDR(LMMAXD,IRMINSO:IRMD,2,NSPIND))
      ALLOCATE(QZLM(LMMAXD,IRMINSO:IRMD,2,NSPIND))
      ALLOCATE(TMAT(0:LMAXD,NSPIND))
c      ALLOCATE(LSM(LMMAXD*2,LMMAXD*2))
      IF (INS .EQ. 1 ) THEN
        ALLOCATE(VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD,NSPIND))
      ELSE
        IRMIN=1
        ALLOCATE(VNSPLL(LMMAXD,LMMAXD,IRMIN:IRMD,NSPIND))
      END IF
      ALLOCATE(HSOFAC(IRMINSO:IRMD))
      ALLOCATE(RDVDR(IRMD))

      IF (KSRA.GE.1) THEN    ! previously this was .GT. which is wrong for kvrel=1
         NSRA = 2
      ELSE
         NSRA = 1
      END IF

c     this was giving a value of 3 for KSRA=2 which was WRONG!
c      NSRA = KSRA + 1

      CALL CINIT(NSPD*LMMAXD*NSPD*LMMAXD,TMATLL)
c
c---> determine wavefunctions and t matrix for spherical averaged pot.
c
      DO ISP=1,NSPIN

        CALL CPLXWB01(PZ(:,:,ISP),FZ(:,:,ISP),QZ(:,:,ISP),SZ(:,:,ISP),
     +            TMAT(:,ISP),ALPHA(:,ISP),E,EK,KSRA,IPAN,
     &            IRCUT,IRWS,VM2Z(:,ISP),
     +            DRDI,R,Z,LMAX,C,DROR,RS,S,IE,IELAST,EF,
     +            REALFLAG,PHASE_SHIFT)

        DET = CONE

      END DO

      IF(TEST('flow    ')) 
     +    write(6,*) "after regular spherical wavefunction"

      IF (INS .EQ. 0 .AND. NSPOH .EQ.1 ) THEN
        
        DO ISP=1,NSPIN
          DO 10 LM1 = 1,LMMAX
            TMATLL(LM1+LMMAX*(ISP-1),LM1+LMMAX*(ISP-1)) = 
     +                                 TMAT(LOFLM(LM1),ISP)
   10     CONTINUE
        END DO


c---> non-spherical input potential, but no spin-orbit coupling

      ELSE IF (INS .EQ. 1 .AND. NSPO .EQ. 1) THEN
c
c---> non spherical input potential
c
        IRC1 = IRCUT(IPAN)
c
c---> determine the lm,lm' dependent potential
c
        CALL VLLNS(IRMIN,IRC1,LMMAX,VNSPLL,VINS,CLEB,ICLEB,IEND,NSPIN)
        

        IF(TEST('flow    ')) 
     +    write(6,*) "after VLLNS (lm,lm' dependent potential)"
c
c---> get wfts of same magnitude by scaling with efac

        DO ISP=1,NSPIN
  
          CALL WFTSCA(DRDI,EFAC,LMAX,PZ(:,:,ISP),QZ(:,:,ISP),IRMIN,IRWS,
     +         IPAN,IRCUT,FZ(:,:,ISP),SZ(:,:,ISP),
     +         NSRA,PZLM(:,IRMIND:IRMD,:,ISP),QZLM(:,IRMIND:IRMD,:,ISP),
     +         PZEKDR(:,IRMIND:IRMD,:,ISP),QZEKDR(:,IRMIND:IRMD,:,ISP),
     +         EK,LOFLM)


        END DO

c
c---> determine the regular non sph. wavefunction
c
        PNS(:,:,:,:)=0d0
        PNS_SO(:,:,:,:)=0d0
        TMATLL=0d0
        ARSP=0d0
        AR=0d0

        DO ISP=1,NSPIN

          PNS(:,:,:,:)=0d0
          AR=0d0
          TMATLL1=0d0

          WRITE(6,*) " before regns"

          CALL REGNS_VOLTERRA(AR(:,:),TMATLL1,EFAC,PNS,
     +             VNSPLL(:,:,:,ISP),ICST,IPAN,IRCUT,
     +             PZLM(:,IRMIND:IRMD,:,ISP),QZLM(:,IRMIND:IRMD,:,ISP),
     +             PZEKDR(:,IRMIND:IRMD,:,ISP),
     +             QZEKDR(:,IRMIND:IRMD,:,ISP),
     +             EK,CMAT,DMAT,NSRA,IRMIN,IRWS,IPAND,LMMAXD)
        
          PNS_SO((ISP-1)*LMMAXD+1:(ISP-1)*LMMAXD+LMMAXD,
     +           (ISP-1)*LMMAXD+1:(ISP-1)*LMMAXD+LMMAXD,IRMIND:IRMD,:)=
     +                     PNS(1:LMMAXD,1:LMMAXD,IRMIND:IRMD,:)
          ARSP((ISP-1)*LMMAXD+1:(ISP-1)*LMMAXD+LMMAXD,
     +           (ISP-1)*LMMAXD+1:(ISP-1)*LMMAXD+LMMAXD)=AR(:,:)

          TMATLL((ISP-1)*LMMAXD+1:(ISP-1)*LMMAXD+LMMAXD,
     +           (ISP-1)*LMMAXD+1:(ISP-1)*LMMAXD+LMMAXD)=TMATLL1(:,:)

        END DO

        DO 20 LM1 = 1,LMMAX
          TMATLL(LM1,LM1) = TMATLL(LM1,LM1) + TMAT(LOFLM(LM1),1)
   20   CONTINUE

        IF (NSPIN ==2 ) THEN
          DO LM1 = LMMAX+1,2*LMMAX
            TMATLL(LM1,LM1) = TMATLL(LM1,LM1) + 
     +                    TMAT(LOFLM(LM1-LMMAX),NSPIN)
          END DO
        END IF 

c---> Spin-orbit coupling

      ELSE IF (NSPO.EQ.2) THEN

c---> calculate the spin-orbit Hamiltonian 

        IF (Z .GT. 0 ) THEN
c ic         WRITE(6,*) "Z before HAM",z
          WRITE(6,*) "SOC"

          CALL HAM_SO(LMAX,LMMAXD,IRMD,IRCUT(IPAN),
     +           VM2Z,R,IE,E,Z,C,
     +           NSRA,LSM,HSOFAC,IRMINSO,IPAN,IRCUT(0:IPAN),NSPIN,RDVDR)

c        WRITE(6,*) "after HAM_SO"

        ELSE !set Spin-orbit coupling to zero for the vacuum

          HSOFAC=0d0
        
        END IF

c       IF (IE==1 ) THEN
!        OPEN (UNIT=34,FILE="HAM_SO",FORM="FORMATTED")

!        DO IR=IRMINSO,IRMD
!          WRITE(34,"((I5),(5e17.9))") IR,R(IR),HSOFAC(IR)
!c             HSOFAC(IR)=0.0d0*HSOFAC(IR)
!c            WRITE(35,"((I5),(5e17.9))") IR,R(IR),HSOFAC(IR)
!        END DO

!        CLOSE(34)
c       END IF

c        NOSOC=.TRUE.
        NOSOC=.FALSE.
        IF (NOSOC == .TRUE. .OR. NSPOH==1) THEN
c        IF (NOSOC == .TRUE.) THEN

          DO IR=IRMINSO,IRMD
            HSOFAC(IR)=0.0d0*HSOFAC(IR)
          END DO

          WRITE(6,*) "HSOFAC set to zero!!!!"

        END IF

        IF(TEST('flow    ')) 
     +    write(6,*) "after spin-orbit Hamiltonian" 

        IRC1 = IRCUT(IPAN)         
c
c---> non-spherical input potential and spin-orbit coupling
c
        IF (INS.EQ.1) THEN 
c
c---> determine the lm,lm' dependent potential

          CALL VLLNS(IRMIN,IRC1,LMMAX,VNSPLL,VINS,CLEB,ICLEB,IEND,NSPIN)
c          VNSPLL=0d0

        ELSE

          VNSPLL=0d0
      
        END IF

c---> get wfts of same magnitude by scaling with efac
c        write(6,*) "before scaling of the WF"

        DO ISP=1,NSPIN 
c        DO LM1=1,IRMD
c        DO LM2=0,LMAXD
c        write(41,*) LM1,LM2,FZ(LM1,LM2,ISP)
c        ENDDO
c        ENDDO

         write(6,*) "befor WFTSCA"
          CALL WFTSCA_SO(DRDI,EFAC,LMAX,PZ(:,:,ISP),QZ(:,:,ISP),IRMINSO,
     +                  IRWS,IPAN,IRCUT,FZ(:,:,ISP),SZ(:,:,ISP),NSRA,
     +                  PZLM(:,:,:,ISP),QZLM(:,:,:,ISP),
     +                  PZEKDR(:,:,:,ISP),QZEKDR(:,:,:,ISP),EK,LOFLM,
     +                  LMMAXD)

        END DO
c calculate the spin-orbit strength
        SOS=0
        IF (SOS.EQ.1) THEN
        IF (I1.EQ.5) THEN
c        IRMT=IRCUT(1)
        IRMT=IRMD
        SQPZ=0d0
        NORMPZ=0d0
        SQPZ1=0d0
        NORMPZ1=0d0
        DO LM1=0,LMAXD
c        DO LM1=1,LMMAXD
c         DO LM1=2,2
         DO IR=IRMINSO,IRMD
          SQPZ(LM1,IR)=PZ(IR,LM1,1)*DCONJG(PZ(IR,LM1,1))
         ENDDO
         CALL SIMP3(SQPZ(LM1,1),NORMPZ(LM1),IRMINSO,IRMT,DRDI)
         DO IR=IRMINSO,IRMD 
          NORMPZLM(LM1,IR)=PZ(IR,LM1,1)/SQRT(NORMPZ(LM1))
         ENDDO
         DO IR=IRMINSO,IRMD
          SQPZ1(LM1,IR)=NORMPZLM(LM1,IR)*DCONJG(NORMPZLM(LM1,IR))
     +                  *RDVDR(IR)
c     +                  *HSOFAC(IR)
         ENDDO
         CALL SIMP3(SQPZ1(LM1,1),NORMPZ1(LM1),IRMINSO,IRMT,DRDI)
        write(6,*) LM1,NORMPZ1(LM1)
       ENDDO
       STOP
       ENDIF
       ENDIF
        write(6,*) "after WFTSCA"

c       write(6,*) "after scaling of the WF"
        IF(TEST('flow    ')) 
     +       write(6,*) "after scaling of the WF"

c---> determine the regular non sph. wavefunction
c
        PNS_SO(:,:,:,:)=0d0

        WRITE(6,*) "before REGNS_SO"
c        HSOFAC=0d0

        CALL REGNS_SO(ARSP,TMATLL,EFAC,PNS_SO,VNSPLL,LSM,HSOFAC,ICST,
     +             IPAN,IRCUT,PZLM,QZLM,PZEKDR,QZEKDR,R,
     +             EK,NSRA,IRMIN,IRMINSO,IRWS,IPAND,LMMAXD,
     +             LMMAXD*NSPO,Z,INS,NSPIN)
        write(6,*) 'after regns spin orbit'

        IF(TEST('flow    ')) 
     +    write(6,*) 'after regns spin orbit'
c calculate t* for perturbed deltaq
        IF (I1.EQ.IMPLAYER) THEN
        IRMT=IRMIN
        SUMR=CZERO
        TSTAR=CZERO
        DO IR=IRMINSO,IRWS
         DO LM1=1,NSPD*LMMAXD
          DO LM2=1,NSPD*LMMAXD
           DO LM3=1,NSPD*LMMAXD
            SUMR(LM1,LM2,IR)=SUMR(LM1,LM2,IR)+
     +       PNS_SO(LM3,LM1,IR,1)*PNS_SO(LM3,LM2,IR,1)
c     +       +(0d0,1d0)*PNS_SO(LM3,LM1,IR,2)*PNS_SO(LM3,LM2,IR,2)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        DO LM1=1,NSPD*LMMAXD
         DO LM2=1,NSPD*LMMAXD
          CALL CSIMP3(SUMR(LM1,LM2,1),TSTAR(LM1,LM2),IRMINSO,IRMT,DRDI)
         ENDDO
        ENDDO
        ENDIF
!       OPEN(UNIT=55,file="TMAT",form="formatted")
!      DO LM1=1,LMMAX*NSPD
!        DO LM2=1,LMMAX*NSPD
!          WRITE(55,"((2I5),(4e17.9))") LM2,LM1,
!     +               TMATLL(LM2,LM1),TMATLL(LM1,LM2)
!        END DO
!      END DO
!      CLOSE(55)

!      OPEN(UNIT=55,file="TMAT_SPHER",form="formatted")
!      DO LM1=1,LMMAX
!        WRITE(55,"((I5),(4e17.9))") LM1,
!     +               TMAT(LOFLM(LM1),1)
!      END DO
!      CLOSE(55)
c       WRITE(6,*) "NSPIN",NSPIN
c       OPEN(UNIT=55,file="PNS_test_in_tmat",form="formatted")
c       WRITE(56,*) "TMATLL"
c       DO IR=1,IRMD
c         DO LM1=1,LMMAX*NSPD
c           DO LM2=1,LMMAX*NSPD
c            WRITE(56,"((2I5),(2e17.9))") LM2,LM1,TMATLL(LM2,LM1)
c             WRITE(55,"((3I5),(4e17.9))") IR,LM2,LM1,
c    +                   PNS_SO(LM2,LM1,IR,1),PNS_SO(LM1,LM2,IR,1)
c           END DO
c         END DO
c          WRITE(55,"((I5),(2e17.9))") LM1,TMAT(LOFLM(LM1),1)
c       END DO
c       CLOSE(55)
c        STOP " TMATLL"

c---> calculate full t-matrix, add the spherical part to the
c     non-spherical contribution calculated in regns_SO           
        DO LM1 = 1,LMMAX
          TMATLL(LM1,LM1) = TMATLL(LM1,LM1) + TMAT(LOFLM(LM1),1)
          TMATLL(LMMAX+LM1,LMMAX+LM1)=TMATLL(LMMAX+LM1,LMMAX+LM1) 
     +                                       + TMAT(LOFLM(LM1),NSPIN)
        END DO
        IF (I1.EQ.IMPLAYER) THEN
        DO LM1=1,NSPD*LMMAXD
         DO LM2=1,NSPD*LMMAXD
c         WRITE(55,'((2I5),(2e17.9))') LM1,LM2,TSTAR(LM1,LM2)
c         WRITE(56,'((2I5),(2e17.9))') LM1,LM2,TMATLL(LM1,LM2)
         ENDDO
        ENDDO
c        STOP
        ENDIF
      END IF         !end of spin-orbit coupling

c      OPEN(UNIT=55,file="PNS_test_in_tmat",form="formatted")
c     OPEN(UNIT=56,file="PNS_TRANS_S",form="formatted")
c     OPEN(UNIT=57,file="PNS_TRANS_P",form="formatted")
c     OPEN(UNIT=58,file="PNS_TRANS_D",form="formatted")
cc      DO IR=IRMIND,IRMD

c     DO IR=IRMINSO,IRMD

c       DO LM1=1,LMMAX*NSPD
c         DO LM2=1,LMMAX*NSPD
c           WRITE(55,"((3I5),(4e17.9))") IR,LM2,LM1,
c    +                 PNS_SO(LM2,LM1,IR,1),PNS_SO(LM1,LM2,IR,1)
c         END DO
c       END DO

c       DO LM1=1,2        
c         DO LM2=1,2
c           WRITE(56,"((3I5),(4e17.9))") IR,LM2,LM1,
c    +                 PNS_SO(LM2,LM1,IR,1),PNS_SO(LM1,LM2,IR,1)
c         END DO
c       END DO

c       DO LM1=3,8        
c         DO LM2=3,8
c           WRITE(57,"((3I5),(4e17.9))") IR,LM2,LM1,
c    +                 PNS_SO(LM2,LM1,IR,1),PNS_SO(LM1,LM2,IR,1)
c         END DO
c       END DO

c       DO LM1=9,18       
c         DO LM2=9,18
c           WRITE(58,"((3I5),(4e17.9))") IR,LM2,LM1,
c    +                 PNS_SO(LM2,LM1,IR,1),PNS_SO(LM1,LM2,IR,1)
c         END DO
c       END DO

c      END DO

c     CLOSE(55)
c     CLOSE(56)
c     CLOSE(57)
c     CLOSE(58)

c      OPEN(UNIT=55,file="TMAT",form="formatted")
c       OPEN(UNIT=56,file="TMAT_ABS",form="formatted")
c      DO LM1=1,LMMAX*NSPD
c       DO LM2=1,LMMAX*NSPD
c         WRITE(55,"((2I5),(4e17.9))") LM2,LM1,
c     +               TMATLL(LM2,LM1),TMATLL(LM1,LM2)
c           WRITE(56,"((2I5),(5e17.9))") LM2,LM1,
c      +               TMATLL(LM2,LM1)-TMATLL(LM1,LM2),
c      +               ABS(REAL(TMATLL(LM2,LM1)-TMATLL(LM1,LM2))),
c      +               ABS(AIMAG(TMATLL(LM2,LM1)-TMATLL(LM1,LM2))),
c      +               ABS(TMATLL(LM2,LM1)-TMATLL(LM1,LM2))
c       END DO
c      END DO
      CLOSE(55)
c       CLOSE(56)

c      STOP "AFTER TMAT " 
      DEALLOCATE(CR)
      DEALLOCATE(DR)
      DEALLOCATE(PNS)
      DEALLOCATE(QNS)
      DEALLOCATE(TMATLL1)
      DEALLOCATE(AR)
      DEALLOCATE(CMAT)
      DeALLOCATE(DMAT)
      DEALLOCATE(EFAC)
      DEALLOCATE(PZEKDR)
      DEALLOCATE(PZLM)
      DEALLOCATE(QZEKDR)
      DeALLOCATE(QZLM)
      DEALLOCATE(TMAT)
c      DEALLOCATE(LSM)
      DEALLOCATE(VNSPLL)
      DEALLOCATE(HSOFAC)
      DEALLOCATE(RDVDR)

      END













