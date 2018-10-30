       SUBROUTINE TMAT_NEWSOLVER(INS,NSPIN,LMAX,R,Z,E,KSRA,
     +                        IRWS,IPAN,IRCUT,IRMIN,C,
     +                        CLEB,ICLEB,IEND,NPAN_LOG,NPAN_EQ,
     +                        NCHEB,R_LOG,
     +                        VINS,VM2Z,PNS_SO,LEFTPNS,
     +                        TMATLL,THETA,PHI,I1,
     +                        IMPLAYER,TSTAR,DRDI)
       IMPLICIT NONE

c-----------------------------------------------------------------
c calculate tmatrix with new solver for SOC, N. H. Long, Juelich, 05.2013
c-----------------------------------------------------------------
       include 'inc.p'
       INTEGER INS,NSPIN,IRMIN,IPAN,IRWS,KSRA,LMAX,IEND,IATOM,IMPLAYER
       INTEGER LMMAXD
       PARAMETER (LMMAXD= (LMAXD+1)**2)
       INTEGER LMMAXSO
       PARAMETER (LMMAXSO=2*LMMAXD)
       INTEGER LMPOTD
       PARAMETER (LMPOTD= (LPOTD+1)**2)
       INTEGER IRMIND
       PARAMETER (IRMIND= IRMD-IRNSD)
       DOUBLE COMPLEX CZERO,CONE
       PARAMETER (CZERO=(0d0,0d0),CONE=(1d0,0d0))
       INTEGER IRCUT(0:IPAN)
       DOUBLE PRECISION C,Z
       DOUBLE COMPLEX E
       DOUBLE PRECISION R(IRMD),CLEB(*)
       INTEGER ICLEB(NCLEB,4)
       DOUBLE PRECISION 
     +   VINS(IRMIND:IRMD,LMPOTD,NSPIND),
     +   VM2Z(IRMD,NSPIND)
       DOUBLE COMPLEX  
     +   PNS_SO(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2),
     +   TMATLL(NSPD*LMMAXD,NSPD*LMMAXD),
     +   TMATTEMP(NSPD*LMMAXD,NSPD*LMMAXD)
       INTEGER I1,IR,IR2,NSRA,USE_SRATRICK,NVEC,LM1,LM2,LM3
       INTEGER NPAN_LOG,NPAN_EQ,NPAN_INST,NCHEB,NPAN_TOT,IRMDNEW
       DOUBLE PRECISION R_LOG,THETA,PHI
       DOUBLE COMPLEX GMATPREFACTOR
       DOUBLE PRECISION, ALLOCATABLE :: RNEW(:),RPAN_INTERVALL(:)
       INTEGER, ALLOCATABLE :: IPAN_INTERVALL(:)
       DOUBLE PRECISION, ALLOCATABLE :: VINSNEW(:,:,:)
       DOUBLE COMPLEX,ALLOCATABLE :: VNSPLL0(:,:,:),VNSPLL(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: HLK(:,:),JLK(:,:),
     +                                HLK2(:,:),JLK2(:,:),
     +                                RLL(:,:,:),SLL(:,:,:),
     +                                RLLLEFT(:,:,:),SLLLEFT(:,:,:),
     +                                RLLTEMP(:,:),PNSTEMP(:,:)
       DOUBLE COMPLEX TMATSPH(2*(LMAX+1))
       DOUBLE COMPLEX TSTAR(NSPD*LMMAXD,NSPD*LMMAXD)
c       DOUBLE COMPLEX, ALLOCATABLE ::  LEFTPNS(:,:,:),SUMR(:,:,:)
c      LEFTPNS is now allocated in main routine
       DOUBLE COMPLEX LEFTPNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2)
       DOUBLE COMPLEX, ALLOCATABLE ::  SUMR(:,:,:)
       DOUBLE PRECISION DRDI(IRMD)
       INTEGER JLK_INDEX(2*LMMAXSO)
       LOGICAL OPT
       EXTERNAL OPT

       IF (KSRA.GE.1) THEN
        NSRA=2
       ELSE
        NSRA=1
       ENDIF
c create new mesh
c data for the new mesh
       NPAN_INST= IPAN-1
       NPAN_TOT= NPAN_LOG+NPAN_EQ+NPAN_INST
       IRMDNEW= NPAN_TOT*(NCHEB+1)
c new mesh
       ALLOCATE(RNEW(IRMDNEW))
       ALLOCATE(RPAN_INTERVALL(0:NPAN_TOT))
       ALLOCATE(IPAN_INTERVALL(0:NPAN_TOT))
       ALLOCATE(VINSNEW(IRMDNEW,LMPOTD,NSPIND))
       CALL CREATE_NEWMESH(INS,NSPIN,R,LMMAXD,IRMIND,IRWS,IRMD,
     +                    LMPOTD,IPAN,IRCUT,IRMIN,VINS,VM2Z,R_LOG,
     +                    NPAN_LOG,NPAN_EQ,NPAN_INST,NCHEB,NPAN_TOT,
     +                    RNEW,RPAN_INTERVALL,IPAN_INTERVALL,VINSNEW)
c set up the non-spherical ll' matrix for potential VLL'
       IF (NSRA.EQ.2) THEN
        USE_SRATRICK=1
       ELSEIF (NSRA.EQ.1) THEN
        USE_SRATRICK=0
       ENDIF
       ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW))
       VNSPLL0=CZERO
       CALL VLLMAT(1,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINSNEW,
     +                 CLEB,ICLEB,IEND,NSPIN,Z,RNEW,USE_SRATRICK,I1)

c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW,RNEW,
     +                    E,Z,C,NSRA,NSPIN,LMPOTD,
     +                    IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,
     +                    NPAN_TOT*(NCHEB+1),VNSPLL0,THETA,PHI,'1')
c extend matrix for the SRA treatment
       IF (NSRA.EQ.2) THEN
        ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW))
        VNSPLL=CZERO
        IF (USE_SRATRICK.EQ.0) THEN
         CALL VLLMATSRA(VNSPLL0,VNSPLL,RNEW,LMMAXSO,IRMDNEW,E,C,
     +                 LMAX,0,'Ref=0')
        ELSEIF (USE_SRATRICK.EQ.1) THEN
         CALL VLLMATSRA(VNSPLL0,VNSPLL,RNEW,LMMAXSO,IRMDNEW,E,C,
     +                 LMAX,0,'Ref=Vsph')
        ENDIF
       ELSE
        ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW))
        VNSPLL=CZERO
        VNSPLL(:,:,:)=VNSPLL0(:,:,:)
       ENDIF

c calculate the source terms in the Lippmann-Schwinger equation
c these are spherical hankel and bessel functions
       ALLOCATE(HLK(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(JLK(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(HLK2(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(JLK2(1:4*(LMAX+1),IRMDNEW))
       HLK=CZERO
       JLK=CZERO
       HLK2=CZERO
       JLK2=CZERO
       GMATPREFACTOR=CZERO
       CALL RLLSLLSOURCETERMS(NSRA,NVEC,E,RNEW,IRMDNEW,LMAX,LMMAXSO,
     +                        1,JLK_INDEX,HLK,JLK,HLK2,JLK2,
     +                        GMATPREFACTOR)
c using spherical potential as reference
       IF (USE_SRATRICK.EQ.1) THEN
       CALL CALCSPH(NSRA,IRMDNEW,LMAX,NSPIN,Z,C,E,LMPOTD,LMMAXSO,
     +              RNEW,VINSNEW,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +              JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,TMATSPH,
     +              USE_SRATRICK)
       ENDIF

c calculate the tmat and wavefunctions
       ALLOCATE(RLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW))
       ALLOCATE(SLL(NVEC*LMMAXSO,LMMAXSO,IRMDNEW))
       RLL=(0d0,0d0)
       SLL=(0d0,0d0)

c right solutions
       TMATLL=CZERO
       CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL,RLL,SLL,TMATLL,NCHEB,
     +             NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +             NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,
     +             '1','1','0',USE_SRATRICK)
       IF (NSRA.EQ.2) THEN
        RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:)=
     +            RLL(LMMAXSO+1:NVEC*LMMAXSO,:,:)/C
        SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:)=
     +            SLL(LMMAXSO+1:NVEC*LMMAXSO,:,:)/C
        
       ENDIF
c add spherical contribution of tmatrix
       
       IF (USE_SRATRICK.EQ.1) THEN
        DO LM1=1,NSPD*LMMAXD
         TMATLL(LM1,LM1)=TMATLL(LM1,LM1)+TMATSPH(JLK_INDEX(LM1))
        ENDDO
       ENDIF

       DEALLOCATE(VNSPLL0)
       DEALLOCATE(VNSPLL)
       DEALLOCATE(HLK)
       DEALLOCATE(JLK)
       DEALLOCATE(HLK2)
       DEALLOCATE(JLK2)
       DEALLOCATE(SLL)
            
c transform radial wavefunction back to old mesh
       ALLOCATE(RLLTEMP(IRMDNEW,LMMAXSO))
       ALLOCATE(PNSTEMP(IRWS,LMMAXSO))
       DO LM1=1,LMMAXSO
        RLLTEMP=CZERO
        PNSTEMP=CZERO
        DO IR=1,IRMDNEW
         DO LM2=1,LMMAXSO
          RLLTEMP(IR,LM2)=RLL(LM1,LM2,IR)
         ENDDO
        ENDDO
        CALL CHEB2OLDGRID(IRWS,IRMDNEW,LMMAXSO,R,NCHEB,NPAN_TOT,
     +                    RPAN_INTERVALL,IPAN_INTERVALL,
     +                    RLLTEMP,PNSTEMP,IRMD)
        DO IR=1,IRWS
         DO LM2=1,LMMAXSO
          PNS_SO(LM1,LM2,IR,1)=PNSTEMP(IR,LM2)
         ENDDO
        ENDDO
       ENDDO ! LM1
c for small component
        IF (NSRA.EQ.2) THEN
         DO LM1=1,LMMAXSO
          RLLTEMP=CZERO
          PNSTEMP=CZERO
          DO IR=1,IRMDNEW
           DO LM2=1,LMMAXSO
            RLLTEMP(IR,LM2)=RLL(LM1+LMMAXSO,LM2,IR)
           ENDDO
          ENDDO
          CALL CHEB2OLDGRID(IRWS,IRMDNEW,LMMAXSO,R,NCHEB,NPAN_TOT,
     +                      RPAN_INTERVALL,IPAN_INTERVALL,
     +                      RLLTEMP,PNSTEMP,IRMD)
          DO IR=1,IRWS
           DO LM2=1,LMMAXSO
            PNS_SO(LM1,LM2,IR,2)=PNSTEMP(IR,LM2)
           ENDDO
          ENDDO
         ENDDO ! LM1
        ENDIF ! NSRA.EQ.2
       DEALLOCATE(RLL)
       DEALLOCATE(RLLTEMP)
       DEALLOCATE(PNSTEMP)
       
c rotate tmatrix and radial wavefunction to global frame
      CALL ROTATEMATRIX(TMATLL,THETA,PHI,LMMAXD,0)
      DO IR=1,IRMD
       CALL ROTATEMATRIX(PNS_SO(1,1,IR,1),THETA,PHI,LMMAXD,0)
       CALL ROTATEMATRIX(PNS_SO(1,1,IR,2),THETA,PHI,LMMAXD,0)
      ENDDO

c calcualte perturb t* using for calculating DP mechanism SOFIELD or
c write out left solutions LEFTSOL
      IF (OPT('SOFIELD ').or.OPT('LEFTSOL ')) THEN
        IF (OPT('SOFIELD').and.I1.NE.IMPLAYER) GO TO 200
c      now allocated in the main routine
c      ALLOCATE(LEFTPNS(NSPD*LMMAXD,NSPD*LMMAXD,IRMD))

c calculate left solution 
c set up the non-spherical ll' matrix for potential VLL'
       ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW))
       VNSPLL0=CZERO
       CALL VLLMAT(1,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINSNEW,
     +                 CLEB,ICLEB,IEND,NSPIN,Z,RNEW,USE_SRATRICK,I1)

c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW,RNEW,
     +                 E,Z,C,NSRA,NSPIN,LMPOTD,
     +                 IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,
     +                 NPAN_TOT*(NCHEB+1),VNSPLL0,THETA,PHI,'transpose')

c extend matrix for the SRA treatment
       IF (NSRA.EQ.2) THEN
        ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW))
        VNSPLL=CZERO
        IF (USE_SRATRICK.EQ.0) THEN
         CALL VLLMATSRA(VNSPLL0,VNSPLL,RNEW,LMMAXSO,IRMDNEW,E,C,
     +                 LMAX,0,'Ref=0')
        ELSEIF (USE_SRATRICK.EQ.1) THEN
         CALL VLLMATSRA(VNSPLL0,VNSPLL,RNEW,LMMAXSO,IRMDNEW,E,C,
     +                 LMAX,0,'Ref=Vsph')
        ENDIF
       ELSE
        ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW))
        VNSPLL=CZERO
        VNSPLL(:,:,:)=VNSPLL0(:,:,:)
       ENDIF

c calculate the source terms in the Lippmann-Schwinger equation
c these are spherical hankel and bessel functions
       ALLOCATE(HLK(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(JLK(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(HLK2(1:4*(LMAX+1),IRMDNEW))
       ALLOCATE(JLK2(1:4*(LMAX+1),IRMDNEW))
       HLK=CZERO
       JLK=CZERO
       HLK2=CZERO
       JLK2=CZERO
       GMATPREFACTOR=CZERO
       CALL RLLSLLSOURCETERMS(NSRA,NVEC,E,RNEW,IRMDNEW,LMAX,LMMAXSO,
     +                        1,JLK_INDEX,HLK,JLK,HLK2,JLK2,
     +                        GMATPREFACTOR)
c using spherical potential as reference
       IF (USE_SRATRICK.EQ.1) THEN
       CALL CALCSPH(NSRA,IRMDNEW,LMAX,NSPIN,Z,C,E,LMPOTD,LMMAXSO,
     +              RNEW,VINSNEW,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +              JLK_INDEX,HLK2,JLK2,HLK,JLK,GMATPREFACTOR,TMATSPH,
     +              USE_SRATRICK)
       ENDIF

c calculate the tmat and wavefunctions
       ALLOCATE(RLLLEFT(NVEC*LMMAXSO,LMMAXSO,IRMDNEW))
       ALLOCATE(SLLLEFT(NVEC*LMMAXSO,LMMAXSO,IRMDNEW))
       RLLLEFT=CZERO
       SLLLEFT=CZERO

c left solutions
       TMATTEMP=CZERO
       CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL,RLLLEFT,SLLLEFT,TMATTEMP,
     +           NCHEB,NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +             NSRA,JLK_INDEX,HLK2,JLK2,HLK,JLK,GMATPREFACTOR,
     +             '1','1','0',USE_SRATRICK)
       IF (NSRA.EQ.2) THEN
        RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:)=
     +            RLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:)/C
        SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:)=
     +            SLLLEFT(LMMAXSO+1:NVEC*LMMAXSO,:,:)/C
        
       ENDIF

       DEALLOCATE(VNSPLL0)
       DEALLOCATE(VNSPLL)
       DEALLOCATE(HLK)
       DEALLOCATE(JLK)
       DEALLOCATE(HLK2)
       DEALLOCATE(JLK2)
       DEALLOCATE(SLLLEFT)

c transform left solution back to old mesh
       ALLOCATE(RLLTEMP(IRMDNEW,LMMAXSO))
       ALLOCATE(PNSTEMP(IRWS,LMMAXSO))
       DO LM1=1,LMMAXSO
        RLLTEMP=CZERO
        PNSTEMP=CZERO
        DO IR=1,IRMDNEW
         DO LM2=1,LMMAXSO
          RLLTEMP(IR,LM2)=RLLLEFT(LM1,LM2,IR)
         ENDDO
        ENDDO
        CALL CHEB2OLDGRID(IRWS,IRMDNEW,LMMAXSO,R,NCHEB,NPAN_TOT,
     +                    RPAN_INTERVALL,IPAN_INTERVALL,
     +                    RLLTEMP,PNSTEMP,IRMD)
        DO IR=1,IRWS
         DO LM2=1,LMMAXSO
          LEFTPNS(LM1,LM2,IR,1)=PNSTEMP(IR,LM2)
         ENDDO
        ENDDO
       ENDDO ! LM1
c for small component
        IF (NSRA.EQ.2) THEN
         DO LM1=1,LMMAXSO
          RLLTEMP=CZERO
          PNSTEMP=CZERO
          DO IR=1,IRMDNEW
           DO LM2=1,LMMAXSO
            RLLTEMP(IR,LM2)=RLLLEFT(LM1+LMMAXSO,LM2,IR)
           ENDDO
          ENDDO
          CALL CHEB2OLDGRID(IRWS,IRMDNEW,LMMAXSO,R,NCHEB,NPAN_TOT,
     +                      RPAN_INTERVALL,IPAN_INTERVALL,
     +                      RLLTEMP,PNSTEMP,IRMD)
          DO IR=1,IRWS
           DO LM2=1,LMMAXSO
            LEFTPNS(LM1,LM2,IR,2)=PNSTEMP(IR,LM2)
           ENDDO
          ENDDO
         ENDDO ! LM1
        ENDIF ! NSRA.EQ.2
       DEALLOCATE(RLLLEFT)
       DEALLOCATE(RLLTEMP)
       DEALLOCATE(PNSTEMP)

      ENDIF!((OPT('SOFIELD ').and.I1.EQ.IMPLAYER).or.OPT('LEFTSOL '))
  200 CONTINUE

      IF (OPT('LEFTSOL ')) THEN
        DO IR=1,IRMD
       CALL ROTATEMATRIX(LEFTPNS(1,1,IR,1),THETA,PHI,LMMAXD,0)
       CALL ROTATEMATRIX(LEFTPNS(1,1,IR,2),THETA,PHI,LMMAXD,0)
      ENDDO 
      ENDIF

c calculate t* on I1=IMPLAYER site (defined in inputcard)
c t*=int{R_leftR_right}
      IF (OPT('SOFIELD ').and.I1.EQ.IMPLAYER) THEN

      ALLOCATE(SUMR(NSPD*LMMAXD,NSPD*LMMAXD,IRMD))
        TSTAR=CZERO
        SUMR=CZERO
        DO IR=1,IRWS
         DO LM1=1,NSPD*LMMAXD
          DO LM2=1,NSPD*LMMAXD
           DO LM3=1,NSPD*LMMAXD
            SUMR(LM1,LM2,IR)=SUMR(LM1,LM2,IR)+
     +        (LEFTPNS(LM3,LM1,IR,1))*PNS_SO(LM3,LM2,IR,1) ! left solution
c     +        PNS_SO(LM3,LM1,IR,1)*PNS_SO(LM3,LM2,IR,1) ! right solution
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        DO LM1=1,NSPD*LMMAXD
         DO LM2=1,NSPD*LMMAXD
c          CALL CSIMP3(SUMR(LM1,LM2,:),TSTAR(LM1,LM2),1,IRMIN,DRDI)
          CALL CSIMP3(SUMR(LM1,LM2,:),TSTAR(LM1,LM2),1,IRCUT(1),DRDI)
c          IF (ABS(TSTAR(LM1,LM2)).GT.1d-10) THEN
c           WRITE(55,'((2I5),(2e17.9))') LM1,LM2,TSTAR(LM1,LM2)
c         ENDIF
         ENDDO
        ENDDO       

c      DEALLOCATE(LEFTPNS)
      DEALLOCATE(SUMR)

      ENDIF!OPT('SOFIELD ').and.I1.EQ.IMPLAYER

      DEALLOCATE(RNEW)
      DEALLOCATE(VINSNEW)
      DEALLOCATE(RPAN_INTERVALL)
      DEALLOCATE(IPAN_INTERVALL)

       END 

