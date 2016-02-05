       SUBROUTINE TMATIMP_NEWSOLVER(INS,NSPIN,LMAX,R,Z,IELAST,E,KSRA,
     +                     IRWS,IPAN,IRCUT,IRMIN,C,CLEB,ICLEB,IEND,
     +                     NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,VINS,VM2Z,
     +                     NATOMIMP,RCLSIMP,ATOMIMP,IHOST,HOSTIMP,RIMP,
     +                     ZIMP,IRWSIMP,IPANIMP,IRCUTIMP,IRMINIMP,
     +                     VINSIMP,VM2ZIMP,DTMTRX)
       IMPLICIT NONE

c-----------------------------------------------------------------
c calculate and write down impurity tmatrix and delta matrix
c N. H. Long, Juelich, 05.2013
c-----------------------------------------------------------------
       include 'inc.p'
       INTEGER INS,NSPIN,IELAST,IRMIN(*),IPAN(*),IRWS(*),KSRA,LMAX,IEND
       INTEGER IRMINIMP(*),IPANIMP(*),IRWSIMP(*),ATOMIMP(NATOMIMPD)
       INTEGER NATOMIMP,NPAN_TOTD,IHOST,HOSTIMP(NATYPD)
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
       INTEGER IRCUT(0:IPAND,*),IRCUTIMP(0:IPAND,*)
       DOUBLE PRECISION C,Z(*),ZIMP(NATOMIMP),RCLSIMP(3,NATOMIMPD)
       DOUBLE COMPLEX E
       DOUBLE PRECISION R(IRMD,*),CLEB(*),RIMP(IRMD,NATOMIMP)
       INTEGER ICLEB(NCLEB,4)
       DOUBLE PRECISION 
     +   VINS(IRMIND:IRMD,LMPOTD,NSPOTD),
     +   VM2Z(IRMD,NSPOTD),
     +   VINSIMP(IRMIND:IRMD,LMPOTD,NSPIN*NATOMIMP),
     +   VM2ZIMP(IRMD,NSPIN*NATOMIMP)
       DOUBLE COMPLEX  
c     +   PNS_SO(NSPD*LMMAXD,NSPD*LMMAXD,IRMD,2),
     +   TMATLL(NSPD*LMMAXD,NSPD*LMMAXD,IHOST),
c     +   TMATLLIMP(NSPD*LMMAXD,NSPD*LMMAXD,NATOMIMP),
c     +   DELTAIMP(NSPD*LMMAXD*NATOMIMP,NSPD*LMMAXD*NATOMIMP),
     +   DTMTRX(NSPD*LMMAXD*NATOMIMP,NSPD*LMMAXD*NATOMIMP)
       INTEGER I1,IR,IR2,NSRA,USE_SRATRICK,NVEC,LM1,LM2,IPOT,ISPIN,I2,
     +         IL1,IL2
       INTEGER NPAN_LOG,NPAN_EQ,NPAN_INST,NCHEB,NPAN_TOT,IRMDNEW
       DOUBLE PRECISION R_LOG,THETA,PHI
       DOUBLE COMPLEX GMATPREFACTOR
       DOUBLE PRECISION, ALLOCATABLE :: RNEW(:),RPAN_INTERVALL(:)
       INTEGER, ALLOCATABLE :: IPAN_INTERVALL(:)
       DOUBLE PRECISION, ALLOCATABLE :: VINSNEW(:,:,:)
       DOUBLE COMPLEX,ALLOCATABLE :: VNSPLL0(:,:,:),VNSPLL(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE :: HLK(:,:),JLK(:,:),
     +                                HLK2(:,:),JLK2(:,:),
     +                                SLL(:,:,:),
     +                                RLLIMP(:,:,:,:),RLLHOST(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: VNSIMP(:,:,:,:),VNSHOST(:,:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: DELTABG(:,:,:),DELTASM(:,:,:)
       DOUBLE COMPLEX, ALLOCATABLE:: RADIALHOST(:,:),RADIALIMP(:,:),
     +                               VLLIMP(:,:),DELTAV(:,:),
     +                               TMATLLIMP(:,:,:),DELTAMTR(:,:,:),
     +                               DELTATMP(:),DELTAIMP(:,:)
       DOUBLE COMPLEX TMATSPH(2*(LMAX+1))
       INTEGER JLK_INDEX(2*LMMAXSO)

       CHARACTER*8 TESTC(16),OPTC(8)
       COMMON /TESTC/TESTC
       COMMON /OPTC/OPTC
       LOGICAL TEST          ! added for test('DltBorn ') by Philipp, 18.01.16
       EXTERNAl TEST         ! added for test('DltBorn ') by Philipp, 18.01.16

       WRITE(6,*) 'in tmatimp'
       IF (KSRA.GE.1) THEN
        NSRA=2
       ELSE
        NSRA=1
       ENDIF
       NPAN_TOTD= NPAN_LOG+NPAN_EQ+IPAND-1
c       NPAN_TOTD= NPAN_LOG+NPAN_EQ+IPAN(1)-1
       ALLOCATE(VNSIMP(NSRA*LMMAXSO,NSRA*LMMAXSO,
     +                 NPAN_TOTD*(NCHEB+1),NATOMIMP))
       ALLOCATE(VNSHOST(NSRA*LMMAXSO,NSRA*LMMAXSO,
     +                  NPAN_TOTD*(NCHEB+1),IHOST))
       ALLOCATE(RLLIMP(NSRA*LMMAXSO,LMMAXSO,
     +                 NPAN_TOTD*(NCHEB+1),NATOMIMP))
       ALLOCATE(RLLHOST(NSRA*LMMAXSO,LMMAXSO,
     +                 NPAN_TOTD*(NCHEB+1),IHOST))
       ALLOCATE(SLL(NSRA*LMMAXSO,LMMAXSO,NPAN_TOTD*(NCHEB+1)))
       ALLOCATE(DELTABG(LMMAXSO,LMMAXSO,NPAN_TOTD*(NCHEB+1)))
       ALLOCATE(DELTASM(LMMAXSO,LMMAXSO,NPAN_TOTD*(NCHEB+1)))
       ALLOCATE(DELTAV(LMMAXSO,LMMAXSO))
       ALLOCATE(TMATLLIMP(LMMAXSO,LMMAXSO,NATOMIMP))
       ALLOCATE(DELTAMTR(LMMAXSO,LMMAXSO,NATOMIMP))
       TMATLL=CZERO
       TMATLLIMP=CZERO
       DELTAMTR=CZERO
       VNSIMP=CZERO
       VNSHOST=CZERO
       RLLHOST=CZERO
       RLLIMP=CZERO

       OPEN(UNIT=11,FILE='nonco_angle.dat',FORM='FORMATTED')
       OPEN(UNIT=12,FILE='nonco_angle_imp.dat',FORM='FORMATTED')

c calculate tmat and radial wavefunctions of host atoms      
      DO I1=1,NATYPD
       READ(11,*) THETA,PHI
       DO I2=1,IHOST
        IF (I1.EQ.HOSTIMP(I2)) THEN
         ISPIN=1
         IPOT=NSPIN*(I1-1)+ISPIN
         WRITE(6,*) 'HOST',I2,I1,IPOT

c create new mesh
c data for the new mesh
       NPAN_INST= IPAN(I1)-1
       NPAN_TOT= NPAN_LOG+NPAN_EQ+NPAN_INST
       IRMDNEW= NPAN_TOT*(NCHEB+1)
c new mesh
       ALLOCATE(RNEW(IRMDNEW))
       ALLOCATE(RPAN_INTERVALL(0:NPAN_TOT))
       ALLOCATE(IPAN_INTERVALL(0:NPAN_TOT))
       ALLOCATE(VINSNEW(IRMDNEW,LMPOTD,NSPIND))
       CALL CREATE_NEWMESH(INS,NSPIN,R(1,I1),LMMAXD,IRMIND,IRWS(I1),
     +                    IRMD,LMPOTD,IPAN(I1),IRCUT(0,I1),IRMIN(I1),
     +                    VINS(IRMIND:IRMD,1:LMPOTD,IPOT:IPOT+NSPIN-1),
     +                    VM2Z(1:IRMD,IPOT:IPOT+NSPIN-1),R_LOG,
     +                    NPAN_LOG,NPAN_EQ,NPAN_INST,NCHEB,NPAN_TOT,
     +                    RNEW,RPAN_INTERVALL,IPAN_INTERVALL,VINSNEW)
c set up the non-spherical ll' matrix for potential VLL'
       USE_SRATRICK=1
       ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW))
       VNSPLL0=CZERO
       CALL VLLMAT(1,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINSNEW,
     +             CLEB,ICLEB,IEND,NSPIN,Z(I1),RNEW,USE_SRATRICK,I1)

c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW,RNEW,
     +                    E,Z(I1),C,NSRA,NSPIN,LMPOTD,
     +                    IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,
     +                    NPAN_TOT*(NCHEB+1),VNSPLL0,THETA,PHI,'1')
c extend matrix for the SRA treatment
       IF (NSRA.EQ.2) THEN
        ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW))
        IF (USE_SRATRICK.EQ.0) THEN
         CALL VLLMATSRA(VNSPLL0,VNSPLL,RNEW,LMMAXSO,IRMDNEW,E,C,
     +                 LMAX,0,'Ref=0')
        ELSEIF (USE_SRATRICK.EQ.1) THEN
         CALL VLLMATSRA(VNSPLL0,VNSPLL,RNEW,LMMAXSO,IRMDNEW,E,C,
     +                 LMAX,0,'Ref=Vsph')
        ENDIF
       ELSE
        ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW))
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
       CALL CALCSPH(NSRA,IRMDNEW,LMAX,NSPIN,Z(I1),C,E,LMPOTD,LMMAXSO,
     +              RNEW,VINSNEW,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +              JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,TMATSPH,
     +              USE_SRATRICK)
       ENDIF

c calculate the tmat and wavefunctions
       SLL=CZERO

c right solutions
       CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL,RLLHOST(1,1,1,I2),
     +             SLL(1,1,1),TMATLL(1,1,I2),NCHEB,
     +             NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +             NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,
     +             '1','1','0',USE_SRATRICK)
       IF (NSRA.EQ.2) THEN
        DO IR=1,IRMDNEW
         DO LM1=1,LMMAXSO
          DO LM2=1,LMMAXSO
           RLLHOST(LM1+LMMAXSO,LM2,IR,I2)=
     +             RLLHOST(LM1+LMMAXSO,LM2,IR,I2)/C
           SLL(LM1+LMMAXSO,LM2,IR)=
     +             SLL(LM1+LMMAXSO,LM2,IR)/C
          ENDDO
         ENDDO
        ENDDO
       ENDIF
c add spherical contribution of tmatrix

       IF (USE_SRATRICK.EQ.1) THEN
        DO LM1=1,NSPD*LMMAXD
         TMATLL(LM1,LM1,I2)=TMATLL(LM1,LM1,I2)+
     &                         TMATSPH(JLK_INDEX(LM1))
        ENDDO
       ENDIF

c rotate tmatrix and radial wavefunction to global frame
       CALL ROTATEMATRIX(TMATLL(1,1,I2),THETA,PHI,LMMAXD,0)

c create SRA potential for host
c set up the non-spherical ll' matrix for potential VLL'
       VNSPLL0=CZERO
       CALL VLLMAT(1,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINSNEW,
     +             CLEB,ICLEB,IEND,NSPIN,Z(I1),RNEW,0,I1)

c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW,RNEW,
     +                    E,Z(I1),C,NSRA,NSPIN,LMPOTD,
     +                    IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,
     +                    NPAN_TOT*(NCHEB+1),VNSPLL0,THETA,PHI,'1')
       DO IR=1,IRMDNEW
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          VNSHOST(LM1,LM2,IR,I2)=VNSPLL0(LM1,LM2,IR)
          VNSHOST(LM1+LMMAXSO,LM2+LMMAXSO,IR,I2)=VNSPLL0(LM1,LM2,IR)
         ENDDO
        ENDDO
       ENDDO

          DEALLOCATE(RNEW)
          DEALLOCATE(RPAN_INTERVALL)
          DEALLOCATE(IPAN_INTERVALL)
          DEALLOCATE(VINSNEW)
          DEALLOCATE(VNSPLL0)
          DEALLOCATE(VNSPLL)
          DEALLOCATE(HLK)
          DEALLOCATE(JLK)
          DEALLOCATE(HLK2)
          DEALLOCATE(JLK2)
         ENDIF
        ENDDO ! I1
       ENDDO ! I2

c calculate tmat and radial wavefunctions of impurity atoms      
       IF (IELAST.EQ.1) THEN
        OPEN(UNIT=20,FILE='DTMTRX',FORM='FORMATTED')
        WRITE(20,'(I5)') NATOMIMP
        DO I1=1,NATOMIMP
         WRITE(20,'(3e17.9)') (RCLSIMP(I2,I1),I2=1,3)
        ENDDO
       ENDIF
       DO I1=1,NATOMIMP
        READ(12,*) THETA,PHI
        ISPIN=1
        IPOT=NSPIN*(I1-1)+ISPIN
        WRITE(6,*) 'IMP',I1,IPOT

c create new mesh
c data for the new mesh
       NPAN_INST= IPANIMP(I1)-1
       NPAN_TOT= NPAN_LOG+NPAN_EQ+NPAN_INST
       IRMDNEW= NPAN_TOT*(NCHEB+1)
c new mesh
       ALLOCATE(RNEW(IRMDNEW))
       ALLOCATE(RPAN_INTERVALL(0:NPAN_TOT))
       ALLOCATE(IPAN_INTERVALL(0:NPAN_TOT))
       ALLOCATE(VINSNEW(IRMDNEW,LMPOTD,NSPIND))
       CALL CREATE_NEWMESH(INS,NSPIN,RIMP(1,I1),LMMAXD,IRMIND,
     +                 IRWSIMP(I1),IRMD,LMPOTD,
     +                 IPANIMP(I1),IRCUTIMP(0,I1),IRMINIMP(I1),
     +                 VINSIMP(IRMIND:IRMD,1:LMPOTD,IPOT:IPOT+NSPIN-1),
     +                 VM2ZIMP(1:IRMD,IPOT:IPOT+NSPIN-1),R_LOG,
     +                 NPAN_LOG,NPAN_EQ,NPAN_INST,NCHEB,NPAN_TOT,
     +                 RNEW,RPAN_INTERVALL,IPAN_INTERVALL,VINSNEW)
c set up the non-spherical ll' matrix for potential VLL'
       USE_SRATRICK=1
       ALLOCATE(VNSPLL0(LMMAXSO,LMMAXSO,IRMDNEW))
       VNSPLL0=CZERO
       CALL VLLMAT(1,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINSNEW,
     +             CLEB,ICLEB,IEND,NSPIN,ZIMP(I1),RNEW,USE_SRATRICK,I1)

c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW,RNEW,
     +                    E,ZIMP(I1),C,NSRA,NSPIN,LMPOTD,
     +                    IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,
     +                    NPAN_TOT*(NCHEB+1),VNSPLL0,THETA,PHI,'1')
c extend matrix for the SRA treatment
       IF (NSRA.EQ.2) THEN
        ALLOCATE(VNSPLL(2*LMMAXSO,2*LMMAXSO,IRMDNEW))
        IF (USE_SRATRICK.EQ.0) THEN
         CALL VLLMATSRA(VNSPLL0,VNSPLL,RNEW,LMMAXSO,IRMDNEW,E,C,
     +                 LMAX,0,'Ref=0')
        ELSEIF (USE_SRATRICK.EQ.1) THEN
         CALL VLLMATSRA(VNSPLL0,VNSPLL,RNEW,LMMAXSO,IRMDNEW,E,C,
     +                 LMAX,0,'Ref=Vsph')
        ENDIF
       ELSE
        ALLOCATE(VNSPLL(LMMAXSO,LMMAXSO,IRMDNEW))
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
       CALL CALCSPH(NSRA,IRMDNEW,LMAX,NSPIN,ZIMP(I1),C,E,LMPOTD,LMMAXSO,
     +              RNEW,VINSNEW,NCHEB,NPAN_TOT,RPAN_INTERVALL,
     +              JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,TMATSPH,
     +              USE_SRATRICK)
       ENDIF

c calculate the tmat and wavefunctions
       SLL=CZERO

c right solutions
       CALL RLLSLL(RPAN_INTERVALL,RNEW,VNSPLL,RLLIMP(1,1,1,I1),
     +             SLL(1,1,1),TMATLLIMP(1,1,I1),NCHEB,
     +             NPAN_TOT,LMMAXSO,NVEC*LMMAXSO,4*(LMAX+1),IRMDNEW,
     +             NSRA,JLK_INDEX,HLK,JLK,HLK2,JLK2,GMATPREFACTOR,
     +             '1','1','0',USE_SRATRICK)
       IF (NSRA.EQ.2) THEN
        DO IR=1,IRMDNEW
         DO LM1=1,LMMAXSO
          DO LM2=1,LMMAXSO
           RLLIMP(LM1+LMMAXSO,LM2,IR,I1)=
     +            RLLIMP(LM1+LMMAXSO,LM2,IR,I1)/C
           SLL(LM1+LMMAXSO,LM2,IR)=
     +            SLL(LM1+LMMAXSO,LM2,IR)/C
          ENDDO
         ENDDO
        ENDDO
       ENDIF
            
c add spherical contribution of tmatrix
       
       IF (USE_SRATRICK.EQ.1) THEN
        DO LM1=1,NSPD*LMMAXD
         TMATLLIMP(LM1,LM1,I1)=TMATLLIMP(LM1,LM1,I1)+
     &                         TMATSPH(JLK_INDEX(LM1))
        ENDDO
       ENDIF

c rotate tmatrix and radial wavefunction to global frame

       CALL ROTATEMATRIX(TMATLLIMP(1,1,I1),THETA,PHI,LMMAXD,0)

c create SRA potential for impurity
c set up the non-spherical ll' matrix for potential VLL'
       VNSPLL0=CZERO
       CALL VLLMAT(1,IRMDNEW,LMMAXD,LMMAXSO,VNSPLL0,VINSNEW,
     +             CLEB,ICLEB,IEND,NSPIN,ZIMP(I1),RNEW,0,I1)

c contruct the spin-orbit coupling hamiltonian and add to potential
       CALL SPINORBIT_HAM(LMAX,LMMAXD,VINSNEW,RNEW,
     +                    E,ZIMP(I1),C,NSRA,NSPIN,LMPOTD,
     +                    IPAN_INTERVALL,RPAN_INTERVALL,NPAN_TOT,NCHEB,
     +                    NPAN_TOT*(NCHEB+1),VNSPLL0,THETA,PHI,'1')
       DO IR=1,IRMDNEW
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          VNSIMP(LM1,LM2,IR,I1)=VNSPLL0(LM1,LM2,IR)
          VNSIMP(LM1+LMMAXSO,LM2+LMMAXSO,IR,I1)=VNSPLL0(LM1,LM2,IR)
         ENDDO
        ENDDO
       ENDDO

c calculate delta_t_imp matrix written in TMATLLIMP

        DO I2=1,IHOST
         IF (ATOMIMP(I1).EQ.HOSTIMP(I2)) THEN
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            TMATLLIMP(LM1,LM2,I1)=TMATLLIMP(LM1,LM2,I1)-
     +                            TMATLL(LM1,LM2,I2)
           ENDDO
          ENDDO
          DO LM1=1,NSRA*LMMAXSO
           DO LM2=1,NSRA*LMMAXSO
            DO IR=1,IRMDNEW
             VNSIMP(LM1,LM2,IR,I1)=VNSIMP(LM1,LM2,IR,I1)-
     +                             VNSHOST(LM1,LM2,IR,I2)
            ENDDO
           ENDDO
          ENDDO
         ENDIF
        ENDDO

c calculate delta matrix
 
       IF (IELAST.EQ.1) THEN
        DELTABG=CZERO
        DELTASM=CZERO
        ALLOCATE(DELTATMP(IRMDNEW))
        ALLOCATE(RADIALHOST(LMMAXSO,LMMAXSO))
        ALLOCATE(RADIALIMP(LMMAXSO,LMMAXSO))
        ALLOCATE(VLLIMP(LMMAXSO,LMMAXSO))
        DELTATMP=CZERO
c big component for SRA
        DO IR=1,IRMDNEW
         RADIALHOST=CZERO
         RADIALIMP=CZERO
         VLLIMP=CZERO 
         DELTAV=CZERO
         DO LM1=1,LMMAXSO
          DO LM2=1,LMMAXSO
           DO I2=1,IHOST
            IF (ATOMIMP(I1).EQ.HOSTIMP(I2)) THEN
             RADIALHOST(LM1,LM2)=RLLHOST(LM1,LM2,IR,I2)
            ENDIF
           ENDDO
           RADIALIMP(LM1,LM2)=RLLIMP(LM1,LM2,IR,I1)
           VLLIMP(LM1,LM2)=VNSIMP(LM1,LM2,IR,I1)
          ENDDO
         ENDDO 
c
         ! set impurity wavefunction to host wavefunction to calculate
         ! Born approximation, i.e. calculate T^Born=c.[Delta.(1+0)].c with Delta=R^*.(V^imp-V^host).R
         ! instead of T=c.[Delta.(1+G^imp.Dt^imp)].c^imp with Delta=R^*.(V^imp-V^host).R^imp
         IF(TEST('DltBorn ')) RADIALIMP = RADIALHOST
c
         CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VLLIMP,
     &              LMMAXSO,RADIALIMP,LMMAXSO,CZERO,DELTAV,LMMAXSO)
c         CALL ZGEMM('T','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
c     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTABG(1,1,IR),LMMAXSO)
         CALL ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTABG(1,1,IR),LMMAXSO)
c small component for SRA         
         IF (NSRA.EQ.2) THEN
         RADIALHOST=CZERO
         RADIALIMP=CZERO
         VLLIMP=CZERO
         DELTAV=CZERO
         DO LM1=1,LMMAXSO
          DO LM2=1,LMMAXSO
           DO I2=1,IHOST
            IF (ATOMIMP(I1).EQ.HOSTIMP(I2)) THEN
             RADIALHOST(LM1,LM2)=RLLHOST(LM1+LMMAXSO,LM2,IR,I2)
            ENDIF
           ENDDO
           RADIALIMP(LM1,LM2)=RLLIMP(LM1+LMMAXSO,LM2,IR,I1)
           VLLIMP(LM1,LM2)=VNSIMP(LM1+LMMAXSO,LM2+LMMAXSO,IR,I1)
          ENDDO
         ENDDO 
         CALL ZGEMM('N','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,VLLIMP,
     &              LMMAXSO,RADIALIMP,LMMAXSO,CZERO,DELTAV,LMMAXSO)
c         CALL ZGEMM('T','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
c     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTASM(1,1,IR),LMMAXSO)
         CALL ZGEMM('C','N',LMMAXSO,LMMAXSO,LMMAXSO,CONE,RADIALHOST,
     &             LMMAXSO,DELTAV,LMMAXSO,CZERO,DELTASM(1,1,IR),LMMAXSO)
          DO LM1=1,LMMAXSO
           DO LM2=1,LMMAXSO
            DELTABG(LM1,LM2,IR)=DELTABG(LM1,LM2,IR)+DELTASM(LM1,LM2,IR)
           ENDDO
          ENDDO
         ENDIF ! NSRA
        ENDDO ! IR
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          DO IR=1,IRMDNEW
           DELTATMP(IR)=DELTABG(LM1,LM2,IR)
          ENDDO
          CALL INTCHEB_CELL(DELTATMP,DELTAMTR(LM1,LM2,I1),
     &            RPAN_INTERVALL,IPAN_INTERVALL,NPAN_TOT,NCHEB,IRMDNEW)
         ENDDO
        ENDDO

       DEALLOCATE(DELTATMP)
       DEALLOCATE(RADIALHOST)
       DEALLOCATE(RADIALIMP)
       DEALLOCATE(VLLIMP)
      
       ENDIF ! IELAST.EQ.1
          
       DEALLOCATE(RNEW)
       DEALLOCATE(RPAN_INTERVALL)
       DEALLOCATE(IPAN_INTERVALL)
       DEALLOCATE(VINSNEW)
       DEALLOCATE(VNSPLL0)
       DEALLOCATE(VNSPLL)
       DEALLOCATE(HLK)
       DEALLOCATE(JLK)
       DEALLOCATE(HLK2)
       DEALLOCATE(JLK2)

      ENDDO ! I1 impurity
      DO I1=1,NATOMIMP
       DO LM1=1,LMMAXSO
        DO LM2=1,LMMAXSO
         IL1=LMMAXSO*(I1-1)+LM1
         IL2=LMMAXSO*(I1-1)+LM2
         DTMTRX(IL1,IL2)=TMATLLIMP(LM1,LM2,I1)
        ENDDO
       ENDDO
      ENDDO

c write down to the file DTMTRX
      IF (IELAST.EQ.1) THEN
       ALLOCATE(DELTAIMP(NSPD*LMMAXD*NATOMIMP,NSPD*LMMAXD*NATOMIMP))
       DELTAIMP=CZERO
       DO I1=1,NATOMIMP
        DO LM1=1,LMMAXSO
         DO LM2=1,LMMAXSO
          IL1=LMMAXSO*(I1-1)+LM1
          IL2=LMMAXSO*(I1-1)+LM2
          DELTAIMP(IL1,IL2)=DELTAMTR(LM1,LM2,I1)
         ENDDO
        ENDDO
       ENDDO
       DO LM1=1,LMMAXSO*NATOMIMP
        DO LM2=1,LMMAXSO*NATOMIMP
          WRITE(20,'((2I5),(4e17.9))') LM2,LM1,DTMTRX(LM2,LM1),
     &                                 DELTAIMP(LM2,LM1)
        ENDDO
       ENDDO
       DEALLOCATE(DELTAIMP)
       WRITE(6,*) 'end of delta t'
       STOP  
      ENDIF ! IELAST.EQ.1




       CLOSE(11)
       CLOSE(12)
       CLOSE(20)
       DEALLOCATE(VNSHOST)
       DEALLOCATE(VNSIMP)
       DEALLOCATE(RLLHOST)
       DEALLOCATE(RLLIMP)
       DEALLOCATE(SLL)
       DEALLOCATE(TMATLLIMP)
       DEALLOCATE(DELTAMTR)
       END 

