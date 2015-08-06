      SUBROUTINE DECITSET(ALAT,BRAVSYS,EZ,IELAST,
     &                    NLBASIS,NRBASIS,FILELEFT,FILERIGHT,
     &                    INS,KVREL,KREL,NSPIN,KMROT,
     &                    VREF,RMTREF,NREF,REFPOT,
     &                    RC,CREL,RREL,
     &                    LEFTTINV,RIGHTTINV,VACFLAG,
     &                    NEMBD1,IEMXD,IRMD,IPAND,
     &                    LMAXD,LMGF0D,LMMAXD,LM2D,NSPIND)
C **********************************************************************
C *                                                                    *
C * This subroutine is thought as an alternative to the < decimaread > *
C * which requires an a priori calculated set of single-site matrices  *
C * over a fixed energy mesh. It is using the potential as written out *
C * in < outpothost > routine and determines the matrix                *
C *                                                                    *
C *            /         \-1   /             \-1                       *
C *            | Delta t |   = |  t    - t   |                         *
C *            \         /     \  sys    ref /                         *
C *                                                                    *
C * for the left and the right host, using the energy mesh as read in  *
C * from the input-file                                                *
C *                                                                    *
C *                                        v.popescu - munich, Dec 04  *
C *                                                                    *
C * Notes: - no charge moments are calculated -- thus this option CAN  *
C *          NOT be used in SCF calculations                           *
C *        - non-spherical case not implemented, neither LDA+U (al-    *
C *          though the interface to regsol in decitmat is supplied)   *
C *        - CPA case not implemented - requires BZ integration        *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalars arguments ..
      INTEGER IEMXD,NEMBD1,LMMAXD,IPAND,NSPIND,IRMD,LMAXD,LM2D,LMGF0D
      INTEGER IELAST,KMROT,NLBASIS,NRBASIS,NREF,INS,KVREL,KREL,NSPIN
      DOUBLE PRECISION ALAT
      CHARACTER*40 FILELEFT,FILERIGHT
C     ..
C     .. Array arguments ..
      INTEGER REFPOT(NEMBD1)
      DOUBLE PRECISION VREF(*),RMTREF(*),BRAVSYS(3,3)
      DOUBLE COMPLEX CREL(LMMAXD,LMMAXD),RREL(LMMAXD,LMMAXD),
     &               RC(LMMAXD,LMMAXD),EZ(IEMXD)
      DOUBLE COMPLEX LEFTTINV(LMMAXD,LMMAXD,NEMBD1,NSPIND,IEMXD),
     &               RIGHTTINV(LMMAXD,LMMAXD,NEMBD1,NSPIND,IEMXD)
      LOGICAL VACFLAG(2)
C     ..
C     .. Local scalars ..
      INTEGER IHOST,I,LL,MM,LNGSTRING,NQHOST,ILHOST
      INTEGER NHOST,ILOOPH
      INTEGER NQ,NT,IQOFF,ITOFF,IE,IH,IQH,IOQ,INFO
      INTEGER IPOT,I1,ISPIN,NSRA,LM1,LM2,IRC1,IREF
      INTEGER NTLEFT,NTRIGHT,NTHOST
C     .. LDA+U
      INTEGER IDOLDAU,LOPT
      DOUBLE PRECISION WLDAUAV
C     ..
      DOUBLE PRECISION EFERMI,RIRC
      DOUBLE COMPLEX ERYD,CARG,CFCTOR
      LOGICAL TEST
      CHARACTER*40 FILEHOST
      CHARACTER*10 SOLVER
C     ..
C     .. Local arrays
      INTEGER KRELH(2),NSPINH(2),INSH(2),IPVT(LMMAXD)
      INTEGER NOQ(NEMBD1),KAOEZ(NEMBD1,NEMBD1),INHOST(2)
      DOUBLE PRECISION BRAVAIS(3,3,2),RBASIS(3,NEMBD1),QMTET(NEMBD1),
     &                 QMPHI(NEMBD1)
      CHARACTER*5 CHHOST(2)
      CHARACTER*9 TXTS(2)
C     ..
C     .. Allocatable local arrays 
CF90--------------------------------------------------------------------
      INTEGER NTMAX
      DOUBLE PRECISION ZAT(:),RWS(:),RMT(:),CONC(:)
      DOUBLE PRECISION RR(:,:),DRDI(:,:),VISP(:,:),DROR(:,:)
      DOUBLE PRECISION SOCSCL(:,:),CSCL(:,:)
      INTEGER IRWS(:),IPAN(:),IQAT(:,:),IRCUT(:,:),LOFLM(:)
      DOUBLE COMPLEX TREFLL(:,:,:),TMATLL(:,:),DHMAT(:,:,:)
      DOUBLE COMPLEX DTREFLL(:,:,:) ! LLY Lloyd
      DOUBLE COMPLEX ALPHAREF(:,:),DALPHAREF(:,:) ! LLY Lloyd Alpha matrix and deriv.
      DOUBLE COMPLEX WN1(:,:)
      DOUBLE PRECISION VTREL(:,:),BTREL(:,:),R2DRDIREL(:,:)
      INTEGER ZREL(:)
      ALLOCATABLE ZAT,RWS,RMT,CONC,RR,DRDI,VISP,DROR,SOCSCL,CSCL
      ALLOCATABLE VTREL,BTREL,R2DRDIREL
      ALLOCATABLE IRWS,IPAN,IQAT,IRCUT,LOFLM,ZREL
      ALLOCATABLE TREFLL,TMATLL,DHMAT,WN1
      ALLOCATABLE DTREFLL,ALPHAREF,DALPHAREF ! LLY
CF90--------------------------------------------------------------------
Comment previous lines and uncomment the following ones for F77 
CF77--------------------------------------------------------------------
CF77      INTEGER NTMAX
CF77      PARAMETER ( NTMAX = 8 )
CF77      DOUBLE PRECISION ZAT(NTMAX),RWS(NTMAX),RMT(NTMAX),CONC(NTMAX)
CF77      DOUBLE PRECISION RR(IRMD,NTMAX),DRDI(IRMD,NTMAX),
CF77     &                 VISP(IRMD,NTMAX*NSPIND),DROR(IRMD,NTMAX)
CF77      INTEGER IRWS(NTMAX),IPAN(NTMAX),IQAT(NEMBD1,NTMAX),
CF77     &        IRCUT(0:IPAND,NTMAX)
CF77      INTEGER LOFLM(LM2D)
CF77      DOUBLE COMPLEX TREFLL(LMMAXD,LMMAXD,NTMAX)
CF77      DOUBLE COMPLEX WN1(KREL*LMGF0D+1-KREL,KREL*LMGF0D+1-KREL)
CF77      DOUBLE COMPLEX TMATLL(LMMAXD,LMMAXD),DHMAT(LMMAXD,LMMAXD,2)
CF77      DOUBLE PRECISION SOCSCL(KREL*LMAXD+1,KREL*NTMAX+(1-KREL))
CF77      DOUBLE PRECISION CSCL(KREL*LMAXD+1,KREL*NTMAX+(1-KREL))
CF77      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL),NTMAX)
CF77      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL),NTMAX)
CF77      DOUBLE PRECISION R2DRDIREL(IRMD*KREL+(1-KREL),NTMAX)
CF77      INTEGER ZREL(NTMAX)
CF77----------------------------------------------------------------
C     .. 
C     .. External subroutines
      EXTERNAL CALCTREF13,CHANGEREP,CINIT,CMATSTR,DECIPOTBAS,
     &         DECIPOTHEAD,DECITMAT,ZAXPY,ZCOPY,ZGETRF,ZGETRI
C     ..
C     .. External Functions ..
      EXTERNAL LNGSTRING,TEST
C     ..
C     .. Data statements
      DATA CHHOST/'LEFT ','RIGHT'/
      DATA TXTS /'spin   UP','spin DOWN'/
C ......................................................................
C
CF77--------------------------------------------------------------------
Cccc      IF ( NREF.GT.NTMAX ) THEN
Cccc         WRITE(6,99001) 'local','NTMAX',NREF
Cccc         STOP
Cccc      END IF
CF77--------------------------------------------------------------------
CF90--------------------------------------------------------------------
      CFCTOR = ALAT/(8.D0*ATAN(1.0D0))           ! = ALAT/(2*PI)
C
      IDOLDAU = 0
      LOPT = -1
      WLDAUAV = 0D0
      ALLOCATE(LOFLM(LM2D),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate LOFLM'
CF90--------------------------------------------------------------------
      WRITE (6,'(5X,A,/,8X,65(1H-))')
     &           'Reading in host potentials'
      VACFLAG(1) = .FALSE.
      VACFLAG(2) = .FALSE.
      NSRA = 1
      IF ( KVREL.GE.1 ) NSRA = 2
      I = 1
      DO LL = 0,2*LMAXD
         DO MM = -LL,LL
            LOFLM(I) = LL
            I = I + 1
         END DO
      END DO
      NTLEFT = 0
      NTRIGHT = 0
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      NHOST = 0
      DO IHOST = 1,2
         FILEHOST = FILELEFT
         NQHOST = NLBASIS
         IF ( IHOST.EQ.2 ) THEN
            FILEHOST = FILERIGHT
            NQHOST = NRBASIS
         END IF
         ILHOST = LNGSTRING(FILEHOST,40)
         CALL DECIPOTHEAD(IHOST,FILEHOST,ILHOST,NQHOST,
     &                    VACFLAG,ALAT,BRAVSYS,NQ,NT,BRAVAIS(1,1,IHOST),
     &                    EFERMI,INSH(IHOST),KRELH(IHOST),NSPINH(IHOST),
     &                    INS,KREL,NSPIN,KMROT)
C
         IF ( .NOT.VACFLAG(IHOST) ) THEN
            NHOST = NHOST + 1
            INHOST(NHOST) = IHOST
            IF ( IHOST.EQ. 1 ) THEN
               NTLEFT = NT
            ELSE
               NTRIGHT = NT
            END IF
         END IF
C
C
CF77--------------------------------------------------------------------
Cccc         IF ( NTLEFT+NTRIGHT.GT.NTMAX ) THEN
Cccc            WRITE(6,99001) 'local','NTMAX',NT
Cccc            STOP
Cccc99001    FORMAT (6X,'Dimension ERROR: please increase the ',A
Cccc     &        ,' parameter',/,6X,A,' to a value >=',I5,/)
Cccc         END IF
CF77--------------------------------------------------------------------
      END DO
C
      IF ( NTLEFT+NTRIGHT.LE.0 ) THEN
         WRITE (6,'(8X,"Vacuum will be considered on both sides",/,
     &              8X,65(1H-))')
         RETURN
      END IF
C
CF90--------------------------------------------------------------------
      NTMAX = NTLEFT+NTRIGHT
      ALLOCATE(ZAT(NTMAX),RWS(NTMAX),RMT(NTMAX),CONC(NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate ZAT/RWS/RMT/CONC'
      ALLOCATE(RR(IRMD,NTMAX),DRDI(IRMD,NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate RR/DRDI'
      ALLOCATE(VISP(IRMD,NTMAX*NSPIND),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate VISP'
      ALLOCATE(IRWS(NTMAX),IPAN(NTMAX),IQAT(NEMBD1,NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate IRWS/IPAN/IQAT'
      ALLOCATE(IRCUT(0:IPAND,NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate IRCUT'
      ALLOCATE(SOCSCL(KREL*LMAXD+1,KREL*NTMAX+(1-KREL)),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate SOCSCL'
      ALLOCATE(CSCL(KREL*LMAXD+1,KREL*NTMAX+(1-KREL)),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate CSCL'
      ALLOCATE(VTREL(IRMD*KREL+(1-KREL),NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate VTREL'
      ALLOCATE(BTREL(IRMD*KREL+(1-KREL),NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate BTREL'
CF90--------------------------------------------------------------------
C
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      DO IHOST = 1,2
         WRITE (6,'(8X,A5," side host: ",$)') CHHOST(IHOST)
         IQOFF = 0
         ITOFF = 0
         NQHOST = NLBASIS
         NTHOST = NTLEFT
         FILEHOST = FILELEFT
         IF ( IHOST.EQ.2 ) THEN
            NQHOST = NRBASIS
            NTHOST = NTRIGHT
            IQOFF = NLBASIS
            ITOFF = NTLEFT
            FILEHOST = FILERIGHT
         END IF
         ILHOST = LNGSTRING(FILEHOST,40)
C
         IF ( FILEHOST(1:7).EQ.'vacuum' ) THEN
            WRITE (6,'(A,/,8X,65(1H-))') 'VACUUM will be used'
         ELSE
            WRITE (6,'(A,/)') FILEHOST(1:ILHOST)
            WRITE (6,99005) KRELH(IHOST),NSPINH(IHOST),INSH(IHOST),
     &                      KMROT,NQHOST,ALAT,EFERMI
            WRITE (6,99006) ((BRAVAIS(LL,MM,IHOST),MM=1,3),LL=1,3)
            CALL DECIPOTBAS(IHOST,IQOFF,ITOFF,NQHOST,NTHOST,
     &                      RBASIS,QMTET,QMPHI,NOQ,KAOEZ,
     &                      ZAT,IQAT,CONC,IRWS,IPAN,IRCUT,
     &                      RR,DRDI,VISP,NSPINH(IHOST),KRELH(IHOST),
     &                      SOLVER,SOCSCL,CSCL,VTREL,BTREL,
     &                      IRMD,IPAND,NEMBD1,NTMAX,NSPIND,LMAXD)
            WRITE (6,'(8X,65(1H-))')
         END IF
      END DO
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
CF90--------------------------------------------------------------------
      ALLOCATE(DROR(IRMD,NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate DROR'
CF90--------------------------------------------------------------------
      ALLOCATE(R2DRDIREL(IRMD*KREL+(1-KREL),NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate R2DRDIREL'
      ALLOCATE(ZREL(NTMAX),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate ZREL'
CF90--------------------------------------------------------------------
C
      WRITE (6,'(/,5X,A,/)') 'Calculating host (Delta_t)^(-1) matrices'
      IF ( KREL.EQ.0 ) THEN
         DO I = 1,NTLEFT+NTRIGHT
            IRC1 = IRCUT(IPAN(I),I)
            DO I1 = 2,IRC1
               DROR(I1,I) = DRDI(I1,I)/RR(I1,I)
            END DO
         END DO
      ELSE
         DO I = 1,NTLEFT+NTRIGHT
            IRC1 = IRCUT(IPAN(I),I)
            DO I1 = 1,IRC1
               R2DRDIREL(I1,I) = RR(I1,I)*RR(I1,I)*DRDI(I1,I)
            END DO
            ZREL(I) = NINT(ZAT(I))
         END DO
      END IF
C
C ******************************************************* energy loop IE
CF90--------------------------------------------------------------------
      ALLOCATE(TREFLL(LMMAXD,LMMAXD,NREF),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate TREFLL'
      ALLOCATE(DTREFLL(LMMAXD,LMMAXD,NREF),STAT=I1)     ! LLY
      IF ( I1.NE.0 ) STOP '    Allocate DTREFLL'        ! LLY
      ALLOCATE(TMATLL(LMMAXD,LMMAXD),DHMAT(LMMAXD,LMMAXD,2),STAT=I1)
      IF ( I1.NE.0 ) STOP '    Allocate TMATLL/DHMAT'
      ALLOCATE( ALPHAREF(0:LMAXD,NREF),DALPHAREF(0:LMAXD,NREF) ,STAT=I1) ! LLY Lloyd Alpha matrix and deriv.
      IF ( I1.NE.0 ) STOP '    Allocate ALPHAREF/DALPHAREF'

CF90--------------------------------------------------------------------
      DO IE = 1,IELAST
         ERYD = EZ(IE)

C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C -> set up t matrices for the reference system
C
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
         IF ( KREL.EQ.0 ) THEN
            DO I1 = 1,NREF
               CALL CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAXD,IH,
     &                   TREFLL(1,1,I1),DTREFLL(1,1,I1),
     &                   ALPHAREF(0,I1),DALPHAREF(0,I1),LMAXD+1,LMGF0D)
            END DO
!          ELSE
! CF90--------------------------------------------------------------------
!             ALLOCATE(WN1(LMGF0D,LMGF0D),STAT=I1)
!             IF ( I1.NE.0 ) STOP '    Allocate WN1'
! CF90--------------------------------------------------------------------
!             DO I1 = 1,NREF
!                CALL CALCTREF13(ERYD,VREF(I1),RMTREF(I1),LMAXD,IH,
!      &                       WN1,DTREFLL(1,1,I1),LMAXD+1,LMGF0D)
! C-----------------------------------------------------------------
! C add second spin-block for relativistic calculation and transform
! C from NREL to REL representation
! C-----------------------------------------------------------------
!                CALL CINIT(LMMAXD*LMMAXD,TMATLL)
!                IF ( LMMAXD.NE.IH*2 ) STOP 'LMMAXD <> IH*2 '
!                DO I=1,IH
!                   TMATLL(I,I) = WN1(I,I)
!                   TMATLL(IH+I,IH+I) = WN1(I,I)
!                END DO
!                CALL CHANGEREP(TMATLL,'RLM>REL',TREFLL(1,1,I1),LMMAXD,
!      &                        LMMAXD,RC,CREL,RREL,'TREFLL',0)
!             END DO
! CF90--------------------------------------------------------------------
!             DEALLOCATE(WN1,STAT=I1)
!             IF ( I1.NE.0 ) STOP '    Deallocate WN1'
! CF90--------------------------------------------------------------------
         END IF
C &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
         DO ILHOST = 1,NHOST
            IHOST = INHOST(ILHOST)
            IQOFF = 0
            ITOFF = 0 
            NQHOST = NLBASIS
            IF ( IHOST.EQ.2 ) THEN
               NQHOST = NRBASIS
               IQOFF = NLBASIS
               ITOFF = NTLEFT
            END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ sites in host
            DO IH = 1,NQHOST
C
C -> assign Delta_t = -t_ref
C
C  Note: REFPOT(1) = REFPOT(NAEZ) (i.e., of the 2D system)
C
               IQH = IQOFF + IH
               IREF = REFPOT(IQH+1)
               DO LM2 = 1,LMMAXD
                  DO LM1 = 1,LMMAXD
                     DHMAT(LM1,LM2,1) = -TREFLL(LM1,LM2,IREF)
                  END DO
               END DO
C
               IF ( NSPINH(IHOST).GT.1 ) THEN
                  DO LM2 = 1,LMMAXD
                     CALL ZCOPY(LMMAXD,DHMAT(1,LM2,1),1,
     &                                 DHMAT(1,LM2,2),1)
                  END DO
               END IF
C ====================================================== spins and atoms
               DO ISPIN = 1,NSPINH(IHOST)
C ----------------------------------------------------------------------
                  DO IOQ = 1,NOQ(IQH)
                     I1 = KAOEZ(IOQ,IQH) + ITOFF
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C -> calculate t_sys for the atom I1 located on site IH
C
                     IPOT = (I1-1) * NSPINH(IHOST) + ISPIN
                     IRC1 = IRCUT(IPAN(I1),I1)
                     RIRC = RR(IRC1,I1)
C     
                     CALL DECITMAT(ERYD,ZAT(I1),IPAN(I1),RR(1,I1),
     &                             DROR(1,I1),VISP(1,IPOT),IRCUT(0,I1),
     &                             RIRC,KREL,NSRA,INS,TMATLL,LOFLM,
     &                             IDOLDAU,LOPT,WLDAUAV,
     &                             SOLVER,SOCSCL(1,KREL*I1+(1-KREL)),
     &                             CSCL(1,KREL*I1+(1-KREL)),ZREL(I1),
     &                             VTREL(1,I1),BTREL(1,I1),
     &                             DRDI(1,I1),R2DRDIREL(1,I1),
     &                             IPAND,IRMD,LMAXD,LMAXD+1,LM2D,LMMAXD)
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tmat calculated
C
C -> Delta_t = Delta_t + CONC(I1)*t_mat(I1)
C                  
                     CARG = CONC(I1)
                     DO LM2 = 1,LMMAXD
                        CALL ZAXPY(LMMAXD,CARG,TMATLL(1,LM2),1,
     &                             DHMAT(1,LM2,ISPIN),1)
                     END DO
                  END DO
C ----------------------------------------------------------------------
C tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
                  IF ( TEST('tmat    ') ) THEN
                     WRITE (6,*)
                     WRITE (6,99002) 
     &                    '      ---> Delta_t  matrix for site: ',IQH
                     IF ( KREL.EQ.0 ) WRITE (6,99003) TXTS(ISPIN)
                     WRITE (6,99004) ', energy: ',ERYD
                     CALL CMATSTR(' ',1,DHMAT(1,1,ISPIN),LMMAXD,LMMAXD,
     &                            2*KREL+1,2*KREL+1,0,1D-8,6)
                     WRITE (6,*)
                  END IF
C tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
C
C --> inversion 
C
                  CALL ZGETRF(LMMAXD,LMMAXD,DHMAT(1,1,ISPIN),LMMAXD,
     &                        IPVT,INFO)
                  CALL ZGETRI(LMMAXD,DHMAT(1,1,ISPIN),LMMAXD,IPVT,
     &                        TMATLL,LMMAXD*LMMAXD,INFO)
               END DO
C ======================================================================
C 
C --> scaling the host t-matrices to p.u. 
C 
               IF ( IHOST.EQ.1 ) THEN
                  DO ISPIN = 1,NSPINH(IHOST)
                     DO LM2 = 1,LMMAXD
                        DO LM1 = 1,LMMAXD
                           LEFTTINV(LM1,LM2,IH,ISPIN,IE) = CFCTOR *
     &                                              DHMAT(LM1,LM2,ISPIN)
                        END DO
                     END DO
                  END DO
               ELSE
                  DO ISPIN = 1,NSPINH(IHOST)
                     DO LM2 = 1,LMMAXD
                        DO LM1 = 1,LMMAXD
                           RIGHTTINV(LM1,LM2,IH,ISPIN,IE) = CFCTOR *
     &                                              DHMAT(LM1,LM2,ISPIN)
                        END DO
                     END DO
                  END DO
               END IF
            END DO
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END DO
C **********************************************************************
CF90--------------------------------------------------------------------
      DEALLOCATE(ZAT,RWS,RMT,CONC,RR,DRDI,VISP,STAT=I1)
      IF ( I1.NE.0 ) STOP '   Deallocate ZAT/RWS/RMT/.../VISP'
      DEALLOCATE(IRWS,IPAN,IQAT,IRCUT,LOFLM,STAT=I1)
      IF ( I1.NE.0 ) STOP '   Deallocate IRWS/IPAN/IQAT/IRCUT/LOFLM'
      DEALLOCATE(TREFLL,TMATLL,DHMAT,STAT=I1)
      IF ( I1.NE.0 ) STOP '   Deallocate TREFLL/TMATLL/DHMAT'
      DEALLOCATE(SOCSCL,CSCL,VTREL,BTREL,STAT=I1)
      IF ( I1.NE.0 ) STOP '   Deallocate SOCSCL/CSCL/VTREL/BTREL'
      DEALLOCATE(ALPHAREF,DALPHAREF ,STAT=I1)
      IF ( I1.NE.0 ) STOP '   Deallocate ALPHAREF/DALPHAREF'
      IF ( KREL.EQ.0 ) THEN
         DEALLOCATE(DROR,STAT=I1)
         IF ( I1.NE.0 ) STOP '   Deallocate DROR'
      ELSE
         DEALLOCATE(R2DRDIREL,ZREL,STAT=I1)
         IF ( I1.NE.0 ) STOP '   Deallocate R2DRDIREL/ZREL'
      END IF
CF90--------------------------------------------------------------------
99002 FORMAT (A,I3,$)
99003 FORMAT (', ',A,$)
99004 FORMAT (A,2F10.6)
99005 FORMAT(10X,'KREL= ',I1,' NSPIN= ',I1,' INS= ',I1,' KMROT= ',I1,/,
     &       10X,'NAEZ=',I3,' ALAT= ',F9.6,' EFERMI= ',F9.6)
99006 FORMAT(10X,'BRAVAIS '/10X,3F8.4/10X,3F8.4/10X,3F8.4/10X,'RBASIS')
      END
