C*==deciopt.f    processed by SPAG 6.05Rc at 19:30 on  5 Jun 2004
      SUBROUTINE DECIOPT(ALAT,INS,KREL,KVREL,KMROT,NSPIN,NAEZ,LMMAX,
     &                   BRAVAIS,TK,NPOL,NPNT1,NPNT2,NPNT3,
     &                   EZ,IELAST,KAOEZ,SCFSTEPS,
     &                   LEFTTINVLL,RIGHTTINVLL,VACFLAG,NLBASIS,NRBASIS,
     &                   CMOMHOST,VREF,RMTREF,NREF,REFPOT,RC,CREL,RREL,
     &                   LMAXD,LMGF0D,LMMAXD,LM2D,NEMBD1,IEMXD,NSPIND,
     &                   LMPOTD,NATYPD,IRMD,IPAND,KORBIT)
C **********************************************************************
C *                                                                    *
C * This routine treats the DECIMATION case setting up the single-site *
C * (Delta t)^(-1) matrices and the charge moments of the host(s).     *
C *                                                                    *
C * This is realised in two ways:                                      *
C *      - either reading in the matrices (and moments - if SCFSTEPS   *
C *        is greater than 1) as written out in a previous (bulk) run  *
C *        -- DECIFILES token points to the files containing the nece- *
C *        ssary information                                           *
C *      - or reading in the self-consistent potential for each host   *
C *        and effectively calculating the matrices; the potential     *
C *        must have the specific format set in < OUTPOTHOST > routine *
C *        and the DECIPOTS token should point to the corresponding    *
C *        potential file(s)                                           *
C *                                                                    *
C * Notes:                                                             *
C *        - DECIFILES token is considered by default and sought first *
C *                          is requiring the same energy mesh for the *
C *                          system as for the host                    *
C *        - DECIPOTS token is not restrictive in this sense           *
C *                         however, it does not calculate charge mo-  *
C *                         ments -- does not work in SCF mode         *
C *                         is not dealing with CPA host so far        *
C *                         is not dealing with NON-SPHERICAL poten-   *
C *                         tials so far                               *
C *                                                                    *
C *                                     v.popescu - munich, Dec 04     *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar arguments
      INTEGER LMMAXD,NEMBD1,IEMXD,NSPIND,LMPOTD,NATYPD,IPAND,IRMD,LMAXD
      INTEGER LM2D,NREF,LMGF0D,KORBIT           ! ruess: for tmat newsolver (SOC)
      INTEGER INS,KREL,KMROT,NSPIN,NAEZ,LMMAX,NPOL,NPNT1,NPNT2,NPNT3
      INTEGER IELAST,NLBASIS,NRBASIS,KVREL,SCFSTEPS
      DOUBLE PRECISION ALAT,TK
C     ..
C     .. Array arguments
      INTEGER KAOEZ(NATYPD,*),REFPOT(NEMBD1)
      DOUBLE PRECISION BRAVAIS(3,3),CMOMHOST(LMPOTD,*)
      DOUBLE PRECISION VREF(*),RMTREF(*)
      DOUBLE COMPLEX LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPIND,IEMXD),
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPIND,IEMXD)
      DOUBLE COMPLEX CREL(LMMAXD,LMMAXD),RREL(LMMAXD,LMMAXD),
     &               RC(LMMAXD,LMMAXD)
      DOUBLE COMPLEX EZ(IEMXD)
      LOGICAL VACFLAG(2)
C     ..
C     .. Local scalars
      INTEGER IERROR,IL,IE,ISPIN,NSPINSO           ! ruess: for tmat new solver
      DOUBLE COMPLEX CFCTOR
      CHARACTER*40 FILELEFT,FILERIGHT
      CHARACTER*256 UIO ! NCOLIO=256
C     ..                                  ! ruess: for NEWSOSOL running option
C     .. External Functions ..
      LOGICAL OPT
      EXTERNAL OPT
C     
C ======================================================================
      WRITE (1337,'(79(1H=))')
      WRITE (1337,'(15X,A,/,79(1H=),/)')
     &               'DECIOPT: reading left/right host decimation files'
      IL = 1
      IERROR = 0
      CALL IOINPUT('DECIFILES       ',UIO,IL,7,IERROR)
C :::::::::::::::::::::::::::::::::::::::::::::::: decifiles (tmatrices)
      IF ( IERROR.EQ.0 ) THEN
         READ (UNIT=UIO,FMT='(A40)') FILELEFT
         CALL IOINPUT('DECIFILES       ',UIO,IL+1,7,IERROR)
         READ (UNIT=UIO,FMT='(A40)') FILERIGHT
C ----------------------------------------------------------------------
C
C --> first call to read the header ( IE = 0 )
C
         IE = 0
         CALL DECIMAREAD(EZ,TK,NPNT1,NPNT2,NPNT3,NPOL,NSPIN,
     &        LEFTTINVLL(1,1,1,1,1),RIGHTTINVLL(1,1,1,1,1),VACFLAG,
     &        IE,NLBASIS,NRBASIS,NAEZ,KAOEZ,KMROT,
     &        INS,NSPIN,LMMAX,IELAST,FILELEFT,FILERIGHT,
     &        KREL,NATYPD,LMMAXD,NEMBD1,KORBIT)       ! ruess: pass KORBIT to decimaread for newsolver
C     
C --> get the left and right host Delta_t matrices
C
         CFCTOR = ALAT/(8.D0*ATAN(1.0D0)) ! = ALAT/(2*PI)
         NSPINSO = NSPIN
         IF (OPT('NEWSOSOL')) NSPINSO = 1 ! ruess: only combined l-s index for newsolver
         DO ISPIN = 1,NSPINSO
            DO IE = 1,IELAST
               CALL DECIMAREAD(EZ,TK,NPNT1,NPNT2,NPNT3,NPOL,ISPIN,
     &                         LEFTTINVLL(1,1,1,ISPIN,IE),
     &                         RIGHTTINVLL(1,1,1,ISPIN,IE),VACFLAG,
     &                         IE,NLBASIS,NRBASIS,NAEZ,KAOEZ,KMROT,INS,
     &                         NSPIN,LMMAX,IELAST,FILELEFT,FILERIGHT,
     &                         KREL,NATYPD,LMMAXD,NEMBD1,KORBIT)
C
C --> host matrices have been written out in true units
C     they are used in p.u. units (see kloopz) --> convert them here
C
               CALL ZSCAL(LMMAXD*LMMAXD*NEMBD1,CFCTOR,
     &                          LEFTTINVLL(1,1,1,ISPIN,IE),1)
               CALL ZSCAL(LMMAXD*LMMAXD*NEMBD1,CFCTOR,
     &                          RIGHTTINVLL(1,1,1,ISPIN,IE),1)
            END DO
         END DO
C
C --> get the left and right host charge moments
C     ( not needed in single-iteration mode calculations )
C
cfivos        IF ( SCFSTEPS.GT.1 ) 
         CALL CMOMSREAD(NLBASIS,NRBASIS,NAEZ,
     &                                       CMOMHOST,VACFLAG,KAOEZ,
     &                                       NATYPD,NEMBD1,LMPOTD)
         CLOSE (37)
         CLOSE (38)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ELSE 
C :::::::::::::::::::::::::::::::::::::::::::::::: decipots (calc tmats)
         IERROR = 0
         CALL IOINPUT('DECIPOTS        ',UIO,IL,7,IERROR)
         IF ( IERROR.NE.0 ) THEN
            WRITE(6,99010)
            STOP
         END IF
         READ (UNIT=UIO,FMT='(A40)') FILELEFT
         CALL IOINPUT('DECIPOTS        ',UIO,IL+1,7,IERROR)
         READ (UNIT=UIO,FMT='(A40)') FILERIGHT
         CALL DECITSET(ALAT,BRAVAIS,EZ,IELAST,
     &                 NLBASIS,NRBASIS,FILELEFT,FILERIGHT,
     &                 INS,KVREL,KREL,NSPIN,KMROT,
     &                 VREF,RMTREF,NREF,REFPOT,RC,CREL,RREL,
     &                 LEFTTINVLL,RIGHTTINVLL,VACFLAG,
     &                 NEMBD1,IEMXD,IRMD,IPAND,
     &                 LMAXD,LMGF0D,LMMAXD,LM2D,NSPIND)
      END IF
C ======================================================================
C
99010 FORMAT (/,6X,'ERROR : Missing decimation files (t-mat or pot)',/,
     &        14X,'Please use one of the tokens DECIFILES/DECIPOTS',/)
      END
