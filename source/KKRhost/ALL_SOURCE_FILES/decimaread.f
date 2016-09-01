C*==decimaread.f    processed by SPAG 6.05Rc at 19:11 on  5 Jun 2004
      SUBROUTINE DECIMAREAD(EZ,TK,NPTP1,NPTP2,NPTP3,NPOL,ISPIN,
     &                      LEFTTINVLL,RIGHTTINVLL,VACFLAG,IENERGY,
     &                      NLBASIS,NRBASIS,NAEZ,KAOEZ,KMROT,INS,NSPIN,
     &                      LMMAX,IELAST,FILELEFT,FILERIGHT,
     &                      KREL,NATYPD,LMMAXD,NEMBD1)!,KORBIT)
C **********************************************************************
C *                                                                    *
C * This subroutine reads in the t-matrices of the left                *
C * and right host for the decimation method.                          *
C *                                                                    *
C * The t-matrices are writen out in kloopz1  (option 'deci-out')      *
C *                                                                    *
C * The host files contain the CMOMHOST data neeed to create the       *
C * interatomic potential in subroutine < vinterface >                 *
C * This is going to be read in < cmomsread >, the decimation files    *
C * are for this reason not rewinded.                                  *
C *                                                                    *
C * In case of 'vacuum' setting on one of the sides, no energy points  *
C * are read in and the VACLAG is set to TRUE.                         *
C *                                                                    *
C * A call of this routine with IENERGY = 0 means that we are not in   *
C * the energy loop - only the header of each decimation file read in  *
C *                                                                    *
C * IENERGY <> 0 reads in the matrices/energy at IENERGY -> returned   *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
      INTEGER KREL,NATYPD,NEMBD1,LMMAXD!,KORBIT
C     ..
C     .. Scalar arguments
      INTEGER NPTP1,NPTP2,NPTP3,NPOL,ISPIN,IENERGY
      INTEGER NLBASIS,NRBASIS,NAEZ,KMROT,INS,NSPIN,LMMAX,IELAST
      DOUBLE PRECISION TK
      CHARACTER*40 FILELEFT,FILERIGHT  
C     ..
C     .. Array arguments
      INTEGER KAOEZ(NATYPD,*)
      DOUBLE COMPLEX EZ(*)
      DOUBLE COMPLEX LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1),
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1)
      LOGICAL VACFLAG(2)
C     ..
C     .. Local variables
      INTEGER IOS,NPT1L,NPT2L,NPT3L,NPOLL,INSL
      INTEGER NSPINL,NAEZL,LMMAXL,KRELL,KMROTL
      INTEGER IEL,IH,LM1,LM2,I,IDUM,IH1,IHOST,ILHOST,NATHOST
      DOUBLE PRECISION ALATL,TEMPL,E2L,QMT,QMP
      DOUBLE COMPLEX EL
      DOUBLE COMPLEX W1(LMMAXD,LMMAXD)
      CHARACTER*80 BANER1
      CHARACTER*40 FILEHOST
      CHARACTER*5 CHHOST(2),STR5
      DOUBLE PRECISION BRAVAISL(3,3),RBASISL(3,NEMBD1)
C     ..
C     .. External Functions
      INTEGER LNGSTRING
      EXTERNAL LNGSTRING
C     ..
C     .. External Subroutines
      EXTERNAL ZCOPY
C     ..
C     .. Data statements
      DATA CHHOST/'LEFT ','RIGHT'/
C     ..
C ========================================================== IENERGY = 0
      IF ( IENERGY.LT.1 ) THEN
         WRITE (1337,'(5X,A,/,8X,65(1H-))')
     &           'Reading in host Delta_t matrices'
         VACFLAG(1) = .FALSE.
         VACFLAG(2) = .FALSE.
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
         DO IHOST = 1,2
            FILEHOST = FILELEFT
            NATHOST = NLBASIS
            IF ( IHOST.EQ.2 ) THEN
               FILEHOST = FILERIGHT
               NATHOST = NRBASIS
            END IF
            ILHOST = LNGSTRING(FILEHOST,40)
C
            WRITE (1337,'(8X,A5," side host: ",$)') CHHOST(IHOST)
C ----------------------------------------------------------------------
            IF ( FILEHOST(1:7).EQ.'vacuum' ) THEN
               WRITE (1337,'(A)') 'VACUUM will be used'
               VACFLAG(IHOST) = .TRUE.
C ----------------------------------------------------------------------
            ELSE
               WRITE (1337,'(A,/)') FILEHOST(1:ILHOST)
               OPEN (36+IHOST,FILE=FILEHOST,STATUS='OLD',IOSTAT=IOS)
C ......................................................................
               IF ( IOS.GT.0 ) THEN
                  WRITE (6,'(/,5X,2A)') 'ERROR: Can not open host file '
     &                                  ,FILEHOST(1:ILHOST)
                  WRITE (6,'(12X,A,A)') 
     &                           'vacuum    entry should be used to set'
     &                           ,' a vacuum-host on this side'
                  STOP '       < DECIMAREAD > '
               END IF
C ......................................................................
               DO IH = 1,3
                  READ (36+IHOST,99007) BANER1
               END DO
               READ (36+IHOST,99009) ALATL,NSPINL,NAEZL,LMMAXL,INSL,
     &                               KRELL,KMROTL
C ......................................................................
               IF ( (KRELL.NE.KREL) .OR. (KMROTL.NE.KMROT) .OR. 
     &              (INSL.NE.INS) .OR. (NSPINL.NE.NSPIN) .OR. 
     &              (LMMAXL.NE.LMMAX) .OR. (NAEZL.NE.NATHOST) ) THEN
                  WRITE (6,'(/,5X,2A)') 'ERROR: ',
     &                       'host not compatible with your input/sytem'
                  WRITE (6,'(14X,6(A6),/,8X,42(1H-))') '  KREL',
     &                   ' KMROT','   INS',' NSPIN',' LMMAX',' BASIS'
                  WRITE (6,'(8X,A6,6I6)') 'syst: ',KREL,KMROT,INS,NSPIN,
     &                   LMMAX,NATHOST
                  WRITE (6,'(8X,A6,6I6,/)') 'host: ',KRELL,KMROTL,INSL,
     &                   NSPINL,LMMAXL,NAEZL
                  STOP '       < DECIMAREAD > '
               END IF
C ......................................................................
               WRITE (1337,99001) ALATL,NSPINL,NAEZL,LMMAXL,INSL,KRELL,
     &                         KMROTL
C
               READ (36+IHOST,99007) BANER1
               READ (36+IHOST,99003) BRAVAISL
               WRITE (1337,99002) BRAVAISL
               READ (36+IHOST,99007) BANER1
               IH = LNGSTRING(BANER1,80)
               WRITE (1337,99008) BANER1(1:IH)
C ......................................................................
               IF ( KREL.EQ.0 ) THEN
                  DO IH = 1,NAEZL
                     READ (36+IHOST,99003) (RBASISL(I,IH),I=1,3)
                     WRITE (1337,99003) (RBASISL(I,IH),I=1,3)
                  END DO
               ELSE
                  DO IH = 1,NAEZL
                     READ (36+IHOST,99014) (RBASISL(I,IH),I=1,3),QMT,QMP
                     WRITE (1337,99004) (RBASISL(I,IH),I=1,3),QMT,QMP
                  END DO
               END IF
C ......................................................................
               READ (36+IHOST,99010) E2L,TEMPL
               READ (36+IHOST,99011) NPT1L,NPT2L,NPT3L,NPOLL
C ......................................................................
               IF ( (NPT1L.NE.NPTP1) .OR. (NPT2L.NE.NPTP2) .OR. 
     &              (NPT3L.NE.NPTP3) .OR. (NPOLL.NE.NPOL) .OR. 
     &              (ABS(TEMPL-TK).GT.1.D-4) ) THEN
                  WRITE (6,'(/,5X,2A)') 'ERROR: Host energy mesh',
     &                                 ' not compatible with your input'
                  WRITE (6,'(14X,5(A6),/,8X,40(1H-))') '   NE1',
     &                   '   NE2','   NE3','  NPOL','  TEMP'
                  WRITE (6,'(8X,A6,4I6,2X,F10.6)') 'input:',NPTP1,NPTP2,
     &                   NPTP3,NPOL,TK
                  WRITE (6,'(8X,A6,4I6,2X,F10.6,/)') ' host:',NPT1L,
     &                   NPT2L,NPT3L,NPOLL,TEMPL
                  STOP '       < DECIMAREAD > '
               END IF
C ......................................................................
               WRITE (1337,99005) E2L,TEMPL
               WRITE (1337,99006) NPT1L,NPT2L,NPT3L,NPOLL
            END IF
C ----------------------------------------------------------------------
            WRITE (1337,'(8X,65(1H-))')
         END DO
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
C                            Headers read in
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         RETURN
      END IF
C ========================================================== IENERGY = 0
C
C
C ======================================================================
C                  READ IN Delta_t matrices for IENERGY > 0
C ======================================================================
      IF ( (ISPIN.NE.2) .OR. (NSPIN.NE.1) ) THEN
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
         DO IHOST = 1,2
            NATHOST = NLBASIS
            IF ( IHOST.EQ.2 ) NATHOST = NRBASIS
C ----------------------------------------------------------------------
C --> read in Delta_t if it is not a vacuum host
C     attention: new format Dec. 2004
C
            IF ( .NOT.VACFLAG(IHOST) ) THEN
C ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ HOST ATOM LOOP
               DO IH = 1,NATHOST
                  CALL CINIT(LMMAXD*LMMAXD,W1)
C
 9910             CONTINUE
                  READ (36+IHOST,'(A5)',END=300) STR5
                  IF ( STR5.NE.'*****' ) GOTO 9910
C
                  READ (36+IHOST,99012) IEL,EL,IDUM
C ......................................................................
                  IF ( (ABS(EL-EZ(IENERGY)).GT.1.D-4) .OR. 
     &                 (IEL.NE.IENERGY) ) THEN
                     WRITE (6,'(/,5X,2A)') 'ERROR: Host energy mesh',
     &                                 ' not compatible with your input'
                     WRITE (6,'(14X,2(A6),/,8X,43(1H-))') '    IE',
     &                                 ' E(IE)'
                     WRITE (6,'(8X,A6,I6,1X,2(1X,F10.6))') 'input:',
     &                    IENERGY,EZ(IENERGY)
                     WRITE (6,'(8X,A6,I6,1X,2(1X,F10.6),/)') ' host:',
     &                    IEL,EL
                     STOP '       < DECIMAREAD > '
                  END IF
C ......................................................................
                  IF ( IDUM.NE.IH ) THEN
                     WRITE (6,'(/,5X,2A,/)') 'ERROR reading host file',
     &                      ' basis indexing wrong'
                     STOP '       < DECIMAREAD > '
                  END IF
C ......................................................................
C
C --> map the t-matrix corectly
C
                  IH1 = KAOEZ(1,NAEZ+(IHOST-1)*NLBASIS+IH)
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 9920             CONTINUE
                  READ(36+IHOST,99013) LM1,LM2,EL
                  IF ( (LM1+LM2).NE.0 ) THEN
                     W1(LM1,LM2) = EL
                     IF ( (LM1+LM2).LT.2*LMMAXD ) GOTO 9920
                  END IF
                  
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
                  IF ( IHOST.EQ.1 ) THEN
                     DO LM1 = 1,LMMAXD  !(KREL+KORBIT+1)*LMMAX
!                         CALL ZCOPY((KREL+KORBIT+1)*LMMAX,W1(1,LM1),1,
                        CALL ZCOPY(LMMAXD,W1(1,LM1),1,
     &                             LEFTTINVLL(1,LM1,IH1),1)
                     END DO
                  ELSE
                     DO LM1 = 1,LMMAXD  !(KREL+KORBIT+1)*LMMAX
!                         CALL ZCOPY((KREL+KORBIT+1)*LMMAX,W1(1,LM1),1,
                        CALL ZCOPY(LMMAXD,W1(1,LM1),1,
     &                             RIGHTTINVLL(1,LM1,IH1),1)
                     END DO
                  END IF
               END DO
C ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            END IF              ! not vacuum
C ----------------------------------------------------------------------
         END DO
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      END IF                    ! ispin.ne.2 .or. .nspin.ne.1
C ======================================================================
C
      IF ( IENERGY.EQ.IELAST ) WRITE (1337,*)
      RETURN
 300  CONTINUE
      STOP '        Error reading hostfile'
C
99001 FORMAT (10X,'ALAT=',F9.6,' NSPIN=',I2,'  NAEZ=',I3,' LMMAX=',I3,
     &        ' INS=',I1,' KREL=',I1,' KMROT=',I1)
99002 FORMAT (10X,'BRAVAIS '/10X,3F8.4/10X,3F8.4/10X,3F8.4)
99003 FORMAT (10X,3F8.4)
99004 FORMAT (10X,3F8.4,2F9.4)
99005 FORMAT (10X,'EF=',F10.6,' TEMP=',F10.4,' Kelvin')
99006 FORMAT (10X,'N1=',I3,' N2=',I3,' N3=',I3,' NPOL=',I3)
C
99007 FORMAT (A80)
99008 FORMAT (10X,A)
99009 FORMAT (5X,F9.6,7X,I2,7X,I3,7X,I3,5X,I1,6X,I1,7X,I1)
99010 FORMAT (3X,F10.6,6X,F10.4)
99011 FORMAT (3X,I3,4X,I3,4X,I3,6X,I3)
99012 FORMAT (7X,I5,2D16.8,6X,I3)
99013 FORMAT (2I5,1P,2D22.14)
99014 FORMAT (3F8.4,2F9.4)
C
      END
