      SUBROUTINE GFSHELLS(ICC,NATOMIMP,NSH1,NSH2,
     &                    IJTABSYM,IJTABSH,IJTABCALC,
     &                    IOFGIJ,JOFGIJ,NOFGIJ,ISH,JSH,
     &                    NSHELL,NAEZ,NATYP,NOQ,RBASIS,BRAVAIS,
     &                    IFILIMP,RATOM,RCLSIMP,
     &                    NSYMAT,ISYMINDEX,RSYMAT,
     &                    KAOEZ,ATOMIMP,
     &                    ROTNAME,IPRINT,HOSTIMP,LMAXD,LMMAXD,
     &                    NAEZD,NATYPD,NATOMIMPD,NEMBD,NSHELD)
C **********************************************************************
C *                                                                    *
C * This subroutine constructs mainly the index arrays                 *
C * NSHELL, NSH1, NSH2 -- NSHELL(0) number of different GF blocks that *
C * have to be calculated, NSH1(I),NSH2(I) the sites connected for     *
C * the block I, I = 1,NSHELL(0)                                       *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
      INTEGER  LMAXD,LMMAXD,NAEZD,NATYPD,NATOMIMPD,NEMBD,NSHELD
C     ..
C     .. Scalar arguments
      INTEGER ICC,NAEZ,NATOMIMP,NATYP,NSYMAT,IPRINT,NOFGIJ
      CHARACTER*40 IFILIMP
C     ..
C     .. Array arguments
      CHARACTER*10 ROTNAME(64)
C     ..
      INTEGER ATOMIMP(NATOMIMPD),HOSTIMP(0:NATYPD)
      INTEGER ISYMINDEX(*),KAOEZ(NATYPD,NAEZD+NEMBD)
      INTEGER NOQ(NAEZD),NSH1(*),NSH2(*),NSHELL(0:NSHELD)
      INTEGER ISH(NSHELD,*),JSH(NSHELD,*)
      INTEGER IJTABSYM(*),IJTABSH(*),IJTABCALC(*),IOFGIJ(*),JOFGIJ(*)
C     ..
      DOUBLE PRECISION BRAVAIS(3,3),RATOM(3,NSHELD)
      DOUBLE PRECISION RBASIS(3,*),RCLSIMP(3,NATOMIMPD)
      DOUBLE PRECISION RSYMAT(64,3,3)
C     .. 
C     .. Local scalars
      INTEGER NB,I,J,POS,II,IO,NS,IN,NDIM,NSIZE
      CHARACTER*9 STR9
      LOGICAL LSURF,OPT
C     ..
C     .. External subroutines
      EXTERNAL IMPCHECK,IMPCOEFS,SHELLGEN2K,OPT
C
      WRITE (1337,99000)
      NSIZE = NATOMIMPD*LMMAXD
C
C **********************************************************************
C     
C --> construction of ratom, nsh1 and nsh2 for a self-consistent
C     calculation
C     
      IF ( .not. OPT('VIRATOMS') ) THEN
       NSHELL(0) = NATYP
      ELSE
       NSHELL(0) = NAEZ
      END IF !( .not. OPT('VIRATOMS') ) THEN

       IF ( NSHELL(0).GT.NSHELD ) THEN
          WRITE(6,99001) 'NSHELD',NSHELL(0)
          STOP 
       END IF
      
       DO I=1,NSHELL(0)
          RATOM(1,I) = 0.0D0
          RATOM(2,I) = 0.0D0
          RATOM(3,I) = 0.0D0
          NSHELL(I) = 0
          
          DO J=1,NAEZ
             DO IO=1,NOQ(J)
                IF (KAOEZ(IO,J).EQ.I) THEN
                   NSHELL(I) = NSHELL(I) + 1
                   IF (NSHELL(I).EQ.1) THEN
                      NSH1(I) = J
                      NSH2(I) = J                    
                   END IF
                END IF
             END DO
          END DO
          IF ( OPT('VIRATOMS') ) THEN
            NSHELL(I)=1
            NSH1(I) = I
            NSH2(I) = I                    
          END IF


C          
          IF (NSHELL(I).EQ.0) THEN
             WRITE(6,99002)
             STOP
          END IF
       END DO
C
       IF ( ICC.EQ.0 ) THEN
          WRITE(1337,99003) NSHELL(0)
          RETURN
       END IF
C
C      end of simple SCF-calculation part. 
C **********************************************************************
C
Check if we are in surface mode
C
      LSURF = .FALSE.
      IF ( BRAVAIS(1,3).EQ.0D0 .AND. BRAVAIS(2,3).EQ.0D0 .AND.
     &     BRAVAIS(3,3).EQ.0D0 ) LSURF = .TRUE. 
      NDIM = 3
      IF (LSURF) NDIM = 2
C     
C **********************************************************************
C      NATOMIMP=0   ! BUG: This initialization breaks the shell generation for
C                   ! ICC=-1, which is set by option XCPL.  B. Zimmermann

       IF (ICC.LT.0) THEN
C
C --->  ICC.LT.1 all shells are (should be) prepared
C
          WRITE(1337,99011) NATOMIMP
       ELSE
C
C --> read-in the cluster coordinates from an external file
C
          REWIND 25
          READ (25,FMT=*) NATOMIMP
C
          IF (NATOMIMP.GT.NATOMIMPD ) THEN
             WRITE(6,99001) 'NATOMIMPD',NATOMIMP
             STOP
          END IF
          WRITE(1337,99004) IFILIMP,NATOMIMP
C
          DO I=1,NATOMIMP
             READ (25,FMT=*) (RCLSIMP(J,I),J=1,3),ATOMIMP(I)
             ATOMIMP(I) = ATOMIMP(I) + ICC - 1
          ENDDO
       END IF
C **********************************************************************
C          
       CALL IMPCHECK(ATOMIMP,NATOMIMP,NAEZ,RCLSIMP,RBASIS,BRAVAIS,NDIM)
C
C **********************************************************************
       IF ( ICC.GT.0 ) THEN
          WRITE(1337,99005) 
C
C --> set up the number of all (I,J)-pairs to be looked for, 
C     avoid considering again the diagonal elements
C
          NOFGIJ = 0 
          DO I = 1,NATOMIMP
             NB = (I-1)*NATOMIMP
             DO J = 1,NATOMIMP
                IJTABCALC(NB+J) = 0
             END DO
             IF ( ATOMIMP(I).GE.0 ) THEN
                DO J = 1,NATOMIMP
                   IF ( ( ATOMIMP(J).GE.0 ).AND.( I.NE.J ) ) THEN
                      NOFGIJ = NOFGIJ + 1
                      IF ( NOFGIJ.GT.NATOMIMPD*NATOMIMPD ) THEN
                         WRITE(6,99001) 'NATOMIMPD',NOFGIJ/NATOMIMP
                         STOP 
                      END IF
                      IOFGIJ(NOFGIJ) = I
                      JOFGIJ(NOFGIJ) = J
                      IJTABCALC(NB+J) = 1
                   END IF
                END DO
             END IF
          END DO
       END IF
C **********************************************************************
C
       CALL SHELLGEN2K(ICC,NATOMIMP,RCLSIMP(1,1),ATOMIMP(1),
     +                 NOFGIJ,IOFGIJ,JOFGIJ,
     +                 NSYMAT,RSYMAT,ISYMINDEX,ROTNAME,
     +                 NSHELL,RATOM(1,1),NSH1,NSH2,ISH,JSH,
     +                 IJTABSYM,IJTABSH,IJTABCALC,2,NSHELD)
C    +                 IJTABSYM,IJTABSH,IJTABCALC,IPRINT,NSHELD)
C
C **********************************************************************
C
C --> now write out the impurity.coefs file for the impurity calculation
C                                                         n.papanikolaou
C
       IF ( ICC.GT.0 .or. OPT('KKRFLEX ')) 
     &      CALL IMPCOEFS(NATOMIMP,NAEZ,ATOMIMP,RCLSIMP,NSHELL,
     &                    NSH1,NSH2,RATOM,NSYMAT,ISYMINDEX,ROTNAME,
     &                    HOSTIMP,NATYPD,LMAXD,NSHELD,NSIZE)
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
       WRITE(1337,99003) NSHELL(0)
       WRITE(1337,99006) 
       NB = MAX(NATYP,NAEZ)
       DO NS = 1,NSHELL(0)
          IF ( NS.EQ.NB+1 ) WRITE(1337,99012) 
          IF ( NS.LE.NB ) THEN
             CALL SETPAIRSTR(NSH1(NS),NSH2(NS),STR9)
             WRITE(1337,99007) 
     +            NS,NSH1(NS),NSH2(NS),(RATOM(II,NS),II=1,3),
     +            SQRT(RATOM(1,NS)**2+RATOM(2,NS)**2+RATOM(3,NS)**2),
     +            STR9
          ELSE
             WRITE(1337,99008) 
     +            NS,NSH1(NS),NSH2(NS),(RATOM(II,NS),II=1,3),
     +            SQRT(RATOM(1,NS)**2+RATOM(2,NS)**2+RATOM(3,NS)**2)
             IO = MIN(2,NSHELL(NS))
             DO I = 1,IO
                CALL SETPAIRSTR(ISH(NS,I),JSH(NS,I),STR9)
                WRITE(1337,'(A9,$)') STR9
             END DO
             WRITE(1337,*)
             POS = (NSHELL(NS)+1)/2
             DO I = 2,POS
                IO = (I-1)*2
                IN = MIN(2,NSHELL(NS)-IO)
                WRITE(1337,99009)
                DO J = 1,IN
                   CALL SETPAIRSTR(ISH(NS,IO+J),JSH(NS,IO+J),STR9)
                   WRITE(1337,'(A9,$)') STR9
                END DO
                WRITE(1337,*)
             END DO
          END IF
       END DO
       WRITE(1337,'(6X,72(1H-))') 
       NB = 0
       DO NS=1,NSHELL(0)
          NB = NB + NSHELL(NS)     
       END DO
       WRITE(1337,99010) NB
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
C ----------------------------------------------------------------------
99000 FORMAT(5X,"< GFSHELLS > : setting up indices of the GF blocks",/)
99001 FORMAT(6X,"Dimension ERROR: please increase the global parameter",
     &     /,6X,A," to a value >=",I5,/)
99002 FORMAT(6X,"ERROR: there are some inconsistencies in your input",/,
     &      13X,"not all atoms defined by NATYP have been found",/)
99003 FORMAT(8X,"Different shells for GF calculation : ",I3,/)
99004 FORMAT(8X,"Reading in cluster impurity sites from file",/,
     &      12X,"file name        : ",A,/,12X,"atoms in cluster : ",I3)
99005 FORMAT(8X,
     &   "Preparing indexing for impurity GF",/,11X,
     &   "- unsymmetrised GF is written out (v. 20.09.2001)",/,11X,
     &   "- files that will be created: impurity.coefs",/,41X,
     &   "intercell_ref",/,41X,"green",/)
99006 FORMAT(6X,72(1H-),/,6X,"shell|"," IQ ",
     &      " JQ"," | ",10X,"vec R_IJ ",11X,"R_IJ   | equiv. pairs",/,
     &      6X,72(1H-))
99007 FORMAT(5X,I5," |",I3,1X,I3," | ",3F9.4,F9.5,1X,"|",A9)
99008 FORMAT(5X,I5," |",I3,1X,I3," | ",3F9.4,F9.5,1X,"|",$)
99009 FORMAT(5X,5X," |",7X," | ",27X,9X,1X,"|",$)
99010 FORMAT(8X,"Number of block elements to be calculated : ",I3,/)
99011 FORMAT(8X,"Setting pairs for task-defined cluster sites ",
     &       "and connections",/,
     &      12X,"atoms in cluster : ",I3)
99012 FORMAT(6X,72(1H:),/,
     &     22X,"(impurity) cluster related data/indexing",/,
     &     6X,72(1H:))

       END
C **********************************************************************
      SUBROUTINE SETPAIRSTR(I,J,STR9)
      IMPLICIT NONE
      CHARACTER*9 STR9,STRD
      INTEGER I,J,L,LSTR
      CHARACTER*20 FMT1
C     ..
      FMT1 = '("(",I'
      FMT1 = FMT1(1:6)//'1'
      LSTR = 4
      IF ( I.GE.10 ) THEN
         FMT1 = FMT1(1:6)//'2'
         LSTR = LSTR + 1
         IF (I.GE.100 ) THEN
            FMT1 = FMT1(1:6)//'3'
            LSTR = LSTR + 1
         END IF
      END IF
      FMT1 = FMT1(1:7)//',",",I'
      FMT1 = FMT1(1:13)//'1'
      LSTR = LSTR + 1
      IF ( J.GE.10 ) THEN
         FMT1 = FMT1(1:13)//'2'
         LSTR = LSTR + 1
         IF ( J.GE.100 ) THEN
            FMT1 = FMT1(1:13)//'3'
            LSTR = LSTR + 1
         END IF
      END IF
      FMT1 = FMT1(1:14)//',")")'
      WRITE(STRD,FMT1) I,J
      DO L = 1,9-LSTR
         STR9(L:L) = ' '
      END DO
      STR9 = STR9(1:9-LSTR)//STRD(1:LSTR)
      END 
C **********************************************************************
