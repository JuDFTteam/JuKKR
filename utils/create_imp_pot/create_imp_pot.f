      PROGRAM CREATE_IMP_POT
      INTEGER NATYPD,NLINED
      PARAMETER (NATYPD=400,NLINED=7000)
C     
      INTEGER NSPIN,IAT,ILIN,IOS
      INTEGER, ALLOCATABLE:: LINESIMP(:),LINES(:)
      INTEGER,ALLOCATABLE :: ORDER(:),ORDERBLK(:),ORDERBLK2(:),IMAP(:)
      CHARACTER*90 TEXT 
      CHARACTER*90 TXT(:,:)
      CHARACTER*90  TXTIMP(:,:)
      CHARACTER*50 POT
      CHARACTER SP
      LOGICAL GO_ON,READORDER,MAP
      ALLOCATABLE TXTIMP,TXT
                                ! READ(5,*) POT
                                ! READ(5,*) 
      READORDER=.FALSE.
      NSPIN = 2
      CALL GETARG(1,POT)
      WRITE(6,*) 'FILENAME OF ORIGINAL POTENTIAL (MAX 50 CHAR)'
      WRITE(6,*) POT
      OPEN(10,FILE=POT,FORM='FORMATTED',IOSTAT=IOS)
      IF (IOS.NE.0) THEN
          WRITE(6,*) ' I CANNOT FIND FILE >',POT
          STOP
      END IF
      CALL GETARG(2,POT)
      MAP=.FALSE.

      IF(POT(1:4).EQ.'-map') THEN
          MAP=.TRUE.
          CALL GETARG(3,POT)
          OPEN(UNIT=60,FILE=POT,IOSTAT=IOS2)
          IAT = 0
          ILIN = 0
          IOS2= 0
          NLINIMP=0
          DO WHILE (IOS2.EQ.0)
              READ(60,9000,END=99) TEXT
              IF (TEXT(37:40).EQ.'exc:') THEN
                  IAT = IAT + 1
                  ILIN = 1 
              ELSE
                  ILIN = ILIN + 1
              END IF
              NLINIMP=MAX(NLINIMP,ILIN)
          END DO
 99       WRITE(6,*) ' **** END OF THE MAPED FILE ****** '
          INATIMP=IAT/2
          ALLOCATE(LINESIMP(0:2*IAT))
          ALLOCATE (TXTIMP(2*IAT,NLINIMP))
          WRITE(6,*) "FOUND: IATOMS,LINES",INATIMP,NLINIMP
          REWIND(60)
          IAT = 0
          ILIN = 0
          IOS2= 0
          DO WHILE (IOS2.EQ.0)
              READ(60,9000,END=98) TEXT
              IF (TEXT(37:40).EQ.'exc:') THEN
                  IAT = IAT + 1
                  ILIN = 1 
              ELSE
                  ILIN = ILIN + 1
                  LINESIMP(IAT) = ILIN
              END IF
              TXTIMP(IAT,ILIN) = TEXT
          END DO
 98       WRITE(6,*) ' **** END OF FILE '
      END IF

      OPEN(9,FILE='composed.potent',FORM='FORMATTED',IOSTAT=IOS)

      IAT = 0
      ILIN = 0
      IOS= 0
      NLINHOST=0
      DO WHILE (IOS.EQ.0)
          READ(10,9000,END=97) TEXT
          IF (TEXT(37:40).EQ.'exc:') THEN
              IAT = IAT + 1
              ILIN = 1 
          ELSE      
              ILIN = ILIN + 1
          END IF
          NLINHOST=MAX(NLINHOST,ILIN)
      END DO 
      
 97   WRITE(6,*) ' **** END OF HOST FILE *************'

      NSPIN1 = 2
      IF (MOD(IAT,2).EQ.0) THEN
          NSPIN1=2  
          IAT = IAT/NSPIN1
      END IF
      ALLOCATE(LINES(0:2*IAT))
      ALLOCATE (TXT(2*IAT,NLINHOST))
      INATHOST=IAT
      WRITE(6,*) ' NO OF ATOMS: ',IAT,' NO OF SPINS : ',NSPIN1
      REWIND(10)

      IAT = 0
      ILIN = 0
      IOS= 0
      DO WHILE (IOS.EQ.0)
          READ(10,9000,END=100) TEXT
          IF (TEXT(37:40).EQ.'exc:') THEN
              IAT = IAT + 1
              ILIN = 1 
          ELSE      
              ILIN = ILIN + 1
              LINES(IAT) = ILIN
          END IF
          TXT(IAT,ILIN) = TEXT   
      END DO 
 100  CONTINUE
      CLOSE(10)

      WRITE(6,*) ' GIVE ORDER OF ATOMS TO BE WRITTEN OUT '
      I = 1
      INQUIRE (FILE='impurity.coefs',EXIST=READORDER)
      IF(READORDER) WRITE(6,*) 'ORDER READ IN FROM IMPURITY.COEFS'
      IF(.NOT. READORDER) STOP "MISSING IMPURITY.COEFS"
      IF(READORDER) THEN
C     NOW READ IN THE IMPURITY.COEFS FILE FOR THE IMPURITY CALCULATION
          OPEN(58,FILE='impurity.coefs',FORM='FORMATTED') 
          READ(58,8001) IDUM,NATOMIMP,IDUM,IDUM,
     &         (IDUM,I=1,NATOMIMP)
          READ(58,*)
          ALLOCATE(ORDERBLK2(NATOMIMP),ORDERBLK(NATOMIMP))
          DO I=1,NATOMIMP
              READ(58,8002) 
     &             RDUM,RDUM,RDUM,ORDERBLK2(I),IDUM,IDUM,
     &             RDUM,ORDERBLK(I)
          END DO
          READ(58,8001) IDUM,IDUM
          READ(58,8001) IDUM
          READ(58,8005) NHOST
          ALLOCATE(ORDER(NHOST+NATOMIMP))
          READ(58,8006) (ORDER(I),I=1,NHOST)
          WRITE(*,9999) 
          DO I=1,NATOMIMP
              ORDER(I+NHOST)=ORDERBLK(I)
          END DO
          N=NATOMIMP+NHOST
          NSPIN2=NSPIN1
          ALLOCATE (IMAP(N))
          IMAP=0
          WRITE(*,9999) 
          CLOSE(58)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     WRITE OUT INPUT FILE FOR IMPURITY PROGRAM     
C     CHECK ALL PARAMETERS YOU NEED. SOME STANDARD (MY JM)
C     SETTING ARE USED HERE (FULLY RELAT. ..ETC)
C     
          OPEN(UNIT=58,STATUS='SCRATCH')
          WRITE(58,*) TXT(1,5)
          REWIND(58)
          READ(58,*) IRM
          CLOSE (58)
          OPEN(UNIT=58,STATUS='SCRATCH')
          WRITE(58,*) TXT(1,4)
          REWIND(58)
          READ(58,*) RDUM,RDUM,VBCC
          CLOSE (58)

 9019     FORMAT (A20)
 9120     FORMAT (8I5)
 9121     FORMAT (8I5,A50)
 9122     FORMAT (4I5,A50)
 9123     FORMAT (7I5,A50)
 9124     FORMAT (3I5,A50)
 9125     FORMAT (3F10.6,A40)
 9126     FORMAT (5I5,A40)
          WRITE(*,*) "WRITING IMPUT FILE INFORMATION INTO INPUT.TMP"
          OPEN(UNIT=58,FILE="INPUT.TMP")
C     NAME OF GREENS FUNC FILE HAS TO BE GREEN
          WRITE(58,9019) 'GREEN               '  
C     POTENTIAL NAME (YOU SHOULD ADD POTENTIAL OF YOUR IMPURITY)
          WRITE(58,9019) 'COMPOSED.POTENT     '  
C     !IMPORTANT CREATED FROM TB-KKR
          WRITE(58,9019) 'IMPURITY.COEFS      ' 
C     !SHAPE FUNCTIONS FOR FP
          WRITE(58,9019) 'SHAPEFUN            ' 

          WRITE(58,9122) 0,12,0,0
     &         ,'                         IRESIST,IFILE,IPE,IWTMAT' 
          WRITE(58,9122) NSPIN2,IRM,0,0
     &         ,'                         NSPIN,IRM,INS,ICST'
          WRITE(58,9123) 3,3,2,0,0,0,2,
     &         '        KCOR,KVREL,KWS,KHYP,KHFELD,KFSR,KXC'
          WRITE(58,9123) 0,0,0,0,0,0,0,
     &         '        KTE,KPRE,KSPH,KEFG,KVMAD,KF,IGGA'
          WRITE(58,9122) NHOST,NATOMIMP,0,0
     &         ,'                    NATREF,NATPER,ICUT,KSHAPE' 
          WRITE(58,9120) (ORDERBLK2(I),I=1,NATOMIMP)
          WRITE(58,9120) (1,I=1,N)
          WRITE(58,9124) 100,4,1
     &         ,'ITCLST,IMIX,IEF(FOR REAL E)'
          WRITE(58,9125) 0.001,1.000,0.000001
     &         ,'STRMIXING,FCM,QBOUND'
          WRITE(58,9125) 0.0010,0.0,0.0
     &         ,'BRYMIX,EF(1),EF(2)'
          WRITE(58,9125) 0.0,VBCC,0.00
     &         ,'HFIELD,VBCC,VCONST'
          WRITE(58,9126) 0,0,0,0,0
     &         ,'KLATR,KSYMM,LKONV,KSYMMAT,KESYM'
          CLOSE(58)
          WRITE(*,9999) 
      END IF
      WRITE(6,9999)
      WRITE(6,*) ' ATOM COPIED FROM ATOM '
C     
C     NOW WRITE OUT
C     
C     
      IF(MAP) THEN
          WRITE(*, 9999)
          WRITE(*,*) " YOU USED OPTION -map:"
          WRITE(*,*)
     &         "WHICH ATOMS SHOULD BE READ IN FROM SECOND POTENTIAL?"
          WRITE(*,*) "FORMAT: IT FORM IMPURITY_COEFS, IT FROM SECNDPOT"
          WRITE(*,*) "IT - 0: STOP"
          DO I1=1,INATIMP
              READ(*,*) IT,IIMP
              WRITE(*,*) IT,IIMP
              IF(IIMP.GT.INATIMP) STOP
     &             "IT FROM SECOND POT .GT. NUMBER OF ATOMS READ IN"
              IF(IT.GT.NATOMIMP) STOP
     &             "IT FORM IMPURITY_COEFS .GT. NUMBER OF ATOMS READ IN"
              IF(IT.EQ.0) GOTO 88
              IMAP(IT)=IIMP
          END DO
 88       CONTINUE
      END IF

      IF (NSPIN1.EQ.1) THEN
          STOP "NSPIN EQ 1 NOT YET IMPLEMENTED"
          IF (NSPIN2.EQ.1) THEN
              DO I=1,N
                  ID = ORDER(I)
                  DO I1=1,LINES(ID)
                      WRITE(9,9000) TXT(ID,I1)
                  END DO  
              END DO 
          ELSE
              DO I=1,N
                  ID = ORDER(I)
                  DO I1=1,LINES(ID)
                      WRITE(9,9000) TXT(ID,I1)
                  END DO
                  DO I1=1,LINES(ID)
                      WRITE(9,9000) TXT(ID,I1)
                  END DO  
              END DO 
              
          END IF
      ELSE 

          IF (NSPIN2.EQ.1) THEN
              DO I=1,N
                  ID = 2*ORDER(I)-1
                  DO I1=1,LINES(ID)
                      WRITE(9,9000) TXT(ID,I1)
                  END DO
              END DO
          ELSE
C     FIRST WRITE OUT HOST 
              DO I=1,NHOST
                  ID = 2*ORDER(I)-1
                  DO I1=1,LINES(ID)
                      WRITE(9,9000) TXT(ID,I1)
                  END DO
                  DO I1=1,LINES(ID+1)
                      WRITE(9,9000) TXT(ID+1,I1)
                  END DO
              END DO  
C     WRITE OUT IMP  
              DO I=NHOST+1,N
                  IF(IMAP(I-NHOST).EQ.0) THEN
                      ID = 2*ORDER(I)-1
                      DO I1=1,LINES(ID)
                          WRITE(9,9000) TXT(ID,I1)
                      END DO
                      DO I1=1,LINES(ID+1)
                          WRITE(9,9000) TXT(ID+1,I1)
                      END DO
                  ELSE
                      ID = 2*IMAP(I-NHOST)-1
                      DO I1=1,LINESIMP(ID)
                          WRITE(9,9000) TXTIMP(ID,I1)
                      END DO
                      DO I1=1,LINESIMP(ID+1)
                          WRITE(9,9000) TXTIMP(ID+1,I1)
                      END DO
                  END IF
              END DO  
          END IF
          
          
      END IF 
      write(6,*) ' Potential written in file : composed.potent'         
      WRITE(*,*) " IMPUT FILE INFORMATION INTO input.tmp"
      close(11)         
 9000 format(a80)
 8001 FORMAT(11I5)
 8002 FORMAT (3F12.8,I4,I4,I5,F10.6,I5)
 8003 FORMAT (5A10)
 8004 FORMAT ('Shells in the reduced format')
 8005 FORMAT ('Host order, no of host sites: ',I5)
 8006 FORMAT (12I4)
 8007 FORMAT (I5,'    Symmetries for the Bulk')
 8010 FORMAT ('     Position of Impurity            Host Imp Shell',
     &     '   Dist     Host id in Bulk')
 8030 FORMAT(3I5,3F12.8,I8,F10.5)
 8040 FORMAT(' >>> Some host atoms are missing in the impurity <<<'/
     &     '     INDEXING WILL BE CHANGED:  PLEASE CHECK ALSO  ') 
 8050 FORMAT(' ******************************************'/
     &     ' *  Prepearing indexing for impurity GF   *'/
     &     ' *  Unsymmetrized GF is written out in    *'/
     &     ' *  this version                          *'/
     &     ' *  Producing:   green,                   *'/
     &     ' *               impurity.coefs,          *'/
     &     ' *               intercell_ref            *'/
     &     ' *                            v. 20.09.01 *'/
     &     ' ******************************************')
 9030 FORMAT(I3,I7,I5,F14.3,2F11.3,I8,F8.3)
 9999 FORMAT(/,79('*'),/)
      END

      
      









