C*==vinterface.f    processed by SPAG 6.05Rc at 16:57 on  1 Feb 2002
      SUBROUTINE VINTERFACE(CMOM,CMINST,LPOT,NSPIN,NLAYERS,NATYP,V,ZAT,
     &                      R,IRWS,IRCUT,IPAN,KSHAPE,NOQ,KAOEZ,IQAT,
     &                      CONC,CATOM,ICC,HOSTIMP,
     &                      NLBASIS,NLEFT,NRBASIS,NRIGHT,
     &                      CMOMHOST,CHRGNT,VINTERS)
C **********************************************************************
C
C  This is calculating the intra-atomic contibution of the potential in 
C  the case of an interface taking into account the bulk potential on 
C  the two sides
C  it uses the structure dependent matrices AVMAD which are calculated
C  once in the subroutine MADELUNG2D and saved in the DA-file
C  avmad.unformatted                                         ( may 2004)
C
C  For each site in a layer the summation in all other layers is split
C  into three parts: within the slab, over the NLEFT*NLBASIS left host
C  sites and over the NRIGHT*NRBASIS right host sites, the last two
C  steps only in case of decimation run
C
C ----------------------------------------------------------------------
C
C  adopted for the case of more atoms on the same site, summation is
C  done over the occupants of that site, the charge is weighted with
C  the appropriate concentration of the occupant  v.popescu feb. 2002
C
C ----------------------------------------------------------------------
C     
C  impurity-program adopted feb. 2004 (according to n.papanikalou)
C     
C **********************************************************************
      IMPLICIT NONE
      INCLUDE 'inc.p'
C     ..
C     .. PARAMETER definitions
      INTEGER LMPOTD
      PARAMETER (LMPOTD=(LPOTD+1)**2)
C     ..
C     .. Scalar arguments ..
      INTEGER LPOT,NSPIN,NLAYERS,NATYP,KSHAPE,ICC
      INTEGER NLBASIS,NLEFT,NRBASIS,NRIGHT
      DOUBLE PRECISION CHRGNT
C     ..
C     .. Array arguments ..
      DOUBLE PRECISION CMOM(LMPOTD,*),CMINST(LMPOTD,*),V(IRMD,LMPOTD,*),
     &                 ZAT(*),R(IRMD,*),CONC(NATYPD),CATOM(NATYPD),
     &                 CMOMHOST(LMPOTD,*)
      INTEGER IRWS(*),IRCUT(0:IPAND,*),IPAN(*),
     &        NOQ(NAEZD),KAOEZ(NATYPD,NAEZD+NEMBD),IQAT(NATYPD),
     &        HOSTIMP(0:NATYPD)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION AC(LMPOTD),AVMAD(LMPOTD,LMPOTD),
     &                 CHARGE(2),CM(LMPOTD),MONOPOL(NAEZD),
     &                 VINTERS(LMPOTD,NAEZD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,CM1,FPI
      INTEGER I,IATOM,IB,IH1,ILAY1,ILAY2,IO2,
     &        IPOT,IRS1,ISPIN,IT1,IT2,L,LM,LM2,LMPOT,M
      LOGICAL OPT,TEST
      INTEGER LRECAMAD,IREC,NLEFTOFF,NRIGHTOFF,NLEFTALL,NRIGHTALL
      INTEGER ILEFT,IRIGHT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
C     .. External Functions/Subroutines
      EXTERNAL OPT,TEST
C     ..
      IF(TEST('flow    ')) WRITE (1337,*) '>>>>>> Vinterface'
C
      LRECAMAD = WLENGTH*2*LMPOTD*LMPOTD
C
      OPEN (69,ACCESS='direct',RECL=LRECAMAD,FILE='avmad.unformatted',
     +     FORM='unformatted')
C
      WRITE(1337,FMT=99001)
      WRITE(1337,FMT=99002)
C
      PI = 4.D0*ATAN(1.D0)
      FPI = 4.D0*PI
      LMPOT = (LPOT+1)**2
C
      IF ( OPT('DECIMATE') ) THEN
C
C     setup the charges to put in the ghost layers in the case of 
C     decimation technique to achieve charge neutrality
C
         CHARGE(1) = -CHRGNT/(2.D0*SQRT(FPI))
         CHARGE(2) = -CHRGNT/(2.D0*SQRT(FPI))
C
         NLEFTOFF = NLAYERS * NLAYERS                     ! record offsets
         NRIGHTOFF = NLEFTOFF + NLAYERS * NLEFT * NLBASIS ! left and right
         NLEFTALL = NLEFT * NLBASIS
         NRIGHTALL = NRIGHT * NRBASIS
      END IF
C **********************************************************************
C                   START CALCULATION IN THE LAYERS
C **********************************************************************
C ********************************************** loop over atoms in slab
C
      DO IT1 = 1,NATYP
C ========================================== take a site occupied by IT1
         ILAY1 = IQAT(IT1)
C
         IF ( KSHAPE.NE.0 ) THEN
            IRS1 = IRCUT(IPAN(IT1),IT1)
         ELSE
            IRS1 = IRWS(IT1)
         END IF
C
         DO LM = 1,LMPOT
            AC(LM) = 0.D0
         END DO
C
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                     1.  Summation in all layers in the slab
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         DO ILAY2 = 1,NLAYERS
            IREC = ILAY2 + NLAYERS*(ILAY1-1) 
            READ(69,REC=IREC) AVMAD

C
C --> Keep the monopole term -- Hoshino is doing 
C                              (SMONOPOL(I) -SMONOPOL(0))
C
            IF ( ILAY1.EQ.ILAY2 ) MONOPOL(ILAY1) = AVMAD(1,1)
C
C ================================ loop over all occupants of site ILAY2
            DO IO2 = 1,NOQ(ILAY2)
               IT2 = KAOEZ(IO2,ILAY2)
C
               DO LM = 1,LMPOT
                  CM(LM) = CMOM(LM,IT2)
C
C---> add contribution of interstial in case of shapes
C
                  IF ( KSHAPE.NE.0 ) CM(LM) = CM(LM) + CMINST(LM,IT2)
               END DO
               CM(1) = CM(1) - ZAT(IT2)/SQRT(FPI)
C
               DO LM = 1,LMPOT
                  DO LM2 = 1,LMPOT
                     AC(LM) = AC(LM) + AVMAD(LM,LM2)*CM(LM2)*CONC(IT2)
                  END DO
               END DO
            END DO
C ================================ loop over all occupants of site ILAY2
         END DO                 ! ILAY2 loop in all interface planes
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         DO ILAY2 = 1,NLAYERS
C ================================ loop over all occupants of site ILAY2
            DO IO2 = 1,NOQ(ILAY2)
               IT2 = KAOEZ(IO2,ILAY2)
C
               CM1 = CMOM(1,IT2)
               IF ( KSHAPE.NE.0 ) CM1 = CM1 + CMINST(1,IT2)
C
               CM1 = CM1 - ZAT(IT2)/SQRT(FPI)
               AC(1) = AC(1) - MONOPOL(ILAY1)*CM1*CONC(IT2)
C
            END DO           
C ================================ loop over all occupants of site ILAY2
         END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C --> Correction: charge neutrality is imposed (see P. Lang)
C
C ######################################################################
         IF ( OPT('DECIMATE') ) THEN
C 
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                     2.  Summation in the LEFT bulk side
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C +++++++++++++++++++++++++++++++++ loop over all occupants of LEFT host
            ILEFT = 0
            DO IH1 = 1,NLEFT
               DO IB = 1,NLBASIS
                  ILEFT = ILEFT + 1
                  IREC = ILEFT + NLEFTALL*(ILAY1-1) + NLEFTOFF
                  READ(69,REC=IREC) AVMAD
C
                  IATOM = IB
                  DO LM = 1,LMPOT
                     DO LM2 = 1,LMPOT
                        AC(LM) = AC(LM) 
     &                       + AVMAD(LM,LM2)*CMOMHOST(LM2,IATOM)
                     END DO
                  END DO
C
                  IF ( (IH1.EQ.1) .AND. (IB.EQ.1) ) AC(1) = AC(1) + 
     &                 (AVMAD(1,1)-MONOPOL(ILAY1))*CHARGE(1)
               END DO           
            END DO              
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF ( ILEFT.NE.NLEFTALL) THEN
               WRITE(6,*) ' < VINTERFACE > : index error ',
     &           'ILEFT <> NLEFT*NLBASIS'
               STOP
            END IF
C 
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C                     3.  Summation in the RIGHT bulk side
C        vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
C ++++++++++++++++++++++++++++++++ loop over all occupants of RIGHT host
            IRIGHT = 0
            DO IH1 = 1,NRIGHT
               DO IB = 1,NRBASIS
                  IRIGHT = IRIGHT + 1
                  IREC = IRIGHT + NRIGHTALL*(ILAY1-1) + NRIGHTOFF
                  READ(69,REC=IREC) AVMAD
C
                  IATOM = NLBASIS + IB 
                  DO LM = 1,LMPOT
                     DO LM2 = 1,LMPOT
                        AC(LM) = AC(LM) 
     &                       + AVMAD(LM,LM2)*CMOMHOST(LM2,IATOM)
                     END DO
                  END DO
C
                  IF ( (IH1.EQ.1) .AND. (IB.EQ.1) ) AC(1) = AC(1) + 
     &                 (AVMAD(1,1)-MONOPOL(ILAY1))*CHARGE(2)
               END DO           
            END DO              
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            IF ( IRIGHT.NE.NRIGHTALL ) THEN
               WRITE(6,*) ' < VINTERFACE > : index error ',
     &                'IRIGHT <> NRIGHT*NRBASIS'
               STOP
            END IF
         END IF                 ! (OPT(DECIMATE)
C ######################################################################
C
         WRITE (1337,FMT=99003) IT1,(CATOM(IT1)-ZAT(IT1)),
     &                       (AC(1)/SQRT(4.D0*PI)),
     &                       (AC(3)/SQRT(4.D0*PI))
C ++++++++++++++++++++++++++++++++++++++++++ loop over spins of atom IT1
         DO ISPIN = 1,NSPIN
C
C---> determine the right potential number
C
            IPOT = NSPIN*(IT1-1) + ISPIN
C
C---> in the case of l=0 : r(1)**l is not defined
C
            V(1,1,IPOT) = V(1,1,IPOT) + AC(1)
C
            DO L = 0,LPOT
               DO M = -L,L
                  LM = L*L + L + M + 1
                  DO I = 2,IRS1
                     V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,IT1))**L*AC(LM)
                  END DO
               END DO
            END DO
         END DO                 
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C this part (ICC.GT.0) should be presumably reconsidered for impurity 
C calculation in host-CPA case
C ----------------------------------------------------------------------
         IF ( ICC.GT.0  .or. OPT('KKRFLEX ')) THEN
            DO L = 0,LPOT
               DO M = -L,L
                  LM = L*L + L + M + 1
                  VINTERS(LM,ILAY1) = AC(LM)
               END DO
            END DO
         END IF
C ----------------------------------------------------------------------
      END DO                    
C **********************************************************************
C
      CLOSE (69)
      WRITE(1337,'(15X,45(1H-),/)')
      WRITE(1337,'(79(1H=))')
      IF ( ICC.EQ.0 .and. OPT('KKRFLEX ')==.false.) RETURN
C
C ######################################################################
C
C Now Prepare output for Impurity calculation
C
      OPEN (91,FILE='intercell_ref',STATUS='unknown',FORM='formatted')
      WRITE(1337,*) 
      WRITE(1337,*) '                     ',
     &     'Writing intercell potential for impurity'
      WRITE(1337,'(/,20X,55(1H-))')
      WRITE(1337,99004) HOSTIMP(0),LMPOT
      WRITE(1337,'(20X,55(1H-),/,35X,"  i host lm  Vint")') 
      DO I=1,HOSTIMP(0)
         WRITE(1337,*)
         LM = 1
         WRITE(1337,'(35X,I4,I4,I3,1X,F10.6)') I, HOSTIMP(I),
     &        LM,VINTERS(LM,HOSTIMP(I))
         DO LM=2,9
            WRITE (1337,'(43X,I3,1X,F10.6)') LM,VINTERS(LM,HOSTIMP(I))
         END DO
         WRITE(1337,'(20X,55(1H-))')
      END DO
      WRITE(1337,'(79(1H=),/)')
      
      WRITE(91,99005) HOSTIMP(0),LMPOT 
      DO I=1,HOSTIMP(0)
         WRITE(91,99006) (VINTERS(LM,HOSTIMP(I)),LM=1,LMPOT)
      END DO
      CLOSE(91)  
C
      RETURN
C
99001 FORMAT (79(1H=),/,25X,' INTERFACE MADELUNG POTENTIALS ')
99002 FORMAT (/,15X,' ATOM ','  Delta_Q  ','   MONOPOLE       DIPOLE',
     &        /,15X,45(1H-))
99003 FORMAT (15X,i4,2X,F10.6,1X,1P,D13.6,1X,1P,D13.6)
99004 FORMAT (22X,I4,' host atoms, LMPOT = ',I2,' output up to LM = 9')
99005 FORMAT (3I6)
99006 FORMAT (4D20.10)      
      END
