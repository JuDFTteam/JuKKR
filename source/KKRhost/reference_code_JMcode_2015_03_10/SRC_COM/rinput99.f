






C*** Substituted by rinput 13 containing defaults, 25.9.2013




      SUBROUTINE RINPUT99(ALAT,RBASIS,ABASIS,BBASIS,CBASIS,CLS,NCLS,
     +           E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,
     +           EBOTSEMI,EMUSEMI,TKSEMI,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,
     +           FSEMICORE,ESHIFT,
     +           ITCLST,NSTEPS,IMIX,MIXING,QBOUND,FCM,ITDBRY,
     +           IRNS,NTCELL,NAEZ,NEMB,KAOEZ,IRM,Z,
     +           NINEQ,NREF,NTCELLR,
     +           ICST,IFILE,IGF,INS,INSREF,IPE,IPF,IPFE,
     +           KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,KSHAPE,KTE,
     +           KFG,KVMAD,KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT, 
     +           NATYP,NSPIN,
     +           LMXC,TXC,KSCOEF,ICC,REFPOT,
     +           IPOTOU,IPRCOR,IRNUMX,ISHIFT,ITCCOR,
     +           INTERVX,INTERVY,INTERVZ,
     +           HFIELD,COMPLX,
     +           KMT,MTFAC,VBC,VCONST,LINIPOL,INIPOL,IXIPOL,LRHOSYM,
     +           MMIN,MMAX,SINN,SOUT,RIN,ROUT,M2,I12,I13,I19,I25,I40,
     &           NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,    
     &           TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY,RMTREF,KFORCE,
     &           KMROT,QMTET,QMPHI,NCPA,ICPA,ITCPAMAX,CPATOL,NAT,   
     &           NOQ,IQAT,CONC,SOLVER,SOCSCL,CSCL,KREL,
     &           LOPT,UEFF,JEFF,EREFLDAU,KREADLDAU,
     &           LMAXD,LPOTD,NSPIND,NAEZD,NATYPD,NEMBD,NPRINCD,
     &           IRMD,IRNSD,NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,IVSHIFT)
      IMPLICIT NONE
C     ..
C     .. Parameters
      DOUBLE PRECISION CVLIGHT
      PARAMETER (CVLIGHT=274.0720442D0)
C     ..
C     .. Scalar arguments ..
      INTEGER  KREL,LMAXD,LPOTD,NSPIND,NAEZD,NATYPD,NEMBD,NPRINCD,
     &         IRMD,IRNSD

C     .. Local Arrays ..
      CHARACTER*4 TSPIN(3)
      CHARACTER*8 TKWS(3)
      CHARACTER*43 TINS(0:3),TKCOR(0:3),TVREL(0:2)
      CHARACTER*2  SOCII(-2:-1)
C NOTE of VP : there should be some crosscheck of competing options
C              e.g., XCPL and CONDUCT cannot be done simultaneously
C              neither SOC1 and SOC2 manipulation etc.
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Array Arguments ..
      INTEGER IRNS(*),KFG(4,*),LMXC(*),NTCELL(*),CLS(*),REFPOT(*)
      INTEGER INIPOL(*),IXIPOL(*),NTCELLR(*)
      DOUBLE PRECISION Z(*),MTFAC(*),VBC(*),RBASIS(3,*),RMTREF(*)
      DOUBLE PRECISION TRIGHT(3,NEMBD+1),TLEFT(3,NEMBD+1)
      DOUBLE PRECISION ZPERLEFT(3),ZPERIGHT(3)
C     variables for spin-orbit/speed of light scaling
      DOUBLE PRECISION SOCSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION CSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      CHARACTER*24 TXC(4)
      CHARACTER*256 UIO !NCOLIO=256
      CHARACTER*10 SOLVER
      CHARACTER*40 I12,I13,I19,I25,I40
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION E1,E2,ESHIFT,FCM,HFIELD,MIXING,QBOUND,TK,
     +       VCONST,ABASIS,BBASIS,CBASIS,RCUTZ,RCUTXY
      INTEGER ICC,ICST,IFILE,IGF,IMIX,INS,INSREF,
     +        IPE,IPF,IPFE,IPOTOU,IPRCOR,
     +        IRM,IRNUMX,ISHIFT,ITCCOR,ITCLST,ITDBRY,KCOR,
     +        KEFG,KFROZN,KHFELD,KHYP,KPRE,KSCOEF,KSHAPE,KTE,KVMAD,
     +        KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT,KFORCE,
     +        NATYP,NPNT1,NPNT2,NPNT3,NPOL,NSPIN,
     +        NPAN_LOG,NPAN_EQ,NCHEB
      INTEGER NPOLSEMI,N1SEMI,N2SEMI,N3SEMI
      DOUBLE PRECISION FSEMICORE,EBOTSEMI,EMUSEMI,TKSEMI,R_LOG
      INTEGER NSTEPS,KMT,NAEZ,NEMB
      INTEGER NINEQ,NEMBZ
      DOUBLE PRECISION ALAT
      INTEGER MMIN,MMAX,SINN,SOUT,RIN,ROUT
      INTEGER INTERVX,INTERVY,INTERVZ,NREF,NCLS
      INTEGER NLBASIS,NRBASIS,NLEFT,NRIGHT          ! new1
      LOGICAL LINIPOL,LRHOSYM,COMPLX,LINTERFACE     ! new1
C----------------------------------------------------------------
C     CPA variables. Routine has been modified to look for
C     the token ATOMINFOC and only afterwards, if not found, for the
C     old token ATOMINFO. The only necessary extra information 
C     required is the site IQAT(*,IATOM) on which the atom IATOM 
C     is located and the occupancy (concentration) CONC(IATOM). 
C     The rest of CPA variables are deduced from these two. 
C     The tolerance for the CPA-cycle and the number of CPA iterations
C     can be modified adding the token <CPAINFO> in the input file.
C
      INTEGER NCPA,ICPA(NAEZD)   ! ncpa = 0/1 CPA flag
                                 ! icpa = 0/1 site-dependent CPA flag
      INTEGER ITCPAMAX           ! max. number of CPA iterations
      REAL*8  CPATOL             ! convergency tolerance for CPA-cycle
      INTEGER NAT(NATYPD)        ! number of diff. sites occupied 
                                 ! by a given atom type
      INTEGER NOQ(NAEZD)         ! number of diff. atom types located
                                 ! on a given site
      INTEGER IQAT(NAEZD,NATYPD) ! the site on which an atom is located
      INTEGER KAOEZ(NATYPD,NAEZD+NEMBD) 
                                 ! atom types located at a given site
      REAL*8 CONC(NATYPD)        ! concentration of a given atom 

      CHARACTER*3 CPAFLAG(0:1)
      REAL*8 SUM
      INTEGER IO,IA,IQ,IPRINT
C-----------------------------------------------------------------------
C     Variables storing the magnetization direction information.
C     QMTET/QMPHI(NAEZD) give the angles to which the magnetic moment
C     on a given site is rotated against the z-axis. Default values
C     0.0 and 0.0, i.e., magnetic moment parallel to the z-axis.
C     The angles are read in after the token RBASISANG is found 
C     (sought in input file prior to RBASIS token)
C
C   *  KMROT                                                           *
C   *  0: no rotation of the magnetisation                             *
C   *  1: individual rotation of the magnetisation for every site      *
C   ( see also the routine < FINDGROUP > and ff)
C
      INTEGER KMROT
      REAL*8 QMTET(NAEZD),QMPHI(NAEZD)
C ----------------------------------------------------------------------
C LDA+U
      INTEGER LOPT(NATYPD),KREADLDAU
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD)
C LDA+U
C ----------------------------------------------------------------------
C IVSHIFT test option
      INTEGER IVSHIFT
      LOGICAL TEST
      EXTERNAL TEST
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BRYMIX,STRMIX,TX,TY,TZ,DVEC(10)
      INTEGER I,IL,J,IER,I1,II,IR,M2,IDOSEMICORE
      CHARACTER*43 TSHAPE
      DOUBLE PRECISION SOCSCALE,CTLSCALE
      INTEGER IMANSOC(NATYPD),NASOC,ISP(NATYPD)
      LOGICAL OPT,MANSOC,MANCTL
c
      CHARACTER*8 TESTC(16),OPTC(8),VERSION
      COMMON /TESTC/TESTC
      COMMON /OPTC/OPTC
C     ..
C     .. Data statements ..
      DATA VERSION /'Feb 2005'/
      DATA TSPIN/'non-','    ','    '/
      DATA TSHAPE/' exact cell treatment (shape correction)  '/
      DATA TVREL/
     +     ' non relativistic calculation              ',
     +     ' s.r.a. calculation                        ',
     +     ' fully relativistic calculation            '/
      DATA TKCOR/
     +     ' frozen core approximation                 ',
     +     ' core relaxation s.r.a.                    ',
     +     ' core relaxation nonsra                    ',
     +     ' core relaxation                           '/
      DATA TINS/' spherical averaged input potential        ',
     +     ' non spherical input potential for cluster ',
     +     ' non spherical input potential for cluster ',
     +     ' non spherical input potential             '/
      DATA TKWS/' full mt','   ws   ',' full ws'/
C
      DATA CPAFLAG/' NO','YES'/
      DATA SOCII/'xy','zz'/ 
C     ..
c
c------------ array set up and definition of input parameter -----------
c
      TXC(1) = ' Morruzi,Janak,Williams '
      TXC(2) = ' von Barth,Hedin        '
      TXC(3) = ' Vosko,Wilk,Nusair      '
      TXC(4) = ' GGA PW91               '

      IPRINT = 0
      WRITE (6,2004) VERSION
C
C read RUNNING options
C
      CALL IoInput('RUNOPT    ',UIO,1,7,IER)
                   READ (UNIT=UIO,FMT=980)(OPTC(I),I=1,8)
C
C read TEST options
C
      CALL IoInput('TESTOPT   ',UIO,1,7,IER)
                     READ(UNIT=UIO,FMT=980)(TESTC(i),i=1,8)
      CALL IoInput('TESTOPT   ',UIO,2,7,IER)
                     READ(UNIT=UIO,FMT=980)(TESTC(8+i),i=1,8)

      IL=1
      CALL IoInput('NSTEPS    ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NSTEPS
      ITCLST = NSTEPS
C------------ although NSPIND is fixed to 1 in REL mode,
C             NSPIN should be used as 1 or 2 at this stage 
C             to indicate a non- or spin-polarised potential 
C             that has to be read in. NSPIN is set to 1 before
C             being passed to the subsequent programs.
C             < TESTDIM > has been accordingly modified
C
      CALL IoInput('NSPIN     ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NSPIN
      WRITE(6,2011) NSTEPS
      WRITE(6,2104)
      WRITE(6,2010) NSPIN
      WRITE(6,2104)
      CALL IoInput('NATYP     ',UIO,1,7,IER)
                        READ (UNIT=UIO,FMT=*) NATYP     
C
      CALL IoInput('NAEZ      ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) naez
      IF (NATYP.GT.NATYPD) THEN
               WRITE(6,*) ' set NATYPD to at least ',natyp
               STOP ' IN < RINPUT99 > '
      END IF

      IF (NAEZ.GT.NAEZD) THEN
               WRITE(6,*) ' set NAEZD to at least ',naez
               STOP ' in < RINPUT99 > '
      END IF


      DO I=1,NAEZ
         ICPA(I) = 0
         NOQ(I) = 0
      END DO
C
c
c --->  set CPA parameters, search inputfile for update
c
      CPATOL = 1D-4
      ITCPAMAX = 20
      IER = 0
      CALL IOINPUT('CPAINFO   ',UIO,1,7,IER)
      IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) CPATOL,ITCPAMAX
      
      NCPA = 0
      DO I=1,NATYP
          NAT(I) = 0
      END DO

      DO I=1,NATYP
         IER = 0
         CALL IoInput('ATOMINFOC ',UIO,I+1,7,IER)
         IA = 1
         IF ( IER.EQ.0 ) THEN
                           READ (UNIT=UIO,FMT=*)    Z(I),
     +                        LMXC(I),
     +                        (KFG(J,I),J=1,4),
     +                        J,
     +                        IER,
     +                        NTCELL(I),
     +                        MTFAC(I),
     +                        IRNS(I),
     +                        RMTREF(IER),IQAT(IA,I),CONC(I)
            IQ = IQAT(IA,I)
            REFPOT(IQ) = IER
            CLS(IQ) = J
            NAT(I) = NAT(I) + 1
            NOQ(IQ) = NOQ(IQ) + 1
            IF ( NOQ(IQ) .GT. 1 ) THEN
                ICPA(IQ) = 1
                NCPA = 1
            END IF
            KAOEZ(NOQ(IQ),IQ) = I
         ELSE
            CALL IoInput('ATOMINFO  ',UIO,I+1,7,IER)
                           READ (UNIT=UIO,FMT=*)    Z(I),
     +                        LMXC(I),
     +                        (KFG(J,I),J=1,4),
     +                        J,
     +                        REFPOT(I),
     +                        NTCELL(I),
     +                        MTFAC(I),
     +                        IRNS(I),
     +                        RMTREF(REFPOT(I))
            IQAT(1,I) = I
            CLS(I) = J
            CONC(I) = 1D0
            NOQ(I) = 1
            NAT(I) = 1
            KAOEZ(1,I) = I
         END IF

      END DO

      IF ( NCPA .NE. 0 ) THEN 
          DO IQ=1,NAEZ
              SUM = 0D0
              DO IO=1,NOQ(IQ)
                  SUM = SUM + CONC(KAOEZ(IO,IQ))
              END DO  
              IF ( ABS(SUM-1D0).GT.1D-6) THEN
                  WRITE(6,*) ' SITE ', IQ, ' CONCENTRATION <> 1.0 !'
                  WRITE(6,*) ' CHECK YOUR <ATOMINFO-CPA> INPUT '
                  STOP       ' IN <RINPUT99>'
              END IF
          END DO
      END IF

         CALL IoInput('KMT       ',UIO,1,7,IER)
                            READ (UNIT=UIO,FMT=*) KMT
      WRITE(6,2028) NATYP
      WRITE(6,2104)
      WRITE(6,1029) (
     +     Z(I),
     +     LMXC(I),
     +     (KFG(J,I),J=1,4),
     +     CLS(IQAT(1,I)),
     +     REFPOT(IQAT(1,I)),
     +     NTCELL(I),
     +     MTFAC(I),
     +     IRNS(I),
     +     IQAT(1,I),CONC(I),I=1,NATYP)
      WRITE(6,2108)
      WRITE(6,2029) KMT
      WRITE(6,2104)
c
c---> read input
c
      IL=1
      CALL IoInput('LMAX      ',UIO,0,7,IER)
                    READ (UNIT=UIO,FMT=*) LMAX
      CALL IoInput('EMIN      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) E1
      CALL IoInput('EMAX      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) E2

      CALL IoInput('TEMPR     ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) TK

      CALL IoInput('NPOL      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPOL
      CALL IoInput('NPT1      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPNT1
      CALL IoInput('NPT2      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPNT2
      CALL IoInput('NPT3      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPNT3
C
C -> semicore 
C
C initialise variables
      IDOSEMICORE = 0
      EBOTSEMI = E1
      EMUSEMI = EBOTSEMI
      NPOLSEMI = 0
      N1SEMI = 0
      N2SEMI = 0
      N3SEMI = 0
      FSEMICORE = 1D0
C
      IER = 0
      IF ( OPT('SEMICORE') ) THEN
         CALL IoInput('EBOTSEMI  ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                       READ (UNIT=UIO,FMT=*) EBOTSEMI
         CALL IoInput('EMUSEMI   ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) EMUSEMI
C
C -> EMUSEMI < EBOT
         IF ( EMUSEMI.GE.E1 ) GOTO 99800
         CALL IoInput('TKSEMI    ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) TKSEMI

         CALL IoInput('NPOLSEMI  ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) NPOLSEMI
         CALL IoInput('N1SEMI    ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) N1SEMI
         CALL IoInput('N2SEMI    ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) N2SEMI
         CALL IoInput('N3SEMI    ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) N3SEMI
         CALL IoInput('FSEMICORE ',UIO,1,7,IER)
         IF ( IER.NE.0 ) GOTO 99800
                      READ (UNIT=UIO,FMT=*) FSEMICORE
         IDOSEMICORE = 1
99800    CONTINUE
         IF ( IDOSEMICORE.EQ.0 ) THEN 
            WRITE (6,*)
            WRITE (6,*) ' WARNING: SEMICORE used',
     &           ' with incomplete/incorrect contour description'
            WRITE (6,*) ' Running option SEMICORE will be ignored'
            WRITE (6,*)
            DO I=1,8 
               IF (OPTC(I)(1:8).EQ.'SEMICORE') OPTC(I)='        '
            END DO
         END IF
      END IF
C      
      CALL IoInput('IRNUMX    ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) irnumx
      CALL IoInput('ITCCOR    ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) itccor
      CALL IoInput('IPRCOR    ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) iprcor

      CALL IoInput('IFILE     ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ifile
      CALL IoInput('IPE       ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ipe
      CALL IoInput('ISHIFT    ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ishift
      IF (OPT('NEWSOSOL')) THEN
      CALL IoInput('NPAN_LOG  ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) NPAN_LOG
      CALL IoInput('NPAN_EQ   ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) NPAN_EQ
      CALL IoInput('NCHEB     ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) NCHEB
      CALL IoInput('R_LOG     ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) R_LOG
      ENDIF
      ESHIFT = 0.D0
      IF ( OPT('rigid-ef') .OR. OPT('DECIMATE') ) THEN
        ISHIFT = 2
        WRITE(6,*) ' Rigid Fermi Energy, ISHIFT is set to ',ISHIFT
      END IF
C
      CALL IoInput('KSHAPE    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) KSHAPE
      IF ( (KREL.EQ.1).AND.(KSHAPE.NE.0) ) THEN
          WRITE(6,*) 
     &          ' WARNING : KSHAPE set to ZERO for REL case'
          KSHAPE = 0
       END IF
C
      CALL IoInput('IRM       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) irm
      CALL IoInput('INS       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) ins
      CALL IoInput('ICST      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) icst
      CALL IoInput('INSREF    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) insref
     

      CALL IoInput('KCOR      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kcor
      CALL IoInput('KVREL     ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kvrel
      CALL IoInput('KWS       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kws
      CALL IoInput('KHYPERF   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) khyp  
      CALL IoInput('KHFIELD   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) khfeld
      CALL IoInput('KEXCOR    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kxc    
c -------------------------------------------------
      CALL IoInput('KTE       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kte
      CALL IoInput('KPRE      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kpre
      CALL IoInput('KEFG      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kefg
      CALL IoInput('KVMAD     ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kvmad
      CALL IoInput('KSCOEF    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kscoef
c
      CALL IoInput('IMIX      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) imix
      CALL IoInput('IPOTOUT   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) ipotou
      CALL IoInput('IGREENFUN ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) igf
      CALL IoInput('ICC       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) icc
      CALL IoInput('ITDBRY    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) itdbry
      CALL IoInput('STRMIX    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) strmix
      CALL IoInput('FCM       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) fcm
      CALL IoInput('QBOUND    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) qbound
      CALL IoInput('BRYMIX    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) brymix
      CALL IoInput('HFIELD    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) hfield
      CALL IoInput('VCONST    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) vconst
      IF (TEST('atptshft')) THEN
        write(*,*) 'READ IN IVSHIFT'
      CALL IoInput('IVSHIFT   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) ivshift
      ENDIF
C----------------------------------------------------------------------
C
C --> determination of properties at Fermi level
C
      IF ( OPT('GF-EF   ') ) THEN
         IGF = 1
         IF (NPOL.GT.0) NPOL = 0
         IF (NPOL.LT.0) THEN 
            NPNT1 = 0
            NPNT3 = 0
         END IF
         NPNT2 = 1
      END IF
C
      IF (OPT('DOS-EF  ')) THEN
        NPOL = 0
        NPNT2 = 1
      END IF
C ----------------------------------------------------------------------
      IF ( ICC.NE.0 .AND. IGF.EQ.0 ) IGF = 1
      IF ( ICC.EQ.0 .AND. IGF.NE.0 ) ICC = -1
      IF ( ICC.EQ.0 .AND. KSCOEF.NE.0 ) THEN
         KSCOEF = 0
         WRITE(6,*) ' WARNING : electron density calculation.'
         WRITE(6,*) '           ',
     &              'No input of shell structure. KSCOEF set to 0'
      END IF
C ---------------------------------------------------------------------
C     
C     Force calculation 18.5.2000
C
      KFORCE = 0
      IF (INS.GT.0) THEN      
         IER = 0
         CALL IoInput('KFORCE    ',UIO,1,7,IER)
         IF (IER.EQ.0) READ (UNIT=UIO,FMT=*) KFORCE
      END IF

      KFROZN = KCOR
      IF (KCOR.EQ.0) KCOR = 2
c ------------------------------------------------------------------------
      WRITE (6,9210) LMAX
      WRITE (6,9301)
      WRITE (6,9220) E1,E2,TK
      WRITE (6,9302)
      WRITE (6,9230) NPOL,NPNT1,NPNT2,NPNT3
      WRITE (6,9304)
      WRITE (6,9240) IRNUMX,ITCCOR,IPRCOR
      WRITE (6,9303)
      WRITE (6,9250) IFILE,IPE,ISHIFT,ESHIFT
      WRITE (6,9305)
      WRITE (6,9260) KSHAPE,IRM,INS,ICST,INSREF
      WRITE (6,9309)
      WRITE (6,9270) KCOR,KVREL,KWS,KHYP,KHFELD,KXC
      WRITE (6,9306)
      WRITE (6,9330) KTE,KPRE,KEFG,KVMAD,KSCOEF
      WRITE (6,9309)
      WRITE (6,9290) IMIX,IPOTOU,IGF,ICC
      WRITE (6,9304)
      WRITE (6,9300) ITDBRY
      WRITE (6,9307)
      WRITE (6,9310) STRMIX,FCM,QBOUND
      WRITE (6,9302)
      WRITE (6,9320) BRYMIX
      WRITE (6,9308)
      WRITE (6,9280) HFIELD,VCONST
c ------------------------------------------------------------------------

c ------------------------------------------------------------------------
c      DO 10 IP = 1,NATYP
c        READ (5,FMT=9010) LMXC(IP), (KFG(IL,IP),IL=1,4)
c        WRITE (6,FMT=9010) LMXC(IP), (KFG(IL,IP),IL=1,4)
c   10 CONTINUE
c
      IF (KSHAPE.NE.0) KWS = 2
c
      IPF = 6
      IPFE = IPF + 3
C
      IF (OPT('SEARCHEF')) THEN
         IMIX=0
         MIXING=0.0d0
         STRMIX=MIXING
         ITDBRY=1
         QBOUND=1.0d-10
         WRITE(6,'(1X,A)') 'Option SEARCHEF used overriding INPUT for'
         WRITE(6,'(1X,A)')
     &        'IMIX,MIX,QBOUND,ITDBRY: 0, 0.0, 1E-10, 1'
         WRITE(6,*)
      ENDIF

      IF (QBOUND.LT.1.D-15) QBOUND = 1.D-4
      IF (IMIX.GT.2) THEN
        FCM = 1.0D0
        MIXING = BRYMIX
      ELSE
        MIXING = STRMIX
      END IF
c
      IF (IMIX.GE.6) WRITE (6,FMT=9110) (IMIX-5),ITDBRY - 1
c
      WRITE (6,FMT=9090) MIXING,QBOUND
c--------------------------------------------------------
      WRITE (6,FMT=9091) CPAFLAG(NCPA)
      IF (NCPA.NE.0) WRITE(6,9092) ITCPAMAX,CPATOL
c--------------------------------------------------------
c
      LMMAX = (LMAX+1)**2
      LPOT  = MIN(2*LMAX,LPOTD)
      LMPOT = (LPOT+1)* (LPOT+1)
c
      WRITE (6,FMT=9020) LMAX,LMAXD,NATYP,NATYPD,IRM,IRMD,NSPIN,NSPIND

c      IF (LMAX.GT.LMAXD .OR. NATYP.GT.NATYPD .OR. IRM.GT.IRMD .OR.
c     +    NSPIN.GT.NSPIND) CALL RCSTOP('18      ')

      IF (INS.GT.0) THEN
        WRITE (6,FMT=9130)
        WRITE (6,FMT=9140)
        DO 20 I = 1,NATYP
          WRITE (6,FMT=9150) I,IRNS(I),IRNSD

          IF (IRNS(I).GT.IRNSD) CALL RCSTOP('19      ')

   20   CONTINUE

        IF (LMAX.NE.LMAXD) THEN
          WRITE (6,FMT=9120)

          CALL RCSTOP('20      ')

        END IF

      END IF


      WRITE (6,FMT=9130)
c
c
c
c
      IF (KHFELD.EQ.1) WRITE (6,FMT=9030) HFIELD
      if (kvrel.le.1 ) then 
          WRITE (6,FMT=9050) TSPIN(NSPIN)
      else
          write (6,fmt=9050) tspin(nspin+1)
      end if
      WRITE (6,FMT=9170) TVREL(KVREL)
      WRITE (6,FMT=9170) TKCOR(KFROZN)
      IF (KSHAPE.EQ.0) THEN
        WRITE (6,FMT=9070) TKWS(KWS+1)

      ELSE
        WRITE (6,FMT=9170) TSHAPE
      END IF

      WRITE (6,FMT=9100) TXC(KXC+1)
      IF (INS.GT.0) WRITE (6,FMT=9160) TINS(INS),ICST
      WRITE (6,FMT=9080)

c
      VBC(1) = VCONST
      VBC(2) = VBC(1)

      CALL IoInput('LINIPOL   ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) linipol
c      READ(7,1003) LINIPOL
      IF (LINIPOL) THEN
c        READ (7,1000) INIPOL
        CALL IoInput('XINIPOL   ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) (inipol(I),I=1,natyp) 
      ELSE
        DO I=1,NATYP
          INIPOL(I) = 0
        END DO
      END IF

      write (6,2021) (inipol(i),i=1,natyp) 
      write(6,2103)

      CALL IoInput('LRHOSYM   ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) lrhosym

      IF ( (NCPA.NE.0).AND.LRHOSYM ) THEN
         WRITE(6,*)
     &        ' WARNING : CHARGE SYMMETRISATION NOT ALLOWED FOR CPA '
         WRITE(6,*) '          YOUR SETTING IN INPUT FILE IS OVERRIDDEN'
         LRHOSYM = .FALSE.
      END IF
           

      IF (LRHOSYM) THEN

        CALL IoInput('IXIPOL    ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) (ixipol(I),I=1,natyp) 
        write (6,2022) (ixipol(i),i=1,natyp) 
        write (6,2103)
        DO I=1,NATYP
          IF ( IXIPOL(I).NE.0 .AND. 
     +         ABS(IXIPOL(ABS(IXIPOL(I)))).NE.I) THEN
            write(6,*) 'Error in IXIPOL at atom ',I,'.'
            stop 'IXIPOL'
          END IF
        END DO
      ELSE
        DO I=1,NATYP
          IXIPOL(I) = 0
        END DO
        write (6,2022) (ixipol(i),i=1,natyp) 
        write (6,2103)
      END IF
      CALL IoInput('NAEZ      ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) naez
      NEMB = 0
      NEMBZ = 0
      write(6,2023) naez,nemb,nembz
      write(6,2110)
      IF(NAEZ.GT.NAEZD) THEN
        write(6,*) 'Please, increase the parameter naezd (',naezd,
     +       ') in inc.p to',naez
        STOP 'ERROR in NAEZD.'
      ENDIF

       COMPLX=.TRUE.

       NINEQ = NAEZ
c
C
C----------------------------------------------------------------------
C
      KMROT = 0

      DO I=1,NAEZ
c
c --->  atoms equivalent by inversional symmetry
c
        QMTET(I)=0D0
        QMPHI(I)=0D0
        IER = 0 
          CALL IoInput('RBASISANG ',UIO,I,7,IER)

          IF( IER.EQ.0 ) THEN
             READ (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3),
     &             QMTET(I),QMPHI(I)
          ELSE
             CALL IoInput('RBASIS    ',UIO,I,7,IER)
             READ (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3)
          END IF

        IF( ABS(QMTET(I)) .GT. 1D-6 ) KMROT = 1
        IF( ABS(QMPHI(I)) .GT. 1D-6 ) KMROT = 1

      ENDDO                         ! I=1,NAEZ
      CALL IDREALS(RBASIS(1,1),3*NAEZ,IPRINT)
c ---------------------------------------------------------
c -                   Add virtual sites                   -
c----------------------------------------------------------
c ----- add virtual positions I25 file (usually called scoef)
c ----- virtual positions are site with a vanshing t-matrix
c ----- and a muffin tin potential of 0 Ryd. It is used to calculate
c ----- the Greens function for sites which are used to place atoms
c ----- impurity KKR code.
c----------------------------------------------------------
c ----- The subroutine does the following :
c -----  * read the scoef file
c -----  * checks if the real space vectors are lattice positions
c -----  * if not it creates a new basis vector
c -----  * adds the basis vector to the list of basis vectors
c -----  * increases NATYP
c -----  * sets IQAT, CLS, CONC, NOQ, NAT, KAOEZ for CPA=0
c -----    calculations







c
c --- > Read Left and Right host and set up the embeding positions
c
c
      LINTERFACE = .FALSE.
      NRIGHT=  10    
      NLBASIS=  1
      NLEFT=  10    
      NRBASIS=  1
C
      CALL IoInput('INTERFACE ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) LINTERFACE
C----------------------------------------------------------------------
      IF (LINTERFACE) THEN

         WRITE(6,9410)
         CALL IoInput('NRIGHTHO  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) NRIGHT
         CALL IoInput('NLEFTHOS  ',UIO,IL,7,IER)
               READ (UNIT=UIO,FMT=*) NLEFT
         CALL IoInput('NLBASIS   ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) NLBASIS
         CALL IoInput('NRBASIS   ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) NRBASIS
c Information is enough to define NEMB
         NEMB = NLBASIS + NRBASIS
         IF(NEMB.GT.NEMBD) THEN
           write(6,*) 'Please, increase the parameter nembd (',nembd,
     &            ') in inc.p to',nemb
            STOP 'ERROR in NEMBD.'
         ENDIF

         DO I=1,NLBASIS
            CALL IoInput('LEFTBASIS ',UIO,I,7,IER)
                READ (UNIT=UIO,FMT=*) (TLEFT(I1,I),I1=1,3),II,IR
            KAOEZ(1,NAEZ+I) = II    ! changed 1.11.99
            REFPOT(NAEZ+I) = IR    
         END DO
         CALL IDREALS(TLEFT,3*(NEMBD+1),IPRINT)
C
         DO I=1,NRBASIS
         CALL IoInput('RIGHBASIS ',UIO,I,7,IER)
              READ (UNIT=UIO,FMT=*) (TRIGHT(I1,I),I1=1,3),II,IR
         KAOEZ(1,NAEZ+NLBASIS+I) = II  ! changed 1.11.99
         REFPOT(NAEZ+NLBASIS+I) = IR  
         END DO
         CALL IDREALS(TRIGHT,3*(NEMBD+1),IPRINT)
c
c Put The additional atoms in the "embeding" positions
c
         DO I=1,NLBASIS
           DO I1=1,3
             RBASIS(I1,NAEZ+I) = TLEFT(I1,I)
           END DO
         END DO
         DO I=1,NRBASIS
           DO I1=1,3
             RBASIS(I1,NAEZ+NLBASIS+I) = TRIGHT(I1,I)
           END DO
         END DO 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c In RBASIS we have first the basis atoms or the interface
c atoms then the left host then the right host the host
c goes in the NEMB positions 
c
c IN CASE OF CPA the host is treated as an effective 
C CPA medium, that is, there is only one kind of atom
C occupying a crystallographic site.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         CALL IoInput('ZPERIODL  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) (ZPERLEFT(i1),I1=1,3)
         CALL IDREALS(ZPERLEFT(1),3,IPRINT)
C
         CALL IoInput('ZPERIODR  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) (ZPERIGHT(i1),I1=1,3) 
         CALL IDREALS(ZPERIGHT(1),3,IPRINT)
C
         WRITE(6,9430) NLEFT,NLBASIS
         WRITE(6,9440) NRIGHT,NRBASIS
         WRITE(6,9450) (ZPERLEFT(i1),I1=1,3)
         WRITE(6,9460) (ZPERIGHT(i1),I1=1,3)
         WRITE(6,9465)
         WRITE(6,9470)
         DO I=NLEFT,1,-1
            DO I1=NLBASIS,1,-1
            tx = TLEFT(1,i1) + (I-1)*ZPERLEFT(1)
            ty = TLEFT(2,i1) + (I-1)*ZPERLEFT(2)
            tz = TLEFT(3,i1) + (I-1)*ZPERLEFT(3)
            WRITE(6,9420) (I-1)*NLBASIS+i1, tx,ty,tz
            END DO 
         END DO
          WRITE(6,9475)
         DO I=1,NAEZ
            WRITE(6,9420) I, (RBASIS(I1,I),I1=1,3)
         END DO
          WRITE(6,9480)
          DO I=1,NRIGHT
            DO I1=1,NRBASIS
            tx = TRIGHT(1,i1) + (I-1)*ZPERIGHT(1)
            ty = TRIGHT(2,i1) + (I-1)*ZPERIGHT(2)
            tz = TRIGHT(3,i1) + (I-1)*ZPERIGHT(3) 
            WRITE(6,9420) (I-1)*NRBASIS+i1,tx,ty,tz
            END DO 
         END DO  

      END IF                  ! LINTERFACE      
C 
C----------------------------------------------------------------------
C 
      CALL IoInput('RCLUSTZ   ',UIO,IL,7,IER)
                READ (UNIT=UIO,FMT=*) RCUTZ
      CALL IoInput('RCLUSTXY  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) RCUTXY
      WRITE(6,*) 'Parameters used for the cluster calculation'
      if (abs(rcutz-rcutxy).lt.1.d-4) then
      write(6,*) 'Clusters inside spheres with radius R = ',rcutz
      else
      write(6,*) 'Clusters inside cylinders with '
      write(6,*) 'Rz = ',rcutz,' Rxy = ',rcutxy
      end if
      write(6,2104)
      write(6,2018)                 ! rbasis
      write(6,2101)
      do i=1,naez
        write(6,2025) i,(rbasis(j,i),j=1,3),
     &       QMTET(I),QMPHI(I),ICPA(I),NOQ(I),(KAOEZ(J,I),J=1,NOQ(I))
      enddo         
c-------------------------------------------------------------

      if (nemb.gt.0) write(6,*) 
      write(6,2031) ((rbasis(j,i),j=1,3),i,refpot(i),
     +     i=naez+1,naez+nemb)

c
      CALL IoInput('BASISCALE ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) (DVEC(I),I=1,3)
      CALL IDREALS(DVEC(1),3,IPRINT)
      ABASIS = DVEC(1)
      BBASIS = DVEC(2) 
      CBASIS = DVEC(3)

      CALL IoInput('ALATBASIS ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) ALAT
      CALL IoInput('BZDIVIDE  ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) INTERVX,INTERVY,INTERVZ
      
      WRITE(6,2019) ABASIS,BBASIS,CBASIS 
      WRITE(6,2107)
      WRITE(6,2014) ALAT
      WRITE(6,2104)
      WRITE(6,2015) INTERVX,INTERVY,INTERVZ 
      WRITE(6,2102)
c ------------------------------------------------------------------------
      IF ( .not. OPT('VIRATOMS') ) THEN
        DO I=1,NAEZ
          DO IO=1,NOQ(I)
            IF (KAOEZ(IO,I).LT.1) STOP 'Error in KAOEZ'
          END DO
        ENDDO
      END IF

      NCLS = 0
      NREF = 0
      do i=1,nref 
        NTCELLR(I) = 0 
      enddo
c
      DO I=1,NATYP 
c       NCLS = MAX(NCLS,CLS(I))  this is wrong: cls has dimension naezd, not natypd (fivos)
        NCLS = MAX(NCLS,CLS(IQAT(1,I))) 
      ENDDO
      IF ( OPT('VIRATOMS') ) THEN
c        NCLS=NAEZ
c        NINEQ=NAEZ
      END IF


      DO I=1,NAEZ
        NREF = MAX(NREF,REFPOT(I)) 
      ENDDO
      DO I=1,NAEZ
        IF (NTCELLR(REFPOT(I)).EQ.0) NTCELLR(REFPOT(I)) = NTCELL(I)
      ENDDO
c
      WRITE(6,2016) NCLS,NREF,NINEQ
      WRITE(6,2110)
      WRITE(6,2032) (NTCELLR(I),I=1,NREF)
      WRITE(6,2103)

c ------------------------------------------------------------------------
      CALL IoInput('MMIN      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*)  MMIN
      CALL IoInput('MMAX      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*)  MMAX
      CALL IoInput('SRINOUT   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*)  SINN,SOUT,RIN,ROUT
      M2 = LMMAX*NAEZ
      WRITE(6,2013) M2,MMIN,MMAX,SINN,SOUT,RIN,ROUT
      WRITE(6,2111)
C
Check for DECIMATE consistency
C
      IF (OPT('DECIMATE')) THEN
         IF ( MOD(NPRINCD,NLBASIS).NE.0 ) THEN
            WRITE(6,*) ' Decimation cannot continue '
            WRITE(6,*) 'NPRINCD=',NPRINCD,' NLBASIS=',NLBASIS
            STOP
         END IF
         IF ( MOD(NPRINCD,NRBASIS).NE.0 )  THEN
            WRITE(6,*) ' Decimation cannot continue '
            WRITE(6,*) 'NPRINCD=',NPRINCD,' NRBASIS=',NRBASIS
            STOP
         END IF
      END IF
C
Check for ITERMDIR consistency -- if KMROT=0 suppress it
C
      IF ( (OPT('ITERMDIR')).AND.(KMROT.EQ.0) ) THEN
         WRITE (6,*)
         WRITE (6,*)
     &        ' WARNING: ITERMDIR running option used with collinear/',
     &        'parallel Oz starting'
         WRITE (6,*)  
     &        '          system (KMROT = 0 ). Please check token',
     &        ' RBASISANG in your input'
         WRITE (6,*) ' Running option ITERMDIR will be ignored'
         WRITE (6,*)
         DO I=1,8 
            IF (OPTC(I)(1:8).EQ.'ITERMDIR') OPTC(I)='        '
         END DO
      END IF
C
Check for XCPL consistency 
C
      MANCTL = ( KMROT.EQ.0 ).AND.( KREL.EQ.0 ).AND.( NPOL.NE.0 ).AND.
     &         ( NSPIN.GT.1 )
      IF ( (OPT('XCPL    ') ).AND.( .NOT.MANCTL ) ) THEN
         WRITE (6,*)
         WRITE (6,*)
     &        ' WARNING: XCPL running option requires collinear ',
     &        'magnetic systems, complex'
         WRITE (6,*)  
     &        '          energy contour (NPOL<>0) in a NON/SCALAR',
     &        ' relativistic mode (KREL=0)'
         WRITE (6,*) ' Running option XCPL will be ignored'
         WRITE (6,*)
         DO I=1,8 
            IF (OPTC(I)(1:8).EQ.'XCPL    ') OPTC(I)='        '
         END DO
      END IF
      IF ( ( OPT('XCPL    ') ).OR.( OPT('CONDUCT ') ) ) ICC = -1
C
Check for LDA+U consistency -- if INS=0 suppress it
      IF ( (OPT('LDA+U   ')).AND.(INS.EQ.0) ) THEN
         WRITE (6,*)
         WRITE (6,*)
     &        ' WARNING: LDA+U should be used only in NON-SPHERICAL',
     &        ' case (INS=1) '
         WRITE (6,*) ' Running option LDA+U will be ignored'
         WRITE (6,*)
         DO I=1,8 
            IF (OPTC(I)(1:8).EQ.'LDA+U   ') OPTC(I)='        '
         END DO
      END IF
C
      WRITE(6,62) (OPTC(I),I=1,8)                                              
 62   FORMAT(79('-')/' EXECUTION OPTIONS:'/1X,A8,7('//',A8)/79('-'))
      WRITE(6,52) (TESTC(I),I=1,16)                                    
 52   FORMAT(79('-')/' TEST OPTIONS:'/2(1X,A8,7('//',A8)/)/79('-'))
 980  FORMAT(8A8)
C
C ---------------------------------------------------------------------
C Initialise SOLVER, SOC and CTL parameters in REL case
C
C 
      IF (KREL.EQ.1) THEN
         SOLVER='BS        '
C
         IER=0
         CALL IoInput('SOLVER    ',UIO,0,7,IER)
         IF (IER.EQ.0) THEN
              READ (UNIT=UIO,FMT=*) SOLVER
              IF ( SOLVER(1:2) .EQ. 'BS' ) THEN
                 SOLVER = 'BS        '
              ELSE
                 IF( SOLVER .NE. 'ABM-OP    ' ) SOLVER='ABM-OP    '
              END IF
         END IF
C
         DO I=1,NATYP
            DO IL=1,LMAXD+1
               SOCSCL(IL,I) = 1D0
               CSCL(IL,I) = CVLIGHT
            END DO
         END DO
         MANSOC=.FALSE.
         MANCTL=.FALSE.
C
C ============================================================= SOC-MAN
C
         IF (OPT('SOC     ')) THEN
            IER=0
            CALL IOInput('SOCSCALE  ',UIO,0,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) SOCSCALE
               IF (SOCSCALE.GT.-2.5D0) THEN 
                  IF (SOCSCALE.GE.0.0D0) THEN           ! SOC-I
                     SOLVER='ABM-SOC   '
                     MANSOC=.TRUE.
                  ELSE                                  ! SOC-II
                     SOLVER       = 'ABM-SOC-II'
                     MANSOC=.TRUE.
                     DO I=1,NATYP
                        DO IL=1,LMAXD+1
                           SOCSCL(IL,I) = SOCSCALE
                        END DO
                     END DO
                     WRITE(6,99010) SOCII(NINT(SOCSCALE))
                  END IF
               ELSE
                  WRITE(6,99001) '< SOC >'
                  WRITE(6,99003)
               END IF
            ELSE
               WRITE(6,99002) '< SOC >'
               WRITE(6,99003)
            END IF
C
            IF ( MANSOC .AND. (SOCSCALE.GE.0D0) ) THEN
               DO I=1,NATYP
                  IMANSOC(I) = 1
               END DO
C
C ---> now look for a possible include/exclude list (SOCLIST= +/- NASOC)
C ---> if SOCLIST is not found, ALL the atoms will have SOC modified with
C ---> SOCSCALE (+NASOC=only NASOC atoms, -NASOC=all but these NASOC atoms)
C      Note that this is allowed only for SOC-I manipulation
C              
               IER=0
               CALL IOInput('SOCLIST   ',UIO,0,7,IER)
               IF (IER.EQ.0) THEN
                  READ(UNIT=UIO,FMT=*) NASOC,(ISP(I),I=1,ABS(NASOC))
                  
                  IF (NASOC.NE.0) THEN
                     IF (NASOC.LT.0) THEN ! exclude this atoms
                        DO I=1,-NASOC
                           IMANSOC(ISP(I)) = 0
                        END DO
                     ELSE
                        DO I=1,NATYP
                           IMANSOC(I) = 0
                        END DO
                        DO I=1,NASOC
                           IMANSOC(ISP(I)) = 1
                        END DO
                     END IF
                  END IF
               END IF
C     
               WRITE(6,2100)
               DO I=1,NATYP
                  IF (IMANSOC(I).EQ.1) THEN
                     DO IL=1,LMAXD+1
                        SOCSCL(IL,I)=SOCSCALE
                     END DO
                  END IF
               END DO
               WRITE(6,99004)
               IF (NASOC.EQ.0) WRITE(6,99005)
               IF (NASOC.GT.0) THEN
                  WRITE(6,99006)
                  WRITE(6,99008) (ISP(I),I=1,NASOC)
               END IF
               IF (NASOC.LT.0) THEN
                  WRITE(6,99007)
                  WRITE(6,99008) (ISP(I),I=1,ABS(NASOC))
               END IF
               WRITE(6,99009) SOCSCALE
               WRITE(6,2100)
            END IF
         END IF
C
C ============================================================= SOC-MAN
C
         WRITE(6,'('' SOLVER used for the DIRAC equation : '',2X,A)')
     &        SOLVER
         WRITE(6,2100)
C
C ============================================================= CTL-MAN
C
         IF (OPT('CSCALE  ')) THEN
            IER=0
            CALL IOInput('CTLSCALE  ',UIO,0,7,IER)
            IF (IER.EQ.0) THEN
               READ (UNIT=UIO,FMT=*) CTLSCALE
               IF (CTLSCALE.GE.1D-12) THEN 
                  MANCTL=.TRUE.
               ELSE
                  WRITE(6,99001) '< CSCALE >'
                  WRITE(6,99011)
               END IF
            ELSE
               WRITE(6,99002) '< CSCALE >'
               WRITE(6,99011)
            END IF
C
            IF (MANCTL) THEN
               DO I=1,NATYP
                  DO IL=1,LMAXD+1
                     CSCL(IL,I)=CSCL(IL,I)/CTLSCALE**0.5D0
                  END DO
               END DO
               WRITE(6,99012)
               WRITE(6,99005)
               WRITE(6,99009) 1D0/CTLSCALE**0.5D0
            END IF
            WRITE(6,2100)
         END IF
C
C ============================================================= CTL-MAN
      END IF
C ================================================================ LDA+U
C
C
C -> Initialise UEFF,JEFF,LOPT,EREFLDAU for all atoms
C
      DO I=1,NATYP
         LOPT(I) = -1           !  not perform lda+u (default)
         UEFF(I) = 0.D0
         JEFF(I) = 0.D0
         EREFLDAU(I) = 0.1D0
      ENDDO

      IF (OPT('LDA+U   ')) THEN
C 
C -> get number of atoms for lda+u:
C
         IER = 0
         CALL IoInput('NAT_LDAU  ',UIO,1,7,IER)
         IF ( IER.NE.0 ) THEN 
            NASOC = NATYP
         ELSE 
            READ (UNIT=UIO,FMT=*) NASOC
            IF ( NASOC.GT.NATYP ) STOP ' main0: NAT_LDAU > NATYP'
         END IF
C
C -> read in UEFF,JEFF,LOPT,EREFLDAU for the desired atoms
C
         IL = 0
         DO I=1,NASOC
            IER = 0
            CALL IoInput('LDAU_PARA ',UIO,I,7,IER)
            IF ( IER.EQ.0 ) THEN
               READ (UNIT=UIO,FMT=*) 
     &              I1,LOPT(I1),UEFF(I1),JEFF(I1),EREFLDAU(I1)
               IL = IL + 1
            END IF
         ENDDO
         IF ( IL.NE.NASOC ) THEN
            WRITE(6,*) ' ERROR: LDA+U invoked for ',NASOC,' atoms'
            WRITE(6,*) '        Some (all) parameters are missing',
     &           ' in the input-file'
            STOP 
         END IF
         KREADLDAU = 0
         IER = 0
         CALL IoInput('KREADLDAU ',UIO,1,7,IER)
         IF ( IER.EQ.0 ) READ (UNIT=UIO,FMT=*) KREADLDAU
      END IF
C ======================================================================

C ---------------------------------------------------------------------

      IL=1
      CALL IoInput('FILES     ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I12
      CALL IoInput('FILES     ',UIO,IL+1,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I13
      CALL IoInput('FILES     ',UIO,IL+2,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I40
      CALL IoInput('FILES     ',UIO,IL+3,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I19
      CALL IoInput('FILES     ',UIO,IL+4,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I25
      IF ( OPT('VIRATOMS') ) THEN

!      CALL ADDVIRATOMS(LINTERFACE, naez, naezd, natypd, nembd,rbasis,
!     +                 NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,TLEFT,TRIGHT,
!     +                 NLBASIS,NRBASIS,ABASIS,BBASIS,CBASIS,NEMB, 
!!    +  &                  lmxc, kfg, refpot, ntcell, mtfac, irns, 
!!    +  &                  rmtref, IQAT, cls, conc, noq, nat,
!     +                     refpot,KAOEZ,noq)

      END IF
C
      WRITE(6,*) 'I12="',I12,'"'
      WRITE(6,*) 'I13="',I13,'"'
      WRITE(6,*) 'I40="',I40,'"'
      WRITE(6,*) 'I19="',I19,'"'
      WRITE(6,*) 'I25="',I25,'"'
      WRITE(6,2100) 
      WRITE(6,2040) KMROT
      WRITE(6,2110)
      WRITE(6,*) ' >>>>>>>>> RINPUT99 EXITS NOW <<<<<<<<<< '
      RETURN
C *********************************************Input-End ********
 1029 FORMAT((F4.0,I4,4x,4I1,3I4,F8.4,I4,I5,1x,f5.3))
C ------------------------------------------------------------------------
 2004 FORMAT( /79(1H*)/
     &     '*',77X,'*'/
     &     '*',10X,'Screened Korringa-Kohn-Rostoker ',
     &             'Electronic Structure Code',10X,'*'/
     &     '*',27X,'for Bulk and Interfaces',27X,'*'/
     &     '*',77X,'*'/
     &     '*',4X,'Juelich-Munich 2001 - 2005',25X,
     &     'Version : ',A8,4X,'*'/
     &     '*',77X,'*'/79(1H*))
 2010 FORMAT(' NSPIN '/I4)
 2011 FORMAT(' NSTEPS'/I4)
 2013 FORMAT('      M2    MMIN    MMAX    SINN',
     +       '    SOUT     RIN    ROUT'/7I8)
 2014 FORMAT('          ALAT = ',F15.8)
 2015 FORMAT('   INTERVX   INTERVY   INTERVZ'/3I10)
 2016 FORMAT('    NCLS    NREF   NINEQ'/,3I8)
 2018 FORMAT(' RBASIS'/,
     &     'SITE                BASIS VECTORS                 ',
     &     'THETA   PHI CPA OCC KAOEZ')
 2019 FORMAT('         ABASIS         BBASIS         CBASIS'/3F15.8)
 2021 FORMAT(' INIPOL'/,(10I4))
 2022 FORMAT(' IXIPOL'/,(10I4))
 2023 FORMAT('    NAEZ    NEMB   NEMBZ'/,3I8)
 2025 FORMAT((i4,3F15.8,2F6.1,2(1x,I3),4I3))
 2028 FORMAT(' NATYP '/,I4/,
     &     '   Z lmx     KFG cls pot ntc  MTFAC irns SITE  CONC')
 2029 FORMAT(' KMT   '/,I4)
 2031 FORMAT((3F15.8,2I6))
 2032 FORMAT(' NTCELLR'/,(10I4))
 2040 FORMAT(' KMROT'/,4I8)
C ------------------------------------------------------------------------
 2100 FORMAT(79(1H-))
 2101 format(   3(1H-),1H+  , 3(14(1H-),1H+),  30(1H-))
 2102 format( 3(9(1H-),1H+) ,49(1H-))
 2103 FORMAT(10(3(1H-),1H+) ,39(1H-))
 2104 format(   3(1H-),1H+  ,75(1H-))
 2107 format( 3(14(1H-),1H+),34(1H-))
 2108 format( 2(3(1H-),1H+),  7(1H-),1H+,      3(3(1H-),1H+),
     +          7(1H-),1H+,   3(1H-),1H+,      39(1H-))
 2110 format( 3(7(1H-),1H+) ,55(1H-))
 2111 format( 7(7(1H-),1H+) ,23(1H-))
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x,
     +       35 ('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x,
     +       'natyp  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6,
     +       ')',/,40x,'nspin  : (',i6,',',i6,')',/)
 9030 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=',
     +       f8.5)
 9050 FORMAT (20x,a4,'spin polarized calculation')
 9070 FORMAT (1x,20x,' calculation with',a8,'-potential')
 9080 FORMAT (1x,79 ('*'))
 9090 FORMAT (' mixing factor used           :',f15.6,/,
     +        ' convergence quality required :',1p,d15.2)
 9091 FORMAT (' make use of CPA algorithm    :',1x,a14)
 9092 FORMAT ('         max. iterations      :',i15,/,
     +        '         req. CPA convergency :',1p,d15.2)
 9100 FORMAT (1x,20x,a24,'exchange-correlation potential')
 9110 FORMAT (/,20x,'broyden"s method # :',i3,
     +       ' is used up to iteration-      ',/,20x,'depth :',i3,
     +       '  then jacobian is fixed and potential      ',/,20x,
     +       'is updated using that jacobian')
 9120 FORMAT (13x,' in case of calculating non - spherical wavefcts ',
     +       'the parameter lmaxd has to be set equal lmax ')
 9130 FORMAT (/)
 9140 FORMAT (20x,'full potential calculation ',
     +       '- cut off of non spherical potential',/,' >',/)
 9150 FORMAT (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9160 FORMAT (21x,a43,/,21x,' using',i3,'-th. born approximation ')
 9170 FORMAT (21x,a43)
 9210 FORMAT (' lmax'/,i4)
 9220 FORMAT ('          E1          E2          TK'/,3f12.6)
 9230 FORMAT ('   NPOL  NPNT1  NPNT2  NPNT3'/,4i7)
 9240 FORMAT (' IRNUMX ITCCOR IPRCOR'/,3i7)
 9250 FORMAT ('  IFILE    IPE ISHIFT ESHIFT'/,3i7,f12.6)
 9260 FORMAT (' KSHAPE    IRM    INS   ICST INSREF'/,5i7)
 9270 FORMAT ('   KCOR  KVREL    KWS   KHYP KHFELD    KXC'/,6i7)
 9280 FORMAT (' external magnetic hfield     :',f15.4/,
     +        ' VCONST                       :',f15.6)
 9290 FORMAT ('   IMIX IPOTOU    IGF    ICC'/,4i7)
 9300 FORMAT (' ITDBRY'/,i7)
 9310 FORMAT ('      STRMIX        FCM       QBOUND'/,3f12.6)
 9320 FORMAT ('      BRYMIX'/,f12.6)
 9330 FORMAT ('    KTE   KPRE   KEFG  KVMAD KSCOEF'/,5i7)
 9301 format(   3(1H-),1H+  ,75(1H-))
 9302 format( 3(11(1H-),1H+),43(1H-))
 9303 format(3(6(1H-),1H+) ,58(1H-))
 9304 format(4(6(1H-),1H+) ,51(1H-))
 9305 format(3(6(1H-),1H+),11(1H-),1H+ ,46(1H-))
 9306 format(6(6(1H-),1H+) ,37(1H-))
 9307 format(6(1H-),1H+,72(1H-))
 9308 format(11(1H-),1H+,67(1H-))
 9309 format(5(6(1H-),1H+) ,44(1H-))
 9410 format('*** SLAB - INTERFACE CALCULATION ***'/)
 9420 format(I5,3F12.6)
 9430 format('Number of LEFT  Host Layers : ',I5,' with ',I5,' basis')
 9440 format('Number of RIGHT Host Layers : ',I5,' with ',I5,' basis')
 9450 format('Left  side periodicity : ',3F10.5)
 9460 format('Right side periodicity : ',3F10.5)
 9465 format('    Geommetry used : '/,
     &       ' ATOM       TX          TY          TZ ')
 9470 format('--------------- Left  Host -------------- ')
 9475 format('---------------   S L A B  -------------- ')
 9480 format('--------------- Right Host -------------- ')
99001 FORMAT(/,1X,
     &     "WARNING: Option ",A," used with an INVALID ",
     &     "scaling parameter.")
99002 FORMAT(/,1X,
     &     "WARNING: Option ",A," found but NO value given for the",
     &     " scaling parameter.")
99003 FORMAT(15X,'++++++++++   SOC option will be IGNORED   ++++++++++',
     &     /,1X,'Please use SOCSCALE= XXX (real>-2.5) in the inputcard',
     &     ' to make your option valid ',/)
99004 FORMAT(1X,'The SOC will be SCALED',$)
99005 FORMAT(' for ALL the atoms in the unit cell.')
99006 FORMAT(' for the FOLLOWING atoms in the unit cell :')
99007 FORMAT(' for all the atoms in the unit cell EXCLUDING :')
99008 FORMAT(1X,6(2X,I3))
99009 FORMAT(1X,'Scaling factor = ',1P,D9.2)
99010 FORMAT(1X,'The SOC is manipulated',' -- part of the SOC kept: ',A)
99011 FORMAT(15X,'+++++++++  CSCALE option will be IGNORED  ++++++++++',
     &     /,1X,'Please use CTLSCALE= X (real>=1D-12) in the inputcard',
     &     ' to make your option valid ',/)
99012 FORMAT(1X,'The CLIGHT will be SCALED',$)
      END
