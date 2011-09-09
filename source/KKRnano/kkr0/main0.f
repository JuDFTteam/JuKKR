      PROGRAM MAIN0
C
C Explanation of most variables follows below
C
C     ALAT                     : lattice constant (in a.u.)
C     ABASIS,BBASIS,CBASIS,    : scaling factors for rbasis
C     E1,E2,                   : energies needed in EMESHT
C     HFIELD                   : external magnetic field, for
C                              : initial potential shift in
C                              : spin polarised case
C     TK,                      : temperature
C     VCONST,                  : potential shift
C     BRAVAIS(3,3),            : bravais lattice vectors
C     RECBV(3,3),              : reciprocal basis vectors
C     CLEB(NCLEB,2),           : GAUNT coefficients (GAUNT)
C     RMTREF(NREFD),           : muffin-tin radius of reference system
C     RBASIS(3,NAEZD),         : position of atoms in the unit cell
C                              : in units of bravais vectors
C     RCLS(3,NACLSD,NCLSD),    : real space position of atom in cluster
C     RR(3,0:NRD)              : set of real space vectors (in a.u.)
C     VBC(2),                  : potential constants
C     WG(LASSLD),              : integr. weights for Legendre polynomials
C     YRG(LASSLD,0:LASSLD)     : spherical harmonics (GAUNT2)
C     ZAT(NAEZD)               : nuclear charge
C     INTERVX,INTERVY,INTERVZ, : number of intervals in x,y,z-direction
C                              : for k-net in IB of the BZ
C     ICST,                    : number of Born approximation
C     IEND,                    : number of nonzero gaunt coeffizients
C     IFILE,                   : unit specifier for potential card
C     IPE,IPF,IPFE,            : not real used, IPFE should be 0
C     KHFELD,                  : 0,1: no / yes external magnetic field
C     KVREL,                   : 0,1 : non / scalar relat. calculation
C     LMAX,                    : maximum l component in
C                              : wave function expansion
C     LPOT,                    : maximum l component in
C                              : potential expansion
C     NAEZ,                    : number of atoms in unit cell
C     NCLS,                    : number of reference clusters
C     NPNT1,NPNT2,NPNT3,       : number of E points (EMESHT)
C     NPOL,                    : number of Matsubara Pols (EMESHT)
C     NR,                      : number of real space vectors rr
C     NREF,                    : number of diff. ref. potentials
C     NSPIN,                   : counter for spin directions
C     IGUESS                   : 0,1 : no / yes (sc) initial guess, set
C                              : IGUESSD to 1 in inc.p if needed
C     BCP                      : 0,1 : no / yes bc-preconditioning, set
C                              : BCPD to 1 in inc.p if needed
C     QBOUND                   : exit condition for self-consistent 
C                              : iteration 
C     QMRBOUND                 : exit condition for QMR iterations
C     SCFSTEPS                 : number of scf iterations
C     INIPOL(NAEZD),           : initial spin polarisation
C
C     ATOM(NACLSD,NAEZD),      : atom at site in cluster
C     CLS(NAEZD),              : cluster around atom
C     NACLS(NCLSD),            : number of atoms in cluster
C     EZOA(NACLSD,NAEZD),      : EZ of atom at site in cluster
C     ICLEB(NCLEB,3),          : pointer array
C     RMT(NAEZD)               : muffin-tin radius of true system
C     RMTNEW(NAEZD)            : adapted muffin-tin radius 
C     RWS(NAEZD)               : Wigner Seitz radius
C     ICST                     : the regular non spherical wavefunctions, the
C                              : alpha matrix and the t-matrix in the ICST-th. born approx
C     IMT(NAEZD),              : r point at MT radius
C     IPAN(NAEZD),             : number of panels in non-MT-region
C     IRC(NAEZD),              : r point for potential cutting
C     IRCUT(0:IPAND,NAEZD),    : r points of panel borders
C     IRMIN(NAEZD),            : max r for spherical treatment
C     IRNS(NAEZD)              : number r points for non spher. treatm.
C     IRWS(NAEZD),             : r point at WS radius
C     LMSP(NAEZD,LMXSPD)       : 0,1 : non/-vanishing lm=(l,m) component
C                              : of non-spherical potential
C     LLMSP(NAEZD,NFUND)       : lm=(l,m) of 'nfund'th nonvanishing
C                              : component of non-spherical pot.
C     JEND(LMPOTD,             : pointer array for icleb()
C     LOFLM(LM2D),             : l of lm=(l,m) (GAUNT)
C     NTCELL(NAEZD),           : index for WS cell
C     REFPOT(NAEZD)            : ref. pot. card  at position
C
C     A(NAEZD),B(NAEZD)        : contants for exponential r mesh
C     R(IRMD,NAEZD)            : radial mesh ( in units a Bohr)
C     DRDI(IRMD,NAEZD)         : derivative dr/di
C     THETAS(IRID,NFUND,NCELLD): shape function
C                              :         ( 0 outer space
C                              : THETA = (
C                              :         ( 1 inside WS cell
C                              : in spherical harmonics expansion
C     LCORE(20,NPOTD)          : angular momentum of core states
C     NCORE(NPOTD)             : number of core states
C ----------------------------------------------------------------------

      USE kkr0_interfaces
      USE inputcard_reader

      IMPLICIT NONE
C
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'
C     .. Parameters ..
      INTEGER   LASSLD,LMMAXD,LMPOTD,LMXSPD,LM2D,LRECPOT,
     &          MAXMSHD,NPOTD,NSYMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      PARAMETER (NPOTD= NSPIND*NAEZD)
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      PARAMETER (LASSLD=4*LMAXD)
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      PARAMETER (NSYMAXD=48)
      PARAMETER (MAXMSHD=8)
      PARAMETER (LRECPOT=8*(LMPOTD*(IRNSD+1)+IRMD+20))
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALAT,ABASIS,BBASIS,CBASIS,RMAX,GMAX,
     &                 E1,E2,E2IN,EFERMI,FCM,
     &                 HFIELD,MIXING,PI,QBOUND,
     &                 RCUTXY,RCUTZ,RCUTJIJ,RCUTTRC,
     &                 TK,VCONST,VOLUME0,QMRBOUND
      INTEGER I,I1,ICST,IE,IELAST,IEND,IFILE,
     &        IMIX,INTERVX,INTERVY,INTERVZ,
     &        IER,IPE,IPF,IPFE,IRM,ISHIFT,
     &        ITDBRY,IFUN,KHFELD,KPRE,
     &        KTE,KVMAD,KVREL,KXC,LM,LMAX,
     &        LMPOT,LPOT,MAXMESH,NAEZ,
     &        NCLS,
     &        NPNT1,NPNT2,NPNT3,NPOL,NR,NREF,
     &        NSPIN,NSRA,NSYMAT,IGUESS,BCP
      CHARACTER(len=40) I13,I19
      CHARACTER(len=80) UIO
      LOGICAL       JIJ,LDAU
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX DEZ(IEMXD),DSYMLL(LMMAXD,LMMAXD,NSYMAXD),EZ(IEMXD),
     &               WEZ(IEMXD)
      DOUBLE PRECISION A(NAEZD),B(NAEZD),BRAVAIS(3,3),CLEB(NCLEB,2),
     &                 DRDI(IRMD,NAEZD),
     &                 R(IRMD,NAEZD),
     &                 RBASIS(3,NAEZD),RCLS(3,NACLSD,NCLSD),
     &                 RECBV(3,3),
     &                 RMT(NAEZD),RMTNEW(NAEZD),
     &                 RMTREF(NREFD),RR(3,0:NRD),
     &                 RWS(NAEZD),
     &                 VBC(2),
     &                 THETAS(IRID,NFUND,NCELLD),
     &                 VREF(NAEZD),WG(LASSLD),
     &                 YRG(LASSLD,0:LASSLD,0:LASSLD),
     &                 ZAT(NAEZD)
      INTEGER ATOM(NACLSD,NAEZD),CLS(NAEZD),
     &        EZOA(NACLSD,NAEZD),
     &        ICLEB(NCLEB,3),IFUNM(LMXSPD,NAEZD),
     &        IMT(NAEZD),INIPOL(NAEZD),IPAN(NAEZD),IRC(NAEZD),
     &        IRCUT(0:IPAND,NAEZD),IRMIN(NAEZD),IRNS(NAEZD),
     &        IRWS(NAEZD),ISYMINDEX(NSYMAXD),ITITLE(20,NPOTD),
     &        JEND(LMPOTD,0:LMAXD,0:LMAXD),
     &        LCORE(20,NPOTD),LOFLM(LM2D),
     &        LLMSP(NFUND,NAEZD),LMSP(LMXSPD,NAEZD),
     &        NACLS(NCLSD),NCORE(NPOTD),NFU(NAEZD),
     &        NTCELL(NAEZD),REFPOT(NAEZD)
      INTEGER NUMN0(NAEZD),INDN0(NAEZD,NACLSD)
      INTEGER KMESH(IEMXD)
      INTEGER SCFSTEPS,KFORCE
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD)
      DOUBLE PRECISION GSH(NGSHD),EREF
      LOGICAL EVREF

      LOGICAL LCARTESIAN
C
C
C-----------------------------------------------------------------------
C     ..
C     .. External Functions ..
      LOGICAL OPT,TEST
C     EXTERNAL OPT,TEST
C     ..
C     .. External Subroutines ..
C      EXTERNAL BZKINT0,CLSGEN99,EMESHT,GAUNT,GAUNT2,
C     +         LATTIX99,RINPUT99,SCALEVEC,STARTB1,TESTDIM,SHAPE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DABS,DBLE,DIMAG,LOG,MAX,SQRT
C     ..
C

      TYPE (InputcardParams) :: input_params
      TYPE (InputcardArrays) :: input_arrays

C===================================================================
C next short section used to search for file 'VREF'
C if present - DP-value given in that file will be used as EREF
C===================================================================
      INQUIRE(FILE='VREF',EXIST=EVREF)
      IF (EVREF) THEN
        OPEN(87,FILE='VREF',FORM='formatted')
        READ(87,*) EREF
        CLOSE(87)
      ENDIF
      DO I1 = 1,NAEZD
        VREF(I1) = 8.D0
        IF (EVREF) VREF(I1) = EREF 
      END DO
C===================================================================
C===================================================================
      PI = 4.0D0*ATAN(1.0D0)
      EFERMI = 0.0d0

      CALL createInputcardArrays(input_arrays, NAEZD, NREFD)

C
      CALL readInput(input_params, input_arrays)

      CALL RINPUT99(ALAT,RBASIS,ABASIS,BBASIS,CBASIS,CLS,NCLS,
     &              E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,
     &              SCFSTEPS,IMIX,MIXING,QBOUND,FCM,
     &              ITDBRY,IRNS,NTCELL,NAEZ,IRM,ZAT,
     &              NREF,ICST,IFILE,IPE,IPF,IPFE,
     &              KHFELD,KPRE,KTE,
     &              KVMAD,KVREL,KXC,LMAX,LMPOT,LPOT,
     &              NSPIN,REFPOT,
     &              ISHIFT,INTERVX,INTERVY,INTERVZ,
     &              HFIELD,VBC,VCONST,INIPOL,
     &              I13,I19,
     &              RCUTZ,RCUTXY,RCUTJIJ,JIJ,RCUTTRC,
     &              LDAU,
     &              RMTREF,KFORCE,
     &              IGUESS,BCP,QMRBOUND,LCARTESIAN,RMAX,GMAX)


C     IF(NAEZD.GT.100) IPE=0

C
      E2IN = E2
      NSRA = 1
      IF (KVREL.GE.1) NSRA = 2
C
      CALL TESTDIM(NSPIN,NAEZ,LMAX,IRM,NREF,
     &             IRNS,NCLS)
C
      OPEN (19,FILE=I19,STATUS='old',FORM='formatted')
      OPEN (IFILE,FILE=I13,STATUS='old',
     &                         FORM='formatted')
C
C
      CALL STARTB1(IFILE,IPF,IPFE,IPE,KVREL,KHFELD,LMAX,
     &             1,NAEZ,
     &             RMTNEW,RMT,ITITLE,HFIELD,IMT,IRC,VCONST,
     &             IRNS,LPOT,NSPIN,IRMIN,NTCELL,IRCUT,IPAN,
     &             THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN,LRECPOT,
     &             VBC,RWS,LCORE,NCORE,DRDI,
     &             R,ZAT,A,B,IRWS,INIPOL,1)

      CLOSE(IFILE)
      CLOSE(19)
C
C ----------------------------------------------------------------------
C update Fermi energy, adjust energy window according to running options
C
      IF ( NPOL.EQ.0 ) EFERMI = E2IN
C     a test if E2IN is changed after call to STARTB1
      IF ( DABS(E2IN-E2).GT.1D-10 .AND. NPOL.NE.0 ) E2 = E2IN
C
C --> set up energy contour
C
      CALL EMESHT(EZ,DEZ,IELAST,E1,E2,E2IN,TK,
     &            NPOL,NPNT1,NPNT2,NPNT3,IEMXD)
      DO IE = 1,IELAST
        WEZ(IE) = -2.D0/PI*DEZ(IE)
      END DO

C
      CALL GAUNT2(WG,YRG)
      CALL GAUNT(LMAX,LPOT,WG,YRG,CLEB,LOFLM,ICLEB,IEND,JEND)
C      OPEN(56,FILE='gaunt',FORM='formatted')
C      WRITE(56,FMT='(I10)') IEND
C      DO I = 1,IEND
C      WRITE(56,FMT='(3I5,1P,D25.17)')
C     +            ICLEB(I,1),ICLEB(I,2),ICLEB(I,3),CLEB(I,1)
C      END DO
C      CLOSE(56)
C
C --> setup of GAUNT coefficients C(l,m;l',m';l'',m'') for all 
C     nonvanishing (l'',m'')-components of the shape functions THETAS
C
      CALL SHAPE(LPOT,NAEZ,GSH,ILM,IMAXSH,LMSP,NTCELL,WG,YRG)
C      OPEN(56,FILE='gaunt_shape',FORM='formatted')
C      WRITE(56,FMT='(I10)') IMAXSH(LMPOTD)
C      DO I = 1,IMAXSH(LMPOTD)
C      WRITE(56,FMT='(3I5,1P,D25.17)')
C     +            ILM(I,1),ILM(I,2),ILM(I,3),GSH(I)
C      END DO
C      CLOSE(56)
C
C ================================================ deal with the lattice
      CALL LATTIX99(ALAT,BRAVAIS,
     &              RECBV,VOLUME0,RR,NR,NRD)
C
      CALL SCALEVEC(RBASIS,ABASIS,BBASIS,CBASIS,
     &              NAEZ,BRAVAIS,LCARTESIAN)
C ======================================================================
C
C ======================================================================



      CALL CLSGEN99(NAEZ,RR,NR,RBASIS,CLS,NACLS,REFPOT,ATOM,
     &              EZOA,
     &              RCLS,RCUTZ,RCUTXY,
     &              NUMN0,INDN0)

cxxcpl test dimensions for Jij-calculation ..
      CALL CLSJIJ0(NAEZ,RR,NR,RBASIS,RCUTJIJ,JIJ)
cxxcpl .

CGC2 ===================================================================
CGC2  generate cluster 2 in order to truncate GLLKE, GLLH, etc.
CGC2 ===================================================================
C
      IF (TRC.EQ.1) THEN
        CALL CLSGEN_TRC(NAEZ,RR,NR,RBASIS,
     &                  RCUTTRC,ALAT)
      ENDIF
CGC2 ===================================================================
CGC2 ===================================================================
C
C
C
C ======================================================================
C     setting up kpoints
C ======================================================================
      CALL BZKINT0(NAEZ,
     &             RBASIS,BRAVAIS,RECBV,
     &             NSYMAT,ISYMINDEX,
     &             DSYMLL,
     &             INTERVX,INTERVY,INTERVZ,
     &             IELAST,EZ,KMESH,MAXMESH,MAXMSHD)
C ======================================================================
C ======================================================================
C
C
C ======================================================================
C =        write out information for the other program parts           =
C ======================================================================
C
      IF (IGUESS.GT.IGUESSD) THEN
        STOP ' ERROR: to activate initial guess, set IGUESSD=1'
      ENDIF
      IF (BCP.GT.BCPD) THEN
        STOP ' ERROR: to activate bc-preconditioning, set BCPD=1'
      ENDIF
C
C     Conversion of RMAX and GMAX to units of ALAT
      RMAX = RMAX*ALAT
      GMAX = GMAX/ALAT
C
      CALL TESTDIMLAT(ALAT,BRAVAIS,RECBV,RMAX,GMAX)
C
      CALL WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,E1,E2,EZ,WEZ,
     &              BRAVAIS,RECBV,VOLUME0,RMAX,GMAX,
     &              EFERMI,VBC,
     &              SCFSTEPS,LCORE,NCORE,
     &              NSRA,NAEZ,NREF,NSPIN,LMAX,
     &              NCLS,ICST,IPAN,IRCUT,ALAT,ZAT,R,DRDI,
     &              REFPOT,RMTREF,VREF,IEND,JEND,CLEB,ICLEB,
     &              ATOM,CLS,RCLS,NACLS,LOFLM,
     &              RBASIS,RR,EZOA,
     &              KMESH,MAXMESH,NSYMAT,
     &              DSYMLL,
     &              A,B,IFUNM,ITITLE,
     &              LMSP,NTCELL,THETAS,
     &              LPOT,LMPOT,
     &              IMIX,MIXING,QBOUND,FCM,KPRE,KTE,
     &              KVMAD,KXC,ISHIFT,KFORCE,
     &              LLMSP,IMT,IRC,IRMIN,IRNS,IRWS,RWS,RMT,NFU,
     &              GSH,ILM,IMAXSH,
     &              IEMXD,IRMD,LMPOTD,NPOTD,NAEZD,
     &              LMMAXD,IPAND,NREFD,LMAXD,
     &              NCLEB,NACLSD,NCLSD,LM2D,NRD,
     &              NSYMAXD,IRID,NFUND,
     &              NCELLD,LMXSPD,NGSHD,NUMN0,INDN0,
     &              IGUESS,BCP,QMRBOUND,
     &              NR,RCUTJIJ,JIJ,LDAU,ISYMINDEX)
C ======================================================================
C

      CALL destroyInputCardArrays(input_arrays)

      END
