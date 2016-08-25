      SUBROUTINE RINPUT99(ALAT,RBASIS,ABASIS,BBASIS,CBASIS,
     +           ALATC,BLATC,CLATC,LATT,CLS,NCLS,
     +           E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,ESHIFT,
     +           DELTA_K_FS,E_FS_IN,                  !new by Swantje
     +           ITCLST,NSTEPS,IMIX,MIXING,QBOUND,FCM,ITDBRY,
     +           IRNS,RMTREF,NTCELL,NAEZ,NEMB,KAOEZ,EQINV,IRM,Z,
     +           NINEQ,NREF,NTCELLR,
     +           ICST,IFILE,IGF,INS,INSREF,IPE,IPF,IPFE,
     +           KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,KSHAPE,KTE,
     +           KFG,KVMAD,KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT, 
     +           NATYP,NSPIN,NSPO,NSPOH,NCL_IMP,IMPLAYER,ILAYERS,
     +           LMXC,TXC,KSCOEF,ICC,REFPOT,
     +           IPOTOU,IPRCOR,IRNUMX,ISHIFT,ITCCOR,
     +           MD,INTERVX,INTERVY,INTERVZ,
     +           HFIELD,COMPLX,
     +           KMT,MTFAC,VBC,VCONST,LINIPOL,INIPOL,IXIPOL,LRHOSYM,
     +           MMIN,MMAX,SINN,SOUT,RIN,ROUT,M2,I12,I13,I19,I25,I40,
     &           NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,    ! new1
     &           TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY, ! new1
     &           NPAN_LOG,NPAN_EQ,NCHEB,R_LOG)                  
c     &       TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY,IRMINSO)                  ! new1
      implicit none
      include 'inc.p'
      include 'inc.cls'
C     .. Local Arrays ..
      CHARACTER*4 TSPIN(2)
      CHARACTER*8 TKWS(3)
      CHARACTER*43 TINS(0:3),TKCOR(0:3),TVREL(0:2)
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Array Arguments ..
      INTEGER IRNS(*),KFG(4,*),LMXC(*),NTCELL(*),CLS(*),REFPOT(*)
      INTEGER EQINV(*),INIPOL(*),IXIPOL(*),NTCELLR(*),KAOEZ(*)
      DOUBLE PRECISION Z(*),MTFAC(*),VBC(*),RBASIS(3,*),RMTREF(*)
      DOUBLE PRECISION TRIGHT(3,*),TLEFT(3,*),ZPERLEFT(3),ZPERIGHT(3) ! new1
      CHARACTER*24 TXC(5)
      CHARACTER*80 UIO
      CHARACTER*40 I13,I19,I25,I40,I12
c      CHARACTER*40 I12,I13,I19,I25,I40
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,E1,E2,ESHIFT,FCM,HFIELD,MIXING,QBOUND,TK,
     +       VCONST,ABASIS,BBASIS,CBASIS,RCUTZ,RCUTXY,
     +       DELTA_K_FS,E_FS_IN
      INTEGER ICC,ICST,IFILE,IGF,IMIX,INS,INSREF,
     +        IPE,IPF,IPFE,IPOTOU,IPRCOR,
     +        IRM,IRNUMX,ISHIFT,ITCCOR,ITCLST,ITDBRY,KCOR,
     +        KEFG,KFROZN,KHFELD,KHYP,KPRE,KSCOEF,KSHAPE,KTE,KVMAD,
     +        KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT,MD,IMPLAYER,ILAYERS,
     +        NATYP,NPNT1,NPNT2,NPNT3,NPOL,NSPIN,NSPO,NCL_IMP,NSPOH,
     +        NPAN_LOG,NPAN_EQ,NCHEB
c     +        IRMINSO
      INTEGER NSTEPS,KMT,NAEZ,NEMB
      INTEGER NINEQ,LATT,NEMBZ,NZ,CENTEROFINV(3)
      DOUBLE PRECISION ALATC,BLATC,CLATC,R_LOG
      INTEGER MMIN,MMAX,SINN,SOUT,RIN,ROUT
      INTEGER INTERVX,INTERVY,INTERVZ,NREF,NCLS
      INTEGER NLBASIS,NRBASIS,NLEFT,NRIGHT          ! new1
      LOGICAL LINIPOL,LRHOSYM,COMPLX,LINTERFACE     ! new1

      LOGICAL OPT           ! added for VELFERMI by Swantje
      EXTERNAl OPT          ! added for VELFERMI by Swantje


C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BRYMIX,STRMIX,E3,TX,TY,TZ
      INTEGER I,IL,IP,J,IER,I1,IR,IC,II,M2
      CHARACTER*43 TSHAPE
c
      CHARACTER*8 TESTC(16),OPTC(8)
      COMMON /TESTC/TESTC
      COMMON /OPTC/OPTC
C     ..
C     .. Data statements ..
      DATA TSPIN/'non-','    '/
      DATA TSHAPE/' exact cell treatment (shape correction)  '/
      DATA TVREL/
     +     ' non relativistic calculation              ',
     +     ' s.r.a. calculation                        ',
     +     ' s.r.a. calculation                        '/
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
C     ..
c
c------------ array set up and definition of input parameter -----------
c
      TXC(1) = ' Morruzi,Janak,Williams '
      TXC(2) = ' von Barth,Hedin        '
      TXC(3) = ' Vosko,Wilk,Nusair      '
      TXC(4) = ' GGA PW91               '
      TXC(5) = ' GGA PW91               '


c      READ (7,1002) (HEAD(I) ,I=1,80)
      IL=1
      CALL IoInput('NSTEPS    ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NSTEPS
      ITCLST = NSTEPS
      CALL IoInput('NSPIN     ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NSPIN
      CALL IoInput('NSPO      ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NSPO 
      CALL IoInput('NSPOH     ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NSPOH 
      CALL IoInput('NCL_IMP   ',UIO,1,7,IER)
                         READ (UNIT=UIO,FMT=*) NCL_IMP
c      CALL IoInput('IRMINSO   ',UIO,1,7,IER)
c                         READ (UNIT=UIO,FMT=*) IRMINSO
c      write(6,2004) dat,vers,head
      write(6,2011) nsteps
      write(6,2104)
      write(6,2010) nspin
      write(6,2104)
      write(6,2009) nspo 
      write(6,2104)
      write(6,2003) nspoh
      write(6,2104)
      write(6,2008) ncl_IMP
      write(6,2104)
      CALL IoInput('NATYP     ',UIO,1,7,IER)
                        READ (UNIT=UIO,FMT=*) NATYP     
      DO I=1,NATYP
         CALL IoInput('ATOMINFO  ',UIO,I+1,7,IER)
                           READ (UNIT=UIO,FMT=*)    Z(I),
     +                        LMXC(I),
     +                       (KFG(J,I),J=1,4),
     +                        CLS(I),
     +                        REFPOT(I),
     +                        NTCELL(I),
     +                        MTFAC(I),
     +                        IRNS(I),
     +                        RMTREF(REFPOT(I))
      END DO
         CALL IoInput('KMT       ',UIO,1,7,IER)
                            READ (UNIT=UIO,FMT=*) KMT
      write(6,2028) natyp
      write(6,2104)
      write(6,1029) (
     +     z(i),
     +     lmxc(i),
     +     (kfg(j,i),j=1,4),
     +     cls(i),
     +     refpot(i),
     +     ntcell(i),
     +     mtfac(i),
     +     irns(i),
     +     rmtref(i),i=1,natyp)
      write(6,2108)
      write(6,2029) kmt
      write(6,2104)
c
c---> read input
c
      IL=1
      CALL IoInput('LMAX      ',UIO,0,7,IER)
                    READ (UNIT=UIO,FMT=*) LMAX
c      WRITE(6,*) "LMAX"
      CALL IoInput('EMIN      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) E1
c      WRITE(6,*) "EMIN"
      CALL IoInput('EMAX      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) E2

      CALL IoInput('TEMPR     ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) TK
c      WRITE(6,*) "TEMP"
      CALL IoInput('NPOL      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPOL
      CALL IoInput('NPT1      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPNT1
      CALL IoInput('NPT2      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) NPNT2
      CALL IoInput('NPT3      ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) npnt3

c      WRITE(6,*) "NPT3"
c new by Swantje for the calculation of the fermi velocity.


c       IF (OPT('VELFERMI') .OR. OPT('TRIANGLE')) then
c       IF (OPT('VELFERMI')) then
c         CALL IoInput('E_F_IN    ',UIO,1,7,IER)
c                    READ (UNIT=UIO,FMT=*) E_FS_IN
c          write(6,*) "E_FS_IN", E_FS_IN
cc         write(6,*) "OPT VELFERMI"
c        CALL IoInput('D_K_FS    ',UIO,1,7,IER)
c                    READ (UNIT=UIO,FMT=*) DELTA_K_FS
c       END IF

      CALL IoInput('IRNUMX    ',UIO,1,7,IER)
                    READ (UNIT=UIO,FMT=*) irnumx
      CALL IoInput('ITCCOR    ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) itccor
      CALL IoInput('IPRCOR    ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) iprcor

c      WRITE(6,*) "IPRCOR"
      CALL IoInput('IFILE     ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ifile
      CALL IoInput('IPE       ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ipe
      CALL IoInput('ISHIFT    ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) ishift
      CALL IoInput('NPAN_LOG  ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) NPAN_LOG
      CALL IoInput('NPAN_EQ   ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) NPAN_EQ
      CALL IoInput('NCHEB     ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) NCHEB
      CALL IoInput('R_LOG     ',UIO,1,7,IER)
                     READ (UNIT=UIO,FMT=*) R_LOG
c      WRITE(6,*) "ISHIFT"
c      CALL IoInput('ESHIFT    ',UIO,1,7,IER)
c                     READ (UNIT=UIO,FMT=*) eshift
      ESHIFT = 0.d0
      CALL IoInput('KSHAPE    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kshape
c      WRITE(6,*) "KSHAPE"
      CALL IoInput('IRM       ',UIO,1,7,IER)
c      WRITE(6,*) "IRM"
                      READ (UNIT=UIO,FMT=*) irm
      CALL IoInput('INS       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) ins
c     WRITE(6,*) "INS"
      CALL IoInput('ICST      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) icst
c     WRITE(6,*) "ICST  "
      CALL IoInput('INSREF    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) insref
c     WRITE(6,*) "INSREF"
     

      CALL IoInput('KCOR      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kcor
      CALL IoInput('KVREL     ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kvrel
      CALL IoInput('KWS       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kws
c     WRITE(6,*) "KWS   "
      CALL IoInput('KHYPERF   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) khyp  
      CALL IoInput('KHFELD    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) khfeld
      CALL IoInput('KEXCOR    ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) kxc    
c     WRITE(6,*) "KEXCOR"
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
c     WRITE(6,*) "KSCOEF"
c
      CALL IoInput('IMIX      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) imix
      CALL IoInput('IPOTOUT   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) ipotou
      CALL IoInput('IGREENFUN ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) igf
      CALL IoInput('ICC       ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) icc
c     WRITE(6,*) "ICC   "
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
c      CALL IoInput('NITER     ',UIO,1,7,IER)
c                      READ (UNIT=UIO,FMT=*) niter 
c     WRITE(6,*) "VCONST"
      CALL IoInput('IMPLAYER  ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) implayer
      CALL IoInput('ILAYERS   ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) ilayers

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
c
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
      IF (KHFELD.EQ.1) WRITE (6,FMT=9030) KHFELD
      WRITE (6,FMT=9050) TSPIN(NSPIN)
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

      E3 = E1
      IF (DABS(TK).GT.1.0D-10 .AND. TK.LT.10.0D0) E3 = TK

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
      write (6,*) "test"
      CALL IoInput('NAEZ      ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) naez
      CALL IoInput('NEMB      ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) nemb 
      CALL IoInput('NEMBZ     ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) nembz 
      NEMB = 0
      NEMBZ = 0
      WRITE(6,*) 'NEMB NEMBZ SET TO ZERO'
      write(6,2023) naez,nemb,nembz
      write(6,2110)
      IF(NAEZ.GT.NAEZD) THEN
        write(6,*) 'Please, increase the parameter naezd (',naezd,
     +       ') in inc.p to',naez
        STOP 'ERROR in NAEZD.'
      ENDIF
      IF(NEMB.GT.NEMBD) THEN
        write(6,*) 'Please, increase the parameter nembd (',nembd,
     +       ') in inc.p to',nemb
        STOP 'ERROR in NEMBD.'
      ENDIF

c      CALL IoInput('NZ        ',UIO,IL,7,IER)
c                      READ (UNIT=UIO,FMT=*) nz
       nz= 0
c      CALL IoInput('LCOMPLEX  ',UIO,IL,7,IER)
c                      READ (UNIT=UIO,FMT=*) complx
       COMPLX=.TRUE.
c      READ(7,1003) COMPLX

      NINEQ = (NAEZ+NZ)/2
      IF (COMPLX) NINEQ = NAEZ
c
      write(6,2024) nz
      write(6,2104)
      write(6,2017) complx
      write(6,2104)
c
      IF (.NOT. COMPLX) THEN
      CALL IoInput('CENTROFIN ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) (centerofinv(i),i=1,3)
c        READ(7,1001) CENTEROFINV
        write(6,2026) (centerofinv(i),i=1,3)
        write(6,2107)
      END IF
c
c!      CALL IoInput('ATOMKIND  ',UIO,1,7,IER)
c!                      READ (UNIT=UIO,FMT=*) KAOEZ(1),
c!     +     (KAOEZ(I),I=2,NINEQ*KAOEZ(1)/MAX(1,KAOEZ(1)))



      IF (KAOEZ(1).eq.0) THEN
        DO I=1,NINEQ
          KAOEZ(I) = I
        END DO
      END IF
c
      DO I=1,NAEZ
c
c --->  atoms equivalent by inversional symmetry
c
        IF (COMPLX) THEN
          CALL IoInput('RBASIS    ',UIO,I,7,IER)
                      READ (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3)
          EQINV(I) = I
        ELSE                        ! (COMPLX)
          IF(I.LE.NINEQ) THEN
            CALL IoInput('RBASIS    ',UIO,I,7,IER)
                      READ (UNIT=UIO,FMT=*) (RBASIS(J,I), J=1,3)
            EQINV(I) = I
          ELSE
            I1 = NAEZ+1+NZ/2-I  ! atom id. by invers. symmetry
            DO J=1,3
              RBASIS(J,I)=2.D0*CENTEROFINV(J)-RBASIS(J,I1)
            ENDDO
            EQINV(I) = I1
            KAOEZ(I) = KAOEZ(I1)
          ENDIF
        END IF                      ! (COMPLX)
      ENDDO                         ! I=1,NAEZ
c
c
c ---> read embedding positions for 1D and 2D systems
c
      IC = (NEMB+NEMBZ)/2
      DO I=1,NEMB
        IF(I.LE.IC) THEN
c          READ(7,1006) (RBASIS(J,NAEZ+I),J=1,3),II
          CALL IoInput('REMBED    ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) (RBASIS(J,NAEZ+I),J=1,3),II
          REFPOT(NATYP+I) = II
          KAOEZ(NAEZ+I) = NATYP + I
        ELSE
          I1 = NEMB+1+NEMBZ/2-I        ! pos. id. by invers. symmetry
          DO J=1,3
            RBASIS(J,NAEZ+I)=2.D0*CENTEROFINV(J)-RBASIS(J,NAEZ+I1)
          ENDDO
          KAOEZ(NAEZ+I) = NATYP + I1
        END IF
      END DO
c
c --- > Read Left and Right host and set up the embending positions
c
c
      CALL IoInput('INTERFACE ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) LINTERFACE
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
c Information is enought to define NEMB
         NEMB = NLBASIS + NRBASIS
         IF(NEMB.GT.NEMBD) THEN
            write(6,*) 'Please, increase the parameter nembd (',nembd,
     +           ') in inc.p to',nemb
            STOP 'ERROR in NEMBD.'
         ENDIF
         DO I=1,NLBASIS
         CALL IoInput('LEFTBASIS ',UIO,I,7,IER)
              READ (UNIT=UIO,FMT=*) (TLEFT(I1,I),I1=1,3),II,IR
         KAOEZ(NATYP+I) = II    ! changed 1.11.99
         REFPOT(NAEZ+I) = IR
         write(6,*) 'this must be 1',KAOEZ(NATYP+I),NATYP+I
         END DO
         DO I=1,NRBASIS
         CALL IoInput('RIGHBASIS ',UIO,I,7,IER)
              READ (UNIT=UIO,FMT=*) (TRIGHT(I1,I),I1=1,3),II,IR
         KAOEZ(NATYP+NLBASIS+I) = II  ! changed 1.11.99
         REFPOT(NAEZ+NLBASIS+I) = IR
         END DO
c
c Put The additional atoms in the "embending" positions
c
         DO I=1,NLBASIS
           !KAOEZ(NAEZ+I) = NATYP + I  ! 2.11.99
           DO I1=1,3
           RBASIS(I1,NAEZ+I) = TLEFT(I1,I)
           END DO
         END DO
         DO I=1,NRBASIS
           !KAOEZ(NAEZ+NLBASIS+I) = NATYP + NLBASIS + I !2.11.99
           DO I1=1,3
           RBASIS(I1,NAEZ+NLBASIS+I) = TRIGHT(I1,I)
           END DO
         END DO 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c In RBASIS we have first the basis atoms or the interface
c atoms then the left host then the right host the host
c goes in the NEMB positions 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         CALL IoInput('ZPERIODL  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) (ZPERLEFT(i1),I1=1,3)
         CALL IoInput('ZPERIODR  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) (ZPERIGHT(i1),I1=1,3) 
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

c      added by Swantje, 08.02.2008
c      set NLBASIS and NRBASIS to zero, in case there is no slab
c      calculation or decimation used
       ELSE
               NLBASIS=0
               NRBASIS=0

       END IF ! LINTERFACE      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      CALL IoInput('RCLUSTZ   ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) RCUTZ
      CALL IoInput('RCLUSTXY  ',UIO,IL,7,IER)
              READ (UNIT=UIO,FMT=*) RCUTXY
      WRITE(6,*) 'Parameters Used for the cluster calculation'
      if (abs(rcutz-rcutxy).lt.1.d-4) then
      write(6,*) 'Clusters inside spheres with radious R = ',rcutz
      else
      write(6,*) 'Clusters inside cylindels with '
      write(6,*) 'Rz = ',rcutz,' Rxy = ',rcutxy
      end if
      write(6,2018)                 ! rbasis
      write(6,2025) ((rbasis(j,i),j=1,3),i,i=1,naez)
      if (nemb.gt.0) write(6,*) 
      write(6,2031) ((rbasis(j,i),j=1,3),i,refpot(kaoez(i)),
     +     i=naez+1,naez+nemb)

      write(6,2101)
      write(6,2012) (kaoez(i),i=1,naez)
      write(6,2103)
      write(6,2030) (eqinv(I),i=1,naez)
      write(6,2103)
c
      CALL IoInput('BASISCALE ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) ABASIS,BBASIS,CBASIS 
      CALL IoInput('LATTICE   ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) LATT
      CALL IoInput('ALATBASIS ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) ALATC,BLATC,CLATC
      CALL IoInput('BZDIVIDE  ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT=*) INTERVX,INTERVY,INTERVZ
      

      IF (LATT.LE.6 .OR. LATT.EQ.10 ) 
     +     BLATC = ALATC
      IF (LATT.LT.3) CLATC = ALATC

      write(6,2019) abasis,bbasis,cbasis 
      write(6,2107)
      write(6,2027) latt
      write(6,2104)
      write(6,2014) alatc,blatc,clatc
      write(6,2107)
      write(6,2015) intervx,intervy,intervz 
      write(6,2102)
c ------------------------------------------------------------------------
      do i=1,nineq
        if (kaoez(i).lt.1) STOP 'Error in KAOEZ'
      enddo
      NCLS = 0
      NREF = 0
      do i=1,nref 
        NTCELLR(I) = 0 
      enddo
c
      DO I=1,NATYP 
        NCLS = MAX(NCLS,CLS(I)) 
      ENDDO
      DO I=1,NATYP 
        NREF = MAX(NREF,REFPOT(I)) 
      ENDDO
      DO I=1,NATYP 
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
c      CALL IoInput('RINOUT   ',UIO,1,7,IER)
c                      READ (UNIT=UIO,FMT=*)  RIN,ROUT


c      READ (7,1005) MMIN,MMAX
c      READ (7,1005) SINN,SOUT,RIN,ROUT
      M2 = LMMAX*NAEZ
      write(6,2013) m2,mmin,mmax,sinn,sout,rin,rout
      write(6,2111)
c
c read test options
c
           
      CALL IoInput('TESTOPT   ',UIO,1,7,IER)
                     READ(UNIT=UIO,FMT=980)(TESTC(i),i=1,8)
      CALL IoInput('TESTOPT   ',UIO,2,7,IER)
                     READ(UNIT=UIO,FMT=980)(TESTC(8+i),i=1,8)
      WRITE(6,52) (TESTC(I),I=1,16)                              
 52   FORMAT(79('-')/' TEST OPTIONS:'/2(1X,A8,7('//',A8)/)/79('-'))
c
c
      CALL IoInput('RUNOPT    ',UIO,1,7,IER)
                   READ (UNIT=UIO,FMT=980)(OPTC(i),i=1,8)
      WRITE(6,62) (OPTC(i),i=1,8)                                 
 62   FORMAT(79('-')/' EXECUTION OPTIONS:'/1X,A8,7('//',A8)/79('-'))
 980  FORMAT(8A8)
c      CALL RTEST(7)
c      CALL ROPT(7)
c
      IL=1

      IF (OPT('TREF_NUM') .EQ. .TRUE. ) THEN
        CALL IoInput('FILES     ',UIO,IL,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I12
      END IF
      CALL IoInput('FILES     ',UIO,IL+1,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I13
      CALL IoInput('FILES     ',UIO,IL+2,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I40
      CALL IoInput('FILES     ',UIO,IL+3,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I19
      CALL IoInput('FILES     ',UIO,IL+4,7,IER)
                      READ (UNIT=UIO,FMT='(A40)')  I25
c      READ (7,FMT='(A40)') I12
c      READ (7,FMT='(A40)') I13
c      READ (7,FMT='(A40)') I40
c      READ (7,FMT='(A40)') I19
c      READ (7,FMT='(A40)') I25

c      write(6,*) 'I12="',I12,'"'
      write(6,*) 'I13="',I13,'"'
      write(6,*) 'I40="',I40,'"'
      write(6,*) 'I19="',I19,'"'
      write(6,*) 'I25="',I25,'"'
      write(6,2100) 
c      WRITE(6,2040) IFULLD,ISPARD,ISLABD
      WRITE(6,2110)
      WRITE(6,*) ' >>>>>>>>> RINPUT99 EXITS NOW <<<<<<<<<< '
      RETURN
C *********************************************Input-End ********
 1000 FORMAT(10I4)
 1001 FORMAT(3F12.9)
 1002 FORMAT(80A1)
 1003 FORMAT(10L4)
 1004 FORMAT(I4)
 1005 FORMAT(4I4)
 1006 FORMAT(3F12.9,I4)
 1029 FORMAT((F4.0,I4,4x,4I1,3I4,F8.4,I4,F8.4))
C ------------------------------------------------------------------------
 2001 FORMAT(/,' ATOM',I3/)
 2002 FORMAT(/,1X,80(A1)///' INPUT DATA :  (NBG=1: LDA; NBG=2: LSDA)')
 2004 FORMAT( 79(1H=)/
     +     'I',77X,'I'/
     +     'I',4X,' Screened Korringa-Kohn-Rostoker ',
     +             'Electronic Structure Code',10X,'I'/
     &     'I',4X,'for Bulk and Interfaces',20X,'I'/
     +     'I',77X,'I',/,'I',4X,'        ',A8,57X,'I'/
     +     'I',77X,'I',/,'I',4X,'Version ',A8,57X,'I'/
     +     'I',77X,'I'/
     +     'I',77X,'I',/,79(1H=),/,/,1X,80A1 )
 2008 FORMAT(' NCLIMP'/I4)
 2009 FORMAT(' NSPO  '/I4)
 2003 FORMAT(' NSPOH '/I4)
 2010 FORMAT(' NSPIN '/I4)
 2011 FORMAT(' NSTEPS'/I4)
 2012 FORMAT(' KAOEZ '/,(10I4))
 2013 FORMAT('      M2    MMIN    MMAX    SINN',
     +       '    SOUT     RIN    ROUT'/7I8)
 2014 FORMAT('          ALATC          BLATC          CLATC'/3F15.8)
 2015 FORMAT('   INTERVX   INTERVY   INTERVZ'/3I10)
 2016 FORMAT('    NCLS    NREF   NINEQ'/,3I8)
 2017 FORMAT(' COMPLX'/(3X,L1))
 2018 FORMAT(' RBASIS')
 2019 FORMAT('         ABASIS         BBASIS         CBASIS'/3F15.8)
 2020 FORMAT(' NPRINCD  NLAYER'/,2I8)
 2021 FORMAT(' INIPOL'/,(10I4))
 2022 FORMAT(' IXIPOL'/,(10I4))
 2023 FORMAT('    NAEZ    NEMB   NEMBZ'/,3I8)
 2024 FORMAT(' NZ    '/,I4)
 2025 FORMAT((3F15.8,I6))
 2026 FORMAT(' COFINV'/,3F15.8)
 2027 FORMAT(' LATT  '/,I4)
 2028 FORMAT(' NATYP '/,I4/,'   Z lmx     KFG cls pot ntc  MTFAC irns
     +  rmtref')
 2029 FORMAT(' KMT   '/,I4)
 2030 FORMAT(' EQINV '/,(10I4))
 2031 FORMAT((3F15.8,2I6))
 2032 FORMAT(' NTCELLR'/,(10I4))
 2040 FORMAT('  IFULLD  ISPARD  ISLABD'/,4I8)
 2050 FORMAT(' Dimension and Input Data CHECK')
C ------------------------------------------------------------------------
 2100 FORMAT(79(1H-))
 2101 format(   3(1H-),1H+  , 3(14(1H-),1H+),  30(1H-))
 2102 format( 3(9(1H-),1H+) ,49(1H-))
 2103 FORMAT(10(3(1H-),1H+) ,39(1H-))
 2104 format(   3(1H-),1H+  ,75(1H-))
 2105 format( 3(3(1H-),1H+) ,67(1H-))
 2106 format( 5(3(1H-),1H+) ,59(1H-))
 2107 format( 3(14(1H-),1H+),34(1H-))
 2108 format( 2(3(1H-),1H+),  7(1H-),1H+,      3(3(1H-),1H+),
     +          7(1H-),1H+,   3(1H-),1H+,      39(1H-))
 2109 format( 5(7(1H-),1H+) ,39(1H-))
 2110 format( 3(7(1H-),1H+) ,55(1H-))
 2111 format( 7(7(1H-),1H+) ,23(1H-))
 2112 format( 2(7(1H-),1H+) ,63(1H-))
 9000 FORMAT (I2,3X,4I5)
 9010 FORMAT (1X,I1,1X,4I1)
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x,
     +       35 ('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x,
     +       'natyp  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6,
     +       ')',/,40x,'nspin  : (',i6,',',i6,')',/)
 9030 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=',
     +       f8.5)
 9040 FORMAT (3f12.7)
 9050 FORMAT (20x,a4,'spin polarized calculation')
 9060 FORMAT (8i4)
 9070 FORMAT (1x,20x,' calculation with',a8,'-potential')
 9080 FORMAT (1x,79 ('*'))
 9090 FORMAT (' mixing factor used           :',f15.6,/,
     +        ' convergence quality required :',1p,d15.2)
 9100 FORMAT (1x,20x,a24,'exchange-correlation potential')
 9110 FORMAT (/,20x,'broyden"s method # :',i3,
     +       ' is used up to iteration-      ',/,20x,'depth :',i3,
     +       '  then jacobian is fixed and potential      ',/,20x,
     +       'is updated using that jacobian')
 9120 FORMAT (13x,' in case of calculating non - spherical wavefcts ',
     +       'the parameter lmaxd has to be set equal lmax ! ')
 9130 FORMAT (/)
 9140 FORMAT (20x,'full potential calculation ',
     +       '- cut off of non spherical potential',/,' >',/)
 9150 FORMAT (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9160 FORMAT (21x,a43,/,21x,' using',i3,'-th. born approximation ')
 9170 FORMAT (21x,a43)
 9180 FORMAT (2i5)
 9190 FORMAT (3f12.7,/,4i4)
 9200 FORMAT (3i4,1f12.7)
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
      END





