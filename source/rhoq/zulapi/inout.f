c ************************************************************************
      SUBROUTINE GREFREAD(NATOM,GINP,ITMAT)
c ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      include 'inc.cls'
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER ITMAT,NATOM
      DOUBLE COMPLEX EZ ! added 18.5.2001
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX GINP(LMAXSQ*NACLSD,*)
c     .. local scalars
      INTEGER M,N
c      .. external statement
      LOGICAL TEST
      EXTERNAL RCSTOP,TEST
      SAVE
C     ..
      READ (ITMAT) N
      IF (N.NE.NATOM) THEN 
        write(6,*) 'natom .ne. n',natom,n
        CALL RCSTOP('GREFREAD')
      END IF
      READ (ITMAT) EZ,((GINP(N,M),M=1,LMAXSQ),N=1,LMAXSQ*NACLSD)
      write(6,*) 'Reading GINPBIAS   EZ : ',EZ
c ez added 18.5.2001     
      RETURN

      END
c ************************************************************************
      SUBROUTINE RINPUT1(ALAT,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,ESHIFT,
     +                  FCM,
     +                  HFIELD,MIXING,QBOUND,VCONST,KFG,LMXC,IRNS,
     +                  NTCELL,
     +                  ICST,IFILE,IGF,IMIX,INS,INSREF,IPE,IPF,IPFE,
     +                  IPOTOU,
     +                  IPRCOR,IRM,IRNUMX,ISHIFT,ITCCOR,ITCLST,ITDBRY,
     +                  KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,KSHAPE,KTE,
     +                  KVMAD,KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT,MD,
     +                  NATYP,NSPIN,NATYPD,NSPIND,IEMXD,IRMD,IRNSD,
     +                  LMAXD,LPOTD,TXC,KSCOEF,ICC)
c ************************************************************************
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
      INTEGER IRNS(*),KFG(4,*),LMXC(*),NTCELL(*)
      CHARACTER*24 TXC(3)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,E1,E2,ESHIFT,FCM,HFIELD,MIXING,QBOUND,TK,
     +       VCONST
      INTEGER ICC,ICST,IEMXD,IFILE,IGF,IMIX,INS,INSREF,
     +        IPE,IPF,IPFE,IPOTOU,IPRCOR,
     +        IRM,IRMD,IRNSD,IRNUMX,ISHIFT,ITCCOR,ITCLST,ITDBRY,KCOR,
     +        KEFG,KFROZN,KHFELD,KHYP,KPRE,KSCOEF,KSHAPE,KTE,KVMAD,
     +        KVREL,KWS,KXC,LMAX,LMAXD,LMMAX,LMPOT,LPOT,LPOTD,MD,
     +        NATYP,NATYPD,NPNT1,NPNT2,NPNT3,NPOL,NSPIN,NSPIND
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BRYMIX,STRMIX
      INTEGER I,IL,IP
      CHARACTER*43 TSHAPE
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
c
c---> read input
c
      read (7,FMT=9180) LMAX
      read (7,FMT=9190) E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3
      read (7,FMT=9060) IRNUMX,ITCCOR,IPRCOR
      read (7,FMT=9200) IFILE,IPE,ISHIFT,ESHIFT
      read (7,FMT=9060) KSHAPE,IRM,INS,ICST,INSREF
      read (7,FMT=9060) KCOR,KVREL,KWS,KHYP,KHFELD,KXC
      KFROZN = KCOR
      IF (KCOR.EQ.0) KCOR = 2
      read (7,FMT=9060) KTE,KPRE,KEFG,KVMAD,KSCOEF
c      read (7,FMT=9060) NATYP,KSHAPE
c      read (7,FMT=9060) (NTCELL(I),I=1,NATYP)
c      read (7,FMT=9060) (IRNS(I),I=1,NATYP)
      read (7,FMT=9060) IMIX,IPOTOU,IGF,ICC
      read (7,FMT=9060) ITDBRY
      read (7,FMT=9040) STRMIX,FCM,QBOUND
      read (7,FMT=9040) BRYMIX
      read (7,FMT=9040) HFIELD,VCONST
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
      IF (KHFELD.EQ.1) WRITE (6,FMT=9030) HFIELD
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

      RETURN


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
      END                           ! RINPUT1
c ************************************************************************
      SUBROUTINE RITES(IFILE,NATPS,NATYP,NSPIN,Z,ALAT,RMT,RMTNEW,RWS,
     +                 ITITLE,R,DRDI,VM2Z,IRWS,A,B,TXC,KXC,INS,IRNS,
     +                 LPOT,VINS,QBOUND,IRC,KSHAPE,EFERMI,VBC,ECORE,
     +                 LCORE,NCORE)
c ************************************************************************
c      this subroutine stores in 'ifile' the necessary results
c      (potentials e.t.c.) to start self-consistency iterations
c
c      modified for the full potential case - if ins .gt. 0 there
c       is written a different potential card
c       if the sum of absolute values of an lm component of vins (non
c       spher. potential) is less than the given rms error qbound this
c       component will not be stored .
c
c        (see to subroutine start , where most of the arrays are
c         described)
c
c                            modified by b. drittler  aug. 1988
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.p'
c      INTEGER IRMD,IRNSD,LPOTD
c      PARAMETER (IRMD=424,IRNSD=148,LPOTD=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,QBOUND
      INTEGER IFILE,INS,KSHAPE,KXC,LPOT,NATPS,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),ECORE(20,*),EFERMI, !(2), 22.5,2000
     +                 R(IRMD,*),RMT(*),RMTNEW(*),RWS(*),VBC(2),
     +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
      INTEGER IRC(*),IRNS(*),IRWS(*),ITITLE(20,*),LCORE(20,*),NCORE(*)
      CHARACTER*24 TXC(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,RMAX,RMT1,RMTNW1,RV,SIGN,SUM,Z1
      INTEGER I,ICORE,IH,INEW,IP,IR,IRMIN,IRNS1,IS,ISAVE,J,LM,LMNR,LMPOT,
     +        NCORE1,NR
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRADI(IRMD),ECORE1(20),RA(IRMD),VM2ZA(IRMD)
      INTEGER LCORE1(20)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
c -------------------------------------------------------------------
      ISAVE = 1
      INEW  = 1
c
c
      LMPOT = (LPOT+1)* (LPOT+1)
      DO 60 IH = 1,NATYP
        DO 50 IS = 1,NSPIN
          IF (IS.EQ.NSPIN) THEN
            SIGN = 1.0D0

          ELSE
            SIGN = -1.0D0
          END IF
          IP = NSPIN* (IH-1) + IS

          RMT1 = RMT(IH)
          RMTNW1 = RMTNEW(IH)
          Z1 = Z(IH)
          RMAX = RWS(IH)
          IF (KSHAPE.EQ.0) THEN
            NR = IRWS(IH)

          ELSE
            NR = IRC(IH)
          END IF

          IRNS1 = IRNS(IH)
          IRMIN = NR - IRNS1
          A1 = A(IH)
          B1 = B(IH)
          NCORE1 = NCORE(IP)
c
          DO 10 J = 1,NR
            RA(J) = R(J,IH)
            DRADI(J) = DRDI(J,IH)
c
c--->       store only lm=1 component of the potential
c
            VM2ZA(J) = VM2Z(J,IP)
   10     CONTINUE
c
          IF (NCORE1.GE.1) THEN
c
            DO 20 J = 1,NCORE1
              LCORE1(J) = LCORE(J,IP)
              ECORE1(J) = ECORE(J,IP)
   20       CONTINUE
          END IF
c
c
          WRITE (IFILE,FMT=9000) (ITITLE(I,IP),I=1,7),TXC(KXC+1)
          WRITE (IFILE,FMT=9010) RMT1,ALAT,RMTNW1
          WRITE (IFILE,FMT=9020) Z1,RMAX,EFERMI,VBC(IS)
          WRITE (IFILE,FMT=9030) NR,A1,B1,NCORE1,INEW
          IF (NCORE1.GE.1) WRITE (IFILE,FMT=9040) (LCORE1(ICORE),
     +        ECORE1(ICORE),ICORE=1,NCORE1)
c
          IF (INS.EQ.0 .OR. (IH.LT.NATPS.AND.INS.LE.2)) THEN
c
c--->       store only the spherically averaged potential 
c           (in mt or as - case) 
c           this is done always for the host
c
            IF (INEW.EQ.0) THEN
              WRITE (IFILE,FMT=9050)
     +             (RA(IR),DRADI(IR),VM2ZA(IR),IR=1,NR)
            ELSE
              WRITE (IFILE,FMT=9051) (VM2ZA(IR),IR=1,NR)
            END IF

          ELSE
c
c--->     store the full potential , but the non spherical contribution
c         only from irns1 up to irws1 ;
c         remember that the lm = 1 contribution is multiplied
c         by a factor 1/sqrt(4 pi)
c
            WRITE (IFILE,FMT=9060) NR,IRNS1,LMPOT,ISAVE
            WRITE (IFILE,FMT=9070) (VM2ZA(IR),IR=1,NR)
            IF (LPOT.GT.0) THEN
              LMNR = 1
              DO 40 LM = 2,LMPOT
                SUM = 0.0D0
                DO 30 IR = IRMIN,NR
                  RV = VINS(IR,LM,IP)*RA(IR)
                  SUM = SUM + RV*RV*DRADI(IR)
   30           CONTINUE

                IF (SQRT(SUM).GT.QBOUND) THEN
                  LMNR = LMNR + 1
                  WRITE (IFILE,FMT=9060) LM
                  WRITE (IFILE,FMT=9070) (VINS(IR,LM,IP),IR=IRMIN,NR)
                END IF

   40         CONTINUE
c
c--->         write a one to mark the end
c
              IF (LMNR.LT.LMPOT) WRITE (IFILE,FMT=9060) ISAVE
            END IF

          END IF

   50   CONTINUE
   60 CONTINUE



 9000 FORMAT (7a4,6x,'  exc:',a24,3x,a10)
 9010 FORMAT (3f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f15.10)
 9030 FORMAT (i4,/,2d15.8,/,2i2)
 9040 FORMAT (i5,1p,d20.11)
 9050 FORMAT (1p,2d15.6,1p,d15.8)
 9051 FORMAT (1p,4d20.12)
 9060 FORMAT (10i5)
 9070 FORMAT (1p,4d20.13)
      END
c ************************************************************************
      SUBROUTINE STARTB1(IFILE,IPF,IPFE,IPE,KVREL,KWS,KHFELD,LMAX,
     +                  NBEG,NEND,
     +                  ALAT,RMTNEW,RMT,ITITLE,HFIELD,IMT,IRC,VCONST,
     +                  INS,IRNS,LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELL,
     +                  IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP,EFERMI,
     +                  VBC,C,DROR,RS,S,VM2Z,RWS,ECORE,LCORE,NCORE,DRDI,
     +                  R,Z,A,B,IRWS,INIPOL,IINFO)
c ************************************************************************
c   reads the input potentials
c
c    units :       ry - units for energy
c                  the lattice constant and all other lengths
c                                           given in bohr units
c                  the planck constant h/2pi=1
c                  the electron charge e=sqrt(2)
c                  the electron mass m=1/2
c                  the speed of light c = 2/alpha = 274.0720442
c                      with the fein structure constant alpha
c
c    remember that the input potentials do not include the electro-
c             static contribution of the nucleus of the cell itself
c             this has to be added explicitly !
c
c   as input is used: lmax=maximum angular momentum
c                    nbeg .. nend=number of different atoms
c
c
c     in case of shape corrections this routine  reads from unit 19
c     a suitable radial  mesh 'xrn',its derivate 'drn' and the shape
c     functions 'thetas' .          thus, the region from the muffin
c     tin to the circumscribed  sphere radii is divided  into  'npan'
c     pannels, each one containing 'nm(ipan)' points in order to take
c     care of the  discontinuities of the shape-function  derivative.
c     at the output one obtains :
c            llmsp (icell,ifun)       = integer array giving the com-
c                                       posite  index  lm=l*(l+1)+m+1
c                                       of the ifun-th shape function
c            lmsp  (icell,lm)         = (0,1)  if the lm-th component
c                                       is vanishing or not
c            nfu   (icell)            = number  of   shape   function
c                                       components in cell 'icell'
c
c
c     modified for bandstructure code
c
c                                 b.drittler nov. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      IMPLICIT NONE
      include 'inc.p'
c      INTEGER NATYPD,NSPIND
c      PARAMETER (NATYPD=1,NSPIND=2)
c      INTEGER IRMD,IRNSD,LMAXD,LPOTD
c      PARAMETER (IRMD=424,IRNSD=148,LMAXD=4,LPOTD=8)
c      INTEGER NFUND,IRID
c      PARAMETER (NFUND=24,IRID=75)
c      INTEGER NCELLD,IPAND
c      PARAMETER (NCELLD=1,IPAND=4)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND,INSLPD
      PARAMETER (IRMIND=IRMD-IRNSD,INSLPD= (IRNSD+1)*LMPOTD)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,C,EFERMI,HFIELD,VBC(*),VCONST
      INTEGER IFILE,IINFO,INS,IPE,IPF,IPFE,
     +        KHFELD,KSHAPE,KVREL,KWS,
     +        LMAX,LPOT,
     +        NBEG,NEND,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),DROR(IRMD,*),ECORE(20,*),
     +                 R(IRMD,*),RMT(*),RMTNEW(*),RS(IRMD,0:LMAXD,*),
     +                 RWS(*),S(0:LMAXD,*),THETAS(IRID,NFUND,*),
     +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
      INTEGER IFUNM(NATYPD,*),IMT(*),INIPOL(*),IPAN(*),
     +        IRC(*),IRCUT(0:IPAND,*),
     +        IRMIN(*),IRNS(*),IRWS(*),ITITLE(20,*),
     +        LCORE(20,*),LLMSP(NATYPD,*),LMSP(NATYPD,*),
     +        NCORE(*),NFU(*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,EA,EFNEW,S1,Z1,DUMMY
      INTEGER I,IA,ICELL,ICORE,IFUN,IH,IMT1,INEW,IO,IPAN1,IR,IRC1,IRI,
     +        IRMINM,IRMINP,IRNS1P,IRT1P,IRWS1,ISAVE,ISPIN,ISUM,
     +        J,
     +        L,LM,LM1,LMPOT,LMPOTP,
     +        N,NCELL,NFUN,NR
      LOGICAL TEST
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRN(IRID,NCELLD),SCALE(NCELLD),U(IRMD),
     +                 XRN(IRID,NCELLD)
      INTEGER MESHN(NCELLD),NM(IPAND,NCELLD),NPAN(NCELLD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CALRMT,POTCUT,RCSTOP,RINIT,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ANINT,EXP,LOG,MAX,MOD,REAL,SQRT
C     ..
C     .. Save statement ..
      SAVE
      INTEGER ISHAPE
      DATA ISHAPE / 0 /
C     ..
c-----------------------------------------------------------------------
c
c ---> output of radial mesh information
c
      IO = 0
      IF (IINFO.NE.0 .AND. TEST('RMESH   ')) IO = 1
c
c---> set speed of light
c
      C = 274.0720442D0
      CALL RINIT(INSLPD*(NEND-NBEG+1),VINS(IRMIND,1,NBEG))
c-----------------------------------------------------------------------
c
c---> read radial mesh information of the shape functions and
c     shape functions THETAS in the first iteration - if needed
c
      IF ((KSHAPE.NE.0) .AND. (ISHAPE.EQ.0)) THEN
        ISHAPE = 1
        READ (19,FMT=9000) NCELL
        WRITE (6,FMT=*) '  ncell : ',NCELL,NCELLD
c
        IF(NCELL.GT.NCELLD) THEN
          WRITE(6,*) 'Please, change the parameter ncelld (',NCELLD,
     +         ') in inc.p to',NCELL
          STOP 'STARTB - NCELLD'
        ENDIF
c
        READ (19,FMT=9010) (SCALE(ICELL),ICELL=1,NCELL)
        DO 30 ICELL = 1,NCELL
          READ (19,FMT=9000) NPAN(ICELL),MESHN(ICELL)
c
          IF(NPAN(ICELL)+1.GT.IPAND) THEN
            WRITE(6,*) 'Please, change the parameter ipand (',IPAND,
     +           ') in inc.p to',NPAN(ICELL)+1
            STOP 'STARTB - IPAND'
          ENDIF
c
          IF(MESHN(ICELL).GT.IRID) THEN
            WRITE(6,*) 'Please, change the parameter irid (',IRID,
     +           ') in inc.p to',MESHN(ICELL)
            STOP 'STARTB - IRID'
          ENDIF
c
          READ (19,FMT=9000) (NM(IPAN1,ICELL),IPAN1=2,NPAN(ICELL)+1)
          READ (19,FMT=9010) (XRN(IR,ICELL),DRN(IR,ICELL),IR=1,
     +      MESHN(ICELL))

          READ (19,FMT=9000) NFU(ICELL)
          NFUN = NFU(ICELL)
          WRITE (6,FMT=*) '  nfun  : ',NFUN,NFUND
c
          IF(NFUN.GT.NFUND) THEN
            WRITE(6,*) 'Please, change the parameter nfund (',NFUND,
     +           ') in inc.p to',NFUN
            STOP 'STARTB - NFUND'
          ENDIF
c
          DO 10 LM = 1,LMXSPD
            LMSP(ICELL,LM) = 0
   10     CONTINUE

          DO 20 IFUN = 1,NFUN
            READ (19,FMT=9000) LM
            IF (LM.LE.LMXSPD) THEN
              LLMSP(ICELL,IFUN) = LM
              LMSP(ICELL,LM) = 1
              IFUNM(ICELL,LM) = IFUN
              READ (19,FMT=9010) (THETAS(N,IFUN,ICELL),N=1,MESHN(ICELL))
            ELSE
              READ (19,FMT=9010) (DUMMY,N=1,MESHN(ICELL))
            END IF
   20     CONTINUE

   30   CONTINUE
      END IF                        ! ((KSHAPE.NE.0) .AND. (IFILE.NE.0))
c-----------------------------------------------------------------------
c
      LMPOT = (LPOT+1)* (LPOT+1)
      DO 150 IH = NBEG,NEND
        DO 140 ISPIN = 1,NSPIN
          I = NSPIN* (IH-1) + ISPIN

          IF (IFILE.NE.0) THEN
            IRCUT(0,IH) = 0
            IF (INS.NE.0) THEN
c p.z.            IF (KSHAPE.NE.0) THEN
              
              ICELL = NTCELL(IH)
              IPAN(IH) = 1 + NPAN(ICELL) 
           ELSE

              IPAN(IH) = 1
            END IF
c
c---> read title of potential card
c
c            WRITE(6,*)"before title of potential card"
            READ (IFILE,FMT=9020) (ITITLE(IA,I),IA=1,20)
            IF (IINFO.NE.0) THEN 
              IF (INS.EQ.0) THEN
                WRITE (6,FMT=9080) (ITITLE(IA,I),IA=1,20)
              ELSE
                WRITE (6,FMT=9081) (ITITLE(IA,I),IA=1,20)
              END IF
            END IF
c
c            WRITE(6,*)"after title of potential card"
c---  >read muffin-tin radius , lattice constant and new muffin radius
c      (new mt radius is adapted to the given radial mesh)
c
            READ (IFILE,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
c
c---> read nuclear charge , lmax of the core states ,
c     wigner seitz radius , fermi energy and energy difference
c     between electrostatic zero and muffin tin zero
c
c            READ (IFILE,FMT=9040) Z(IH),RWS(IH),EFNEW,VBC(ISPIN)
            READ (IFILE,*) Z(IH),RWS(IH),EFNEW,VBC(ISPIN)
c
c---> if efermi .eq. 0 use value from in5
c
            IF (EFNEW.NE.0.0D0 .AND. I.EQ.1) EFERMI = EFNEW
c
c---> read : number of radial mesh points
c     (in case of ws input-potential: last mesh point corresponds
c     to ws-radius, in case of shape-corrected input-potential
c     last mesh point of the exponential mesh corresponds to
c     mt-radius/nevertheless this point is always in the array
c     irws(ih)),number of points for the radial non-muffin-tin
c     mesh  needed for shape functions, the constants a and b
c     for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
c     the no. of different core states and some other stuff
c
c            READ (IFILE,FMT=9050) IRWS(IH),A(IH),B(IH),NCORE(I),INEW
            READ (IFILE,FMT=9050) IRWS(IH)
c            WRITE(6,*)"before IRWS"
            READ (IFILE,FMT=9051) A(IH),B(IH),NCORE(I),INEW
c            WRITE(*,*)IRWS(IH)
c
            NR = IRWS(IH)
            IF (NR.GT.IRMD) THEN
              write(6,*) 'Increase parameter IRMD in ''inc.p''',
     +             ' to a value .ge. ',NR,' (= IRWS(',IH,')).'
              STOP 'STARTB1 - IRWS'
            END IF
c
c---> read the different core states : l and energy
c
            IF (NCORE(I).GE.1) READ (IFILE,FMT=9070) (LCORE(ICORE,I),
     +          ECORE(ICORE,I),ICORE=1,NCORE(I))
c
            IF (INS.LT.1) THEN
c
c--->  read radial mesh points, its derivative, the spherically averaged
c      charge density and the input potential without the nuclear pot.
c
              IF (INEW.EQ.0) THEN
                READ (IFILE,FMT=9060) (R(IR,IH),DRDI(IR,IH),VM2Z(IR,I),
     +               IR=1,NR)
c               DO IR=1,NR
c                 WRITE(120,"((4e17.9))") R(IR,IH),DRDI(IR,IH)
c                END DO
              ELSE
                READ (IFILE,FMT=*) (VM2Z(IR,I),IR=1,NR)
              END IF

            ELSE                    ! (INS.LT.1)
c
c--->  read full potential - the non spherical contribution from irmin
c      to irt - remember that the lm = 1 contribution is multiplied by
c      1/sqrt(4 pi)
c
              READ (IFILE,FMT=9090) IRT1P,IRNS1P,LMPOTP,ISAVE
              IRMINP = IRT1P - IRNS1P
              IRMINM = MAX(IRMINP,IRMIND)
              READ (IFILE,FMT=9100) (VM2Z(IR,I),IR=1,NR)
              IF (LMPOTP.GT.1) THEN
                LM1 = 2
                DO 50 LM = 2,LMPOTP
                  IF (LM1.NE.1) THEN
                    IF (ISAVE.EQ.1) THEN
                      READ (IFILE,FMT=9090) LM1

                    ELSE
                      LM1 = LM
                    END IF

                    IF (LM1.GT.1) THEN

                      READ (IFILE,FMT=9100) (U(IR),IR=IRMINP,NR)
                
                      IF (LM1.LE.LMPOT) THEN
                        DO 40 IR = IRMINM,NR
                          VINS(IR,LM1,I) = U(IR)
   40                   CONTINUE
                      END IF

                    END IF

                  END IF

   50           CONTINUE

              END IF

            END IF                  ! (INS.LT.1)
c

            IRWS1 = IRWS(IH)
c
c---> redefine new mt-radius in case of shape corrections
c
            IF (INS.NE.0) THEN
c p.z.      IF (KSHAPE.NE.0) THEN
              RMTNEW(IH) = SCALE(ICELL)*ALAT*XRN(1,ICELL)
              IMT1 = ANINT(LOG(RMTNEW(IH)/B(IH)+1.0D0)/A(IH)) + 1
c
c---> for proper core treatment imt must be odd
c     shift potential by one mesh point if imt is even
c
              IF (MOD(IMT1,2).EQ.0) THEN
                IMT1 = IMT1 + 1
                DO 60 IR = IMT1,2,-1
                  VM2Z(IR,I) = VM2Z(IR-1,I)
   60           CONTINUE
              END IF
c
              IMT(IH) = IMT1
              B(IH) = RMTNEW(IH)/ (EXP(A(IH)*REAL(IMT1-1))-1.0D0)
           ELSE                 ! (INS.NE.0)
             IF (MOD(IRWS1,2).EQ.0) THEN
                 write(6,*) '************************************ '
                 write(6,*) ' For proper core treatment radial    '
                 write(6,*) ' mesh points must be odd.            '
                 write(6,*) ' Please change mesh points in the    '
                 write(6,*) ' potential card. PROGRAM STOPS       '
                 write(6,*) '************************************ '
                 STOP 'VINTERS1 stop 2'
             END IF  
           END IF               ! (INS.NE.0)
c
c---> generate radial mesh - potential only is stored in potential card
c     INEW = 1
c     p. zahn, jan. 99
c
            A1 = A(IH)
            B1 = B(IH)
            R(1,IH) = 0.0D0
            DRDI(1,IH) = A1*B1
            DO 70 IR = 2,IRWS1
              EA = EXP(A1*REAL(IR-1))
              R(IR,IH) = B1* (EA-1.0D0)
              DRDI(IR,IH) = A1*B1*EA
              DROR(IR,IH) = A1/ (1.0D0-1.0D0/EA)
   70       CONTINUE
c
c---> fill cell-type depending mesh points in the non-muffin-tin-region
c
            IF (INS.NE.0) THEN
c p.z.      IF (KSHAPE.NE.0) THEN
              DO 80 IRI = 1,MESHN(ICELL)
                IR = IRI + IMT1
                R(IR,IH) = SCALE(ICELL)*ALAT*XRN(IRI,ICELL)
                DRDI(IR,IH) = SCALE(ICELL)*ALAT*DRN(IRI,ICELL)
                DROR(IR,IH) = DRDI(IR,IH)/R(IR,IH)
   80         CONTINUE
            END IF

            RWS(IH) = R(IRWS1,IH)
c
c---> kshape.eq.0 : calculate new rmt adapted to exp. mesh
c
            CALL CALRMT(IPF,IPFE,IPE,IMT(IH),Z(IH),RMT(IH),RWS(IH),
     +                  RMTNEW(IH),ALAT,DRDI(1,IH),A(IH),B(IH),IRWS1,
     +                  R(1,IH),IO,INS)
c p.z. +                  R(1,IH),IO,KSHAPE)
c
            IF (INS.GT.0) THEN
c p.z.            IF (KSHAPE.GT.0) THEN
              IRCUT(1,IH) = IMT(IH)
              ISUM = IMT(IH)
             
              DO 90 IPAN1 = 2,IPAN(ih)
                ISUM = ISUM + NM(IPAN1,ICELL)
                IRCUT(IPAN1,IH) = ISUM
                
   90         CONTINUE
              NR = ISUM

            ELSE                    ! (KSHAPE.GT.0)

              NR = IRWS(IH)
              IF (KWS.GE.1) THEN
                IRCUT(1,IH) = IRWS1

              ELSE
                IRCUT(1,IH) = IMT(IH)
              END IF

            END IF                  ! (KSHAPE.GT.0)
c
            IRC(IH) = IRCUT(IPAN(IH),IH) 
            
            !!!!!IRC(IH) = NR
c
c---> fill array irmin in case of full potential
c
            IF (INS.NE.0) IRMIN(IH) = NR - IRNS(IH)
c
c---> generate arrays for the calculation of the wave functions
c
            Z1 = Z(IH)
            DO 110 L = 0,LMAX
              IF (KVREL.GE.1) THEN
                S1 = SQRT(REAL(L*L+L+1)-4.0D0*Z1*Z1/ (C*C))
                IF (Z1.EQ.0.0D0) S1 = REAL(L)

              ELSE

                S1 = REAL(L)
              END IF

              S(L,IH) = S1
              RS(1,L,IH) = 0.0D0
              DO 100 IR = 2,NR
                RS(IR,L,IH) = R(IR,IH)**S1
  100         CONTINUE
  110       CONTINUE                ! L = 0,LMAX
c
c---> cut input potential at rmt if given only at exponential mesh
c
            IF (KSHAPE.EQ.1) THEN
              IMT1 = IMT(IH)
              IRC1 = IRCUT(IPAN(IH),IH)
              CALL POTCUT(IMT1,IRC1,INS,LMPOT,R(1,IH),VM2Z(1,I),
     +                    VINS(IRMIND,1,I),Z(IH))
            END IF
c
c--->  first iteration : shift all potentials (only for test purpose)
c
            DO 120 J = 1,NR
              VM2Z(J,I) = VM2Z(J,I) + VCONST
  120       CONTINUE

          END IF                    ! (IFILE.NE.0)

c
          IF (KSHAPE.EQ.0 .AND. KWS.EQ.0) THEN
c
c---> in case of a mt calculation cut potential at mt radius
c
            IMT1 = IMT(IH)
            IRWS1 = IRWS(IH)
            CALL POTCUT(IMT1,IRWS1,INS,LMPOT,R(1,IH),VM2Z(1,I),
     +                  VINS(IRMIND,1,I),Z(IH))

          END IF                    ! KSHAPE.EQ.0 .AND. KWS.EQ.0
c
! ruess: move this out and replace by bshift_ns in main.f to apply also in NS
! part of the potential
!          IF (KHFELD.EQ.1 ) THEN
!c          IF (KHFELD.EQ.1 .AND. NSPIN.EQ.2) THEN
!c
!c--->       maybe apply a magnetic field
!c
!            write(6,*) 'atom',ih,'spin',ispin,'shifted by',
!     +           -REAL(2*ISPIN-3)*HFIELD!*INIPOL(IH)
!            DO 130 J = 1,IRCUT(IPAN(IH),IH)
!              VM2Z(J,I) = VM2Z(J,I) - REAL(2*ISPIN-3)*HFIELD!*INIPOL(IH)
!  130       CONTINUE
!          END IF

  140   CONTINUE                    ! ISPIN = 1,NSPIN

  150 CONTINUE                      ! IH = NBEG,NEND

      RETURN


 9000 FORMAT (16i5)
 9010 FORMAT (4d20.12)
 9020 FORMAT (20a4)
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i4)
 9051 FORMAT (2d15.8,/,2i2)
 9060 FORMAT (1p,2d15.6,1p,d15.8)
 9061 FORMAT (1p,5d15.8)
 9070 FORMAT (i5,1p,d20.11)
c 9080 FORMAT (10x,20a4)
 9080 FORMAT (' < ',20a4)
 9081 FORMAT (' <#',20a4)
 9090 FORMAT (10i5)
 9100 FORMAT (1p,4d20.13)
      END                           ! STARTB1
c ************************************************************************
      SUBROUTINE TMREAD(TMATLL,ITMAT)
c ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMSQ
      PARAMETER (LMSQ= (LMAXD+1)**4)
C     ..
C     .. Scalar Arguments ..
      INTEGER ITMAT
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX TMATLL(LMSQ)
C     ..
      READ (ITMAT) TMATLL
      RETURN

      END
c ************************************************************************
      SUBROUTINE TMWRIT(TMATLL,ITMAT)
c ************************************************************************
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMSQ
      PARAMETER (LMSQ= (LMAXD+1)**4)
C     ..
C     .. Scalar Arguments ..
      INTEGER ITMAT
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX TMATLL(LMSQ)
C     ..
      WRITE (ITMAT) TMATLL
      RETURN

      END
c ************************************************************************
      SUBROUTINE WFWRIT(PNS,QNS,ALPHA,DET,AR,CR,PZ,QZ,FZ,SZ,IOWFCT,INS)
c ************************************************************************
C     .. Parameters ..
      include 'inc.p'
c      INTEGER IRMD,IRNSD,LMAXD,LMX
c      PARAMETER (IRMD=424,IRNSD=148,LMAXD=4,LMX=LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=(LMAXD+1)**2)
      INTEGER IRLLD
      PARAMETER (IRLLD= (IRNSD+1)*LMMAXD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DET
      INTEGER INS,IOWFCT
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),
     +               FZ(IRMD,0:LMAXD),PNS(IRLLD),PZ(IRMD,0:LMAXD),
     +               QNS(IRLLD),QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD)
C     ..
c
c---> store wavefunctions and matrices
c
      IF (INS.NE.0) WRITE (IOWFCT) PNS,QNS
      WRITE (IOWFCT) ALPHA,DET,AR,CR,PZ,QZ,FZ,SZ
      RETURN
      END                           ! WFWRIT
c ************************************************************************
      SUBROUTINE WFREAD(PNS,QNS,ALPHA,DET,AR,CR,PZ,QZ,FZ,SZ,IOWFCT,INS)
c ************************************************************************
C     .. Parameters ..
      include 'inc.p'
c      INTEGER IRMD,IRNSD,LMAXD,LMX
c      PARAMETER (IRMD=424,IRNSD=148,LMAXD=4,LMX=LMAXD+1)
      INTEGER LMX
      PARAMETER (LMX=LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
      INTEGER IRLLD
      PARAMETER (IRLLD= (IRNSD+1)*LMMAXD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DET
      INTEGER INS,IOWFCT
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),
     +               FZ(IRMD,0:LMAXD),PNS(IRLLD),PZ(IRMD,0:LMAXD),
     +               QNS(IRLLD),QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD)
C     ..
c
c---> read wavefunctions and matrices
c
      IF (INS.NE.0) READ (IOWFCT) PNS,QNS
      READ (IOWFCT) ALPHA,DET,AR,CR,PZ,QZ,FZ,SZ
      RETURN
      END                           ! WFREAD
c ************************************************************************
c    EOF inout.f
c ************************************************************************
