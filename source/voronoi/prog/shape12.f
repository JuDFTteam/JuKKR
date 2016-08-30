C=====================================================================
  
      SUBROUTINE SHAPE(NPOI8,AFACE8,BFACE8,CFACE8,DFACE8,
     &                 TOLVDIST,   ! Max. tolerance for distance of two vertices
     &                 TOLEULER,   ! Used in calculation of Euler angles, subr. EULER
     &                 NMIN,       ! Min. number of points in panel
     &                 NVERTICES8,XVERT8,YVERT8,ZVERT8,NFACE8,LMAX8,
     &                 DLT8,KEYPAN8,NM8,ICLUSTER,
     &                 NCELL_S,    ! number of cell
     &                 SCALE_S,    ! scaling factor
     &                 NPAN_S,     ! panels for the shape
     &                 MESHN_S,    ! number of mesh points 
     &                 NM_S,       ! array, mesh for each panel
     &                 XRN_S,      ! radial mesh
     &                 DRN_S,      ! drdi for radial mesh
     &                 NFUN_S,     ! number of nonzero shape-functions
     &                 LMIFUN_S,   ! lm of the ifun shapefunction
     &                 THETAS_S)   ! Thetas(r,ifun)
      implicit none
C----------------------------------------------------------------------
C
C                       S H A P E   P R O G R A M 
C
C        F O R  A R B I T R A R Y  V O R O N O I  P O L Y H E D R A
C
C
C                                                           N.Stefanou
C----------------------------------------------------------------------
C     In order to improve the efficency of the code and to
C     check more, this program was changed in summer 1998 by N.Stefanou
C     ATTENTION: BUG is removed! 
C                In the old version IPAN=-1 is wrong and has to be
C                set to IPAN=0. 
C
C     THIS PROGRAM CALCULATES THE ANGULAR MOMENTUM COMPONENTS OF THE 
C     SHAPE FUNCTION FOR AN ARBITRARY VORONOI POLYHEDRON.
C     A REAL SPHERICAL HARMONIC BASIS IS USED FOR THE DECOMPOSITION.
C     ON INPUT WE GIVE :
C
C        LMAX          :  MAXIMUM ANGULAR MOMENTUM
C
C        DLT           :  DEFINES THE STEP FOR GAUSS-LEGENDRE CALC.
C  TIME (DEC)  DLT    TOTAL ENERGY (BCC-TEST CdSb in Ge 12 Shells)
C  1104S      0.002  -.59583341\
C   315S      0.005  -.59583341 \
C    51S      0.050  -.59583341  --> NO CHANGES ALSO IN
C    40S      0.100  -.59583341 /    CHARGES OR FORCES
C    38S      0.200  -.59583341/
C    38S      0.300  -.59583354---> CHARGES DIFFER IN 10NTH DIGIT
C    35S      0.400  -.59583673---> CHARGES DIFFER IN 7NNTH DIGIT
C                                   FORCES DIFFER IN 5TH DIGIT
C    --->TO BE AT THE SAVE SIDE USE 0.05 OR 0.1 (SHOULD BE QUITE GOOD)
C        SIMILAR RESULTS WERE HELD FOR Cu in Fe NN relaxation.
C
C
C        NFACE         :  NUMBER OF FACES OF THE POLYHEDRON
C        KEYPAN        :  KEY TO DEFINE  THE  RADIAL  MESH.  IF KEYPAN=
C                         THE DEFAULT  RADIAL  DIVISION  OF PANNELS GIVE
C                         IN DATA STATEMENT IS USED.OTHERWISE THE  NUMBE
C                         OF  RADIAL  MESH  POINTS PER PANNEL  (NM(IPAN)
C                         IS READ IN INPUT
C                      ** IN THIS VERSION THE MESH IS DETERMINED
C                         BY SUBROUTINE MESH0.
C        Z(I)          :  COEFFICIENTS OF THE EQUATION OF A FACE
C                         Z(1)*X + Z(2)*Y + Z(3)*Z  =  1
C        NVERT         :  NUMBER OF VERTICES OF A FACE
C        V(I,IVERT)    :  COORDINATES OF THE VERTICES OF A FACE
C        NEWSCH(IFACE) :  INTEGER   PARAMETER TO CALCULATE   (=1)
C                         THE CONTRIBUTION OF THE CORRESPONDING
C                         PYRAMID TO THE SHAPE FUNCTIONS  .  IF
C                         NEWSCH.NE.1 THE CONTRIBUTION IS TAKEN
C                         EQUAL TO THAT OF THE PREVIOUS PYRAMID
C
C
C     IN ORDER TO SAVE MEMORY WE STORE IN LOCAL TEMPORARY FILES IN  UNIT
C     30+1 , 30+2 , ... , 30+NFACE THE TRANSFORMATION MATRICES ASSOCIATE
C     WITH THE ROTATION OF EACH PYRAMID. THE TEMPORARY DIRECT ACCESS FIL
C     IN UNIT 10 CONTAINS THE CALCULATED COMPONENTS OF THE SHAPE FUNCTIO
C
C                  ...........I N P U T  C A R D...(Bcc/fcc)
C
C                                       if not (Bcc/fcc) change main prg
C bcc                          <----- Gives the lattice parameters
C    16    1   0.05000                lmax,nkey,division
C                                     LMAX=4*LMAX(KKR), 
C                                     NKEY is not used, 
C                                     DIVISION is DLT      
C   125    0                          number of mesh points,keypan
C                                     Number of mesh points 
C                                     used for the radial mesh (Depends
C                                     on the number of pannels).
C                                     If keypan is 1, 
C                                     then the radial mesh division is 
C                                     taken from the input
C   -3.30000 -3.30000  -3.30000      relaxation percent
C
C    63   32   30    7   21   15   15   15   15   15 \     
C    15   17   15   15   23   15   15    0    0    0  \ This is the
C     0   17   15   15   23   15   15    0    0    0  / radial mesh info
C     0    0    0    0    0    0    0    0    0    0 /
C                  .........................................
C
C
C     THE DEFINITION OF REAL SPHERICAL HARMONICS IS NOT THE STANDARD  ON
C     REFERED IN THE PAPER:
C     N.STEFANOU,H.AKAI AND R.ZELLER,COMPUTER PHYS.COMMUN. 60 (1990) 231
C     IF YOU WANT TO HAVE ANGULAR MOMENTUM COMPONENTS IN THE STANDARD BA
C     SIS CHANGE THE FOLLOWING STATEMENTS IN THE ROUTINES :
C     IN CCOEF      ISI=1                        ---->    ISI=1-2*MOD(M,
C     IN DREAL      IF(MOD(M+MP),2).EQ.0) D=-D   ---->    DELETE THE LIN
C
C
C-----------------------------------------------------------------------
C
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'
      INTEGER    ICD,ICED,IBMAXD,ISUMD,MESHND
      INTEGER    NVTOTD
      PARAMETER (ICD=1729,ICED=((LMAXD1+1)*(LMAXD1+2))/2)
c
c     PARAMETER (IBMAXD=(LMAXD1+1)*(LMAXD1+1),ISUMD=23426)
c
      PARAMETER (IBMAXD=(LMAXD1+1)*(LMAXD1+1),ISUMD=100000)
      PARAMETER (MESHND=IRID)
      PARAMETER (NVTOTD=NFACED*NVERTD)
C
C     .. SCALAR VARIABLES ..
C
      INTEGER     ICE,IFACE,IMAX,IP,IPAN,IPMAX,IREC,IS,ISU,ISUM,IS0
      INTEGER     ITEMP,ITET,IV,IVERT,IVTOT,K,KEYPAN,K0,L,LMAX,M,MESHN
      INTEGER     NFUN,LM0,NCELL,NPOI,NMIN
      INTEGER    MO,MP,N,NFACE,NPAN,NTET,NVERT,NVTOT,I,IB,IBM,IBMAX,IC
      REAL*8      ARG1,ARG2,DLT,FK,FL,FPISQ,R,RAP,RDOWN,SQ3O3,COA
      REAL*8      SCALE,A1,A2,A3,A4
      REAL*8      TOLVDIST,TOLEULER
      LOGICAL    KHCP
C
C     .. ARRAY VARIABLES ..
C
      INTEGER   ISW(IBMAXD),NM(NPAND),LOFM(IBMAXD),MOFM(IBMAXD)
      INTEGER   NEWSCH(NFACED)
      INTEGER   NVERTICES(NFACED)
      REAL*8     CL(ICD)
      REAL*8     V(3,NVERTD),Z(3),CRT(NPAND),DMATL(ISUMD)
      REAL*8     XRN(MESHND),DRN(MESHND),BB(MESHND),C(ICED),B(IBMAXD)
      REAL*8     S3(-LMAXD1:LMAXD1,0:LMAXD1),RUPSQ(NVTOTD)
      REAL*8     S(-LMAXD1:LMAXD1,0:LMAXD1),SUM(0:LMAXD1,2)
      REAL*8     S1(-LMAXD1:LMAXD1,0:LMAXD1),S2(-LMAXD1:LMAXD1,0:LMAXD1)
      REAL*8    AFACE(NFACED),BFACE(NFACED),CFACE(NFACED),DFACE(NFACED)
      REAL*8    XVERT(NVERTD,NFACED),YVERT(NVERTD,NFACED),
     &          ZVERT(NVERTD,NFACED)
 
      integer NPOI8,NVERTICES8(NFACED),NFACE8,LMAX8,KEYPAN8,NM8(NPAND)
      integer ICLUSTER,icount
      REAL*8           AFACE8(NFACED),BFACE8(NFACED),CFACE8(NFACED),
     &                 DFACE8(NFACED)
      REAL*8           XVERT8(NVERTD,NFACED),YVERT8(NVERTD,NFACED),
     &          ZVERT8(NVERTD,NFACED)
      REAL*8           dlt8
      integer j,i1
      logical test
c
c  Basic output
c
      integer NCELL_S,NPAN_S,MESHN_S,NFUN_S 
      integer NM_S(NPAND),LMIFUN_S(IBMAXD)
      REAL*8           XRN_S(MESHND),DRN_S(MESHND),
     &                 THETAS_S(MESHND,IBMAXD) 
      REAL*8           scale_s
C
C     .. SCALARS IN COMMON ..
C
      REAL*8       PI
C
C     .. ARRAYS IN COMMON ..
C
      INTEGER     ISIGNU(NVTOTD),NTT(NFACED)
      REAL*8      ALPHA(NFACED),BETA(NFACED),GAMMA(NFACED),R0(NFACED)
      REAL*8      FA(NVTOTD),FB(NVTOTD),FD(NVTOTD),RD(NVTOTD)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC ABS,ACOS,DMAX1,DMIN1,ATAN,DFLOAT,SQRT
C
C     .. EXTERNAL ROUTINES ..
C
      EXTERNAL CCOEF,CRIT,DREAL,MESH,PINTG,TEST
C
C     .. COMMON BLOCKS ..
C
      COMMON /ANGLES/ PI,ALPHA,BETA,GAMMA
      COMMON /TETRA/  FA,FB,FD,R0,RD,ISIGNU,NTT
C
C     .. DATA STATEMENTS ..
C
c$      DATA NM/NPAND*25/
C-----------------------------------------------------------------------
      PI=4.D0*DATAN(1.D0)
      FPISQ=DSQRT(4.D0*PI)
      KHCP=.false.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      NPOI = NPOI8
      do i=1,nfaced
         NVERTICES(i) = NVERTICES8(i) !(NFACED)
      end do
      NFACE = NFACE8
      LMAX = LMAX8
      KEYPAN = KEYPAN8
      do i=1,npand
         NM(i) = NM8(i)         ! (NPAND)
      end do
      do i=1,nfaced 
         AFACE(i) = AFACE8(i)   ! (NFACED)
         BFACE(i) = BFACE8(i)   ! (NFACED)
         CFACE(i) = CFACE8(i)   ! (NFACED)
         DFACE(i) = DFACE8(i)   ! (NFACED)
         do i1=1,nvertd
            
            XVERT(i1,i) = XVERT8(i1,i) ! (NVERTD,NFACED)
            YVERT(i1,i) = YVERT8(i1,i) ! (NVERTD,NFACED)
            ZVERT(i1,i) = ZVERT8(i1,i) ! (NVERTD,NFACED)
         end do
      end do 
      DLT = dlt8

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C-----------------------------------------
C  THIS CALL DOES SOME GEOMETRICAL TESTS (N.STEFANOU 98)
      CALL POLCHK(NFACE,NVERTICES,XVERT,YVERT,ZVERT,TOLVDIST)
C**********   FOR HCP CASE ONLY    *********
      if (khcp) then
      COA  =SQRT(8.D0/3.D0)
      COA  =2.D0
      SQ3O3=SQRT(3.D0)/3.D0
      end if
C**********   FOR HCP CASE  ONLY   *********
c$      READ(7,104) NFACE,LMAX,KEYPAN,DLT
      IBMAX=(LMAX+1)*(LMAX+1)
      ISUM=0
      DO 19 L=0,LMAX
      IS0=(2*L+1)*(2*L+1)
   19 ISUM=ISUM+IS0
      IF(ISUM.GT.ISUMD . OR . LMAX.GT.LMAXD1) GO TO 200
      IPAN=0 
      IVTOT=0
C.......................................................................
C     S T O R A G E            I N    C O M M O N        B L O C K S
C     C A L C U L A T I O N    O F    R O T A T I O N    M A T R I C E S
C.......................................................................
      DO 1 IFACE=1,NFACE
      ITEMP=30+IFACE
c$      READ(7,103) A1,A2,A3,A4,NVERT,NEWSCH(IFACE)       !1.3.2000
      A1 = AFACE(IFACE)
      A2 = BFACE(IFACE)
      A3 = CFACE(IFACE)
      A4 = DFACE(IFACE)
      NVERT = NVERTICES(IFACE)
      NEWSCH(IFACE) = 1         ! THIS is ALWAYS ONE!
c -----------------
      Z(1)=A1/A4
      Z(2)=A2/A4
      Z(3)=A3/A4

C************    FOR HCP CASE ONLY (TO REMOVE OTHERWISE)   ************
      if (khcp) then
      Z(1)=Z(1)*SQ3O3
      Z(3)=Z(3)*8.D0/COA/3.D0
      end if
C************    FOR HCP CASE ONLY (TO REMOVE OTHERWISE)   ************

      DO 2 IVERT=1,NVERT
c$      READ(7,100) (V(I,IVERT),I=1,3)           ! 1.3.2000
        
        V(1,IVERT) = XVERT(IVERT,IFACE)  
        V(2,IVERT) = YVERT(IVERT,IFACE)        
        V(3,IVERT) = ZVERT(IVERT,IFACE)

C************    FOR HCP CASE ONLY (TO REMOVE OTHERWISE)   ************
      if (khcp) then
      V(1,IVERT)=V(1,IVERT)*SQ3O3
      V(3,IVERT)=V(3,IVERT)*COA
      end if
C************    FOR HCP CASE ONLY (TO REMOVE OTHERWISE)   ************

    2 CONTINUE
      CALL CRIT(IFACE,NVERT,V,Z,IPAN,IVTOT,TOLEULER,TOLVDIST,CRT)
      IF  (TEST('verb0   ')) THEN
         WRITE(6,105) IFACE,NTT(IFACE)
         IF(NEWSCH(IFACE).EQ.0) WRITE(6,110)
      ENDIF
      CALL DREAL(LMAX,ALPHA(IFACE),BETA(IFACE),GAMMA(IFACE),ITEMP)
      REWIND ITEMP
    1 CONTINUE

C.......................................................................
C     D E F I N I T I O N    O F    T H E    S U I T A B L E    M E S H
C.......................................................................
      NVTOT=IVTOT
      NPAN=IPAN
c$      IF(KEYPAN.NE.0) THEN 
c$      READ(7,106) (NM(IPAN),IPAN=1,NPAN-1)
c$      IF (NM(NPAN-1).LT.6) THEN
c$      WRITE(6,*) 'Check number of points for each panel'
c$      STOP
c$      END IF
c$      END IF
      CALL MESH(CRT,NPAN,NM,XRN,DRN,MESHN,NPOI,KEYPAN,NMIN)
      IF (TEST('verb0   ')) THEN 
      WRITE(6,102)
      DO 3 IV=1,NVTOT
      WRITE(6,101) IV,FA(IV)/PI,FB(IV)/PI,FD(IV)/PI,RD(IV),ISIGNU(IV)
    3 CONTINUE
      END IF
C.......................................................................
C     E X P A N S I O N    C O E F F I C I E N T S
C.......................................................................
      CALL CCOEF(LMAX,CL,C)
      IVTOT=0
      DO 21 IFACE=1,NFACE
      NTET=NTT(IFACE)
      DO 21 ITET=1,NTET
      IVTOT=IVTOT+1
      RUPSQ(IVTOT)=SQRT((RD(IVTOT)-R0(IFACE))*(RD(IVTOT)+R0(IFACE)))
   21 CONTINUE
      DO 27 IBM=1,IBMAX
   27 ISW(IBM)=0
      OPEN(11,STATUS='SCRATCH',ACCESS='DIRECT',
     *        FORM='UNFORMATTED',RECL=80)
C.......................................................................
C     L O O P    O V E R    R A D I A L    M E S H    P O I N T S
C.......................................................................
      DO 12 N=1,MESHN
      R=XRN(N)
      DO 9 IBM=1,IBMAX
    9 B(IBM)=0.D0
      IVTOT=0
C.......................................................................
C     L O O P    O V E R    P Y R A M I D S
C.......................................................................
      DO 13 IFACE=1,NFACE
      NTET=NTT(IFACE)
      ITEMP=30+IFACE
      IF(R.GT.R0(IFACE))  GO TO 31
      IVTOT=IVTOT+NTET
      DO 33 I=0,LMAX
   33 S(0,I) =0.D0
      DO 34 M=1,LMAX
      DO 34 I=0,LMAX-M
      S(-M,I)=0.D0
   34 S( M,I)=0.D0
      GO TO 13
   31 IF(NEWSCH(IFACE).EQ.1) GO TO 35
      IVTOT=IVTOT+NTET
      GO TO 32
   35 ARG1=R0(IFACE)/R
      RDOWN=SQRT((R-R0(IFACE))*(R+R0(IFACE)))
      DO 18 I=0,LMAX
   18 S(0,I) =0.D0
      DO 4 M=1,LMAX
      DO 4 I=0,LMAX-M
      S(-M,I)=0.D0
    4 S( M,I)=0.D0
C.......................................................................
C     L O O P     O V E R     T E T R A H E D R A
C.......................................................................
      DO 14 ITET=1,NTET
      IVTOT=IVTOT+1
      IF(R.LE.RD(IVTOT))      T H E N
      CALL PINTG(FA(IVTOT),FB(IVTOT),DLT,S1,LMAX,ISIGNU(IVTOT),
     *           ARG1,FD(IVTOT),0)
      DO 22 I=0,LMAX
   22 S(0,I)=S(0,I)+S1(0,I)
      DO 10 M=1,LMAX
      DO 10 I=0,LMAX-M
      S(-M,I)=S(-M,I)+S1(-M,I)
   10 S( M,I)=S( M,I)+S1( M,I)
                              E L S E
      RAP =RUPSQ(IVTOT)/RDOWN
      ARG2=RUPSQ(IVTOT)/R0(IFACE)
      FK=FD(IVTOT)-ACOS(RAP)
      FL=FD(IVTOT)+ACOS(RAP)
C CRAY AMAX1
      FK=DMAX1(FA(IVTOT),FK)
      FL=DMAX1(FA(IVTOT),FL)
      FK=DMIN1(FB(IVTOT),FK)
      FL=DMIN1(FB(IVTOT),FL)
      CALL PINTG(FA(IVTOT),FK,DLT,S1,LMAX,ISIGNU(IVTOT),
     *           ARG1,FD(IVTOT),0)
      CALL PINTG(FK       ,FL,DLT,S2,LMAX,ISIGNU(IVTOT),
     *           ARG2,FD(IVTOT),1)
      CALL PINTG(FL,FB(IVTOT),DLT,S3,LMAX,ISIGNU(IVTOT),
     *           ARG1,FD(IVTOT),0)
      DO 23 I=0,LMAX
   23 S(0,I)=S(0,I)+S1(0,I)+S2(0,I)+S3(0,I)
      DO 20 M=1,LMAX
      DO 20 I=0,LMAX-M
      S(-M,I)=S(-M,I)+S1(-M,I)+S2(-M,I)+S3(-M,I)
   20 S( M,I)=S( M,I)+S1( M,I)+S2( M,I)+S3( M,I)
                              E N D   I F
   14 CONTINUE
C.......................................................................
C     I N T E G R A L   E X P A N S I O N        B A C K - R O T A T I O
C.......................................................................
   32 CONTINUE
      IB=0
      IC=0
      ICE=0
      READ(ITEMP) (DMATL(ISU),ISU=1,ISUM)
      ISU=0
      DO 5 L=0,LMAX
      IB=IB+L+1
      DO 6 MP=L,1,-1
      SUM(MP,1)=0.D0
      SUM(MP,2)=0.D0
      ICE=ICE+1
      K0=(L+MP+1)/2
      DO 7 K=L,K0,-1
      IS=2*K-L-MP
      IC=IC+1
C CRAY FLOAT
      SUM(MP,2)=SUM(MP,2)+CL(IC)*S(-MP,IS)
    7 SUM(MP,1)=SUM(MP,1)+CL(IC)*S( MP,IS)
      SUM(MP,2)=SUM(MP,2)*C(ICE)
      SUM(MP,1)=SUM(MP,1)*C(ICE)
    6 CONTINUE
      SUM(0,1)=0.D0
      ICE=ICE+1
      K0=(L+1)/2
      DO 24 K=L,K0,-1
      IS=2*K-L
      IC=IC+1
C CRAY FLOAT
   24 SUM(0,1)=SUM(0,1)+CL(IC)*S(0,IS)
      SUM(0,1)=SUM(0,1)*C(ICE)
      IMAX=1
      M=0
    8 CONTINUE
      DO 11 I=1,IMAX
      MO=(3-2*I)*M
      IBM=IB+MO
      LOFM(IBM)=L
      MOFM(IBM)=MO
      IPMAX=1
      MP=0
   16 CONTINUE
      DO 17 IP=1,IPMAX
      ISU=ISU+1
      B(IBM)=B(IBM)+SUM(MP,IP)*DMATL(ISU)
   17 CONTINUE
      IPMAX=2
      MP=MP+1
      IF(MP.LE.L) GO TO 16
   11 CONTINUE
      IMAX=2
      M=M+1
      IF(M.LE.L) GO TO 8
      IB=IB+L
    5 CONTINUE
      REWIND ITEMP
   13 CONTINUE
C.......................................................................
C     D E F I N E S   A N D    S A V E S   S H A P E    F U N C T I O N
C.......................................................................
      B(1)=FPISQ-B(1)/FPISQ
      DO 15 IBM=2,IBMAX
      B(IBM)=-B(IBM)/FPISQ
   15 CONTINUE
      DO 25 IBM=1,IBMAX
C     write(6,*) ibm,b(ibm)
      IF(ABS(B(IBM)).GT.1.D-6) ISW(IBM)=1
      IREC=(IBM-1)*MESHN+N
      WRITE(11,REC=IREC) B(IBM)
   25 CONTINUE
   12 CONTINUE
      NFUN=0
      DO 36 IBM=1,IBMAX
      IF(ISW(IBM).EQ.1)  NFUN=NFUN+1
   36 CONTINUE
C THIS IS FOR DIFFERENT SHELLS...
        NCELL=1
        SCALE=1.d0
C
        WRITE(9,106) NCELL
        WRITE(9,111) SCALE
        WRITE(9,106) NPAN-1,MESHN
        WRITE(9,106) (NM(IPAN),IPAN=1,NPAN-1)
        WRITE(9,111) (XRN(N),DRN(N),N=1,MESHN)
        WRITE(9,106) NFUN
        NCELL_S = NCELL
        SCALE_S = SCALE
        NPAN_S = NPAN-1
        MESHN_S = MESHN
        DO IPAN=1,NPAN-1
        NM_S(IPAN) = NM(IPAN)
        END DO
        DO N=1,MESHN
        XRN_S(N) = XRN(N)
        DRN_S(N) = DRN(N)
        END DO
        NFUN_S = NFUN 
        
C     WRITE(9,106) MESHN,NFUN
C     WRITE(9,108) (XRN(N),N=1,MESHN)
      ICOUNT = 0 
      DO 28 IBM=1,IBMAX
C here is;l;lsd

      IF(ISW(IBM).EQ.0) GO TO 28
 
      ICOUNT = ICOUNT + 1
!      WRITE(6,109) LOFM(IBM),MOFM(IBM)
        LM0=LOFM(IBM)*LOFM(IBM)+LOFM(IBM)+MOFM(IBM)+1
        WRITE(9,106)LM0
        LMIFUN_S(ICOUNT) = LM0
      DO 29 N=1,MESHN
      IREC=(IBM-1)*MESHN+N
      READ(11,REC=IREC) BB(N)
      THETAS_S(N,ICOUNT) = BB(N)
   29 CONTINUE
        
        WRITE(9,111)(BB(N),N=1,MESHN)
!      WRITE(6,108) (BB(N),N=1,MESHN)
   28 CONTINUE
      CLOSE(11)
      RETURN
C     STOP
  200 WRITE(6,107) ISUM,ISUMD,LMAX,LMAXD1
      STOP
  100 FORMAT(4F16.8)
  101 FORMAT(I10,4F10.4,I10)
  102 FORMAT(//15X,'FA/PI',5X,'FB/PI',5X,'FD/PI',6X,'RD',8X,'ISIGNU'/)
  103 FORMAT(4F16.8,2I5)
  104 FORMAT(3I5,F10.5)
  105 FORMAT(/10X,I3,'-TH PYRAMID SUBDIVIDED IN ',I3,' TETRAHEDRA')
  106 FORMAT(16I5)
  107 FORMAT(23X,'FROM MAIN : ISUM=',I7,'  GREATER THAN DIMENSIONED',I7/
     *       23X,'       OR   LMAX=',I7,'  GREATER THAN DIMENSIONED',I7)
  108 FORMAT(5D14.7)
  109 FORMAT(/3X,'ANGULAR MOMENTUM : (',I3,',',I4,')'/3X,29('-')/)
  110 FORMAT(11X,'BUT IS IDENTICAL TO A PREVIOUS ONE.')
  111 FORMAT(4D20.12)
      END
      SUBROUTINE ROTATE(V,VZ,IFACE,NVERT)
      implicit none
C-----------------------------------------------------------------------
C     THIS ROUTINE PERFORMS THE ROTATION OF NVERT VECTORS THROUGH THE
C     EULER ANGLES: ALPHA(IFACE),BETA(IFACE),GAMMA(IFACE).
C     V (I,IVERT) : INPUT   VECTORS
C     VZ(I,IVERT) : ROTATED VECTORS
C-----------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'
!     INTEGER NFACED,NVERTD
!     PARAMETER(NFACED=200,NVERTD=250)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   IFACE,NVERT
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 V(3,NVERTD),VZ(3,NVERTD)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I,J,IVERT
      REAL*8    SA,SB,SG,CA,CB,CG
C
C     .. LOCAL ARRAYS ..
C
      REAL*8 A(3,3)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC COS,SIN
C
C     .. SCALARS IN COMMON ..
C
      REAL*8 PI
C
C     .. ARRAYS IN COMMON ..
C
      REAL*8 ALPHA(NFACED),BETA(NFACED),GAMMA(NFACED)
C
C     .. COMMON BLOCKS ..
C
       COMMON /ANGLES/ PI,ALPHA,BETA,GAMMA
C-----------------------------------------------------------------------
      CA=DCOS(ALPHA(IFACE))
      SA=DSIN(ALPHA(IFACE))
      CB=DCOS(BETA(IFACE))
      SB=DSIN(BETA(IFACE))
      CG=DCOS(GAMMA(IFACE))
      SG=DSIN(GAMMA(IFACE))
      A(1,1)=CA*CB*CG-SA*SG
      A(2,1)=SA*CB*CG+CA*SG
      A(3,1)=-SB*CG
      A(1,2)=-CA*CB*SG-SA*CG
      A(2,2)=-SA*CB*SG+CA*CG
      A(3,2)=SB*SG
      A(1,3)=CA*SB
      A(2,3)=SA*SB
      A(3,3)=CB
      DO 1 IVERT=1,NVERT
      DO 2 I=1,3
      VZ(I,IVERT)=0D0
      DO 2 J=1,3
    2 VZ(I,IVERT)=VZ(I,IVERT)+A(I,J)*V(J,IVERT)
    1 CONTINUE
      RETURN
      END
      SUBROUTINE EULER(Z,XX,IFACE,TOLEULER)
      implicit none
C-----------------------------------------------------------------------
C     GIVEN TWO DISTINCT POINTS (Z(1),Z(2),Z(3)) AND (XX(1),XX(2),XX(3))
C     THIS ROUTINE DEFINES  A LOCAL COORDINATE  SYSTEM WITH THE  Z- AXIS
C     PASSING THROUGH (Z(1),Z(2),Z(3))  AND THE X- AXIS PARALLEL TO  THE
C     VECTOR : (XX(1)-Z(1),XX(2)-Z(2),XX(3)-Z(3)).
C     THE EULER ANGLES ROTATING THIS LOCAL COORDINATE SYSTEM BACK TO THE
C     ORIGINAL FRAME OF REFERENCE ARE CALCULATED  AND STORED  IN COMMON.
C-----------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'
!     INTEGER NFACED
!     PARAMETER(NFACED=200)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   IFACE
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8    XX(3),Z(3)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I
      REAL*8    RX,RZ,S,P,RZP,SA,CA,SG,CG
      REAL*8    TOLEULER   ! introduced by Phivos (05.2008) to account for inaccuracies.
      ! Earlier, 1.D-5 was hard-coded at the places in this subr. where TOL is used
C     DATA TOLEULER /1.D-10/
C
C     .. LOCAL ARRAYS ..
C
      REAL*8    X(3),Y(3)
C
C     .. SCALARS IN COMMON ..
C
      REAL*8    PI
C
C     .. ARRAYS IN COMMON ..
C
      REAL*8    ALPHA(NFACED),BETA(NFACED),GAMMA(NFACED)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC SQRT,ACOS,ABS,ATAN2
C
C     .. COMMON BLOCKS ..
C
      COMMON /ANGLES/ PI,ALPHA,BETA,GAMMA
C-----------------------------------------------------------------------
      IF(IFACE.GT.NFACED) GO TO 20
      DO 4 I=1,3
    4 X(I)=XX(I)-Z(I)
      RX=DSQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
      RZ=DSQRT(Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3))
      IF (RX.LT.1D-6 .OR. RZ.LT.1D-6)  GO TO 30
      S=X(1)*Z(1)+X(2)*Z(2)+X(3)*Z(3)
      IF (S.GT.1D-6)                    GO TO 30
      P=DSQRT(Z(1)*Z(1)+Z(2)*Z(2))
      DO 1 I=1,3
      X(I)=X(I)/RX
    1 Z(I)=Z(I)/RZ
      ALPHA(IFACE)=0D0
      BETA(IFACE) =DACOS(Z(3))
      IF (P .LT. TOLEULER) GO TO 10
      RZP=RZ/P
      Y(1)=Z(2)*X(3)-Z(3)*X(2)
      Y(2)=Z(3)*X(1)-Z(1)*X(3)
      Y(3)=Z(1)*X(2)-Z(2)*X(1)
      SA=Y(3)*RZP
      CA=X(3)*RZP
      IF(DABS(SA).LT.TOLEULER .AND. DABS(CA+1.D0).LT.TOLEULER)  T H E N
      ALPHA(IFACE)=PI
                                                      E L S E
      ALPHA(IFACE)=2D0*DATAN2(SA,CA+1D0)
                                                      E N D    I F
      SG= Z(2)*RZP
      CG=-Z(1)*RZP
      IF(DABS(SG).LT.TOLEULER .AND. DABS(CG+1D0).LT.TOLEULER)  T H E N
      GAMMA(IFACE)=PI
                                                      E L S E
      GAMMA(IFACE)=2D0*DATAN2(SG,CG+1D0)
                                                      E N D    I F
      DO 2 I=1,3
    2 Z(I)=Z(I)*RZ
      RETURN
   10 SG=-Z(3)*X(2)
      CG= Z(3)*X(1)
      IF(DABS(SG).LT.TOLEULER .AND. DABS(CG+1D0).LT.TOLEULER)  T H E N
      GAMMA(IFACE)=PI
                                                      E L S E
      GAMMA(IFACE)=2D0*DATAN2(SG,CG+1D0)
                                                      E N D    I F
      DO 3 I=1,3
    3 Z(I)=Z(I)*RZ
      RETURN
   20 WRITE(6,100) IFACE,NFACED
      STOP
   30 WRITE(6,101) (X(I),I=1,3),(Z(I),I=1,3)
      STOP
  100 FORMAT(//13X,'NUMBER OF FACES:',I5,' GREATER THAN DIMENSIONED',I5)
  101 FORMAT(/13X,'FROM EULER,ILLEGAL VECTORS:'/13X,2(' (',3E13.6,' )'))
      END
      SUBROUTINE PERP(R0,R1,R2,RD,TOLVDIST,INSIDE)
      implicit none
C-----------------------------------------------------------------------
C     GIVEN  TWO  DISTINCT  POINTS   R1 , R2, THIS  ROUTINE CALCULATES
C     THE COORDINATES  OF THE FOOT OF  THE  PERPENDICULAR FROM A POINT
C     R0 TO THE LINE JOINING R1   AND  R2. THE LOGICAL VARIABLE INSIDE
C     GIVES THE ADDITIONAL INFORMATION WHETHER THE FOOT OF THE PERPEN-
C     DICULAR LIES WITHIN THE SEGMENT OR NOT.
C-----------------------------------------------------------------------
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 R0(3),R1(3),R2(3),RD(3)
      REAL*8 TOLVDIST
C
C     .. LOGICAL ARGUMENTS ..
C
      LOGICAL INSIDE
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I
      REAL*8    DX,DY,DZ,S,D,DA,DB,DC,D1,D2,CO
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC DMAX1
C---------------------------------------------------------------------
      DX=R2(1)-R1(1)
      DY=R2(2)-R1(2)
      DZ=R2(3)-R1(3)
      S=R0(1)*DX+R0(2)*DY+R0(3)*DZ
      D=DX*DX+DY*DY+DZ*DZ
      IF(SQRT(D).LT.TOLVDIST) GO TO 100
      DA=S*DX+DY*(R1(1)*R2(2)-R1(2)*R2(1))+DZ*(R1(1)*R2(3)-R1(3)*R2(1))
      DB=S*DY+DZ*(R1(2)*R2(3)-R1(3)*R2(2))+DX*(R1(2)*R2(1)-R1(1)*R2(2))
      DC=S*DZ+DX*(R1(3)*R2(1)-R1(1)*R2(3))+DY*(R1(3)*R2(2)-R1(2)*R2(3))
      RD(1)=DA/D
      RD(2)=DB/D
      RD(3)=DC/D
      D1=(RD(1)-R1(1))**2+(RD(2)-R1(2))**2+(RD(3)-R1(3))**2
      D2=(RD(1)-R2(1))**2+(RD(2)-R2(2))**2+(RD(3)-R2(3))**2
C CRAY AMAX1
      CO=D-DMAX1(D1,D2)
      INSIDE=.FALSE.
      IF ( CO.GT.TOLVDIST)  INSIDE=.TRUE.
      RETURN
  100 WRITE(6,200) (R1(I),I=1,3),(R2(I),I=1,3)
      STOP
  200 FORMAT(///33X,'FROM PERP:   IDENTICAL POINTS'/33X,2('(',3E14.6,')'
     *,3X))
      END
      SUBROUTINE CRIT(IFACE,NVERT,V,Z,IPAN,IVTOT,TOLEULER,TOLVDIST,CRT)
      implicit none
C-----------------------------------------------------------------------
C     THIS ROUTINE CALCULATES THE CRITICAL POINTS 'CRT' OF THE SHAPE
C     FUNCTIONS DUE TO THE FACE: Z(1)*X + Z(2)*Y + Z(3)*Z = 1
C     THE FACE IS ROTATED THROUGH THE APPROPRIATE EULER ANGLES TO BE
C     PERPENDICULAR TO THE Z-AXIS. A FURTHER SUBDIVISION OF THE CEN-
C     TRAL PYRAMID INTO ELEMENTARY TETRAHEDRA IS PERFORMED. THE  NE-
C     CESSARY QUANTITIES FOR THE CALCULATION ARE STORED IN COMMON.
C-----------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'
      INTEGER NVTOTD
!     PARAMETER (NVERTD=250,NFACED=200)
      PARAMETER (NVTOTD=NFACED*NVERTD)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   IFACE,NVERT,IPAN,IVTOT
      REAL*8 TOLEULER,TOLVDIST
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 V(3,NVERTD),Z(3),CRT(NPAND)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I,IX,ICORN,IVERT,INEW,IP,IVERTP,IVERT1,IBACK
      REAL*8    ARG,A1,A2,A3,CF1,CF2,CF3,CO,CRRT,DD,DOWN,D1,D2
      REAL*8    FF,F1,F2,OMEGA,RDD,S,SF1,SF2,SF3,UP,XJ,YJ,ZMOD2
      REAL*8    ZVMOD
C
C     .. LOCAL ARRAYS ..
C
      INTEGER   IN(NVERTD)
      REAL*8    VZ(3,NVERTD),RDV(3),ORIGIN(3)
C
C     .. LOCAL LOGICAL ..
C
      LOGICAL INSIDE,TEST
C
C     .. SCALARS IN COMMON ..
C
      REAL*8 PI
C
C     .. ARRAYS IN COMMON ..
C
      INTEGER   ISIGNU(NVTOTD),NTT(NFACED)
      REAL*8    RD(NVTOTD),R0(NFACED),FA(NVTOTD),FB(NVTOTD),FD(NVTOTD)
      REAL*8    ALPHA(NFACED),BETA(NFACED),GAMMA(NFACED)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC ABS,ACOS,DMAX1,DMIN1,ATAN2,SIGN,SQRT
C
C     .. EXTERNAL ROUTINES ..
C
      EXTERNAL EULER,PERP,ROTATE
C
C     .. COMMON BLOCKS ..
C
      COMMON/ANGLES/  PI,ALPHA,BETA,GAMMA
      COMMON/TETRA/ FA,FB,FD,R0,RD,ISIGNU,NTT
C
C     .. DATA STATEMENTS ..
C
      DATA ORIGIN/3*0.D0/
      PI=4.D0*DATAN(1.D0)  ! added 1.3.2012, fivos
C-----------------------------------------------------------------------
      NTT(IFACE)=0
      IF (TEST('verb0   ')) WRITE(6,203) IFACE,(Z(I),I=1,3)
      ZMOD2=Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3)
      IF(ZMOD2.LE.1D-6) GO TO 100
      S=2D0*PI
      DO 1 I=1,3
    1 Z(I)=Z(I)/ZMOD2
      IX=1
      ZVMOD=DSQRT((V(1,1)-Z(1))**2+(V(2,1)-Z(2))**2+(V(3,1)-Z(3))**2)
      IF(ZVMOD.LT.1D-6)  IX=2 
      CALL EULER(Z,V(1,IX),IFACE,TOLEULER)
      IF (TEST('verb0   '))
     &     WRITE(6,204) ALPHA(IFACE)/PI,BETA(IFACE)/PI,GAMMA(IFACE)/PI
      CALL ROTATE(V,VZ,IFACE,NVERT)
      R0(IFACE)=1D0/DSQRT(ZMOD2)
      ICORN=0
      IF (TEST('verb0   ')) WRITE(6,207)
      DO 2 IVERT =1,NVERT
      IF(DABS(R0(IFACE)-VZ(3,IVERT)).GT.1D-6) GO TO 101
C.......................................................................
C     D I S T A N C E S   O F   V E R T I C E S   F R O M   C E N T E R
C.......................................................................
      CRRT=DSQRT(VZ(1,IVERT)**2+VZ(2,IVERT)**2+VZ(3,IVERT)**2)
      INEW=1
      DO 3 IP=1,IPAN
      IF(DABS(CRRT-CRT(IP)).LT.1D-6) INEW=0
    3 CONTINUE
      IF(INEW.EQ.1)          T H E N
      IPAN=IPAN+1
      IF(IPAN.GT.NPAND) GO TO 102
      CRT(IPAN)=CRRT
                             E N D   I F
      IVERTP=IVERT+1
      IF(IVERT.EQ.NVERT) IVERTP=1
C.......................................................................
C     D I S T A N C E S   O F   E D G E S   F R O M   C E N T E R
C.......................................................................
      CALL PERP(ORIGIN,VZ(1,IVERT),VZ(1,IVERTP),RDV,TOLVDIST,INSIDE)
      RDD=DSQRT(RDV(1)*RDV(1)+RDV(2)*RDV(2)+RDV(3)*RDV(3))
      IF(INSIDE)             T H E N
      INEW=1
      DO 4 IP=1,IPAN
      IF(DABS(RDD-CRT(IP)).LT.1D-4) INEW=0
    4 CONTINUE
      IF(INEW.EQ.1)          T H E N
      IPAN=IPAN+1
      IF(IPAN.GT.NPAND) GO TO 102
      CRT(IPAN)=RDD
                             E N D   I F
                             E N D    I F
      A1=DSQRT(VZ(1,IVERT )*VZ(1,IVERT )+VZ(2,IVERT )*VZ(2,IVERT ))
      A2=DSQRT(VZ(1,IVERTP)*VZ(1,IVERTP)+VZ(2,IVERTP)*VZ(2,IVERTP))
      DOWN=A1*A2
      UP=VZ(1,IVERT)*VZ(1,IVERTP)+VZ(2,IVERT)*VZ(2,IVERTP)
      IF(DOWN.GT.1D-6)      T H E N
      ARG=UP/DOWN
      IF(DABS(ARG).GE.1D0) ARG=SIGN(1D0,ARG)
      OMEGA=DACOS(ARG)
      S=S-OMEGA
      IF(DABS(OMEGA-PI).GT.1D-6)                 T H E N
C.......................................................................
C     S U B D I V I S I O N    I N T O    T E T R A H E D R A
C.......................................................................
      NTT(IFACE)=NTT(IFACE)+1
      IVTOT=IVTOT+1
      IF (TEST('verb0   ')) WRITE(6,205) IVTOT,IVERT,(VZ(I,IVERT),I=1,3)
      IF (TEST('verb0   ')) WRITE(6,206)     IVERTP,(VZ(I,IVERTP),I=1,3)
      A3=DSQRT(RDV(1)*RDV(1)+RDV(2)*RDV(2))
      RD    (IVTOT)=RDD
      ISIGNU(IVTOT)=1
      CF1=VZ(1,IVERT )/A1
      CF2=VZ(1,IVERTP)/A2
      SF1=VZ(2,IVERT )/A1
      SF2=VZ(2,IVERTP)/A2
      CF3=RDV(1)/A3
      SF3=RDV(2)/A3
      IF(DABS(SF1).LT.TOLEULER .AND. DABS(CF1+1D0).LT.TOLEULER)  T H E N
      F1=PI
                                                           E L S E
      F1=2D0*DATAN2(SF1,CF1+1D0)
                                                           E N D    I F
      IF(DABS(SF2).LT.TOLEULER .AND. DABS(CF2+1D0).LT.TOLEULER)  T H E N
      F2=PI
                                                           E L S E
      F2=2D0*DATAN2(SF2,CF2+1D0)
                                                           E N D    I F
      IF(DABS(SF3).LT.TOLEULER .AND. DABS(CF3+1D0).LT.TOLEULER)  T H E N
      FD(IVTOT)=PI
                                                           E L S E
      FD(IVTOT)=2D0*DATAN2(SF3,CF3+1D0)
                                                           E N D    I F
C CRAY AMIN1
      FA(IVTOT)=DMIN1(F1,F2)
      FB(IVTOT)=DMAX1(F1,F2)
      IF((FB(IVTOT)-FA(IVTOT)).GT.PI)                 T H E N
      FF=FA(IVTOT)+2D0*PI
      FA(IVTOT)=FB(IVTOT)
      FB(IVTOT)=FF
                                                      E N D   I F
      IF((FA(IVTOT)-FD(IVTOT)).GT.PI)  FD(IVTOT)= 2D0*PI+FD(IVTOT)
      IF((FD(IVTOT)-FA(IVTOT)).GT.PI)  FD(IVTOT)=-2D0*PI+FD(IVTOT)
                                                 E N D    I F
                             E L S E
      ICORN=1
                             E N D    I F
    2 CONTINUE
C.......................................................................
C     F O O T   O F   T H E    P E R P E N D I C U L A R   TO    T H E
C     F A C E   O U T S I D E   O R  I N S I D E   T H E   P O L Y G O N
C.......................................................................
      IF(S.LT.1D-06.OR.ICORN.EQ.1)               T H E N
      INEW=1
      DO 5 IP=1,IPAN
      IF(DABS(R0(IFACE)-CRT(IP)).LT.1D-4) INEW=0
    5 CONTINUE
      IF(INEW.EQ.1)          T H E N
      IPAN=IPAN+1
      IF(IPAN.GT.NPAND) GO TO 102
      CRT(IPAN)=R0(IFACE)
                             E N D   I F
                                                 E L S E
      DO 6 IVERT1=1,NVERT
      IN(IVERT1)=0
      DO 7 IVERT=1,NVERT
      IVERTP=IVERT+1
      IF(IVERT.EQ.NVERT) IVERTP=1
      IF(IVERT.EQ.IVERT1.OR.IVERTP.EQ.IVERT1)    GO TO 7
      DOWN=VZ(2,IVERT1)*(VZ(1,IVERTP)-VZ(1,IVERT))
     *    -VZ(1,IVERT1)*(VZ(2,IVERTP)-VZ(2,IVERT))
      IF(DABS(DOWN).LE.1D-06)                     GO TO 7
      UP  =VZ(1,IVERT1)*(VZ(2,IVERT)*(VZ(1,IVERTP)+VZ(1,IVERT))
     *    -              VZ(1,IVERT)*(VZ(2,IVERTP)+VZ(2,IVERT)))
      XJ=UP/DOWN
      YJ=XJ*VZ(2,IVERT1)/VZ(1,IVERT1)
      DD=(VZ(1,IVERTP)-VZ(1,IVERT))**2+(VZ(2,IVERTP)-VZ(2,IVERT))**2
      D1=(XJ-VZ(1,IVERT ))**2+(YJ-VZ(2,IVERT ))**2
      D2=(XJ-VZ(1,IVERTP))**2+(YJ-VZ(2,IVERTP))**2
C CRAY AMAX1
      CO=DD-DMAX1(D1,D2)
      IF(CO.GT.1D-06)        T H E N
      IN(IVERT1)=1
      GO TO 6
                             E N D   I F
    7 CONTINUE
    6 CONTINUE
      IBACK=IVTOT-NVERT
      DO 8 IVERT=1,NVERT
      IBACK=IBACK+1
      IVERTP=IVERT+1
      IF(IVERT.EQ.NVERT) IVERTP=1
      IF(IN(IVERT).EQ.0.AND.IN(IVERTP).EQ.0)     ISIGNU(IBACK)=-1
    8 CONTINUE
                                                 E N D   I F
      RETURN
  100 WRITE(6,200) IFACE,(Z(I),I=1,3)
      STOP
  101 WRITE(6,201) IFACE,R0(IFACE),((VZ(I,IVERT),I=1,3),IVERT=1,NVERT)
      STOP
  102 WRITE(6,202) IPAN,NPAND
      STOP
  200 FORMAT(//13X,'FATAL ERROR FROM CRIT: THE',I3,'-TH FACE OF THE POLY
     *HEDRON PASSES THROUGH THE CENTER'/13X,'(',3E14.7,' )')
  201 FORMAT(//13X,'FATAL ERROR FROM CRIT: THE VERTICES OF THE',I3,'-TH
     *ROTATED POLYGON DO NOT LIE ON THE PLANE:',E13.6,' *Z = 1'/30(/13X,
     *3E13.6))
  202 FORMAT(//13X,'ERROR FROM CRIT: NUMBER OF PANNELS=',I5,' GREATER TH
     *AN DIMENSIONED=',I5)
  203 FORMAT(//80('*')/3X,'FACE:',I3,' EQUATION:',F10.4,'*X +',F10.4,
     *'*Y +',F10.4,'*Z  =  1')
  204 FORMAT(3X,'ROTATION ANGLES  :',3(F10.4,4X)/)
  205 FORMAT(I5,'       VZ(',I2,')  =  (',3F10.4,' )')
  206 FORMAT(5X,'       VZ(',I2,')  =  (',3F10.4,' )')
  207 FORMAT(/'TETRAHEDRON',14X,'COORDINATES'/11('*'),14X,11('*')/)
      END
      SUBROUTINE MESH(CRT,NPAN,NM,XRN,DRN,MESHN,NPOI,KEYPAN,NMIN)
      implicit none
C-----------------------------------------------------------------------
C     THIS ROUTINE DEFINES A UNIQUE SUITABLE RADIAL MESH 'XRN,DRN' OF
C     'MESHN' POINTS,DISTRIBUTED INTO 'NPAN' PANNELS  DEFINED  BY THE
C     CRITICAL POINTS 'CRT'
C-----------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'
c      INTEGER NPAND,MESHND
c      PARAMETER (NPAND=80,MESHND=1000)
       INTEGER MESHND
       PARAMETER (MESHND=IRID)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   NPAN,MESHN,NPOI,KEYPAN,NMIN
C
C     .. ARRAY ARGUMENTS ..
C
      INTEGER   NM(NPAND)
      REAL*8    CRT(NPAND),XRN(MESHND),DRN(MESHND)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   IORD,IPAN,N1,N2,K
      REAL*8    C,D
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC DFLOAT
C-----------------------------------------------------------------------
      DO 10 IORD=1,NPAN
      C=CRT(IORD)
      DO 20 IPAN=NPAN,IORD,-1
      IF(CRT(IPAN).GT.C)   GO TO 20
      C=CRT(IPAN)
      CRT(IPAN)=CRT(IORD)
      CRT(IORD)=C
   20 CONTINUE
   10 CONTINUE
C
C     CALCULATE AN APPROPRIATE MESH
C
      IF (KEYPAN.EQ.0) THEN
      CALL MESH0(CRT,NM,NPAN,NPOI,NMIN)
      END IF
C
      WRITE(6,103)
      N2=0
      DO 50 IPAN = 1,NPAN-1
      WRITE(6,104) IPAN,CRT(IPAN),CRT(IPAN+1),NM(IPAN)
      N1 = N2 + 1
      N2 = N2 + NM(IPAN)
      IF (MESHND.GE.N2)      T  H  E  N
C CRAY FLOAT
      C = (CRT(IPAN+1)-CRT(IPAN))/DFLOAT(N2-N1)
      D = CRT(IPAN) - C*DFLOAT(N1)
      DO 60 K = N1,N2
      XRN(K) = C*DFLOAT(K) + D
      DRN(K) = C
   60 CONTINUE
                             E  L  S  E
      GO TO 70
                             E  N  D    I  F
   50 CONTINUE
      WRITE(6,105)
      MESHN = N2
c      WRITE(6,101)(K,DRN(K),XRN(K),K=1,MESHN)
      RETURN
   70 WRITE (6,102) MESHND
      STOP
  101 FORMAT (/'    NEW MESH  K,DRN,XRN'/(1H ,I3,2F12.7))
  102 FORMAT ('   *** FROM MESH  :    NXR=',I4,' IS TOO SMALL')
  103 FORMAT(/50('-')/'I',13X,'SUITABLE RADIAL MESH',15X,'I'/'I',13X,
     *20('*'),15X,'I'/'I',3X,'IPAN',7X,'FROM',7X,'TO',13X,'POINTS  I'/
     *'I',48X,'I')
  104 FORMAT('I',2X,I5,2E24.16,I10,'   I')
  105 FORMAT(50('-'))
      END
      SUBROUTINE PINTG(X1,X2,DLT,S,LMAX,ISI,ARG,FD,ITYPE)
      implicit none
C-----------------------------------------------------------------------
C     THIS ROUTINE  ACCOMPLISHES THE  FI-INTEGRATION  OF REAL  SPHERICAL
C     HARMONICS BY THE REPEATED SIMPSON'S METHOD , OR ANALYTICALLY ACCOR
C     DING TO THE VALUE OF ITYPE. THE OBTAINED RESULTS HAVE TO BE MULTI-
C     PLIED BY THE APPROPRIATE EXPANSION COEFFICIENTS.
C-----------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'
c      INTEGER LMAXD,NDIM
c      PARAMETER (LMAXD=25,NDIM=1000)
      INTEGER NDIM
      PARAMETER (NDIM=1000)
C
C     .. SCALAR ARGUMENTS ..
C
      REAL*8 X1,X2,DLT,ARG,FD
      INTEGER   LMAX,ISI,ITYPE
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 S(-LMAXD1:LMAXD1,0:LMAXD1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I,M,N,K
      REAL*8    X,THETA,W
C
C     .. LOCAL ARRAYS ..
C
      REAL*8    XX(NDIM),WW(NDIM)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC DFLOAT,ACOS,ATAN,COS,IABS
C
C     .. EXTERNAL ROUTINE ..
C
      EXTERNAL RECUR,RECUR0,GAULEG
C-----------------------------------------------------------------------
      IF(LMAX.LE.LMAXD1) GO TO 1
      WRITE(6,200)LMAX,LMAXD1
  200 FORMAT(3X,'FROM PINTG: LMAX=',I4,' GREATER THAN DIMENSIONED',I4)
      STOP
    1 CONTINUE
      DO 6 I=0,LMAX
    6 S( 0,I)=0.D0
      DO 5 M=1,LMAX
      DO 5 I=0,LMAX-M
      S(-M,I)=0.D0
    5 S( M,I)=0.D0
      IF(ITYPE.NE.0)      GO TO 10
      THETA=ACOS(ARG)
      CALL RECUR0(LMAX,X1,THETA,-DFLOAT(ISI),S)
      CALL RECUR0(LMAX,X2,THETA, DFLOAT(ISI),S)
      RETURN
C                         E N D    I F
   10 CONTINUE
      N=(X2-X1)/DLT+3
      IF(N.GT.NDIM) STOP 'INCREASE NDIM'
      CALL GAULEG(X1,X2,XX,WW,N)
      DO 2 K=1,N
      X=XX(K)
      W=DFLOAT(ISI)*WW(K)
      THETA=ATAN(ARG/COS(X-FD))
      CALL RECUR(LMAX,X,THETA,W,S)
    2 CONTINUE
      RETURN
      END
      SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT NONE
C     ----------------------------------------------------------------
C     GINEN THE LOWER AND UPPER LIMITS OF INTEGRATION  X1 AND X2, AND
C     GIVEN N, THIS SUBROUTINE RETURNS THE  ARRAYS X(1:N) AND  W(1:N)
C     OF LENGTH N, CONTAINING THE ABSCISSAS AND WEIGHTS OF THE  GAUSS
C     LEGENDRE N-POINT QUADRATURE FORMULA (NUMERICAL RECIPES,2ND ED.).
C     ----------------------------------------------------------------
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   N
      REAL*8    X1,X2
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8    X(1),W(1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I,J,M
      REAL*8    P1,P2,P3,PP,XL,XM,Z,Z1,PI314
C     ----------------------------------------------------------------
      PI314 = 4.D0*DATAN(1.D0)

      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
      Z=DCOS(PI314*(I-.25D0)/(N+.5D0))
    1 CONTINUE
      P1=1.D0
      P2=0.D0
      DO 11 J=1,N
      P3=P2
      P2=P1
      P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
   11 CONTINUE
      PP=N*(Z*P1-P2)/(Z*Z-1.D0)
      Z1=Z
      Z=Z1-P1/PP
      IF(ABS(Z-Z1).GT.3.D-14) GO TO 1
      X(I)=XM-XL*Z
      X(N+1-I)=XM+XL*Z
      W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
      W(N+1-I)=W(I)
   12 CONTINUE
      RETURN
      END
      SUBROUTINE RECUR(LMAX,X,THETA,FAC,S)
      implicit none
C-----------------------------------------------------------------------
C     THIS ROUTINE IS USED TO PERFORM THE FI-INTEGRATION OF REAL SPHE-
C     RICAL HARMONICS .THE THETA-INTEGRATION IS PERFORMED ANALYTICALLY
C     USING RECURRENCE RELATIONS.
C-----------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
       include 'inc.geometry'
c      INTEGER LMAXD
c      PARAMETER (LMAXD=25)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   LMAX
      REAL*8 X,THETA,FAC
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 S(-LMAXD1:LMAXD1,0:LMAXD1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   M,I
      REAL*8 OL0,OL,EL0,EL,C1,C2,SS,CC
      REAL*8 C01(LMAXD1),C02(LMAXD1),SSA(LMAXD1+2),CCA(LMAXD1+2)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC COS,SIN,DFLOAT
C-----------------------------------------------------------------------
      SS=SIN(THETA)
      CC=COS(THETA)
      DO 13 I=1,LMAX
      C01(I)=FAC*SIN(DFLOAT(I)*X)
   13 C02(I)=FAC*COS(DFLOAT(I)*X)
      DO 11 I=1,LMAX+2
      SSA(I)=SS**I
      CCA(I)=CC**I
   11 CONTINUE
      OL0=(THETA-SS*CC)/2.D0
      EL0=0D0
      DO 1 M=1,LMAX,2
      OL=OL0
      EL=EL0
      C1=C01(M)
      C2=C02(M)
      DO 2 I=0,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
    2 OL=(DFLOAT(I+1)*OL+(SSA(M+2))*(CCA(I+1)))/DFLOAT(I+M+3)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+3)
      OL0=(DFLOAT(M+2)*OL0-(SSA(M+2))*CC)/DFLOAT(M+3)
    1 CONTINUE
      OL0=SSA(3)/3.D0
      EL0=0D0
      DO 3 M=1,LMAX,2
      OL=OL0
      EL=EL0
      C1=C01(M)
      C2=C02(M)
      DO 4 I=1,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SSA(M+2)*CCA(I+1))/DFLOAT(I+M+3)
    4 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-(SSA(M+2))*CCA(2))/DFLOAT(M+4)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+4)
    3 CONTINUE
      OL0=-CC
      EL0=-1D0
      OL=OL0
      EL=EL0
      C2=FAC
      DO 9 I=0,LMAX,2
      S(0,I)=S(0,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SSA(2)*CCA(I+1))/DFLOAT(I+3)
    9 EL= DFLOAT(I+1)*EL/DFLOAT(I+3)
      OL0=(2.D0*OL0-SSA(2)*CC)/3.D0
      EL0= 2.D0*EL0/3.D0
      DO 5 M=2,LMAX,2
      OL=OL0
      EL=EL0
      C1=C01(M)
      C2=C02(M)
      DO 6 I=0,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SSA(M+2)*CCA(I+1))/DFLOAT(I+M+3)
    6 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-SSA(M+2)*CC)/DFLOAT(M+3)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+3)
    5 CONTINUE
      OL0=-CCA(2)/2D0
      EL0=-0.5D0
      OL=OL0
      EL=EL0
      C2=FAC
      DO 10 I=1,LMAX,2
      S(0,I)=S(0,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SS*SS*CCA(I+1))/DFLOAT(I+3)
   10 EL= DFLOAT(I+1)*EL/DFLOAT(I+3)
      OL0=(2.D0*OL0-SSA(2)*CCA(2))/4.D0
      EL0= 2.D0*EL0/4.D0
      DO 7 M=2,LMAX,2
      OL=OL0
      EL=EL0
      C1=C01(M)
      C2=C02(M)
      DO 8 I=1,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SSA(M+2)*CCA(I+1))/DFLOAT(I+M+3)
    8 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-SSA(M+2)*CCA(2))/DFLOAT(M+4)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+4)
    7 CONTINUE
      RETURN
      END
      SUBROUTINE RECUR0(LMAX,X,THETA,FAC,S)
      implicit none
C-----------------------------------------------------------------------
C     THIS ROUTINE IS USED TO PERFORM  THE  FI - INTEGRATION  OF REAL SP
C     RICAL HARMONICS ANALYTICALLY.  THE  THETA-INTEGRATION  IS   PERFOR
C     ALSO ANALYTICALLY USING RECURRENCE RELATIONS.(THETA IS FI-INDEPEND
C-----------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
       include 'inc.geometry'
c      INTEGER LMAXD
c      PARAMETER (LMAXD=25)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   LMAX
      REAL*8 X,THETA,FAC
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8 S(-LMAXD1:LMAXD1,0:LMAXD1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   M,I
      REAL*8 OL0,OL,EL0,EL,C1,C2,SS,CC
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC COS,SIN,DFLOAT
C-----------------------------------------------------------------------
      SS=SIN(THETA)
      CC=COS(THETA)
      OL0=(THETA-SS*CC)/2D0
      EL0=0D0
      DO 1 M=1,LMAX,2
      OL=OL0
      EL=EL0
      C1=-FAC*COS(DFLOAT(M)*X)/DFLOAT(M)
      C2= FAC*SIN(DFLOAT(M)*X)/DFLOAT(M)
      DO 2 I=0,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
    2 OL=(DFLOAT(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DFLOAT(I+M+3)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+3)
      OL0=(DFLOAT(M+2)*OL0-(SS**(M+2))*CC)/DFLOAT(M+3)
    1 CONTINUE
      OL0=SS*SS*SS/3D0
      EL0=0D0
      DO 3 M=1,LMAX,2
      OL=OL0
      EL=EL0
      C1=-FAC*COS(DFLOAT(M)*X)/DFLOAT(M)
      C2= FAC*SIN(DFLOAT(M)*X)/DFLOAT(M)
      DO 4 I=1,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DFLOAT(I+M+3)
    4 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-(SS**(M+2))*CC*CC)/DFLOAT(M+4)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+4)
    3 CONTINUE
      OL0=-CC
      EL0=-1D0
      OL=OL0
      EL=EL0
      C2=FAC*X
      DO 9 I=0,LMAX,2
      S(0,I)=S(0,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SS*SS*(CC**(I+1)))/DFLOAT(I+3)
    9 EL= DFLOAT(I+1)*EL/DFLOAT(I+3)
      OL0=(2.D0*OL0-SS*SS*CC)/3.D0
      EL0= 2.D0*EL0/3.D0
      DO 5 M=2,LMAX,2
      OL=OL0
      EL=EL0
      C1=-FAC*COS(DFLOAT(M)*X)/DFLOAT(M)
      C2 =FAC*SIN(DFLOAT(M)*X)/DFLOAT(M)
      DO 6 I=0,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DFLOAT(I+M+3)
    6 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-(SS**(M+2))*CC)/DFLOAT(M+3)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+3)
    5 CONTINUE
      OL0=-CC*CC/2D0
      EL0=-0.5D0
      OL=OL0
      EL=EL0
      C2=FAC*X
      DO 10 I=1,LMAX,2
      S(0,I)=S(0,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+SS*SS*(CC**(I+1)))/DFLOAT(I+3)
   10 EL= DFLOAT(I+1)*EL/DFLOAT(I+3)
      OL0=(2.D0*OL0-SS*SS*CC*CC)/4.D0
      EL0= 2.D0*EL0/4.D0
      DO 7 M=2,LMAX,2
      OL=OL0
      EL=EL0
      C1=-FAC*COS(DFLOAT(M)*X)/DFLOAT(M)
      C2= FAC*SIN(DFLOAT(M)*X)/DFLOAT(M)
      DO 8 I=1,LMAX-M,2
      S(-M,I)=S(-M,I)+C1*(OL-EL)
      S( M,I)=S( M,I)+C2*(OL-EL)
      OL=(DFLOAT(I+1)*OL+(SS**(M+2))*(CC**(I+1)))/DFLOAT(I+M+3)
    8 EL= DFLOAT(I+1)*EL/DFLOAT(I+M+3)
      OL0=(DFLOAT(M+2)*OL0-(SS**(M+2))*CC*CC)/DFLOAT(M+4)
      EL0= DFLOAT(M+2)*EL0/DFLOAT(M+4)
    7 CONTINUE
      RETURN
      END
      SUBROUTINE REDUCE(NMBR,IFMX,IFI,IEXP)
      implicit none
C-----------------------------------------------------------------------
C     THIS ROUTINE REDUCES A POSITIVE INTEGER   INPUT NUMBER 'NMBR'
C     TO A PRODUCT  OF  FIRST  NUMBERS 'IFI' , AT POWERS  'IEXP'.
C-----------------------------------------------------------------------
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   NMBR,IFMX
C
C     .. ARRAY  ARGUMENTS ..
C
      INTEGER   IEXP(1),IFI(1)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   I,NMB
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC MOD
C-----------------------------------------------------------------------
      IF(NMBR.LE.0) GO TO 4
      DO 5 I=1,IFMX
    5 IEXP(I)=0
      IF(NMBR.EQ.1) RETURN
      NMB=NMBR
      DO 1 I=1,IFMX
      IEXP(I)=0
    2 IF(MOD(NMB,IFI(I)).NE.0) GO TO 3
      NMB=NMB/IFI(I)
      IEXP(I)=IEXP(I)+1
      GO TO 2
    3 CONTINUE
      IF(NMB.EQ.1) RETURN
    1 CONTINUE
      WRITE(6,100) NMBR
      STOP
    4 WRITE(6,101) NMBR
      STOP
  100 FORMAT(3X,I15,'  CANNOT BE REDUCED IN THE BASIS OF FIRST NUMBERS G
     *IVEN'/20X,'INCREASE THE BASIS OF FIRST NUMBERS')
  101 FORMAT(3X,I15,'  NON POSITIVE NUMBER')
      END
      SUBROUTINE CCOEF(LMAX,CL,COE)
      implicit none
C-----------------------------------------------------------------------
C     THIS ROUTINE CALCULATES THE COEFFICIENTS OF A POLYNOMIAL EXPANSION
C     OF RENORMALIZED LEGENDRE FUNCTIONS IN POWERS OF COSINES.
C     THE POSSIBILITY OF OVERFLOW (HIGH LMAX) IS AVOIDED BY USING FACTO-
C     RIZED FORMS FOR THE NUMBERS.
C-----------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'
      INTEGER ICD,ICED,IFMX,LMA2D
      PARAMETER (ICD=1729,ICED=((LMAXD1+1)*(LMAXD1+2))/2)
      PARAMETER (IFMX=25,LMA2D=LMAXD1/2+1)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   LMAX
C
C     .. ARRAY ARGUMENTS ..
C
      REAL*8    CL(ICD),COE(ICED)
C
C     .. LOCAL SCALARS ..
C
      INTEGER   ICMAX,L,LI,ICE,IC,I1,L2P,M,K,K0,ISI,IRE,IR,IC1,IC2
      INTEGER   LA,LB,IEUPSQ,IEINT,IEMOD
      REAL*8    UP,DOWN,UPSQ
C
C     .. LOCAL ARRAYS ..
C
      INTEGER   IE(IFMX,LMA2D),IED(IFMX),IFI(IFMX)
      INTEGER   L1ST(IFMX),L2ST(IFMX),L1(IFMX),L2(IFMX),JM0(IFMX)
      INTEGER   IEA(IFMX),IEB(IFMX),IL2P(IFMX)
      logical test
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC MOD,SQRT,DFLOAT
C
C     .. EXTERNAL ROUTINES ..
C
      EXTERNAL REDUCE,TEST
C
C     .. DATA STATEMENTS ..
C
      DATA IFI/2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,
     &         61,67,71,73,79,83,89,97/
C-----------------------------------------------------------------------
      ICMAX=0
      DO 9 L=0,LMAX
      LI=L/2+1
    9 ICMAX=ICMAX+(LMAX+1-L)*LI
      IF(LMAX.GT.LMAXD1.OR.ICMAX.GT.ICD) GO TO 100
      if (TEST('verb0   ')) WRITE(6,203) ICMAX
      ICE=0
      IC=1
      L=0
      DO 11 I1=1,IFMX
      L1ST(I1)=0
   11 L2ST(I1)=0
    1 CONTINUE
      L2P=2*L+1
      CALL REDUCE(L2P,IFMX,IFI,IL2P)
      IL2P(1)=IL2P(1)+1
      M=L
      DO 12 I1=1,IFMX
      L1(I1)=L1ST(I1)
      L2(I1)=L2ST(I1)
   12 JM0(I1)=0
    2 CONTINUE
      ICE=ICE+1
      ISI=1
C  THIS IS CHANGED
C     ISI=1-2*MOD(M,2)
C
      K0=(L+M+1)/2
      K=L
      IRE=1
      IC1=IC
      DO 13 I1=1,IFMX
   13 IE(I1,IRE)=JM0(I1)
    3 CONTINUE
      IF((K-1). LT .K0)  GO TO 30
      IRE=IRE+1
      IC=IC+1
      LA=(2*K-L-M)*(2*K-L-M-1)
      LB=2*(2*K-1)*(L-K+1)
      CALL REDUCE(LA,IFMX,IFI,IEA)
      CALL REDUCE(LB,IFMX,IFI,IEB)
      DO 14 I1=1,IFMX
   14 IE(I1,IRE)=IE(I1,IRE-1)+IEA(I1)-IEB(I1)
      K=K-1
      GO TO 3
   30 CONTINUE
      IC2=IC
      DO 4 I1=1,IFMX
      IED(I1)=IE(I1,1)
      DO 5 IR=2,IRE
      IF(IE(I1,IR). LT .IED(I1))  IED(I1)=IE(I1,IR)
    5 CONTINUE
      DO 6 IR=1,IRE
    6 IE(I1,IR)=IE(I1,IR)-IED(I1)
    4 CONTINUE
      IR=0
      DO 7 IC=IC1,IC2
      IR=IR+1
      CL(IC)=1.D0
      DO 8 I1=1,IFMX
    8 CL(IC)=CL(IC)*IFI(I1)**IE(I1,IR)
      CL(IC)=ISI*CL(IC)
      ISI=-ISI
    7 CONTINUE
      IF(M.EQ.0) IL2P(1)=IL2P(1)-1
      UP  =1.D0
      UPSQ=1.D0
      DOWN=1.D0
      DO 17 I1=1,IFMX
      IEUPSQ=2*IED(I1)+IL2P(I1)+L1(I1)
      IEINT =IEUPSQ/2-L2(I1)
      IEMOD =MOD(IEUPSQ,2)
      UPSQ =UPSQ*IFI(I1) **IEMOD
      IF(IEINT.GE.0)                   T H E N
      UP   =UP  *IFI(I1) **IEINT
                                       E L S E
      DOWN =DOWN*IFI(I1) **(-IEINT)
                                       E N D    I F
   17 CONTINUE
      COE(ICE)=SQRT(UPSQ)* UP / DOWN
      if (TEST('SHAPE   ')) 
     &      WRITE(6,201) L,M,UP,UPSQ,DOWN,(CL(IC),IC=IC1,IC2)
      IF(M. EQ .0)  GO TO 20
      LA=L+M
      LB=L-M+1
      CALL REDUCE(LA,IFMX,IFI,IEA)
      CALL REDUCE(LB,IFMX,IFI,IEB)
      DO 15 I1=1,IFMX
      JM0(I1)=JM0(I1)+IEA(I1)-IEB(I1)
   15 L1 (I1)=L1 (I1)-IEA(I1)+IEB(I1)
      M=M-1
      GO TO 2
   20 CONTINUE
      if (TEST('SHAPE   ')) WRITE(6,202)
      IF(L. EQ .LMAX)   GO TO 10
      LA=(2*L+1)*(2*L+2)
      LB=(L+1)*2
      CALL REDUCE(LA,IFMX,IFI,IEA)
      CALL REDUCE(LB,IFMX,IFI,IEB)
      DO 16 I1=1,IFMX
      L1ST(I1)=L1ST(I1)+IEA(I1)
   16 L2ST(I1)=L2ST(I1)+IEB(I1)
      L=L+1
      GO TO 1
   10 CONTINUE
      RETURN
  100 WRITE(6,204) LMAX,LMAXD1,ICMAX,ICD
      STOP
  201 FORMAT(2X,'L=',I2,' M=',I2,F10.3,' *SQRT(',F16.2,')/',F10.3/
     *       2X,'CL  :',6F14.2)
  202 FORMAT(80('*'))
  203 FORMAT(13X,'THERE ARE',I5,'  COEFFICIENTS'/)
  204 FORMAT(13X,'FROM CCOEF: INCONSISTENCY DATA-DIMENSION'/
     *       14X,'LMAX:',2I5/13X,'ICMAX:',2I5)
      END
      SUBROUTINE DREAL(LMAX,ALPHA,BETA,GAMMA,ITEMP)
      implicit none
C------------------------------------------------------------------
C     THIS ROUTINE COMPUTES TRANSFORMATION MATRICES ASSOCIATED TO
C     THE ROTATION THROUGH THE EULER ANGLES ALPHA,BETA,GAMMA  FOR
C     REAL SPHERICAL HARMONICS UP TO QUANTUM NUMBER LMAX. THE RE-
C     SULTS ARE STORED IN THE TEMPORARY FILE IN UNIT ITEMP.
C------------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'
c      INTEGER LMAXD,ISUMD
c     PARAMETER (LMAXD=25,ISUMD=23426)
c      PARAMETER (LMAXD=25,ISUMD=100000)
      INTEGER ISUMD
      PARAMETER (ISUMD=100000)
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   LMAX,ITEMP
      REAL*8    ALPHA,BETA,GAMMA
C
C     .. LOCAL SCALARS ..
C
      INTEGER   L,M,MP,I,IMAX,IP,IPMAX,ISU,ISUM
      REAL*8    SQR2,FAC,FAC1,FAC2,D,D1,D2
C
C     .. LOCAL ARRAYS ..
C
      REAL*8 DMN(LMAXD1+1,LMAXD1+1),DPL(LMAXD1+1,LMAXD1+1),DMATL(ISUMD)
C
C     .. EXTERNAL ROUTINES ..
C
      REAL*8 DROT
      EXTERNAL DROT
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC COS,SIN,MOD,SQRT
C-----------------------------------------------------------------------
      SQR2=SQRT(2.D0)
      ISU=0
      DO 2 L=0,LMAX
      FAC2=1.D0
      M=0
    1 FAC1=1.D0
      MP=0
    3 CONTINUE
      FAC=FAC1*FAC2/2.D0
      D1=DROT(L,MP ,M,BETA)
      D2=DROT(L,MP,-M,BETA)
      I F ( M O D ( M , 2       ) .N E. 0 )  D 2 = - D 2
      DPL(MP+1,M+1)=(D1+D2)*FAC
      DMN(MP+1,M+1)=(D1-D2)*FAC
      I F ( M O D ( M + M P , 2 ) .N E. 0 )  G O  T O  4
      DPL(M+1,MP+1)=DPL(MP+1,M+1)
      DMN(M+1,MP+1)=DMN(MP+1,M+1)
                                             G O  T O  5
    4 DMN(M+1,MP+1)=-DMN(MP+1,M+1)
      DPL(M+1,MP+1)=-DPL(MP+1,M+1)
    5 CONTINUE
      FAC1=SQR2
      MP=1+MP
      IF (MP .LE. M) GO TO 3
      FAC2=SQR2
      M=1+M
      IF (M .LE. L) GO TO 1
      IMAX=1
      M=0
   12 CONTINUE
      DO 13 I=1,IMAX
      IPMAX=1
      MP=0
   11 CONTINUE
      DO 6 IP=1,IPMAX
      I F ( I   . E Q . 2 )                      G O  T O   7
      I F ( I P . E Q . 2 )                      G O  T O  10
      D= COS(MP*ALPHA)*COS(M*GAMMA)*DPL(MP+1,M+1)
     *  -SIN(MP*ALPHA)*SIN(M*GAMMA)*DMN(MP+1,M+1)
                                                 G O  T O   9
    7 I F ( I P . E Q . 2 )                      G O  T O   8
      D=-COS(MP*ALPHA)*SIN(M*GAMMA)*DPL(MP+1,M+1)
     *  -SIN(MP*ALPHA)*COS(M*GAMMA)*DMN(MP+1,M+1)
                                                 G O  T O   9
    8 D=-SIN(MP*ALPHA)*SIN(M*GAMMA)*DPL(MP+1,M+1)
     *  +COS(MP*ALPHA)*COS(M*GAMMA)*DMN(MP+1,M+1)
                                                 G O  T O   9
   10 D= SIN(MP*ALPHA)*COS(M*GAMMA)*DPL(MP+1,M+1)
     *  +COS(MP*ALPHA)*SIN(M*GAMMA)*DMN(MP+1,M+1)
    9 CONTINUE
C THIS IS CHANGED
      IF(MOD(M+MP,2) . NE . 0)  D=-D
C
      ISU=ISU+1
      DMATL(ISU)=D
    6 CONTINUE
      IPMAX=2
      MP=MP+1
      IF(MP.LE.L) GO TO 11
   13 CONTINUE
      IMAX=2
      M=M+1
      IF(M .LE.L) GO TO 12
    2 CONTINUE
      ISUM=ISU
      WRITE(ITEMP) (DMATL(ISU),ISU=1,ISUM)
      RETURN
      END
      FUNCTION DROT(L,MP,M,BETA)
      implicit none
C-----------------------------------------------------------------------
C     CALCULATION OF D COEFFICIENT ACCORDING TO ROSE, ELEMENTARY THEORY
C     ANGULAR MOMENTUM,J.WILEY & SONS ,1957 , EQ. (4.13).
C-----------------------------------------------------------------------
C
C     .. SCALAR ARGUMENTS ..
C
      INTEGER   L,M,MP
      REAL*8    BETA,DROT
C
C     .. LOCAL SCALARS ..
C
      INTEGER   L0,M0,MP0,N1,N2,N3,N4,NN,I,KMIN,KMAX,LTRM,N,K
      REAL*8    SINB,TERM,FF,BETA2,COSB
C
C     .. LOCAL ARRAYS ..
C
      INTEGER        NF(4)
      EQUIVALENCE (N1,NF(1)),(N2,NF(2)),(N3,NF(3)),(N4,NF(4))
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC IABS,ABS,SQRT,COS,SIN,MOD,MIN0,MAX0
C-----------------------------------------------------------------------
      DATA L0,M0,MP0/-1,1,1/
         L0=-1
         M0=1
         MP0=1
   10 IF(L.NE.L0) GO TO 2
                           WRITE(6,300) L0,M0,MP0
  300                      FORMAT(3X,'L0,M0,MP0=',3I3)
      IF((IABS(M ).EQ.IABS(M0).AND.IABS(MP).EQ.IABS(MP0)).OR.
     1   (IABS(MP).EQ.IABS(M0).AND.IABS(M ).EQ.IABS(MP0))) GO TO 1
    2 FF=1.D0
      IF(IABS(M).LE.L.AND.IABS(MP).LE.L) GO TO 3
      WRITE(6,100) L,M,MP
  100 FORMAT('     L=',I5,'    M=',I5,'    MP =',I5)
      RETURN
    3 N1   =L+M
      N2   =L-M
      N3   =L+MP
      N4   =L-MP
      L0=L
      M0=M
      MP0=MP
      DO 4 N=1,4
      NN=NF(N)
      IF(NN.EQ.0) GO TO 4
      DO 5 I=1,NN
    5 FF=FF*I
    4 CONTINUE
      FF=SQRT(FF)
    1 BETA2=BETA/2.D0
      COSB= COS(BETA2)
      SINB=-SIN(BETA2)
      IF(ABS(COSB).LT.1.E-4) GO TO 9
      IF(ABS(SINB).LT.1.E-4) GO TO 11
      KMAX=MIN0(L-MP,L+M)
      KMIN=MAX0(M-MP,0)
      TERM=COSB**(2*L+M-MP-2*KMIN)*SINB**(MP-M+2*KMIN)*FF
      GO TO 12
    9 LTRM=L
      TERM=FF
      IF(SINB.LT.0.D0.AND.MOD(MP-M,2).NE.0) TERM=-TERM
      GO TO 14
   11 LTRM=0
      TERM=FF
      IF(COSB.LT.0.D0.AND.MOD(MP-M,2).NE.0) TERM=-TERM
   14 KMAX=M-MP
      IF(MOD(KMAX,2).NE.0) GO TO 13
      KMAX=LTRM+KMAX/2
      IF(KMAX.LT.MAX0(M-MP,0)) GO TO 13
      IF(KMAX.GT.MIN0(L-MP,L+M)) GO TO 13
      KMIN=KMAX
   12 IF(MOD(KMIN,2).NE.0) TERM=-TERM
      N1   =L-MP-KMIN
      N2   =L+M-KMIN
      N3   =KMIN+MP-M
      N4   =KMIN
      DO 6 N=1,4
      NN=NF(N)
      IF(NN.EQ.0) GO TO 6
      DO 7 I=1,NN
    7 TERM=TERM/I
    6 CONTINUE
      DROT=TERM
      IF(KMIN.EQ.KMAX) RETURN
      KMIN=KMIN+1
      COSB=COSB**2
      SINB=SINB**2
      N3=N3   +1
      DO 8 K=KMIN,KMAX
      TERM=-N1*N2*TERM*SINB/(COSB*K*N3)
      DROT=DROT+TERM
      N1=N1-1
      N2=N2-1
    8 N3=N3+1
      RETURN
   13 DROT=0.D0
      RETURN
      END
      SUBROUTINE INTSIM(FD,RATIO,X1,X2,DLT,XEL)
      implicit none
      REAL*8 FD,RATIO,X1,X2,DLT
      REAL*8 H,X
      INTEGER   I,INC,INB,N,K
      REAL*8 FFF
      EXTERNAL FFF
      REAL*8 XEL(5)
      INTRINSIC DFLOAT

      DO 1 I=1,5
    1 XEL(I)=0.D0
      N=(X2-X1)/DLT + 1
      INC=2*N
      INB=INC-1
      H=(X2-X1)/DFLOAT(INC)
      DO 2 K=1,INB,2
      X=DFLOAT(K)*H+X1
      DO 12 I=1,5
   12 XEL(I)=XEL(I)+4.D0*FFF(I,X,FD,RATIO)
    2 CONTINUE
      DO 3 K=2,INB,2
      X=DFLOAT(K)*H+X1
      DO 13 I=1,5
   13 XEL(I)=XEL(I)+2.D0*FFF(I,X,FD,RATIO)
    3 CONTINUE
      DO 14 I=1,5
   14 XEL(I)=XEL(I)+FFF(I,X1,FD,RATIO)+FFF(I,X2,FD,RATIO)
      DO 15 I=1,5
   15 XEL(I)=H*XEL(I)/3.D0
      RETURN
      END

      REAL*8 FUNCTION FFF(I,X,FD,RATIO)
      implicit none
      REAL*8 X,FD,RATIO
      REAL*8 AEXP,A1,A2,A,C1,C2,C,B1,B
      INTEGER   I
      INTRINSIC SIN,COS,SQRT

      AEXP=3.D0/2.D0
      A1=RATIO*COS(X-FD)*(1.D0-RATIO*RATIO)
      A2=(1.D0-RATIO*RATIO*SIN(X-FD)*SIN(X-FD))**AEXP
      A=A1/A2
      C1=RATIO*COS(X-FD)
      C2=RATIO*RATIO*(1+2.D0*SIN(X-FD)*SIN(X-FD))-3.D0
      C=2.D0+C1*C2/A2
      B1=(1.D0-RATIO*RATIO)**AEXP
      B=B1/A2
      IF(I.EQ.1) FFF=SQRT(15.D0     )*SIN(2.D0*X)*C/6.D0
      IF(I.EQ.2) FFF=-SQRT(5.D0/3.D0)*SIN(     X)*B
      IF(I.EQ.3) FFF=SQRT( 5.D0     )            *A/2.D0
      IF(I.EQ.4) FFF=-SQRT(5.D0/3.D0)*COS(     X)*B
      IF(I.EQ.5) FFF=SQRT(15.D0     )*COS(2.D0*X)*C/6.D0
      RETURN
      END
      SUBROUTINE MESH0(CRT,NM,NPAN,NAPROX,NMIN)
      IMPLICIT NONE
C ***********************************************************
C *  THIS SUBROUTINE CALCULATES AN APPROPRIATE MESH FOR
C *  THE SHAPE FUNCTIONS. MORE THAN NMIN POINTS BETWEEN TWO
C *  CRITICAL POINTS
c *  In case of more dense mesh increase NMIN 
C *
C ***********************************************************
c      INTEGER NPAND,MESHND
c      PARAMETER (NPAND=80,MESHND=1000)
       include 'inc.geometry'
       INTEGER MESHND
       PARAMETER (MESHND=IRID)
C
C
      REAL*8 CRT(NPAND)
      REAL*8 DIST,D1
      INTEGER   NM(NPAND)
      INTEGER   NAPROX,NPAN,NMIN,N,NTOT,I,NT,NA
      INTRINSIC DFLOAT
c     DATA NMIN/3/  ! 7
      IF ((NPAN-1)*NMIN.GE.NAPROX) THEN
      WRITE(6,*) NPAN,NMIN,NAPROX
      STOP ' INCREASE NUMBER OF POINTS'
      END IF
      DIST=ABS(CRT(1)-CRT(NPAN))
      DO 1 I=1,NPAN-1
      D1=ABS(CRT(I)-CRT(I+1))
      NM(I)=NAPROX*D1/DIST
      IF (NM(I).LT.NMIN) THEN
      NM(I)=NMIN
      END IF
    1 CONTINUE
      N = NMIN*(NPAN-1)
      N = NAPROX -N
      IF (N.LE.0) THEN
         WRITE(*,*) NAPROX,NMIN*(NPAN-1),N
         STOP '*** INCREASE NUMBER OF MESH POINTS ***'
      ENDIF
      D1=DFLOAT(N)/DFLOAT(NAPROX)
      NTOT=N
      DO 3 I=1,NPAN-1
      NA = NINT(D1*FLOAT(NM(I)))
      IF (NM(I).GT.NMIN.AND.(NTOT-NA).GT.0) THEN
      NM(I)=NMIN+NA
      NTOT=NTOT-NA
      END IF
  3   CONTINUE
      NM(1) = NM(1) + NTOT 
      RETURN
      END
C=====================================================================
      SUBROUTINE POLCHK(NFACE,NVERTICES,XVERT,YVERT,ZVERT,TOLVDIST)
      IMPLICIT NONE
C     ----------------------------------------------------------------
C     THIS SUBROUTINE READS THE COORDINATES OF THE VERTICES OF EACH
C     (POLYGON)  FACE OF  A CONVEX POLYHEDRON AND  CHECKS  IF THESE 
C     VERTICES ARRANGED  CONSECUTIVELY DEFINE A  POLYGON. THEN  THE 
C     SUBROUTINE  DETERMINES  THE  VERTICES  AND  THE  EDGES OF THE 
C     POLYHEDRON AND CHECKS IF  THE  NUMBER  OF  VERTICES  PLUS THE  
C     NUMBER OF FACES EQUALS THE NUMBER OF EDGES PLUS 2.
C
C     DATA ARE READ FROM FILE IN UNIT 7, WHICH WE FINALLY REWIND
C     ----------------------------------------------------------------
C
C     .. PARAMETER STATEMENTS ..
C
      include 'inc.geometry'

      INTEGER NVRTD,NEDGED
      PARAMETER (NVRTD=500,NEDGED=NVRTD+NFACED-2)
c
c     ...Arrays ......
c
      INTEGER   NVERTICES(NFACED)
      REAL*8    XVERT(NVERTD,NFACED),YVERT(NVERTD,NFACED),
     &          ZVERT(NVERTD,NFACED)
c     ...Scalars....
      REAL*8 TOLVDIST
C
C     .. LOCAL SCALARS ..
C
      INTEGER   IVERT,INEW,IVERTP,IVERTM,IVRT,IEDGE,NVRT,NEDGE
      INTEGER   IFACE,NFACE,NVERT,I
      REAL*8    ARG,A1,A2,DOWN,UP,FISUM,T
      REAL*8    VRTX,VRTY,VRTZ,VRTPX,VRTPY,VRTPZ,VRTMX,VRTMY,VRTMZ
      REAL*8    PI314
C
C     .. LOCAL ARRAYS ..
C
      REAL*8    V1(3,NEDGED),V2(3,NEDGED),V(3,NVERTD),VRT(3,NVRTD)
C
C     .. INTRINSIC FUNCTIONS ..
C
      INTRINSIC ABS,ACOS,SIGN,SQRT
C     ----------------------------------------------------------------
c$      READ(7,100) NFACE,LDUM,KDUM,DDUM
      PI314 = 4.D0*DATAN(1.D0)

      NVRT=0
      NEDGE=0
      DO 10 IFACE=1,NFACE
c$      READ(7,101) DUM1,DUM2,DUM3,DUM4,NVERT
      NVERT = NVERTICES(IFACE)
      FISUM=(NVERT-2)*PI314
      !write(6,*) 'starting ',fisum
      DO 25 IVERT=1,NVERT
c$      READ(7,102) (V(I,IVERT),I=1,3)
        V(1,IVERT) = XVERT(IVERT,IFACE)
        V(2,IVERT) = YVERT(IVERT,IFACE)
        V(3,IVERT) = ZVERT(IVERT,IFACE)
   25 CONTINUE
C
C------> T R E A T M E N T   O F   V E R T I C E S
C
      DO 2 IVERT =1,NVERT
      VRTX=V(1,IVERT)
      VRTY=V(2,IVERT)
      VRTZ=V(3,IVERT)
      INEW=1                          ! Save all different vertices
      DO 13 IVRT=1,NVRT
      T=(VRTX-VRT(1,IVRT))**2+(VRTY-VRT(2,IVRT))**2
     & +(VRTZ-VRT(3,IVRT))**2
      IF(T.LT.TOLVDIST) INEW=0
   13 CONTINUE
      IF(INEW.EQ.1)                  T H E N
      NVRT=NVRT+1
      IF(NVRT.GT.NVRTD) STOP 'INCREASE NVRTD'
      VRT(1,NVRT)=V(1,IVERT)
      VRT(2,NVRT)=V(2,IVERT)
      VRT(3,NVRT)=V(3,IVERT)
                                     E N D   I F
      IVERTP=IVERT+1                  
      IF(IVERT.EQ.NVERT) IVERTP=1
      VRTPX=V(1,IVERTP)
      VRTPY=V(2,IVERTP)
      VRTPZ=V(3,IVERTP)
      IVERTM=IVERT-1
      IF(IVERT.EQ.1) IVERTM=NVERT
      VRTMX=V(1,IVERTM)
      VRTMY=V(2,IVERTM)               ! Check if the  consecutive
      VRTMZ=V(3,IVERTM)               ! vertices define a polygon
      A1=SQRT((VRTPX-VRTX)**2+(VRTPY-VRTY)**2+(VRTPZ-VRTZ)**2)
      A2=SQRT((VRTMX-VRTX)**2+(VRTMY-VRTY)**2+(VRTMZ-VRTZ)**2)
      DOWN=A1*A2
      UP=(VRTPX-VRTX)*(VRTMX-VRTX)+(VRTPY-VRTY)*(VRTMY-VRTY)+
     &   (VRTPZ-VRTZ)*(VRTMZ-VRTZ)
      ! write(6,*)
      ! write(6,*) VRTX,VRTY,VRTZ
      ! write(6,*) VRTMX,VRTMY,VRTMZ
      ! write(6,*) VRTPX,VRTPY,VRTPZ
      ! write(6,*) 'fisum ',ivert,a1,a2,up,arg,ACOS(ARG)
      IF(DOWN.GE.TOLVDIST)   T H E N
      ARG=UP/DOWN
      IF(ABS(ARG).GE.1.D0) ARG=SIGN(1.D0,ARG)
      FISUM=FISUM-ACOS(ARG)
      
                             E L S E
      STOP 'IDENTICAL CONSECUTIVE VERTICES'
                             E N D    I F
C
C------> T R E A T M E N T   O F   E D G E S 
C
      INEW=1                          ! Save all different edges
      DO 14 IEDGE=1,NEDGE
      T=(VRTX-V1(1,IEDGE))**2+(VRTY-V1(2,IEDGE))**2 
     & +(VRTZ-V1(3,IEDGE))**2
      IF(T.LT.TOLVDIST)         THEN
      T=(VRTPX-V2(1,IEDGE))**2+(VRTPY-V2(2,IEDGE))**2 
     & +(VRTPZ-V2(3,IEDGE))**2
      IF(T.LT.TOLVDIST) INEW=0
                             ELSE
      T=(VRTX-V2(1,IEDGE))**2+(VRTY-V2(2,IEDGE))**2 
     & +(VRTZ-V2(3,IEDGE))**2
      IF(T.LT.TOLVDIST) THEN
      T=(VRTPX-V1(1,IEDGE))**2+(VRTPY-V1(2,IEDGE))**2 
     & +(VRTPZ-V1(3,IEDGE))**2
      IF(T.LT.TOLVDIST) INEW=0
                     END IF
                             END IF
   14 CONTINUE
      IF(INEW.EQ.1)                T H E N
      NEDGE=NEDGE+1
      IF(NEDGE.GT.NEDGED) STOP 'INSUFFICIENT NEDGED'
      V1(1,NEDGE)=V(1,IVERT )
      V1(2,NEDGE)=V(2,IVERT )
      V1(3,NEDGE)=V(3,IVERT )
      V2(1,NEDGE)=V(1,IVERTP)
      V2(2,NEDGE)=V(2,IVERTP)
      V2(3,NEDGE)=V(3,IVERTP)
                                   E N D   I F
    2 CONTINUE
      IF(FISUM.GT.1.D-6) THEN
      write(6,*) 'fisum =',fisum
      STOP 'NOT CONSECUTIVE VERTICES OF A POLYGON'
      END IF
   10 CONTINUE
      IF((NVRT+NFACE).NE.(NEDGE+2)) THEN
       WRITE(6,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
       WRITE(6,*) '            WARNING FROM SHAPE      '
!      WRITE(6,*) '   >>  STOP ILLEGAL POLYHEDRON      '
       WRITE(6,*) 'NVRT=',NVRT,' ; NFACE=',NFACE,' ; NEDGE=',NEDGE
c changed 6.10.2000
c       IF((NVRT+NFACE).NE.(NEDGE+2)) STOP 'ILLEGAL POLYHEDRON'
      END IF 
      REWIND 7 
      RETURN
  100 FORMAT(3I5,F10.5)
  101 FORMAT(4F16.8,2I5)
  102 FORMAT(4F16.8)
      END
