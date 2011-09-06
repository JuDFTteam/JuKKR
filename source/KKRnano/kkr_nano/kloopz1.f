      SUBROUTINE KLOOPZ1(GMATN,ALAT,IE,IELAST,ITER,
     +                   NAEZ,NOFKS,VOLBZ,BZKP,VOLCUB,CLS,
     +                   NACLS,RR,EZOA,ATOM,GINP_LOCAL,DGINP,
     +                   NSYMAT,DSYMLL,
     &                   TSST_LOCAL,DTDE_LOCAL,
     &                   NUMN0,INDN0,I2,
     >                   NATRC,ATTRC,EZTRC,NUTRC,INTRC,
     <                   SPRS,
     &                   PRSC,EKM,NOITER,
     &                   EZPRE,QMRBOUND,IGUESS,BCP,CNVFAC,
     >                   NXIJ,XCCPL,IXCP,ZKRXIJ,
     <                   LLY_GRDT,TR_ALPH,GMATXIJ,
     >                   LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE,
     >                   LSMYRANK,LSRANK,LSMPIB,LSMPIC)
C
c **********************************************************************
      IMPLICIT NONE
      include 'mpif.h'
C     .. Parameters ..
      include 'inc.p'
      include 'inc.cls'
C
      INTEGER          NSYMAXD
      PARAMETER        (NSYMAXD=48)
      INTEGER          LMMAXD
      PARAMETER        (LMMAXD= (KREL+1) * (LMAXD+1)**2)
      INTEGER          LMGF0D
      PARAMETER        (LMGF0D= (LMAXD+1)**2)
      INTEGER          LLYALM
      PARAMETER        (LLYALM=LLY*(NAEZD*LMMAXD-1)+1)
      DOUBLE COMPLEX   CONE,CZERO
      PARAMETER        (CONE  = ( 1.0D0,0.0D0),CZERO  = ( 0.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
C     ..
      DOUBLE PRECISION ALAT,VOLBZ,QMRBOUND
      INTEGER          IE,IELAST,ITER,NOFKS,I2,NXIJ,
     +                 NAEZ,IGUESS,BCP,EKM,NOITER
C     ..
C     .. Array Arguments ..
C     ..
      DOUBLE COMPLEX   GMATN(LMMAXD,LMMAXD,IEMXD),
     +                 DGINP(LMGF0D,LMGF0D,NACLSD,NCLSD),
     +                 GINP_LOCAL(LMGF0D,LMGF0D,NACLSD,NCLSD),
     +                 GMATXIJ(LMMAXD,LMMAXD,NXIJD),
     +                 GSXIJ(LMMAXD,LMMAXD,NSYMAXD,NXIJD),
     +                 TSST_LOCAL(LMMAXD,LMMAXD),
     +                 DTDE_LOCAL(LMMAXD,LMMAXD),
     +                 LLY_GRDT,TR_ALPH,
     +                 GLL(LMMAXD,LMMAXD),
     +                 GS(LMMAXD,LMMAXD,NSYMAXD)
      DOUBLE PRECISION RR(3,0:NRD),ZKRXIJ(48,3,NXIJD),
     +                 BZKP(3,KPOIBZ),VOLCUB(KPOIBZ),
     +                 CNVFAC(EKMD)
      INTEGER          ATOM(NACLSD,*),
     +                 CLS(*),
     +                 EZOA(NACLSD,*),
     +                 NACLS(*)
      INTEGER          NUTRC,              ! number of inequivalent atoms in the cluster
     +                 INTRC(NATRCD),      ! pointer to atoms in the unit cell
     +                 NATRC,              ! number of atoms in cluster
     +                 ATTRC(NATRCD),      ! index to atom in elem/cell at site in cluster
     +                 EZTRC(NATRCD)       ! index to bravais lattice  at site in cluster
C     .. 
      INTEGER          NUMN0(NAEZD),INDN0(NAEZD,NACLSD),IXCP(NXIJD)
C     ..
C     .. Local Scalars ..
C     ..
      DOUBLE COMPLEX   TAUVBZ
      DOUBLE PRECISION RFCTOR
      INTEGER          LM1,LM2,NSYMAT,INFO,
     +                 IERR,IU
C     ..
C     .. Local Arrays ..
C     ..
C     effective (site-dependent) Delta_t^(-1) matrix
      DOUBLE COMPLEX   MSSQ(LMMAXD,LMMAXD),
     +                 TMATLL(LMMAXD,LMMAXD,NAEZD),
     +                 W1(LMMAXD,LMMAXD),
     +                 DSYMLL(LMMAXD,LMMAXD,NSYMAXD),
     +                 TPG(LMMAXD,LMMAXD),
     +                 XC(LMMAXD,LMMAXD)
      INTEGER          IPVT(LMMAXD)
C
C----- preconditioning ---------------------------------------------------
      COMPLEX          PRSC(NGUESSD*LMMAXD,EKMD)
      DOUBLE COMPLEX   EZPRE(IEMXD)
      INTEGER          SPRS(NGUESSD*LMMAXD+1,EKMD+1)
C=======================================================================
C-----------------------------------------------------------------------
C     .. MPI ..
C     .. L-MPI
      INTEGER      MYLRANK(LMPID*SMPID*EMPID),
     +             LCOMM(LMPID*SMPID*EMPID),
     +             LGROUP(LMPID*SMPID*EMPID),
     +             LSIZE(LMPID*SMPID*EMPID),
     +             LMPI,LMPIC
C     .. LS-MPI
      INTEGER      LSMYRANK(LMPID,NAEZD*SMPID*EMPID),
     +             LSRANK(LMPID,NAEZD*SMPID*EMPID),
     +             LSMPI,LSMPIB,LSMPIC
C     .. External Subroutines ..
      LOGICAL TEST,OPT,XCCPL
C
      EXTERNAL TEST,OPT,KKRMAT01,
     &         ZCOPY,ZAXPY,ZGETRF,ZGETRS,ZGETRI,ZGEMM,ZSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE
C     ..
C
C     RFCTOR=A/(2*PI) conversion factor to p.u.
      RFCTOR = ALAT/(8.D0*ATAN(1.0D0))           ! = ALAT/(2*PI)
C
C
C
      IF ( TEST('flow    ') ) WRITE (6,*) 
     &       '>>> KLOOPZ1: invert delta_t and do Fourier transformation'
C --> convert inverted delta_t-matrices to p.u.
C
            DO LM2 = 1,LMMAXD
               DO LM1 = 1,LM2
                  TSST_LOCAL(LM1,LM2) = 0.5D0/RFCTOR * 
     &                ( TSST_LOCAL(LM1,LM2) + TSST_LOCAL(LM2,LM1) )
                  TSST_LOCAL(LM2,LM1) = TSST_LOCAL(LM1,LM2) 
               END DO
            END DO
C

      CALL MPI_ALLGATHER(TSST_LOCAL,LMMAXD*LMMAXD,MPI_DOUBLE_COMPLEX,
     +     TMATLL,LMMAXD*LMMAXD,MPI_DOUBLE_COMPLEX,
     +     LCOMM(LMPIC),IERR)

C ---------------------------------------------------------------------

C

      DO LM2 = 1,LMMAXD
        DO LM1 = 1,LMMAXD
           MSSQ(LM1,LM2) =  TMATLL(LM1,LM2,I2)
        END DO
      END DO

C
C ---> inversion 
C
        CALL ZGETRF(LMMAXD,LMMAXD,MSSQ,LMMAXD,IPVT,INFO)
        CALL ZGETRI(LMMAXD,MSSQ,LMMAXD,IPVT,W1,
     &              LMMAXD*LMMAXD,INFO)
C
C
C=======================================================================
C
C
C
      TAUVBZ = 1.D0/VOLBZ

      CALL KKRMAT01(BZKP,NOFKS,GS,VOLCUB,VOLBZ,TMATLL,MSSQ,
     &              IE,IELAST,ITER,
     &              ALAT,NSYMAT,NAEZ,CLS,NACLS,RR,EZOA,ATOM,
     &              GINP_LOCAL,DGINP,
     &              NUMN0,INDN0,I2,
     >              NATRC,ATTRC,EZTRC,NUTRC,INTRC,
     <              SPRS,PRSC,
     >              EKM,NOITER,
     &              EZPRE,QMRBOUND,IGUESS,BCP,CNVFAC,
     &              DTDE_LOCAL,
     <              GSXIJ,
     >              NXIJ,XCCPL,IXCP,ZKRXIJ,
     <              LLY_GRDT,TR_ALPH,
     >              LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE,
     >              LSMYRANK,LSRANK,LSMPIB,LSMPIC)
C
C
C-------------------------------------------------------- SYMMETRISE GLL
C
         DO IU = 1,NSYMAT
C
C --->    GLL = sum(i=1,iumax)(tauvbz * DLL(i) * GS * DLL(i)^T)
C
            IF ( IU.EQ.1 ) THEN
C
C --->    ull(1) is equal to unity matrix
C
               CALL ZCOPY(LMMAXD*LMMAXD,GS(1,1,1),1,GLL,1)
               CALL ZSCAL(LMMAXD*LMMAXD,TAUVBZ,GLL,1)
C
            ELSE
C
C --->      tpg = tauvbz * DLL * GS
C
               CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,TAUVBZ,
     &                    DSYMLL(1,1,IU),LMMAXD,GS(1,1,IU),LMMAXD,
     &                    CZERO,TPG,LMMAXD)
C
C --->    GLL = GLL + TPG * DLL(i)^T
C                           C  ! dsymll might be complex in REL case
C
               CALL ZGEMM('N','C',LMMAXD,LMMAXD,LMMAXD,CONE,TPG,LMMAXD,
     &                    DSYMLL(1,1,IU),LMMAXD,CONE,GLL,LMMAXD)
            END IF
C
         END DO
C-------------------------------------------------------- IU = 1,NSYMAT

C
C
C --->  XC = Delta_t^(-1) * GLL
C
        CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,MSSQ,
     &             LMMAXD,GLL,LMMAXD,CZERO,XC,LMMAXD)
C
C
C --->  GLL = - Delta_t^(-1) - Delta_t^(-1) * GLL * Delta_t^(-1)
C
          CALL ZCOPY(LMMAXD**2,MSSQ,1,GLL,1)
          CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,-CONE,XC,LMMAXD,
     &               MSSQ,LMMAXD,-CONE,GLL,LMMAXD)
C
C
C --->  GMATN = GMATLL = GLL/RFCTOR...............renamed
C
        DO LM1 = 1,LMMAXD
          DO LM2 = 1,LMMAXD
            GMATN(LM2,LM1,IE) = GLL(LM2,LM1)/RFCTOR
          END DO
        END DO
C      ENDIF
C
C================================
       IF (XCCPL) THEN
C================================

       CALL SYMJIJ(
     >             ALAT,TAUVBZ,
     >             NSYMAT,DSYMLL,
     >             NXIJ,IXCP,
     >             TMATLL,MSSQ,
     >             GSXIJ,
     <             GMATXIJ)

C================================
      ENDIF
C================================
C
C
      IF ( TEST('flow    ') ) WRITE (6,*) '<<< KLOOPZ1'
C
      RETURN

      END
