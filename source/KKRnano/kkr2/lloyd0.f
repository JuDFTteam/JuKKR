      SUBROUTINE LLOYD0(EZ,WEZ,CLEB,DRDI,R,IRMIN,
     +                  VINS,VISP,THETAS,ZAT,ICLEB,
     +                  IFUNM1,IPAN,IRCUT,LMSP1,JEND,LOFLM,NTCELL,ICST,
     +                  IELAST,IEND,NAEZ,NSPIN,NSRA,
     <                  WEZRN,RNORM,
     >                  GMATN,
     >                  LLY_GRDT,
     >                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU,
     <                  DMATLDAU,
     >                  LMPIC,MYLRANK,
     >                  LGROUP,LCOMM,LSIZE)
C
      IMPLICIT NONE
C
      INCLUDE 'mpif.h'
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'


C     .. Parameters ..
      INTEGER             LMMAXD,LMPOTD
      PARAMETER          (LMPOTD= (LPOTD+1)**2)
      PARAMETER          (LMMAXD= (LMAXD+1)**2)
      INTEGER             LMAXD1
      PARAMETER          (LMAXD1=LMAXD+1)
      INTEGER             MMAXD
      PARAMETER          (MMAXD=2*LMAXD+1)
      INTEGER             LM2D
      PARAMETER          (LM2D= (2*LMAXD+1)**2)
      INTEGER             LMXSPD
      PARAMETER          (LMXSPD= (2*LPOTD+1)**2)
      INTEGER             NPOTD
      PARAMETER          (NPOTD= (2*KREL + (1-KREL)*NSPIND)*NAEZD)
      INTEGER             IRMIND
      PARAMETER          (IRMIND=IRMD-IRNSD)
      DOUBLE COMPLEX      CZERO
      PARAMETER          (CZERO=(0.0D0,0.0D0))
C     ..
C     .. Scalars ..
      INTEGER            I1,ICELL,ICST,IE,IELAST,IEND,
     +                   ISPIN,L,NAEZ,NSPIN,NSRA,NLDAU
      LOGICAL            LDAU
C     ..
C     .. Arrays ..
      DOUBLE COMPLEX WORK1(4*IEMXD),WORK2(4*IEMXD)
C     ..
      DOUBLE COMPLEX DOS(IEMXD,2),DOS0(IEMXD),DOS1(IEMXD)
      DOUBLE COMPLEX DEN0(0:LMAXD1,IEMXD,NSPIND),WEZRN(IEMXD,2)
      DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND)
      DOUBLE COMPLEX LLY_GRDT(IEMXD,NSPIND),
     +               PHILDAU(IRMD,LMAXD1),
     +               DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
      DOUBLE PRECISION DLOYD,DLOYDINT
      DOUBLE PRECISION DLOC,D0LOC,D0LOCINT,D1LOC,D1LOCINT
      DOUBLE PRECISION RNORM(IEMXD,2),
     +                 CLEB(NCLEB,2),DRDI(IRMD,NAEZD),R(IRMD,NAEZD),
     +                 ESPV(0:LMAXD1,NSPIND)
      DOUBLE PRECISION RHO2N1(IRMD,LMPOTD,2),RHO2N2(IRMD,LMPOTD,2)
      DOUBLE PRECISION VINS(IRMIND:IRMD,LMPOTD,2),VISP(IRMD,2),
     +                 THETAS(IRID,NFUND,NCELLD),
     +                 ZAT(NAEZD),
     +                 WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
      INTEGER ICLEB(NCLEB,3),IFUNM1(LMXSPD,NAEZD),IPAN(NAEZD),
     +        IRCUT(0:IPAND,NAEZD),IRMIN(NAEZD),
     +        LMSP1(LMXSPD,NAEZD),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     +        LOFLM(LM2D),NTCELL(NAEZD),
     +        LLDAU(LMAXD1)
C     ..
C     .. MPI ..
C     .. N-MPI
      INTEGER IERR,MAPBLOCK
C      INTEGER MYRANK,NROFNODES
C      COMMON /MPI/MYRANK,NROFNODES
C     .. L-MPI
      INTEGER      MYLRANK(LMPID*SMPID*EMPID),
     +             LCOMM(LMPID*SMPID*EMPID),
     +             LGROUP(LMPID*SMPID*EMPID),
     +             LSIZE(LMPID*SMPID*EMPID),
     +             LMPI,LMPIC
C
C     .. External Subroutines ..
      EXTERNAL CINIT,RHOVAL0,RHOVAL
      EXTERNAL MPI_ALLREDUCE,MPI_COMM_RANK,MPI_COMM_SIZE,
     +         MPI_FINALIZE,MPI_INIT,ZCOPY,DCOPY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
C
      DO IE=1,IELAST
        RNORM(IE,1)=1.D0
        RNORM(IE,2)=1.D0
        DOS0(IE)       =CZERO
        DOS1(IE)       =CZERO
        DO ISPIN = 1, NSPIN
          DOS(IE,ISPIN)=CZERO
        ENDDO
      ENDDO

C=================================================================
C==  calculate DOS  ==============================================
C=================================================================

      CALL CINIT(IEMXD*(LMAXD+2)*NSPIND,DEN0)

        DO I1 = 1,NAEZ
        IF (MYLRANK(LMPIC).EQ.
     +      MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) THEN
          DO ISPIN = 1,NSPIN
            ICELL = NTCELL(I1)

            CALL RHOVAL(.FALSE.,ICST,IELAST,NSRA,
     &                  ISPIN,NSPIN,
     &                  EZ,WEZ,DRDI(1,I1),R(1,I1),IRMIN(I1),
     &                  VINS(IRMIND,1,ISPIN),VISP(1,ISPIN),
     &                  ZAT(I1),IPAN(I1),IRCUT(0,I1),
     &                  THETAS(1,1,ICELL),IFUNM1(1,ICELL),
     &                  LMSP1(1,ICELL),
     &                  RHO2N1,RHO2N2,DEN0(0,1,ISPIN),ESPV(0,ISPIN),
     &                  CLEB,LOFLM,ICLEB,IEND,JEND,
     &                  GMATN,
     >                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU,
     <                  DMATLDAU)

            DO IE = 1,IELAST
              DO L = 0,LMAXD1
             DOS(IE,ISPIN) = DOS(IE,ISPIN) + WEZ(IE)*DEN0(L,IE,ISPIN)
              END DO
            END DO
          END DO
          DO IE = 1,IELAST
            CALL RHOVAL0(
     &           EZ(IE),WEZ(IE),
     &           DRDI(1,I1),R(1,I1),
     &           IPAN(I1),IRCUT(0,I1),
     &           THETAS(1,1,ICELL),
     &           DOS0(IE),DOS1(IE))
          END DO
        ENDIF
        END DO


C==  allreduce DOS  ==============================================

        CALL ZCOPY(IEMXD,DOS0,1,WORK1,1)
        CALL ZCOPY(IEMXD,DOS1,1,WORK1(IEMXD+1),1)
        CALL ZCOPY(2*IEMXD,DOS,1,WORK1(2*IEMXD+1),1)

        CALL MPI_ALLREDUCE(WORK1,WORK2,4*IEMXD,
     +                     MPI_DOUBLE_COMPLEX,MPI_SUM,LCOMM(LMPIC),
     +                     IERR)

        CALL ZCOPY(IEMXD,WORK2,1,DOS0,1)
        CALL ZCOPY(IEMXD,WORK2(IEMXD+1),1,DOS1,1)
        CALL ZCOPY(2*IEMXD,WORK2(2*IEMXD+1),1,DOS,1)

C ================================================================ NAEZ



C ======================================================================

      DO ISPIN=1,NSPIN

        D0LOCINT=0.0D0
        D1LOCINT=0.0D0
        DLOYDINT=0.0D0
        DO IE=1,IELAST

            DLOYD=DIMAG(WEZ(IE)*LLY_GRDT(IE,ISPIN))

            DLOC=DIMAG(DOS(IE,ISPIN))
            D0LOC=DIMAG(DOS0(IE))
            D1LOC=DIMAG(DOS1(IE))

            RNORM(IE,ISPIN) = (DLOYD+D0LOC)/DLOC

            D0LOCINT=D0LOCINT+D0LOC
            D1LOCINT=D1LOCINT+D1LOC
            DLOYDINT=DLOYDINT+DLOYD

            WEZRN(IE,ISPIN) = WEZ(IE)*RNORM(IE,ISPIN)

        ENDDO

          DO IE=1,IELAST
            WEZRN(IE,ISPIN) = WEZRN(IE,ISPIN)*
     +      (DLOYDINT+D0LOCINT-D1LOCINT)/(DLOYDINT+D0LOCINT)
            RNORM(IE,ISPIN) = RNORM(IE,ISPIN)*
     +      (DLOYDINT+D0LOCINT-D1LOCINT)/(DLOYDINT+D0LOCINT)
          ENDDO

      ENDDO

        RETURN
        END
