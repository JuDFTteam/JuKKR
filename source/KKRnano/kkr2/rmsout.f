C================================================================
C Collects contributions to the rms errors from all sites and
C prints rms errors. Sets new Fermi energy
C Does formatted output of potential if file VFORM exists.
C
C output of a) rms-error
C and       b) (optional) potential in formatted format
C              files are named VPOT.X, X = atom number
C              >> output can be enable providing file VFORM
C
C called by main2
C subroutine called: rites   
C
C================================================================
C
      SUBROUTINE RMSOUT(RMSAVQ,RMSAVM,ITER,E2,EFOLD,
     +                  SCFSTEPS,VBC,QBOUND,NSPIN,NAEZ,
     >                  KXC,LPOT,A,B,IRC,
     >                  VINS,VISP,DRDI,IRNS,R,RWS,RMT,ALAT,
     >                  ECORE,LCORE,NCORE,ZAT,ITITLE,
     >                  LMPIC,MYLRANK,
     >                  LGROUP,LCOMM,LSIZE,
C                       new input parameters after inc.p removal
     &                  irmd, irnsd, prod_lmpid_smpid_empid)
C
      IMPLICIT NONE
C
      INCLUDE 'mpif.h'

      INTEGER irmd
      INTEGER irnsd
      INTEGER prod_lmpid_smpid_empid

C     INTEGER LMPOTD
C     PARAMETER (LMPOTD= (LPOTD+1)**2)
C
      DOUBLE PRECISION QBOUND,RMSAVM,RMSAVQ,RMSAV0,
     +                 E2,EFOLD,EFNEW,ALAT,VBC(2),
     +                 WORK1(2),WORK2(2),
     +                 A(NAEZ),B(NAEZ)

C     DOUBLE PRECISION VINS(IRMIND:IRMD,LMPOTD,2)

      DOUBLE PRECISION VINS(IRMD-IRNSD:IRMD,(LPOT+1)**2,2),
     +                 VISP(IRMD,2),
     +                 DRDI(IRMD,NAEZ),
     +                 ECORE(20,2),
     +                 ZAT(NAEZ),
     +                 R(IRMD,NAEZ),RWS(NAEZ),RMT(NAEZ)

      INTEGER ITER,SCFSTEPS,NSPIN,NAEZ,KXC,LPOT
C
      INTEGER IRNS(NAEZ),IRC(NAEZ),
     +        LCORE(20,NSPIN*NAEZ),
     +        NCORE(NSPIN*NAEZ),
     +        ITITLE(20,NAEZ*NSPIN)

      LOGICAL OPT
C     .. local scalars ..
      DOUBLE PRECISION RMSQ,RMSM
      INTEGER          I1,D1,D10,D100,D1000,OFF(3)
      CHARACTER*12     FNAME
      LOGICAL          VFORM
C     ..
C     .. MPI variables ..
C     .. L-MPI ..
      INTEGER      MYLRANK(prod_lmpid_smpid_empid),
     +             LCOMM(prod_lmpid_smpid_empid),
     +             LGROUP(prod_lmpid_smpid_empid),
     +             LSIZE(prod_lmpid_smpid_empid),
     +             LMPI,LMPIC
C
C     .. N-MPI ..
      INTEGER IERR, MAPBLOCK

      EXTERNAL MPI_REDUCE
      EXTERNAL MAPBLOCK
C     ..

C****************************************************** MPI COLLECT DATA
C
      WORK1(1) = RMSAVQ
      WORK1(2) = RMSAVM

      CALL MPI_ALLREDUCE(WORK1,WORK2,2,
     +                MPI_DOUBLE_PRECISION,MPI_SUM,LCOMM(LMPIC),
     +                IERR)

      RMSQ = SQRT(WORK2(1)/NAEZ)
      RMSM = SQRT(WORK2(2)/NAEZ)

      RMSAV0 = 1.0D10

C****************************************************** MPI COLLECT DATA
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ======================================================================


C ================== MYRANK.EQ.0 =======================================

      IF(MYLRANK(LMPIC).EQ.0) THEN

        WRITE(6,'(79(1H-),/)')
        IF (NSPIN.EQ.2) THEN
          WRITE (6,FMT=9041) ITER,RMSQ,RMSM
        ELSE
          WRITE (6,FMT=9051) ITER,RMSQ
        END IF
        WRITE(6,'(79(1H-))')

 9041 FORMAT ('      ITERATION',I4,' average rms-error : v+ + v- = ',
     +       1p,d11.4,/,39x,' v+ - v- = ',1p,d11.4)
 9051 FORMAT ('      ITERATION',I4,' average rms-error : v+ + v- = ',
     +       1p,d11.4)

C
        IF (ITER.NE.1) RMSAV0 = 1.0d2*MAX(RMSQ,RMSM)
C
      EFNEW = E2
      IF (OPT('rigid-ef')) EFNEW = EFOLD
C
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC CONVERGENCY TESTS
C
      IF (MAX(RMSQ,RMSM).LT.QBOUND) THEN
         WRITE(6,'(17X,A)') '++++++ SCF ITERATION CONVERGED ++++++'
         WRITE(6,'(79(1H*))')
         GO TO 260
C ----------------------------------------------------------------------
      ELSE
C ----------------------------------------------------------------------
         IF (MAX(RMSQ,RMSM).GT.RMSAV0) THEN
            WRITE(6,*) 'ITERATION DIVERGED ---'
            GO TO 260
         END IF
C ----------------------------------------------------------------------
         IF (ITER.GE.SCFSTEPS) THEN
            WRITE(6,'(12X,A)')
     &           '++++++ NUMBER OF SCF STEPS EXHAUSTED ++++++'
            WRITE(6,'(79(1H*))')
            GOTO 260
         END IF
      END IF
C ----------------------------------------------------------------------
C
 260  CONTINUE                  ! jump mark
C
C ======================================================================
C =             write out information for the next iteration           =
C ======================================================================

      OPEN (28,FILE='not.converged',FORM='formatted')
      WRITE (28,'(1P,4D17.10)') EFOLD,VBC
      CLOSE (28)
      END IF

C ============= MYRANK.EQ.0 ============================================
C
C
C ......................................................................
C formatted output                       A.Thiess 09/09
C ......................................................................
C
      VFORM = .FALSE.
      INQUIRE(FILE='VFORM',EXIST=VFORM)
C
      IF (VFORM.AND.(ITER.EQ.SCFSTEPS)) THEN
C
      DO I1 = 1,NAEZ
      IF(MYLRANK(LMPIC).EQ.
     +   MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) THEN
C
C      IRNS(1) = 208 !ART is this necessary?
C
        D1 = mod(I1,10)
        D10 = int( (mod(I1,100) + 0.5)/10 )
        D100 = int( (mod(I1,1000) + 0.5)/100 )
        D1000 = int( (mod(I1,10000) + 0.5)/1000 )
C
        OFF(1) = iachar('1')-1
        OFF(2) = iachar('1')-1
        OFF(3) = iachar('1')-1
C
        IF ( D10.GE.10 ) OFF(1) = iachar('7')
        IF ( D100.GE.100 ) OFF(2) = iachar('7')
        IF ( D1000.GE.1000 ) OFF(3) = iachar('7')
        FNAME='VPOT.'
     +   //achar(D1000+OFF(3))
     +   //achar(D100+OFF(2))
     +   //achar(D10+OFF(1))
     +   //achar(D1+iachar('1')-1)
C
        OPEN(11,FILE=FNAME,FORM='formatted')
C
        CALL RITES(11,I1,NAEZ,NSPIN,ZAT,ALAT,RMT,RMT,RWS,
     +             ITITLE,R,DRDI,VISP,A,B,KXC,IRNS,LPOT,VINS,
     +             QBOUND,IRC,EFNEW,VBC,ECORE,LCORE,NCORE,
     &             irmd, irnsd)
C
        CLOSE (11)
C
      ENDIF
      ENDDO
C
      ENDIF
C ......................................................................
C ......................................................................
C

      END
