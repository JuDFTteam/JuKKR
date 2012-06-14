c*********************************************************************
      SUBROUTINE BRYDBM_com(VISP,V,VINS,LMPOT,
     +                  R,DRDI,ALPHA,IRC,IRMIN,NSPIN,
     +                  IATOM,IMIX,ITER,
     +                  UI2,VI2,WIT,SM1S,FM1S,
     >                  MYLRANK,
     >                  communicator,
C                       new input parameters after inc.p removal
     &                  itdbryd, irmd, irnsd, nspind)
c*********************************************************************
c     imix :
c       3      broyden's            f i r s t  m e t h o d
c       4      broyden's          s e c o n d  m e t h o d
c       5      anderson's     g e n e r a l i z e d   m e t h o d
c
c     implemented here according to notes of s.b.
c     broyden's iteration scheme following the papers of :
c     srivastava, j. phys. , 17 (1984) , pp l317
c     c.g. broyden in math.comput., 19 , pp 577, 1965
c     c.g. broyden in ibid, 21 ,pp 368 ,1967
c     the method has been generalized to include a metric.the
c     definition of the necessary inner products are similar to the
c     discription given in the notes of m.weinert.the algorithm
c     discribed in the paper srivastava  has been simplified
c     ( see notes of s.b.)
c     the files ui,vi are stored on high speed ssd memory.
c     broyden's update treats charge and spin on the same footing
c                  s. bluegel , kfa , may 1987
c     the anderson method (d.g. anderson, j. acm 12, 547 (1964)) has
c     been generalized and reformulated as an improvement of broyden's
c     second method. successive linesearch is replaced by successive
c     search on hyperplanes. ( see notes of s.b. )
c                  s. bluegel , issp , july 1989
c
c     modified for non spherical potential
c                  b. drittler , aug. 1988
c
c     parallelized over the number of atoms 
c                  a. thiess , jun. 2008
c*********************************************************************
      IMPLICIT NONE

      INCLUDE 'mpif.h'

C     array dimension for broyden mixing history
      INTEGER itdbryd
C     number of radial mesh points - total
      INTEGER irmd
C     number of radial mesh points of non-spherical part
      INTEGER irnsd
C     number of spin directions
      INTEGER nspind

C      INTEGER LMPOTD
C      PARAMETER (LMPOTD= (LPOTD+1)**2)
C      INTEGER IRMIND
C      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NTIRD
C      PARAMETER (NTIRD=(IRMD+(IRNSD+1)*(LMPOTD-1))*NSPIND)
      INTEGER ITDTHD
      PARAMETER (ITDTHD=40)
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),R(IRMD,*)
      DOUBLE PRECISION V(IRMD,LMPOT,*)
C     DOUBLE PRECISION VINS(IRMIND:IRMD,LMPOT,*)
      DOUBLE PRECISION VINS(IRMD-IRNSD:IRMD,LMPOT,*)
      DOUBLE PRECISION VISP(IRMD,*)
      INTEGER IRC(*),IRMIN(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ONE,RMIXIV,VMDENO,VMNORM,VOLINV,ZERO,
     +                 SMNORM_LOCAL,SMNORM_GLOBAL,DDOT_LOCAL,
     +                 DDOT_GLOBAL,CMM_LOCAL,CMM_GLOBAL 
      INTEGER IJ,IMAP,IR,IRC1,IRMIN1,ISP,IT,LM,MIT,ITER

      INTEGER IERR
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL BRYSH1,BRYSH2,BRYSH3,DAXPY,DSCAL,RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Save statement ..
      SAVE ZERO,ONE
C     ..
C     .. Local Arrays ..
C     .. these arrays are automatic (dynamically allocated)
C        Fortran arrays
      DOUBLE PRECISION AM(2:ITDBRYD),BM(2:ITDBRYD),
     +                 AM_LOCAL(2:ITDBRYD),BM_LOCAL(2:ITDBRYD)

C     the following arrays have dimension NTIRD
C     PARAMETER (NTIRD=(IRMD+(IRNSD+1)*(LMPOTD-1))*NSPIND)
      DOUBLE PRECISION FM  ((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)
      DOUBLE PRECISION FM1 ((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)
      DOUBLE PRECISION G   ((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)
      DOUBLE PRECISION SM  ((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)
      DOUBLE PRECISION SM1 ((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)
      DOUBLE PRECISION SM1S((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)
      DOUBLE PRECISION FM1S((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)
      DOUBLE PRECISION WORK((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)
      DOUBLE PRECISION VI3 ((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)


      DOUBLE PRECISION UI3((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND)

      DOUBLE PRECISION UI2((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND,
     &                     2:ITDBRYD)
      DOUBLE PRECISION VI2((IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND,
     &                     2:ITDBRYD)
      DOUBLE PRECISION WIT(2:ITDBRYD)

C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER IMIX,LMPOT,IATOM,NSPIN
C     ..
C     .. Data statements ..
      DATA ZERO,ONE/0.0D0,1.0D0/
C     ..
C     .. MPI variables ..

C     .. L-MPI ..
      INTEGER      MYLRANK,
     +             communicator
C
      EXTERNAL MPI_ALLREDUCE
C
C
      NTIRD=(IRMD+(IRNSD+1)*(LMPOT-1))*NSPIND

      IF (MYLRANK.EQ.0) WRITE(*,*) 'ITERATION: ',ITER

      IF (ITDBRYD.GT.ITDTHD .OR. ITDTHD.GT.200) CALL RCSTOP('ITDBRYD  ')
      IF (IMIX.LE.2 .OR. IMIX.GT.5) CALL RCSTOP('IMIXD   ')

      MIT = ITER
      DO WHILE (MIT.GT.ITDBRYD)
        MIT = MIT - ITDBRYD
      ENDDO

C     For the moment the broyden mixing is started from scratch if ITDBRYD 
C     is exceeded - in principle also possible to make the transition
C     continous                                         16/9/2008        


      IF (MYLRANK.EQ.0) THEN
      IF (IMIX.EQ.3) WRITE (6,FMT='('' broyden"s 1st method used '')')
      IF (IMIX.EQ.4) WRITE (6,FMT='('' broyden"s 2nd method used '')')
      IF (IMIX.EQ.5) WRITE (6,FMT='('' gen. anderson method used '')')
      WRITE(6,*) 'ITERATION: ',MIT
      ENDIF

c
      RMIXIV = ONE/ALPHA
c
c---->  the following block is activated only one iteration before
c        broyden iteration scheme is used
c---->  set up of : sm1 = rho(1) ; fm1=fm[1]=f(rho(1)) - rho(1) ;
c                   metric  g := r*r*drdi
c---->  map data of all muffin-tin spheres into one single vector
c
c
      CALL BRYSH3(SM1,VISP,VINS,IRMIN,IRC,IATOM,
     &            NSPIN,IMAP,LMPOT,
     &            IRMD, IRNSD)


      CALL BRYSH1(FM1,V,IRMIN,IRC,IATOM,
     +            NSPIN,IMAP,LMPOT,
     &            irmd)


      IF (IMAP.GT.NTIRD) CALL RCSTOP('NIRDBRY ')

      DO 10 IJ = 1,IMAP
        FM1(IJ) = RMIXIV* (FM1(IJ)-SM1(IJ))
   10 CONTINUE
c
      IJ = 0
      DO 60 ISP = 1,NSPIN
C        DO 50 IA = NATPS,NAEZ
C         !!! get the atom number from somewhere
          IRC1 = IRC(IATOM)
          VOLINV = 3.0D0/(R(IRC1,IATOM)**3)
          DO 20 IR = 1,IRC1
            IJ = IJ + 1
            G(IJ) = VOLINV*R(IR,IATOM)*R(IR,IATOM)*DRDI(IR,IATOM)
   20     CONTINUE
C
          IF (LMPOT.GT.1) THEN
            IRMIN1 = IRMIN(IATOM)
            DO 40 LM = 2,LMPOT
              DO 30 IR = IRMIN1,IRC1
                IJ = IJ + 1
                G(IJ) = VOLINV*R(IR,IATOM)*R(IR,IATOM)*DRDI(IR,IATOM)
   30         CONTINUE
   40       CONTINUE
          END IF

   50   CONTINUE

   60 CONTINUE
c
C=========================================================================
C=========================================================================
C=========================================================================
C=====  For MIT GT 1 activ  ==============================================
C=========================================================================
C=========================================================================
C=========================================================================

      IF (MIT.GT.1) THEN

        DO IJ = 1,IMAP
          SM1(IJ)=SM1S(IJ)
          FM1(IJ)=FM1S(IJ)
        ENDDO      

c
c----> map rho(m) of all mt-spheres into one single vector
c
        CALL BRYSH3(SM,VISP,VINS,IRMIN,IRC,IATOM,
     &              NSPIN,IMAP,LMPOT,
     &              IRMD, IRNSD)

c
c----> map f[m] = f(m) - rho(m) = f(rho(m)) - rho(m) of all mt-spheres
c      into one single vector
c
        CALL BRYSH1(FM,V,IRMIN,IRC,IATOM,
     +              NSPIN,IMAP,LMPOT,
     &              irmd)

        DO 70 IJ = 1,IMAP
          FM(IJ) = RMIXIV* (FM(IJ)-SM(IJ))
   70   CONTINUE
c
c----> calculate  sm = rho(m) - rho(m-1)
c----> calculate dfm = f[m] - f[m-1]
c
        DO 80 IJ = 1,IMAP
          SM1(IJ) = SM(IJ) - SM1(IJ)
          FM1(IJ) = FM(IJ) - FM1(IJ)
   80   CONTINUE

c
c----> loop to generate u[m] = u(ij,mit)
c
        DO 90 IJ = 1,IMAP
          UI3(IJ) = ALPHA*FM1(IJ) + SM1(IJ)
   90   CONTINUE


        DO IT = 2,MIT - 1
          DO IJ = 1,IMAP
          WORK(IJ)=VI2(IJ,IT)
          ENDDO

          AM_LOCAL(IT) = DDOT(IMAP,FM1,1,WORK,1)
        ENDDO

        CALL MPI_ALLREDUCE(AM_LOCAL,AM,(ITDBRYD-1),
     &              MPI_DOUBLE_PRECISION,MPI_SUM,communicator,IERR)

        DO IT = 2,MIT - 1
          DO IJ = 1,IMAP
          WORK(IJ)=UI2(IJ,IT)
          ENDDO
          CALL DAXPY(IMAP,-AM(IT),WORK,1,UI3,1)
        ENDDO

c
c----> print amj = the importance of the history of ui
c
c        WRITE (IPF,FMT='(5x,'' amj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')
c     +    MIT - 1, (AM(IT),IT=2,MIT-1)
c
c

C=========================================================================
C=========================================================================
C=========== BROYD 1 =====================================================
C=========================================================================
C=========================================================================

        IF (IMIX.EQ.3) THEN
c-------->     b r o y d e n ' s   f i r s t   m e t h o d
c
c----> calculate dsmnorm
c
          SMNORM_LOCAL = ZERO
          DO 110 IJ = 1,IMAP
            SMNORM_LOCAL = SMNORM_LOCAL + SM1(IJ)*G(IJ)*SM1(IJ)
  110     CONTINUE

        CALL MPI_ALLREDUCE(SMNORM_LOCAL,SMNORM_GLOBAL,1,
     &              MPI_DOUBLE_PRECISION,MPI_SUM,communicator,IERR)


c
c----> convolute dsm with the metric g
c
          DO 120 IJ = 1,IMAP
            SM1(IJ) = G(IJ)*SM1(IJ)
  120     CONTINUE
c
c----> loop to generate v[m] = v(ij,mit)
c
          DO 130 IJ = 1,IMAP
            VI3(IJ) = ALPHA*SM1(IJ)
  130     CONTINUE


          DO IT = 2,MIT - 1
            DO IJ = 1,IMAP
            WORK(IJ)=UI2(IJ,IT)
            ENDDO
            BM_LOCAL(IT) = DDOT(IMAP,SM1,1,WORK,1)
          ENDDO

        CALL MPI_ALLREDUCE(BM_LOCAL,BM,(ITDBRYD-1),
     &              MPI_DOUBLE_PRECISION,MPI_SUM,communicator,IERR)

          DO IT = 2,MIT - 1
            DO IJ = 1,IMAP
            WORK(IJ)=VI2(IJ,IT)
            ENDDO
            CALL DAXPY(IMAP,-BM(IT),WORK,1,VI3,1)
          ENDDO


c
c----> complete the evaluation of v[m]
c
        DDOT_LOCAL = DDOT(IMAP,SM1,1,UI3,1)

        CALL MPI_ALLREDUCE(DDOT_LOCAL,DDOT_GLOBAL,1,
     &              MPI_DOUBLE_PRECISION,MPI_SUM,communicator,IERR)

        VMDENO = DDOT_GLOBAL - SMNORM_GLOBAL

          IF (ABS(VMDENO).LT.1D-70) CALL RCSTOP('BRY0SN  ')

          CALL DSCAL(IMAP,ONE/VMDENO,VI3,1)
c
c----> print bmj = the importance of the history of vi
c
c          WRITE (IPF,FMT='(5x,'' bmj , ---> j=2,'',i3,/,(9x,1p,7d10.2))'
c     +      ) MIT - 1, (BM(IT),IT=2,MIT-1)
c


C=========================================================================
C=========================================================================
C=========== BROYD 2 =====================================================
C=========================================================================
C=========================================================================

        ELSE IF (IMIX.EQ.4) THEN

c-------->     b r o y d e n ' s   s e c o n d    m e t h o d
c----> calculate v[m] ; convoluted with the metric g

          DO 150 IJ = 1,IMAP
            VI3(IJ) = G(IJ)*FM1(IJ)
  150     CONTINUE

c----> calculate #vm# and normalize v[m]

          DDOT_LOCAL = DDOT(IMAP,VI3,1,FM1,1)

          CALL MPI_ALLREDUCE(DDOT_LOCAL,DDOT_GLOBAL,1,
     &                MPI_DOUBLE_PRECISION,MPI_SUM,communicator,IERR)

          VMNORM = DDOT_GLOBAL

          CALL DSCAL(IMAP,ONE/VMNORM,VI3,1)

C=========================================================================
C=========================================================================
C=========== ANDERSON ====================================================
C=========================================================================
C=========================================================================

        ELSE IF (IMIX.EQ.5) THEN

c-------->     g e n e r a l i z e d   a n d e r s o n   m e t h o d
c
c----> calculate v[m] ; convoluted with the metric g
c
          DO 160 IJ = 1,IMAP
            VI3(IJ) = G(IJ)*FM1(IJ)
  160     CONTINUE
          DO 170 IT = 2,MIT - 1

            CALL DAXPY(IMAP,-AM(IT)*WIT(IT),VI2,1,VI3,1)
  170     CONTINUE

C----> complete the evaluation of v[m]

          DDOT_LOCAL = DDOT(IMAP,FM1,1,VI3,1)

          CALL MPI_ALLREDUCE(DDOT_LOCAL,DDOT_GLOBAL,1,
     &                MPI_DOUBLE_PRECISION,MPI_SUM,communicator,IERR)

          VMDENO = DDOT_GLOBAL

          IF (ABS(VMDENO).LT.1D-70) CALL RCSTOP('BRY1SN  ')

          CALL DSCAL(IMAP,ONE/VMDENO,VI3,1)

C----> save wit(mit) for next iteration

          WIT(MIT) = VMDENO

        END IF

C=========================================================================
C=========================================================================
C============ END MIXING, NOW OUTPUT =====================================
C=========================================================================
C=========================================================================

c
c----> write u3(ij) and v3(ij) on disk
c
        DO IJ = 1,IMAP
           UI2(IJ,MIT)=UI3(IJ)
           VI2(IJ,MIT)=VI3(IJ)
        ENDDO

c
c----> update f[m-1] = f[m]  ; rho(m) = rho(m-1)
c
        DO 180 IJ = 1,IMAP
          FM1(IJ) = FM(IJ)
          SM1(IJ) = SM(IJ)
  180   CONTINUE
c
c----> calculate cmm
c

        CMM_LOCAL = DDOT(IMAP,FM,1,VI3,1)

        CALL MPI_ALLREDUCE(CMM_LOCAL,CMM_GLOBAL,1,
     &              MPI_DOUBLE_PRECISION,MPI_SUM,communicator,IERR)

C           WRITE (IPF,FMT='(5X,'' CMM = '',1P,D12.4)') CMM
c
c----> update rho(m+1)
c

        CALL DAXPY(IMAP,ONE-CMM_GLOBAL,UI3,1,SM,1)

c
c----> map solution back into each mt-sphere
c
        CALL BRYSH2(SM,V,IRMIN,IRC,IATOM,
     &              NSPIN,IMAP,LMPOT,
     &              IRMD)
c

      END IF

C=========================================================================
C=========================================================================
C=========================================================================
C=====  For MIT GT 1 activ  ==============================================
C=========================================================================
C=========================================================================
C=========================================================================

C      MIT = MIT + 1

      DO IJ = 1,IMAP
      SM1S(IJ)=SM1(IJ)
      FM1S(IJ)=FM1(IJ)
      ENDDO

      RETURN

      END
