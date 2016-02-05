c ************************************************************************
      SUBROUTINE GLL95(LSTART,EZ,CLEB,ICLEB,LOFLM,IEND,
     +                 TMATLL,ATOM,KAOEZ,REFPOT,
     +                 RATOM,NATOM,ALAT,OUT,GREF0)
c ************************************************************************
c
c     solution of the DYSON equation for a cluster of potentials
c     (TMATLL) centered at positions RATOM in free space,
c
c     (modified version of GLL91 by P. Zahn, Sept. 95)
c ------------------------------------------------------------------------

C     .. Parameters ..
      include 'inc.p'
      include 'inc.cls'
c
      INTEGER LMAX,NATOMD
c      PARAMETER (LMAX=4,NATOMD=79)
      PARAMETER (LMAX=LMAXD,NATOMD=NACLSD)
      INTEGER LMAXSQ,NGD
      PARAMETER (LMAXSQ= (LMAX+1)**2,NGD=LMAXSQ*NATOMD)
c      INTEGER NCLEB
c      PARAMETER (NCLEB= (LMAX*2+1)**2*LMAXSQ)
      DOUBLE COMPLEX CONE,CZERO,CONEM,CTWO
      PARAMETER (CONE  = (1.D0,0.D0),
     +           CTWO  = (2.D0,0.D0),
     +           CZERO = (0.D0,0.D0),
     +           CONEM = (-1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EZ
      INTEGER IEND,NATOM,OUT
      DOUBLE PRECISION ALAT
      LOGICAL LSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CLEB(*),RATOM(3,*)
c
      DOUBLE COMPLEX GREF0(NGD,*), GREF1(NACLSD*LMAXSQ,LMAXSQ),
     +               TMATLL(LMAXSQ,LMAXSQ,*)
c
      INTEGER 
     +     ICLEB(NCLEB,*),
     +     KAOEZ(*),            ! kind of atom at site in el. cell
     +     LOFLM(*),
     +     REFPOT(*),
     +     ATOM(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,L,LM,LM1,LM2,M,N,N1,N2,NDIM,NLM1,NLM2
      DOUBLE COMPLEX A,B
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX 
     +     GLL(LMAXSQ,LMAXSQ),
     +     GREF(NGD,NGD),
     +     GTREF(NGD,LMAXSQ)
c
      DOUBLE PRECISION RDIFF(3)
c
      INTEGER LF(144)
      DATA LF/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11/
c      INTEGER LF(81)
c      DATA LF/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8/
C     ..
C     .. External Subroutines ..
      LOGICAL TEST,OPT
      EXTERNAL GFREE,GREFSY,OPT,TEST,TMREAD,ZCOPY,ZGEMM
C     ..
C     .. Save statement ..
      SAVE
C     ..
c ------------------------------------------------------------------------
C      IF (LSTART) THEN
C        IF(TEST('flow    ')) write(6,*) 'LSTART=', lstart
C        LM = 0
C        DO 20 L = 1,LMAX + 1
C          DO 10 M = 1,L + L - 1
C            LM = LM + 1
C            LF(LM) = L
C 10       CONTINUE
C 20     CONTINUE
C        LSTART = .FALSE.
C      END IF
c ------------------------------------------------------------------------
      IF(TEST('flow    ')) write(6,*) '>>> GLL95'

      NDIM = LMAXSQ*NATOM
c
c ---> construct free Green's function
c
      DO 90 N1 = 1,NATOM
        DO 80 N2 = 1,NATOM
          DO 30 I = 1,3
c            RDIFF(I) = (RATOM(I,N1) - RATOM(I,N2))*ALAT
c           changed P.Z. 4.7.97
            RDIFF(I) = -(RATOM(I,N1) - RATOM(I,N2))*ALAT
 30       CONTINUE
          IF (N1.NE.N2) THEN
            IF(TEST('N1 N2   ')) 
     +           write(6,*) N1,N2,RDIFF(1),RDIFF(2),RDIFF(3)
            CALL GFREE(RDIFF,EZ,GLL,CLEB,ICLEB,LOFLM,IEND)

c            IF ( (N1 .LE. 5) .AND. (N2 .LE. 5)) then
c              write(434,*) "GFREE"
c              write(434,"((A7),(2I5))") "N1,N2", N1,N2
c            END If

            DO 50 LM2 = 1,LMAXSQ
              NLM2 = (N2-1)*LMAXSQ + LM2
              DO 40 LM1 = 1,LMAXSQ
                NLM1 = (N1-1)*LMAXSQ + LM1

C                IF (TEST('GFREE   ') .AND. 
C     +               N1.EQ.1 .AND.
C     +               abs((GLL(LM1,LM2)-(conem)**(lf(lm1)+lf(lm2))
C     +               *GLL(LM2,LM1))/
C     +               (GLL(LM1,LM2)+(conem)**(lf(lm1)+lf(lm2))*
C     +               GLL(LM2,LM1))).gt.1.d-10 )
C     +               write(6,FMT='(3i4,1p,2d15.6)')
C     +                N2,LM1,LM2,GLL(LM1,LM2)

                IF (TEST('GFREE   ') .AND. 
     +               N1.EQ.1 .AND. LM1.eq.LM2 .AND.
     +               zabs(GLL(LM1,LM2)).gt.1.d-10)
     +               write(6,FMT='(3i4,1p,2d15.6,d17.6,0p,f8.4)')
     +                N2,LM1,LM2,GLL(LM1,LM2),zabs(GLL(LM1,LM2)),
     +               DREAL(GLL(LM1,LM2)/zabs(GLL(LM1,LM2)))
                

                GREF(NLM1,NLM2) = GLL(LM1,LM2)

c                IF ( (N1 .LE. 5) .AND. (N2 .LE. 5)) then
c                  write(434,"((2I5),(2e17.9))") LM2,LM1,GREF(NLM1,NLM2)
c                END If
c                GREF(NLM2,NLM1) = GLL(LM1,LM2)
 40           CONTINUE
 50         CONTINUE
          ELSE

            DO 70 LM2 = 1,LMAXSQ
              NLM2 = (N2-1)*LMAXSQ + LM2
              DO 60 LM1 = 1,LMAXSQ
                NLM1 = (N1-1)*LMAXSQ + LM1
                GREF(NLM1,NLM2) = CZERO

 60           CONTINUE
 70         CONTINUE
          END IF
 80     CONTINUE
 90   CONTINUE
      
      IF(TEST('flow    ')) write(6,*) 'GFREE o.k.'
c ------------------------------------------------------------------------

      CALL ZCOPY(NGD*LMAXSQ,GREF,1,GREF0,1)
        
      DO 100 N2 = 1,NATOM
        NLM2 = (N2-1)*LMAXSQ + 1
        CALL ZGEMM('N','N',NDIM,LMAXSQ,LMAXSQ,
     +       -CONE,GREF(1,NLM2),NGD,
     +       TMATLL(1,1,REFPOT(KAOEZ(ABS(ATOM(N2))))),
     +       LMAXSQ,CZERO,GTREF,NGD)
        CALL ZCOPY(NGD*LMAXSQ,GTREF,1,GREF(1,NLM2),1)
        IF (TEST('REFPOT  '))       
     +       write(6,*) N2,REFPOT(KAOEZ(ABS(ATOM(N2))))
 100  CONTINUE
        
      IF (TEST('WAIT    ')) write(6,*) 'Input I'
      IF (TEST('WAIT    ')) read(5,*) I

      CALL GREFSY(GREF,GREF0,NDIM)
      IF(TEST('flow    ')) write(6,*) 'GREFSY o.k.'
c
c test added 18.7.2000
c
cc      CALL GREFNEW(LSTART,EZ,CLEB,ICLEB,LOFLM,IEND,
cc     +             RATOM,NATOM,ALAT,GREF1)
c
c Now compare
c
c      IF (2.LT.1) THEN
c      DO N1=1,NATOM
c         DO LM1=1,LMAXSQ
c             DO LM2=1,LMAXSQ
c             NLM1 = (N1-1)*LMAXSQ + LM1
c             write(99,8000) nlm1,lm2, GREF0(NLM1,Lm2),
c     &                                GREF1(NLM1,Lm2)
c             END do
c         end do
c      end do 
c      end if
c 8000 format(2I5,4F15.8)
c
c test added 18.7.2000
c
  
      IF (OUT.GT.0) 
     +     WRITE (out) ez,((GREF0(N,M),M=1,LMAXSQ),N=1,LMAXSQ*NACLSD)
c small change 18.5.01 ez is added
c ------------------------------------------------------------------------
 999  RETURN
 9000 format(2f12.6)
      end




