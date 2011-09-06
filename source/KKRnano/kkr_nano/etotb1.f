      SUBROUTINE ETOTB1(ECOU,EFERMI,EPOTIN,ESPC,ESPV,EXC,
     &                  EULDAU,EDCLDAU,LDAU,
     &                  KPRE,LMAX,LPOT,
     &                  LCOREMAX,NSPIN,I1,NAEZ)
C ************************************************************************
C     calculate the total energy of the cluster .
C     gather all energy-parts which are calculated in different
C     subroutines .
C     since the program uses group theory only shell-indices
C     are used instead of atom-indices .
C
C                               b.drittler   may 1987
C
C     modified for supercells with nshell(i) atoms of type i in the
C     unit cell
C                               p.zahn       oct. 95
C
C
C-----------------------------------------------------------------------
C
C  MAIN2
C    | 
C	 +- RESULTS
C    |     |
C    |     +-ETOTB1
C
C-----------------------------------------------------------------------
      IMPLICIT NONE

      INCLUDE 'inc.p'
C
      INTEGER LMAXD1
      PARAMETER (LMAXD1=LMAXD+1)
C
C Dummy arguments
C
      DOUBLE PRECISION EFERMI
      INTEGER KPRE,LMAX,LPOT,NAEZ,NSPIN
      DOUBLE PRECISION ECOU(0:LPOTD),EPOTIN,ESPC(0:3,NSPIND),
     &                 ESPV(0:LMAXD1,NSPIND),EXC(0:LPOTD),
     &                 EULDAU,EDCLDAU
      INTEGER LCOREMAX
C
C Local variables
C
      DOUBLE PRECISION BANDESUM,BANDET,ECOUS,EDC,EFCTOR,ET,ETOT,EXCS
      DOUBLE PRECISION DBLE
      INTEGER I1,IPOT,IS,ISPIN,L
      LOGICAL TEST,LDAU
      CHARACTER*4 TEXTL(0:6)
      CHARACTER*5 TEXTNS
      CHARACTER*13 TEXTS(3)
C     ..
C     .. externals
      EXTERNAL TEST
C     ..
      SAVE ETOT,BANDESUM
C     .. Data statements ..
      DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h =',' i ='/
      DATA TEXTS/' spin down   ',' spin  up    ',' paramagnetic'/
      DATA TEXTNS/' ns ='/
C ------------------------------------------------------------------------
      EFCTOR = 1.0D0/13.6058D0
C
      IF(I1.EQ.1) THEN
        ETOT = 0.0D0
        BANDESUM = 0.0D0
C
        IF ( KPRE.EQ.1 ) WRITE (6,FMT=99001)
      END IF
C
         IF ( KPRE.EQ.1 ) WRITE (6,FMT=99002) I1
C
         EDC = 0.0D0
         ET = 0.0D0
         BANDET = 0.0D0
         ECOUS = 0.0D0
         EXCS = 0.0D0
C
         IS = 0
         IF ( NSPIN.EQ.1 ) IS = IS + 2
         DO ISPIN = 1,NSPIN
            IS = IS + 1
            IPOT = (I1-1)*NSPIN + ISPIN
C
            IF ( KPRE.EQ.1 ) THEN
               WRITE (6,FMT=99003) TEXTS(IS)
               WRITE (6,FMT=99004) (TEXTL(L),ESPC(L,ISPIN),L=0,LCOREMAX)
               WRITE (6,FMT=99005) (TEXTL(L),ESPV(L,ISPIN),L=0,LMAX)
               WRITE (6,FMT=99006) TEXTNS,ESPV(LMAXD1,ISPIN)
            END IF
C
            DO L = 0,LCOREMAX
               ET = ET + ESPC(L,ISPIN)
            END DO
C
            DO L = 0,LMAX
               BANDET = BANDET + ESPV(L,ISPIN)
               ET = ET + ESPV(L,ISPIN)
            END DO
            BANDET = BANDET + ESPV(LMAXD1,ISPIN)
            ET = ET + ESPV(LMAXD1,ISPIN)
         END DO
C
         ET = ET + EULDAU
         BANDET = BANDET + EULDAU
         IF ((LDAU).AND.(EULDAU.NE.0.0D0)) WRITE(6,99019) EULDAU
C
C --->  sum up Coulomb and Ex.-Corel. contribution
C
         DO L = 0,LPOT
            ECOUS = ECOUS + ECOU(L)
            EXCS = EXCS + EXC(L)
         END DO
C
         IF ( KPRE.EQ.1 ) THEN
            WRITE (6,FMT=99007) ET
            WRITE (6,FMT=99008) BANDET
            WRITE (6,FMT=99009) (L,ECOU(L),L=0,LPOT)
            WRITE (6,FMT=99010)
            WRITE (6,FMT=99018) ECOUS
            WRITE (6,FMT=99011) (L,EXC(L),L=0,LPOT)
            WRITE (6,FMT=99010)
            WRITE (6,FMT=99017) EXCS
            WRITE (6,FMT=99015) EPOTIN
         END IF
C
C
         ET = ET + ECOUS + EXCS
         EDC = EDC + ECOUS + EXCS
C
         ET = ET + EPOTIN - EDCLDAU
         EDC = EDC + EPOTIN - EDCLDAU
         IF ((LDAU).AND.(EDCLDAU.NE.0.0D0)) WRITE(6,99020) -EDCLDAU
C
         IF ( KPRE.EQ.1 ) THEN
           WRITE (6,FMT=99016) EDC
           IF ((LDAU).AND.(EULDAU-EDCLDAU).NE.0.0D0)
     +     WRITE (6,FMT=99021) (EULDAU-EDCLDAU)
         ENDIF
C
         IF ( NAEZ.GT.1) WRITE (6,FMT=99012)
     &        I1,ET
C
          ETOT = ETOT + ET
          BANDESUM = BANDESUM + BANDET
C
      IF (I1.EQ.NAEZ) THEN 
      WRITE (6,FMT=99013) BANDESUM
      WRITE (6,FMT=99014) ETOT,ETOT/EFCTOR
      END IF
C
      RETURN
C

99001 FORMAT (32('='),' TOTAL ENERGIES ',31('='),/)
99002 FORMAT (3X,'Total energies atom ',i5,/,3x,23('-'))
99003 FORMAT (5X,'single particle energies ',a13)
99004 FORMAT (7X,'  core   contribution : ',2(a4,f15.8),/,
     &        (31X,2(a4,f15.8)))
99005 FORMAT (7X,'valence  contribution : ',2(a4,f15.8),/,
     &        (31X,2(a4,f15.8)))
99006 FORMAT (7x,'                        ',a4,f15.8)
99007 FORMAT (5X,68('-'),/,5X,
     &        'total contribution of the single particle energies :',1X,
     &        f15.8)
99008 FORMAT (5X,'                              band energy per atom :',
     &        1X,f15.10,/)
99009 FORMAT (5X,'coulomb  contribution : ',2(i3,1X,f15.8),/,
     &        (29X,2(i5,1X,f15.8)))
99010 FORMAT (5X,68('-'))
C
99011 FORMAT (5X,'ex.-cor. contribution : ',2(i3,1X,f15.8),/,
     &        (29X,2(i5,1X,f15.8)))
99012 FORMAT (3x,'Total contribution of atom',i5,' =',f15.8)
99013 FORMAT (5X,'        sum of band energies :',
     &        1X,F20.10,/,3X,70('-'))
99014 FORMAT (/,3X,70('+'),/,15x,'TOTAL ENERGY in ryd. : ',f21.8,/15x,
     &        '                 eV  : ',F21.8,/,3x,70('+'))
99015 FORMAT (5X,'eff. pot. contribution     : ',f15.8)
99016 FORMAT (5X,'total double counting contribution                 :',
     &        1X,f15.8)
99017 FORMAT (5X,'tot. ex.-cor. contribution : ',f15.8,/)
99018 FORMAT (5X,'tot. coulomb contribution : ',f15.8,/)
99019 FORMAT (5X,'LDA+U correction to the single particle energy   :',
     &        F16.8)
99020 FORMAT (5X,'LDA+U double counting contribution               :',
     &        F16.8)
99021 FORMAT (3X,'   including LDA+U correction :',F15.8)
C
      END
