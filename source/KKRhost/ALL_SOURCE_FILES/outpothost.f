C*==outpothost.f    processed by SPAG 6.05Rc at 13:41 on  3 Dec 2004
      SUBROUTINE OUTPOTHOST(ALAT,INS,KREL,KMROT,NSPIN,NAEZ,NATYP,EFERMI,
     &                      BRAVAIS,RBASIS,QMTET,QMPHI,NOQ,KAOEZ,IQAT,
     &                      ZAT,CONC,IPAN,IRCUT,SOLVER,SOC,CTL,IRWS,RMT,
     &                      RWS,RR,DRDI,VISP,IRSHIFT,RMREL,DRDIREL,
     &                      VTREL,BTREL,LMAXD,NATYPD,NAEZD,IPAND,IRMD)
C **********************************************************************
C *                                                                    *
C * writes decimation potential-file  'decimate.pot' to be later used  *
C * for 2D systems with the DECIMATE option. Based on the host         *
C * potentials, the single-site matrices of the host can be calculated *
C * directly on each particular energy-mesh                            *
C *                                        v.popescu - munich, Dec 04  *
C *                                                                    *
C * Note: so far, only SPHERICAL case implemented                      *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..      
C     .. Scalar arguments
      DOUBLE PRECISION ALAT,EFERMI
      INTEGER INS,IPAND,IRMD,KMROT,KREL,LMAXD,NAEZ
      INTEGER NAEZD,NATYP,NATYPD,NSPIN
      CHARACTER*10 SOLVER
C     ..
C     .. Array arguments
      DOUBLE PRECISION BRAVAIS(3,3),BTREL(IRMD*KREL+(1-KREL),*),CONC(*),
     &                 CTL(KREL*LMAXD+1,*),DRDI(IRMD,*),
     &                 DRDIREL(IRMD*KREL+(1-KREL),*),QMPHI(*),QMTET(*),
     &                 RBASIS(3,*),RMREL(IRMD*KREL+(1-KREL),*),RMT(*),
     &                 RR(IRMD,*),RWS(*),SOC(KREL*LMAXD+1,*),
     &                 VISP(IRMD,*),VTREL(IRMD*KREL+(1-KREL),*),ZAT(*)
      INTEGER IPAN(*),IQAT(*),IRCUT(0:IPAND,*),IRSHIFT(*),IRWS(*),
     &        KAOEZ(NATYPD,*),NOQ(NAEZD)
C     ..      
C     .. Locals
      CHARACTER*3 ELEMNAME(0:113)
      INTEGER I,IQ,IR,IS
      INTEGER INT
      CHARACTER*9 TXTREL(2),TXTSPIN(3)
C     ..
      DATA TXTSPIN/'         ','spin UP  ','spin DOWN'/
      DATA TXTREL/'(UP+DN)/2','(UP-DN)/2'/
C     .. 1     2     3     4    5     6     7     8     9     0
      DATA ELEMNAME/'Vac','H  ','He ','Li ','Be','B  ','C  ','N  ',
     &     'O  ','F  ','Ne ','Na ','Mg ','Al ','Si','P  ','S  ','Cl ',
     &     'Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr','Mn ','Fe ','Co ',
     &     'Ni ','Cu ','Zn ','Ga ','Ge ','As ','Se','Br ','Kr ','Rb ',
     &     'Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru','Rh ','Pd ','Ag ',
     &     'Cd ','In ','Sn ','Sb ','Te ','I  ','Xe','Cs ','Ba ','La ',
     &     'Ce ','Pr ','Nd ','Pm ','Sm ','Eu ','Gd','Tb ','Dy ','Ho ',
     &     'Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W ','Re ','Os ','Ir ',
     &     'Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po','At ','Rn ','Fr ',
     &     'Ra ','Ac ','Th ','Pa ','U  ','Np ','Pu','Am ','Cm ','Bk ',
     &     'Cf ','Es ','Fm ','Md ','No ','Lr ','Rf','Db ','Sg ','Bh ',
     &     'Hs ','Mt ','Uun','Uuu','Uub','NoE'/
C     ..
      WRITE (6,'(5X,A,A,/)') '< OUTPOTHOST > : ',
     &                     'creating decimate.pot file - host potential'
      OPEN (37,FILE='decimate.pot',STATUS='unknown')
      WRITE (37,FMT=*) 'Host structure and potential for decimation'
      WRITE (37,FMT=99001)
      WRITE (37,FMT=99005) KREL,INS,NSPIN,KMROT
      WRITE (37,FMT=99006) NAEZ,NATYP,ALAT
      WRITE (37,FMT=99004) EFERMI
      WRITE (37,FMT=99002) BRAVAIS
C ----------------------------------------------------------------------
C here insert whatever is more needed for the structure (BZ,SYM etc),
C ref. system and so on for a full host calculation (1 iteration)
C ----------------------------------------------------------------------
      WRITE (37,FMT=99003)
      DO IQ = 1,NAEZ
         WRITE (37,FMT=99007) IQ,(RBASIS(I,IQ),I=1,3)
      END DO
      WRITE (37,FMT=99008)
      DO IQ = 1,NAEZ
         WRITE (37,FMT=99009) IQ,QMTET(IQ),QMPHI(IQ),NOQ(IQ),
     &                        (KAOEZ(I,IQ),I=1,NOQ(IQ))
      END DO
      IF ( KREL.EQ.1 ) WRITE (37,99012) SOLVER
      WRITE (37,99010)
      DO I = 1,NATYP
         WRITE (37,FMT=99011) I,ZAT(I),IQAT(I),CONC(I),IRWS(I),
     &                        IPAN(I),(IRCUT(IQ,I),IQ=0,IPAN(I))
         IF ( KREL.EQ.1 ) WRITE (37,99013) SOC(1,I),CTL(1,I)
      END DO
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      DO I = 1,NATYP
         WRITE (37,'(80(1H*))')
         IR = INT(ZAT(I))
         WRITE (37,99014) I,ELEMNAME(IR),RMT(I),RWS(I)
         IF ( KREL.EQ.0 ) THEN
            WRITE (37,'(4(A9,11X))') 'R MESH   ','DRDI     ',
     &                               (TXTSPIN(NSPIN+IS-1),IS=1,NSPIN)
            DO IR = 1,IRWS(I)
               WRITE (37,99015) RR(IR,I),DRDI(IR,I),
     &                          (VISP(IR,(I-1)*NATYP+IS),IS=1,NSPIN)
            END DO
         ELSE
            WRITE (37,99016) IRSHIFT(I)
            WRITE (37,'(4(A9,11X))') 'R MESH   ','DRDI     ',
     &                               (TXTREL(IS),IS=1,2)
            DO IR = 1,IRWS(I) - IRSHIFT(I)
               WRITE (37,99015) RMREL(IR,I),DRDIREL(IR,I),VTREL(IR,I),
     &                          BTREL(IR,I)
            END DO
         END IF
      END DO
C
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CLOSE (37)
C     ..
99001 FORMAT ('Vectors in lattice constant units')
99002 FORMAT ('BRAVAIS ',/,3F8.4,/,3F8.4,/,3F8.4)
99003 FORMAT ('RBASIS')
99004 FORMAT ('EFERMI=',F10.6)
99005 FORMAT ('KREL  =',I3,' INS  =',I3,' NSPIN=',I3,' KMROT=',I3)
99006 FORMAT ('NAEZ  =',I3,' NATYP=',I3,' ALAT =',F12.8)
99007 FORMAT ('SITE  :',I3,3F12.8)
99008 FORMAT ('Magnetisation angles, occupancies, types on each site')
99009 FORMAT ('SITE  :',I3,' THETA=',F9.4,' PHI  =',F9.4,' NOQ  =',I3,
     &        ' ITOQ :',8I3)
99010 FORMAT ('ATOMS')
99011 FORMAT ('TYPE  :',I3,' Z    =',F4.0,' IQAT =',I3,' CONC =',F7.4,/,
     &        10X,' IRWS =',I4,' IPAN =',I3,' IRCUT=',6I4)
99012 FORMAT ('SOLVER=',A10)
99013 FORMAT (10X,' SOC  =',F10.6,' CTL  =',D13.6)
99014 FORMAT ('ATOM  :',I3,1X,A3,': mesh and potential data',/,
     &        'RMT   :',F12.8,/,'RWS   :',F12.8)
99015 FORMAT (1P,4D20.12)
99016 FORMAT ('ISHIFT:',I3)
      END
