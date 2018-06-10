!*==outpothost.f    processed by SPAG 6.05Rc at 13:41 on  3 Dec 2004
SUBROUTINE outpothost(alat,ins,krel,kmrot,nspin,naez,natyp,efermi,  &
    bravais,rbasis,qmtet,qmphi,noq,kaoez,iqat,  &
    zat,conc,ipan,ircut,solver,soc,ctl,irws,rmt,  &
    rws,rr,drdi,visp,irshift,rmrel,drdirel,  &
    vtrel,btrel,lmaxd,natypd,naezd,ipand,irmd)
! **********************************************************************
! *                                                                    *
! * writes decimation potential-file  'decimate.pot' to be later used  *
! * for 2D systems with the DECIMATE option. Based on the host         *
! * potentials, the single-site matrices of the host can be calculated *
! * directly on each particular energy-mesh                            *
! *                                        v.popescu - munich, Dec 04  *
! *                                                                    *
! * Note: so far, only SPHERICAL case implemented                      *
! *                                                                    *
! **********************************************************************
      use mod_version_info
      IMPLICIT NONE
!..      
!.. Scalar arguments
      DOUBLE PRECISION ALAT,EFERMI
      INTEGER INS,IPAND,IRMD,KMROT,KREL,LMAXD,NAEZ
      INTEGER NAEZD,NATYP,NATYPD,NSPIN
      CHARACTER*10 SOLVER
!..
!.. Array arguments
DOUBLE PRECISION BRAVAIS(3,3),BTREL(IRMD*KREL+(1-KREL),*),CONC(*), &
                 CTL(KREL*LMAXD+1,*),DRDI(IRMD,*), &
                 DRDIREL(IRMD*KREL+(1-KREL),*),QMPHI(*),QMTET(*), &
                 RBASIS(3,*),RMREL(IRMD*KREL+(1-KREL),*),RMT(*), &
                 RR(IRMD,*),RWS(*),SOC(KREL*LMAXD+1,*), &
                 VISP(IRMD,*),VTREL(IRMD*KREL+(1-KREL),*),ZAT(*)
INTEGER IPAN(*),IQAT(*),IRCUT(0:IPAND,*),IRSHIFT(*),IRWS(*), &
        KAOEZ(NATYPD,*),NOQ(NAEZD)
!..      
!.. Locals
      CHARACTER*3 ELEMNAME(0:113)
      INTEGER I,IQ,IR,IS
      INTEGER INT
      CHARACTER*9 TXTREL(2),TXTSPIN(3)
!..
      DATA TXTSPIN/'         ','spin UP  ','spin DOWN'/
      DATA TXTREL/'(UP+DN)/2','(UP-DN)/2'/
!.. 1     2     3     4    5     6     7     8     9     0
DATA ELEMNAME/'Vac','H  ','He ','Li ','Be','B  ','C  ','N  ',        &      
     'O  ','F  ','Ne ','Na ','Mg ','Al ','Si','P  ','S  ','Cl ',     &
     'Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr','Mn ','Fe ','Co ',     &
     'Ni ','Cu ','Zn ','Ga ','Ge ','As ','Se','Br ','Kr ','Rb ',     &
     'Sr ','Y  ','Zr ','Nb ','Mo ','Tc ','Ru','Rh ','Pd ','Ag ',     &
     'Cd ','In ','Sn ','Sb ','Te ','I  ','Xe','Cs ','Ba ','La ',     &
     'Ce ','Pr ','Nd ','Pm ','Sm ','Eu ','Gd','Tb ','Dy ','Ho ',     &
     'Er ','Tm ','Yb ','Lu ','Hf ','Ta ','W ','Re ','Os ','Ir ',     &
     'Pt ','Au ','Hg ','Tl ','Pb ','Bi ','Po','At ','Rn ','Fr ',     &
     'Ra ','Ac ','Th ','Pa ','U  ','Np ','Pu','Am ','Cm ','Bk ',     &
     'Cf ','Es ','Fm ','Md ','No ','Lr ','Rf','Db ','Sg ','Bh ',     &
     'Hs ','Mt ','Uun','Uuu','Uub','NoE'/ 
!     ..
WRITE (1337,'(5X,A,A,/)') '< OUTPOTHOST > : ',  &
    'creating decimate.pot file - host potential'

OPEN (37,FILE='decimate.pot',STATUS='unknown')
CALL version_print_header(37)
WRITE (37,FMT=*) 'Host structure and potential for decimation'
WRITE (37,FMT=99001)
WRITE (37,FMT=99005) krel,ins,nspin,kmrot
WRITE (37,FMT=99006) naez,natyp,alat
WRITE (37,FMT=99004) efermi
WRITE (37,FMT=99002) bravais
! ----------------------------------------------------------------------
! here insert whatever is more needed for the structure (BZ,SYM etc),
! ref. system and so on for a full host calculation (1 iteration)
! ----------------------------------------------------------------------
WRITE (37,FMT=99003)
DO iq = 1,naez
  WRITE (37,FMT=99007) iq,(rbasis(i,iq),i=1,3)
END DO
WRITE (37,FMT=99008)
DO iq = 1,naez
  WRITE (37,FMT=99009) iq,qmtet(iq),qmphi(iq),noq(iq),  &
      (kaoez(i,iq),i=1,noq(iq))
END DO
IF ( krel == 1 ) WRITE (37,99012) solver
WRITE (37,99010)
DO i = 1,natyp
  WRITE (37,FMT=99011) i,zat(i),iqat(i),conc(i),irws(i),  &
      ipan(i),(ircut(iq,i),iq=0,ipan(i))
  IF ( krel == 1 ) WRITE (37,99013) soc(1,i),ctl(1,i)
END DO
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
DO i = 1,natyp
  WRITE (37,'(80(1H*))')
  ir = INT(zat(i))
  WRITE (37,99014) i,elemname(ir),rmt(i),rws(i)
  IF ( krel == 0 ) THEN
    WRITE (37,'(4(A9,11X))') 'R MESH   ','DRDI     ',  &
        (txtspin(nspin+is-1),is=1,nspin)
    DO ir = 1,irws(i)
      WRITE (37,99015) rr(ir,i),drdi(ir,i),  &
          (visp(ir,(i-1)*natyp+is),is=1,nspin)
    END DO
  ELSE
    WRITE (37,99016) irshift(i)
    WRITE (37,'(4(A9,11X))') 'R MESH   ','DRDI     ', (txtrel(is),is=1,2)
    DO ir = 1,irws(i) - irshift(i)
      WRITE (37,99015) rmrel(ir,i),drdirel(ir,i),vtrel(ir,i), btrel(ir,i)
    END DO
  END IF
END DO

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
CLOSE (37)
!     ..
99001 FORMAT ('Vectors in lattice constant units')
99002 FORMAT ('BRAVAIS ',/,3F8.4,/,3F8.4,/,3F8.4)
99003 FORMAT ('RBASIS')
99004 FORMAT ('EFERMI=',f10.6)
99005 FORMAT ('KREL  =',i3,' INS  =',i3,' NSPIN=',i3,' KMROT=',i3)
99006 FORMAT ('NAEZ  =',i3,' NATYP=',i3,' ALAT =',f12.8)
99007 FORMAT ('SITE  :',i3,3F12.8)
99008 FORMAT ('Magnetisation angles, occupancies, types on each site')
99009 FORMAT ('SITE  :',i3,' THETA=',f9.4,' PHI  =',f9.4,' NOQ  =',i3,  &
    ' ITOQ :',8I3)
99010 FORMAT ('ATOMS')
99011 FORMAT ('TYPE  :',i3,' Z    =',f4.0,' IQAT =',i3,' CONC =',f7.4,/,  &
    10X,' IRWS =',i4,' IPAN =',i3,' IRCUT=',6I4)
99012 FORMAT ('SOLVER=',a10)
99013 FORMAT (10X,' SOC  =',f10.6,' CTL  =',d13.6)
99014 FORMAT ('ATOM  :',i3,1X,a3,': mesh and potential data',/,  &
    'RMT   :',f12.8,/,'RWS   :',f12.8)
99015 FORMAT (1P,4D20.12)
99016 FORMAT ('ISHIFT:',i3)
END SUBROUTINE outpothost
