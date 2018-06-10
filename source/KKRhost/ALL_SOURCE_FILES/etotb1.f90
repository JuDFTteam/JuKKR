! 17.10.95 ***************************************************************
SUBROUTINE etotb1(ecou,epotin,espc,espv,exc,kpre,lmax,lpot,  &
    lcoremax,nspin,natyp,nshell,conc,idoldau, lopt,eu,edcldau)
! ************************************************************************
!     calculate the total energy of the cluster .
!     gather all energy-parts which are calculated in different
!     subroutines .
!     since the program uses group theory only shell-indices
!     are used instead of atom-indices .

!                               b.drittler   may 1987

!     modified for supercells with nshell(i) atoms of type i in the
!     unit cell
!                               p.zahn       oct. 95

!     adopted for more atoms per site (CPA) v.popescu feb. 02
!-----------------------------------------------------------------------
use mod_types, only: t_inc
IMPLICIT NONE
INCLUDE 'inc.p'

! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************

! PARAMETER definitions
      INTEGER LMAXD1,NPOTD
      PARAMETER (LMAXD1=LMAXD+1,NPOTD=(2*KREL+(1-KREL)*NSPIND)*NATYPD)

!..Dummy arguments
INTEGER KPRE,LMAX,LPOT,NATYP,NSPIN,IDOLDAU
DOUBLE PRECISION CONC(NATYPD)
DOUBLE PRECISION ECOU(0:LPOTD,*),EPOTIN(*),ESPC(0:3,NPOTD), &
                 ESPV(0:LMAXD1,NPOTD),EXC(0:LPOTD,*), &
                 EU(*),EDCLDAU(*)
INTEGER LCOREMAX(*),NSHELL(*),LOPT(*)

!..Local variables
      DOUBLE PRECISION BANDESUM,BANDET,ECOUS,EDC,EFCTOR,ET,ETOT,EXCS
      DOUBLE PRECISION ETOTLDAU
      DOUBLE PRECISION DBLE
      INTEGER IATYP,IPOT,IS,ISPIN,L
      LOGICAL TEST
      CHARACTER*4 TEXTL(0:6)
      CHARACTER*5 TEXTNS
      CHARACTER*13 TEXTS(3)
!..
!.. externals
      EXTERNAL TEST
!..
!.. Data statements ..
      DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h =',' i ='/
      DATA TEXTS/' spin down   ',' spin  up    ',' paramagnetic'/
      DATA TEXTNS/' ns ='/
! ------------------------------------------------------------------------
efctor = 1.0D0/13.6058D0

etot = 0.0D0
bandesum = 0.0D0
etotldau = 0.0D0

IF ( kpre == 1 .AND. (t_inc%i_write>0)) WRITE (1337,FMT=99001)

!---> loop over host atoms

DO iatyp = 1,natyp
  
  IF ( kpre == 1 .AND. (t_inc%i_write>0)) WRITE (1337,FMT=99002) iatyp
  
  edc = 0.0D0
  et = 0.0D0
  bandet = 0.0D0
  ecous = 0.0D0
  excs = 0.0D0
  
  is = 0
  IF ( nspin == 1 ) is = is + 2
  DO ispin = 1,nspin
    is = is + 1
    ipot = (iatyp-1)*nspin + ispin
    
    IF ( kpre == 1 .AND. (t_inc%i_write>0)) THEN
      WRITE (1337,FMT=99003) texts(is)
      WRITE (1337,FMT=99004) (textl(l),espc(l,ipot),l=0, lcoremax(iatyp))
      WRITE (1337,FMT=99005) (textl(l),espv(l,ipot),l=0,lmax)
      WRITE (1337,FMT=99006) textns,espv(lmaxd1,ipot)
    END IF
    
    DO l = 0,lcoremax(iatyp)
      et = et + espc(l,ipot)
    END DO
    
    DO l = 0,lmax
      bandet = bandet + espv(l,ipot)
      et = et + espv(l,ipot)
    END DO
    bandet = bandet + espv(lmaxd1,ipot)
    et = et + espv(lmaxd1,ipot)
  END DO
  
! -> LDA+U
  
  et = et + eu(iatyp)
  bandet = bandet + eu(iatyp)
  IF ( kpre == 1 .AND. idoldau == 1 .AND. lopt(iatyp) >= 0  &
      .AND. (t_inc%i_write>0)) WRITE(1337,99019) eu(iatyp)
  
! --->  sum up Coulomb and Ex.-Corel. contribution
  
  DO l = 0,lpot
    ecous = ecous + ecou(l,iatyp)
    excs = excs + exc(l,iatyp)
  END DO
  
  IF ( kpre == 1 .AND. (t_inc%i_write>0)) THEN
    WRITE (1337,FMT=99007) et
    WRITE (1337,FMT=99008) bandet
    WRITE (1337,FMT=99009) (l,ecou(l,iatyp),l=0,lpot)
    WRITE (1337,FMT=99010)
    WRITE (1337,FMT=99018) ecous
    WRITE (1337,FMT=99011) (l,exc(l,iatyp),l=0,lpot)
    WRITE (1337,FMT=99010)
    WRITE (1337,FMT=99017) excs
    WRITE (1337,FMT=99015) epotin(iatyp)
  END IF
  
  IF ( .NOT.(test('NoMadel ')) ) THEN
    
    et = et + ecous + excs
    edc = edc + ecous + excs
    
    et = et + epotin(iatyp) - edcldau(iatyp)
    edc = edc + epotin(iatyp) - edcldau(iatyp)
    
    IF ( kpre == 1 .AND. (t_inc%i_write>0)) THEN
      IF ( idoldau == 1 .AND. lopt(iatyp) >= 0 )  &
          WRITE(1337,99020) -edcldau(iatyp)
      WRITE (1337,FMT=99016) edc
    END IF
    
  END IF
  
  IF ( natyp > 1 .OR. nshell(iatyp) > 1 ) THEN
    IF(t_inc%i_write>0) WRITE (1337,FMT=99012) iatyp,et
    IF ( kpre == 1 .AND. idoldau == 1 .AND. lopt(iatyp) >= 0  &
        .AND. (t_inc%i_write>0)) WRITE(1337,99021) eu(iatyp) - edcldau(iatyp)
    WRITE (1337,FMT=99022)
  END IF
  
  etot = etot + et*DBLE(nshell(iatyp))*conc(iatyp)
  bandesum = bandesum + bandet*DBLE(nshell(iatyp))*conc(iatyp)
  
END DO                        ! IATYP = 1,NATYP

IF(t_inc%i_write>0) WRITE (1337,FMT=99013) bandesum
IF(t_inc%i_write>0) WRITE (1337,FMT=99014) etot,etot/efctor
WRITE (*,FMT=99024) etot


RETURN

99001 FORMAT (32('='),' TOTAL ENERGIES ',31('='),/)
99002 FORMAT (3X,'Total energies atom ',i3,/,3X,23('-'))
99003 FORMAT (5X,'single particle energies ',a13)
99004 FORMAT (7X,'  core   contribution : ',2(a4,f15.8),/, (31X,2(a4,f15.8)))
99005 FORMAT (7X,'valence  contribution : ',2(a4,f15.8),/, (31X,2(a4,f15.8)))
99006 FORMAT (7X,'                        ',a4,f15.8)
99007 FORMAT (5X,68('-'),/,5X,  &
    'total contribution of the single particle energies :',1X, f15.8)
99008 FORMAT (5X,'                              band energy per atom :',  &
    1X,f15.10,/)
99009 FORMAT (5X,'coulomb  contribution : ',2(i3,1X,f15.8),/,  &
    (29X,2(i3,1X,f15.8)))
99010 FORMAT (5X,68('-'))

99011 FORMAT (5X,'ex.-cor. contribution : ',2(i3,1X,f15.8),/,  &
    (29X,2(i3,1X,f15.8)))
99012 FORMAT (/,3X,'Total contribution of atom',i3,' =',f15.8)
99013 FORMAT (5X,'                              sum of band energies :',  &
    1X,f15.10,/,3X,70('-'))
99014 FORMAT (/,3X,70('+'),/,15X,'TOTAL ENERGY in ryd. : ',f25.8,/15X,  &
    '                 eV  : ',f25.8,/,3X,70('+'))
99015 FORMAT (5X,'eff. pot. contribution     : ',f15.8)
99016 FORMAT (5X,'total double counting contribution                 :',  &
    1X,f15.8)
99017 FORMAT (5X,'tot. ex.-cor. contribution : ',f15.8,/)
99018 FORMAT (5X,'tot. coulomb contribution : ',f15.8,/)
99019 FORMAT (/,5X, 'LDA+U correction to the single particle energy     :',  &
    f16.8)
99020 FORMAT (/,5X, 'LDA+U double counting contribution                 :',  &
    f16.8)

99021 FORMAT (3X,'   including LDA+U correction :',f15.8)
99022 FORMAT (3X,70('-'))

99024 FORMAT ('TOTAL ENERGY in ryd. : ',f25.8,/15X )
END SUBROUTINE etotb1
