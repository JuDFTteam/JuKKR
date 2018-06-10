SUBROUTINE scfiterang(itrscf,itoq,fact,mvphi,mvtet,mvgam,  &
    qmphi,qmtet,qmgam,nq,nk,erravang, nqmax,ntmax,nmvecmax,nkmmax)
!   ********************************************************************
!   *                                                                  *
!   * applied to mix the angles specifying                             *
!   * the local frame of reference                                     *
!   *                                                                  *
!   ********************************************************************
use mod_wunfiles, only: t_params
IMPLICIT NONE

! PARAMETER definitions

INTEGER IXTRMAX
PARAMETER (IXTRMAX=4)
DOUBLE PRECISION DANGMAX
PARAMETER (DANGMAX=3D0)

! Dummy arguments

DOUBLE PRECISION ERRAVANG
INTEGER ITRSCF,NK,NKMMAX,NMVECMAX,NQ,NQMAX,NTMAX
DOUBLE PRECISION FACT(0:100),MVGAM(NTMAX,NMVECMAX), &
       MVPHI(NTMAX,NMVECMAX),MVTET(NTMAX,NMVECMAX), &
       QMGAM(NQMAX),QMPHI(NQMAX),QMTET(NQMAX)
INTEGER ITOQ(NTMAX,NQMAX)

! Local variables

DOUBLE PRECISION A,B,C,PHIXTR,TETXTR, D12, D23, D3X
DOUBLE PRECISION DELPHI,DELTET,LASTERR,MIXING,QMGAMMIX, &
       QMPHIMIX,QMTETMIX,WN,WO
INTEGER I,IMV,IPREV,IPRINT,IQ,IT,ITAB,IXTR
DOUBLE PRECISION QMGAMTAB(NQMAX,3),QMPHITAB(NQMAX,3), &
                 QMTETTAB(NQMAX,3)
DOUBLE COMPLEX DROTQ(NKMMAX,NKMMAX,NQMAX)

mixing = 1.0D0
wo = 1.0D0 - mixing
wn = mixing
iprint = 0
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ITERMDIR I/O

qmtet    = t_params%qmtet
qmphi    = t_params%qmphi
qmphitab = t_params%qmphitab
qmtettab = t_params%qmtettab
qmgamtab = t_params%qmgamtab
itab     = t_params%itab
lasterr  = t_params%lasterr
!       READ (67) QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ITAB,LASTERR
!       REWIND (67)

!      ITERMDIR I/O
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

WRITE (1337,FMT=99002)

IF ( itrscf == 0 ) THEN
  
! ----------------------------------------------------- store old angles
! --------------------------------------- dummy for tbkkr, done in main0
  DO iq = 1,nq
    qmphitab(iq,1) = qmphi(iq)
    qmtettab(iq,1) = qmtet(iq)
    qmgamtab(iq,1) = qmgam(iq)
    DO i = 2,3
      qmphitab(iq,i) = 0D0
      qmtettab(iq,i) = 0D0
      qmgamtab(iq,i) = 0D0
    END DO
  END DO
  
  RETURN
  
ELSE
  
!------------------------- set new local frame of reference according to
!-------------------------------------------- orientation of spin moment
  
  imv = 1
  
  DO iq = 1,nq
    
    it = itoq(1,iq)
    
    qmphi(iq) = mvphi(it,imv)
    qmtet(iq) = mvtet(it,imv)
    qmgam(iq) = mvgam(it,imv)
    
  END DO
  
!=======================================================================
  
  IF ( itrscf <= 2 ) THEN
    itab = itrscf
    DO iq = 1,nq
      qmphi(iq) = wo*qmphitab(iq,itab) + wn*qmphi(iq)
      qmtet(iq) = wo*qmtettab(iq,itab) + wn*qmtet(iq)
      qmgam(iq) = wo*qmgamtab(iq,itab) + wn*qmgam(iq)
    END DO
    
    itab = itrscf + 1
    DO iq = 1,nq
      
      qmphitab(iq,itab) = qmphi(iq)
      qmtettab(iq,itab) = qmtet(iq)
      qmgamtab(iq,itab) = qmgam(iq)
    END DO
    
    iprev = itab - 1
    
  ELSE
    IF ( lasterr > 2D0 ) THEN
      ixtr = 2
    ELSE
      ixtr = ixtrmax
    END IF
    
    DO iq = 1,nq
      
      qmphimix = wo*qmphitab(iq,3) + wn*qmphi(iq)
      qmtetmix = wo*qmtettab(iq,3) + wn*qmtet(iq)
      qmgammix = wo*qmgamtab(iq,3) + wn*qmgam(iq)
      
      a = qmphitab(iq,2)
      b = (qmphitab(iq,3)-qmphitab(iq,1))*0.5D0
      c = (qmphitab(iq,3)+qmphitab(iq,1)-2.d0*a)*0.5D0
      phixtr = a + b*DBLE(ixtr) + c*DBLE(ixtr)**2
      IF ( phixtr >= 0.0D0 .AND. phixtr <= 360.0D0 .AND.  &
            ABS(phixtr-qmphitab(iq,3)) < dangmax ) THEN
        qmphi(iq) = phixtr
      ELSE
        qmphi(iq) = qmphimix
      END IF
      
      IF ( qmphi(iq) < 0D0 ) qmphi(iq) = qmphi(iq) + 360D0
      IF ( qmphi(iq) > 360D0 ) qmphi(iq) = qmphi(iq) - 360D0
      
      a = qmtettab(iq,2)
      b = (qmtettab(iq,3)-qmtettab(iq,1))*0.5D0
      c = (qmtettab(iq,3)+qmtettab(iq,1)-2.d0*a)*0.5D0
      tetxtr = a + b*DBLE(ixtr) + c*DBLE(ixtr)**2
      d12=qmtettab(iq,1)-qmtettab(iq,2)
      d23=qmtettab(iq,2)-qmtettab(iq,3)
      d3x=qmtettab(iq,3)-tetxtr
      IF ( tetxtr >= 0.0D0 .AND. tetxtr <= 180.0D0 .AND.  &
            ABS(tetxtr-qmtettab(iq,3)) < dangmax .AND.  &
            (d12*d23 > 0D0) .AND. (d23*d3x > 0D0) ) THEN
        qmtet(iq) = tetxtr
      ELSE
        qmtet(iq) = qmtetmix
      END IF
      
      a = qmgamtab(iq,2)
      b = (qmgamtab(iq,3)-qmgamtab(iq,1))*0.5D0
      c = (qmgamtab(iq,3)+qmgamtab(iq,1)-2.d0*a)*0.5D0
      qmgam(iq) = a + b*DBLE(ixtr) + c*DBLE(ixtr)**2
      
      IF ( iprint > 0 ) WRITE (1337,99001) iq,ixtr,lasterr,  &
          ('PHI',i,qmphitab(iq,i),i=1,3),'PHIMIX',qmphimix,  &
          'PHIXTR',phixtr,'PHINEW',qmphi(iq),  &
          ('TET',i,qmtettab(iq,i),i=1,3),'TETMIX',qmtetmix,  &
          'TETXTR',tetxtr,'TETNEW',qmtet(iq)
      
      DO i = 1,2
        qmphitab(iq,i) = qmphitab(iq,i+1)
        qmtettab(iq,i) = qmtettab(iq,i+1)
        qmgamtab(iq,i) = qmgamtab(iq,i+1)
      END DO
      qmphitab(iq,3) = qmphi(iq)
      qmtettab(iq,3) = qmtet(iq)
      qmgamtab(iq,3) = qmgam(iq)
      
    END DO
    
    iprev = 2
  END IF
  
  WRITE (1337,99004)
  
  erravang = 0D0
  DO iq = 1,nq
    
    delphi = ABS(qmphi(iq)-qmphitab(iq,iprev))
    deltet = ABS(qmtet(iq)-qmtettab(iq,iprev))
    erravang = MAX(erravang,delphi,deltet)
    
    WRITE (1337,99005) iq,qmphitab(iq,iprev),qmtettab(iq,iprev),  &
        qmphi(iq),qmtet(iq),delphi,deltet
    
! --> update the rotation matrices DROTQ for the new angles
    
    CALL calcrotmat(nk,3,qmphi(iq),qmtet(iq),0.0D0, drotq(1,1,iq),fact,nkmmax)
    
  END DO
  
  WRITE (1337,FMT=99003) itrscf,erravang
  WRITE (1337,'(I5,4F10.3,'' #  ANGLES'',/,79(1H+),/)') itrscf,  &
      (qmphi(iq),qmtet(iq),iq=1,MIN(2,nq))
  
  lasterr = erravang
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!      ITERMDIR I/O
  
!          WRITE (67) QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ITAB,LASTERR
!          WRITE (67) DROTQ
  t_params%qmtet    = qmtet
  t_params%qmphi    = qmphi
  t_params%qmphitab = qmphitab
  t_params%qmtettab = qmtettab
  t_params%qmgamtab = qmgamtab
  t_params%itab     = itab
  t_params%lasterr  = lasterr
  t_params%drotq    = drotq
  
!      ITERMDIR I/O
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
END IF
99001 FORMAT (/,5X,'IQ',i3,'  IXTR',i3,'  LAST ERROR',f12.4,/,  &
    3(3(5X,a,5X,'TAB',i2,2X,f12.6,/),3(5X,a,9X,f12.6,/)))
!   ====================================================================
99002 FORMAT (79('+'),/,28X,'ANGLE - mixing scheme')
99003 FORMAT (5X,'iter.',i4,'     max. CHANGE = ',f12.6)
99004 FORMAT (/,5X,' setting new LOCAL frame of reference ',//,5X,  &
    'IQ  old   phi      tet    --> new   phi      tet', '     del   phi      tet')
99005 FORMAT (5X,i2,4X,2F9.4,8X,2F9.4,5X,2F9.4)
END SUBROUTINE scfiterang
