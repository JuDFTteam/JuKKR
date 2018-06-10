SUBROUTINE decipotbas(ihost,iqoff,itoff,nq,nt,rbasis,qmtet,qmphi,  &
    noq,kaoez,zat,iqat,conc,irws,ipan,ircut,rr,  &
    drdi,visp,nspin,krel,solver,socscl,cscl,  &
    vtrel,btrel,irmd,ipand,nembd1,ntmax,nspind, lmaxd)
! **********************************************************************
! *                                                                    *
! * reads in the potential data for the host atoms from the potential  *
! * file 'decimate.pot'                                                *
! *                                        v.popescu - munich, Dec 04  *
! *                                                                    *
! * Note: so far, only  SPHERICAL case implemented                     *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
!..
!.. Arguments
INTEGER IHOST,IPAND,IQOFF,IRMD,ITOFF,KREL,LMAXD,NEMBD1,NQ,NSPIN, &
        NSPIND,NT,NTMAX
CHARACTER*10 SOLVER
DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL),NTMAX),CONC(NTMAX), &
                 CSCL(KREL*LMAXD+1,KREL*NTMAX+(1-KREL)), &
                 DRDI(IRMD,NTMAX),QMPHI(NEMBD1),QMTET(NEMBD1), &
                 RBASIS(3,NEMBD1),RR(IRMD,NTMAX), &
                 SOCSCL(KREL*LMAXD+1,KREL*NTMAX+(1-KREL)), &
                 VISP(IRMD,NTMAX*NSPIND), &
                 VTREL(IRMD*KREL+(1-KREL),NTMAX),ZAT(NTMAX)
INTEGER IPAN(NTMAX),IQAT(NEMBD1,NTMAX),IRCUT(0:IPAND,NTMAX), &
        IRWS(NTMAX),KAOEZ(NEMBD1,NEMBD1),NOQ(NEMBD1)
!..
!.. Locals
      INTEGER I,IH,IHF,IL,IPOT1,IPOT2
      INTEGER NINT
      DOUBLE PRECISION RMT(NTMAX),RWS(NTMAX)
      CHARACTER*3 TXTT(NT)
! ......................................................................
! --> read basis

DO ih = 1,nq
  ihf = ih + iqoff
  READ (36+ihost,99001) il,(rbasis(i,ihf),i=1,3)
  IF ( ih /= il ) STOP ' Inconsistent data '
  WRITE (1337,99001) il,(rbasis(i,ihf),i=1,3)
END DO
READ (36+ihost,*)
WRITE (1337,99003)
DO ih = 1,nq
  ihf = ih + iqoff
  READ (36+ihost,FMT=99002) il,qmtet(ihf),qmphi(ihf),noq(ihf),  &
      (kaoez(i,ihf),i=1,noq(ihf))
  IF ( ih /= il ) STOP ' Inconsistent data '
  WRITE (1337,99004) ih,qmtet(ihf),qmphi(ihf),noq(ihf),  &
      (kaoez(i,ihf),i=1,noq(ihf))
END DO
WRITE (1337,99005)
IF ( krel == 1 ) READ (36+ihost,'(7X,A10)') solver
READ (36+ihost,*)

! --> read atoms

DO ih = 1,nt
  ihf = ih + itoff
  READ (36+ihost,99006) il,zat(ihf),iqat(1,ihf),conc(ihf),  &
      irws(ihf),ipan(ihf), (ircut(i,ihf),i=0,ipan(ihf))
  IF ( ih /= il ) STOP ' Inconsistent data '
  
  IF ( krel == 1 ) THEN
    READ (36+ihost,99007) socscl(1,ihf),cscl(1,ihf)
    DO il = 2,lmaxd + 1
      socscl(il,ihf) = socscl(1,ihf)
      cscl(il,ihf) = cscl(1,ihf)
    END DO
  END IF
  
END DO
DO ih = 1,nt
  ihf = ih + itoff
  ipot1 = (ihf-1)*nspin + 1
  ipot2 = ipot1 + 1
  READ (36+ihost,*)
  READ (36+ihost,99008) il,txtt(ih),rmt(ihf),rws(ihf)
  IF ( ih /= il ) STOP ' Inconsistent data '
  WRITE (1337,99009) ih,txtt(ih),nint(zat(ihf)),conc(ihf), irws(ihf),rws(ihf)
  IF ( krel == 0 ) THEN
    READ (36+ihost,*)
    DO i = 1,irws(ihf)
      READ (36+ihost,99010) rr(i,ihf),drdi(i,ihf), (visp(i,il),il=ipot1,ipot2)
    END DO
  ELSE
    READ (36+ihost,'(7X,I3)') il
    irws(ihf) = irws(ihf) - il
    ircut(ipan(ihf),ihf) = irws(ihf)
    READ (36+ihost,*)
    DO i = 1,irws(ihf)
      READ (36+ihost,99010) rr(i,ihf),drdi(i,ihf),vtrel(i,ihf), btrel(i,ihf)
    END DO
  END IF
END DO

99001 FORMAT (9X,i3,3F12.8)
99002 FORMAT (7X,i3,2(7X,f9.4),7X,i3,7X,8I3)
99003 FORMAT (9X,39('-'),/,9X,'   THETA   ','   PHI   ','OCC',' IT')
99004 FORMAT (9X,i3,2(f9.4),i3,8I3)
99005 FORMAT (9X,39('-'),/,10X,'ATOMS',/,15X,'Z   CONC  IWS    RWS')
99006 FORMAT (7X,i3,7X,f4.0,7X,i3,7X,f7.4,/,17X,i4,7X,i3,7X,6I4)
99007 FORMAT (17X,f10.6,7X,d13.6)
99008 FORMAT (7X,i3,1X,a3,2(/,7X,f12.8))
99009 FORMAT (9X,i3,1X,a3,i3,f7.4,i4,f10.6)
99010 FORMAT (1P,4D20.12)
END SUBROUTINE decipotbas
