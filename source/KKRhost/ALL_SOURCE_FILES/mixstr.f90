! 13.10.95 ***************************************************************
SUBROUTINE mixstr(rmsavq,rmsavm,ins,lpot,lmpot,natref,nshell,  &
    nstart,nend,conc,nspin,itc,rfpi,fpi,ipf,  &
    mixing,fcm,irc,irmin,r,drdi,vons,visp,vins, vspsmo,vspsme,lsmear)
! ************************************************************************
!.. Parameters ..
      include 'inc.p'
      INTEGER LMPOTD,IRMIND
      PARAMETER (LMPOTD= (LPOTD+1)**2, &
                IRMIND=IRMD-IRNSD)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION FCM,FPI,MIXING,RFPI,RMSAVM,RMSAVQ
      INTEGER INS,IPF,ITC,LMPOT,LPOT,NATREF,NEND,NSPIN,NSTART
      INTEGER LSMEAR
!..
!.. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,NATYPD),R(IRMD,NATYPD), &
                       VINS(IRMIND:IRMD,LMPOTD,*),VISP(IRMD,*), &
                       VONS(IRMD,LMPOTD,*),CONC(NATYPD), &
                       VSPSMO(IRMD,NSPOTD),VSPSME(IRMD,NSPOTD)
      INTEGER IRC(NATYPD),IRMIN(NATYPD),NSHELL(0:NSHELD)
!..
!.. Local Scalars ..
      DOUBLE PRECISION FAC,RMSERM,RMSERQ,VMN,VNM,VNP,VOLDM,VOLDP,VPN
      DOUBLE PRECISION NATOM
      INTEGER I,IH,IHP1,IRC1,IRMIN1,J,LM,NP
!..
!.. Intrinsic Functions ..
      INTRINSIC MOD,REAL,SQRT
!     ..
rmsavq = 0.0D0
rmsavm = 0.0D0

!---> final construction of the potentials
!     attention : the spherical averaged potential is the lm=1
!                     component of vons times sqrt(4 pi).

!     first mixing scheme : straight mixing
!---> determination of the root mean sqare error

natom = 0.0D0

DO  np = nstart,nend
  
  i = np - natref
  natom = natom + DBLE(nshell(i))*conc(i)
  
  IF (nspin == 2) THEN
    ih = 2*np - 1
    ihp1 = ih + 1
  ELSE
    ih = np
    ihp1 = ih
  END IF
  
  irc1 = irc(np)
  rmserq = 0.0D0
  rmserm = 0.0D0
  fac = 0.5D0/rfpi
  
  DO  j = 1,irc1
    vnp = fac* (vons(j,1,ih)+vons(j,1,ihp1))
    vnm = fac* (vons(j,1,ih)-vons(j,1,ihp1))
    voldp = 0.5D0* (visp(j,ih)+visp(j,ihp1))
    voldm = 0.5D0* (visp(j,ih)-visp(j,ihp1))
    rmserq = rmserq + 2.0D0*REAL(1+MOD(j,2))*r(j,np)*r(j,np)*  &
        drdi(j,np)* (vnp-voldp)* (vnp-voldp)
    rmserm = rmserm + 2.0D0*REAL(1+MOD(j,2))*r(j,np)*r(j,np)*  &
        drdi(j,np)* (vnm-voldm)* (vnm-voldm)
    vpn = voldp + mixing* (vnp-voldp)
    vmn = voldm + fcm*mixing* (vnm-voldm)
    vons(j,1,ihp1) = vpn - vmn
    vons(j,1,ih) = vpn + vmn
  END DO
  
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  IF ( lsmear >= 3 ) THEN
    DO j = 1,irc1
      vnp   = 0.5D0* (vspsmo(j,ih)+vspsmo(j,ihp1))
      vnm   = 0.5D0* (vspsmo(j,ih)-vspsmo(j,ihp1))
      voldp = 0.5D0* (vspsme(j,ih)+vspsme(j,ihp1))
      voldm = 0.5D0* (vspsme(j,ih)-vspsme(j,ihp1))
      vpn   = voldp + mixing* (vnp-voldp)
      vmn   = voldm + fcm*mixing* (vnm-voldm)
      vspsmo(j,ihp1) = vpn - vmn
      vspsmo(j,ih)   = vpn + vmn
    END DO
  END IF
  
  IF ( (lsmear == 1) .OR. (lsmear == 2) ) THEN
    DO j = 1,irc1
      vspsme(j,ihp1) = vspsmo(j,ihp1)
      vspsme(j,ih)   = vspsmo(j,ih)
    END DO
  END IF
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  
  rmserq = rmserq/ (r(irc1,np)**3)
  rmserm = rmserm/ (r(irc1,np)**3)
  rmsavq = rmsavq + rmserq*nshell(i)*conc(i)
  rmsavm = rmsavm + rmserm*nshell(i)*conc(i)
  
  IF (nspin == 2) THEN
    WRITE (ipf,FMT=9000) i,SQRT(rmserq),SQRT(rmserm)
  ELSE
    WRITE (ipf,FMT=9020) i,SQRT(rmserq)
  END IF
  
  IF (ins /= 0 .AND. lpot > 0) THEN
    
    rmserq = 0.0D0
    rmserm = 0.0D0
    irmin1 = irmin(np)
    DO  lm = 2,lmpot
      DO  j = irmin1,irc1
        vnp = 0.5D0* (vons(j,lm,ih)+vons(j,lm,ihp1))
        vnm = 0.5D0* (vons(j,lm,ih)-vons(j,lm,ihp1))
        voldp = 0.5D0* (vins(j,lm,ih)+vins(j,lm,ihp1))
        voldm = 0.5D0* (vins(j,lm,ih)-vins(j,lm,ihp1))
        rmserq = rmserq + 2.0D0*REAL(1+MOD(j,2))*r(j,np)*r(j,np)*  &
            drdi(j,np)* (vnp-voldp)* (vnp-voldp)
        rmserm = rmserm + 2.0D0*REAL(1+MOD(j,2))*r(j,np)*r(j,np)*  &
            drdi(j,np)* (vnm-voldm)* (vnm-voldm)
        vpn = voldp + mixing* (vnp-voldp)
        vmn = voldm + fcm*mixing* (vnm-voldm)
        vons(j,lm,ihp1) = vpn - vmn
        vons(j,lm,ih) = vpn + vmn
      END DO
    END DO
    rmserq = rmserq/ (r(irc1,np)**3)/fpi
    rmserm = rmserm/ (r(irc1,np)**3)/fpi
    rmsavq = rmsavq + rmserq*nshell(i)*conc(i)
    rmsavm = rmsavm + rmserm*nshell(i)*conc(i)
    
    IF (nspin == 2) THEN
      WRITE (ipf,FMT=9010) i,SQRT(rmserq),SQRT(rmserm)
    ELSE
      WRITE (ipf,FMT=9030) i,SQRT(rmserq)
    END IF
    
  END IF
  
END DO


rmsavq = SQRT(rmsavq/natom)
rmsavm = SQRT(rmsavm/natom)

WRITE(1337,'(79(1H-),/)')
IF (nspin == 2) THEN
  WRITE (ipf,FMT=9040) itc,rmsavq,rmsavm
  WRITE ( 6 ,FMT=9040) itc,rmsavq,rmsavm
ELSE
  WRITE (ipf,FMT=9050) itc,rmsavq
  WRITE ( 6 ,FMT=9050) itc,rmsavq
END IF
WRITE(1337,'(79(1H-))')

9000 FORMAT (5X,' rms-error for atom',i3,1X,':','v+ + v- = ',1P,d11.4,  &
    2X,',',2X,'v+ - v- = ',1P,d11.4)
9010 FORMAT (5X,' rms-error non spherical contribution for atom ',i3,  &
    1X,':','v+ + v- = ',1P,d11.4,02X,',',2X,'v+ - v- = ',1P, d11.4)
9020 FORMAT (5X,' rms-error for atom',i3,1X,':','v+ + v- = ',1P,d11.4)
9030 FORMAT (5X,' rms-error non spherical contribution for atom ',i3,  &
    1X,':','v+ + v- = ',1P,d11.4)
9040 FORMAT ('      ITERATION',i4,' average rms-error : v+ + v- = ',  &
    1P,d11.4,/,39X,' v+ - v- = ',1P,d11.4)
9050 FORMAT ('      ITERATION',i4,' average rms-error : v+ + v- = ', 1P,d11.4)
END SUBROUTINE mixstr
