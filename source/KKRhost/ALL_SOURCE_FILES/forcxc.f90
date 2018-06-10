SUBROUTINE forcxc(flm,flmc,lmax,nspin,nstart,nend,rhoc,v,r,alat,  &
        drdi,irws,natref)
!>>>>>BEWARE!!! RM commented away!!! -->Dipole Tensor is useless
!     SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,ALAT,
!    +                  RM,NSHELL,DRDI,IRWS,NATREF)
!-----------------------------------------------------------------------
!     calculates the force on nucleus m
!     from a given non spherical charge density at the nucleus site r
!     with core correction(exchange contribution)

!-----------------------------------------------------------------------
use mod_types, only: t_inc
IMPLICIT NONE
!     .. Parameters ..
INCLUDE 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER LMAX,NATREF,NEND,NSPIN,NSTART
!..
!.. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),FLM(-1:1,*),FLMC(-1:1,*),R(IRMD,*), &
            RHOC(IRMD,*), &
!     $     RM(3,*),
          V(IRMD,LMPOTD,*)
!      INTEGER IRWS(*),NSHELL(*)
       INTEGER IRWS(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION DV,FAC,PI,RWS,TRP,VINT1,VOL!,DVOL
      INTEGER I,IATYP,IPER,IPOT,IREP,IRWS1,ISPIN,LM,M!,J
!..
!.. Local Arrays ..
      DOUBLE PRECISION F(3,NATYPD),FLMH(-1:1,NATYPD), &
             FLMXC(-1:1,NATYPD),P(NATYPD),V1(IRMD)
!..
!.. External Subroutines ..
      EXTERNAL SIMP3
!..
!.. Save statement ..
      SAVE PI
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,DSQRT
!     ..
pi = 4.d0*ATAN(1.d0)
fac = DSQRT((4.0D0*pi)/3.0D0)
trp = 0.0D0
IF (lmax < 1) THEN
  WRITE (6,FMT=9000)
  STOP
  
END IF

IF(t_inc%i_write>0) WRITE (1337,FMT=9200)
IF(t_inc%i_write>0) WRITE (1337,FMT=9100)
IF(t_inc%i_write>0) WRITE (1337,FMT=9200)

irep = 1
DO  iatyp = nstart,nend
  
  iper = iatyp - natref
  p(iper) = 0.0D0
  IF(t_inc%i_write>0) WRITE (1337,FMT=9400) iper
  
  irws1 = irws(iatyp)
  rws = r(irws1,iatyp)
  vol = 0.25*alat**3
  
!---> determine the right potential numbers
  
  
  DO  m = -1,1
    lm = 2 + m + 1
    
    DO  i = 1,irws1
      v1(i) = 0.0D0
    END DO
    
    DO  ispin = 1,nspin
      
      ipot = nspin* (iatyp-1) + ispin
      
      
      dv = (-3.0D0*v(1,lm,ipot)-10.0D0*v(2,lm,ipot)+  &
          18.0D0*v(3,lm,ipot)-6.0D0*v(4,lm,ipot)+v(5,lm,ipot))/  &
          (12.0D0*drdi(2,iatyp))
      
      v1(2) = rhoc(2,ipot)* (2.0D0*v(2,lm,ipot)/r(2,iatyp)+dv)/  &
          (4.0D0*pi) + v1(2)
      
      DO  i = 3,irws1 - 2
        
        dv = (v(i-2,lm,ipot)-v(i+2,lm,ipot)+  &
            8.0D0* (v(i+1,lm,ipot)-v(i-1,lm,ipot)))/ (12.0D0*drdi(i,iatyp))
        
        v1(i) = rhoc(i,ipot)* (2.0D0*v(i,lm,ipot)/r(i,iatyp)+  &
            dv)/ (4.0D0*pi) + v1(i)
      END DO
      
      dv = (-v(irws1-4,lm,ipot)+6.0D0*v(irws1-3,lm,ipot)-  &
          18.0D0*v(irws1-2,lm,ipot)+10.0D0*v(irws1-1,lm,ipot)+  &
          3.0D0*v(irws1,lm,ipot))/ (12.0D0*drdi(irws1-1,iatyp))
      v1(irws1-1) = rhoc(irws1-1,ipot)*  &
          (2.0D0*v(irws1-1,lm,ipot)/r(irws1-1,iatyp)+  &
          dv)/ (4.0D0*pi) + v1(irws1-1)
      
      dv = (3.0D0*v(irws1-4,lm,ipot)-16.0D0*v(irws1-3,lm,ipot)+  &
          36.0D0*v(irws1-2,lm,ipot)-48.0D0*v(irws1-1,lm,ipot)+  &
          25.0D0*v(irws1,lm,ipot))/ (12.0D0*drdi(irws1,iatyp))
      
      v1(irws1) = rhoc(irws1,ipot)*  &
          (2.0D0*v(irws1,lm,ipot)/r(irws1,iatyp)+dv)/ (4.0D0*pi) + v1(irws1)
    END DO
    
!---> integrate with simpson subroutine
    
    CALL simp3(v1,vint1,1,irws1,drdi(1,iatyp))
    
    flmh(m,iatyp) = flm(m,iatyp) - flmc(m,iatyp)
    flmxc(m,iatyp) = -fac*vint1 - flmc(m,iatyp)
    flm(m,iatyp) = flm(m,iatyp) + flmxc(m,iatyp)
    
    
  END DO
  
  IF(t_inc%i_write>0) THEN
    WRITE (1337,FMT=9600) flmh(1,iatyp),flmc(1,iatyp),  &
        flmxc(1,iatyp),flm(1,iatyp)
    WRITE (1337,FMT=9601) flmh(-1,iatyp),flmc(-1,iatyp),  &
        flmxc(-1,iatyp),flm(-1,iatyp)
    WRITE (1337,FMT=9602) flmh(0,iatyp),flmc(0,iatyp),  &
        flmxc(0,iatyp),flm(0,iatyp)
  END IF
  
  f(1,iatyp) = flm(1,iatyp)
  f(2,iatyp) = flm(-1,iatyp)
  f(3,iatyp) = flm(0,iatyp)
  
  
!         DO 60 J = 1,3
!            P(IPER) = P(IPER) + RM(J,IREP)*NSHELL(IPER)*F(J,IATYP)*ALAT
!   60    CONTINUE
!         TRP = TRP + P(IPER)
  
!         IREP = IREP + NSHELL(IPER)
  
!        write (6,*) '-->Tensor is useless'
!        WRITE (6,FMT=9700) P(IPER)
  
END DO

!     DVOL = TRP/ (3.0D0*VOL)

IF(t_inc%i_write>0) THEN
  WRITE (1337,FMT=9200)
!       WRITE (6,FMT=9101)
!       WRITE (6,FMT=9200)
!       WRITE (6,FMT=9800) DVOL
!       WRITE (6,FMT=9200)
  WRITE (1337,FMT=9102)
  WRITE (1337,FMT=9200)
END IF

9000 FORMAT (13X,'error stop in subroutine force :',  &
    ' the charge density has to contain non spherical',  &
    ' contributions up to l=1 at least ')
!  9101 FORMAT (1x,33 ('-'),' volume change ',33 ('-'),/,34x,
!      +       ' in units Ry/(a(Bohr)**3 ')
9102 FORMAT (1X,81 ('-'))
9100 FORMAT (1X,33 ('-'),' force on the nucleus ',33 ('-'),/,34X,  &
    ' in units Ry/(a(Bohr) ')
9200 FORMAT (1X,'>')
9400 FORMAT (3X,i5,'th shell')
9600 FORMAT (7X,'fhx=',e13.6,2X,'fcx=',e13.6,2X,'fxcx=',e13.6,2X,'fx=',  &
    e13.6,' Ry/(a(Bohr))')
9601 FORMAT (7X,'fhy=',e13.6,2X,'fcy=',e13.6,2X,'fxcy=',e13.6,2X,'fy=',  &
    e13.6,' Ry/(a(Bohr))')
9602 FORMAT (7X,'fhz=',e13.6,2X,'fcz=',e13.6,2X,'fxcz=',e13.6,2X,'fz=',  &
    e13.6,' Ry/(a(Bohr))')
!  9700 FORMAT (10x,'contribution to the trace of the dipol force tensor:'
!      +       ,3x,e12.6,' Ry')
!  9800 FORMAT (7x,' volume change dvol/vol=',2x,e12.6,' Ry/(a(Bohr))**3',
!      +       /,7x,'( notice: has to be divided',
!      +       ' by the bulk modulus of the host)')

END SUBROUTINE forcxc
