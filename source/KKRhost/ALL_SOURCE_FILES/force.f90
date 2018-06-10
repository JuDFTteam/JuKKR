SUBROUTINE force(flm,flmc,lmax,nspin,nstart,nend,rhoc,v,r,drdi,  &
        irws)
!-----------------------------------------------------------------------
!     calculates the force on nucleus m
!     from a given non spherical charge density at the nucleus site r
!     with core correction (coulomb contribution)

!-----------------------------------------------------------------------
IMPLICIT NONE

!     .. Parameters ..
INCLUDE 'inc.p'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
!..
!.. Scalar Arguments ..
      INTEGER LMAX,NEND,NSPIN,NSTART
!..
!.. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),FLM(-1:1,*),FLMC(-1:1,*),R(IRMD,*), &
             RHOC(IRMD,*),V(IRMD,LMPOTD,*)
      INTEGER IRWS(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION DV,FAC,PI,RWS,VINT1
      INTEGER I,IATYP,IPOT,IRWS1,ISPIN,LM,M
!..
!.. Local Arrays ..
      DOUBLE PRECISION FLMH(-1:1,NATYPD),V1(IRMD)
!..
!.. External Subroutines ..
      EXTERNAL SIMP3
!..
!.. Save statement ..
      SAVE PI
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,DSQRT
!..
pi = 4.d0*ATAN(1.d0)
fac = DSQRT((4.0D0*pi)/3.0D0)
IF (lmax < 1) THEN
  WRITE (6,FMT=9000)
  STOP
  
END IF

!---> loop over rep. atoms

DO  iatyp = nstart,nend
  
  
  irws1 = irws(iatyp)
  rws = r(irws1,iatyp)
  
  
  
  DO  m = -1,1
    lm = 2 + m + 1
    
!---> initialize v1
    
    DO  i = 1,irws1
      v1(i) = 0.0D0
    END DO
    
    DO  ispin = 1,nspin
      
!---> determine the right potential numbers
      
      ipot = nspin* (iatyp-1) + ispin
      
!---> determine the derivative of the potential using a 5-point formular
      
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
    
    flmh(m,iatyp) = fac*flm(m,iatyp)
    flmc(m,iatyp) = -fac*vint1
    flm(m,iatyp) = flmh(m,iatyp) + flmc(m,iatyp)
    
    
  END DO
  
  
END DO

9000 FORMAT (13X,'error stop in subroutine force :',  &
    ' the charge density has to contain non spherical',  &
    ' contributions up to l=1 at least ')

END SUBROUTINE force
