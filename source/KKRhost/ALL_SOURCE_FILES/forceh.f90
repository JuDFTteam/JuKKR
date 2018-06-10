SUBROUTINE forceh(cmom,flmh,lmax,nspin,nstart,nend,r2rho,v,r,drdi,  &
        irws,z)
!-----------------------------------------------------------------------
!     calculates the force on nucleus m with hellmann - feynman theorem
!     from a given non spherical charge density at the nucleus site r


!-----------------------------------------------------------------------
IMPLICIT NONE
!.. Parameters ..
include 'inc.p'
INTEGER LMPOTD
PARAMETER (LMPOTD= (LPOTD+1)**2)
!..
!.. Scalar Arguments ..
INTEGER LMAX,NEND,NSPIN,NSTART
!..
!.. Array Arguments ..
DOUBLE PRECISION CMOM(LMPOTD,*),DRDI(IRMD,*),FLMH(-1:1,*), &
       R(IRMD,*), &
       R2RHO(IRMD,LMPOTD,NATYPD,*),V(IRMD,LMPOTD,*),Z(*)
INTEGER IRWS(*)
!..
!.. Local Scalars ..
DOUBLE PRECISION PI,RWS,VINT1
INTEGER I,IATYP,IPOT,IRWS1,LM,M
!..
!.. Local Arrays ..
DOUBLE PRECISION FLM(-1:1,2),V1(IRMD)
!..
!.. External Subroutines ..
EXTERNAL SIMP3
!..
!.. Save statement ..
SAVE PI
!..

!.. Intrinsic Functions ..
INTRINSIC ATAN
!..
pi = 4.d0*ATAN(1.d0)
IF (lmax < 1) THEN
  WRITE (6,FMT=9000)
  STOP
  
END IF

!---> loop over the rep. atoms

DO  iatyp = nstart,nend
  
!---> reading the right Wigner-S. radius
  
  irws1 = irws(iatyp)
  rws = r(irws1,iatyp)
  
!---> determine the right potential numbers
  
  ipot = nspin* (iatyp-1) + 1
  
  DO  m = -1,1
    lm = 2 + m + 1
    
    v1(1) = 0.0D0
    DO  i = 2,irws1
      v1(i) = r2rho(i,lm,iatyp,1)* (r(i,iatyp)** (-2.0D0))
    END DO
    
!---> integrate with simpson subroutine
    
    CALL simp3(v1,vint1,1,irws1,drdi(1,iatyp))
    
    flm(m,1) = 2.0D0*vint1
    
!---> use coulomb potential to determine extra atomic contribution
    
    flm(m,2) = v(irws1,lm,ipot)* (3.0D0/ (4.0D0*pi*rws)) -  &
        2.0D0*cmom(lm,iatyp)/ (rws**3)
    
!---> total Hellman-Feynman force
    
    flmh(m,iatyp) = (flm(m,1)+flm(m,2))*z(iatyp)
  END DO
END DO


9000 FORMAT (13X,'error stop in subroutine force :',  &
    ' the charge density has to contain non spherical',  &
    ' contributions up to l=1 at least ')

END SUBROUTINE forceh
