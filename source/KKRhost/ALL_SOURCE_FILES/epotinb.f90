! 13.10.95 ***************************************************************
SUBROUTINE epotinb(epotin,nspin,natyp,rho2ns,vm2z,r,drdi,ins,  &
    irmin,irws,lpot,vins,ircut,ipan,z)
! ************************************************************************

!     attention : energy zero ---> electro static zero

!                 since input potential and single particle energies
!                 are using muffin tin zero as zero the energy shift
!                 is cancelled in the kinetic energy contribution !


!     calculate the energy of the input potential
!     the energy for the representive atom i is given by

!                               rws
!       epotin(i) = - sqrt(4 pi) {  dr' vm2z(r',i)*rho2ns(r',1,i)
!                                0

!     in case of non spherical input potential one has to add

!                 rirt
!            {  -  {  dr' vins(r',lm,i)rho2ns(r',lm,i)   }
!                 rmin
!                                        (summed over lm)

!     remember : the non spherical part of the input potential is
!                different from zero only between r(irmin) and r(irt)

!             (see notes by b.drittler)

!     attention: vm2z is the spherically averaged input potential ,
!                vins contains the non spherical contribution of the
!                potential and rho2ns(...,1) is the  real charge density
!                times r**2. vins and rho2ns are expanded into spherical
!                harmonics . (see deck rholm or rhons)

!     remember :  in case of shape corrections  the contribution of
!                 the nuclear potential - 2*Z/r has to be explicitly
!                 taken into account between muffin tin sphere and
!                 circum scribed sphere .
!                 only within the muffin tin sphere this term is
!                 analytically cancelled wtih the contribution of
!                 the coulomb potential - see deck ecoulom


!                 modified for non spherical potential and shape correc-
!                  tions

!                               b.drittler   oct. 1989
!-----------------------------------------------------------------------
!.. Parameters ..
      include 'inc.p'
!..
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
!..
!.. Scalar Arguments ..
      INTEGER INS,LPOT,NATYP,NSPIN
!..
!.. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),EPOTIN(*),R(IRMD,*), &
                       RHO2NS(IRMD,LMPOTD,NATYPD,*), &
                       VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IRMIN(*),IRWS(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION PI,R2RHOD,R2RHOU,RFPI,TEMP,ZZOR
      INTEGER I,IATYP,IC,IPAN1,IPOTD,IPOTU,IRC1,IRMIN1,IRS1,L1,LM,M1
!..
!.. Local Arrays ..
      DOUBLE PRECISION ENS(0:LPOTD,NATYPD),ER(IRMD)
      INTEGER IRCUTM(0:IPAND)
!..
!.. External Subroutines ..
      EXTERNAL SIMP3,SIMPK
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
!     ..
pi = 4.0D0*ATAN(1.0D0)
rfpi = SQRT(4.0D0*pi)

DO  iatyp = 1,natyp
  
  ipan1 = ipan(iatyp)
  irc1 = ircut(ipan1,iatyp)
  
  IF (ipan1 > 1) THEN
    irs1 = ircut(1,iatyp)
  ELSE
    irs1 = irws(iatyp)
  END IF
  
  IF (nspin == 1) THEN
    ipotu = iatyp
    ipotd = iatyp
  ELSE
    ipotu = 2*iatyp - 1
    ipotd = 2*iatyp
  END IF
  
  DO  i = 1,irs1
    
!---> calculate charge density times input potential
    
    r2rhou = (rho2ns(i,1,iatyp,1)-rho2ns(i,1,iatyp,nspin))/2.0D0
    r2rhod = (rho2ns(i,1,iatyp,1)+rho2ns(i,1,iatyp,nspin))/2.0D0
    er(i) = - (r2rhou*vm2z(i,ipotu)+r2rhod*vm2z(i,ipotd))*rfpi
  END DO
  
!--->  remember the form of vm2z between mt sphere and rirc
  
  IF (ipan1 > 1) THEN
    DO  i = irs1 + 1,irc1
      r2rhou = (rho2ns(i,1,iatyp,1)-rho2ns(i,1,iatyp,nspin))/2.0D0
      r2rhod = (rho2ns(i,1,iatyp,1)+rho2ns(i,1,iatyp,nspin))/2.0D0
      zzor = 2.0D0*z(iatyp)/r(i,iatyp)
      er(i) = - (r2rhou* (vm2z(i,ipotu)-zzor)+  &
          r2rhod* (vm2z(i,ipotd)-zzor))*rfpi
    END DO
  END IF
  
!--->   now integrate er to get epotin
  
  IF (ipan1 > 1) THEN
    CALL simpk(er,temp,ipan(iatyp),ircut(0,iatyp),drdi(1,iatyp))
  ELSE
    CALL simp3(er,temp,1,irs1,drdi(1,iatyp))
  END IF
  
  epotin(iatyp) = temp
  ens(0,iatyp) = temp
  
!--->   add non spher. contribution in case of non spher. input potential
  
  DO  l1 = 1,lpot
    ens(l1,iatyp) = 0.0D0
  END DO
  
  IF (ins /= 0) THEN
    
    irmin1 = irmin(iatyp)
    IF (irmin1 <= irs1) THEN
      
      ircutm(0) = irmin1 - 1
      DO  ic = 1,ipan1
        ircutm(ic) = ircut(ic,iatyp)
      END DO
      
      DO  l1 = 1,lpot
        
        DO  i = 1,irmd
          er(i) = 0.0D0
        END DO
        
        DO  m1 = -l1,l1
          lm = l1* (l1+1) + m1 + 1
          DO  i = irmin1,irc1
            
!---> calculate charge density times potential
            
            r2rhou = (rho2ns(i,lm,iatyp,1)- rho2ns(i,lm,iatyp,nspin))/2.0D0
            r2rhod = (rho2ns(i,lm,iatyp,1)+ rho2ns(i,lm,iatyp,nspin))/2.0D0
            er(i) = er(i) - r2rhou*vins(i,lm,ipotu) - r2rhod*vins(i,lm,ipotd)
          END DO
        END DO
        CALL simpk(er,temp,ipan1,ircutm,drdi(1,iatyp))
        
        epotin(iatyp) = epotin(iatyp) + temp
        ens(l1,iatyp) = temp
      END DO
      
    END IF                    ! (IRMIN1.LE.IRS1)
    
  END IF                      ! (INS.NE.0)
  
END DO

RETURN
END SUBROUTINE epotinb
