SUBROUTINE vxcgga(exc,kte,lpot,nspin,rho2ns,v,r,drdi,a, &
        irws,ircut,ipan,gsh,ilm,imaxsh,  &
        ifunm,thetas,wtyr,ijend,lmsp,thet,ylm,dylmt1, &
        dylmt2,dylmf1,dylmf2,dylmtf, &
        irmd, irid, nfund, ngshd, ipand, kxc)
!-----------------------------------------------------------------------
!     add the exchange-correlation-potential to the given potential
!     and if total energies should be calculated (kte=1) the exchange-
!     correlation-energies are calculated .
!     use as input the charge density times r**2 (rho2ns(...,1)) and
!     in the spin-polarized case (nspin=2) the spin density times r**2
!     (rho2ns(...,2)) .
!     the density times 4 pi is generated at an angular mesh .
!     the exchange-correlation potential and the exchange-correlation
!     energy are calculated at those mesh points with a subroutine .
!     in the non-spin-polarized case the "spin-density" is
!     set equal zero .
!     after that the exchange-correlation potential and in the case of
!     total energies (kte=1) the exchange-correlation energy are
!     expanded into spherical harmonics .
!     the ex.-cor. potential is added to the given potential .
!     the expansion into spherical harmonics uses the orthogonality
!     of these harmonics . - therefore a gauss-legendre integration
!     for "theta" and a gauss-tschebyscheff integration for "phi"
!     is used .
!     all needed values for the angular mesh and angular integration
!     are generate in the subroutine sphere .

!     the ex.-cor. potential is extrapolated to the origin only
!     for the lm=1 value .

!                               b.drittler   june 1987

!     modified for shape functions
!                                       b. drittler oct. 1989
!     simplified and modified for Paragon X/PS
!                                       R. Zeller Nov. 1993
!-----------------------------------------------------------------------

use XCFunctionals_mod

IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT)            :: exc(0:(lpot+1)**2)
INTEGER, INTENT(IN)                      :: kte
INTEGER, INTENT(IN)                      :: lpot
INTEGER, INTENT(IN)                      :: nspin
DOUBLE PRECISION, INTENT(IN)             :: rho2ns(irmd,(lpot+1)**2,2)
DOUBLE PRECISION, INTENT(OUT)            :: v(irmd,(lpot+1)**2,2)
DOUBLE PRECISION, INTENT(IN)             :: r(irmd)
DOUBLE PRECISION, INTENT(IN OUT)         :: drdi(irmd)
DOUBLE PRECISION, INTENT(IN)             :: a
INTEGER, INTENT(IN)                      :: irws
INTEGER, INTENT(IN)                      :: ircut(0:ipand)
INTEGER, INTENT(IN)                      :: ipan
DOUBLE PRECISION, INTENT(IN)             :: gsh(*)
INTEGER, INTENT(IN)                      :: ilm(ngshd,3)
INTEGER, INTENT(IN)                      :: imaxsh(0:(lpot+1)**2)
INTEGER, INTENT(IN)                      :: ifunm((2*lpot+1)**2)
DOUBLE PRECISION, INTENT(IN OUT)         :: thetas(irid,nfund)
DOUBLE PRECISION, INTENT(IN)             :: wtyr(ijend,*)
INTEGER, INTENT(IN OUT)                  :: ijend
INTEGER, INTENT(IN)                      :: lmsp((2*lpot+1)**2)
DOUBLE PRECISION, INTENT(IN OUT)         :: thet(ijend)
DOUBLE PRECISION, INTENT(IN OUT)         :: ylm(ijend,(lpot+1)**2)
DOUBLE PRECISION, INTENT(IN OUT)         :: dylmt1(ijend,(lpot+1)**2)
DOUBLE PRECISION, INTENT(IN OUT)         :: dylmt2(ijend,(lpot+1)**2)
DOUBLE PRECISION, INTENT(IN OUT)         :: dylmf1(ijend,(lpot+1)**2)
DOUBLE PRECISION, INTENT(IN OUT)         :: dylmf2(ijend,(lpot+1)**2)
DOUBLE PRECISION, INTENT(IN OUT)         :: dylmtf(ijend,(lpot+1)**2)
INTEGER, INTENT(IN)                      :: kxc


!     ..
!     .. Local Scalars ..
DOUBLE PRECISION :: chgden,dx,elmxc,fpi,r1,r2,rpoint,spiden,vlmxc,  &
    vxc1,vxc2,vxc3,zero,zero1,ddot
INTEGER :: ifun,ipan1,ipot,ir,irc0,irc1,irh,irs1,ispin,j,l,l1max,lm,  &
    lm2,lmmax,m,mesh,nspin2,irmd,irid,nfund,ngshd,ipand
!     ..
!     .. Local Arrays ..

DOUBLE PRECISION :: ddrrl(irmd,(lpot+1)**2)
DOUBLE PRECISION :: ddrrul(irmd,(lpot+1)**2)
DOUBLE PRECISION :: drrl(irmd,(lpot+1)**2)
DOUBLE PRECISION :: drrul(irmd,(lpot+1)**2)
DOUBLE PRECISION :: er(irmd,0:(lpot+1)**2)
DOUBLE PRECISION :: estor(irmd,(lpot+1)**2)
DOUBLE PRECISION :: excij(ijend)
DOUBLE PRECISION :: rhol(irmd,2,(lpot+1)**2)
DOUBLE PRECISION :: rholm((lpot+1)**2,2)
DOUBLE PRECISION :: vxc(ijend,2)
DOUBLE PRECISION :: vxcr(2:3,2)
!     ..

!     .. External Functions ..
EXTERNAL ddot
!     ..
!     .. External Subroutines ..
EXTERNAL gradrl,mkxcpe,simpk
!     ..


!     .. Intrinsic Functions ..
INTRINSIC ABS,ATAN,MOD
!     ..

!     .. Data statements ..
DATA zero,zero1/0.d0,1.d-12/

INTEGER :: lmpotd

lmpotd= (lpot+1)**2

!     ..
WRITE (6,FMT=*) ' GGA CALCULATION '
fpi = 16.0D0*ATAN(1.0D0)
lmmax = (lpot+1)* (lpot+1)

!---> loop over given representive atoms

ipan1 = ipan
irc1 = ircut(ipan)
irs1 = ircut(1)
irc0 = 2

DO  ispin = 1,nspin
  vxcr(2,ispin) = 0.0D0
  vxcr(3,ispin) = 0.0D0
END DO

!---> initialize for ex.-cor. energy

IF (kte == 1) THEN
  DO  l = 0,lpot
    exc(l) = 0.0D0
    DO  ir = 1,irc1
      er(ir,l) = 0.0D0
    END DO
  END DO
  
  DO  lm = 1,lmmax
    DO  ir = 1,irc1
      estor(ir,lm) = 0.0D0
    END DO
  END DO
END IF

l1max = lpot + 1
mesh = irws
dx = a

IF (nspin == 2) THEN
  DO  lm = 1,lmmax
    DO  ir = 2,mesh
      r1 = r(ir)
      r2 = r1*r1
      chgden = rho2ns(ir,lm,1)/r2
      spiden = rho2ns(ir,lm,2)/r2
      IF (ABS(chgden) <= zero1) chgden = zero
      IF (ABS(spiden) <= zero1) spiden = zero
      rhol(ir,2,lm) = (chgden+spiden)/2.d0
      rhol(ir,1,lm) = (chgden-spiden)/2.d0
    END DO
    
!       extrapolate
    
    rhol(1,1,lm) = rhol(2,1,lm)
    rhol(1,2,lm) = rhol(2,2,lm)
  END DO
  
ELSE
  
  DO  lm = 1,lmmax
    DO  ir = 2,mesh
      r1 = r(ir)
      r2 = r1*r1
      
      chgden = rho2ns(ir,lm,1)/r2
      IF (ABS(chgden) <= zero1) chgden = zero
      rhol(ir,1,lm) = chgden/2.d0
      rhol(ir,2,lm) = chgden/2.d0
    END DO
    
!       extrapolate
    rhol(1,1,lm) = rhol(2,1,lm)
    rhol(1,2,lm) = rhol(2,2,lm)
  END DO
END IF


CALL gradrl(nspin,mesh,l1max,dx,rhol,r,drdi,ipan1,ipand,ircut,  &
    drrl,ddrrl,drrul,ddrrul,irmd,lmpotd)


!---> loop over radial mesh


DO  ir = irc0,irc1
  rpoint = r(ir)
  
!---> calculate the ex.-cor. potential
  
  nspin2 = 2
  
  DO  ispin = 1,nspin2
    DO  lm = 1,lmmax
      rholm(lm,ispin) = rhol(ir,ispin,lm)
    END DO
  END DO
  
!    only for spin-polarized
  
!        write(6,*) ' before mkxcpe '
IF (kxc==3) THEN
!  WRITE(*,*) 'kxc=3: Using GGA functional PW91'  
  CALL mkxcpe_pw91(nspin2,ir,ijend,l1max,rpoint,rholm,vxc,excij,thet,  &
      ylm,dylmt1,dylmt2,dylmf1,dylmf2,dylmtf,drrl,ddrrl, drrul,ddrrul,irmd,lmpotd)
ELSE IF (kxc==4) THEN
!  WRITE(*,*) 'kxc=4: Using GGA functional PBE'
  CALL mkxcpe_pbe(IR,IJEND,RPOINT,RHOLM,VXC,EXCIJ,YLM,DYLMT1,  &
               DYLMF1,DYLMF2,DYLMTF,DRRL,DDRRL,DRRUL,DDRRUL,  &
               IRMD,LMPOTD,LMMAX,.false.)
ELSE IF (kxc==5) THEN
!  WRITE(*,*) 'kxc=5: Using GGA functional PBEsol'
  CALL mkxcpe_pbe(IR,IJEND,RPOINT,RHOLM,VXC,EXCIJ,YLM,DYLMT1,  &
               DYLMF1,DYLMF2,DYLMTF,DRRL,DDRRL,DRRUL,DDRRUL,  &
               IRMD,LMPOTD,LMMAX,.true.)
ELSE
WRITE(*,*) 'Error in vxcgga: In order to use GGA set kxc=3 (PW91), kxc=4 (PBE) &
or kxc=5 (PBEsol). LDA is used for kxc<3)'
STOP
END IF 
  
  
  
  
!---> expand the ex.-cor. potential into spherical harmonics ,
!       using the orthogonality
  
  DO  ispin = 1,nspin
    
!---> determine the corresponding potential number
    
    ipot = ispin
    DO  lm = 1,lmmax
      vlmxc = ddot(ijend,vxc(1,ispin),1,wtyr(1,lm),1)
      v(ir,lm,ipot) = v(ir,lm,ipot) + vlmxc
      
!---> store the ex.-c. potential of ir=2 and =3 for the extrapolation
      
      IF (lm == 1 .AND. (ir == 2.OR.ir == 3)) vxcr(ir, ispin) = vlmxc
    END DO
  END DO
  
!---> file er in case of total energies
  
  IF (kte == 1) THEN
    
!---> expand ex.-cor. energy into spherical harmonics
!       using the orthogonality
    
    DO  l = 0,lpot
      DO  m = -l,l
        lm = l*l + l + m + 1
        elmxc = ddot(ijend,excij,1,wtyr(1,lm),1)
        
!---> multiply the lm-component of the ex.-cor. energy with the same
!     lm-component of the charge density times r**2 and sum over lm
!     this corresponds to a integration over the angular .
        
        IF (ir > irs1) THEN
          estor(ir,lm) = elmxc
          
        ELSE
          
          er(ir,l) = er(ir,l) + rho2ns(ir,lm,1)*elmxc
        END IF
        
      END DO
      
    END DO
    
  END IF
  
END DO

!---> integrate er in case of total energies to get exc

IF (kte == 1) THEN
  
  DO  l = 0,lpot
    DO  m = -l,l
      lm = l*l + l + m + 1
      
!---> convolute with shape function
      
      DO  j = imaxsh(lm-1) + 1,imaxsh(lm)
        lm2 = ilm(j,2)
        IF (lmsp(ilm(j,3)) > 0) THEN
          ifun = ifunm(ilm(j,3))
          DO  ir = irs1 + 1,irc1
            irh = ir - irs1
            er(ir,l) = er(ir,l) + rho2ns(ir,lm,1)*gsh(j)*  &
                thetas(irh,ifun)*estor(ir,lm2)
          END DO
        END IF
      END DO
    END DO
    CALL simpk(er(1,l),exc(l),ipan1,ircut,drdi)
  END DO
  
END IF

!---> extrapolate ex.-cor potential to the origin only for lm=1

DO  ispin = 1,nspin
  ipot = ispin
  
  vxc2 = vxcr(2,ispin)
  vxc3 = vxcr(3,ispin)
  vxc1 = vxc2 - r(2)* (vxc3-vxc2)/ (r(3)-r(2))
  
  v(1,1,ipot) = v(1,1,ipot) + vxc1
END DO

END SUBROUTINE vxcgga
