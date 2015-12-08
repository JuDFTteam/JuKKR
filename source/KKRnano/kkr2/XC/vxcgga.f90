subroutine vxcgga(exc,kte,lpot,nspin,rho2ns,v,r,drdi,a, &
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
!     simplified and modified for paragon x/ps
!                                       r. zeller nov. 1993
!-----------------------------------------------------------------------

  use xcfunctionals_mod ! todo: fill only-list
  use quadrature_mod, only: simpson

implicit none

double precision, intent(out)            :: exc(0:(lpot+1)**2)
integer, intent(in)                      :: kte
integer, intent(in)                      :: lpot
integer, intent(in)                      :: nspin
double precision, intent(in)             :: rho2ns(irmd,(lpot+1)**2,2)
double precision, intent(out)            :: v(irmd,(lpot+1)**2,2)
double precision, intent(in)             :: r(irmd)
double precision, intent(in)             :: drdi(irmd)
double precision, intent(in)             :: a
integer, intent(in)                      :: irws
integer, intent(in)                      :: ircut(0:ipand)
integer, intent(in)                      :: ipan
double precision, intent(in)             :: gsh(*)
integer, intent(in)                      :: ilm(ngshd,3)
integer, intent(in)                      :: imaxsh(0:(lpot+1)**2)
integer, intent(in)                      :: ifunm((2*lpot+1)**2)
double precision, intent(in)             :: thetas(irid,nfund)
double precision, intent(in)             :: wtyr(ijend,*)
integer, intent(in)                      :: ijend
integer, intent(in)                      :: lmsp((2*lpot+1)**2)
double precision, intent(in out)         :: thet(ijend)
double precision, intent(in out)         :: ylm(ijend,(lpot+1)**2)
double precision, intent(in out)         :: dylmt1(ijend,(lpot+1)**2)
double precision, intent(in out)         :: dylmt2(ijend,(lpot+1)**2)
double precision, intent(in out)         :: dylmf1(ijend,(lpot+1)**2)
double precision, intent(in out)         :: dylmf2(ijend,(lpot+1)**2)
double precision, intent(in out)         :: dylmtf(ijend,(lpot+1)**2)
integer, intent(in)                      :: kxc

!     ..
!     .. local scalars ..
double precision :: chgden,dx,elmxc,fpi,r1,r2,rpoint,spiden,vlmxc,  &
    vxc1,vxc2,vxc3
integer :: ifun,ipan1,ipot,ir,irc0,irc1,irh,irs1,ispin,j,l,l1max,lm,  &
    lm2,lmmax,m,mesh,nspin2,irmd,irid,nfund,ngshd,ipand
!     ..
!     .. local arrays ..

double precision :: ddrrl(irmd,(lpot+1)**2)
double precision :: ddrrul(irmd,(lpot+1)**2)
double precision :: drrl(irmd,(lpot+1)**2)
double precision :: drrul(irmd,(lpot+1)**2)
double precision :: er(irmd,0:(lpot+1)**2)
double precision :: estor(irmd,(lpot+1)**2)
double precision :: excij(ijend)
double precision :: rhol(irmd,2,(lpot+1)**2)
double precision :: rholm((lpot+1)**2,2)
double precision :: vxc(ijend,2)
double precision :: vxcr(2:3,2)
!     ..

double precision, external :: ddot
external :: gradrl, mkxcpe

double precision, parameter :: zero = 0.d0, zero1 = 1.d-12

integer :: lmpotd

lmpotd= (lpot+1)**2

!     ..
write (6,fmt=*) ' GGA CALCULATION '
fpi = 16.d0*atan(1.d0)
lmmax = (lpot+1)* (lpot+1)

!---> loop over given representive atoms

ipan1 = ipan
irc1 = ircut(ipan)
irs1 = ircut(1)
irc0 = 2

do ispin = 1,nspin
  vxcr(2,ispin) = 0.d0
  vxcr(3,ispin) = 0.d0
enddo

!---> initialize for ex.-cor. energy

if (kte == 1) then
  do l = 0,lpot
    exc(l) = 0.d0
    do ir = 1,irc1
      er(ir,l) = 0.d0
    enddo
  enddo
  
  do lm = 1,lmmax
    do ir = 1,irc1
      estor(ir,lm) = 0.d0
    enddo
  enddo
endif

l1max = lpot + 1
mesh = irws
dx = a

if (nspin == 2) then
  do lm = 1,lmmax
    do ir = 2,mesh
      r1 = r(ir)
      r2 = r1*r1
      chgden = rho2ns(ir,lm,1)/r2
      spiden = rho2ns(ir,lm,2)/r2
      if (abs(chgden) <= zero1) chgden = zero
      if (abs(spiden) <= zero1) spiden = zero
      rhol(ir,2,lm) = (chgden+spiden)/2.d0
      rhol(ir,1,lm) = (chgden-spiden)/2.d0
    enddo
    
!       extrapolate
    
    rhol(1,1,lm) = rhol(2,1,lm)
    rhol(1,2,lm) = rhol(2,2,lm)
  enddo
  
else
  
  do lm = 1,lmmax
    do ir = 2,mesh
      r1 = r(ir)
      r2 = r1*r1
      
      chgden = rho2ns(ir,lm,1)/r2
      if (abs(chgden) <= zero1) chgden = zero
      rhol(ir,1,lm) = chgden/2.d0
      rhol(ir,2,lm) = chgden/2.d0
    enddo
    
!       extrapolate
    rhol(1,1,lm) = rhol(2,1,lm)
    rhol(1,2,lm) = rhol(2,2,lm)
  enddo
endif


call gradrl(nspin,mesh,l1max,dx,rhol,r,drdi,ipan1,ipand,ircut,  &
    drrl,ddrrl,drrul,ddrrul,irmd,lmpotd)


!---> loop over radial mesh


do ir = irc0,irc1
  rpoint = r(ir)
  
!---> calculate the ex.-cor. potential
  
  nspin2 = 2
  
  do ispin = 1,nspin2
    do lm = 1,lmmax
      rholm(lm,ispin) = rhol(ir,ispin,lm)
    enddo
  enddo
  
!    only for spin-polarized
  
!        write(6,*) ' before mkxcpe '
if (kxc==3) then
!  write(*,*) 'kxc=3: using gga functional pw91'  
  call mkxcpe_pw91(nspin2,ir,ijend,l1max,rpoint,rholm,vxc,excij,thet,  &
      ylm,dylmt1,dylmt2,dylmf1,dylmf2,dylmtf,drrl,ddrrl, drrul,ddrrul,irmd,lmpotd)
else if (kxc==4) then
!  write(*,*) 'kxc=4: using gga functional pbe'
  call mkxcpe_pbe(ir,ijend,rpoint,rholm,vxc,excij,ylm,dylmt1,  &
               dylmf1,dylmf2,dylmtf,drrl,ddrrl,drrul,ddrrul,  &
               irmd,lmpotd,lmmax,.false.)
else if (kxc==5) then
!  write(*,*) 'kxc=5: using gga functional pbesol'
  call mkxcpe_pbe(ir,ijend,rpoint,rholm,vxc,excij,ylm,dylmt1,  &
               dylmf1,dylmf2,dylmtf,drrl,ddrrl,drrul,ddrrul,  &
               irmd,lmpotd,lmmax,.true.)
else
write(*,*) 'error in vxcgga: in order to use gga set kxc=3 (pw91), kxc=4 (pbe) &
or kxc=5 (pbesol). lda is used for kxc<3)'
stop
endif 
  
  
  
  
!---> expand the ex.-cor. potential into spherical harmonics ,
!       using the orthogonality
  
  do ispin = 1,nspin
    
!---> determine the corresponding potential number
    
    ipot = ispin
    do lm = 1,lmmax
      vlmxc = ddot(ijend,vxc(1,ispin),1,wtyr(1,lm),1)
      v(ir,lm,ipot) = v(ir,lm,ipot) + vlmxc
      
!---> store the ex.-c. potential of ir=2 and =3 for the extrapolation
      
      if (lm == 1 .and. (ir == 2.or.ir == 3)) vxcr(ir, ispin) = vlmxc
    enddo
  enddo
  
!---> file er in case of total energies
  
  if (kte == 1) then
    
!---> expand ex.-cor. energy into spherical harmonics
!       using the orthogonality
    
    do l = 0,lpot
      do m = -l,l
        lm = l*l + l + m + 1
        elmxc = ddot(ijend,excij,1,wtyr(1,lm),1)
        
!---> multiply the lm-component of the ex.-cor. energy with the same
!     lm-component of the charge density times r**2 and sum over lm
!     this corresponds to a integration over the angular .
        
        if (ir > irs1) then
          estor(ir,lm) = elmxc
          
        else
          
          er(ir,l) = er(ir,l) + rho2ns(ir,lm,1)*elmxc
        endif
        
      enddo
      
    enddo
    
  endif
  
enddo

!---> integrate er in case of total energies to get exc

if (kte == 1) then
  
  do l = 0,lpot
    do m = -l,l
      lm = l*l + l + m + 1
      
!---> convolute with shape function
      
      do j = imaxsh(lm-1) + 1,imaxsh(lm)
        lm2 = ilm(j,2)
        if (lmsp(ilm(j,3)) > 0) then
          ifun = ifunm(ilm(j,3))
          do ir = irs1 + 1,irc1
            irh = ir - irs1
            er(ir,l) = er(ir,l) + rho2ns(ir,lm,1)*gsh(j)*  &
                thetas(irh,ifun)*estor(ir,lm2)
          enddo
        endif
      enddo
    enddo
!   call simpk(er(1,l),exc(l),ipan1,ircut,drdi)
    exc(l) = simpson(er(1:,l), ipan1, ircut, drdi)
  enddo
  
endif

!---> extrapolate ex.-cor potential to the origin only for lm=1

do ispin = 1,nspin
  ipot = ispin
  
  vxc2 = vxcr(2,ispin)
  vxc3 = vxcr(3,ispin)
  vxc1 = vxc2 - r(2)* (vxc3-vxc2)/ (r(3)-r(2))
  
  v(1,1,ipot) = v(1,1,ipot) + vxc1
enddo

end subroutine vxcgga
