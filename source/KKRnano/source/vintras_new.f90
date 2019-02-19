! 13.10.95 ***************************************************************
subroutine vintras_new(lpot,nspin,rho2ns,vons,r, &
drdi,ircut,ipan,ilm,ifunm, &
imaxsh,gsh,thetas,lmsp, &
irmd, irid, nfund, ngshd, ipand)
  ! ************************************************************************
  !     calculate the electron-intracell-potentials and the charge-
  !     moments of given charge densities . ( for each spin-direc-
  !     tion the potential is the same in the polarized case . )
  !     initialize the potential v with the electron-intracell-potentials
  !     the intracell-potential is expanded into spherical harmonics .
  !     the lm-term of the intracell-potential of the representive atom i
  !     is given by
  !                    8pi        r      r'** l
  !      v(r,lm,i) =  ----- *  (  s dr' --------   rho2ns(r',lm,i,1)
  !                   2*l+1       0     r **(l+1)
  !
  !                                 rcut    r ** l
  !                               +  s dr' ---------   rho2ns(r',lm,i,1) )
  !                                  r     r' **(l+1)
  !
  !             (see notes by b.drittler and u.klemradt)
  !
  !     attention : rho2ns(...,1) is the real charge density times r**2
  !                 developed into spherical harmonics . (see deck rholm)
  !
  !                               b.drittler   may 1987
  !-----------------------------------------------------------------------
  implicit none
  !     .. Parameters ..

  integer irmd
  integer irid
  integer nfund
  integer ngshd
  integer ipand

  !
  !     INTEGER LMPOTD
  !     PARAMETER (LMPOTD= (LPOTD+1)**2)
  !     ..
  !     .. Scalar Arguments ..
  integer lpot,nspin
  !     ..
  !     .. Array Arguments ..

  double precision drdi(irmd)
  double precision gsh(*)
  double precision r(irmd)
  !     DOUBLE PRECISION RHO2NS(IRMD,LMPOTD)
  double precision rho2ns(irmd,(lpot+1)**2)

  double precision thetas(irid,nfund)
  !     DOUBLE PRECISION VONS(IRMD,LMPOTD,2)
  double precision vons(irmd,(lpot+1)**2,2)

  integer ifunm(*)
  integer ilm(ngshd,3)
  !     INTEGER IMAXSH(0:LMPOTD)
  integer imaxsh(0:(lpot+1)**2)
  integer ipan
  integer ircut(0:ipand)
  integer lmsp(*)

  !     ..
  !     .. Local Scalars ..
  double precision fac,pi,rl
  integer i,iend,ifun,ipot,irc1,irs1,istart,j,l,lm,lm2, &
  lm3,m
  !     ..
  !     .. Local Arrays ..
  !     Fortran 90 automatic arrays
  double precision v1(irmd),v2(irmd),vint1(irmd),vint2(irmd)
  integer ircutm(0:ipand)
  !     ..
  !     .. External Subroutines ..
  external sinwk,soutk
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic atan,real
  !     ..

  v1 = 0.0d0
  v2 = 0.0d0
  vint1 = 0.0d0
  vint2 = 0.0d0
  ircutm = 0

  pi = 4.d0*atan(1.d0)
  !
  irs1 = ircut(1)
  irc1 = ircut(ipan)

  do 10 i = 0,ipan
    ircutm(i) = ircut(i)
10 continue

   !---> determine the right potential numbers
   ipot = nspin

   do 100 l = 0,lpot
     fac = 8.0d0*pi/real(2*l+1)
     do 90 m = -l,l
       lm = l*l + l + m + 1
       !
       !---> set up of the integrands v1 and v2
       !
       v1(1) = 0.0d0
       v2(1) = 0.0d0
       do 20 i = 2,irs1
         rl = r(i)**l
         v1(i) = rho2ns(i,lm)*rl*drdi(i)
         v2(i) = rho2ns(i,lm)/r(i)/rl*drdi(i)
20     continue
       !
       !---> convolute charge density of interstial with shape function
       !

       do 30 i = irs1 + 1,irc1
         v1(i) = 0.0d0
30     continue

       istart = imaxsh(lm-1) + 1
       iend = imaxsh(lm)

       do 50 j = istart,iend
         lm2 = ilm(j,2)
         lm3 = ilm(j,3)

         if (lmsp(lm3).gt.0) then
           ifun = ifunm(lm3)
           do 40 i = irs1 + 1,irc1
             v1(i) = v1(i) + gsh(j)*rho2ns(i,lm2)* &
             thetas(i-irs1,ifun)
40         continue
         end if

50     continue

       do 60 i = irs1 + 1,irc1
         rl = r(i)**l
         v2(i) = v1(i)/r(i)/rl*drdi(i)
         v1(i) = v1(i)*rl*drdi(i)
60     continue
       !
       !---> now integrate v1 and v2
       !
       call soutk(v1,vint1,ipan,ircutm)  ! integrals from 0 to r    (r varies)
       call sinwk(v2,vint2,ipan,ircutm)  ! integrals from r to r_BS (r varies)
       !
       !---> gather all parts

       ! deal with r=0
       if (lm.eq.1) then
         vons(1,lm,ipot) = fac*vint2(1)
       else
         vons(1,lm,ipot) = 0.0d0
       end if

       do 70 i = 2,irc1
         rl = r(i)**l
         vons(i,lm,ipot) = fac* (vint1(i)/r(i)/rl+vint2(i)*rl)
70     continue

       ! intra-cell potential is the same for the other spin direction
       if (nspin.eq.2) then
         do 80 i = 1,irc1
           vons(i,lm,ipot-1) = vons(i,lm,ipot)
80       continue
       end if

90   continue ! loop over M

100 continue  ! loop over L

    return

  end
