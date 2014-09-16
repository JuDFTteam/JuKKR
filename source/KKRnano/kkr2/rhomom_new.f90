! rho2ns = r^2 n_L(r)

! 13.10.95 ***************************************************************
subroutine rhomom_new(cmom,cminst,lpot,rho2ns,r, &
drdi,ircut,ipan,ilm,ifunm, &
imaxsh,gsh,thetas,lmsp, &
irmd, irid, nfund, ipand, ngshd)

  ! ************************************************************************
  !     calculate charge moments of given charge densities
  !
  !                             rcut
  !              cmom(lm,i) =    s dr' r'** l rho2ns(r',lm,i,1)
  !                              0
  !-----------------------------------------------------------------------
  implicit none
  integer irmd
  integer irid
  integer nfund
  integer ipand
  integer ngshd

  !     INTEGER LMPOTD
  !     PARAMETER (LMPOTD= (LPOTD+1)**2) ! = (2*LMAX+1)**2
  !     ..
  !     .. Scalar Arguments ..
  integer lpot
  !     ..
  double precision cminst((lpot+1)**2)
  double precision cmom((lpot+1)**2)
  double precision drdi(irmd)
  double precision gsh(*)
  double precision r(irmd)
  double precision rho2ns(irmd,(lpot+1)**2)
  double precision thetas(irid,nfund)
  integer ifunm(*)
  integer ilm(ngshd,3)
  integer imaxsh(0:(lpot+1)**2)
  integer ipan
  integer ircut(0:ipand)
  integer lmsp(*)
  !     ..
  !     .. Local Scalars ..
  double precision rl
  integer i,iend,ifun,irc1,irs1,istart,j,l,lm,lm2, &
  lm3,m
  !     ..
  !     .. Local Arrays ..
  double precision v1(irmd),vint1(irmd)

  !     .. External Subroutines ..
  external sinwk,soutk
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic atan,real
  !     ..

  v1 = 0.0d0
  vint1 = 0.0d0

  irs1 = ircut(1)
  irc1 = ircut(ipan)

  do 100 l = 0,lpot

    do 90 m = -l,l
      lm = l*l + l + m + 1
      !
      !---> set up of the integrands v1 and v2
      !
      v1(1) = 0.0d0
      do 20 i = 2,irs1
        rl = r(i)**l
        v1(i) = rho2ns(i,lm)*rl*drdi(i)
20    continue
      !
      !---> convolute charge density of interstial with shape function
      !
      do 30 i = irs1 + 1,irc1
        v1(i) = 0.0d0
30    continue

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
40        continue
        end if
50    continue

      do 60 i = irs1 + 1,irc1
        rl = r(i)**l
        v1(i) = v1(i)*rl*drdi(i)
60    continue
      !
      !---> now integrate
      !
      call soutk(v1,vint1,ipan,ircut)

      cmom(lm) = vint1(irs1)
      cminst(lm) = vint1(irc1) - vint1(irs1)

90  continue

100 continue

    return

  end
