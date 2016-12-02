! rho2ns = r^2 n_L(r)

! 13.10.95 ***************************************************************
subroutine rhomom_new(cmom, cminst, lpot, rho2ns, r, drdi, ircut, ipan, ilm, ifunm, imaxsh, gsh, thetas, lmsp, irmd, irid, nfund, ipand, ngshd)

  ! ************************************************************************
  !     calculate charge moments of given charge densities
  !
  !                             rcut
  !              cmom(lm,i) =    s dr' r'** l rho2ns(r',lm,i,1)
  !                              0
  !-----------------------------------------------------------------------
  implicit none
  integer, intent(in) :: irmd
  integer, intent(in) :: irid
  integer, intent(in) :: nfund
  integer, intent(in) :: ipand
  integer, intent(in) :: ngshd

  !     INTEGER LMPOTD
  !     PARAMETER (LMPOTD= (LPOTD+1)**2) ! = (2*LMAX+1)**2
  !     ..
  !     .. Scalar Arguments ..
  integer, intent(in) :: lpot
  !     ..
  double precision, intent(out) :: cminst((lpot+1)**2)
  double precision, intent(out) :: cmom((lpot+1)**2)
  double precision, intent(in) :: drdi(irmd)
  double precision, intent(in) :: gsh(*)
  double precision, intent(in) :: r(irmd)
  double precision, intent(in) :: rho2ns(irmd,(lpot+1)**2)
  double precision, intent(in) :: thetas(irid,nfund)
  integer, intent(in) :: ifunm(*)
  integer, intent(in) :: ilm(ngshd,3)
  integer, intent(in) :: imaxsh(0:(lpot+1)**2)
  integer, intent(in) :: ipan
  integer, intent(in) :: ircut(0:ipand)
  integer, intent(in) :: lmsp(*)
  !     ..
  !     .. Locals ..
  double precision :: rl
  integer :: i,iend,ifun,irc1,irs1,istart,j,l,lm,lm2,lm3,m
  double precision :: v1(irmd), vint1(irmd)

  external :: sinwk, soutk

  v1 = 0.0d0
  vint1 = 0.0d0

  irs1 = ircut(1)
  irc1 = ircut(ipan)

  do l = 0, lpot
    do m = -l, l
      lm = l*l + l + m + 1
      !
      !---> set up of the integrands v1 and v2
      !
      v1(1) = 0.0d0
      do i = 2, irs1
        rl = r(i)**l
        v1(i) = rho2ns(i,lm)*rl*drdi(i)
      enddo ! i
      !
      !---> convolute charge density of interstial with shape function
      !
      v1(irs1+1:irc1) = 0.0d0

      istart = imaxsh(lm-1) + 1
      iend = imaxsh(lm)

      do j = istart, iend
        lm2 = ilm(j,2)
        lm3 = ilm(j,3)
        if (lmsp(lm3) > 0) then
          ifun = ifunm(lm3)
          do i = irs1 + 1,irc1
            v1(i) = v1(i) + gsh(j)*rho2ns(i,lm2)*thetas(i-irs1,ifun)
          enddo ! i
        endif
      enddo ! j

      do i = irs1 + 1, irc1
        rl = r(i)**l
        v1(i) = v1(i)*rl*drdi(i)
      enddo ! i
      !
      !---> now integrate
      !
      call soutk(v1,vint1,ipan,ircut)

      cmom(lm) = vint1(irs1)
      cminst(lm) = vint1(irc1) - vint1(irs1)

      enddo ! m
    enddo ! l

  endsubroutine rhomom_new
