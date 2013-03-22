subroutine irwns(cr,dr,efac,qns,vnspll,icst,ipan,ircut,nsra, &
pzlm,qzlm,pzekdr,qzekdr,cder,cmat,dder,dmat, &
irmind,irmd,ipand,lmmaxd)
  implicit none
  !-----------------------------------------------------------------------
  !     determines the irregular non spherical wavefunctions in the n-th.
  !       born approximation ( n given by input parameter icst ) .
  !
  !
  !     using the wave functions pz and qz ( regular and irregular
  !       solution ) of the spherically averaged potential , the ir-
  !       regular wavefunction qns is determined by
  !
  !           qns(ir,lm1,lm2) = cr(ir,lm1,lm2)*pz(ir,l1)
  !
  !                                   + dr(ir,lm1,lm2)*qz(ir,l1)
  !
  !      the matrices cr and dr are determined by integral equations
  !        containing qns and only the non spherical contributions of
  !        the potential , stored in vinspll . these integral equations
  !        are solved iteratively with born approximation up to given n.
  !
  !     the original way of writing the cr and dr matrices in the equa-
  !        tion above caused numerical troubles . therefore here are used
  !        rescaled cr and dr matrices (compare subroutine wftsca):
  !
  !              ~
  !              cr(ir,lm1,lm2) = sqrt(e)**(l1+l2)
  !                             * cr(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
  !
  !              ~
  !              dr(ir,lm1,lm2) = sqrt(e)**(l2-l1)
  !                             * dr(ir,lm1,lm2)*((2*l1-1)!!/(2*l2-1)!!)
  !
  !     attention :  the sign of the dr matrix is changed to reduce the
  !     ===========  number of floating point operations
  !
  !     modified for the use of shape functions
  !
  !                              (see notes by b.drittler)
  !
  !                                b.drittler   mar.  1989
  !-----------------------------------------------------------------------
  !     modified by R. Zeller      Aug. 1994
  !-----------------------------------------------------------------------
  !     .. Parameters ..
  double complex CONE
  parameter (CONE= (1.d0,0.d0))
  !     ..
  !     .. Scalar Arguments ..
  integer icst,ipan,ipand,irmd,irmind,lmmaxd,nsra
  !     ..
  !     .. Array Arguments ..
  double complex cder(lmmaxd,lmmaxd,irmind:irmd), &
  cmat(lmmaxd,lmmaxd,irmind:irmd),cr(lmmaxd,lmmaxd), &
  dder(lmmaxd,lmmaxd,irmind:irmd), &
  dmat(lmmaxd,lmmaxd,irmind:irmd),dr(lmmaxd,lmmaxd), &
  efac(lmmaxd),pzekdr(lmmaxd,irmind:irmd,2), &
  pzlm(lmmaxd,irmind:irmd,2), &
  qns(lmmaxd,lmmaxd,irmind:irmd,2), &
  qzekdr(lmmaxd,irmind:irmd,2), &
  qzlm(lmmaxd,irmind:irmd,2)

  double precision vnspll(lmmaxd,lmmaxd,irmind:irmd)
  integer ircut(0:ipand)
  !     ..
  !     .. Local Scalars ..
  double complex efac2
  integer i,ir,irc1,j,lm1,lm2
  !     ..
  !     .. External Subroutines ..
  external csinwd,wfint,wfint0
  !     ..
  irc1 = ircut(ipan)

  do 70 i = 0,icst
    !---> set up integrands for i-th born approximation
    if (i.eq.0) then
      call wfint0(cder,dder,qzlm,qzekdr,pzekdr,vnspll,nsra,irmind, &
      irmd,lmmaxd)
    else
      call wfint(qns,cder,dder,qzekdr,pzekdr,vnspll,nsra,irmind, &
      irmd,lmmaxd)
    end if

    !---> call integration subroutines
    call csinwd(cder,cmat,lmmaxd**2,irmind,irmd,ipan,ircut)
    call csinwd(dder,dmat,lmmaxd**2,irmind,irmd,ipan,ircut)
    do 20 ir = irmind,irc1
      do 10 lm2 = 1,lmmaxd
        dmat(lm2,lm2,ir) = dmat(lm2,lm2,ir) - CONE
10    continue
20  continue

    !---> calculate non sph. wft. in i-th born approximation
    do 60 j = 1,nsra
      do 50 ir = irmind,irc1
        do 40 lm1 = 1,lmmaxd
          do 30 lm2 = 1,lmmaxd
            qns(lm1,lm2,ir,j) = cmat(lm1,lm2,ir)*pzlm(lm1,ir,j) - &
            dmat(lm1,lm2,ir)*qzlm(lm1,ir,j)
30        continue
40      continue
50    continue
60  continue
70 continue

   do 90 lm2 = 1,lmmaxd
     !---> store c - and d - matrix
     do 80 lm1 = 1,lmmaxd
       cr(lm1,lm2) = cmat(lm1,lm2,irmind)
       dr(lm1,lm2) = -dmat(lm1,lm2,irmind)
80   continue
90 continue

   !---> rescale with efac
   do 130 j = 1,nsra
     do 120 lm2 = 1,lmmaxd
       efac2 = 1.d0/efac(lm2)
       do 110 ir = irmind,irc1
         do 100 lm1 = 1,lmmaxd
           qns(lm1,lm2,ir,j) = qns(lm1,lm2,ir,j)*efac2
100      continue
110    continue
120  continue
130 continue

  end
