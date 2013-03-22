subroutine rhoout(cden,df,gmat,ek,pns,qns,rho2ns,thetas,ifunm, &
ipan1,imt1,lmsp,cdenns,nsra,cleb,icleb,iend, &
lmaxd, irmd, irnsd, irid, nfund, ncleb)
  !-----------------------------------------------------------------------
  !
  !     calculates the charge density from r(irmin) to r(irc)
  !      in case of a non spherical input potential .
  !
  !     fills the array cden for the complex density of states
  !
  !     attention : the gaunt coeffients are stored in index array
  !                   (see subroutine gaunt)
  !
  !     the structured part of the greens-function (gmat) is symmetric in
  !       its lm-indices , therefore only one half of the matrix is
  !       calculated in the subroutine for the back-symmetrisation .
  !       the gaunt coeffients are symmetric too (since the are calculated
  !       using the real spherical harmonics) . that is why the lm2- and
  !       the lm02- loops are only only going up to lm1 or lm01 and the
  !       summands are multiplied by a factor of 2 in the case of lm1 .ne.
  !       lm2 or lm01 .ne. lm02 .
  !
  !             (see notes by b.drittler)
  !
  !                               b.drittler   aug. 1988
  !-----------------------------------------------------------------------
  !     .. Parameters ..

  !     INTEGER LMMAXD
  !     PARAMETER (LMMAXD= (LMAXD+1)**2)
  !     INTEGER LMPOTD
  !     PARAMETER (LMPOTD= (LPOTD+1)**2) ! = (2*LMAXD+1)**2
  !     INTEGER IRMIND
  !     PARAMETER (IRMIND=IRMD-IRNSD)
  !     ..
  !     .. Scalar Arguments ..

  implicit none

  integer lmaxd
  integer irmd
  integer ncleb
  integer irnsd
  integer irid
  integer nfund

  double complex df,ek
  integer iend,imt1,ipan1,nsra
  !     ..
  !     .. Array Arguments ..
  !     DOUBLE COMPLEX CDEN(IRMD,0:*),CDENNS(*),GMAT(LMMAXD,LMMAXD),
  !    +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
  !    +               QNSI(LMMAXD,LMMAXD),
  !    +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
  !     DOUBLE PRECISION CLEB(*),RHO2NS(IRMD,LMPOTD),THETAS(IRID,NFUND)

  double complex cden(irmd,0:*)
  double complex cdenns(*)
  double complex gmat((lmaxd+1)**2,(lmaxd+1)**2)
  double complex pns((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd,2)
  double complex qnsi((lmaxd+1)**2,(lmaxd+1)**2)
  double complex qns((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd,2)
  double precision cleb(*)
  double precision rho2ns(irmd,(2*lmaxd+1)**2)
  double precision thetas(irid,nfund)

  integer icleb(ncleb,3),ifunm(*),lmsp(*)

  !     .. Local Scalars ..
  double complex cltdf,cone,czero
  double precision c0ll
  integer i,ifun,ir,j,l1,lm1,lm2,lm3,m1
  !     ..
  !     .. Local Arrays ..
  !     DOUBLE COMPLEX WR(LMMAXD,LMMAXD,IRMIND:IRMD)
  double complex wr((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
  !     ..
  !     .. External Subroutines ..
  external zgemm
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic atan,dimag,sqrt
  !     ..
  !     .. Save statement ..
  !     SAVE
  !     ..
  !     .. Data statements ..
  data czero/ (0.0d0,0.0d0)/
  data cone/ (1.0d0,0.0d0)/

  integer irmind
  integer lmmaxd

  lmmaxd= (lmaxd+1)**2
  irmind=irmd-irnsd

  !     C0LL = 1/sqrt(4*pi)
  c0ll = 1.0d0/sqrt(16.0d0*atan(1.0d0))
  !
  !
  !---> initialize array for complex charge density
  !
  do l1 = 0,lmaxd
    do i = 1,irmd
      cden(i,l1) = czero
    enddo
  enddo
  !------------------------------------------------------------------
  !
  !---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
  !                                      summed over lm3
  !---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
  !                                               summed over lm3
  do 50 ir = irmind + 1,irmd
    do lm1=1,lmmaxd
      do lm2=1,lmmaxd
        qnsi(lm1,lm2)=qns(lm1,lm2,ir,1)
      enddo
    enddo

    call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,pns(1,1,ir,1), &
    lmmaxd,gmat,lmmaxd,ek,qnsi,lmmaxd)

    call zgemm('N','T',lmmaxd,lmmaxd,lmmaxd,cone,pns(1,1,ir,1), &
    lmmaxd,qnsi,lmmaxd,czero,wr(1,1,ir),lmmaxd)

    if (nsra.eq.2) then
      do lm1=1,lmmaxd
        do lm2=1,lmmaxd
          qnsi(lm1,lm2)=qns(lm1,lm2,ir,2)
        enddo
      enddo

      call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,pns(1,1,ir,2), &
      lmmaxd,gmat,lmmaxd,ek,qnsi,lmmaxd)

      call zgemm('N','T',lmmaxd,lmmaxd,lmmaxd,cone,pns(1,1,ir,2), &
      lmmaxd,qnsi,lmmaxd,cone,wr(1,1,ir),lmmaxd)

    end if

    do lm1 = 1,lmmaxd
      do lm2 = 1,lm1 - 1
        wr(lm1,lm2,ir) = wr(lm1,lm2,ir) + wr(lm2,lm1,ir)
      enddo
    enddo
50 CONTINUe
   !
   !---> first calculate only the spherically symmetric contribution
   !
   do 100 l1 = 0,lmaxd
     do 70 m1 = -l1,l1
       lm1 = l1* (l1+1) + m1 + 1
       do 60 ir = irmind + 1,irmd
         !
         !---> fill array for complex density of states
         !
         cden(ir,l1) = cden(ir,l1) + wr(lm1,lm1,ir)
60     continue
70   continue
     !
     !---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
     !
     do 80 ir = irmind + 1,irmd
       rho2ns(ir,1) = rho2ns(ir,1) + c0ll*dimag(cden(ir,l1)*df)
80   continue
     !
     !
     !
     if (ipan1.gt.1) then
       do 90 i = imt1 + 1,irmd
         cden(i,l1) = cden(i,l1)*thetas(i-imt1,1)*c0ll
90     continue
     end if

100 continue


    if (ipan1.gt.1) then
      do 110 i = 1,irmd
        cdenns(i) = 0.0d0
110   continue
    end if

    do 140 j = 1,iend
      lm1 = icleb(j,1)
      lm2 = icleb(j,2)
      lm3 = icleb(j,3)
      cltdf = df*cleb(j)
      !
      !---> calculate the non spherically symmetric contribution
      !
      do 120 ir = irmind + 1,irmd
        rho2ns(ir,lm3) = rho2ns(ir,lm3) + dimag(cltdf*wr(lm1,lm2,ir))
120   continue

      if (ipan1.gt.1 .and. lmsp(lm3).gt.0) then
        !       IF (IPAN1.GT.1) THEN
        ifun = ifunm(lm3)
        do 130 i = imt1 + 1,irmd
          cdenns(i) = cdenns(i) + cleb(j)*wr(lm1,lm2,i)* &
          thetas(i-imt1,ifun)
130     continue

      end if

140 continue


  end
