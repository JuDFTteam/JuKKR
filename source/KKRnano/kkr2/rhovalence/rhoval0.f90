subroutine rhoval0(ez,wez,drdi,r,ipan,ircut, &
thetas,dos0,dos1, &
lmaxd, irmd, irid, ipand, nfund)
  !
  implicit none
  !
  !     .. Parameters ..

  integer lmaxd
  integer irmd
  integer irid
  integer ipand
  integer nfund

  double complex CONE,CZERO,CI
  parameter ( CONE=(1.d0,0.d0),CZERO=(0.d0,0.d0),CI=(0.d0,1.d0) )
  !     ..
  !     .. Scalar Arguments ..
  integer ipan
  double complex ez,wez,dos0,dos1
  !     ..
  !     .. Array Arguments ..
  double precision drdi(irmd), &
  r(irmd), &
  thetas(irid,nfund)
  integer ircut(0:ipand)
  !     ..
  !     .. Local Scalars ..
  double complex ek,ciek,denl
  double precision c0ll
  integer ir,l,l1,imt1
  !     ..
  !     .. Local Arrays ..
  double complex pz(irmd,0:lmaxd),qz(irmd,0:lmaxd)
  double complex bessjw(0:lmaxd+1)
  double complex bessyw(0:lmaxd+1)
  double complex hankws(0:lmaxd+1)
  double complex cden0(irmd,0:lmaxd+1)
  double complex cden1(irmd,0:lmaxd+1)
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic atan,dble,sqrt

  integer lmaxd1

  lmaxd1 = lmaxd+1
  !     ..
  !
  ek = sqrt(ez)
  c0ll = 1.0d0/sqrt(16.0d0*atan(1.0d0))
  ciek=(0.0d0,1.0d0)*ek
  !
  ! initialize arrays ...
  !
  do ir = 1, irmd
    do l1 = 0, lmaxd1
      cden0(ir,l1) = CZERO
      cden1(ir,l1) = CZERO
    enddo
  enddo
  !
  !
  !=======================================================================
  do ir = 2,irmd
    call beshan(hankws,bessjw,bessyw,r(ir)*ek,lmaxd1)
    do l = 0,lmaxd
      pz(ir,l) = bessjw(l)*r(ir)
      qz(ir,l) = (bessyw(l) - CI*bessjw(l))*r(ir)
    enddo
  enddo

  imt1=ircut(1)
  do l1 = 0,lmaxd1
    cden0(1,l1) = (0.0d0,0.0d0)
    cden1(1,l1) = (0.0d0,0.0d0)
  end do

  do ir = 2,irmd
    cden0(ir,0) = ek*pz(ir,0)*qz(ir,0)
    cden1(ir,0) = ek*pz(ir,0)**2*(0.d0,-1.d0)
    cden1(ir,lmaxd1) = ciek*r(ir)**2
  end do

  do l1 = 1,lmaxd
    do ir = 2,irmd
      cden0(ir,l1) = ek*pz(ir,l1)*qz(ir,l1)*(l1+l1+1)
      cden1(ir,l1) = ek*pz(ir,l1)**2*(0.d0,-1.d0)*(l1+l1+1)
    end do  
  end do

  do l1 = 0,lmaxd1
    if (ipan.gt.1) then
      do ir = imt1 + 1,irmd
        cden0(ir,l1) = cden0(ir,l1)*thetas(ir-imt1,1)*c0ll
        cden1(ir,l1) = cden1(ir,l1)*thetas(ir-imt1,1)*c0ll
      end do
    end if
    call csimpk(cden0(1,l1),denl,ipan,ircut,drdi)
    dos0 = dos0 + wez*denl
    call csimpk(cden1(1,l1),denl,ipan,ircut,drdi)
    dos1 = dos1 + wez*denl
  end do

end
