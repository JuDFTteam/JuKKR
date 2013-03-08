!>    @param print_info  0=print some info  1=print nothing
subroutine madelung3d(lpot,yrg,wg,alat, &
     rmax,gmax,bravais,recbv, &
     lmxspd,lassld,lpotd,lmpotd, &
     nmaxd,ishld, &
     lmpot,cleb,icleb,iend, &
     nclebd,loflm,dfac, &
     ngmax,nrmax,nsg,nsr,nshlg,nshlr,gn,rm, &
     print_info)
  ! **********************************************************************
  ! *                                                                    *
  ! * This subroutine calculates the Madelung potential coefficients     *
  ! *                                                                    *
  ! **********************************************************************
  implicit none
  !     ..
  !     .. Scalar Arguments ..
  integer print_info

  integer lpot
  integer lmxspd,lassld,lpotd,lmpotd,nmaxd,ishld
  double precision alat,rmax,gmax
  !     ..
  !     .. Array Arguments ..
  double precision yrg(lassld,0:lassld,0:lassld),wg(lassld)
  double precision bravais(3,3),recbv(3,3)
  !     ..
  !     .. Local Scalars ..
  integer i,iend,iprint,l,m,l1,l2,lmpot
  integer nclebd
  integer ngmax,nrmax,nshlg,nshlr

  double precision pi,fpi

  !     ..
  !     .. Local Arrays ..
  !     .. Attention: Dimension LMXSPD*LMPOTD appears sometimes as NCLEB1
  double precision cleb(lmxspd*lmpotd)
  double precision gn(3,nmaxd),rm(3,nmaxd)
  double precision dfac(0:lpotd,0:lpotd)
  integer icleb(lmxspd*lmpotd,3)
  integer nsg(ishld),nsr(ishld)
  integer loflm(lmxspd)
  !     ..
  !     .. External subroutines
  external lattice3d,madelgaunt
  ! ......................................................................
  iprint = 0
  nclebd = lmxspd*lmpotd
  pi = 4.0d0*atan(1.0d0)
  fpi = 4.0d0*pi
  lmpot = (lpot+1)**2
  i = 1
  !
  ! --> determine the l-value for given lm
  !
  do l = 0,2*lpot
     do m = -l,l
        loflm(i) = l
        i = i + 1
     end do
  end do
  !
  ! --> calculate:                             (2*(l+l')-1)!!
  !                 dfac(l,l') = 4pi**2 *  ----------------------
  !                                        (2*l+1)!! * (2*l'+1)!!
  !
  dfac(0,0) = fpi*fpi
  do l1 = 1,lpot
     dfac(l1,0) = dfac(l1-1,0)*dble(2*l1-1)/dble(2*l1+1)
     dfac(0,l1) = dfac(l1,0)
     do l2 = 1,l1
        dfac(l1,l2) = dfac(l1,l2-1)*dble(2*(l1+l2)-1)/dble(2*l2+1)
        dfac(l2,l1) = dfac(l1,l2)
     end do
  end do
  !
  ! --> calculate the gaunt coefficients
  !
  call madelgaunt(lpot,yrg,wg,cleb,icleb,iend,lassld,nclebd)
  !
  ! ======================================================================
  call lattice3d(alat,bravais,recbv,ngmax,nrmax,nshlg,nshlr,nsg,nsr, &
       rmax,gmax,gn,rm,iprint,nmaxd,ishld, print_info)
  !
  ! ======================================================================
  !
end subroutine madelung3d
