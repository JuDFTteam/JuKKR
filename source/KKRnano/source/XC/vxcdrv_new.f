      subroutine vxcdrv_new(exc,kte,kxc,lpot,nspin,rho2ns,vons,
     +                  r,drdi,a,irws,ircut,ipan,gsh,ilm,
     +                  imaxsh,ifunm,thetas,lmsp,
     &                  irmd, irid, nfund, ngshd, ipand)
      implicit none

      integer, intent(in) :: irmd, irid, nfund, ngshd, ipand


c     fix: hardcoded
      integer, parameter :: ijd = 434
c     ..
c     .. scalar arguments ..
      integer, intent(in) :: kte, kxc, lpot, nspin
c     ..
c     .. array arguments ..
c     double precision :: a(naezd),drdi(irmd,*),exc(0:lpotd),gsh(*),
c    +                 r(irmd,*),rho2ns(irmd,lmpotd,2),
c    +                 thetas(irid,nfund,*),vons(irmd,lmpotd,2)
c     integer ifunm(*),ilm(ngshd,3),imaxsh(0:lmpotd),ipan(*),
c    +        ircut(0:ipand,*),irws(*),lmsp(*)

      double precision :: a
      double precision :: drdi(irmd)
      double precision :: exc(0:lpot)
      double precision :: gsh(*)
      double precision :: r(irmd)
      double precision :: rho2ns(irmd,(lpot+1)**2,2)
      double precision :: thetas(irid,nfund)
      double precision :: vons(irmd,(lpot+1)**2,2)
      integer ifunm(*)
      integer ilm(ngshd,3)
      integer imaxsh(0:(lpot+1)**2)
      integer ipan
      integer ircut(0:ipand)
      integer irws
      integer lmsp(*)
c     ..
c     .. external subroutines ..
      external :: sphere_gga, sphere_nogga, vxcgga, vxclm
c     ..
c     .. local arrays .. fortran 90 automatic arrays
c     double precision :: dylmf1(ijd,lmpotd),dylmf2(ijd,lmpotd),
c    +                 dylmt1(ijd,lmpotd),dylmt2(ijd,lmpotd),
c    +                 dylmtf(ijd,lmpotd),
c    +                 rij(ijd,3),thet(ijd),wtyr(ijd,lmpotd),
c    +                 ylm(ijd,lmpotd),yr(ijd,lmpotd)
c     integer ifunmiat(lmxspd)

      double precision :: dylmf1(ijd,(lpot+1)**2)
      double precision :: dylmf2(ijd,(lpot+1)**2)
      double precision :: dylmt1(ijd,(lpot+1)**2)
      double precision :: dylmt2(ijd,(lpot+1)**2)
      double precision :: dylmtf(ijd,(lpot+1)**2)
      double precision :: rij(ijd,3)
      double precision :: thet(ijd)
      double precision :: wtyr(ijd,(lpot+1)**2)
      double precision :: ylm(ijd,(lpot+1)**2)
      double precision :: yr(ijd,(lpot+1)**2)
c     ..
c     .. local scalars ..

      integer lmpotd

      lmpotd = (lpot+1)**2
c     ..
      if (kxc.lt.3) then
        call sphere_nogga(lpot,yr,wtyr,rij,ijd)
      else
        call sphere_gga(lpot,yr,wtyr,rij,ijd,lmpotd,thet,ylm,dylmt1,
     +                  dylmt2,dylmf1,dylmf2,dylmtf)
      end if

        if (kxc.lt.3) then

          call vxclm(exc,kte,kxc,lpot,nspin,rho2ns,
     +               vons,r,drdi,
     +               ircut,ipan,
     +               gsh,ilm,imaxsh,ifunm,thetas,
     +               yr,wtyr,ijd,lmsp,
     &               irmd, irid, nfund, ngshd, ipand)
        else
c
c gga ex-cor potential
c
          call vxcgga(exc,kte,lpot,nspin,rho2ns,
     +                vons,r,drdi,a,
     +                irws,ircut,ipan,
     +                gsh,ilm,imaxsh,ifunm,thetas,
     +                wtyr,ijd,lmsp,thet,ylm,dylmt1,dylmt2,
     +                dylmf1,dylmf2,dylmtf,
     &                irmd, irid, nfund, ngshd, ipand, kxc)
        end if
      end
