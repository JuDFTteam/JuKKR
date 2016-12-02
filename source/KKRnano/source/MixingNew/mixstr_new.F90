! ATTENTION: the spherical part of the potential is divided by sqrt(4*pi) here

! 13.10.95 ***************************************************************
      subroutine mixstr_new(rmsavq, rmsavm, lmpot, nspin, mixing, fcm, irc1, irmin1, r, drdi, vons, visp, vins, irmd, irnsd)
      implicit none

      integer, intent(in) :: irmd
      integer, intent(in) :: irnsd
      double precision, intent(in) :: fcm, mixing
      double precision, intent(out) :: rmsavm, rmsavq
      integer, intent(in) :: lmpot, nspin
      double precision, intent(in) :: drdi(irmd), r(irmd), vins((irmd-irnsd):irmd,lmpot,2), visp(irmd,2)
      double precision, intent(inout) :: vons(irmd,lmpot,2)
!     ..
!     .. local scalars ..
      double precision :: fac, rmserm, rmserq, vmn, vnm, vnp, voldm, voldp, vpn, pi, fpi, rfpi
      integer :: ih, ihp1, irc1, irmin1, j, lm
!     ..
!
!---> final construction of the potentials
!     attention : the spherical averaged potential is the lm=1 component of vons times sqrt(4 pi).
!
!     first mixing scheme : straight mixing
!---> determination of the root mean sqare error
!
      rmsavq = 0.d0
      rmsavm = 0.d0

      pi = 4.d0*atan(1.d0)
      fpi = 4.d0*pi
      rfpi = sqrt(fpi)

        if (nspin == 2) then
          ih = 1
          ihp1 = ih + 1
        else
          ih = 1
          ihp1 = ih
        endif ! nspin == 2

!       ok, that means that for
!       nspin=1 (or not 2)  -> ih=1 and ihp1=1
!       nspin=2             -> ih=1 and ihp1=2
!       ihp1 means ih plus 1

        rmserq = 0.d0
        rmserm = 0.d0
        fac = 0.5d0/rfpi
!
        do j = 1, irc1
          vnp = fac*(vons(j,1,ih) + vons(j,1,ihp1))
          vnm = fac*(vons(j,1,ih) - vons(j,1,ihp1))
          voldp = 0.5d0* (visp(j,ih)+visp(j,ihp1))
          voldm = 0.5d0* (visp(j,ih)-visp(j,ihp1))
          rmserq = rmserq + 2.d0*real(1+mod(j,2))*r(j)*r(j)*drdi(j)*(vnp - voldp)*(vnp - voldp)
          rmserm = rmserm + 2.d0*real(1+mod(j,2))*r(j)*r(j)*drdi(j)*(vnm - voldm)*(vnm - voldm)
          vpn = voldp + mixing*(vnp - voldp)
          vmn = voldm + fcm*mixing*(vnm - voldm)
          vons(j,1,ihp1) = vpn - vmn
          vons(j,1,ih) = vpn + vmn
        enddo ! j
!
        rmserq = rmserq/(r(irc1)**3)
        rmserm = rmserm/(r(irc1)**3)
        rmsavq = rmsavq + rmserq
        rmsavm = rmsavm + rmserm

        if (lmpot > 1) then

          rmserq = 0.d0
          rmserm = 0.d0

          do lm = 2,lmpot
            do j = irmin1,irc1
              vnp = 0.5d0*(vons(j,lm,ih) + vons(j,lm,ihp1))
              vnm = 0.5d0*(vons(j,lm,ih) - vons(j,lm,ihp1))
              voldp = 0.5d0*(vins(j,lm,ih) + vins(j,lm,ihp1))
              voldm = 0.5d0*(vins(j,lm,ih) - vins(j,lm,ihp1))
              rmserq = rmserq + 2.d0*real(1+mod(j,2))*r(j)*r(j)*drdi(j)*(vnp - voldp)*(vnp - voldp)
              rmserm = rmserm + 2.d0*real(1+mod(j,2))*r(j)*r(j)*drdi(j)*(vnm - voldm)*(vnm - voldm)
              vpn = voldp + mixing*(vnp - voldp)
              vmn = voldm + fcm*mixing*(vnm - voldm)
              vons(j,lm,ihp1) = vpn - vmn
              vons(j,lm,ih) = vpn + vmn
            enddo ! j
          enddo ! lm
          rmserq = rmserq/(r(irc1)**3)/fpi
          rmserm = rmserm/(r(irc1)**3)/fpi
          rmsavq = rmsavq + rmserq
          rmsavm = rmsavm + rmserm

        endif ! lmpot > 1

      endsubroutine ! mixstr_new
