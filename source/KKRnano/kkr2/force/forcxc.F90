  subroutine forcxc(flm, flmc, lpot, nspin, rhoc,v,r, drdi, irws, irmd)
      use Quadrature_mod, only: simpson
      use Constants_mod, only: pi
      implicit none
!-----------------------------------------------------------------------
!     calculates the force on nucleus m
!     from a given non spherical charge density at the nucleus site r
!     with core correction(exchange contribution)
!-----------------------------------------------------------------------
!      integer, parameter :: lmpotd = (lpotd+1)**2
      integer, intent(in) :: lpot, nspin, irmd, irws
      double precision, intent(in) :: drdi(irmd), flmc(-1:1), r(irmd), rhoc(irmd,*), v(irmd,(lpot+1)**2,2)
      double precision, intent(inout) :: flm(-1:1)

!     .. locals ..
      double precision :: dv, fac, rws, vint1
      integer :: i, ipot, irws1, ispin, lm, m
      double precision :: flmxc(-1:1), v1(irmd), tail_cor(-1:1)
      
      fac = sqrt((4.d0*pi)/3.d0)

      if (lpot < 1) then
        write(6, fmt="('error stop in subroutine force : the charge density has to contain non spherical contributions up to l=1 at least')")
        stop
      endif

        tail_cor(:) = 0.d0

        irws1 = irws
        rws = r(irws1)
 
        do m = -1, 1
          lm = 2 + m + 1

          v1(1:irws1) = 0.d0

          do ispin = 1, nspin
              ipot = ispin

              i = 2
                dv = (-3.d0*v(i-1,lm,ipot) -10.d0*v(i,lm,ipot) +18.d0*v(i+1,lm,ipot) -6.d0*v(i+2,lm,ipot) + v(i+3,lm,ipot))/(12.d0*drdi(2))
                v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)
                
              do i = 3, irws1 - 2
                dv = (v(i-2,lm,ipot) - v(i+2,lm,ipot) + 8.d0*(v(i+1,lm,ipot) - v(i-1,lm,ipot)))/(12.d0*drdi(i))
                v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)
              enddo ! i

              i = irws1 - 1
                dv = (-v(i-3,lm,ipot) +6.d0*v(i-2,lm,ipot) -18.d0*v(i-1,lm,ipot) +10.d0*v(i,lm,ipot) +3.d0*v(i+1,lm,ipot))/(12.d0*drdi(irws1-1))
                v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)

              i = irws1
                dv = (3.d0*v(i-4,lm,ipot) -16.d0*v(i-3,lm,ipot) +36.d0*v(i-2,lm,ipot) -48.d0*v(i-1,lm,ipot) +25.d0*v(i,lm,ipot))/(12.d0*drdi(irws1))
                v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)

!debug        tail correction
            tail_cor(m) = tail_cor(m) + fac*rhoc(irws1,ipot)/(4.d0*pi)*v(irws1,lm,ipot)
!debug

          enddo ! ispin
!
!---> integrate with simpson subroutine
!
          vint1 = simpson(v1, 1, irws1, drdi)

          flmxc(m) = -fac*vint1 - flmc(m)
          flm(m) = flm(m) + flmxc(m)

       enddo ! m

!debug tail correction
!      flm = flm + tail_cor
!debug

!     result: total force in flm

!  8999 format(a,1x,i5,a,f5.0,a,3(f8.4))
!  9000 FORMAT (13x,'error stop in subroutine force : the charge density has to contain non spherical contributions up to l=1 at least ')
!  9100 FORMAT (1x,33 ('-'),' force on the nucleus ',33 ('-'),/,34x,' in units Ry/(a(Bohr) ')
!  9600 FORMAT (7x,'fhx=',e12.6,2x,'fcx=',e12.6,2x,'fxcx=',e12.6,2x,'fx=',e12.6,' Ry/(a(Bohr))')
!  9601 FORMAT (7x,'fhy=',e12.6,2x,'fcy=',e12.6,2x,'fxcy=',e12.6,2x,'fy=',e12.6,' Ry/(a(Bohr))')
!  9602 FORMAT (7x,'fhz=',e12.6,2x,'fcz=',e12.6,2x,'fxcz=',e12.6,2x,'fz=',e12.6,' Ry/(a(Bohr))')
  endsubroutine ! forcxc
