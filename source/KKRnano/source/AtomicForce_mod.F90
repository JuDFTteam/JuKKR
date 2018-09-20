!-------------------------------------------------------------------------------
!> Summary: Forces onto the atomic positions
!> Author: Rudolf Zeller, Paul F Baumeister, et al.
!> Category: KKRnano, forces, physical-observables, potential
!-------------------------------------------------------------------------------
module AtomicForce_mod
implicit none
  private
  
  public :: force, forceh, forcxc

  contains

  subroutine force(flm, flmc, lpot, nspin, rhoc, v, r, drdi, irws, irmd)
  use Quadrature_mod, only: Simpson
  use Constants_mod, only: pi
!-----------------------------------------------------------------------
!     calculates the force on nucleus m
!     from a given non spherical charge density at the nucleus site r
!     with core correction (Coulomb contribution)
!-----------------------------------------------------------------------
    integer, intent(in) :: lpot, nspin
    integer, intent(in) :: irws, irmd
    double precision, intent(in) :: drdi(irmd), r(irmd), rhoc(irmd,*), v(irmd,(lpot+1)**2,2)
    double precision, intent(out) :: flm(-1:1), flmc(-1:1)

    !     .. local vars ..
    integer :: i, ipot, ispin, lm, m
    double precision :: dv, fac, rws, vint1
    double precision :: flmh(-1:1), v1(irmd)

    fac = dsqrt((4.d0*pi)/3.d0)
    if (lpot < 1) then
      write(*, *) "error stop in subroutine force : the charge density has to contain non spherical contributions up to l=1 at least"
      flm  = 0.d0
      flmc = 0.d0
      return
    endif

    rws = r(irws)

    do m = -1, 1
      lm = 2 + m + 1  ! only lm = 2,3,4 needed (l = 1, m =-1,0,1)
!
!---> initialize v1
      v1(1:irws) = 0.d0
!
      do ispin = 1, nspin
!
!---> determine the right potential numbers
!
          ipot = ispin
!
!---> determine the derivative of the potential using a 5-point formula
!
        i = 2
          dv = (-3.d0*v(i-1,lm,ipot) -10.d0*v(i,lm,ipot) +18.d0*v(i+1,lm,ipot) -6.d0*v(i+2,lm,ipot) +1.d0*v(i+3,lm,ipot))/(12.d0*drdi(2))
          v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)

        do i = 3, irws - 2
          dv = (v(i-2,lm,ipot) - v(i+2,lm,ipot) + 8.d0*(v(i+1,lm,ipot) - v(i-1,lm,ipot)))/(12.d0*drdi(i))
          v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)
        enddo ! i

        i = irws - 1
          dv = (-v(i-3,lm,ipot) +6.d0*v(i-2,lm,ipot) -18.d0*v(i-1,lm,ipot) +10.d0*v(i,lm,ipot) +3.d0*v(i+1,lm,ipot))/(12.d0*drdi(i))
          v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)

        i = irws
          dv = (3.d0*v(i-4,lm,ipot) -16.d0*v(i-3,lm,ipot) +36.d0*v(i-2,lm,ipot) -48.d0*v(i-1,lm,ipot) +25.d0*v(i,lm,ipot))/(12.d0*drdi(i))
          v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)

      enddo ! ispin
!
!---> integrate with Simpson subroutine
!
      vint1 = Simpson(v1, 1, irws, drdi)

      flmh(m) = fac*flm(m)
      flmc(m) = -fac*vint1
      flm(m) = flmh(m) + flmc(m)
    enddo ! m

  endsubroutine ! force

!     calculates the Hellmann-Feynman force.
!
!     attention: a factor sqrt(4*pi/3) is missing here, this is corrected in routine 'force'
  subroutine forceh(flmh, lpot, rho2ns, v, r, drdi, irws, z, irmd)
  use Quadrature_mod, only: Simpson
  use Constants_mod, only: pi
!-----------------------------------------------------------------------
!     calculates the force on nucleus m with Hellmann - Feynman theorem
!     from a given non spherical charge density at the nucleus site r
!-----------------------------------------------------------------------
    integer, intent(in) :: irmd, lpot, irws
    double precision, intent(in) :: drdi(irmd), r(irmd), rho2ns(irmd,(lpot+1)**2), v(irmd,(lpot+1)**2,2), z
    double precision, intent(out) :: flmh(-1:1)
    
    !  .. local vars ..
    double precision :: rws, vint1
    integer :: i, ipot, lm, m
    double precision :: flm(-1:1,2), v1(irmd)
    
    if (lpot < 1) then
      write(*, *) "error stop in subroutine force : the charge density has to contain non spherical contributions up to l=1 at least"
      flmh = 0.d0
      return
    endif
!
!
!---> reading the right wigner-s. radius
!
    rws = r(irws)
!
!---> determine the right potential numbers
!
    ipot = 1

    do m = -1, 1
      lm = 2 + m + 1
!
      v1(1) = 0.d0
      do i = 2, irws
          v1(i) = rho2ns(i,lm)*(r(i)**(-2.d0) - r(i)/rws**3)
      enddo ! i
!
!---> integrate with Simpson subroutine
!
      vint1 = Simpson(v1, 1, irws, drdi)
!
      flm(m,1) = 2.d0*vint1
!
!---> use Coulomb potential to determine extra atomic contribution
!
      flm(m,2) = v(irws,lm,ipot)*(3.d0/(4.d0*pi*rws))
!
!---> total Hellmann-Feynman force
!
      flmh(m) = (flm(m,1) + flm(m,2))*z
    enddo ! m

  endsubroutine ! forceh

  
  subroutine forcxc(flm, flmc, lpot, nspin, rhoc,v,r, drdi, irws, irmd)
    use Quadrature_mod, only: Simpson
    use Constants_mod, only: pi
!-----------------------------------------------------------------------
!     calculates the force on nucleus m
!     from a given non spherical charge density at the nucleus site r
!     with core correction(exchange contribution)
!-----------------------------------------------------------------------
!   integer, parameter :: lmpotd = (lpotd+1)**2
    integer, intent(in) :: lpot, nspin, irmd, irws
    double precision, intent(in) :: drdi(irmd), flmc(-1:1), r(irmd), rhoc(irmd,*), v(irmd,(lpot+1)**2,2)
    double precision, intent(inout) :: flm(-1:1)

!     .. locals ..
    double precision :: dv, fac, rws, vint1
    integer :: i, ipot, ispin, lm, m
    double precision :: flmxc(-1:1), v1(irmd), tail_cor(-1:1)
    
    fac = sqrt((4.d0*pi)/3.d0)

    if (lpot < 1) then
      write(*, *) "error stop in subroutine force : the charge density has to contain non spherical contributions up to l=1 at least"
      return
    endif

    tail_cor(:) = 0.d0

    rws = r(irws)

    do m = -1, 1
      lm = 2 + m + 1

      v1(1:irws) = 0.d0

      do ispin = 1, nspin
          ipot = ispin

          i = 2
            dv = (-3.d0*v(i-1,lm,ipot) -10.d0*v(i,lm,ipot) +18.d0*v(i+1,lm,ipot) -6.d0*v(i+2,lm,ipot) + v(i+3,lm,ipot))/(12.d0*drdi(2))
            v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)
            
          do i = 3, irws - 2
            dv = (v(i-2,lm,ipot) - v(i+2,lm,ipot) + 8.d0*(v(i+1,lm,ipot) - v(i-1,lm,ipot)))/(12.d0*drdi(i))
            v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)
          enddo ! i

          i = irws - 1
            dv = (-v(i-3,lm,ipot) +6.d0*v(i-2,lm,ipot) -18.d0*v(i-1,lm,ipot) +10.d0*v(i,lm,ipot) +3.d0*v(i+1,lm,ipot))/(12.d0*drdi(irws-1))
            v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)

          i = irws
            dv = (3.d0*v(i-4,lm,ipot) -16.d0*v(i-3,lm,ipot) +36.d0*v(i-2,lm,ipot) -48.d0*v(i-1,lm,ipot) +25.d0*v(i,lm,ipot))/(12.d0*drdi(irws))
            v1(i) = rhoc(i,ipot)*(2.d0*v(i,lm,ipot)/r(i) + dv)/(4.d0*pi) + v1(i)

!debug        tail correction
        tail_cor(m) = tail_cor(m) + fac*rhoc(irws,ipot)/(4.d0*pi)*v(irws,lm,ipot)
!debug

      enddo ! ispin
!
!---> integrate with Simpson subroutine
!
      vint1 = Simpson(v1, 1, irws, drdi)

      flmxc(m) = -fac*vint1 - flmc(m)
      flm(m) = flm(m) + flmxc(m)

    enddo ! m

!debug tail correction
!      flm = flm + tail_cor
!debug

!     result: total force in flm

!  8999 format(a,1x,i5,a,f5.0,a,3(f8.4))
!  9000 format (13x,'error stop in subroutine force : the charge density has to contain non spherical contributions up to l=1 at least ')
!  9100 format (1x,33 ('-'),' force on the nucleus ',33 ('-'),/,34x,' in units ry/(a(bohr) ')
!  9600 format (7x,'fhx=',e12.6,2x,'fcx=',e12.6,2x,'fxcx=',e12.6,2x,'fx=',e12.6,' ry/(a(bohr))')
!  9601 format (7x,'fhy=',e12.6,2x,'fcy=',e12.6,2x,'fxcy=',e12.6,2x,'fy=',e12.6,' ry/(a(bohr))')
!  9602 format (7x,'fhz=',e12.6,2x,'fcz=',e12.6,2x,'fxcz=',e12.6,2x,'fz=',e12.6,' ry/(a(bohr))')
  endsubroutine ! forcxc

endmodule ! AtomicForce_mod
