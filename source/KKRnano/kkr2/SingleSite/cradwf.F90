      subroutine cradwf(e,ek,nsra,alpha,ipan,ircut,cvlight,rs,s,pz,fz,qz,sz,tmat,vm2z,drdi,r,z,ldau,nldau,lldau,wmldauav,ldaucut, lmaxd, irmd, ipand)
      implicit none
!-----------------------------------------------------------------------
!  subroutine for radial wave functions of spherical potentials
!
!             the generalized phase shifts are calculated by
!             a wronski relation :
!
!                 alpha(z,l) =-sqrt(z)*wronski{hl(r;z),rl(r;z)}; r->0
!
!             where hl is the free hankel function and rl the regular
!             solution . using the analytical behaviour of rl at the
!             origin (rl = alphal * r**(l+1)  ; r->0),
!             the generalized phase shifts can be calculated
!             directly with the renormalization alphal .
!                                           b.drittler nov.1987
!-----------------------------------------------------------------------

      integer, intent(in) :: lmaxd, irmd, ipand
      double complex, intent(in) :: e, ek
      double precision, intent(in) :: cvlight, z
      integer, intent(in) :: ipan,nsra,nldau
      logical, intent(in) :: ldau
      double complex alpha(0:lmaxd),fz(irmd,0:lmaxd),pz(irmd,0:lmaxd), qz(irmd,0:lmaxd),sz(irmd,0:lmaxd),tmat(0:lmaxd)
      double precision   drdi(irmd),r(irmd), rs(irmd,0:lmaxd),s(0:lmaxd), vm2z(irmd), ldaucut(irmd), wmldauav(lmaxd+1)
      integer            ircut(0:ipand), lldau(lmaxd+1)
      
      
      double complex, parameter :: ci=(0.d0,1.d0), zero=(0.0d0,0.0d0)
      double complex :: alphal,arg,bl,eklfac,hl,pn,qf,slope,tl,tlsqez, value,w,x,y
      double precision :: rirc,rsirc,s1
      integer :: ir,irc1,l,n, lmaxp1
      double complex :: bessjw(0:lmaxd+1),bessyw(0:lmaxd+1),dlogdp(0:lmaxd)
      double complex :: hamf(irmd,0:lmaxd),hankws(0:lmaxd+1),mass(irmd)
      double precision :: dror(irmd)

      lmaxp1 = lmaxd+1

! initialisations
      dror = 0.d0
      hamf = zero
      mass = zero
! initialisations

      irc1 = ircut(ipan)
      do ir = 2,irc1
        dror(ir) = drdi(ir)/r(ir)
      end do
      rirc = r(irc1)
      arg = rirc*ek
      call beshan(hankws,bessjw,bessyw,arg,lmaxp1)
!
!    attention : contrary to abramowitz and stegun and
!                the bessel functions of third kind ( hankel functions)
!                are defined as:      hl(l) = nl(l) - i * jl(l)
!
      do l = 0,lmaxp1
        hankws(l) = bessyw(l) - ci*bessjw(l)
      end do
!
!---> calculate regular wavefunctions
!
      call regsol(cvlight,e,nsra,dlogdp,fz,hamf,mass,pz,dror,r,s,vm2z, &
                    z,ipan,ircut,irmd,ipand,lmaxd, &
                    ldau,nldau,lldau,wmldauav,ldaucut)
!
      eklfac = ek
!
      do l = 0, lmaxd
        s1 = s(l)
        rsirc = rs(irc1,l)
        eklfac = eklfac/ek*dble(2*l+1)
!
!---> determine t - matrix
!
        pn = pz(irc1,l)*rsirc
        n = l+1
        qf = dble(l)/rirc
        hl = hankws(l)
        bl = bessjw(l)
        x = qf*hl - ek*hankws(n)
        y = qf*bl - ek*bessjw(n)
        w = dlogdp(l)
        tlsqez = (bl*w - y)/(x - hl*w)
        tl = tlsqez/ek
        tmat(l) = tl
!
!---> determine the renormalization
!
        alphal = (bl + hl*tlsqez)*rirc/pn
!
!---> determine the alpha matrix
!
        alpha(l) = alphal*eklfac

        pz(2:irc1,l) = pz(2:irc1,l)*alphal
        fz(2:irc1,l) = fz(2:irc1,l)*alphal
!
        value = hl*rirc*rsirc
        slope = dble(l+1)*hl - rirc*ek*hankws(l+1)
        slope = (slope*rsirc + s1/rirc*value)
        qz(irc1,l) = value
        sz(irc1,l) = (slope*rirc - (s1 + 1.0d0)*value)/mass(irc1)*dror(irc1)
     enddo ! l
!
!---> calculate irregular wavefunctions
!
      call irwsol(ek,fz,hamf,mass,pz,qz,sz,dror,s,ipan,ircut,irmd, ipand,lmaxd)

      do l = 0, lmaxd
        if (nsra == 2) then
          pz(2:irc1,l) = pz(2:irc1,l)*rs(2:irc1,l)
          qz(2:irc1,l) = qz(2:irc1,l)/rs(2:irc1,l)
          fz(2:irc1,l) = fz(2:irc1,l)*rs(2:irc1,l)/cvlight
          sz(2:irc1,l) = sz(2:irc1,l)/rs(2:irc1,l)/cvlight
        else
          pz(2:irc1,l) = pz(2:irc1,l)*rs(2:irc1,l)
          qz(2:irc1,l) = qz(2:irc1,l)/rs(2:irc1,l)
          fz(2:irc1,l) = zero
          sz(2:irc1,l) = zero
        endif ! nsra
      enddo ! l 

      end
