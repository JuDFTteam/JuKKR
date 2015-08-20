#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

subroutine regsol(cvlight,e,nsra,dlogdp,fz,hamf,mass,pz,dror,r, &
s,vm2z,z,ipan,ircut,irmd,ipand,lmaxd, &
ldau,nldau,lldau,wmldauav,ldaucut)
  !
  implicit none
  !-----------------------------------------------------------------------
  !  calculates the regular solution of the schroedinger equation or
  !    in semi relativistic approximation for a spherically averaged
  !    potential and given energy . to archieve greater presion the
  !    leading power r**s ( in schroedinger case s = l , in case of sra
  !    s = sqrt( (l*l+l-1) - 4*z*z/c/c ) ) is analytically separated
  !    from the wavefunction .
  !
  !  the t - matrix has to be determined at the mt radius in case of
  !    a mt calculation or at the ws radius in case of a ws calcu-
  !    lation . therefore the logarithmic derivative is calculated
  !    at that point (=ircut(ipan) )
  !
  !  the differential equation is solved with a 5 point adams - bashforth
  !    and adams - moulton predictor corrector method integrating
  !    outwards and extended for potentials with kinks
  !
  !                                               b.drittler   nov 1989
  !-----------------------------------------------------------------------

  ! Arguments
  double precision :: cvlight
  double complex :: e
  integer :: nsra
  double complex, dimension(0:lmaxd) :: dlogdp
  double complex, dimension(irmd,0:lmaxd) :: fz
  double complex, dimension(irmd,0:lmaxd) :: hamf
  double complex, dimension(irmd) :: mass
  double complex, dimension(irmd,0:lmaxd) :: pz
  double precision, dimension(irmd) :: dror
  double precision, dimension(irmd) :: r
  double precision, dimension(0:lmaxd) :: s
  double precision, dimension(irmd) :: vm2z
  double precision :: z
  integer :: ipan
  integer, dimension(0:ipand) :: ircut
  integer :: irmd
  integer :: ipand
  integer :: lmaxd
  logical :: ldau
  integer :: nldau
  integer, dimension(lmaxd+1) :: lldau
  double precision, dimension(lmaxd+1) :: wmldauav
  double precision, dimension(irmd) :: ldaucut

  ! Local variables of regsol
  double complex :: dfd0
  double complex :: dpd0
  double complex :: fip0
  double complex :: fip1
  double complex :: hamf1
  double complex :: k1f
  double complex :: k1p
  double complex :: k2f
  double complex :: k2p
  double complex :: k3f
  double complex :: k3p
  double complex :: k4f
  double complex :: k4p
  double complex :: mass1
  double complex :: pip0
  double complex :: pip1
  double complex :: vme
  double complex :: vmefac
  double complex :: vmetr1
  double precision :: dror1
  double precision :: drsm1
  double precision :: drsp1
  double precision :: s1
  double precision :: sm1
  double precision :: sp1
  double precision :: srafac
  double precision :: zocsq
  integer :: ip
  integer :: ir
  integer :: irc
  integer :: ire
  integer :: irs
  integer :: irsp1
  integer :: j
  integer :: k
  integer :: L
  integer :: ildau
  double complex, dimension(-1:4) :: a
  double complex, dimension(0:4) :: b
  double complex, dimension(-4:0) :: dfdi
  double complex, dimension(-4:0) :: dpdi
  double complex, dimension(irmd) :: hamfldau

  !     .. Intrinsic Functions ..
  intrinsic cmplx,dble
  !     ..
  CHECKASSERT(irmd > 6)

  if (nsra == 2) then
    !
    !---> in case of sra  srafac = 1/c - otherwise srafac = 0
    !
    srafac = 1.d0/cvlight

  else
    srafac = 0.d0
  end if
  !
  irc = ircut(ipan)
  !
  do ir = 2, irc
    vmetr1 = (vm2z(ir)-e)*r(ir) - 2.d0*z
    hamf(ir,0) = vmetr1*dror(ir)
    mass(ir) = r(ir) - srafac*srafac*vmetr1
  enddo ! ir
   !
   do L = 1, lmaxd
     do ir = 7, irc
       hamf(ir,L) = dble(L*L+L)/mass(ir)*dror(ir) + hamf(ir,0)
     enddo ! ir
   enddo ! L
   !
   !-----------------------------------------------------------------------
   ! LDA+U
   !
   !  Account for potential shift in case of LDA+U (averaged over m)
   !  by adding the average WLDAUAV to the spherical part of the
   !  potential.
   !
   if (ldau) then
     do ildau = 1, nldau
       if (lldau(ildau) >= 0) then
         s1 = dble(lldau(ildau)*lldau(ildau) + lldau(ildau))
         do ir = 2,irc
           vmetr1 = (vm2z(ir) - e + wmldauav(ildau)*ldaucut(ir))*r(ir)
           hamfldau(ir) = (vmetr1 - 2.d0*z)*dror(ir)
         enddo ! ir
         do ir = 7, irc
           hamf(ir,lldau(ildau)) = s1/mass(ir)*dror(ir) + hamfldau(ir)
         enddo ! ir
       endif
     enddo
   endif
   !
   ! LDA+U
   !-----------------------------------------------------------------------
   !

   do ir = 2, irc
     mass(ir) = mass(ir)*dror(ir)
   enddo ! ir
   
   do L = 0, lmaxd
     !
     s1 = s(L)
     sm1 = s1 - 1.d0
     sp1 = s1 + 1.d0
     !
     !---> loop over the number of kinks
     !
     do ip = 1, ipan
       !
       if (ip == 1) then
         irs = 2
         ire = ircut(1)

         !
         !---> initial values
         !
         vme = vm2z(2) - e
         vmefac = 1.d0 - vme*srafac*srafac
         if (nsra == 2 .and. z > 0.d0) then

           zocsq = -2.d0*z*z/(cvlight*cvlight)
           a(-1) = 0.d0
           a(0) = 1.d0
           b(0) = cmplx(sm1*cvlight*cvlight/(2*z), 0.d0)
           do j = 1,3
             a(j) = (0.d0,0.d0)
             b(j) = (0.d0,0.d0)
           enddo ! j

         else

           a(0) = 0.d0
           b(0) = dble(L)/vmefac
           a(1) = 1.d0
           do j = 2, 4
             a(j) = (vme*vmefac*a(j-2) - 2.d0*z*a(j-1))/dble((j-1)*(j+2*L))
             b(j-1) = dble(L+j-1)*a(j)/vmefac
           enddo ! j

         endif
         !
         k = -4
         !
         !---> power series near origin
         !
         do ir = 2, 6
           pip0 = a(3)
           dpd0 = 3.d0*a(3)
           fip0 = b(3)
           dfd0 = 3.d0*b(3)
           do j = 2, 0, -1
             pip0 = a(j) + pip0*r(ir)
             dpd0 = dble(j)*a(j) + dpd0*r(ir)
             fip0 = b(j) + fip0*r(ir)
             dfd0 = dble(j)*b(j) + dfd0*r(ir)
           enddo ! j
           !
           pz(ir,L) = pip0
           fz(ir,L) = fip0
           dpdi(k) = dpd0*dror(ir)
           dfdi(k) = dfd0*dror(ir)
           !
           k = k + 1
         enddo ! ir

       else ! panels 2, 3, ...
         !
         !---> runge kutta step to restart algorithm
         !
         irs = ircut(ip-1) + 1
         ire = ircut(ip)
         CHECKASSERT((ire - irs) >= 4) ! minimum of 5 points per panel!
         irsp1 = irs + 1
         pip0 = pz(irs,L)
         fip0 = fz(irs,L)
         drsp1 = dror(irs)*sp1
         drsm1 = dror(irs)*sm1
         dpdi(-4) = mass(irs)*  fip0 - drsm1*pip0
         dfdi(-4) = hamf(irs,L)*pip0 - drsp1*fip0
         !
         !---> first step - 4 point runge kutta with interpolation
         !
         k1p = dpdi(-4)
         k1f = dfdi(-4)
         !
         dror1 = (3.d0*dror(irs+3)   - 15.d0*dror(irs+2)   + 45.d0*dror(irsp1)   + 15.d0*dror(irs)  )/48.d0
         mass1 = (3.d0*mass(irs+3)   - 15.d0*mass(irs+2)   + 45.d0*mass(irsp1)   + 15.d0*mass(irs)  )/48.d0
         hamf1 = (3.d0*hamf(irs+3,L) - 15.d0*hamf(irs+2,L) + 45.d0*hamf(irsp1,L) + 15.d0*hamf(irs,L))/48.d0
         drsp1 = dror1*sp1
         drsm1 = dror1*sm1
         k2p = mass1*(fip0 + 0.5d0*k1f) - drsm1*(pip0 + 0.5d0*k1p)
         k2f = hamf1*(pip0 + 0.5d0*k1p) - drsp1*(fip0 + 0.5d0*k1f)
         k3p = mass1*(fip0 + 0.5d0*k2f) - drsm1*(pip0 + 0.5d0*k2p)
         k3f = hamf1*(pip0 + 0.5d0*k2p) - drsp1*(fip0 + 0.5d0*k2f)
         !
         drsp1 = dror(irsp1)*sp1
         drsm1 = dror(irsp1)*sm1
         k4p = mass(irsp1)  *(fip0 + k3f) - drsm1*(pip0 + k3p)
         k4f = hamf(irsp1,L)*(pip0 + k3p) - drsp1*(fip0 + k3f)
         pip0 = pip0 + (k1p + 2.d0*(k2p + k3p) + k4p)/6.d0
         fip0 = fip0 + (k1f + 2.d0*(k2f + k3f) + k4f)/6.d0
         !
         pz(irsp1,L) = pip0
         fz(irsp1,L) = fip0
         dpdi(-3) = mass(irsp1)*  fip0 - drsm1*pip0
         dfdi(-3) = hamf(irsp1,L)*pip0 - drsp1*fip0
         !
         k = -2
         !
         !---> 4 point runge kutta with h = i+2 - i
         !
         do ir = irs + 2, irs + 4
           pip0 = pz(ir-2,L)
           fip0 = fz(ir-2,L)
           k1p = dpdi(k-2)
           k1f = dfdi(k-2)
           k2p = mass(ir-1)* (fip0+k1f) - drsm1* (pip0+k1p)
           k2f = hamf(ir-1,L)* (pip0+k1p) - drsp1* (fip0+k1f)
           k3p = mass(ir-1)* (fip0+k2f) - drsm1* (pip0+k2p)
           k3f = hamf(ir-1,L)* (pip0+k2p) - drsp1* (fip0+k2f)
           !
           drsp1 = dror(ir)*sp1
           drsm1 = dror(ir)*sm1
           !
           k4p = mass(ir)*  (fip0 + 2.d0*k3f) - drsm1*(pip0 + 2.d0*k3p)
           k4f = hamf(ir,L)*(pip0 + 2.d0*k3p) - drsp1*(fip0 + 2.d0*k3f)
           pip0 = pip0 + (k1p + 2.d0*(k2p + k3p) + k4p)/3.d0
           fip0 = fip0 + (k1f + 2.d0*(k2f + k3f) + k4f)/3.d0
           !
           pz(ir,L) = pip0
           fz(ir,L) = fip0
           dpdi(k) = mass(ir)*  fip0 - drsm1*pip0
           dfdi(k) = hamf(ir,L)*pip0 - drsp1*fip0
           k = k + 1
         enddo ! ir
       endif

       ! if >5 points in panel, calculate those points with 5-point formula
       do ir = irs + 5, ire
         drsp1 = dror(ir)*sp1
         drsm1 = dror(ir)*sm1
         !
         !---> predictor : 5 point adams - bashforth
         !
         pip1 = pip0 + (1901.d0*dpdi(0) - 2774.d0*dpdi(-1) + 2616.d0*dpdi(-2) - 1274.d0*dpdi(-3) + 251.d0*dpdi(-4))/720.d0
         fip1 = fip0 + (1901.d0*dfdi(0) - 2774.d0*dfdi(-1) + 2616.d0*dfdi(-2) - 1274.d0*dfdi(-3) + 251.d0*dfdi(-4))/720.d0
         !
         dpdi(-4) = dpdi(-3)
         dpdi(-3) = dpdi(-2)
         dpdi(-2) = dpdi(-1)
         dpdi(-1) = dpdi( 0)
         dfdi(-4) = dfdi(-3)
         dfdi(-3) = dfdi(-2)
         dfdi(-2) = dfdi(-1)
         dfdi(-1) = dfdi( 0)
         !
         dpdi( 0) = mass(ir)*  fip1 - drsm1*pip1
         dfdi( 0) = hamf(ir,L)*pip1 - drsp1*fip1
         !
         !---> corrector : 5 point adams - moulton
         !
         pip0 = pip0 + (251.d0*dpdi(0) + 646.d0*dpdi(-1) - 264.d0*dpdi(-2) + 106.d0*dpdi(-3) - 19.d0*dpdi(-4))/720.d0
         fip0 = fip0 + (251.d0*dfdi(0) + 646.d0*dfdi(-1) - 264.d0*dfdi(-2) + 106.d0*dfdi(-3) - 19.d0*dfdi(-4))/720.d0
         !
         pz(ir,L) = pip0
         fz(ir,L) = fip0
         dpdi(0) = mass(ir)*  fip0 - drsm1*pip0
         dfdi(0) = hamf(ir,L)*pip0 - drsp1*fip0
       enddo ! ir
       !
       !---> remember that the r - mesh contains the kinks two times
       !     store the values of pz and fz to restart the algorithm
       !
       if (ip /= ipan) then
         pz(ire+1,L) = pip0
         fz(ire+1,L) = fip0
       end if

    enddo ! ip ! end loop over panels

     !
     !---> logarithmic derivate of real wavefunction ( r**s *pz / r)
     !
     dlogdp(L) = (dpdi(0)/(pip0*dror(irc)) + sm1)/r(irc)
  enddo ! L ! end loop over L

end
