subroutine rhoin(ar, cden, cr, df, gmat, ek, rho2ns, irc1, nsra, efac, pz, fz, &
  qz, sz, cleb, icleb, jend, iend, ekl, cdenlm) ! lm-dos
!-----------------------------------------------------------------------
!
!     calculates the charge density inside r(irmin) in case
!      of a non spherical input potential .
!
!     fills the array cden for the complex density of states
!
!      the non spher. wavefunctions are approximated in that region
!       in the following way :
!
!           the regular one (ir < irmin = irws-irns) :
!
!              pns(ir,lm1,lm2) = pz(ir,l1) * ar(lm1,lm2)
!
!          where pz is the regular wavefct of the spherically symmetric
!          part of the potential and ar the alpha matrix .
!          (see subroutine regns)
!
!
!           the irregular one (ir < irmin) :
!
!              qns(ir,lm1,lm2) = pz(ir,l1) * cr(lm1,lm2)
!                                    + qz(ir,l1) * dr(lm1,lm2)
!
!          where pz is the regular and qz is the irregular
!          wavefct of the spherically symmetric part of the
!          potential and cr , dr the matrices calculated
!          at the point irmin .  (see subroutine irwns)
!
!     attention : the gaunt coeffients which are used here
!                 are ordered in a special way !   (see subroutine
!                 gaunt)
!
!                 remember that the matrices ar,cr,dr are rescaled !
!                 (see subroutines irwns and regns)
!
!                 arrays rho2ns and cden are initialize in subroutine
!                 rhoout .
!
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
  implicit none
  include 'inc.p'
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
!
!..
!.. Scalar Arguments ..
!..
!.. Array Arguments ..
  integer :: lmmaxd
  integer :: lmpotd
  parameter (lmmaxd=(krel+1)*(lmaxd+1)**2)
  parameter (lmpotd=(lpotd+1)**2)
! lm-dos
!..
  double complex :: df, ek
  integer :: iend, irc1, nsra
!.. Local Scalars ..
!..
  double complex :: ar(lmmaxd, *), cden(irmd, 0:lmaxd), cr(lmmaxd, *), &
    efac(*), ekl(0:lmaxd), fz(irmd, 0:lmaxd), gmat(lmmaxd, lmmaxd), &
    pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd), &
    cdenlm(irmd, lmmaxd) !.. Local Arrays ..
  double precision :: cleb(*), rho2ns(irmd, lmpotd)
  integer :: icleb(ncleb, 4), jend(lmpotd, 0:lmaxd, 0:lmaxd)
!..
!.. External Functions ..
  double complex :: czero, efac1, efac2, ffz, gmatl, ppz, v1, v2
  double precision :: c0ll
  integer :: i, ir, j, j0, j1, l, l1, l2, lm1, lm2, lm3, lm3max, ln2, ln3, m
!..
!.. Intrinsic Functions ..
  double complex :: vr(lmmaxd, lmmaxd), wf(irmd, 0:lmaxd, 0:lmaxd), &
    wr(lmmaxd, lmmaxd)
!..
!.. Save statement ..
  double complex :: zdotu
  external :: zdotu
!..
!.. Data statements ..
  intrinsic :: atan, dimag, sqrt
!..
!
  save :: czero
!C0LL = 1/sqrt(4*pi)
!
  data czero/(0.0d0, 0.0d0)/

!
!---> set up array wr(lm1,lm2)
  c0ll = 1.0d0/sqrt(16.0d0*atan(1.0d0))
!        use first vr
!
  lm3max = icleb(iend, 3)

!
!---> using symmetry of structural green function
!
  do lm2 = 1, lmmaxd
    ln2 = lm2
    v2 = efac(lm2)*efac(lm2)*gmat(ln2, ln2)
    do lm1 = 1, lmmaxd
      vr(lm1, lm2) = ek*cr(lm1, lm2) + v2*ar(lm1, lm2)
    end do
  end do
!
!
!---> using symmetry of gaunt coeffients
!
  do lm2 = 2, lmmaxd
    ln2 = lm2
    efac2 = 2.0d0*efac(lm2)
    do lm3 = 1, lm2 - 1
      ln3 = lm3
      v1 = efac2*gmat(ln3, ln2)*efac(lm3)
      do lm1 = 1, lmmaxd
        vr(lm1, lm2) = vr(lm1, lm2) + v1*ar(lm1, lm3)
      end do
    end do
  end do
!
  do lm1 = 1, lmmaxd
    efac1 = efac(lm1)
    wr(lm1, lm1) = zdotu(lmmaxd, ar(lm1,1), lmmaxd, vr(lm1,1), lmmaxd)/ &
      (efac1*efac1)
    do lm2 = 1, lm1 - 1
!---> set up array wf(l1,l2) = pz(l1)*pz(l2)
!

      efac2 = efac(lm2)
      wr(lm1, lm2) = (zdotu(lmmaxd,ar(lm1,1),lmmaxd,vr(lm2, &
        1),lmmaxd)+zdotu(lmmaxd,ar(lm2,1),lmmaxd,vr(lm1, &
        1),lmmaxd))/(efac1*efac2)
    end do
  end do
!
!---> first calculate only the spherically symmetric contribution
!     remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
  if (nsra==2) then
    do l1 = 0, lmaxd
      do l2 = 0, l1
        do ir = 2, irc1
          wf(ir, l1, l2) = pz(ir, l1)*pz(ir, l2) + fz(ir, l1)*fz(ir, l2)
        end do
      end do
    end do
!
  else
    do l1 = 0, lmaxd
      do l2 = 0, l1
        do ir = 2, irc1
          wf(ir, l1, l2) = pz(ir, l1)*pz(ir, l2)
        end do
      end do
    end do
  end if
!
! Implicit integration over energies


  do l = 0, lmaxd
    gmatl = czero
    do m = -l, l
      lm1 = l*(l+1) + m + 1
      gmatl = gmatl + wr(lm1, lm1)
    end do
!lm-dos
    if (nsra==2) then
      do i = 2, irc1
        ppz = pz(i, l)
        ffz = fz(i, l)
        cden(i, l) = ppz*(gmatl*ppz+ekl(l)*qz(i,l)) + &
          ffz*(gmatl*ffz+ekl(l)*sz(i,l))
        rho2ns(i, 1) = rho2ns(i, 1) + c0ll*dimag(df*cden(i,l)) !lm-dos    
!lm-dos
!lm-dos
        do m = -l, l !lm-dos
          lm1 = l*(l+1) + m + 1 
          cdenlm(i, lm1) = ppz*(wr(lm1,lm1)*ppz+ek*qz(i,l)) + & 
            ffz*(wr(lm1,lm1)*ffz+ek*sz(i,l)) ! Implicit integration over energies
        end do 
!lm-dos
      end do
!lm-dos    
    else
      do i = 2, irc1
        ppz = pz(i, l)
        cden(i, l) = ppz*(gmatl*ppz+ekl(l)*qz(i,l))
        rho2ns(i, 1) = rho2ns(i, 1) + c0ll*dimag(df*cden(i,l)) !lm-dos
!lm-dos
        do m = -l, l 
          lm1 = l*(l+1) + m + 1 
          cdenlm(i, lm1) = ppz*(wr(lm1,lm1)*ppz+ek*qz(i,l)) !
        end do !---> calculate the non spherically symmetric contribution
!        to speed up the pointer jend generated in gaunt is used
      end do
    end if
!        remember that the wavefunctions are l and not lm dependent
  end do
!
!
!
!
!
  j0 = 1
!
  do lm3 = 2, lm3max
    do l1 = 0, lmaxd
      do l2 = 0, l1
!---> sum over m1,m2 for fixed lm3,l1,l2
        j1 = jend(lm3, l1, l2)
!
        if (j1/=0) then

          gmatl = czero
!
!

          do j = j0, j1
            lm1 = icleb(j, 1)
            lm2 = icleb(j, 2)
            gmatl = gmatl + cleb(j)*wr(lm1, lm2)
          end do


          j0 = j1 + 1

          gmatl = df*gmatl
          do i = 2, irc1
            rho2ns(i, lm3) = rho2ns(i, lm3) + dimag(gmatl*wf(i,l1,l2))
          end do
        end if
! lm-dos
      end do
!-----------------------------------------------------------------------
    end do
!
  end do
!     calculates the charge density inside r(irmin) in case
end subroutine
