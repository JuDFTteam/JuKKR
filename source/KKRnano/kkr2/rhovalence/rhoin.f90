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
  !       summands are multiplied by a factor of 2 in the case of lm1  /= 
  !       lm2 or lm01  /=  lm02 .
  !
  !             (see notes by b.drittler)
  !
  !                               b.drittler   aug. 1988
  !-----------------------------------------------------------------------
subroutine rhoin(ar, cden, cr, df, gmat, ek, rho2ns, irc1, nsra, efac, pz,  &
                fz, qz, sz, cleb, icleb, jend, iend, ekl, &
                lmax, irmd, ncleb)
  implicit none
  integer, intent(in) :: lmax, irmd, ncleb
  double complex, intent(in) :: df, ek
  integer, intent(in) :: iend, irc1, nsra
  double complex, intent(in) :: ar((lmax+1)**2,*)
  double complex, intent(out) :: cden(irmd,0:lmax)
  double complex, intent(in) :: cr((lmax+1)**2,*)
  double complex, intent(in) :: efac(*)
  double complex, intent(in) :: ekl(0:lmax)
  double complex, intent(in) :: fz(irmd,0:lmax)
  double complex, intent(in) :: gmat((lmax+1)**2,(lmax+1)**2)
  double complex, intent(in) :: pz(irmd,0:lmax)
  double complex, intent(in) :: qz(irmd,0:lmax)
  double complex, intent(in) :: sz(irmd,0:lmax)
  double precision, intent(in) :: cleb(*)
  double precision, intent(inout) :: rho2ns(irmd,(2*lmax+1)**2)
  integer, intent(in) :: icleb(ncleb,3)
  integer, intent(in) :: jend((2*lmax+1)**2,0:lmax,0:lmax)

  double complex, external :: zdotu
  double complex, parameter :: zero = (0.d0, 0.d0)
  double complex :: efac1,efac2,ffz,gmatl,ppz,v1,v2
  double precision :: c0ll
  integer :: i,ir,j,j0,j1,l,l1,l2,lm1,lm2,lm3,lm3max,ln2,ln3,m
  double complex :: vr((lmax+1)**2, (lmax+1)**2)
  double complex :: wf(irmd,0:lmax,0:lmax)
  double complex :: wr((lmax+1)**2, (lmax+1)**2)
  integer :: lmmaxd, lmaxd
  
  lmaxd = lmax
  lmmaxd= (lmaxd+1)**2

  c0ll = 1.d0/sqrt(16.d0*atan(1.d0))
  !

  lm3max = icleb(iend,3)
  !
  !---> set up array wr(lm1,lm2)
  !        use first vr
  !
  do 20 lm2 = 1,lmmaxd
    ln2 = lm2
    v2 = efac(lm2)*efac(lm2)*gmat(ln2,ln2)
    do 10 lm1 = 1,lmmaxd
      vr(lm1,lm2) = ek*cr(lm1,lm2) + v2*ar(lm1,lm2)
10  continue
20 continue

   !
   !---> using symmetry of structural green function
   !
   do 50 lm2 = 2,lmmaxd
     ln2 = lm2
     efac2 = 2.d0*efac(lm2)
     do 40 lm3 = 1,lm2 - 1
       ln3 = lm3
       v1 = efac2*gmat(ln3,ln2)*efac(lm3)
       do 30 lm1 = 1,lmmaxd
         vr(lm1,lm2) = vr(lm1,lm2) + v1*ar(lm1,lm3)
30     continue
40   continue
50 continue
   !

   do 70 lm1 = 1,lmmaxd
     efac1 = efac(lm1)
     wr(lm1,lm1) = zdotu(lmmaxd,ar(lm1,1),lmmaxd,vr(lm1,1),lmmaxd)/(efac1*efac1)
     do 60 lm2 = 1,lm1 - 1
       !
       !---> using symmetry of gaunt coeffients
       !
       efac2 = efac(lm2)
       wr(lm1,lm2) = (zdotu(lmmaxd,ar(lm1,1),lmmaxd,vr(lm2,1), lmmaxd) + zdotu(lmmaxd,ar(lm2,1),lmmaxd,vr(lm1,1),lmmaxd))/(efac1*efac2)
60   continue
70 continue

   !
   !---> set up array wf(l1,l2) = pz(l1)*pz(l2)
   !
   if (nsra == 2) then
     do 100 l1 = 0,lmaxd
       do 90 l2 = 0,l1
         do 80 ir = 2,irc1
           wf(ir,l1,l2) = pz(ir,l1)*pz(ir,l2) + fz(ir,l1)*fz(ir,l2)
80       continue
90     continue
100  continue

   else
     do 130 l1 = 0,lmaxd
       do 120 l2 = 0,l1
         do 110 ir = 2,irc1
           wf(ir,l1,l2) = pz(ir,l1)*pz(ir,l2)
110      continue
120    continue
130  continue
   end if
   !
   !---> first calculate only the spherically symmetric contribution
   !     remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
   !

   do 170 l = 0,lmaxd
     gmatl = zero
     do 140 m = -l,l
       lm1 = l* (l+1) + m + 1
       gmatl = gmatl + wr(lm1,lm1)
140  continue
     !
     if (nsra == 2) then
       do 150 i = 2,irc1
         ppz = pz(i,l)
         ffz = fz(i,l)
         cden(i,l) = ppz*(gmatl*ppz+ekl(l)*qz(i,l)) + ffz*(gmatl*ffz+ekl(l)*sz(i,l))
         rho2ns(i,1) = rho2ns(i,1) + c0ll*aimag(df*cden(i,l))
150    continue

     else
       do 160 i = 2,irc1
         ppz = pz(i,l)
         cden(i,l) = ppz*(gmatl*ppz+ekl(l)*qz(i,l))
         rho2ns(i,1) = rho2ns(i,1) + c0ll*aimag(df*cden(i,l))
160    continue
     end if

170 continue

    !
    !---> calculate the non spherically symmetric contribution
    !        to speed up the pointer jend generated in gaunt is used
    !        remember that the wavefunctions are l and not lm dependent
    !
    j0 = 1
    !
    do 220 lm3 = 2,lm3max
      do 210 l1 = 0,lmaxd
        do 200 l2 = 0,l1

          !
          j1 = jend(lm3,l1,l2)
          !
          if (j1 /= 0) then
            !
            gmatl = zero
            !
            !---> sum over m1,m2 for fixed lm3,l1,l2
            !
            do 180 j = j0,j1
              lm1 = icleb(j,1)
              lm2 = icleb(j,2)
              gmatl = gmatl + cleb(j)*wr(lm1,lm2)
180         continue

            !
            j0 = j1 + 1
            !
            gmatl = df*gmatl
            do 190 i = 2,irc1
              rho2ns(i,lm3) = rho2ns(i,lm3) + aimag(gmatl*wf(i,l1,l2))
190         continue
          end if

200     continue
210   continue
220 continue

endsubroutine ! rhoin
