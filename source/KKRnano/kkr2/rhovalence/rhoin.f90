subroutine rhoin(ar,cden,cr,df,gmat,ek,rho2ns,irc1,nsra,efac,pz, &
fz,qz,sz,cleb,icleb,jend,iend,ekl, &
lmax, irmd, ncleb)
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

  !     INTEGER LMMAXD
  !     INTEGER LMPOTD
  !     parameter (lmmaxd= (lmaxd+1)**2)
  !     PARAMETER (LMPOTD= (LPOTD+1)**2) ! = (2*LMAX+1)**2
  !     ..
  integer lmax
  integer irmd
  integer ncleb

  !     .. Scalar Arguments ..
  double complex df,ek
  integer iend,irc1,nsra
  !     ..
  !     .. Array Arguments ..
  !     DOUBLE COMPLEX AR(LMMAXD,*),CDEN(IRMD,0:LMAXD),CR(LMMAXD,*),
  !    +               EFAC(*),EKL(0:LMAXD),FZ(IRMD,0:LMAXD),
  !    +               GMAT(LMMAXD,LMMAXD),PZ(IRMD,0:LMAXD),
  !    +               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD)
  !     DOUBLE PRECISION CLEB(*),RHO2NS(IRMD,LMPOTD)
  !     INTEGER ICLEB(NCLEB,3),JEND(LMPOTD,0:LMAXD,0:LMAXD)

  double complex ar((lmax+1)**2,*)
  double complex cden(irmd,0:lmax)
  double complex cr((lmax+1)**2,*)
  double complex efac(*)
  double complex ekl(0:lmax)
  double complex fz(irmd,0:lmax)
  double complex gmat((lmax+1)**2,(lmax+1)**2)
  double complex pz(irmd,0:lmax)
  double complex qz(irmd,0:lmax)
  double complex sz(irmd,0:lmax)
  double precision cleb(*)
  double precision rho2ns(irmd,(2*lmax+1)**2)
  integer icleb(ncleb,3)
  integer jend((2*lmax+1)**2,0:lmax,0:lmax)

  !     ..
  !     .. Local Scalars ..
  double complex czero,efac1,efac2,ffz,gmatl,ppz,v1,v2
  double precision c0ll
  integer i,ir,j,j0,j1,l,l1,l2,lm1,lm2,lm3,lm3max,ln2,ln3,m
  !     ..
  !     .. Local Arrays ..
  !     DOUBLE COMPLEX VR(LMMAXD,LMMAXD),WF(IRMD,0:LMAXD,0:LMAXD),
  !    +               WR(LMMAXD,LMMAXD)
  double complex vr((lmax+1)**2, (lmax+1)**2)
  double complex wf(irmd,0:lmax,0:lmax)
  double complex wr((lmax+1)**2, (lmax+1)**2)
  !     ..
  !     .. External Functions ..
  double complex zdotu
  external zdotu
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic atan,dimag,sqrt
  !     ..
  !     .. Save statement ..
  save czero
  !     ..
  !     .. Data statements ..
  data czero/ (0.0d0,0.0d0)/
  !     ..
  !
  integer lmmaxd
  integer lmaxd
  lmaxd = lmax
  lmmaxd= (lmaxd+1)**2

  !     C0LL = 1/sqrt(4*pi)
  c0ll = 1.0d0/sqrt(16.0d0*atan(1.0d0))
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
     efac2 = 2.0d0*efac(lm2)
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
     wr(lm1,lm1) = zdotu(lmmaxd,ar(lm1,1),lmmaxd,vr(lm1,1),lmmaxd)/ &
     (efac1*efac1)
     do 60 lm2 = 1,lm1 - 1
       !
       !---> using symmetry of gaunt coeffients
       !
       efac2 = efac(lm2)
       wr(lm1,lm2) = (zdotu(lmmaxd,ar(lm1,1),lmmaxd,vr(lm2,1), &
       lmmaxd)+zdotu(lmmaxd,ar(lm2,1),lmmaxd,vr(lm1,1), &
       lmmaxd))/ (efac1*efac2)
60   continue
70 continue

   !
   !---> set up array wf(l1,l2) = pz(l1)*pz(l2)
   !
   if (nsra.eq.2) then
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
     gmatl = czero
     do 140 m = -l,l
       lm1 = l* (l+1) + m + 1
       gmatl = gmatl + wr(lm1,lm1)
140  continue
     !
     if (nsra.eq.2) then
       do 150 i = 2,irc1
         ppz = pz(i,l)
         ffz = fz(i,l)
         cden(i,l) = ppz* (gmatl*ppz+ekl(l)*qz(i,l)) + &
         ffz* (gmatl*ffz+ekl(l)*sz(i,l))
         rho2ns(i,1) = rho2ns(i,1) + c0ll*dimag(df*cden(i,l))
150    continue

     else
       do 160 i = 2,irc1
         ppz = pz(i,l)
         cden(i,l) = ppz* (gmatl*ppz+ekl(l)*qz(i,l))
         rho2ns(i,1) = rho2ns(i,1) + c0ll*dimag(df*cden(i,l))
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
          if (j1.ne.0) then
            !
            gmatl = czero
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
              rho2ns(i,lm3) = rho2ns(i,lm3) + dimag(gmatl*wf(i,l1,l2))
190         continue
          end if

200     continue
210   continue
220 continue

  end
