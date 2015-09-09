module ValenceDensity_mod
  implicit none
  private
  public :: rhoval, rhoval0
  
  double complex, parameter :: cone=(1.d0,0.d0), zero=(0.d0,0.d0), ci=(0.d0,1.d0)
  double precision, parameter :: cvlight=274.0720442d0
  
  contains

  subroutine irwns(cr, dr, efac, qns, vnspll, icst, ipan, ircut, nsra,  &
                    pzlm, qzlm, pzekdr, qzekdr, cder, cmat, dder, dmat,  &
                    irmind, irmd, ipand, lmmaxd)
    use SingleSiteHelpers_mod, only: csinwd, wfint, wfint0
    !-----------------------------------------------------------------------
    !     determines the irregular non spherical wavefunctions in the n-th.
    !       born approximation ( n given by input parameter icst ) .
    !
    !
    !     using the wave functions pz and qz ( regular and irregular
    !       solution ) of the spherically averaged potential , the ir-
    !       regular wavefunction qns is determined by
    !
    !           qns(ir,lm1,lm2) = cr(ir,lm1,lm2)*pz(ir,l1)
    !
    !                                   + dr(ir,lm1,lm2)*qz(ir,l1)
    !
    !      the matrices cr and dr are determined by integral equations
    !        containing qns and only the non spherical contributions of
    !        the potential , stored in vinspll . these integral equations
    !        are solved iteratively with born approximation up to given n.
    !
    !     the original way of writing the cr and dr matrices in the equa-
    !        tion above caused numerical troubles . therefore here are used
    !        rescaled cr and dr matrices (compare subroutine wftsca):
    !
    !              ~
    !              cr(ir,lm1,lm2) = sqrt(e)**(l1+l2)
    !                             * cr(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
    !
    !              ~
    !              dr(ir,lm1,lm2) = sqrt(e)**(l2-l1)
    !                             * dr(ir,lm1,lm2)*((2*l1-1)!!/(2*l2-1)!!)
    !
    !     attention :  the sign of the dr matrix is changed to reduce the
    !     ===========  number of floating point operations
    !
    !     modified for the use of shape functions
    !
    !                              (see notes by b.drittler)
    !
    !                                b.drittler   mar.  1989
    !-----------------------------------------------------------------------
    !     modified by r. zeller      aug. 1994
    !-----------------------------------------------------------------------
    integer, intent(in) :: icst,ipan,ipand,irmd,irmind,lmmaxd,nsra
    double complex, intent(out)   :: cder(lmmaxd,lmmaxd,irmind:irmd), dder(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(inout) :: cmat(lmmaxd,lmmaxd,irmind:irmd), dmat(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(out)   :: cr(lmmaxd,lmmaxd), dr(lmmaxd,lmmaxd)
    double complex, intent(inout) :: qns(lmmaxd,lmmaxd,irmind:irmd,2)
    double complex, intent(in)    :: efac(lmmaxd)
    double complex, intent(in)    :: pzlm(lmmaxd,irmind:irmd,2)
    double complex, intent(in)    :: pzekdr(lmmaxd,irmind:irmd,2), qzekdr(lmmaxd,irmind:irmd,2)
    double complex, intent(in)    :: qzlm(lmmaxd,irmind:irmd,2)

    double precision, intent(in) :: vnspll(lmmaxd,lmmaxd,irmind:irmd)
    integer, intent(in) :: ircut(0:ipand)
    
    double complex :: efac2
    integer :: i, ir, irc1,j, lm
    
    irc1 = ircut(ipan)

    do i = 0, icst
      !---> set up integrands for i-th born approximation
      if (i == 0) then
        call wfint0(cder,dder,qzlm,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
      else
        call wfint(qns,cder,dder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
      endif ! first born iteration

      !---> call integration subroutines
      call csinwd(cder,cmat,lmmaxd**2,irmind,irmd,ipan,ircut)
      call csinwd(dder,dmat,lmmaxd**2,irmind,irmd,ipan,ircut)
      do ir = irmind, irc1
        do lm = 1, lmmaxd
          dmat(lm,lm,ir) = dmat(lm,lm,ir) - cone
        enddo ! lm
      enddo ! ir

      !---> calculate non sph. wft. in i-th born approximation
      do j = 1, nsra
        do ir = irmind, irc1
          do lm = 1, lmmaxd
            qns(:,lm,ir,j) = cmat(:,lm,ir)*pzlm(:,ir,j) - dmat(:,lm,ir)*qzlm(:,ir,j)
          enddo ! lm
        enddo ! ir
      enddo ! j
      
    enddo ! i

    !---> store c- and d- matrix
    cr(:,:) =  cmat(:,:,irmind)
    dr(:,:) = -dmat(:,:,irmind)

    !---> rescale with efac
    do j = 1, nsra
      do lm = 1, lmmaxd
        efac2 = 1.d0/efac(lm)
        do ir = irmind, irc1
          qns(:,lm,ir,j) = qns(:,lm,ir,j)*efac2
        enddo ! ir
      enddo ! lm
    enddo ! j

  endsubroutine irwns
  
  
  
subroutine pnsqns(ar, cr, dr, drdi, ek, icst, pz, qz, fz, sz, pns, qns, nsra, &
  vins, ipan, ircut, cleb, icleb, iend, loflm, lkonv, ispin, ldau, nldau, lldau, &
  wmldau, wmldauav, ldaucut, lmaxd, nspind, irmd, irnsd, ipand, ncleb)
  use singlesite_mod, only: regns
  use SingleSiteHelpers_mod, only: vllns, wftsca

  integer, intent(in) :: lmaxd, nspind, irmd, irnsd, ipand, ncleb
  double complex, intent(in) :: ek
  integer, intent(in) :: icst,iend,ipan,lkonv,nsra,nldau,ispin
  logical, intent(in) :: ldau
  double complex, intent(out) :: ar((lmaxd+1)**2,(lmaxd+1)**2)
  double complex, intent(out) :: cr((lmaxd+1)**2,(lmaxd+1)**2), dr((lmaxd+1)**2,(lmaxd+1)**2)
  double complex, intent(in) :: fz(irmd,0:lmaxd)
  double complex, intent(in) :: sz(irmd,0:lmaxd)
  double complex, intent(inout) :: pns((lmaxd+1)**2,(lmaxd+1)**2,(irmd-irnsd):irmd,2)
  double complex, intent(inout) :: qns((lmaxd+1)**2,(lmaxd+1)**2,(irmd-irnsd):irmd,2)
  double complex, intent(in) :: pz(irmd,0:lmaxd)
  double complex, intent(in) :: qz(irmd,0:lmaxd)

  double precision, intent(in) :: cleb(ncleb,2)
  double precision, intent(in) :: drdi(irmd)
  double precision, intent(in) :: vins(irmd-irnsd:irmd,(2*lmaxd+1)**2)
  double precision, intent(in) :: wmldau(2*lmaxd+1,2*lmaxd+1,lmaxd+1,nspind)
  double precision, intent(in) :: ldaucut(irmd)

  integer, intent(in) :: icleb(ncleb,3), ircut(0:ipand), loflm(*)
  integer, intent(in) :: lldau(lmaxd+1)

  integer            ir,irc1,lm1,lm2,lmmkonv,m1,m2, lmlo,lmhi,ildau
  double complex     cmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
  double complex     dmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
  double complex     efac((lmaxd+1)**2)
  double complex     pzekdr((lmaxd+1)**2,irmd-irnsd:irmd,2)
  double complex     pzlm((lmaxd+1)**2,(irmd-irnsd):irmd,2)
  double complex     qzekdr((lmaxd+1)**2,(irmd-irnsd):irmd,2)
  double complex     qzlm((lmaxd+1)**2,(irmd-irnsd):irmd,2)
  double complex     tmatll((lmaxd+1)**2,(lmaxd+1)**2)
  double precision :: vnspll((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)

  double precision   wmldauav(lmaxd+1)
  integer             lmmaxd
  integer             irmind

  irmind = irmd-irnsd
  lmmaxd = (lmaxd+1)**2

  irc1 = ircut(ipan)

  ! calculate v_ll' potential
  ! v_{ll'} = \sum_{l''} c_{l l' l''} v_{l''}
  call vllns(vnspll,vins,cleb,icleb,iend, lmaxd, irmd, irnsd, ncleb)


  if (lkonv /= lmaxd) then
    lmmkonv = (lkonv+1)**2
    do lm1 = 1, lmmaxd
      do lm2 = lmmkonv + 1, lmmaxd
        do ir = irmind, irmd
          vnspll(lm2,lm1,ir) = 0.d0
          vnspll(lm1,lm2,ir) = 0.d0
        enddo ! ir
      enddo
    enddo
  else
    lmmkonv = lmmaxd
  endif

  !-----------------------------------------------------------------------
  ! lda+u
  ! add wldau to non-spherical porential vins in case of lda+u
  ! use the average wldau (=wldauav) and the deviation
  ! of wldau from this. use the deviation in the born series
  ! for the non-spherical wavefunction, while the average is
  ! used for the spherical wavefunction.
  !
  if (ldau) then
    do ildau = 1, nldau
      if (lldau(ildau) >= 0) then
        !
        lmlo = lldau(ildau)**2 + 1
        lmhi = (lldau(ildau)+1)**2
        !
        do ir = irmind,irmd
          !
          ! -> first add wldau to all elements ...
          !
          do lm2 = lmlo,lmhi
            m2 = lm2 - lmlo + 1
            do lm1 = lmlo,lmhi
              m1 = lm1 - lmlo + 1
              vnspll(lm1,lm2,ir) = vnspll(lm1,lm2,ir) + wmldau(m1,m2,ildau,ispin)*ldaucut(ir)
            enddo ! lm1
            !
            ! ... and then subtract average from diag. elements
            !
            vnspll(lm2,lm2,ir) = vnspll(lm2,lm2,ir) - wmldauav(ildau)*ldaucut(ir)
          enddo ! lm2
        enddo ! ir
      endif ! lldau
    enddo ! ildau
  endif ! ldau
  !
  ! lda+u
  !-----------------------------------------------------------------------


  !
  !---> get wfts of same magnitude by scaling with efac
  !
  call wftsca(drdi,efac,pz,qz,fz,sz,nsra,pzlm,qzlm,pzekdr,qzekdr, ek,loflm,irmind,irmd,lmaxd,lmmaxd)
  !
  !---> determine the irregular non sph. wavefunction
  !
  call irwns(cr,dr,efac,qns,vnspll,icst,ipan,ircut,nsra,pzlm,qzlm, pzekdr,qzekdr,qns(1,1,irmind,1),cmat, qns(1,1,irmind,2),dmat,irmind,irmd,ipand,lmmaxd)
  !
  !---> determine the regular non sph. wavefunction
  !
  call regns(ar,tmatll,efac,pns,vnspll,icst,ipan,ircut,pzlm,qzlm, pzekdr,qzekdr,ek,pns(1,1,irmind,1),cmat, pns(1,1,irmind,2),dmat,nsra,irmind,irmd,ipand,lmmaxd)

endsubroutine pnsqns



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
  !              qns(ir,lm1,lm2) = pz(ir,l1) * cr(lm1,lm2) + qz(ir,l1) * dr(lm1,lm2)
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
                fz, qz, sz, cleb, icleb, jend, iend, ekl, lmax, irmd, ncleb)
  
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

  double complex, external :: zdotu ! from blas
  double complex, parameter :: zero = (0.d0, 0.d0)
  double complex :: efac1,efac2,ffz,gmatl,ppz,v1,v2
  double precision :: c0ll
  integer :: ir,j,j0,j1,l,l1,l2,lm1,lm2,lm3,lm3max,ln2,ln3,m
  double complex :: vr((lmax+1)**2,(lmax+1)**2)
  double complex :: wf(irmd,0:lmax,0:lmax)
  double complex :: wr((lmax+1)**2,(lmax+1)**2)
  integer :: lmmaxd, lmaxd
  
  lmaxd = lmax
  lmmaxd = (lmaxd+1)**2

  !
  !---> set up array wr(lm1,lm2), use first vr
  !
  do lm2 = 1, lmmaxd
    ln2 = lm2
    v2 = efac(lm2)*efac(lm2)*gmat(ln2,ln2)
    vr(1:lmmaxd,lm2) = ek*cr(1:lmmaxd,lm2) + v2*ar(1:lmmaxd,lm2)
  enddo ! lm2

  !
  !---> using symmetry of structural green function
  !
  do lm2 = 2, lmmaxd
    ln2 = lm2
    efac2 = 2.d0*efac(lm2)
    do lm3 = 1, lm2-1
      ln3 = lm3
      v1 = efac2*gmat(ln3,ln2)*efac(lm3)
      vr(1:lmmaxd,lm2) = vr(1:lmmaxd,lm2) + v1*ar(1:lmmaxd,lm3)
    enddo ! lm3
  enddo ! lm2

  do lm1 = 1, lmmaxd
    efac1 = efac(lm1)
    wr(lm1,lm1) = zdotu(lmmaxd,ar(lm1,1),lmmaxd,vr(lm1,1),lmmaxd)/(efac1*efac1)
    do lm2 = 1, lm1-1
      !
      !---> using symmetry of gaunt coeffients
      !
      efac2 = efac(lm2)
      wr(lm1,lm2) = (zdotu(lmmaxd,ar(lm1,1),lmmaxd,vr(lm2,1), lmmaxd) + zdotu(lmmaxd,ar(lm2,1),lmmaxd,vr(lm1,1),lmmaxd))/(efac1*efac2)
    enddo ! lm2
  enddo ! lm1

  !
  !---> set up array wf(l1,l2) = pz(l1)*pz(l2)
  !
  if (nsra == 2) then
    do l1 = 0, lmaxd
      do l2 = 0, l1
        wf(2:irc1,l1,l2) = pz(2:irc1,l1)*pz(2:irc1,l2) + fz(2:irc1,l1)*fz(2:irc1,l2)
      enddo ! l2
    enddo ! l1
  else
    do l1 = 0, lmaxd
      do l2 = 0, l1
        wf(2:irc1,l1,l2) = pz(2:irc1,l1)*pz(2:irc1,l2)
      enddo ! l2
    enddo ! l1
  endif
  !
  !---> first calculate only the spherically symmetric contribution
  !     remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
  !
  c0ll = 1.d0/sqrt(16.d0*atan(1.d0))

  do l = 0, lmaxd
  
    gmatl = zero
    do m = -l, l
      lm1 = l*(l+1)+m+1
      gmatl = gmatl + wr(lm1,lm1)
    enddo ! m
    !
    if (nsra == 2) then
      do ir = 2, irc1
        ppz = pz(ir,l)
        ffz = fz(ir,l)
        cden(ir,l) = ppz*(gmatl*ppz + ekl(l)*qz(ir,l)) + ffz*(gmatl*ffz + ekl(l)*sz(ir,l))
        rho2ns(ir,1) = rho2ns(ir,1) + c0ll*dimag(df*cden(ir,l))
      enddo ! ir

    else
      do ir = 2, irc1
        ppz = pz(ir,l)
        cden(ir,l) = ppz*(gmatl*ppz + ekl(l)*qz(ir,l))
        rho2ns(ir,1) = rho2ns(ir,1) + c0ll*dimag(df*cden(ir,l))
      enddo ! ir
    endif

  enddo ! l

  !
  !---> calculate the non spherically symmetric contribution
  !        to speed up the pointer jend generated in gaunt is used
  !        remember that the wavefunctions are l and not lm dependent
  !
  lm3max = icleb(iend,3)
  j0 = 1
  !
  do lm3 = 2, lm3max
    do l1 = 0, lmaxd
      do l2 = 0, l1
        !
        j1 = jend(lm3,l1,l2)
        if (j1 /= 0) then
          !
          gmatl = zero
          !
          !---> sum over m1,m2 for fixed lm3,l1,l2
          !
          do j = j0,j1
            lm1 = icleb(j,1)
            lm2 = icleb(j,2)
            gmatl = gmatl + cleb(j)*wr(lm1,lm2)
          enddo ! j
          !
          j0 = j1 + 1
          !
          gmatl = df*gmatl
          rho2ns(2:irc1,lm3) = rho2ns(2:irc1,lm3) + dimag(gmatl*wf(2:irc1,l1,l2))
        endif

      enddo ! l2
    enddo ! l1
  enddo ! lm3

endsubroutine ! rhoin




subroutine rhons(den, df, drdi, gmat, ek, rho2ns, ipan, ircut, thetas, &
  ifunm, lmsp, nsra, qns, pns, ar, cr, pz, fz, qz, sz, cleb, &
  icleb, jend, iend, ekl, lmax, irmd, irnsd, irid, ipand, nfund, ncleb)
  !-----------------------------------------------------------------------
  !
  !     the charge density is developed in spherical harmonics :
  !
  !             rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)
  !
  !          rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
  !                                                           unit sphere)
  !     in the case of spin-polarization :
  !       the spin density is developed in spherical harmonics :
  !
  !            sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)
  !
  !         sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
  !                                                           unit sphere)
  !     n(r,e) is developed in
  !
  !        n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }
  !
  !     therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
  !     has to be used to calculate the lm-contribution of the charge
  !     density .
  !
  !
  !     calculate the valence density of states , in the spin-polarized
  !      case spin dependent .
  !     recognize that the density of states is always complex also in
  !      the case of "real-energy-integation" (ief>0) since in that case
  !      the energy integration is done parallel to the real energy axis
  !      but not on the real energy axis .
  !     in the last energy-spin loop the l-contribution of the valence
  !      charge is calculated .
  !
  !                               b.drittler   aug. 1988
  !
  !     modified for the use of shape functions
  !
  !     attention : irmin + 3 has to be less then imt
  !                 if shape functions are used
  !
  !                               b.drittler   july 1989
  !-----------------------------------------------------------------------
  use Quadrature_mod, only: simpson
  
  integer lmax
  integer irmd
  integer ncleb
  integer irnsd
  integer irid
  integer ipand
  integer nfund

  !     integer irmind
  !     parameter (irmind=irmd-irnsd)
  !     integer lmpotd,lmmaxd
  !     parameter (lmpotd= (lpotd+1)**2) ! lmpotd= (2*lmax+1)**2
  !     parameter (lmmaxd= (lmaxd+1)**2)
  !     integer lmaxd1
  !     parameter (lmaxd1= lmaxd+1)
  !     ..
  !     .. scalar arguments ..
  double complex df,ek
  integer iend,ipan,nsra
  !     ..
  !     .. array arguments ..
  !     double complex ar(lmmaxd,lmmaxd),cr(lmmaxd,lmmaxd),den(0:lmaxd1),
  !    +               ekl(0:lmaxd),fz(irmd,0:lmaxd),gmat(lmmaxd,lmmaxd),
  !    +               pns(lmmaxd,lmmaxd,irmind:irmd,2),pz(irmd,0:lmaxd),
  !    +               qns(lmmaxd,lmmaxd,irmind:irmd,2),qz(irmd,0:lmaxd),
  !    +               sz(irmd,0:lmaxd)
  !     double precision cleb(*),drdi(irmd),rho2ns(irmd,lmpotd),
  !    +                 thetas(irid,nfund)
  !     integer icleb(ncleb,3),ifunm(*),ircut(0:ipand),
  !    +        jend(lmpotd,0:lmaxd,0:lmaxd),lmsp(*)

  double complex ar((lmax+1)**2,(lmax+1)**2)
  double complex cr((lmax+1)**2,(lmax+1)**2)
  double complex den(0:lmax+1)
  double complex ekl(0:lmax)
  double complex fz(irmd,0:lmax)
  double complex gmat((lmax+1)**2,(lmax+1)**2)
  double complex pns((lmax+1)**2,(lmax+1)**2,irmd-irnsd:irmd,2)
  double complex pz(irmd,0:lmax)
  double complex qns((lmax+1)**2,(lmax+1)**2,irmd-irnsd:irmd,2)
  double complex qz(irmd,0:lmax)
  double complex sz(irmd,0:lmax)

  !     double precision cleb(*),drdi(irmd),rho2ns(irmd,lmpotd),
  !    +                 thetas(irid,nfund)
  !     integer icleb(ncleb,3),ifunm(*),ircut(0:ipand),
  !    +        jend(lmpotd,0:lmaxd,0:lmaxd),lmsp(*)

  double precision cleb(*)
  double precision drdi(irmd)
  double precision rho2ns(irmd,(2*lmax+1)**2)
  double precision thetas(irid,nfund)

  integer icleb(ncleb,3)
  integer ifunm(*)
  integer ircut(0:ipand)
  integer jend((2*lmax+1)**2,0:lmax,0:lmax)
  integer lmsp(*)
  
  double complex v1
  integer imt1,l,lm,m
  double complex cden(irmd,0:lmax)
  double complex cdenns(irmd)
  double complex efac((lmax+1)**2)
  
  integer irmind
  irmind=irmd-irnsd
  !
  !---> set up efac(lm) = sqrt(e))**l/(2l - 1)!!
  !
  efac(1) = 1.d0
  v1 = 1.d0
  do l = 1, lmax
    v1 = v1*ek/dble(2*l-1)
    do m = -l, l
      lm = l* (l+1) + m + 1
      efac(lm) = v1
    enddo ! m
  enddo ! l
   !
   imt1 = ircut(1)
   !

   call rhoout(cden,df,gmat,ek,pns,qns,rho2ns,thetas,ifunm,ipan, imt1,lmsp,cdenns,nsra,cleb,icleb,iend, lmax, irmd, irnsd, irid, nfund, ncleb)

   call rhoin(ar,cden,cr,df,gmat,ek,rho2ns,irmind,nsra,efac,pz,fz, qz,sz,cleb,icleb,jend,iend,ekl, lmax, irmd, ncleb)

   !
   !---> calculate complex density of states
   !
   do l = 0, lmax
     !
     !---> call integration subroutine
     !
     den(l) = simpson(cden(1:,l), ipan, ircut, drdi)
   enddo ! l

   if (ipan > 1) den(lmax+1) = simpson(cdenns, ipan, ircut, drdi)

 endsubroutine
 
 
 
 
subroutine rhoout(cden, df, gmat, ek, pns, qns, rho2ns, thetas, ifunm,  &
                  ipan1, imt1, lmsp, cdenns, nsra, cleb, icleb, iend,  &
                  lmaxd, irmd, irnsd, irid, nfund, ncleb)
  !-----------------------------------------------------------------------
  !
  !     calculates the charge density from r(irmin) to r(irc)
  !      in case of a non spherical input potential .
  !
  !     fills the array cden for the complex density of states
  !
  !     attention : the gaunt coeffients are stored in index array
  !                   (see subroutine gaunt)
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
  

  integer, intent(in) :: lmaxd, irmd, ncleb, irnsd, irid, nfund
  double complex, intent(in) :: df, ek
  integer, intent(in) :: iend, imt1, ipan1, nsra
  double complex, intent(out) :: cden(irmd,0:*)
  double complex, intent(inout) :: cdenns(*)
  double complex, intent(in) :: gmat((lmaxd+1)**2,(lmaxd+1)**2)
  double complex, intent(in) :: pns((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd,2)
  double complex, intent(in) :: qns((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd,2)
  double precision, intent(in) :: cleb(*)
  double precision, intent(inout) :: rho2ns(irmd,(2*lmaxd+1)**2)
  double precision, intent(in) :: thetas(irid,nfund)
  integer, intent(in) :: icleb(ncleb,3), ifunm(*), lmsp(*)

  external :: zgemm ! from blas
  double complex :: cltdf
  double precision :: c0ll
  integer :: ifun, ir, j, l1, lm1, lm2, lm3, m1
  integer :: lmmaxd, irmind
  double complex :: qnsi((lmaxd+1)**2,(lmaxd+1)**2)
  double complex :: wr((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)

  lmmaxd = (lmaxd+1)**2
  irmind = irmd-irnsd

  !------------------------------------------------------------------
  !
  !---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) } summed over lm3
  !---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) } summed over lm3
  do ir = irmind+1, irmd
  
    qnsi(1:lmmaxd,1:lmmaxd) = qns(1:lmmaxd,1:lmmaxd,ir,1)

    call zgemm('n','n',lmmaxd,lmmaxd,lmmaxd,cone,pns(:,1,ir,1),lmmaxd,gmat,lmmaxd,ek,qnsi,lmmaxd)
    call zgemm('n','t',lmmaxd,lmmaxd,lmmaxd,cone,pns(:,1,ir,1),lmmaxd,qnsi,lmmaxd,zero,wr(:,1,ir),lmmaxd)

    if (nsra == 2) then
      qnsi(1:lmmaxd,1:lmmaxd) = qns(1:lmmaxd,1:lmmaxd,ir,2)

      call zgemm('n','n',lmmaxd,lmmaxd,lmmaxd,cone,pns(:,1,ir,2),lmmaxd,gmat,lmmaxd,ek,qnsi,lmmaxd)
      call zgemm('n','t',lmmaxd,lmmaxd,lmmaxd,cone,pns(:,1,ir,2),lmmaxd,qnsi,lmmaxd,cone,wr(:,1,ir),lmmaxd)

    endif ! nsra

    do lm1 = 1, lmmaxd
      do lm2 = 1, lm1-1
        wr(lm1,lm2,ir) = wr(lm1,lm2,ir) + wr(lm2,lm1,ir) ! symmetrize
      enddo ! lm2
    enddo ! lm1
    
  enddo ! ir
  
  c0ll = 1.d0/sqrt(16.d0*atan(1.d0))
  
  cden(1:irmd,0:lmaxd) = zero !---> initialize array for complex charge density
  !
  !---> first calculate only the spherically symmetric contribution
  !
  do l1 = 0, lmaxd
    do m1 = -l1, l1
      lm1 = l1*(l1+1)+m1+1
      do ir = irmind+1, irmd
        cden(ir,l1) = cden(ir,l1) + wr(lm1,lm1,ir) !---> fill array for complex density of states
      enddo ! ir
    enddo ! m1
    !
    !---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi) == c0ll
    !
    rho2ns(irmind+1:irmd,1) = rho2ns(irmind+1:irmd,1) + c0ll*aimag(cden(irmind+1:irmd,l1)*df)
    !
    if (ipan1 > 1) cden(imt1+1:irmd,l1) = cden(imt1+1:irmd,l1)*thetas(1:irmd-imt1,1)*c0ll

  enddo ! l1

  if (ipan1 > 1) cdenns(1:irmd) = 0.d0

  do j = 1, iend
    lm1 = icleb(j,1)
    lm2 = icleb(j,2)
    lm3 = icleb(j,3)
    cltdf = df*cleb(j)
    !
    !---> calculate the non spherically symmetric contribution
    !
    do ir = irmind+1, irmd
      rho2ns(ir,lm3) = rho2ns(ir,lm3) + aimag(cltdf*wr(lm1,lm2,ir))
    enddo ! ir

    if (ipan1 > 1 .and. lmsp(lm3) > 0) then
    
      ifun = ifunm(lm3)
      do ir = imt1+1, irmd
        cdenns(ir) = cdenns(ir) + cleb(j)*thetas(ir-imt1,ifun)*wr(lm1,lm2,ir)
      enddo ! ir

    endif ! ipan1 ...

  enddo ! j

endsubroutine rhoout



  subroutine rhoval(ldorhoef,icst,ielast,nsra, &
                    ispin,nspin, &
                    ez,wez,drdi,r,irmin, &
                    vins,visp,zat,ipan,ircut, &
                    thetas,ifunm,lmsp,rho2ns,r2nef,den, &
                    espv,cleb,loflm,icleb,iend,jend, &
                    gmatn, &                                 ! input
                    ldau,nldau,lldau,phildau,wmldau, &       ! input
                    dmatldau, &                              ! output
                    iemxd, lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)
    use SingleSiteHelpers_mod, only: cradwf, wfmesh

    integer, intent(in) :: iemxd, lmaxd, irmd, ncleb, irnsd, irid, ipand, nfund
    double precision, intent(in) :: zat
    integer, intent(in) :: icst, ielast, iend, ipan, ispin, nspin, nsra, irmin, nldau
    logical, intent(in) :: ldorhoef, ldau

    double complex, intent(in) :: den(0:lmaxd+1,iemxd)
    double complex, intent(in) :: ez(iemxd)
    double complex, intent(in) :: wez(iemxd)
    double complex, intent(in) :: phildau(irmd,lmaxd+1)
    double complex, intent(in) :: dmatldau(2*lmaxd+1,2*lmaxd+1,nspin,lmaxd+1)
    double complex, intent(in) :: gmatn((lmaxd+1)**2,(lmaxd+1)**2,iemxd,nspin)

    double precision, intent(in) :: cleb(ncleb,2)
    double precision, intent(in) :: drdi(irmd)
    double precision, intent(out) :: espv(0:lmaxd+1,1)
    double precision, intent(in) :: r(irmd)
    double precision, intent(out) :: rho2ns(irmd,(2*lmaxd+1)**2,2)
    double precision, intent(out) :: r2nef(irmd,(2*lmaxd+1)**2,2)
    double precision, intent(in) :: thetas(irid,nfund)
    double precision, intent(in) :: vins(irmd-irnsd:irmd,(2*lmaxd+1)**2)
    double precision, intent(in) :: visp(irmd)
    double precision, intent(in) :: wmldau(2*lmaxd+1,2*lmaxd+1,lmaxd+1,nspin)

    integer, intent(in) :: icleb(ncleb,3)
    integer, intent(in) :: ifunm((4*lmaxd+1)**2)
    integer, intent(in) :: ircut(0:ipand)
    integer, intent(in) :: jend((2*lmaxd+1)**2,0:lmaxd,0:lmaxd)
    integer, intent(in) :: lmsp((4*lmaxd+1)**2)
    integer, intent(in) :: loflm((2*lmaxd+1)**2)
    integer, intent(in) :: lldau(lmaxd+1)

    external :: daxpy, dscal ! from blas
    
    double complex :: df, eryd, ek
    integer :: idm, ie, l, lmlo, lmhi, mmax, im, ildau

    ! the following arrays are local
    double complex :: alpha(0:lmaxd)
    double complex :: ar((lmaxd+1)**2,(lmaxd+1)**2)
    double complex :: dr((lmaxd+1)**2,(lmaxd+1)**2)
    double complex :: cr((lmaxd+1)**2,(lmaxd+1)**2)
    double complex :: ekl(0:lmaxd)
    double complex :: fz(irmd,0:lmaxd)
    double complex :: gmatll((lmaxd+1)**2,(lmaxd+1)**2)

    double complex :: qz(irmd,0:lmaxd)
    double complex :: sz(irmd,0:lmaxd)
    double complex :: tmat(0:lmaxd)
    double complex :: pz(irmd,0:lmaxd)

    double precision :: rs(irmd,0:lmaxd)
    double precision :: s(0:lmaxd)
    double precision :: ldaucut(irmd)
    double precision :: wmldauav(lmaxd+1)
    double complex :: dendum(0:lmaxd+1)

    ! dynamically allocate large arrays
    double complex, allocatable :: pns(:,:,:,:), qns(:,:,:,:) ! dims: ((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd,2)

    integer :: memory_stat
    integer :: lmmaxd, lmaxd1, nspind, lmpotd

    lmpotd = (2*lmaxd+1)**2
    lmmaxd= (lmaxd+1)**2
    lmaxd1= lmaxd+1
    nspind = nspin

    allocate(pns(lmmaxd,lmmaxd,irmd-irnsd:irmd,2), qns(lmmaxd,lmmaxd,irmd-irnsd:irmd,2), stat=memory_stat)
    if (memory_stat /= 0) then
      write(*,*) "rhoval: fatal error, failure to allocate memory, probably out of memory."
      stop
    endif

    !-----------------------------------------------------------------------
    ! initialise local variables to be on the safe side
    ek = zero
    pns = zero
    qns = zero
    alpha = zero
    ar = zero
    dr = zero
    cr = zero
    ekl = zero
    fz = zero
    !gmatll = zero initialised further down
    qz = zero
    sz = zero
    tmat = zero
    pz = zero
    rs = 0.d0
    s = 0.d0
    ldaucut = 0.d0
    wmldauav = 0.d0
    dendum = zero

    !-----------------------------------------------------------------------
    ! ldau

    if (ldau) then
      do ildau = 1, nldau
        wmldauav(ildau) = 0.d0
        lmlo = lldau(ildau)**2 + 1
        lmhi = (lldau(ildau)+1)**2
        mmax = lmhi - lmlo + 1
        do im = 1, mmax
          wmldauav(ildau) = wmldauav(ildau) + wmldau(im,im,ildau,ispin)
        enddo ! im
        wmldauav(ildau) = wmldauav(ildau)/dble(mmax)
      enddo ! ildau
      
      ! -> note: application if wldau makes the potential discontinuous.
      !    a cutoff can be used if necessary to make the potential continuous
      !    for example (array bounds should be adjusted):
        
    !    if(test('cutoff  ')) then
    !      do ir = 1,irmd
    !        ldaucut(ir) = ( 1.d0 + dexp( 20.d0*(r(ir)-r(349)) ) ) * &
    !        ( 1.d0 + dexp( 20.d0*(r(276)-r(ir)) ) )
    !        ldaucut(ir) = 1d0/ldaucut(ir)
    !      enddo
    !    else
    !      do ir = 1,irmd
    !        ldaucut(ir) = 1.d0
    !      enddo
    !    endif
    endif

    ! ldau
    !-----------------------------------------------------------------------

    rho2ns(:,:,ispin) = 0.d0
    r2nef(:,:,ispin) = 0.d0

    espv = 0.d0

    do ie = 1, ielast

      gmatll(:,:) = gmatn(:,:,ie,ispin)

      eryd = ez(ie)
      df = wez(ie)/dble(nspin)
      
      !=======================================================================
      call wfmesh(eryd,ek,cvlight,nsra,zat,r,s,rs,ircut(ipan),irmd,lmaxd)

      call cradwf(eryd,ek,nsra,alpha,ipan,ircut,cvlight,rs,s, pz,fz,qz,sz,tmat,visp,drdi,r,zat, &
                  ldau,nldau,lldau,wmldauav,ldaucut, lmaxd, irmd, ipand)
      !-----------------------------------------------------------------------
      ! non-spherical
      
      call pnsqns(ar,cr,dr,drdi,ek,icst,pz,qz,fz,sz, pns,qns,nsra,vins,ipan,ircut, &
                  cleb,icleb,iend,loflm,lmaxd,ispin, ldau,nldau,lldau, &
                  wmldau,wmldauav,ldaucut, lmaxd, nspind, irmd, irnsd, ipand, ncleb)

      do l = 0, lmaxd
        ekl(l) = ek*dble(2*l+1)
      enddo ! l
      !-----------------------------------------------------------------------
      call rhons(den(0,ie),df,drdi,gmatll,ek, rho2ns(1,1,ispin),ipan,ircut,thetas,ifunm,lmsp, &
                 nsra,qns,pns,ar,cr,pz,fz,qz,sz,cleb(1,1),icleb, jend,iend,ekl, lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)
      !-----------------------------------------------------------------------


      do l = 0, lmaxd1
        espv(l,1) = espv(l,1) + dimag(eryd*den(l,ie)*df)
      enddo ! l
      
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !     get the charge at the fermi energy (ielast)
      !     call with the energy weight cone --> not overwrite df
      !          with the dummy dendum       --> not overwrite den
      
      if ( (ie == ielast) .and. (ldorhoef) ) then
        call rhons(dendum,cone,drdi,gmatll,ek, r2nef(1,1,ispin),ipan,ircut,thetas,ifunm,lmsp, nsra,qns,pns,ar,cr,pz,fz,qz,sz, & 
                   cleb(1,1),icleb, jend,iend,ekl, lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)
      endif

      
      if (ldau .and. nldau >= 1) then
        call ldaudmat(df,pz,qz,pns,qns,ar,cr,dr,gmatll, ipan,ircut,drdi,ek, irmin,lldau,phildau,nldau, dmatldau,ispin, lmaxd, nspind, irmd, irnsd, ipand)
      endif

    enddo ! ie

    ! this should really be separated into another routine
    if (ispin == 2) then
      idm = irmd*lmpotd
      call dscal(idm,2.d0,rho2ns(1,1,1),1)
      call daxpy(idm,-.5d0,rho2ns(1,1,1),1,rho2ns(1,1,2),1)
      call daxpy(idm,1.d0,rho2ns(1,1,2),1,rho2ns(1,1,1),1)
      
      ! --> do the same at the fermi energy
      
      call dscal(idm,2.d0,r2nef(1,1,1),1)
      call daxpy(idm,-.5d0,r2nef(1,1,1),1,r2nef(1,1,2),1)
      call daxpy(idm,1.d0,r2nef(1,1,2),1,r2nef(1,1,1),1)
    endif ! ispin == 2

    deallocate(pns, qns)
    
  endsubroutine rhoval
  

  subroutine rhoval0(ez, wez, drdi, r, ipan, ircut, thetas, dos0, dos1, lmaxd, irmd, irid, ipand, nfund)
    use SingleSiteHelpers_mod, only: beshan
    use Quadrature_mod, only: simpson
    
    integer, intent(in) :: lmaxd, irmd, irid, ipand, nfund, ipan
    double complex, intent(in) :: ez,wez
    double complex, intent(inout) :: dos0, dos1
    double precision, intent(in) :: drdi(irmd), r(irmd), thetas(irid,nfund)
    integer, intent(in) :: ircut(0:ipand)

    double complex ek,ciek
    double precision c0ll
    integer ir,l,l1,imt1
    double complex pz(irmd,0:lmaxd),qz(irmd,0:lmaxd)
    double complex bessjw(0:lmaxd+1)
    double complex bessyw(0:lmaxd+1)
    double complex hankws(0:lmaxd+1)
    double complex cden0(irmd,0:lmaxd+1)
    double complex cden1(irmd,0:lmaxd+1)

    integer lmaxd1

    lmaxd1 = lmaxd+1

    ek = sqrt(ez)
    c0ll = 1.d0/sqrt(16.d0*atan(1.d0))
    ciek = (0.d0,1.d0)*ek

    ! initialize arrays ...
    cden0 = zero
    cden1 = zero

    do ir = 2, irmd
      call beshan(hankws,bessjw,bessyw,r(ir)*ek,lmaxd1)
      do l = 0, lmaxd
        pz(ir,l) = bessjw(l)*r(ir)
        qz(ir,l) = (bessyw(l) - ci*bessjw(l))*r(ir)
      enddo ! l
    enddo ! ir

    imt1 = ircut(1)
    do l1 = 0, lmaxd1
      cden0(1,l1) = (0.d0,0.d0)
      cden1(1,l1) = (0.d0,0.d0)
    enddo ! l1

    do ir = 2, irmd
      cden0(ir,0) = ek*pz(ir,0)*qz(ir,0)
      cden1(ir,0) = ek*pz(ir,0)**2*(0.d0,-1.d0)
      cden1(ir,lmaxd1) = ciek*r(ir)**2
    enddo ! ir

    do l1 = 1, lmaxd
      do ir = 2, irmd
        cden0(ir,l1) = ek*pz(ir,l1)*qz(ir,l1)*(2*l1+1)
        cden1(ir,l1) = ek*pz(ir,l1)**2*(0.d0,-1.d0)*(2*l1+1)
      enddo ! ir
    enddo ! l1

    do l1 = 0, lmaxd1
      if (ipan > 1) then
        do ir = imt1 + 1, irmd
          cden0(ir,l1) = cden0(ir,l1)*thetas(ir-imt1,1)*c0ll
          cden1(ir,l1) = cden1(ir,l1)*thetas(ir-imt1,1)*c0ll
        enddo ! ir
      endif ! ipan > 1
      dos0 = dos0 + wez*simpson(cden0(1:,l1), ipan, ircut, drdi)
      dos1 = dos1 + wez*simpson(cden1(1:,l1), ipan, ircut, drdi)
    enddo ! l1

  endsubroutine ! rhoval0

endmodule ValenceDensity_mod