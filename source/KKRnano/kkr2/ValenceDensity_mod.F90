module ValenceDensity_mod
  implicit none
  private
  public :: rhoval, rhoval0
  
  double complex, parameter :: cone=(1.d0,0.d0), zero=(0.d0,0.d0)
  
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
  !     modified by R. Zeller      Aug. 1994
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
    endif ! first Born iteration

    !---> call integration subroutines
    call csinwd(cder,cmat,lmmaxd**2,irmind,irmd,ipan,ircut)
    call csinwd(dder,dmat,lmmaxd**2,irmind,irmd,ipan,ircut)
    do ir = irmind, irc1
      do lm = 1, lmmaxd
        dmat(lm,lm,ir) = dmat(lm,lm,ir) - CONE
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

  !---> store c - and d - matrix
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

  endsubroutine
  
  
  
subroutine pnsqns(ar,cr,dr,drdi,ek,icst,pz,qz,fz,sz,pns,qns,nsra, &
vins,ipan,ircut,cleb,icleb,iend,loflm,lkonv, &
ispin,ldau,nldau,lldau, &
wmldau,wmldauav,ldaucut, &
lmaxd, nspind, irmd, irnsd, ipand, ncleb)
  use SingleSite_mod, only: regns
  use SingleSiteHelpers_mod, only: vllns, wftsca
  

  integer lmaxd
  integer nspind
  integer irmd
  integer irnsd
  integer ipand
  integer ncleb
  !     ..
  !     .. Scalar Arguments ..
  double complex     ek
  integer            icst,iend,ipan,lkonv,nsra,nldau,ispin
  logical            ldau
  !     ..
  !     .. Array Arguments ..

  double complex     ar((lmaxd+1)**2,(lmaxd+1)**2)
  double complex     cr((lmaxd+1)**2,(lmaxd+1)**2)
  double complex     fz(irmd,0:lmaxd)

  !     DOUBLE COMPLEX     PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
  double complex     pns((lmaxd+1)**2,(lmaxd+1)**2, &
  (irmd-irnsd):irmd,2)

  double complex     pz(irmd,0:lmaxd)
  !     DOUBLE COMPLEX     QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
  double complex     qns((lmaxd+1)**2,(lmaxd+1)**2, &
  (irmd-irnsd):irmd,2)

  double complex     qz(irmd,0:lmaxd)
  double complex     sz(irmd,0:lmaxd)

  double precision   cleb(ncleb,2)
  double precision   drdi(irmd)
  !     DOUBLE PRECISION   VINS(IRMIND:IRMD,LMPOTD)
  double precision   vins(irmd-irnsd:irmd,(2*lmaxd+1)**2)
  !     DOUBLE PRECISION   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  double precision   wmldau(2*lmaxd+1,2*lmaxd+1,lmaxd+1,nspind)
  double precision   ldaucut(irmd)

  integer            icleb(ncleb,3),ircut(0:ipand),loflm(*)
  integer            lldau(lmaxd+1)
  !     ..
  !     .. Local Scalars ..
  integer            i,ir,irc1,lm1,lm2,lmmkonv,m1,m2, &
  lmlo,lmhi,ildau
  !     ..
  !     .. Local Arrays ..
  !     DOUBLE COMPLEX     CMAT(LMMAXD,LMMAXD,IRMIND:IRMD)
  double complex     cmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
  !     DOUBLE COMPLEX     DMAT(LMMAXD,LMMAXD,IRMIND:IRMD)
  double complex     dmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
  !     DOUBLE COMPLEX     DR(LMMAXD,LMMAXD)
  double complex     dr((lmaxd+1)**2, (lmaxd+1)**2)
  !     DOUBLE COMPLEX     EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2)
  double complex     efac((lmaxd+1)**2)
  double complex     pzekdr((lmaxd+1)**2,irmd-irnsd:irmd,2)
  !     DOUBLE COMPLEX     PZLM(LMMAXD,IRMIND:IRMD,2)
  double complex     pzlm((lmaxd+1)**2,(irmd-irnsd):irmd,2)
  !     DOUBLE COMPLEX     QZEKDR(LMMAXD,IRMIND:IRMD,2)
  double complex     qzekdr((lmaxd+1)**2,(irmd-irnsd):irmd,2)
  !     DOUBLE COMPLEX     QZLM(LMMAXD,IRMIND:IRMD,2)
  double complex     qzlm((lmaxd+1)**2, (irmd-irnsd):irmd,2)
  !     DOUBLE COMPLEX     TMATLL(LMMAXD,LMMAXD)
  double complex     tmatll((lmaxd+1)**2, (lmaxd+1)**2)
  !     DOUBLE PRECISION   VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
  double precision   vnspll((lmaxd+1)**2, (lmaxd+1)**2, &
  irmd-irnsd:irmd)

  double precision   wmldauav(lmaxd+1)
  

  integer             lmmaxd
  integer             irmind

  irmind= irmd-irnsd
  lmmaxd= (lmaxd+1)**2

  irc1 = ircut(ipan)
  !

  ! calculate V_LL' potential
  ! V_{LL'} = \sum_{L''} C_{L L' L''} V_{L''}
  call vllns(vnspll,vins,cleb,icleb,iend, &
  lmaxd, irmd, irnsd, ncleb)


  if (lkonv /= lmaxd) then
    lmmkonv = (lkonv+1)* (lkonv+1)
    do lm1 = 1,lmmaxd
      do lm2 = lmmkonv + 1,lmmaxd
        do i = irmind,irmd
          vnspll(lm2,lm1,i) = 0.0d0
          vnspll(lm1,lm2,i) = 0.0d0
        end do
      end do
    end do
  else
    lmmkonv = lmmaxd
  end if

  !-----------------------------------------------------------------------
  ! LDA+U
  ! Add WLDAU to non-spherical porential VINS in case of LDA+U
  ! Use the average wldau (=wldauav) and the deviation
  ! of wldau from this. Use the deviation in the Born series
  ! for the non-spherical wavefunction, while the average is
  ! used for the spherical wavefunction.
  !
  if (ldau) then
    do ildau=1,nldau
      if (lldau(ildau) >= 0) then
        !
        lmlo = lldau(ildau)*lldau(ildau) + 1
        lmhi = (lldau(ildau)+1)*(lldau(ildau)+1)
        !
        do ir = irmind,irmd
          !
          ! -> First add wldau to all elements ...
          !
          do lm2 = lmlo,lmhi
            m2 = lm2 - lmlo + 1
            do lm1 = lmlo,lmhi
              m1 = lm1 - lmlo + 1
              vnspll(lm1,lm2,ir) =vnspll(lm1,lm2,ir) &
              +wmldau(m1,m2,ildau,ispin)*ldaucut(ir)
            enddo
            !
            ! ... and then subtract average from diag. elements
            !
            vnspll(lm2,lm2,ir) =  vnspll(lm2,lm2,ir) &
            - wmldauav(ildau) * ldaucut(ir)
          enddo
        enddo
      endif
    enddo
  endif
  !
  ! LDA+U
  !-----------------------------------------------------------------------


  !
  !---> get wfts of same magnitude by scaling with efac
  !
  call wftsca(drdi,efac,pz,qz,fz,sz,nsra,pzlm,qzlm,pzekdr,qzekdr, &
  ek,loflm,irmind,irmd,lmaxd,lmmaxd)
  !
  !---> determine the irregular non sph. wavefunction
  !
  call irwns(cr,dr,efac,qns,vnspll,icst,ipan,ircut,nsra,pzlm,qzlm, &
  pzekdr,qzekdr,qns(1,1,irmind,1),cmat, &
  qns(1,1,irmind,2),dmat,irmind,irmd,ipand,lmmaxd)
  !
  !---> determine the regular non sph. wavefunction
  !
  call regns(ar,tmatll,efac,pns,vnspll,icst,ipan,ircut,pzlm,qzlm, &
  pzekdr,qzekdr,ek,pns(1,1,irmind,1),cmat, &
  pns(1,1,irmind,2),dmat,nsra,irmind,irmd,ipand,lmmaxd)
  !

endsubroutine
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

  double complex, external :: zdotu ! from BLAS
  double complex, parameter :: zero = (0.d0, 0.d0)
  double complex :: efac1,efac2,ffz,gmatl,ppz,v1,v2
  double precision :: c0ll
  integer :: ir,j,j0,j1,l,l1,l2,lm1,lm2,lm3,lm3max,ln2,ln3,m
  double complex :: vr((lmax+1)**2,(lmax+1)**2)
  double complex :: wf(irmd,0:lmax,0:lmax)
  double complex :: wr((lmax+1)**2,(lmax+1)**2)
  integer :: lmmaxd, lmaxd
  
  lmaxd = lmax
  lmmaxd= (lmaxd+1)**2


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
  end if
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
    end if

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
        end if

      enddo ! l2
    enddo ! l1
  enddo ! lm3

endsubroutine ! rhoin




subroutine rhons(den,df,drdi,gmat,ek,rho2ns,ipan,ircut,thetas, &
ifunm,lmsp,nsra,qns,pns,ar,cr,pz,fz,qz,sz,cleb, &
icleb,jend,iend,ekl, &
lmax, irmd, irnsd, irid, ipand, nfund, ncleb)
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
  

  integer lmax
  integer irmd
  integer ncleb
  integer irnsd
  integer irid
  integer ipand
  integer nfund

  !     INTEGER IRMIND
  !     PARAMETER (IRMIND=IRMD-IRNSD)
  !     INTEGER LMPOTD,LMMAXD
  !     PARAMETER (LMPOTD= (LPOTD+1)**2) ! LMPOTD= (2*LMAX+1)**2
  !     PARAMETER (LMMAXD= (LMAXD+1)**2)
  !     INTEGER LMAXD1
  !     PARAMETER (LMAXD1= LMAXD+1)
  !     ..
  !     .. Scalar Arguments ..
  double complex df,ek
  integer iend,ipan,nsra
  !     ..
  !     .. Array Arguments ..
  !     DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),DEN(0:LMAXD1),
  !    +               EKL(0:LMAXD),FZ(IRMD,0:LMAXD),GMAT(LMMAXD,LMMAXD),
  !    +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),PZ(IRMD,0:LMAXD),
  !    +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),QZ(IRMD,0:LMAXD),
  !    +               SZ(IRMD,0:LMAXD)
  !     DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
  !    +                 THETAS(IRID,NFUND)
  !     INTEGER ICLEB(NCLEB,3),IFUNM(*),IRCUT(0:IPAND),
  !    +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)

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

  !     DOUBLE PRECISION CLEB(*),DRDI(IRMD),RHO2NS(IRMD,LMPOTD),
  !    +                 THETAS(IRID,NFUND)
  !     INTEGER ICLEB(NCLEB,3),IFUNM(*),IRCUT(0:IPAND),
  !    +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(*)

  double precision cleb(*)
  double precision drdi(irmd)
  double precision rho2ns(irmd,(2*lmax+1)**2)
  double precision thetas(irid,nfund)

  integer icleb(ncleb,3)
  integer ifunm(*)
  integer ircut(0:ipand)
  integer jend((2*lmax+1)**2,0:lmax,0:lmax)
  integer lmsp(*)
  
  external :: csimpk
  double complex denns,v1
  integer imt1,l,lm,m
  double complex cden(irmd,0:lmax)
  double complex cdenns(irmd)
  double complex efac((lmax+1)**2)
  
  integer irmind
  irmind=irmd-irnsd
  !
  !---> set up efac(lm) = sqrt(e))**l/(2l - 1)!!
  !
  efac(1) = 1.0d0
  v1 = 1.0d0
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

   call rhoout(cden,df,gmat,ek,pns,qns,rho2ns,thetas,ifunm,ipan, &
   imt1,lmsp,cdenns,nsra,cleb,icleb,iend, &
   lmax, irmd, irnsd, irid, nfund, ncleb)

   !
   call rhoin(ar,cden,cr,df,gmat,ek,rho2ns,irmind,nsra,efac,pz,fz, &
   qz,sz,cleb,icleb,jend,iend,ekl, &
   lmax, irmd, ncleb)

   !
   !---> calculate complex density of states
   !
   do l = 0,lmax
     !
     !---> call integration subroutine
     !
     call csimpk(cden(1,l),den(l),ipan,ircut,drdi)
   enddo ! l

   if (ipan > 1) then
     call csimpk(cdenns,denns,ipan,ircut,drdi)
     den(lmax+1) = denns
   end if

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

  external :: zgemm ! from BLAS
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

    call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,pns(:,1,ir,1),lmmaxd,gmat,lmmaxd,ek,qnsi,lmmaxd)
    call zgemm('N','T',lmmaxd,lmmaxd,lmmaxd,cone,pns(:,1,ir,1),lmmaxd,qnsi,lmmaxd,zero,wr(:,1,ir),lmmaxd)

    if (nsra == 2) then
      qnsi(1:lmmaxd,1:lmmaxd) = qns(1:lmmaxd,1:lmmaxd,ir,2)

      call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,pns(:,1,ir,2),lmmaxd,gmat,lmmaxd,ek,qnsi,lmmaxd)
      call zgemm('N','T',lmmaxd,lmmaxd,lmmaxd,cone,pns(:,1,ir,2),lmmaxd,qnsi,lmmaxd,cone,wr(:,1,ir),lmmaxd)

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

end subroutine rhoout



subroutine RHOVAL(LDORHOEF,ICST,IELAST,NSRA, &
                  ISPIN,NSPIN, &
                  EZ,WEZ,DRDI,R,IRMIN, &
                  VINS,VISP,ZAT,IPAN,IRCUT, &
                  THETAS,IFUNM,LMSP,RHO2NS,R2NEF,DEN, &
                  ESPV,CLEB,LOFLM,ICLEB,IEND,JEND, &
                  GMATN, &                                 ! input
                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &       ! input
                  DMATLDAU, &                              ! output
                  ! new parameters after inc.p removal
                  iemxd, &
                  lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)
  use SingleSiteHelpers_mod, only: CRADWF, WFMESH

  integer :: iemxd
  integer :: lmaxd
  integer :: irmd
  integer :: ncleb
  integer :: irnsd
  integer :: irid
  integer :: ipand
  integer :: nfund

  double precision, parameter :: CVLIGHT=274.0720442D0
  !     ..
  !     .. Scalar Arguments ..
  double precision ::   ZAT
  integer ::            ICST,IELAST,IEND,IPAN,ISPIN,NSPIN,NSRA, &
  IRMIN,NLDAU
  logical ::            LDORHOEF,LDAU
  !     ..
  !     .. Array Arguments ..
  !     DOUBLE COMPLEX     DEN(0:LMAXD1,IEMXD),EZ(IEMXD),
  !    +                   WEZ(IEMXD),
  !    +                   PHILDAU(IRMD,LMAXD1),
  !    +                   DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  !     DOUBLE PRECISION   CLEB(NCLEB,2),DRDI(IRMD),
  !    +                   ESPV(0:LMAXD1,1),
  !    +                   R(IRMD),RHO2NS(IRMD,LMPOTD,2),
  !    +                   R2NEF(IRMD,LMPOTD,2),   ! at fermi energy
  !    +                   THETAS(IRID,NFUND),VINS(IRMIND:IRMD,LMPOTD),
  !    +                   VISP(IRMD),
  !    +                   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  !     INTEGER            ICLEB(NCLEB,3),IFUNM(LMXSPD),IRCUT(0:IPAND),
  !    +                   JEND(LMPOTD,0:LMAXD,0:LMAXD),
  !    +                   LMSP(LMXSPD),LOFLM(LM2D),
  !    +                   LLDAU(LMAXD1)

  double complex     DEN(0:LMAXD+1,IEMXD)
  double complex     EZ(IEMXD)
  double complex     WEZ(IEMXD)
  double complex     PHILDAU(IRMD,LMAXD+1)
  double complex     DMATLDAU(2*LMAXD+1,2*LMAXD+1,NSPIN,LMAXD+1)

  double complex    GMATN((LMAXD+1)**2,(LMAXD+1)**2,IEMXD,NSPIN)

  double precision ::   CLEB(NCLEB,2)
  double precision ::   DRDI(IRMD)
  double precision ::   ESPV(0:LMAXD+1,1)
  double precision ::   R(IRMD)
  double precision ::   RHO2NS(IRMD,(2*LMAXD+1)**2,2)
  double precision ::   R2NEF(IRMD,(2*LMAXD+1)**2,2)
  double precision ::   THETAS(IRID,NFUND)
  double precision ::   VINS(IRMD-IRNSD:IRMD,(2*LMAXD+1)**2)
  double precision ::   VISP(IRMD)
  double precision ::   WMLDAU(2*LMAXD+1,2*LMAXD+1,LMAXD+1,NSPIN)

  integer ::            ICLEB(NCLEB,3)
  integer ::            IFUNM((4*LMAXD+1)**2)
  integer ::            IRCUT(0:IPAND)
  integer ::            JEND((2*LMAXD+1)**2,0:LMAXD,0:LMAXD)
  integer ::            LMSP((4*LMAXD+1)**2)
  integer ::            LOFLM((2*LMAXD+1)**2)
  integer ::            LLDAU(LMAXD+1)

  !     .. Local Scalars ..
  double complex     DF,ERYD,EK
  integer ::            IDIM,IE,IR,L,LM1,LM2, &
  LMLO,LMHI,MMAX,IM,ILDAU
  !     ..
  !     .. Local Arrays ..
  !     DOUBLE COMPLEX     ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),
  !    +                   DR(LMMAXD,LMMAXD),
  !    +                   CR(LMMAXD,LMMAXD),
  !    +                   EKL(0:LMAXD),FZ(IRMD,0:LMAXD),
  !    +                   GMATLL(LMMAXD,LMMAXD),
  !    +                   GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND),
  !    +                   PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
  !    +                   PZ(IRMD,0:LMAXD),
  !    +                   QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
  !    +                   QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),
  !    +                   TMAT(0:LMAXD)
  !     DOUBLE PRECISION   RS(IRMD,0:LMAXD),S(0:LMAXD),
  !    +                   LDAUCUT(IRMD),
  !    +                   WMLDAUAV(LMAXD1)
  !     DOUBLE COMPLEX     DENDUM(0:LMAXD1)

  ! The following arrays are local
  double complex    ALPHA(0:LMAXD)
  double complex    AR((LMAXD+1)**2,(LMAXD+1)**2)
  double complex    DR((LMAXD+1)**2,(LMAXD+1)**2)
  double complex    CR((LMAXD+1)**2,(LMAXD+1)**2)
  double complex    EKL(0:LMAXD)
  double complex    FZ(IRMD,0:LMAXD)
  double complex    GMATLL((LMAXD+1)**2,(LMAXD+1)**2)

  double complex    QZ(IRMD,0:LMAXD)
  double complex    SZ(IRMD,0:LMAXD)
  double complex    TMAT(0:LMAXD)
  double complex    PZ(IRMD,0:LMAXD)

  double precision ::   RS(IRMD,0:LMAXD)
  double precision ::   S(0:LMAXD)
  double precision ::   LDAUCUT(IRMD)
  double precision ::   WMLDAUAV(LMAXD+1)

  double complex     DENDUM(0:LMAXD+1)

  ! dynamically allocate large arrays
  ! DOUBLE COMPLEX    PNS((LMAXD+1)**2,(LMAXD+1)**2,IRMD-IRNSD:IRMD,2)
  ! DOUBLE COMPLEX    QNS((LMAXD+1)**2,(LMAXD+1)**2,IRMD-IRNSD:IRMD,2)
  double complex, dimension(:,:,:,:), allocatable :: PNS
  double complex, dimension(:,:,:,:), allocatable :: QNS

  external :: DAXPY,DSCAL ! from BLAS

  integer :: memory_stat
  logical :: memory_fail

  integer ::             LMMAXD
  integer ::             LMAXD1
  integer ::             NSPIND
  integer ::             LMPOTD

  LMPOTD = (2*LMAXD+1)**2
  LMMAXD= (LMAXD+1)**2
  LMAXD1= LMAXD+1
  NSPIND = NSPIN

  !------------------- Array allocations ---------------------------------
  memory_stat = 0
  memory_fail = .false.

  allocate(PNS(LMMAXD,LMMAXD,IRMD-IRNSD:IRMD,2), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(QNS(LMMAXD,LMMAXD,IRMD-IRNSD:IRMD,2), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (memory_fail .eqv. .true.) then
    write(*,*) "RHOVAL: FATAL Error, failure to allocate memory."
    write(*,*) "        Probably out of memory."
    stop
  end if

  !-----------------------------------------------------------------------
  ! Initialise local variables to be on the safe side
  EK = zero
  PNS = zero
  QNS = zero
  ALPHA = zero
  AR = zero
  DR = zero
  CR = zero
  EKL = zero
  FZ = zero
  !GMATLL = zero initialised further down
  QZ = zero
  SZ = zero
  TMAT = zero
  PZ = zero
  RS = 0.0d0
  S = 0.0d0
  LDAUCUT = 0.0d0
  WMLDAUAV = 0.0d0
  DENDUM = zero

  !-----------------------------------------------------------------------
  ! LDAU

  if (LDAU) then
    do ILDAU=1,NLDAU
      WMLDAUAV(ILDAU) = 0.D0
      LMLO = LLDAU(ILDAU)*LLDAU(ILDAU) + 1
      LMHI = (LLDAU(ILDAU)+1)*(LLDAU(ILDAU)+1)
      MMAX = LMHI - LMLO + 1
      do IM = 1,MMAX
        WMLDAUAV(ILDAU)=WMLDAUAV(ILDAU)+WMLDAU(IM,IM,ILDAU,ISPIN)
      enddo
      WMLDAUAV(ILDAU) = WMLDAUAV(ILDAU)/DBLE(MMAX)
    enddo
    
    ! -> Note: Application if WLDAU makes the potential discontinuous.
    !    A cutoff can be used if necessary to make the potential continuous
    !    for example (array bounds should be adjusted):
    
!    if(TEST('CUTOFF  ')) then
!      do IR = 1,IRMD
!        LDAUCUT(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) * &
!        ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
!        LDAUCUT(IR) = 1D0/LDAUCUT(IR)
!      enddo
!    else
!      do IR = 1,IRMD
!        LDAUCUT(IR) = 1.D0
!      enddo
!    endif
  endif

  ! LDAU
  !-----------------------------------------------------------------------



  do LM1 = 1,LMPOTD
    do IR = 1,IRMD
      RHO2NS(IR,LM1,ISPIN) = 0.0D0
      R2NEF(IR,LM1,ISPIN) = 0.0D0
    end do
  end do

  ESPV = 0.0D0

  do IE = 1,IELAST

    do LM2 = 1,LMMAXD
      do LM1 = 1,LMMAXD
        GMATLL(LM1,LM2) = GMATN(LM1,LM2,IE,ISPIN)
      end do
    end do

    ERYD = EZ(IE)
    DF = WEZ(IE)/DBLE(NSPIN)
    
    !=======================================================================
    call WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,R,S,RS,IRCUT(IPAN),IRMD,LMAXD)

    call CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,S, &
                PZ,FZ,QZ,SZ,TMAT,VISP,DRDI,R,ZAT, &
                LDAU,NLDAU,LLDAU,WMLDAUAV,LDAUCUT, &
                lmaxd, irmd, ipand)
    !-----------------------------------------------------------------------
    ! non-spherical
    
    call PNSQNS(AR,CR,DR,DRDI,EK,ICST,PZ,QZ,FZ,SZ, &
                PNS,QNS,NSRA,VINS,IPAN,IRCUT, &
                CLEB,ICLEB,IEND,LOFLM,LMAXD,ISPIN, &
                LDAU,NLDAU,LLDAU, &
                WMLDAU,WMLDAUAV,LDAUCUT, &
                lmaxd, nspind, irmd, irnsd, ipand, ncleb)


    do L = 0,LMAXD
      EKL(L) = EK*DBLE(2*L+1)
    end do
    !-----------------------------------------------------------------------
    call RHONS(DEN(0,IE),DF,DRDI,GMATLL,EK, &
               RHO2NS(1,1,ISPIN),IPAN,IRCUT,THETAS,IFUNM,LMSP, &
               NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB, &
               JEND,IEND,EKL, &
               lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)

    !-----------------------------------------------------------------------


    do L = 0,LMAXD1
      ESPV(L,1) = ESPV(L,1) + DIMAG(ERYD*DEN(L,IE)*DF)
    end do
    
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     get the charge at the Fermi energy (IELAST)
    !     call with the energy weight CONE --> not overwrite DF
    !          with the dummy DENDUM       --> not overwrite DEN
    
    if ( (IE == IELAST) .and. (LDORHOEF) ) then
      call RHONS(DENDUM,CONE,DRDI,GMATLL,EK, &
                 R2NEF(1,1,ISPIN),IPAN,IRCUT,THETAS,IFUNM,LMSP, &
                 NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB, &
                 JEND,IEND,EKL, &
                 lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)
    end if

    
    if (LDAU .and. NLDAU >= 1) then
      call LDAUDMAT(DF,PZ,QZ,PNS,QNS,AR,CR,DR,GMATLL, &
                    IPAN,IRCUT,DRDI,EK, &
                    IRMIN,LLDAU,PHILDAU,NLDAU, &
                    DMATLDAU,ISPIN, &
                    lmaxd, nspind, irmd, irnsd, ipand)
        
    endif


  end do

  ! this should really be separated into another routine
  if (ISPIN == 2) then
    IDIM = IRMD*LMPOTD
    call DSCAL(IDIM,2.D0,RHO2NS(1,1,1),1)
    call DAXPY(IDIM,-0.5D0,RHO2NS(1,1,1),1,RHO2NS(1,1,2),1)
    call DAXPY(IDIM,1.0D0,RHO2NS(1,1,2),1,RHO2NS(1,1,1),1)
    
    ! --> do the same at the Fermi energy
    
    call DSCAL(IDIM,2.D0,R2NEF(1,1,1),1)
    call DAXPY(IDIM,-0.5D0,R2NEF(1,1,1),1,R2NEF(1,1,2),1)
    call DAXPY(IDIM,1.0D0,R2NEF(1,1,2),1,R2NEF(1,1,1),1)
  end if

  !------------------- Array deallocations ---------------------------------
  deallocate(PNS)
  deallocate(QNS)
!-----------------------------------------------------------------------

end subroutine RHOVAL


subroutine rhoval0(ez,wez,drdi,r,ipan,ircut, &
thetas,dos0,dos1, &
lmaxd, irmd, irid, ipand, nfund)
  use SingleSiteHelpers_mod, only: beshan
  
  !
  !     .. Parameters ..

  integer lmaxd
  integer irmd
  integer irid
  integer ipand
  integer nfund

  double complex CONE,zero,CI
  parameter ( CONE=(1.d0,0.d0),zero=(0.d0,0.d0),CI=(0.d0,1.d0) )
  !     ..
  !     .. Scalar Arguments ..
  integer ipan
  double complex ez,wez,dos0,dos1
  !     ..
  !     .. Array Arguments ..
  double precision drdi(irmd), &
  r(irmd), &
  thetas(irid,nfund)
  integer ircut(0:ipand)
  !     ..
  !     .. Local Scalars ..
  double complex ek,ciek,denl
  double precision c0ll
  integer ir,l,l1,imt1
  !     ..
  !     .. Local Arrays ..
  double complex pz(irmd,0:lmaxd),qz(irmd,0:lmaxd)
  double complex bessjw(0:lmaxd+1)
  double complex bessyw(0:lmaxd+1)
  double complex hankws(0:lmaxd+1)
  double complex cden0(irmd,0:lmaxd+1)
  double complex cden1(irmd,0:lmaxd+1)

  integer lmaxd1

  lmaxd1 = lmaxd+1
  !     ..
  !
  ek = sqrt(ez)
  c0ll = 1.0d0/sqrt(16.0d0*atan(1.0d0))
  ciek=(0.0d0,1.0d0)*ek
  !
  ! initialize arrays ...
  !
  do ir = 1, irmd
    do l1 = 0, lmaxd1
      cden0(ir,l1) = zero
      cden1(ir,l1) = zero
    enddo
  enddo
  !
  !
  !=======================================================================
  do ir = 2,irmd
    call beshan(hankws,bessjw,bessyw,r(ir)*ek,lmaxd1)
    do l = 0,lmaxd
      pz(ir,l) = bessjw(l)*r(ir)
      qz(ir,l) = (bessyw(l) - CI*bessjw(l))*r(ir)
    enddo
  enddo

  imt1=ircut(1)
  do l1 = 0,lmaxd1
    cden0(1,l1) = (0.0d0,0.0d0)
    cden1(1,l1) = (0.0d0,0.0d0)
  end do

  do ir = 2,irmd
    cden0(ir,0) = ek*pz(ir,0)*qz(ir,0)
    cden1(ir,0) = ek*pz(ir,0)**2*(0.d0,-1.d0)
    cden1(ir,lmaxd1) = ciek*r(ir)**2
  end do

  do l1 = 1,lmaxd
    do ir = 2,irmd
      cden0(ir,l1) = ek*pz(ir,l1)*qz(ir,l1)*(l1+l1+1)
      cden1(ir,l1) = ek*pz(ir,l1)**2*(0.d0,-1.d0)*(l1+l1+1)
    end do  
  end do

  do l1 = 0,lmaxd1
    if (ipan > 1) then
      do ir = imt1 + 1,irmd
        cden0(ir,l1) = cden0(ir,l1)*thetas(ir-imt1,1)*c0ll
        cden1(ir,l1) = cden1(ir,l1)*thetas(ir-imt1,1)*c0ll
      end do
    end if
    call csimpk(cden0(1,l1),denl,ipan,ircut,drdi)
    dos0 = dos0 + wez*denl
    call csimpk(cden1(1,l1),denl,ipan,ircut,drdi)
    dos1 = dos1 + wez*denl
  end do

endsubroutine

endmodule ValenceDensity_mod