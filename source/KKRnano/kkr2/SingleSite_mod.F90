module SingleSite_mod
  implicit none
  
  private
  public :: calcdtmat, calcdtmat_deltaez, calctmat, regns


  double complex, parameter :: cone=(1.d0,0.d0), zero=(0.d0,0.d0), ci=(0.d0,1.d0)
  double precision, parameter :: pi=4.d0*atan(1.d0)

  contains
  
  subroutine calcdtmat(ldau,nldau,icst, nsra,ez,dz, drdi,r,vins, visp,zat,ipan, ircut,cleb,loflm,icleb,iend, dtde,tr_alph,lmax, lldau,wmldau_ispin, ncleb, ipand, irmd, irnsd, method)
    integer, intent(in) :: ncleb, ipand, irmd, irnsd, lmax,iend
    double complex, intent(in) :: dz
    double complex, intent(in) :: ez
    double precision, intent(in) :: cleb(ncleb,2)
    double precision, intent(in) :: vins((irmd-irnsd):irmd,(2*lmax+1)**2), visp(irmd), wmldau_ispin(2*lmax+1, 2*lmax+1, lmax + 1)
    double complex, intent(out) :: tr_alph
    double complex, intent(out) :: dtde((lmax+1)**2,(lmax+1)**2)
    double precision, intent(in) :: drdi(irmd)
    double precision, intent(in) :: r(irmd)
    double precision, intent(in) :: zat
    integer, intent(in) :: ipan,nldau
    integer, intent(in) :: ircut(0:ipand),lldau(lmax + 1)
    integer, intent(in) :: icleb(ncleb,3),loflm((2*lmax+1)**2)
    integer, intent(in) :: icst,nsra
    logical, intent(in) :: ldau
    integer, intent(in) :: method
    
    double complex :: tmatn1((lmax+1)**2,(lmax+1)**2), tmatn2((lmax+1)**2,(lmax+1)**2)
    double complex :: tr_alph1, tr_alph2
    double complex :: ez1, ez2

    ! interpolation points for difference quotient
    ez1 = ez + dz !.. perpendicular to the contour
    ez2 = ez - dz !.. parallel to the contour

    call calctmat(ldau,nldau,icst, nsra,ez1, drdi,r,vins, visp,zat,ipan, ircut,cleb,loflm,icleb,iend, tmatn1,tr_alph1,lmax, lldau,wmldau_ispin, ncleb, ipand, irmd, irnsd, method)
    call calctmat(ldau,nldau,icst, nsra,ez2, drdi,r,vins, visp,zat,ipan, ircut,cleb,loflm,icleb,iend, tmatn2,tr_alph2,lmax, lldau,wmldau_ispin, ncleb, ipand, irmd, irnsd, method)

    ! dt(e)
    ! -----
    !  de
    dtde(:,:) = (tmatn1(:,:) - tmatn2(:,:))*(0.5d0/dz)
    tr_alph = -(tr_alph1 - tr_alph2)*(0.5d0/dz)
  
  endsubroutine calcdtmat

  !------------------------------------------------------------------------------
  !> get \delta e(z) for difference quotient to calculate the derivative of
  !> \delta t_ref.
  !> @param[out] dz delta e(z)
  subroutine calcdtmat_deltaez(dz, ie, npnt1, npnt2, npnt3, tk)
    double complex, intent(out) :: dz
    integer, intent(in) :: ie, npnt1, npnt2, npnt3
    double precision, intent(in) :: tk
    
    double precision, parameter :: kb=6.333659d-6

    dz = zero

    if (ie <= npnt1 .or. ie > (npnt1+npnt2+npnt3)) then
      dz = dcmplx(0.01d0*pi*kb*tk, 0.d0)
    else
      dz = dcmplx(0.d0, 0.01d0*pi*kb*tk)
    endif
    
  endsubroutine calcdtmat_deltaez


  subroutine calctmat(ldau,nldau,icst, nsra,ez, drdi,r,vins,visp,zat,ipan, ircut,cleb,loflm,icleb,iend, tmatn,tr_alph,lmax, lldau,wmldau_ispin, ncleb, ipand, irmd, irnsd, method)
    use SingleSiteHelpers_mod, only: wfmesh, cradwf
    integer, intent(in) :: ncleb, ipand, irmd, irnsd, icst, iend, ipan, nsra, lmax, nldau
    double precision, intent(in) :: zat
    double complex, intent(in) :: ez
    double complex, intent(out) :: tr_alph
    logical, intent(in) :: ldau
    double complex, intent(out) :: tmatn((lmax+1)**2,(lmax+1)**2)
    double precision, intent(in) :: cleb(ncleb,2), drdi(irmd), r(irmd)
    double precision, intent(in) :: vins((irmd-irnsd):irmd,(2*lmax+1)**2)
    double precision, intent(in) :: visp(irmd)
    double precision, intent(in) :: wmldau_ispin(2*lmax+1,2*lmax+1,lmax+1)
    integer, intent(in) :: icleb(ncleb,3), ircut(0:ipand)
    integer, intent(in) :: loflm((2*lmax+1)**2)
    integer, intent(in) :: lldau(lmax+1)
    integer, intent(in) :: method

    double precision, parameter :: cvlight = 274.0720442d0
    double complex   :: eryd, ek, det
    integer          :: l, lmlo, lmhi, mmax, im, ildau
    double complex   :: alpha(0:lmax), fz(irmd,0:lmax)
    double complex   :: pns((lmax+1)**2,(lmax+1)**2,(irmd-irnsd):irmd,2)
    double complex   :: pz(irmd,0:lmax), qz(irmd,0:lmax), sz(irmd,0:lmax), tmat(0:lmax)
    double precision :: rs(irmd,0:lmax), s(0:lmax), ldaucut(irmd), wmldauav(lmax+1)
    integer          :: lmmaxd

    lmmaxd = (lmax+1)**2

  !------------------------------------------------------------------------------
  ! initialisation of local variables to be on the safe side
    ek = zero
    det = zero
    alpha = zero
    fz = zero
    pns = zero
    pz = zero
    qz = zero
    sz = zero
    tmat = zero
    rs = 0.d0
    s = 0.d0
    ldaucut = 0.d0
    wmldauav = 0.d0
  ! -----------------------------------------------------------------------------

    !ldau

    if (ldau) then

      do ildau = 1, nldau
        wmldauav(ildau) = 0.d0
        lmlo = lldau(ildau)**2 + 1
        lmhi = (lldau(ildau)+1)**2
        mmax = lmhi - lmlo + 1
        do im = 1, mmax
          wmldauav(ildau) = wmldauav(ildau) + wmldau_ispin(im,im,ildau)
        enddo ! im
        wmldauav(ildau) = wmldauav(ildau)/dble(mmax)
      enddo ! ildau
      
      ! -> note: application if wldau makes the potential discontinuous.
      !    a cutoff can be used if necessary to make the potential continuous
      !    for example (array bounds should be adjusted):
      !       if(test('cutoff  ')) then
      !         do ir = 1,irmd
      !           ldaucut(ir) = ( 1.d0 + dexp( 20.d0*(r(ir)-r(349)) ) ) * ( 1.d0 + dexp( 20.d0*(r(276)-r(ir)) ) )
      !           ldaucut(ir) = 1d0/ldaucut(ir)
      !         enddo
      !       else
      !         do ir = 1,irmd
      !           ldaucut(ir) = 1.d0
      !         enddo
      !       endif
    endif

    !ldau
    tmatn(:,:) = zero

    eryd = ez

    call wfmesh(eryd,ek,cvlight,nsra,zat,r,s,rs,ircut(ipan), irmd,lmax)

    call cradwf(eryd,ek,nsra,alpha,ipan,ircut,cvlight,rs,s, pz,fz,qz,sz,tmat,visp,drdi,r,zat, ldau,nldau,lldau,wmldauav,ldaucut, lmax, irmd, ipand)

    call pnstmat(drdi,ek,icst,pz,qz,fz,sz,pns, tmatn, vins,ipan,ircut,nsra,cleb,icleb,iend,loflm, tmat,det,lmax, & 
         ldau,nldau,lldau, wmldau_ispin,wmldauav,ldaucut, lmax, irmd, irnsd, ipand, ncleb, method)

    tr_alph = log(det)
    do l = 0, lmax
      tr_alph = tr_alph + (2.d0*l + 1.d0)*log(alpha(l))
    enddo ! l

  endsubroutine ! calctmat

  
  
   
      
      
      
!>    construct the t-matrix from radial solutions and potential
  subroutine pnstmat(drdi,ek,icst,pz,qz,fz,sz,pns,tmatll,vins, ipan,ircut,nsra,cleb,icleb,iend,loflm,tmat, &
                    det,lkonv, ldau,nldau,lldau, wmldau_ispin,wmldauav,ldaucut, lmaxd, irmd, irnsd, ipand, ncleb, method)
    use SingleSiteHelpers_mod, only: vllns, wftsca
    integer, intent(in) :: lmaxd, irmd, irnsd, ipand, ncleb
    double complex, intent(in) :: ek
    double complex, intent(out) :: det
    integer, intent(in) :: icst,iend,ipan,lkonv,nsra,nldau
    logical, intent(in) :: ldau
    double complex, intent(in) :: fz(irmd,0:lmaxd)
    double complex, intent(inout) :: pns((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd,2)
    double complex, intent(in) :: pz(irmd,0:lmaxd)
    double complex, intent(in) :: qz(irmd,0:lmaxd)
    double complex, intent(in) :: sz(irmd,0:lmaxd)
    double complex, intent(in) :: tmat(0:lmaxd)
    double complex, intent(out) :: tmatll((lmaxd+1)**2,(lmaxd+1)**2)
    double precision, intent(in) :: cleb(ncleb,2)
    double precision, intent(in) :: drdi(irmd)
    double precision, intent(in) :: vins(irmd-irnsd:irmd,(2*lmaxd+1)**2)
    double precision, intent(in) :: wmldau_ispin(2*lmaxd+1,2*lmaxd+1,lmaxd+1)
    double precision, intent(in) :: ldaucut(irmd)
    double precision, intent(in) :: wmldauav(lmaxd+1)
    integer, intent(in) :: icleb(ncleb,3), ircut(0:ipand), loflm(*), lldau(lmaxd+1)
    integer, intent(in) :: method !< method for single site solver, Volterra or Fredholm
    
    external :: zgetrf ! from BLAS
    double complex  :: ar((lmaxd+1)**2,(lmaxd+1)**2)
    double complex  :: cmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
    double complex  :: dmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
    double complex  :: efac((lmaxd+1)**2)
    double complex  :: pzekdr((lmaxd+1)**2,irmd-irnsd:irmd,2)
    double complex  :: pzlm((lmaxd+1)**2, irmd-irnsd:irmd,2)
    double complex  :: qzekdr((lmaxd+1)**2,irmd-irnsd:irmd,2)
    double complex  :: qzlm((lmaxd+1)**2,irmd-irnsd:irmd,2)
    double precision :: vnspll((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)

    integer :: i,irc1,lm1,lm2,lmmkonv,ir,m1,m2,lmlo,lmhi,ildau
    integer :: ipvt((lmaxd+1)**2)
    integer :: info, irmind, lmmaxd

!     initialisation of some local arrays
    pzlm = zero
    qzlm = zero
    pzekdr = zero
    qzekdr = zero
    ar = zero
!     end initialisation

    irmind = irmd-irnsd
    lmmaxd = (lmaxd+1)**2

    irc1 = ircut(ipan)

    call vllns(vnspll,vins,cleb,icleb,iend, lmaxd, irmd, irnsd, ncleb)

    if (lkonv /= lmaxd) then
      lmmkonv = (lkonv+1)**2
      do lm1 = 1, lmmaxd
        do lm2 = lmmkonv + 1, lmmaxd
          do i = irmind, irmd
            vnspll(lm2,lm1,i) = 0.d0
            vnspll(lm1,lm2,i) = 0.d0
          enddo
        enddo
      enddo
    else ! lkonv /= lmaxd
      lmmkonv = lmmaxd
    endif


!-----------------------------------------------------------------------
! lda+u
! add wldau to non-spherical porential vins in case of lda+u
! use the average wldau (=wldauav) and calculate the deviation
! of wldau from this. use the deviation in the born series
! for the non-spherical wavefunction, while the average is
! used for the spherical wavefunction.
!
    if (ldau) then
      do ildau = 1, nldau
        if (lldau(ildau) >= 0) then

          lmlo = lldau(ildau)*lldau(ildau) + 1
          lmhi = (lldau(ildau)+1)*(lldau(ildau)+1)
          do ir = irmind, irmd
            do lm2 = lmlo, lmhi
              m2 = lm2 - lmlo + 1
              do lm1 = lmlo, lmhi
                m1 = lm1 - lmlo + 1
                vnspll(lm1,lm2,ir) = vnspll(lm1,lm2,ir) + wmldau_ispin(m1,m2,ildau)*ldaucut(ir) ! -> first add wldau to all elements ...
              enddo ! lm1
              vnspll(lm2,lm2,ir) = vnspll(lm2,lm2,ir) - wmldauav(ildau) * ldaucut(ir) ! ... and then subtract average from diag. elements
            enddo ! lm2
          enddo ! ir
          
        endif ! (lldau(ildau) >= 0
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
!---> determine the regular non sph. wavefunction
!
    call regns(ar,tmatll,efac,pns,vnspll,icst,ipan,ircut,pzlm,qzlm, pzekdr,qzekdr,ek,pns(1,1,irmind,1),cmat, pns(1,1,irmind,2),dmat,nsra,irmind,irmd,ipand,lmmaxd, method)

    do lm1 = 1, lmmkonv
      tmatll(lm1,lm1) = tmatll(lm1,lm1) + tmat(loflm(lm1))
    enddo ! lm1
    
    call zgetrf(lmmaxd,lmmaxd,ar,lmmaxd,ipvt,info)
    
    det = cone
    do lm1 = 1, lmmaxd
      if (ipvt(lm1) /= lm1) det = -det
      det = ar(lm1,lm1)*det
    enddo ! lm1

  endsubroutine ! pnstmat
      
      
!-----------------------------------------------------------------------
!     determines the regular non spherical wavefunctions , the
!       alpha matrix and the t - matrix in the n-th. born appro-
!       ximation ( n given by input parameter icst )
!
!
!     using the wave functions pz and qz ( regular and irregular
!       solution ) of the spherically averaged potential , the
!       regular wavefunction pns is determined by
!
!           pns(ir,lm1,lm2) = ar(ir,lm1,lm2)*pz(ir,l1)
!                                   + br(ir,lm1,lm2)*qz(ir,l1)
!
!      the matrices ar and br are determined by integral equations
!        containing pns and only the non spherical contributions of
!        the potential , stored in vinspll . these integral equations
!        are  solved iteratively with born approximation up to given n.
!
!     the original way of writing the cr and dr matrices in the equa-
!        tions above caused numerical troubles . therefore here are used
!        rescaled ar and br matrices :
!
!              ~
!              ar(ir,lm1,lm2) = sqrt(e)**(l1-l2) * ar(ir,lm1,lm2)*((2*l2-1)!!/(2*l1-1)!!)
!
!              ~
!              br(ir,lm1,lm2) = sqrt(e)**(-l1-l2) * br(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
!
!     for lloyd's formular is only the determinant of the alpha -
!        matrix is needed which is identical with the determinant
!        of the rescaled ar - matrix at the innerst point .
!
!     the non spherical t - matrix is the br matrix at r(irc)
!
!     modified for the use of shape functions
!
!                              (see notes by b.drittler)
!
!                                b.drittler   mar.  1989
!-----------------------------------------------------------------------
!     modified by r. zeller      aug. 1994
!-----------------------------------------------------------------------
!     added volterra equation by m. ogura      jan. 2006
!     Fredholm: true -> use Fredholm equation
!               false -> volterra equation
!-----------------------------------------------------------------------
  subroutine regns(ar, br, efac, pns, vnspll, icst, ipan, ircut, pzlm,  &
                  qzlm, pzekdr, qzekdr, ek, ader, amat, bder, bmat, nsra,  &
                  irmind, irmd, ipand, lmmaxd, method)
    use SingleSiteHelpers_mod, only: csinwd, csout, wfint, wfint0, zgeinv1
    double complex, intent(in) :: ek
    integer, intent(in) ::  icst,ipan,ipand,irmd,irmind,lmmaxd,nsra
    double complex, intent(out) :: ader(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(out) :: ar(lmmaxd,lmmaxd)
    double complex, intent(out) :: amat(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(out) :: bder(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(out) :: bmat(lmmaxd,lmmaxd,irmind:irmd)
    double complex, intent(out) :: br(lmmaxd,lmmaxd)
    double complex, intent(in) :: efac(*)
    double complex, intent(inout) :: pns(lmmaxd,lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: pzekdr(lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: pzlm(lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: qzekdr(lmmaxd,irmind:irmd,2)
    double complex, intent(in) :: qzlm(lmmaxd,irmind:irmd,2)
    double precision, intent(in) ::  vnspll(lmmaxd,lmmaxd,irmind:irmd)
    integer, intent(in) :: ircut(0:ipand)
    integer, intent(in) :: method !< 0: Fredholm (default), 1: Volterra

    external :: zgemm ! from BLAS
    double precision :: err
    integer :: i, ir, irc1, j, lm
    double complex :: pns0(lmmaxd,lmmaxd,irmind:irmd,2), pns1(lmmaxd,lmmaxd,irmind:irmd)
    integer :: ipiv(lmmaxd)
    logical :: Volterra ! false ==> Fredholm equation, true ==> Volterra equation
    
    Volterra = (method /= 0)
    
    irc1 = ircut(ipan)
  
    do i = 0, icst
    
!---> set up integrands for i-th born approximation
      if (i == 0) then
        call wfint0(ader,bder,pzlm,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
      else  ! i
        call wfint(pns,ader,bder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
      endif ! i

!---> call integration subroutines
      if (Volterra) then
        call csout (ader,amat,lmmaxd**2,irmind,irmd,ipan,ircut)
      else ! Fredholm equation
        call csinwd(ader,amat,lmmaxd**2,irmind,irmd,ipan,ircut)
      endif
      call csout(bder,bmat,lmmaxd**2,irmind,irmd,ipan,ircut)
      
      do ir = irmind, irc1
        if (Volterra) amat(:,:,ir) = -amat(:,:,ir)
        do lm = 1, lmmaxd
          amat(lm,lm,ir) = cone + amat(lm,lm,ir)
        enddo ! lm
      enddo ! ir
      
!---> calculate non sph. wft. in i-th born approximation
      do j = 1, nsra
        do ir = irmind, irc1
          do lm = 1, lmmaxd
            pns(:,lm,ir,j) = amat(:,lm,ir)*pzlm(:,ir,j) + bmat(:,lm,ir)*qzlm(:,ir,j)
          enddo ! lm
        enddo ! ir
      enddo ! j
      
      if (i == icst-1) pns0(:,:,irmind:irc1,1:nsra) = pns(:,:,irmind:irc1,1:nsra) ! store a copy of pns in pns0 for the convergence check 

    enddo ! i
  
!-----------------------------------------------------------------------
! check convergence of pns after the last iteration comparing to the previous iteration
    pns0(:,:,irmind:irc1,1:nsra) = pns0(:,:,irmind:irc1,1:nsra) - pns(:,:,irmind:irc1,1:nsra)
      
    do j = 1, nsra
      call csout(pns0(1,1,irmind,j), pns1, lmmaxd**2, irmind, irmd, ipan, ircut)
      err = maxval(abs(pns1(:,:,irc1)))
      ! convergence check
      if (err > 1d-3) then
        if (Volterra) then
          write(*,*)'regns.f: Volterra equation does not converge -> consider &
          increasing icst in input.conf or switching to Fredholm'
        else
          write(*,*)'regns.f: Fredholm equation does not converge -> consider &
          increasing icst in input.conf or switching to Volterra'
        endif
        stop 'error 1 in regns.f'
      endif
    enddo ! j
! end convergence check      
!-----------------------------------------------------------------------
    
  
    if (Volterra) then
!-----------------------------------------------------------------------
! only Volterra equation
      
      call zgeinv1(amat(1,1,irc1), ar, br, ipiv, lmmaxd) ! invert the last amat

      do ir = irmind, irc1
      
        call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,amat(1,1,ir),lmmaxd,ar,lmmaxd,zero,ader(1,1,ir),lmmaxd)
        call zgemm('N','N',lmmaxd,lmmaxd,lmmaxd,cone,bmat(1,1,ir),lmmaxd,ar,lmmaxd,zero,bder(1,1,ir),lmmaxd)

        ! overwrite amat and bmat, keep the values in ader and bder (since these are intent(out)
        amat(:,:,ir) = ader(:,:,ir)
        bmat(:,:,ir) = bder(:,:,ir)
      enddo ! ir

      ! create the final solution pns from amat and bmat
      do j = 1, nsra
        do ir = irmind, irc1
          do lm = 1, lmmaxd
            pns(:,lm,ir,j) = efac(lm)*(amat(:,lm,ir)*pzlm(:,ir,j) + bmat(:,lm,ir)*qzlm(:,ir,j))
          enddo ! lm
        enddo ! ir
      enddo ! j
      
!-----------------------------------------------------------------------
    else  ! Volterra
    
!---> rescale with efac (Fredholm only)
      do j = 1, nsra
        do ir = irmind, irc1
          do lm = 1, lmmaxd
            pns(:,lm,ir,j) = pns(:,lm,ir,j)*efac(lm)
          enddo ! lm
        enddo ! ir
      enddo ! j
    
    endif ! Volterra

!---> store alpha and t - matrix
    do lm = 1, lmmaxd
      ar(:,lm) = amat(:,lm,irmind)
      br(:,lm) = bmat(:,lm,irc1)*efac(1:lmmaxd)*efac(lm)/ek !---> t-matrix
    enddo ! lm
      

  endsubroutine ! regns

endmodule SingleSite_mod
