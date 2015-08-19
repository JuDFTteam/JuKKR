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
!              ar(ir,lm1,lm2) = sqrt(e)**(l1-l2)
!                             * ar(ir,lm1,lm2)*((2*l2-1)!!/(2*l1-1)!!)
!
!              ~
!              br(ir,lm1,lm2) = sqrt(e)**(-l1-l2)
!                             * br(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
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
!     Fredholm: true -> use Fredholmholm equation
!           false -> volterra equation
!-----------------------------------------------------------------------
      subroutine regns(ar, br, efac, pns, vnspll, icst, ipan, ircut, pzlm,  &
                       qzlm, pzekdr, qzekdr, ek, ader, amat, bder, bmat, nsra,  &
                       irmind, irmd, ipand, lmmaxd)
      implicit none

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

      external :: csinwd, csout, wfint, wfint0, zgeinv1
      double complex efac2!efac1,
      double precision err
      integer i,ir,irc1,j,lm1,lm2,lm3
      double complex pns0(lmmaxd,lmmaxd,irmind:irmd,2), pns1(lmmaxd,lmmaxd,irmind:irmd)
      integer ipiv(lmmaxd)
      double complex, parameter :: cone=(1.0d0,0.0d0), zero=(0.d0,0.d0)
      logical, parameter :: Fredholm = .false.
      
      irc1 = ircut(ipan)
    if (Fredholm) then
!-----------------------------------------------------------------------
! begin Fredholm equation
   
      do i = 0, icst
!---> set up integrands for i-th born approximation
        if (i == 0) then
          call wfint0(ader,bder,pzlm,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
        else  ! i
          call wfint(pns,ader,bder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
        endif ! i
!---> call integration subroutines
        call csinwd(ader,amat,lmmaxd**2,irmind,irmd,ipan,ircut)
        call csout(bder,bmat,lmmaxd**2,irmind,irmd,ipan,ircut)
        do ir = irmind,irc1
          do lm2 = 1,lmmaxd
            amat(lm2,lm2,ir) = cone + amat(lm2,lm2,ir)
          enddo ! lm2
        enddo ! ir
!---> calculate non sph. wft. in i-th born approximation
        do j = 1,nsra
          do ir = irmind,irc1
              do lm2 = 1,lmmaxd
                pns(:,lm2,ir,j) = amat(:,lm2,ir)*pzlm(:,ir,j) + bmat(:,lm2,ir)*qzlm(:,ir,j)
              enddo ! lm2
          enddo ! ir
        enddo ! j
!-----------------------------------------------------------------------
! check convergence
      do j = 1,nsra
      do ir = irmind,irc1
      do lm1 = 1,lmmaxd
      do lm2 = 1,lmmaxd
        pns0(lm1,lm2,ir,j) = pns0(lm1,lm2,ir,j)-pns(lm1,lm2,ir,j)
      enddo ! lm2
      enddo ! lm1
      enddo ! ir
      enddo ! j
      
      err=0.d0
      
      do j=1,nsra
      call csout(pns0(1,1,irmind,j),pns1,lmmaxd**2,irmind,irmd,ipan, ircut)
      do lm1=1,lmmaxd
      do lm2=1,lmmaxd
        err=max(err,abs(pns1(lm1,lm2,irc1)))
      enddo ! lm2
      enddo ! lm1
      enddo ! j

      ! convergence check
      if(i == icst .and. err > 1d-3) then
        write(*,*)'regns.f: Fredholmholm equation does not converge'
        stop 'error 1 in regns.f'
      endif

      do j = 1,nsra
      do ir = irmind,irc1
      do lm1 = 1,lmmaxd
      do lm2 = 1,lmmaxd
        pns0(lm1,lm2,ir,j) = pns(lm1,lm2,ir,j)
      enddo ! lm2
      enddo ! lm1
      enddo ! ir
      enddo ! j
 
      enddo ! i
   
! end Fredholmholm equation
!-----------------------------------------------------------------------
    else ! Fredholm
!-----------------------------------------------------------------------
! begin Volterra equation
      do i = 0,icst
!---> set up integrands for i-th born approximation
        if (i == 0) then
          call wfint0(ader,bder,pzlm,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
        else
          call wfint(pns,ader,bder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
        end if
!---> call integration subroutines
        call csout(ader,amat,lmmaxd**2,irmind,irmd,ipan,ircut)
        call csout(bder,bmat,lmmaxd**2,irmind,irmd,ipan,ircut)
        do ir = irmind,irc1
          do lm2 = 1,lmmaxd
          do lm1 = 1,lmmaxd
            if(lm1 == lm2)then
              amat(lm1,lm2,ir) = cone - amat(lm1,lm2,ir)
            else
              amat(lm1,lm2,ir) = - amat(lm1,lm2,ir)
            endif
          enddo ! lm1
          enddo ! lm2
        enddo ! ir
!---> calculate non sph. wft. in i-th born approximation
        do j = 1,nsra
          do ir = irmind,irc1
            do lm2 = 1,lmmaxd
              do lm1 = 1,lmmaxd
                pns(lm1,lm2,ir,j) = amat(lm1,lm2,ir)*pzlm(lm1,ir,j) + bmat(lm1,lm2,ir)*qzlm(lm1,ir,j)
          enddo ! lm1
        enddo ! lm2
      enddo ! ir
      enddo ! j
!-----------------------------------------------------------------------
! check convergence
       do j = 1,nsra
       do ir = irmind,irc1
       do lm2 = 1,lmmaxd
       do lm1 = 1,lmmaxd
         pns0(lm1,lm2,ir,j) = pns0(lm1,lm2,ir,j) - pns(lm1,lm2,ir,j)
        enddo ! lm1
        enddo ! lm2
        enddo ! ir
        enddo ! j
  
       err = 0.d0
       do j=1,nsra
       call csout(pns0(1,1,irmind,j),pns1,lmmaxd**2,irmind,irmd,ipan, ircut)
       do lm2=1,lmmaxd
       do lm1=1,lmmaxd
        err=max(err,abs(pns1(lm1,lm2,irc1)))
      enddo ! lm1
      enddo ! lm2
      enddo ! j
  

      ! convergence check
      if(i == icst.and.err > 1d-3) then
      write(*,*)'regns.f: Volterra equation does not converge'
      stop 'error 1 in regns.f'
      endif

       do j = 1,nsra
       do ir = irmind,irc1
       do lm2 = 1,lmmaxd
       do lm1 = 1,lmmaxd
        pns0(lm1,lm2,ir,j) = pns(lm1,lm2,ir,j)
      enddo ! lm1
      enddo ! lm2
      enddo ! ir
      enddo ! j
  
!-----------------------------------------------------------------------
      enddo ! i
      call zgeinv1(amat(1,1,irc1),ar,br,ipiv,lmmaxd)
      do ir=irmind,irc1
      do lm2=1,lmmaxd
      do lm1=1,lmmaxd
        ader(lm1,lm2,ir)= zero
        bder(lm1,lm2,ir)= zero
      enddo ! lm1
      enddo ! lm2
      enddo ! ir

  
      do ir=irmind,irc1
      do lm2=1,lmmaxd
      do lm3=1,lmmaxd
      do lm1=1,lmmaxd
        ader(lm1,lm2,ir) = ader(lm1,lm2,ir) + amat(lm1,lm3,ir)*ar(lm3,lm2)
        bder(lm1,lm2,ir) = bder(lm1,lm2,ir) + bmat(lm1,lm3,ir)*ar(lm3,lm2)
      enddo ! lm1
      enddo ! lm3
      enddo ! lm2
      enddo ! ir

      do ir=irmind,irc1
      do lm2=1,lmmaxd
      do lm1=1,lmmaxd
        amat(lm1,lm2,ir) = ader(lm1,lm2,ir)
        bmat(lm1,lm2,ir) = bder(lm1,lm2,ir)
      enddo ! lm1
      enddo ! lm2
      enddo ! ir

      do j = 1,nsra
      do ir = irmind,irc1
      do lm2 = 1,lmmaxd
      do lm1 = 1,lmmaxd
        pns(lm1,lm2,ir,j) = amat(lm1,lm2,ir)*pzlm(lm1,ir,j) + bmat(lm1,lm2,ir)*qzlm(lm1,ir,j)
   enddo ! lm1
   enddo ! lm2
   enddo ! ir
   enddo ! j
! end Volterra equation
!-----------------------------------------------------------------------
  endif ! Fredholm
  
      do lm2 = 1,lmmaxd
        efac2 = efac(lm2)
!---> store alpha and t - matrix
        ar(:,lm2) = amat(:,lm2,irmind)
        br(:,lm2) = bmat(:,lm2,irc1)*efac(1:lmmaxd)*efac2/ek !---> t-matrix
      enddo ! lm2
!---> rescale with efac
      do j = 1,nsra
         do ir = irmind,irc1
            do lm2 = 1,lmmaxd
               efac2 = efac(lm2)
               pns(:,lm2,ir,j) = pns(:,lm2,ir,j)*efac2
         enddo ! lm2
      enddo ! ir
   enddo ! j
 
      endsubroutine ! regns
! ************************************************************************

      subroutine zgeinv1(a,u,aux,ipiv,n)
! ************************************************************************
!   - inverts a general double complex matrix a,
!   - the result is return in u,
!   - input matrix a is returned unchanged,
!   - aux is a auxiliary matrix,
!   - a,u and aux are of dimension (n,n),
! ------------------------------------------------------------------------
      implicit none
      integer, intent(in) :: n
      double complex, intent(in) :: a(n,*)
      double complex, intent(out) :: u(n,*)
      double complex, intent(inout) :: aux(n,*)
      integer, intent(inout) :: ipiv(*)
      
      double complex, parameter :: cone=(1.d0, 0.d0), zero=(0.d0, 0.d0)
      integer :: i, info
      external :: zcopy, zgetrs, zgetrf
! ------------------------------------------------------------------------
      u(:,1:n) = zero
      do i = 1, n
        u(i,i) = cone
      enddo ! i

      call zcopy(n*n,a,1,aux,1)
      call zgetrf(n,n,aux,n,ipiv,info)
      call zgetrs('n',n,n,aux,n,ipiv,u,n,info)
      
      endsubroutine
