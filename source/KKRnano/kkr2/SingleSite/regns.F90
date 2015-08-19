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
!     Fredholm: true -> use Fredholm equation
!               false -> volterra equation
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
      double precision :: err
      integer :: i, ir, irc1, j,lm2,lm3, lm
      double complex pns0(lmmaxd,lmmaxd,irmind:irmd,2), pns1(lmmaxd,lmmaxd,irmind:irmd)
      integer ipiv(lmmaxd)
      double complex, parameter :: cone=(1.0d0,0.0d0), zero=(0.d0,0.d0)
      logical, parameter :: Volterra = .true. ! false ==> Fredholm equation
      
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
        call csout(pns0(1,1,irmind,j),pns1,lmmaxd**2,irmind,irmd,ipan, ircut)
        err = maxval(abs(pns1(:,:,irc1)))
        ! convergence check
        if (err > 1d-3) then
          if (Volterra) then
            write(*,*)'regns.f: Volterra equation does not converge'
          else
            write(*,*)'regns.f: Fredholm equation does not converge'
          endif
          stop 'error 1 in regns.f'
        endif
      enddo ! j
! end convergence check      
!-----------------------------------------------------------------------
      
    
      if (Volterra) then
!-----------------------------------------------------------------------
! only Volterra equation
        
        call zgeinv1(amat(1,1,irc1),ar,br,ipiv,lmmaxd) ! invert the last amat

        do ir = irmind, irc1
#define _USE_ZGEMM_in_REGNS_
#ifdef  _USE_ZGEMM_in_REGNS_
!!!! excerpt from zgemm.f, case TRANSA='N', TRANSB='N'
! ! ! ! ! ! ! *
! ! ! ! ! ! ! *           Form  C := alpha*A*B + beta*C.
! ! ! ! ! ! ! *
! ! ! ! ! ! !               DO 90 J = 1,N
! ! ! ! ! ! !                   IF (BETA.EQ.ZERO) THEN
! ! ! ! ! ! !                       DO 50 I = 1,M
! ! ! ! ! ! !                           C(I,J) = ZERO
! ! ! ! ! ! !    50                 CONTINUE
! ! ! ! ! ! !                   ELSE IF (BETA.NE.ONE) THEN
! ! ! ! ! ! !                       DO 60 I = 1,M
! ! ! ! ! ! !                           C(I,J) = BETA*C(I,J)
! ! ! ! ! ! !    60                 CONTINUE
! ! ! ! ! ! !                   END IF
! ! ! ! ! ! !                   DO 80 L = 1,K
! ! ! ! ! ! !                       IF (B(L,J).NE.ZERO) THEN
! ! ! ! ! ! !                           TEMP = ALPHA*B(L,J)
! ! ! ! ! ! !                           DO 70 I = 1,M
! ! ! ! ! ! !                               C(I,J) = C(I,J) + TEMP*A(I,L)
! ! ! ! ! ! !    70                     CONTINUE
! ! ! ! ! ! !                       END IF
! ! ! ! ! ! !    80             CONTINUE
! ! ! ! ! ! !    90         CONTINUE
!!!!! identify lm1=I, lm2=J, lm3=L, ALPHA=cone, BETA=zero, M,N,K=lmmaxd, A=amat, B=ar
          call ZGEMM('N','N',lmmaxd,lmmaxd,lmmaxd,cone,amat(1,1,ir),lmmaxd,ar,lmmaxd,zero,ader(1,1,ir),lmmaxd)
          call ZGEMM('N','N',lmmaxd,lmmaxd,lmmaxd,cone,bmat(1,1,ir),lmmaxd,ar,lmmaxd,zero,bder(1,1,ir),lmmaxd)
#else
          ! the following 8 lines should be written as two zgemm calls
          ader(:,:,ir) = zero
          bder(:,:,ir) = zero
          do lm2 = 1, lmmaxd
            do lm3 = 1, lmmaxd
              ader(:,lm2,ir) = ader(:,lm2,ir) + amat(:,lm3,ir)*ar(lm3,lm2)
              bder(:,lm2,ir) = bder(:,lm2,ir) + bmat(:,lm3,ir)*ar(lm3,lm2)
            enddo ! lm3
          enddo ! lm2
          ! end zgemm
#endif          
          ! overwrite amat and bmat, keep the valus in ader and bder (since these are intent(out)
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
