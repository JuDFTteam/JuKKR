      subroutine wfint(qns,cder,dder,qzekdr,pzekdr,vnspll,nsra,irmind, irmd,lmmaxd)
      implicit none
!-----------------------------------------------------------------------
!      determines the integrands cder, dder or ader, bder in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinspll.
!
!      r. zeller      aug. 1994
!-----------------------------------------------------------------------
!     .. scalar arguments ..
      integer, intent(in) :: irmd,irmind,lmmaxd,nsra
!     ..
!     .. array arguments ..
      double complex, intent(out) :: cder(lmmaxd,lmmaxd,irmind:irmd)
      double complex, intent(out) :: dder(lmmaxd,lmmaxd,irmind:irmd)
      double complex, intent(in) :: pzekdr(lmmaxd,irmind:irmd,2)
      double complex, intent(in) :: qns(lmmaxd,lmmaxd,irmind:irmd,2)
      double complex, intent(in) :: qzekdr(lmmaxd,irmind:irmd,2)
      double precision, intent(in) :: vnspll(lmmaxd,lmmaxd,irmind:irmd)

#define _USE_ZGEMM_in_WFINT_      
#ifndef _USE_ZGEMM_in_WFINT_     

      external :: dgemm
      integer :: ir, lm
      double precision :: qnsi(lmmaxd,lmmaxd), qnsr(lmmaxd,lmmaxd)
      double precision :: vtqnsi(lmmaxd,lmmaxd), vtqnsr(lmmaxd,lmmaxd)
      
      do ir = irmind, irmd
      
        qnsr =  dble(qns(:,:,ir,1))
        qnsi = dimag(qns(:,:,ir,1))
        
        call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsr,lmmaxd,0.d0,vtqnsr,lmmaxd)
        call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsi,lmmaxd,0.d0,vtqnsi,lmmaxd)
        
        do lm = 1, lmmaxd
          cder(:,lm,ir) = qzekdr(:,ir,1)*dcmplx(vtqnsr(:,lm), vtqnsi(:,lm))
          dder(:,lm,ir) = pzekdr(:,ir,1)*dcmplx(vtqnsr(:,lm), vtqnsi(:,lm))
        enddo ! lm
   
        if (nsra == 2) then
          qnsr =  dble(qns(:,:,ir,2))
          qnsi = dimag(qns(:,:,ir,2))
   
          call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsr,lmmaxd,0.d0,vtqnsr,lmmaxd)
          call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsi,lmmaxd,0.d0,vtqnsi,lmmaxd)
          do lm = 1, lmmaxd
            cder(:,lm,ir) = cder(:,lm,ir) + qzekdr(:,ir,2)*dcmplx(vtqnsr(:,lm), vtqnsi(:,lm))
            dder(:,lm,ir) = dder(:,lm,ir) + pzekdr(:,ir,2)*dcmplx(vtqnsr(:,lm), vtqnsi(:,lm))
          enddo ! lm
        endif ! nsra

      enddo ! ir
#else
! new version
      external :: zgemm
      integer :: ir, lm
      double complex :: cvnspll(lmmaxd,lmmaxd), vtqns(lmmaxd,lmmaxd)
      double complex, parameter :: zero=(0.d0,0.d0), cone=(1.d0,0.d0)
      
      do ir = irmind, irmd
      
        cvnspll = dcmplx(vnspll(:,:,ir), 0.d0) ! convert the potential to complex

        call zgemm('n','n',lmmaxd,lmmaxd,lmmaxd,cone,cvnspll,lmmaxd,qns(1,1,ir,1),lmmaxd,zero,vtqns,lmmaxd)
        
        do lm = 1, lmmaxd
          cder(:,lm,ir) = qzekdr(:,ir,1)*vtqns(:,lm)
          dder(:,lm,ir) = pzekdr(:,ir,1)*vtqns(:,lm)
        enddo ! lm
   
        if (nsra == 2) then
          call zgemm('n','n',lmmaxd,lmmaxd,lmmaxd,cone,cvnspll,lmmaxd,qns(1,1,ir,2),lmmaxd,zero,vtqns,lmmaxd)
   
          do lm = 1, lmmaxd
            cder(:,lm,ir) = cder(:,lm,ir) + qzekdr(:,ir,2)*vtqns(:,lm)
            dder(:,lm,ir) = dder(:,lm,ir) + pzekdr(:,ir,2)*vtqns(:,lm)
          enddo ! lm
        endif ! nsra

      enddo ! ir
#endif
      endsubroutine 

      
      subroutine wfint0(cder,dder,qzlm,qzekdr,pzekdr,vnspll,nsra, irmind,irmd,lmmaxd)
      implicit none
!-----------------------------------------------------------------------
!      determines the integrands cder, dder or ader, bder in the
!        integral equations for the non-spherical wavefunctions from
!        the non-spherical contributions of the potential vinspll.
!        (this subroutine is used in zeroth order born approximation,
!         otherwise subroutine wfint must be used)
!      r. zeller      aug. 1994
!-----------------------------------------------------------------------
      integer, intent(in) :: irmd,irmind,lmmaxd,nsra
      double complex, intent(out) :: cder(lmmaxd,lmmaxd,irmind:irmd)
      double complex, intent(out) :: dder(lmmaxd,lmmaxd,irmind:irmd)
      double complex, intent(in) :: pzekdr(lmmaxd,irmind:irmd,2)
      double complex, intent(in) :: qzekdr(lmmaxd,irmind:irmd,2)
      double complex, intent(in) :: qzlm(lmmaxd,irmind:irmd,2)
      double precision, intent(in) :: vnspll(lmmaxd,lmmaxd,irmind:irmd)
      
      integer :: ir, lm

      if (nsra == 2) then
        do ir = irmind, irmd
          do lm = 1,lmmaxd
            cder(:,lm,ir) = qzekdr(:,ir,1)*vnspll(:,lm,ir)*qzlm(lm,ir,1) + qzekdr(:,ir,2)*vnspll(:,lm,ir)*qzlm(lm,ir,2)
            dder(:,lm,ir) = pzekdr(:,ir,1)*vnspll(:,lm,ir)*qzlm(lm,ir,1) + pzekdr(:,ir,2)*vnspll(:,lm,ir)*qzlm(lm,ir,2)
          enddo ! lm
        enddo ! ir
      else
        do ir = irmind, irmd
          do lm = 1,lmmaxd
            cder(:,lm,ir) = qzekdr(:,ir,1)*vnspll(:,lm,ir)*qzlm(lm,ir,1)
            dder(:,lm,ir) = pzekdr(:,ir,1)*vnspll(:,lm,ir)*qzlm(lm,ir,1)
          enddo ! lm
        enddo ! ir
      endif
      
      endsubroutine
      