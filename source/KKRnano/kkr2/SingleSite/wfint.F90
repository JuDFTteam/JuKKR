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
      
      endsubroutine 
