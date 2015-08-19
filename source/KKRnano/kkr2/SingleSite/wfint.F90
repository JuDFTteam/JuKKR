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
      integer :: ir, lm1, lm2
      double precision :: qnsi(lmmaxd,lmmaxd), qnsr(lmmaxd,lmmaxd)
      double precision :: vtqnsi(lmmaxd,lmmaxd), vtqnsr(lmmaxd,lmmaxd)
      
      do 90 ir = irmind,irmd
        do 20 lm2 = 1,lmmaxd
          do 10 lm1 = 1,lmmaxd
            qnsr(lm1,lm2) =  dble(qns(lm1,lm2,ir,1))
            qnsi(lm1,lm2) = dimag(qns(lm1,lm2,ir,1))
   10     continue
   20   continue
        call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsr,lmmaxd,0.d0,vtqnsr,lmmaxd)
        call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsi,lmmaxd,0.d0,vtqnsi,lmmaxd)
        do 40 lm1 = 1,lmmaxd
          do 30 lm2 = 1,lmmaxd
            cder(lm1,lm2,ir) = qzekdr(lm1,ir,1)*dcmplx(vtqnsr(lm1,lm2),vtqnsi(lm1,lm2))
            dder(lm1,lm2,ir) = pzekdr(lm1,ir,1)*dcmplx(vtqnsr(lm1,lm2),vtqnsi(lm1,lm2))
   30     continue
   40   continue
        if (nsra == 2) then
          do 60 lm2 = 1,lmmaxd
            do 50 lm1 = 1,lmmaxd
              qnsr(lm1,lm2) =  dble(qns(lm1,lm2,ir,2))
              qnsi(lm1,lm2) = dimag(qns(lm1,lm2,ir,2))
   50       continue
   60     continue
          call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsr,lmmaxd,0.d0,vtqnsr,lmmaxd)
          call dgemm('n','n',lmmaxd,lmmaxd,lmmaxd,1.d0,vnspll(1,1,ir),lmmaxd,qnsi,lmmaxd,0.d0,vtqnsi,lmmaxd)
          do 80 lm2 = 1,lmmaxd
            do 70 lm1 = 1,lmmaxd
              cder(lm1,lm2,ir) = cder(lm1,lm2,ir) + qzekdr(lm1,ir,2)*dcmplx(vtqnsr(lm1,lm2),vtqnsi(lm1,lm2))
              dder(lm1,lm2,ir) = dder(lm1,lm2,ir) + pzekdr(lm1,ir,2)*dcmplx(vtqnsr(lm1,lm2),vtqnsi(lm1,lm2))
   70       continue
   80     continue
        end if

   90 continue
      endsubroutine 
