!>    @param[out] efac
!>    @param[in]  pz
!>    @param[in]  qz
!>    @param[in]  fz
!>    @param[in]  sz
!>    @param[in]  nsra
!>    @param[out] pzlm
!>    @param[out] qzlm
!>    @param[out] pzekdr
!>    @param[out] qzekdr
!>    @param[in]  ek
      subroutine wftsca(drdi, efac, pz, qz, fz, sz, nsra, pzlm, qzlm,  &
        pzekdr, qzekdr, ek, loflm, irmind, irmd, lmaxd, lmmaxd)
      implicit none
!-----------------------------------------------------------------------
!                 r. zeller      oct. 1993
!-----------------------------------------------------------------------
      double complex, intent(in) :: ek
      integer, intent(in) :: irmd, irmind, lmaxd, lmmaxd, nsra
      double complex, intent(out) :: efac(lmmaxd)
      double complex, intent(in) :: fz(irmd,0:lmaxd), pz(irmd,0:lmaxd)
      double complex, intent(out) :: pzekdr(lmmaxd,irmind:irmd,2)
      double complex, intent(out) :: pzlm(lmmaxd,irmind:irmd,2)
      double complex, intent(in) :: qz(irmd,0:lmaxd)
      double complex, intent(out) :: qzekdr(lmmaxd,irmind:irmd,2)
      double complex, intent(out) :: qzlm(lmmaxd,irmind:irmd,2)
      double complex, intent(in) :: sz(irmd,0:lmaxd)
      double precision, intent(in) :: drdi(*)
      integer, intent(in) :: loflm(*)
      
      double complex, parameter :: cone=(1.d0,0.d0)
      double complex :: efac1, v1
      integer :: ir, j, l, l1, lm1, m

!
!---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!
!
      efac(1) = cone
      v1 = cone
      do l = 1, lmaxd
        v1 = v1*ek/dble(2*l-1)
        do m = -l, l
          lm1 = l*(l+1)+m+1
          efac(lm1) = v1
        enddo ! m
      enddo ! l
!
!---> get wfts of same magnitude by scaling with efac
!
      do lm1 = 1, lmmaxd
        l1 = loflm(lm1)
        efac1 = efac(lm1)
        do ir = irmind, irmd
          pzlm(lm1,ir,1) = pz(ir,l1)/efac1
          qzlm(lm1,ir,1) = qz(ir,l1)*efac1
        enddo ! ir
        if (nsra == 2) then
          do ir = irmind, irmd
            pzlm(lm1,ir,2) = fz(ir,l1)/efac1
            qzlm(lm1,ir,2) = sz(ir,l1)*efac1
          enddo ! ir
        endif ! nsra == 2

        do j = 1, nsra
          do ir = irmind, irmd
            pzekdr(lm1,ir,j) = pzlm(lm1,ir,j)*ek*drdi(ir)
            qzekdr(lm1,ir,j) = qzlm(lm1,ir,j)*ek*drdi(ir)
          enddo ! ir
        enddo ! j
   
      enddo ! lm1

      endsubroutine ! wftsca
