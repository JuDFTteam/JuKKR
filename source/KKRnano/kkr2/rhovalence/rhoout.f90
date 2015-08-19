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
  !       summands are multiplied by a factor of 2 in the case of lm1 .ne.
  !       lm2 or lm01 .ne. lm02 .
  !
  !             (see notes by b.drittler)
  !
  !                               b.drittler   aug. 1988
  !-----------------------------------------------------------------------
  implicit none

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

  external :: zgemm
  double complex, parameter :: cone=(1.d0,0.d0), zero=(0.d0,0.d0)
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
