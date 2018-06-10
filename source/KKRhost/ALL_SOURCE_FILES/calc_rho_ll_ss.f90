subroutine calc_rho_ll_ss(lmmax, rll, ircut, ipan, icell, thetas, cleb, icleb, &
  iend, ifunm, lmsp, irws, drdi, dens)

  implicit none

  include 'inc.p'
  integer :: lmmaxd
  parameter (lmmaxd=(lmaxd+1)**2)
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  integer :: irmind
  parameter (irmind=irmd-irnsd)
!.. Array Arguments ..
! non-sph. eigen states of single pot 
  integer :: iend, lmmax, irws ! derivative dr/di

! local variables
  double complex :: rll(irmd, lmmaxd, lmmaxd), & !,DENS(:,:,:)
    dens
  double precision :: cleb(*), thetas(irid, nfund, *), drdi(irmd) !                                  RGES_W(:,:,:,:), &
  integer :: icleb(ncleb, 4), ifunm(natypd, lmpotd), lmsp(natypd, *), &
    ircut(0:ipand), ipan, icell, ifun
!                                  DENS_GESAMT(:,:), &
!                                  DENS_GESAMT_I1(:,:,:)
  double precision :: c0ll
  double complex :: clt
  double complex, allocatable :: rsp(:), rges(:) !     ..
!  ---> first calculate only the spherically symmetric contribution
!       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
!       multiplied with the shape functions...
  integer :: lm1p, lm2p, lm3p, ir, j, i
  integer :: ircutm(0:ipand)

!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)


!      WRITE(6,*) "In rho ll"


  allocate (rges(irmd))
  allocate (rsp(irmd))


!WRITE(56,"((I5),(2e17.9))") IR,THETAS(IR-IRCUT(1),1,ICELL)
  c0ll = 1.0d0/dsqrt(16.0d0*datan(1.0d0))
  rsp = 0d0
  rges = 0d0

  do lm1p = 1, lmmax
    do ir = 1, irmd
      rsp(ir) = rsp(ir) + rll(ir, lm1p, lm1p)
    end do
  end do
!      STOP " "
  do ir = 1, ircut(ipan)
    rges(ir) = rsp(ir)
  end do
!WRITE(6,*) "IRCUT(1)",IRCUT(1)
  if (ipan>1) then
    do ir = ircut(1) + 1, ircut(ipan)
!WRITE(6,*) "IRCUT(IPAN)",IRCUT(IPAN)
      rges(ir) = rsp(ir)*c0ll*thetas(ir-ircut(1), 1, icell)
    end do
  end if
!WRITE(6,*) "IRCUT(IPAN)-IRMIND",IRCUT(IPAN)-IRMIND
!WRITE(6,*) "IRMIND",IRMIND


!---> calculate the non spherically symmetric contribution

!WRITE(156,*) "IFUN",IFUN
  do j = 1, iend
    lm1p = icleb(j, 1)
    lm2p = icleb(j, 2)
    lm3p = icleb(j, 3)
    clt = cleb(j)



    if (ipan>1 .and. lmsp(icell,lm3p)>0) then
      ifun = ifunm(icell, lm3p)

      if (lm1p==lm2p) then
        do ir = ircut(1) + 1, ircut(ipan)
          rges(ir) = rges(ir) + rll(ir, lm2p, lm1p)*cleb(j)*thetas(ir-ircut(1) &
            , ifun, icell)
        end do
      else
        do ir = ircut(1) + 1, ircut(ipan)
          rges(ir) = rges(ir) + cleb(j)*thetas(ir-ircut(1), ifun, icell)*(rll( &
            ir,lm2p,lm1p)+rll(ir,lm1p,lm2p))
        end do
      end if
    end if

  end do

  if (ipan==1) then
    ircutm(0) = 0
    ircutm(1) = irws
  else
    do i = 0, ipan
      ircutm(i) = ircut(i)
    end do
  end if

  call csimpk(rges(:), dens, ipan, ircutm, drdi)
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
  deallocate (rges)
  deallocate (rsp)
! SET ACCORDING TO lmax VALUE OF INPUTCARD
end subroutine
