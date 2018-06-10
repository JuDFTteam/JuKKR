subroutine calc_rho_ll_ss_lmdos(rll, ircut, ipan, icell, thetas, cleb, icleb, &
  iend, ifunm, lmsp, irws, drdi, dens, lmdos)
  implicit none

  include 'inc.p'
  integer :: lmmaxd
  parameter (lmmaxd=(lmaxd+1)**2)
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  integer :: irmind
  parameter (irmind=irmd-irnsd)
! non-sph. eigen states of single pot 
! derivative dr/di
  integer :: iend, irws, lmdos


  double complex :: rll(irmd, lmmaxd, lmmaxd), & ! local variables
    dens
  double precision :: cleb(*), thetas(irid, nfund, *), drdi(irmd) 
  integer :: icleb(ncleb, 4), ifunm(natypd, lmpotd), lmsp(natypd, *), &
    ircut(0:ipand), ipan, icell, ifun

!     ..
!  ---> first calculate only the spherically symmetric contribution
!       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
  double precision :: c0ll
  double complex :: clt
  double complex, allocatable :: rsp(:), rges(:)
  integer :: lm1p, lm2p, lm3p, ir, j, i
  integer :: ircutm(0:ipand)
!       multiplied with the shape functions...

!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)


!      WRITE(6,*) "In rho ll"


  allocate (rges(irmd))
  allocate (rsp(irmd))



  c0ll = 1.0d0/dsqrt(16.0d0*datan(1.0d0))
  rsp = 0d0
  rges = 0d0
!      STOP " "
  do ir = 1, irmd
    rsp(ir) = rsp(ir) + rll(ir, lmdos, lmdos)
  end do
!      WRITE(6,*) "IRCUT(1)",IRCUT(1)
  do ir = 1, ircut(ipan)
    rges(ir) = rsp(ir)
  end do
!      WRITE(6,*) "IRCUT(IPAN)",IRCUT(IPAN)
  if (ipan>1) then
    do ir = ircut(1) + 1, ircut(ipan)
      rges(ir) = rsp(ir)*c0ll*thetas(ir-ircut(1), 1, icell)
    end do
  end if
!      WRITE(6,*) "IRCUT(IPAN)-IRMIND",IRCUT(IPAN)-IRMIND
!      WRITE(6,*) "IRMIND",IRMIND


!---> calculate the non spherically symmetric contribution

!          WRITE(156,*) "IFUN",IFUN
  do j = 1, iend
    lm1p = icleb(j, 1)
    lm2p = icleb(j, 2)
    lm3p = icleb(j, 3)
    clt = cleb(j)
!          ELSE
!            DO IR = IRCUT(1)+1,IRCUT(IPAN)
!              RGES(IR) = RGES(IR)+
    if (ipan>1 .and. lmsp(icell,lm3p)>0) then
      ifun = ifunm(icell, lm3p)
!     +             CLEB(J)*THETAS(IR-IRCUT(1),IFUN,ICELL)*
      if (lm1p==lm2p .and. lm1p==lmdos) then
        do ir = ircut(1) + 1, ircut(ipan)
          rges(ir) = rges(ir) + rll(ir, lm2p, lm1p)*cleb(j)*thetas(ir-ircut(1) &
            , ifun, icell)
        end do
!     +          (RLL(IR,LM2P,LM1P)+RLL(IR,LM1P,LM2P))
!            END DO




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
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
  call csimpk(rges(:), dens, ipan, ircutm, drdi)
! SET ACCORDING TO lmax VALUE OF INPUTCARD
  deallocate (rges)
  deallocate (rsp)
!      PARAMETER ( NRD = 20000, KPOIBZ = 32000 )
end subroutine
