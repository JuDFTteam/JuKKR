    Subroutine calc_rho_ll_ss_lmdos(rll, ircut, ipan, icell, thetas, cleb, &
      icleb, iend, ifunm, lmsp, irws, drdi, dens, lmdos)
      Use mod_datatypes, Only: dp
      Implicit None

      Include 'inc.p'
      Integer :: lmmaxd
      Parameter (lmmaxd=(lmaxd+1)**2)
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
      Integer :: irmind
      Parameter (irmind=irmd-irnsd)
! non-sph. eigen states of single pot 
! derivative dr/di
      Integer :: iend, irws, lmdos


! local variables
      Complex (Kind=dp) :: rll(irmd, lmmaxd, lmmaxd), dens
      Real (Kind=dp) :: cleb(*), thetas(irid, nfund, *), drdi(irmd)
      Integer :: icleb(ncleb, 4), ifunm(natypd, lmpotd), lmsp(natypd, *), &
        ircut(0:ipand), ipan, icell, ifun

!     ..
!  ---> first calculate only the spherically symmetric contribution
!       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
      Real (Kind=dp) :: c0ll
      Complex (Kind=dp) :: clt
      Complex (Kind=dp), Allocatable :: rsp(:), rges(:)
      Integer :: lm1p, lm2p, lm3p, ir, j, i
      Integer :: ircutm(0:ipand)
!       multiplied with the shape functions...

!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)


!      WRITE(6,*) "In rho ll"


      Allocate (rges(irmd))
      Allocate (rsp(irmd))



      c0ll = 1.0E0_dp/sqrt(16.0E0_dp*atan(1.0E0_dp))
      rsp = 0E0_dp
      rges = 0E0_dp
!      STOP " "
      Do ir = 1, irmd
        rsp(ir) = rsp(ir) + rll(ir, lmdos, lmdos)
      End Do
!      WRITE(6,*) "IRCUT(1)",IRCUT(1)
      Do ir = 1, ircut(ipan)
        rges(ir) = rsp(ir)
      End Do
!      WRITE(6,*) "IRCUT(IPAN)",IRCUT(IPAN)
      If (ipan>1) Then
        Do ir = ircut(1) + 1, ircut(ipan)
          rges(ir) = rsp(ir)*c0ll*thetas(ir-ircut(1), 1, icell)
        End Do
      End If
!      WRITE(6,*) "IRCUT(IPAN)-IRMIND",IRCUT(IPAN)-IRMIND
!      WRITE(6,*) "IRMIND",IRMIND


!---> calculate the non spherically symmetric contribution

!          WRITE(156,*) "IFUN",IFUN
      Do j = 1, iend
        lm1p = icleb(j, 1)
        lm2p = icleb(j, 2)
        lm3p = icleb(j, 3)
        clt = cleb(j)
!          ELSE
!            DO IR = IRCUT(1)+1,IRCUT(IPAN)
!              RGES(IR) = RGES(IR)+
        If (ipan>1 .And. lmsp(icell,lm3p)>0) Then
          ifun = ifunm(icell, lm3p)
!     +             CLEB(J)*THETAS(IR-IRCUT(1),IFUN,ICELL)*
          If (lm1p==lm2p .And. lm1p==lmdos) Then
            Do ir = ircut(1) + 1, ircut(ipan)
              rges(ir) = rges(ir) + rll(ir, lm2p, lm1p)*cleb(j)*thetas(ir- &
                ircut(1), ifun, icell)
            End Do
!     +          (RLL(IR,LM2P,LM1P)+RLL(IR,LM1P,LM2P))
!            END DO




          End If
        End If

      End Do

      If (ipan==1) Then
        ircutm(0) = 0
        ircutm(1) = irws
      Else
        Do i = 0, ipan
          ircutm(i) = ircut(i)
        End Do
      End If
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
      Call csimpk(rges(:), dens, ipan, ircutm, drdi)
! SET ACCORDING TO lmax VALUE OF INPUTCARD
      Deallocate (rges)
      Deallocate (rsp)
!      PARAMETER ( NRD = 20000, KPOIBZ = 32000 )
    End Subroutine
