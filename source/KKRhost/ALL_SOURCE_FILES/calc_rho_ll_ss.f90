    Subroutine calc_rho_ll_ss(lmmax, rll, ircut, ipan, icell, thetas, cleb, &
      icleb, iend, ifunm, lmsp, irws, drdi, dens)
      Use mod_datatypes, Only: dp

      Implicit None

      Include 'inc.p'
      Integer :: lmmaxd
      Parameter (lmmaxd=(lmaxd+1)**2)
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
      Integer :: irmind
      Parameter (irmind=irmd-irnsd)
!.. Array Arguments ..
! non-sph. eigen states of single pot 
      Integer :: iend, lmmax, irws ! derivative dr/di

! local variables
!,DENS(:,:,:)
      Complex (Kind=dp) :: rll(irmd, lmmaxd, lmmaxd), dens
      Real (Kind=dp) :: cleb(*), thetas(irid, nfund, *), drdi(irmd) !                                  RGES_W(:,:,:,:), &
      Integer :: icleb(ncleb, 4), ifunm(natypd, lmpotd), lmsp(natypd, *), &
        ircut(0:ipand), ipan, icell, ifun
!                                  DENS_GESAMT(:,:), &
!                                  DENS_GESAMT_I1(:,:,:)
      Real (Kind=dp) :: c0ll
      Complex (Kind=dp) :: clt
      Complex (Kind=dp), Allocatable :: rsp(:), rges(:) !     ..
!  ---> first calculate only the spherically symmetric contribution
!       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
!       multiplied with the shape functions...
      Integer :: lm1p, lm2p, lm3p, ir, j, i
      Integer :: ircutm(0:ipand)

!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)


!      WRITE(6,*) "In rho ll"


      Allocate (rges(irmd))
      Allocate (rsp(irmd))


!WRITE(56,"((I5),(2e17.9))") IR,THETAS(IR-IRCUT(1),1,ICELL)
      c0ll = 1.0E0_dp/sqrt(16.0E0_dp*atan(1.0E0_dp))
      rsp = 0E0_dp
      rges = 0E0_dp

      Do lm1p = 1, lmmax
        Do ir = 1, irmd
          rsp(ir) = rsp(ir) + rll(ir, lm1p, lm1p)
        End Do
      End Do
!      STOP " "
      Do ir = 1, ircut(ipan)
        rges(ir) = rsp(ir)
      End Do
!WRITE(6,*) "IRCUT(1)",IRCUT(1)
      If (ipan>1) Then
        Do ir = ircut(1) + 1, ircut(ipan)
!WRITE(6,*) "IRCUT(IPAN)",IRCUT(IPAN)
          rges(ir) = rsp(ir)*c0ll*thetas(ir-ircut(1), 1, icell)
        End Do
      End If
!WRITE(6,*) "IRCUT(IPAN)-IRMIND",IRCUT(IPAN)-IRMIND
!WRITE(6,*) "IRMIND",IRMIND


!---> calculate the non spherically symmetric contribution

!WRITE(156,*) "IFUN",IFUN
      Do j = 1, iend
        lm1p = icleb(j, 1)
        lm2p = icleb(j, 2)
        lm3p = icleb(j, 3)
        clt = cleb(j)



        If (ipan>1 .And. lmsp(icell,lm3p)>0) Then
          ifun = ifunm(icell, lm3p)

          If (lm1p==lm2p) Then
            Do ir = ircut(1) + 1, ircut(ipan)
              rges(ir) = rges(ir) + rll(ir, lm2p, lm1p)*cleb(j)*thetas(ir- &
                ircut(1), ifun, icell)
            End Do
          Else
            Do ir = ircut(1) + 1, ircut(ipan)
              rges(ir) = rges(ir) + cleb(j)*thetas(ir-ircut(1), ifun, icell)*( &
                rll(ir,lm2p,lm1p)+rll(ir,lm1p,lm2p))
            End Do
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

      Call csimpk(rges(:), dens, ipan, ircutm, drdi)
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
      Deallocate (rges)
      Deallocate (rsp)
! SET ACCORDING TO lmax VALUE OF INPUTCARD
    End Subroutine
