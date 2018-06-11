!     This subroutine computes a matrix that is the basis for constructing
!     the KKR representation of the torque operator. It is adapted from the
!     CALC_RHO_LL_SS subroutine, but the spin dependent part, i.e., the exhange
!     field, replaces the shape function in the integration.
!
!                                     Guillaume Geranton, September 2014
    Subroutine calc_torq_ll_ss(lmmax, rll, ircut, ipan, icell, cleb, icleb, &
      iend, ifunm, lmsp, irws, drdi, dens, visp, nspin, iatom, vins, irmin)
      Use mod_datatypes, Only: dp

      Implicit None

      Include 'inc.p'
      Integer :: lmmaxd
      Parameter (lmmaxd=(lmaxd+1)**2)
      Integer :: lmpotd
      Parameter (lmpotd=(lpotd+1)**2)
      Integer :: irmind
      Parameter (irmind=irmd-irnsd)
!     .. Array Arguments ..
! non-sph. eigen states of single pot  &
      Integer :: iend, lmmax, irws, nspin, iatom, irmin ! derivative dr/di  &
!              spherical part of the potential  &
! non-sph. part of the potential
      Complex (Kind=dp) :: rll(irmd, lmmaxd, lmmaxd), dens
! local variables
      Real (Kind=dp) :: cleb(*), drdi(irmd), visp(irmd, *), &
        vins(irmind:irmd, lmpotd, *)
      Integer :: icleb(ncleb, 4), ifunm(natypd, lmpotd), lmsp(natypd, *), &
        ircut(0:ipand), ipan, icell, ifun

!     ..
!  ---> first calculate only the spherically symmetric contribution
      Real (Kind=dp) :: c0ll
      Complex (Kind=dp) :: clt
      Complex (Kind=dp), Allocatable :: rsp(:), rges(:)
      Integer :: lm1p, lm2p, lm3p, ir, j, i
      Integer :: ircutm(0:ipand)
!       (for all points r; if r>r_MT (or IR> IRMIN),the density has to
      External :: test
      Logical :: test
!       multiplied with the shape functions...

!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)



!     Compute spherical contribution to the torque (LM=1)
!     Sph. potential has to be multiplied by sqrt(4 PI) !
      Allocate (rges(irmd))
      Allocate (rsp(irmd))

      c0ll = 1.0E0_dp/sqrt(16.0E0_dp*atan(1.0E0_dp))
      rsp = 0E0_dp
      rges = 0E0_dp

! cut contributions from outside the MT if recquired

      Do lm1p = 1, lmmax
        Do ir = 1, irmd
          rsp(ir) = rsp(ir) + rll(ir, lm1p, lm1p)*c0ll*(-1)*(visp(ir,nspin*( &
            iatom-1)+2)-visp(ir,nspin*(iatom-1)+1))*0.5_dp* &
            sqrt(16.0E0_dp*atan(1.0E0_dp))
! always >= 2 here
        End Do
      End Do

      Do ir = 1, irmd
!--->   calculate the non spherically symmetric contribution
        If (test('ONLYMT  ') .And. (ir>ircut(1))) Then
          rges(ir) = 0
        Else
          rges(ir) = rsp(ir)
        End If
      End Do
!             DO 150 IR = IRCUT(1)+1,IRCUT(IPAN)
      If (.Not. test('ONLYSPH ')) Then
        Do j = 1, iend
          lm1p = icleb(j, 1)
          lm2p = icleb(j, 2)
          lm3p = icleb(j, 3) !             DO IR = IRCUT(1)+1,IRCUT(IPAN)
          clt = cleb(j)


          If (ipan>1 .And. lmsp(icell,lm3p)>0) Then
            ifun = ifunm(icell, lm3p)
            If (lm1p==lm2p) Then

              Do ir = irmin, ircut(ipan)
                rges(ir) = rges(ir) + rll(ir, lm2p, lm1p)*cleb(j)*(-1)*(vins( &
                  ir,lm3p,nspin*(iatom-1)+2)-vins(ir,lm3p,nspin*(iatom-1)+1))* &
                  0.5_dp
              End Do
            Else

              Do ir = irmin, ircut(ipan)
                rges(ir) = rges(ir) + cleb(j)*(-1)*(vins(ir,lm3p,nspin*(iatom- &
                  1)+2)-vins(ir,lm3p,nspin*(iatom-1)+1))*0.5_dp* &
                  (rll(ir,lm2p,lm1p)+rll(ir,lm1p,lm2p))
              End Do
            End If
          End If

        End Do
      End If
!     This subroutine computes a matrix that is the basis for constructing
      If (ipan==1) Then
        ircutm(0) = 0
        ircutm(1) = irws
      Else
        Do i = 0, ipan
          ircutm(i) = ircut(i)
        End Do
      End If
!     the KKR representation of the torque operator. It is adapted from the
      Call csimpk(rges(:), dens, ipan, ircutm, drdi)
!     CALC_RHO_LL_SS subroutine, but the spin dependent part, i.e., the exhange
      Deallocate (rges)
      Deallocate (rsp)
!     field, replaces the shape function in the integration.
    End Subroutine
