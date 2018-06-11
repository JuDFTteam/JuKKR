    Subroutine spinorbit_ham(lmax, lmmaxd, vins, rnew, e, z, c, socscale, &
      nspin, lmpotd, theta, phi, ipan_intervall, rpan_intervall, npan_tot, &
      ncheb, irmdnew, nrmaxd, vnspll, vnspll1, mode)
      Use mod_datatypes, Only: dp
      Implicit None

      Integer :: lmax, lmmaxd, nspin, npan_tot, ncheb, irmdnew, nrmaxd
      Integer :: lmpotd
      Real (Kind=dp) :: c, z
      Complex (Kind=dp) :: e
      Real (Kind=dp) :: socscale
      Real (Kind=dp) :: vins(irmdnew, lmpotd, nspin), rnew(nrmaxd), &
        rpan_intervall(0:npan_tot)
      Complex (Kind=dp) :: vnspll(2*lmmaxd, 2*lmmaxd, irmdnew)
      Complex (Kind=dp) :: vnspll1(2*lmmaxd, 2*lmmaxd, irmdnew)
      Integer :: ipan_intervall(0:npan_tot)
      Real (Kind=dp) :: vr(irmdnew), dvdr(irmdnew)
      Real (Kind=dp) :: rmass(irmdnew), hsofac(irmdnew)
      Real (Kind=dp) :: rnucl, atn, widthfac, phi, theta
      Integer :: ir, ip, lm1, lm2, ispin, irmin, irmax, ncoll
      Complex (Kind=dp) :: lsmh(2*lmmaxd, 2*lmmaxd), temp
      Real (Kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)
      Character (Len=*) :: mode
      Logical :: test, opt
      External :: test, opt

      vnspll1 = (0E0_dp, 0E0_dp)
      vr = 0E0_dp
      Do ispin = 1, nspin
        Do ir = 1, ipan_intervall(npan_tot)
          vr(ir) = vr(ir) + vins(ir, 1, ispin)/nspin
        End Do
      End Do
! derivative of potential
      dvdr = 0E0_dp
      Call getclambdacinv(ncheb, clambdacinv)
      Do ip = 1, npan_tot
        irmin = ipan_intervall(ip-1) + 1
        irmax = ipan_intervall(ip)
        widthfac = 2E0_dp/(rpan_intervall(ip)-rpan_intervall(ip-1))
        Call dgemv('N', ncheb+1, ncheb+1, 1E0_dp, clambdacinv, ncheb+1, &
          vr(irmin:irmax), 1, 0E0_dp, dvdr(irmin:irmax), 1)
        dvdr(irmin:irmax) = dvdr(irmin:irmax)*widthfac
      End Do
! core potential
      If (z>24E0_dp) Then
        atn = -16.1532921_dp + 2.70335346_dp*z
      Else
        atn = 0.03467714_dp + 2.04820786_dp*z
      End If
      rnucl = 1.2E0_dp/0.529177E0_dp*atn**(1._dp/3E0_dp)*1.E-5_dp

      Do ir = 1, ipan_intervall(npan_tot)
        If (rnew(ir)<=rnucl) Then
!        DVDR(IR)=DVDR(IR)+2d0*Z*RNEW(IR)/RNUCL**3d0
        Else
!        DVDR(IR)=DVDR(IR)+2d0*Z/RNEW(IR)**2d0
        End If
        dvdr(ir) = dvdr(ir) + 2E0_dp*z/rnew(ir)**2E0_dp
      End Do
! contruct LS matrix

      Call spin_orbit_compl(lmax, lmmaxd, lsmh)

! roate LS matrix
      ncoll = 1
      If (ncoll==1) Then
        Call rotatematrix(lsmh, theta, phi, lmmaxd, 1)
      End If

      If (mode=='transpose') Then
        Do lm1 = 1, 2*lmmaxd
          Do lm2 = 1, lm1 - 1
            temp = lsmh(lm2, lm1)
            lsmh(lm2, lm1) = lsmh(lm1, lm2)
            lsmh(lm1, lm2) = temp
          End Do
        End Do
      Else If (mode=='1') Then
      End If
! contruct prefactor of spin-orbit hamiltonian

      hsofac = 0E0_dp
      If (test('NOSOC   ') .Or. z<1E-6_dp) Then
        Do ir = 1, irmdnew
          Do lm1 = 1, 2*lmmaxd
            Do lm2 = 1, 2*lmmaxd
              vnspll1(lm1, lm2, ir) = vnspll(lm1, lm2, ir)
            End Do
          End Do
        End Do
      Else
        Do ir = 1, irmdnew
          rmass(ir) = 0.5E0_dp - 0.5E0_dp/c**2*((vr(ir)-real(e))-2E0_dp*z/rnew &
            (ir))
          hsofac(ir) = socscale/(2E0_dp*rmass(ir)**2*c**2*rnew(ir))*dvdr(ir)

! add to potential
          Do lm1 = 1, 2*lmmaxd
            Do lm2 = 1, 2*lmmaxd
              vnspll1(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + &
                hsofac(ir)*lsmh(lm1, lm2)
            End Do
          End Do
        End Do
      End If
    End Subroutine
