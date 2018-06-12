subroutine spinorbit_ham(lmax, lmmaxd, vins, rnew, e, z, c, socscale, nspin, &
  lmpotd, theta, phi, ipan_intervall, rpan_intervall, npan_tot, ncheb, &
  irmdnew, nrmaxd, vnspll, vnspll1, mode)
  use :: mod_datatypes, only: dp
  implicit none

  integer :: lmax, lmmaxd, nspin, npan_tot, ncheb, irmdnew, nrmaxd
  integer :: lmpotd
  real (kind=dp) :: c, z
  complex (kind=dp) :: e
  real (kind=dp) :: socscale
  real (kind=dp) :: vins(irmdnew, lmpotd, nspin), rnew(nrmaxd), &
    rpan_intervall(0:npan_tot)
  complex (kind=dp) :: vnspll(2*lmmaxd, 2*lmmaxd, irmdnew)
  complex (kind=dp) :: vnspll1(2*lmmaxd, 2*lmmaxd, irmdnew)
  integer :: ipan_intervall(0:npan_tot)
  real (kind=dp) :: vr(irmdnew), dvdr(irmdnew)
  real (kind=dp) :: rmass(irmdnew), hsofac(irmdnew)
  real (kind=dp) :: rnucl, atn, widthfac, phi, theta
  integer :: ir, ip, lm1, lm2, ispin, irmin, irmax, ncoll
  complex (kind=dp) :: lsmh(2*lmmaxd, 2*lmmaxd), temp
  real (kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)
  character (len=*) :: mode
  logical :: test, opt
  external :: test, opt

  vnspll1 = (0e0_dp, 0e0_dp)
  vr = 0e0_dp
  do ispin = 1, nspin
    do ir = 1, ipan_intervall(npan_tot)
      vr(ir) = vr(ir) + vins(ir, 1, ispin)/nspin
    end do
  end do
  ! derivative of potential
  dvdr = 0e0_dp
  call getclambdacinv(ncheb, clambdacinv)
  do ip = 1, npan_tot
    irmin = ipan_intervall(ip-1) + 1
    irmax = ipan_intervall(ip)
    widthfac = 2e0_dp/(rpan_intervall(ip)-rpan_intervall(ip-1))
    call dgemv('N', ncheb+1, ncheb+1, 1e0_dp, clambdacinv, ncheb+1, &
      vr(irmin:irmax), 1, 0e0_dp, dvdr(irmin:irmax), 1)
    dvdr(irmin:irmax) = dvdr(irmin:irmax)*widthfac
  end do
  ! core potential
  if (z>24e0_dp) then
    atn = -16.1532921_dp + 2.70335346_dp*z
  else
    atn = 0.03467714_dp + 2.04820786_dp*z
  end if
  rnucl = 1.2e0_dp/0.529177e0_dp*atn**(1._dp/3e0_dp)*1.e-5_dp

  do ir = 1, ipan_intervall(npan_tot)
    if (rnew(ir)<=rnucl) then
      ! DVDR(IR)=DVDR(IR)+2d0*Z*RNEW(IR)/RNUCL**3d0
    else
      ! DVDR(IR)=DVDR(IR)+2d0*Z/RNEW(IR)**2d0
    end if
    dvdr(ir) = dvdr(ir) + 2e0_dp*z/rnew(ir)**2e0_dp
  end do
  ! contruct LS matrix

  call spin_orbit_compl(lmax, lmmaxd, lsmh)

  ! roate LS matrix
  ncoll = 1
  if (ncoll==1) then
    call rotatematrix(lsmh, theta, phi, lmmaxd, 1)
  end if

  if (mode=='transpose') then
    do lm1 = 1, 2*lmmaxd
      do lm2 = 1, lm1 - 1
        temp = lsmh(lm2, lm1)
        lsmh(lm2, lm1) = lsmh(lm1, lm2)
        lsmh(lm1, lm2) = temp
      end do
    end do
  else if (mode=='1') then
  end if
  ! contruct prefactor of spin-orbit hamiltonian

  hsofac = 0e0_dp
  if (test('NOSOC   ') .or. z<1e-6_dp) then
    do ir = 1, irmdnew
      do lm1 = 1, 2*lmmaxd
        do lm2 = 1, 2*lmmaxd
          vnspll1(lm1, lm2, ir) = vnspll(lm1, lm2, ir)
        end do
      end do
    end do
  else
    do ir = 1, irmdnew
      rmass(ir) = 0.5e0_dp - 0.5e0_dp/c**2*((vr(ir)-real(e))-2e0_dp*z/rnew(ir) &
        )
      hsofac(ir) = socscale/(2e0_dp*rmass(ir)**2*c**2*rnew(ir))*dvdr(ir)

      ! add to potential
      do lm1 = 1, 2*lmmaxd
        do lm2 = 1, 2*lmmaxd
          vnspll1(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + &
            hsofac(ir)*lsmh(lm1, lm2)
        end do
      end do
    end do
  end if
end subroutine spinorbit_ham
