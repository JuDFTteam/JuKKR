subroutine spinorbit_ham(lmax, lmmaxd, vins, rnew, e, z, c, socscale, nspin, &
  lmpotd, theta, phi, ipan_intervall, rpan_intervall, npan_tot, ncheb, &
  irmdnew, nrmaxd, vnspll, vnspll1, mode)
  implicit none

  integer :: lmax, lmmaxd, nspin, npan_tot, ncheb, irmdnew, nrmaxd
  integer :: lmpotd
  double precision :: c, z
  double complex :: e
  double precision :: socscale
  double precision :: vins(irmdnew, lmpotd, nspin), rnew(nrmaxd), &
    rpan_intervall(0:npan_tot)
  double complex :: vnspll(2*lmmaxd, 2*lmmaxd, irmdnew)
  double complex :: vnspll1(2*lmmaxd, 2*lmmaxd, irmdnew)
  integer :: ipan_intervall(0:npan_tot)
  double precision :: vr(irmdnew), dvdr(irmdnew)
  double precision :: rmass(irmdnew), hsofac(irmdnew)
  double precision :: rnucl, atn, widthfac, phi, theta
  integer :: ir, ip, lm1, lm2, ispin, irmin, irmax, ncoll
  double complex :: lsmh(2*lmmaxd, 2*lmmaxd), temp
  double precision :: clambdacinv(0:ncheb, 0:ncheb)
  character (len=*) :: mode
  logical :: test, opt
  external :: test, opt

  vnspll1 = (0d0, 0d0)
  vr = 0d0
  do ispin = 1, nspin
    do ir = 1, ipan_intervall(npan_tot)
      vr(ir) = vr(ir) + vins(ir, 1, ispin)/nspin
    end do
  end do
! derivative of potential
  dvdr = 0d0
  call getclambdacinv(ncheb, clambdacinv)
  do ip = 1, npan_tot
    irmin = ipan_intervall(ip-1) + 1
    irmax = ipan_intervall(ip)
    widthfac = 2d0/(rpan_intervall(ip)-rpan_intervall(ip-1))
    call dgemv('N', ncheb+1, ncheb+1, 1d0, clambdacinv, ncheb+1, &
      vr(irmin:irmax), 1, 0d0, dvdr(irmin:irmax), 1)
    dvdr(irmin:irmax) = dvdr(irmin:irmax)*widthfac
  end do
! core potential
  if (z>24d0) then
    atn = -16.1532921 + 2.70335346*z
  else
    atn = 0.03467714 + 2.04820786*z
  end if
  rnucl = 1.2d0/0.529177d0*atn**(1./3d0)*1.d-5

  do ir = 1, ipan_intervall(npan_tot)
    if (rnew(ir)<=rnucl) then
!        DVDR(IR)=DVDR(IR)+2d0*Z*RNEW(IR)/RNUCL**3d0
    else
!        DVDR(IR)=DVDR(IR)+2d0*Z/RNEW(IR)**2d0
    end if
    dvdr(ir) = dvdr(ir) + 2d0*z/rnew(ir)**2d0
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

  hsofac = 0d0
  if (test('NOSOC   ') .or. z<1d-6) then
    do ir = 1, irmdnew
      do lm1 = 1, 2*lmmaxd
        do lm2 = 1, 2*lmmaxd
          vnspll1(lm1, lm2, ir) = vnspll(lm1, lm2, ir)
        end do
      end do
    end do
  else
    do ir = 1, irmdnew
      rmass(ir) = 0.5d0 - 0.5d0/c**2*((vr(ir)-real(e))-2d0*z/rnew(ir))
      hsofac(ir) = socscale/(2d0*rmass(ir)**2*c**2*rnew(ir))*dvdr(ir)

! add to potential
      do lm1 = 1, 2*lmmaxd
        do lm2 = 1, 2*lmmaxd
          vnspll1(lm1, lm2, ir) = vnspll(lm1, lm2, ir) + &
            hsofac(ir)*lsmh(lm1, lm2)
        end do
      end do
    end do
  end if
end subroutine


