module mod_spinorbit_ham

contains

!< subroutine that constructs SOC potential from radial derivative of vins and adds this to vnspll (output is vnspll1=vnspll+V_SOC)
subroutine spinorbit_ham(lmax, lmmaxd, vins, rnew, eryd, Zat, cvlight, socscale, nspin, &
  lmpotd, theta, phi, ipan_intervall, rpan_intervall, npan_tot, ncheb, &
  irmdnew, nrmaxd, vnspll, vnspll1, mode)
  use mod_datatypes, only: dp
  use mod_cheb, only: getclambdacinv
  use mod_spin_orbit_compl
  use mod_rotatespinframe, only: rotatematrix
  implicit none

  !inputs
  integer, intent(in) :: lmax                              ! < l_max cutoff
  integer, intent(in) :: lmmaxd                            ! < (l_max+1)^2  maximal number in combined L=(l,m) index (L_max)
  integer, intent(in) :: nspin                             ! < number of spin channels
  integer, intent(in) :: lmpotd                            ! < L_max cutoff of potential (from Gaunt coefficients <= 4 l_max)
  integer, intent(in) :: npan_tot                          ! < total number of Chebychev panels of radial mesh
  integer, intent(in) :: ncheb                             ! < number of Chebychev polynomials (radial points per panel)
  integer, intent(in) :: irmdnew                           ! < 
  integer, intent(in) :: nrmaxd                            ! < maximal number of radial points (NPAN_TOT*NCHEB)
  real (kind=dp), intent(in) :: cvlight                    ! < speed of light
  real (kind=dp), intent(in) :: Zat                        ! < atom charge
  complex (kind=dp), intent(in) :: eryd                    ! < complex energy
  real (kind=dp), intent(in) :: socscale                   ! < scaling factor for SOC strength
  real (kind=dp), intent(in) :: vins(irmdnew,lmpotd,nspin) ! < non-sperical input potential in (l,m) basis, separately spin-polarized
  real (kind=dp), intent(in) :: rnew(nrmaxd)               ! < radial points of Chebychev mesh
  real (kind=dp), intent(in) :: rpan_intervall(0:npan_tot) ! < 
  integer, intent(in) :: ipan_intervall(0:npan_tot)        ! < 
  complex (kind=dp), intent(in) :: vnspll(2*lmmaxd, 2*lmmaxd, irmdnew) ! < input potential in (l,m,s) basis
  character (len=*), intent(in) :: mode                    ! < either '1' or 'transpose', depending whether SOC potential is constructed for right or left solution

  !outputs
  complex (kind=dp), intent(out) :: vnspll1(2*lmmaxd, 2*lmmaxd, irmdnew) ! < output potential (sum of input + V_SOC) in (l,m,s) basis

  !locals
  real (kind=dp) :: vr(irmdnew)
  real (kind=dp) :: dvdr(irmdnew)
  real (kind=dp) :: rmass(irmdnew)
  real (kind=dp) :: hsofac(irmdnew)
  !real (kind=dp) :: rnucl, atn
  real (kind=dp) :: widthfac, phi, theta
  integer :: ir, ip, lm1, lm2, ispin, irmin, irmax, ncoll
  complex (kind=dp) :: lsmh(2*lmmaxd, 2*lmmaxd), temp
  real (kind=dp) :: clambdacinv(0:ncheb, 0:ncheb)
  logical :: test, opt
  external :: test, opt

  ! fill radial potential (used to construct radial derivative of potential)
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
  

  !! core potential
  !if (Zat>24e0_dp) then
  !  atn = -16.1532921_dp + 2.70335346_dp*Zat
  !else
  !  atn = 0.03467714_dp + 2.04820786_dp*Zat
  !end if
  !rnucl = 1.2e0_dp/0.529177e0_dp*atn**(1._dp/3e0_dp)*1.e-5_dp

  ! add core potential
  do ir = 1, ipan_intervall(npan_tot)
    !if (rnew(ir)<=rnucl) then
      ! DVDR(IR)=DVDR(IR)+2d0*Zat*RNEW(IR)/RNUCL**3d0
    !else
      ! DVDR(IR)=DVDR(IR)+2d0*Zat/RNEW(IR)**2d0
    !end if
    dvdr(ir) = dvdr(ir) + 2e0_dp*Zat/rnew(ir)**2e0_dp
  end do


  ! contruct LS matrix (output: lsmh)
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
    ! do nothing
  end if


  ! contruct prefactor of spin-orbit hamiltonian
  hsofac = 0e0_dp
  vnspll1 = (0e0_dp, 0e0_dp)
  if (test('NOSOC   ') .or. Zat<1e-6_dp) then
    vnspll1(1:2*lmmaxd, 1:2*lmmaxd, 1:irmdnew) = vnspll(1:2*lmmaxd, 1:2*lmmaxd, 1:irmdnew)
  else
    do ir = 1, irmdnew
      rmass(ir) = 0.5e0_dp - 0.5e0_dp/cvlight**2*((vr(ir)-real(eryd))-2e0_dp*Zat/rnew(ir))
      hsofac(ir) = socscale/(2e0_dp*rmass(ir)**2*cvlight**2*rnew(ir))*dvdr(ir)

      ! and add to potential
      vnspll1(1:2*lmmaxd, 1:2*lmmaxd, ir) = vnspll(1:2*lmmaxd, 1:2*lmmaxd, ir) + hsofac(ir)*lsmh(1:2*lmmaxd, 1:2*lmmaxd)
    end do
  end if

  ! now output vnspll1 is sum of input and SOC potential (eventually scaled with SOCFAC)

end subroutine spinorbit_ham

end module mod_spinorbit_ham
