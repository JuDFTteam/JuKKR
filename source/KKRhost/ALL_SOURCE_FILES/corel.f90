subroutine corel(nsra, ipr, ip, rhoc, v, ecore, lcore, ncore, drdi, z, qc, a, &
  b, is, nspin, nr, rmax, irmd)
  ! -----------------------------------------------------------------------
  ! subroutine for core states
  ! -----------------------------------------------------------------------
  ! lmxc = lmaxcore = (0,1,2,...), .e.g, argon core : lmxc = 1
  ! krypton core : lmxc = 2
  ! kfg = configuration of core, e.g., argon core: 3300=3s,3p,0d
  ! krypton core: 4430=4s,4p,3d
  ! xenon core: 5540=5s,5p,4d
  ! -----------------------------------------------------------------------
  use :: mod_types, only: t_inc
  use :: mod_datatypes, only: dp
  implicit none
  ! .. Parameters ..
  integer :: nitmax, irnumx
  parameter (nitmax=40, irnumx=10)
  real (kind=dp) :: zero
  parameter (zero=0.0d0)
  ! ..
  ! .. Scalar Arguments ..
  real (kind=dp) :: a, b, qc, rmax, z
  integer :: ip, ipr, irmd, is, ncore, nr, nspin, nsra
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: drdi(*), ecore(*), rhoc(*), v(*)
  integer :: lcore(*)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: e, e1, e2, ediff, ei, slope, sum, tol, value, wgt
  integer :: ic, in, inuc, ir, l, lmp1, lmxc, lp1, nc, nmax, nn, nre
  logical :: vlnc
  ! ..
  ! .. Local Arrays ..
  real (kind=dp) :: f(irmd), g(irmd), rho(irmd)
  integer :: kfg(4)
  character (len=4) :: spn(2), text(5)
  ! ..
  ! .. External Subroutines ..
  external :: intcor, simp3
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: dble, real
  ! ..
  ! .. Save statement ..
  save :: spn, text
  ! ..
  ! .. Data statements ..
  data spn, text/'down', 'up  ', 's   ', 'p   ', 'd   ', 'f   ', 'g   '/
  ! ..
  vlnc = .false.
  value = 1.d-8
  slope = -1.d-8
  e2 = 50.0d0

  do ic = 1, 4
    kfg(ic) = 0
  end do
  do ic = 1, ncore
    if (lcore(ic)==0) kfg(1) = kfg(1) + 1
    if (lcore(ic)==1) kfg(2) = kfg(2) + 1
    if (lcore(ic)==2) kfg(3) = kfg(3) + 1
    if (lcore(ic)==3) kfg(4) = kfg(4) + 1
  end do
  if (kfg(2)/=0) kfg(2) = kfg(2) + 1
  if (kfg(3)/=0) kfg(3) = kfg(3) + 2
  if (kfg(4)/=0) kfg(4) = kfg(4) + 3
  lmxc = 0
  if (kfg(2)/=0) lmxc = 1
  if (kfg(3)/=0) lmxc = 2
  if (kfg(4)/=0) lmxc = 3

  tol = 1.0d-12*(z*z+1.d0)
  lmp1 = lmxc + 1
  nc = 0
  inuc = -irnumx

  do ir = 1, irmd
    rhoc(ir) = zero
    rho(ir) = zero
  end do

  do lp1 = 1, lmp1
    l = lp1 - 1
    e1 = (-5.d0-((z+1.d0)/dble(lp1))**2)*1.5d0 - 50.d0
    nmax = kfg(lp1)
    if (nmax/=0) then
      do in = lp1, nmax
        nn = in - lp1
        nc = nc + 1
        inuc = inuc + irnumx
        e = ecore(nc)
        ei = ecore(nc)
        if ((t_inc%i_write>0) .and. (ipr/=0)) write (1337, fmt=100) in, &
          text(lp1), nn, spn(is), ip, e
        call intcor(e1, e2, rho, g, f, v, value, slope, l, nn, e, sum, nre, &
          vlnc, a, b, z, rmax, nr, tol, irmd, ipr, nitmax, nsra)
        ediff = e - ei
        ecore(nc) = e
        wgt = real(l+l+1)/sum*2.d0/real(nspin)
        if ((t_inc%i_write>0) .and. (ipr/=0)) write (1337, fmt=110) ei, ediff, &
          e

        ! ---> sum up contributions to total core charge

        do ir = 2, nre
          rhoc(ir) = rhoc(ir) + rho(ir)*wgt
          rho(ir) = zero
        end do
      end do
    end if

  end do
  if (nc*irnumx>150 .or. irnumx>10) stop 'corel'

  ! ---> integrate core density to get core charge

  call simp3(rhoc, qc, 1, nr, drdi)

100 format (1x, 90('*'), /, '  n = ', i1, '  l = ', a4, '   nnode = ', i1, &
    '  spin=', a4, i5, 'th cell', '    einput = ', 1p, d16.8)
110 format (1x, '  einput =', 1p, d16.8, '   eout - ein =', 1p, d16.8, &
    '   eoutput = ', 1p, d16.8)
end subroutine corel
