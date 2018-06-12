subroutine gradrl(nspin, mesh, l1max, dx, rhol, rv, drdi, ipan, ipand, ircut, &
  drrl, ddrrl, drrul, ddrrul, irmd, lmpotd)
  ! ------------------------------------------------------------------
  ! gradient of rl with rl defined by charge density=sum(rl*ylm).
  ! mesh,l1max: max of mesh and l+1.
  ! IRMD,LMPOTD: maxima of corresponding dimension parameters.
  ! drrl=d(rl)/dr, ddrrl=d(drrl)/dr, drrul=d(rl-up)/dr,
  ! ztal: zeta for each l-component necessary to get down-components.
  ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------
  use :: mod_types, only: t_inc
  use :: mod_datatypes, only: dp
  implicit none
  ! .. Parameters ..
  real (kind=dp) :: zero, zero1
  parameter (zero=0.d0, zero1=1.d-12)
  ! ..
  ! .. Scalar Arguments ..
  real (kind=dp) :: dx
  integer :: ipan, ipand, irmd, l1max, lmpotd, mesh, nspin
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: ddrrl(irmd, lmpotd), ddrrul(irmd, lmpotd), drdi(irmd), &
    drrl(irmd, lmpotd), drrul(irmd, lmpotd), rhol(irmd, 2, lmpotd), rv(irmd)
  integer :: ircut(0:ipand)
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: chgden, pi, r2, s4, spiden
  integer :: i1, ien, ip, ir, ist, llmax
  ! ..
  ! .. Local Arrays ..
  real (kind=dp) :: drdi2(irmd), rl1(irmd), rl1udm(irmd), ztal(irmd)
  ! ..
  ! .. External Subroutines ..
  external :: gradr
  ! ..
  ! .. Intrinsic Functions ..
  intrinsic :: abs, acos, sqrt
  ! ..
  ! ------------------------------------------------------------------
  pi = cos(-1.d0)
  s4 = sqrt(4.d0*pi)
  llmax = l1max*l1max

  do ip = 1, ipan
    ist = ircut(ip-1) + 1
    ien = ircut(ip)
    if (t_inc%i_write>0) write (1337, fmt=130) ip, ist, ien
    if (ip==1) then
      do ir = ist, ien
        drdi2(ir) = dx
      end do
    else
      do ir = ist, ien
        drdi2(ir) = zero
      end do
    end if
  end do


  do i1 = 1, llmax

    if (nspin==1) go to 100

    do ir = 2, mesh
      r2 = rv(ir)*rv(ir)
      chgden = rhol(ir, 1, i1) + rhol(ir, 2, i1)
      spiden = rhol(ir, 2, i1) - rhol(ir, 1, i1)
      if (abs(chgden)>=zero1) then
        rl1(ir) = chgden
        ztal(ir) = spiden/chgden
      else
        rl1(ir) = zero
        ztal(ir) = zero
      end if
    end do

    go to 110


100 continue
    do ir = 2, mesh
      r2 = rv(ir)*rv(ir)
      rl1(ir) = rhol(ir, 1, i1) + rhol(ir, 2, i1)
      ztal(ir) = zero
    end do

110 continue

    rl1(1) = rl1(2)
    ztal(1) = ztal(2)


    do ip = 1, ipan
      ist = ircut(ip-1) + 1
      ien = ircut(ip)

      call gradr(nspin, ist, ien, 1.d0, drdi, drdi2, rl1, ztal, drrl(1,i1), &
        ddrrl(1,i1), drrul(1,i1), ddrrul(1,i1), rl1udm, irmd)

      if (ip==1) then
        do ir = 1, 4
          drrl(ir, i1) = drrl(5, i1)
          ddrrl(ir, i1) = ddrrl(5, i1)
          drrul(ir, i1) = drrul(5, i1)
          ddrrul(ir, i1) = ddrrul(5, i1)
        end do
      end if

      if (nspin==1) then
        do ir = ist, ien
          drrul(ir, i1) = drrl(ir, i1)/2.d0
          ddrrul(ir, i1) = ddrrl(ir, i1)/2.d0
        end do
      end if

    end do

  end do

  return
120 format (1x, ' l1max=', i5, ' mesh=', i5, 'nspi=', i5, ' ipan=', i5)
130 format (1x, '  ip ist ien', 3i5)
end subroutine gradrl
