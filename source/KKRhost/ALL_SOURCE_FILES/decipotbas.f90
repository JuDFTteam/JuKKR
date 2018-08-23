module mod_decipotbas

contains

subroutine decipotbas(ihost, iqoff, itoff, nq, nt, rbasis, qmtet, qmphi, noq, &
  kaoez, zat, iqat, conc, irws, ipan, ircut, rr, drdi, visp, nspin, krel, &
  solver, socscl, cscl, vtrel, btrel, irmd, ipand, nembd1, ntmax, nspind, &
  lmaxd)
  use :: mod_datatypes, only: dp
  ! **********************************************************************
  ! *                                                                    *
  ! * reads in the potential data for the host atoms from the potential  *
  ! * file 'decimate.pot'                                                *
  ! *                                        v.popescu - munich, Dec 04  *
  ! *                                                                    *
  ! * Note: so far, only  SPHERICAL case implemented                     *
  ! *                                                                    *
  ! **********************************************************************
  implicit none
  ! ..
  ! .. Arguments
  integer :: ihost, ipand, iqoff, irmd, itoff, krel, lmaxd, nembd1, nq, nspin, &
    nspind, nt, ntmax
  character (len=10) :: solver
  real (kind=dp) :: btrel(irmd*krel+(1-krel), ntmax), conc(ntmax), &
    cscl(krel*lmaxd+1, krel*ntmax+(1-krel)), drdi(irmd, ntmax), qmphi(nembd1), &
    qmtet(nembd1), rbasis(3, nembd1), rr(irmd, ntmax), &
    socscl(krel*lmaxd+1, krel*ntmax+(1-krel)), visp(irmd, ntmax*nspind), &
    vtrel(irmd*krel+(1-krel), ntmax), zat(ntmax)
  integer :: ipan(ntmax), iqat(nembd1, ntmax), ircut(0:ipand, ntmax), &
    irws(ntmax), kaoez(nembd1, nembd1), noq(nembd1)
  ! ..
  ! .. Locals
  integer :: i, ih, ihf, il, ipot1, ipot2
  integer :: nint
  real (kind=dp) :: rmt(ntmax), rws(ntmax)
  character (len=3) :: txtt(nt)
  ! ......................................................................
  ! --> read basis

  do ih = 1, nq
    ihf = ih + iqoff
    read (36+ihost, 100) il, (rbasis(i,ihf), i=1, 3)
    if (ih/=il) stop ' Inconsistent data '
    write (1337, 100) il, (rbasis(i,ihf), i=1, 3)
  end do
  read (36+ihost, *)
  write (1337, 120)
  do ih = 1, nq
    ihf = ih + iqoff
    read (36+ihost, fmt=110) il, qmtet(ihf), qmphi(ihf), noq(ihf), &
      (kaoez(i,ihf), i=1, noq(ihf))
    if (ih/=il) stop ' Inconsistent data '
    write (1337, 130) ih, qmtet(ihf), qmphi(ihf), noq(ihf), &
      (kaoez(i,ihf), i=1, noq(ihf))
  end do
  write (1337, 140)
  if (krel==1) read (36+ihost, '(7X,A10)') solver
  read (36+ihost, *)

  ! --> read atoms

  do ih = 1, nt
    ihf = ih + itoff
    read (36+ihost, 150) il, zat(ihf), iqat(1, ihf), conc(ihf), irws(ihf), &
      ipan(ihf), (ircut(i,ihf), i=0, ipan(ihf))
    if (ih/=il) stop ' Inconsistent data '

    if (krel==1) then
      read (36+ihost, 160) socscl(1, ihf), cscl(1, ihf)
      do il = 2, lmaxd + 1
        socscl(il, ihf) = socscl(1, ihf)
        cscl(il, ihf) = cscl(1, ihf)
      end do
    end if

  end do
  do ih = 1, nt
    ihf = ih + itoff
    ipot1 = (ihf-1)*nspin + 1
    ipot2 = ipot1 + 1
    read (36+ihost, *)
    read (36+ihost, 170) il, txtt(ih), rmt(ihf), rws(ihf)
    if (ih/=il) stop ' Inconsistent data '
    write (1337, 180) ih, txtt(ih), nint(zat(ihf)), conc(ihf), irws(ihf), &
      rws(ihf)
    if (krel==0) then
      read (36+ihost, *)
      do i = 1, irws(ihf)
        read (36+ihost, 190) rr(i, ihf), drdi(i, ihf), &
          (visp(i,il), il=ipot1, ipot2)
      end do
    else
      read (36+ihost, '(7X,I3)') il
      irws(ihf) = irws(ihf) - il
      ircut(ipan(ihf), ihf) = irws(ihf)
      read (36+ihost, *)
      do i = 1, irws(ihf)
        read (36+ihost, 190) rr(i, ihf), drdi(i, ihf), vtrel(i, ihf), &
          btrel(i, ihf)
      end do
    end if
  end do

100 format (9x, i3, 3f12.8)
110 format (7x, i3, 2(7x,f9.4), 7x, i3, 7x, 8i3)
120 format (9x, 39('-'), /, 9x, '   THETA   ', '   PHI   ', 'OCC', ' IT')
130 format (9x, i3, 2(f9.4), i3, 8i3)
140 format (9x, 39('-'), /, 10x, 'ATOMS', /, 15x, 'Z   CONC  IWS    RWS')
150 format (7x, i3, 7x, f4.0, 7x, i3, 7x, f7.4, /, 17x, i4, 7x, i3, 7x, 6i4)
160 format (17x, f10.6, 7x, d13.6)
170 format (7x, i3, 1x, a3, 2(/,7x,f12.8))
180 format (9x, i3, 1x, a3, i3, f7.4, i4, f10.6)
190 format (1p, 4d20.12)
end subroutine decipotbas

end module mod_decipotbas
