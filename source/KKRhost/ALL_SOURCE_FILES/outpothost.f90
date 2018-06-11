!*==outpothost.f    processed by SPAG 6.05Rc at 13:41 on  3 Dec 2004
subroutine outpothost(alat, ins, krel, kmrot, nspin, naez, natyp, efermi, &
  bravais, rbasis, qmtet, qmphi, noq, kaoez, iqat, zat, conc, ipan, ircut, &
  solver, soc, ctl, irws, rmt, rws, rr, drdi, visp, irshift, rmrel, drdirel, &
  vtrel, btrel, lmaxd, natypd, naezd, ipand, irmd)
! **********************************************************************
! *                                                                    *
! * writes decimation potential-file  'decimate.pot' to be later used  *
! * for 2D systems with the DECIMATE option. Based on the host         *
! * potentials, the single-site matrices of the host can be calculated *
! * directly on each particular energy-mesh                            *
! *                                        v.popescu - munich, Dec 04  *
! *                                                                    *
! * Note: so far, only SPHERICAL case implemented                      *
! *                                                                    *
! **********************************************************************
  use :: mod_version_info
  implicit none
!..      
!.. Scalar arguments
  double precision :: alat, efermi
  integer :: ins, ipand, irmd, kmrot, krel, lmaxd, naez
  integer :: naezd, natyp, natypd, nspin
  character (len=10) :: solver
!..
!.. Array arguments
  double precision :: bravais(3, 3), btrel(irmd*krel+(1-krel), *), conc(*), &
    ctl(krel*lmaxd+1, *), drdi(irmd, *), drdirel(irmd*krel+(1-krel), *), &
    qmphi(*), qmtet(*), rbasis(3, *), rmrel(irmd*krel+(1-krel), *), rmt(*), &
    rr(irmd, *), rws(*), soc(krel*lmaxd+1, *), visp(irmd, *), &
    vtrel(irmd*krel+(1-krel), *), zat(*)
  integer :: ipan(*), iqat(*), ircut(0:ipand, *), irshift(*), irws(*), &
    kaoez(natypd, *), noq(naezd)
!..      
!.. Locals
  character (len=3) :: elemname(0:113)
  integer :: i, iq, ir, is
  integer :: int
  character (len=9) :: txtrel(2), txtspin(3)
!..
  data txtspin/'         ', 'spin UP  ', 'spin DOWN'/
  data txtrel/'(UP+DN)/2', '(UP-DN)/2'/
!.. 1     2     3     4    5     6     7     8     9     0
  data elemname/'Vac', 'H  ', 'He ', 'Li ', 'Be', 'B  ', 'C  ', 'N  ', 'O  ', &
    'F  ', 'Ne ', 'Na ', 'Mg ', 'Al ', 'Si', 'P  ', 'S  ', 'Cl ', 'Ar ', &
    'K  ', 'Ca ', 'Sc ', 'Ti ', 'V  ', 'Cr', 'Mn ', 'Fe ', 'Co ', 'Ni ', &
    'Cu ', 'Zn ', 'Ga ', 'Ge ', 'As ', 'Se', 'Br ', 'Kr ', 'Rb ', 'Sr ', &
    'Y  ', 'Zr ', 'Nb ', 'Mo ', 'Tc ', 'Ru', 'Rh ', 'Pd ', 'Ag ', 'Cd ', &
    'In ', 'Sn ', 'Sb ', 'Te ', 'I  ', 'Xe', 'Cs ', 'Ba ', 'La ', 'Ce ', &
    'Pr ', 'Nd ', 'Pm ', 'Sm ', 'Eu ', 'Gd', 'Tb ', 'Dy ', 'Ho ', 'Er ', &
    'Tm ', 'Yb ', 'Lu ', 'Hf ', 'Ta ', 'W ', 'Re ', 'Os ', 'Ir ', 'Pt ', &
    'Au ', 'Hg ', 'Tl ', 'Pb ', 'Bi ', 'Po', 'At ', 'Rn ', 'Fr ', 'Ra ', &
    'Ac ', 'Th ', 'Pa ', 'U  ', 'Np ', 'Pu', 'Am ', 'Cm ', 'Bk ', 'Cf ', &
    'Es ', 'Fm ', 'Md ', 'No ', 'Lr ', 'Rf', 'Db ', 'Sg ', 'Bh ', 'Hs ', &
    'Mt ', 'Uun', 'Uuu', 'Uub', 'NoE'/
!     ..
  write (1337, '(5X,A,A,/)') '< OUTPOTHOST > : ', &
    'creating decimate.pot file - host potential'

  open (37, file='decimate.pot', status='unknown')
  call version_print_header(37)
  write (37, fmt=*) 'Host structure and potential for decimation'
  write (37, fmt=100)
  write (37, fmt=140) krel, ins, nspin, kmrot
  write (37, fmt=150) naez, natyp, alat
  write (37, fmt=130) efermi
  write (37, fmt=110) bravais
! ----------------------------------------------------------------------
! here insert whatever is more needed for the structure (BZ,SYM etc),
! ref. system and so on for a full host calculation (1 iteration)
! ----------------------------------------------------------------------
  write (37, fmt=120)
  do iq = 1, naez
    write (37, fmt=160) iq, (rbasis(i,iq), i=1, 3)
  end do
  write (37, fmt=170)
  do iq = 1, naez
    write (37, fmt=180) iq, qmtet(iq), qmphi(iq), noq(iq), &
      (kaoez(i,iq), i=1, noq(iq))
  end do
  if (krel==1) write (37, 210) solver
  write (37, 190)
  do i = 1, natyp
    write (37, fmt=200) i, zat(i), iqat(i), conc(i), irws(i), ipan(i), &
      (ircut(iq,i), iq=0, ipan(i))
    if (krel==1) write (37, 220) soc(1, i), ctl(1, i)
  end do
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  do i = 1, natyp
    write (37, '(80("*"))')
    ir = int(zat(i))
    write (37, 230) i, elemname(ir), rmt(i), rws(i)
    if (krel==0) then
      write (37, '(4(A9,11X))') 'R MESH   ', 'DRDI     ', &
        (txtspin(nspin+is-1), is=1, nspin)
      do ir = 1, irws(i)
        write (37, 240) rr(ir, i), drdi(ir, i), (visp(ir,(i-1)*natyp+is), is=1 &
          , nspin)
      end do
    else
      write (37, 250) irshift(i)
      write (37, '(4(A9,11X))') 'R MESH   ', 'DRDI     ', &
        (txtrel(is), is=1, 2)
      do ir = 1, irws(i) - irshift(i)
        write (37, 240) rmrel(ir, i), drdirel(ir, i), vtrel(ir, i), &
          btrel(ir, i)
      end do
    end if
  end do

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  close (37)
!     ..
100 format ('Vectors in lattice constant units')
110 format ('BRAVAIS ', /, 3f8.4, /, 3f8.4, /, 3f8.4)
120 format ('RBASIS')
130 format ('EFERMI=', f10.6)
140 format ('KREL  =', i3, ' INS  =', i3, ' NSPIN=', i3, ' KMROT=', i3)
150 format ('NAEZ  =', i3, ' NATYP=', i3, ' ALAT =', f12.8)
160 format ('SITE  :', i3, 3f12.8)
170 format ('Magnetisation angles, occupancies, types on each site')
180 format ('SITE  :', i3, ' THETA=', f9.4, ' PHI  =', f9.4, ' NOQ  =', i3, &
    ' ITOQ :', 8i3)
190 format ('ATOMS')
200 format ('TYPE  :', i3, ' Z    =', f4.0, ' IQAT =', i3, ' CONC =', f7.4, /, &
    10x, ' IRWS =', i4, ' IPAN =', i3, ' IRCUT=', 6i4)
210 format ('SOLVER=', a10)
220 format (10x, ' SOC  =', f10.6, ' CTL  =', d13.6)
230 format ('ATOM  :', i3, 1x, a3, ': mesh and potential data', /, 'RMT   :', &
    f12.8, /, 'RWS   :', f12.8)
240 format (1p, 4d20.12)
250 format ('ISHIFT:', i3)
end subroutine
