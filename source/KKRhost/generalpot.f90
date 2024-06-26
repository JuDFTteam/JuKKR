!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_generalpot
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Writes potential out interpolated to a generalized radial mesh
  !> Author: 
  !> Category: KKRhost, potential, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> The subroutine writes out the potential cards
  !> in a standard r-mesh that can be read in and
  !> interpolated to a different r-mesh from subroutine
  !> start No shape function information is needed
  !> and all nessecery data are stored in the potential
  !> card.
  !>                                      ver. 18.5.2000
  !>
  !> @warning Using the general potential and consecutive read-in does not result in a converged
  !> potential since the convolution with the shape functions is then done twice 
  !> (which does only not make a difference in the limit l_max-> infinity).
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine generalpot(ifile, natps, natyp, nspin, z, alat, rmt, rmtnew, rws, r, drdi, vm2z, irws, a, b, ins, irns, lpot, vins, qbound, irc, kshape, efermi, vbc, ecore, lcore, &
    ncore, lmpotd, irmd, irmind)

    use mod_splint, only: splint_real
    use mod_spline, only: spline_real

    implicit none

    integer :: lmpotd, irmd, irmind
    real (kind=dp) :: alat, qbound
    integer :: ifile, ins, kshape, lpot, natps, natyp, nspin

    real (kind=dp) :: a(*), b(*), drdi(irmd, *), ecore(20, *), efermi, r(irmd, *), rmt(*), rmtnew(*), rws(*), vbc(2), vins(irmind:irmd, lmpotd, *), vm2z(irmd, *), z(*)
    integer :: irc(*), irns(*), irws(*), lcore(20, *), ncore(*)

    real (kind=dp) :: a1, b1, rmax, rmt1, rmtnw1, rv, sum, z1, parsum, parsumderiv, r0, rinter, dr, maxa
    integer :: i, icore, ih, ip, ir, irmin, irns1, is, isave, j, lm, lmnr, lmpot, ncore1, nr, nz1, nr_u, irmin_u, irns_u, imt1, lm1, irnstot

    real (kind=dp) :: dradi(irmd), ecore1(20), ra(irmd), vm2za(irmd), rr_u(irmd)
    real (kind=dp) :: vm2zb(irmd), vm2z_u(irmd), vins_u(irmind:irmd, lmpotd), vinsa(irmind:irmd, lmpotd), vinsb(irmind:irmd, lmpotd)
    integer :: lcore1(20)
    character (len=4) :: elemname(0:113)

    intrinsic :: sqrt
    ! 1      2      3      4      5      6      7      8      9
    data elemname/'VAC', 'H   ', 'He  ', 'Li  ', 'Be  ', 'B   ', 'C   ', 'N   ', 'O   ', 'F   ', 'Ne  ', 'Na  ', 'Mg  ', 'Al  ', 'Si  ', 'P   ', 'S   ', 'Cl  ', 'Ar  ', 'K   ', &
      'Ca  ', 'Sc  ', 'Ti  ', 'V   ', 'Cr  ', 'Mn  ', 'Fe  ', 'Co  ', 'Ni  ', 'Cu  ', 'Zn  ', 'Ga  ', 'Ge  ', 'As  ', 'Se  ', 'Br  ', 'Kr  ', 'Rb  ', 'Sr  ', 'Y   ', 'Zr  ', &
      'Nb  ', 'Mo  ', 'Tc  ', 'Ru  ', 'Rh  ', 'Pd  ', 'Ag  ', 'Cd  ', 'In  ', 'Sn  ', 'Sb  ', 'Te  ', 'I   ', 'Xe  ', 'Cs  ', 'Ba  ', 'La  ', 'Ce  ', 'Pr  ', 'Nd  ', 'Pm  ', &
      'Sm  ', 'Eu  ', 'Gd  ', 'Tb  ', 'Dy  ', 'Ho  ', 'Er  ', 'Tm  ', 'Yb  ', 'Lu  ', 'Hf  ', 'Ta  ', 'W   ', 'Re  ', 'Os  ', 'Ir  ', 'Pt  ', 'Au  ', 'Hg  ', 'Tl  ', 'Pb  ', &
      'Bi  ', 'Po  ', 'At  ', 'Rn  ', 'Fr  ', 'Ra  ', 'Ac  ', 'Th  ', 'Pa  ', 'U   ', 'Np  ', 'Pu  ', 'Am  ', 'Cm  ', 'Bk  ', 'Cf  ', 'Es  ', 'Fm  ', 'Md  ', 'No  ', 'Lr  ', &
      'Rf  ', 'Db  ', 'Sg  ', 'Bh  ', 'Hs  ', 'Mt  ', 'Uun ', 'Uuu ', 'Uub ', 'NoE '/


    isave = 1
    lmpot = (lpot+1)*(lpot+1)

    do ih = 1, natyp
      do is = 1, nspin
        do lm = 1, lmpotd
          do ir = irmind, irmd
            vinsa(ir, lm) = 0.e0_dp
            vinsb(ir, lm) = 0.e0_dp
          end do
        end do
        ip = nspin*(ih-1) + is
        rmt1 = rmt(ih)
        rmtnw1 = rmtnew(ih)
        z1 = z(ih)
        rmax = rws(ih)
        if (kshape==0) then
          nr = irws(ih)
          irns1 = 0
        else
          nr = irc(ih)
          irns1 = irns(ih)
        end if

        irmin = nr - irns1
        a1 = a(ih)
        b1 = b(ih)
        ncore1 = ncore(ip)

        do j = 1, nr
          ra(j) = r(j, ih)
          dradi(j) = drdi(j, ih)
          vm2za(j) = vm2z(j, ip)
        end do
        do lm1 = 1, lmpot
          do j = irmind, irmd
            vinsa(j, lm1) = vins(j, lm1, ip)
          end do
        end do

        if (ncore1>=1) then

          do j = 1, ncore1
            lcore1(j) = lcore(j, ip)
            ecore1(j) = ecore(j, ip)
          end do
        end if

        ! Generate uniform mesh RUNI

        nr_u = nr
        irns_u = irns1
        irmin_u = nr_u
        if (ins>0) irmin_u = nr_u - irns_u

        if (ins==0) then
          do i = 1, nr_u
            rr_u(i) = ra(i)
          end do
          imt1 = 0
        else
          imt1 = nint(log(rmtnw1/b1+1.0e0_dp)/a1) + 1
          do i = 1, imt1
            rr_u(i) = ra(i)
          end do
          rinter = rmax - rmtnw1
          dr = rinter/real(nr-imt1, kind=dp)
          do i = 1, nr - imt1
            rr_u(imt1+i) = rr_u(imt1) + dr*real(i, kind=dp)
          end do
          call doubleraus1(nr, irmin, lmpot, ra, dradi, vm2za, vinsa, irmd, irmind, lmpotd)

          ! After this sub the arrays are rearanged and nr is not
          ! the same anymore in the case of FP. If ins.eq.0 there is
          ! no nead for doubleraus1. IRMIN should remain the same

        end if
        ! ----------------------------------------------------------------
        ! Now the new mesh is generated

        ! test
        ! write(6,*) nr_u,imt1,irns_u
        ! do i=1,nr_u
        ! write(6,*) i,ra(i),rr_u(i)
        ! end do
        ! test

        maxa = 1.e35_dp
        call spline_real(irmd, ra, vm2za, nr, maxa, maxa, vm2zb)
        if (ins>0) then
          do lm1 = 1, lmpot
            irnstot = nr - irmin   ! nr has changed irmin is the same
            ! write(6,*) ' Testing ',nr,irmin,irnstot,irmind
            call spline_real(irmd-irmind, ra(irmind), vinsa(irmind,lm1), irnstot, maxa, maxa, vinsb(irmind,lm1))
          end do                   ! LM1
        end if

        ! OK with spline

        do ir = 1, nr_u
          r0 = rr_u(ir)
          call splint_real(ra, vm2za, vm2zb, nr, r0, parsum, parsumderiv)
          vm2z_u(ir) = parsum
        end do
        if (ins>0) then
          ! IRNSTOT = NR_U - IRMIN_U
          do lm1 = 1, lmpot
            do ir = irmin_u, nr_u
              r0 = rr_u(ir)
              call splint_real(ra(irmind), vinsa(irmind,lm1), vinsb(irmind,lm1), irnstot, r0, parsum, parsumderiv)
              vins_u(ir, lm1) = parsum
            end do
          end do
        end if
        ! write(6,*) ' All interpolation ok now write'
        ! --------------------------------------------------------------
        write (ifile, fmt=100)
        nz1 = nint(z1)
        if (nspin==1) then
          write (ifile, fmt=110) elemname(nz1), z1
        else if (is==1) then
          write (ifile, fmt=130) elemname(nz1), z1
        else if (is==2) then
          write (ifile, fmt=120) elemname(nz1), z1
        end if
        write (ifile, fmt=140)
        ! write (ifile,*) ALAT,RMAX,RMTNW1,RMT1
        write (ifile, fmt=150) alat, rmax, rmtnw1, rmt1
        write (ifile, fmt=160) nr_u, imt1, irns1
        write (ifile, fmt=170) a1, b1
        write (ifile, fmt=180) efermi, vbc(is)
        write (ifile, fmt=190) ncore1, lmpot
        if (ncore1>=1) write (ifile, fmt=240)(lcore1(icore), ecore1(icore), icore=1, ncore1)

        if (ins==0 .or. (ih<natps .and. ins<=2)) then

          ! ---  >       store only the spherically averaged potential
          ! (in mt or as - case)
          ! this is done always for the host

          write (ifile, fmt=260)(vm2z_u(ir), ir=1, nr_u)
        else

          ! ---  >     store the full potential , but the non spherical
          ! contribution
          ! only from irns1 up to irws1 ;
          ! remember that the lm = 1 contribution is multiplied
          ! by a factor 1/sqrt(4 pi)

          write (ifile, fmt=270) nr_u, irns1, lmpot, isave
          write (ifile, fmt=280)(vm2z_u(ir), ir=1, nr_u)
          if (lpot>0) then
            lmnr = 1
            do lm = 2, lmpot
              sum = 0.0e0_dp
              do ir = irmin, nr_u
                rv = vins_u(ir, lm)*rr_u(ir)
                sum = sum + rv*rv*dradi(ir)
              end do

              if (sqrt(sum)>qbound) then
                lmnr = lmnr + 1
                write (ifile, fmt=270) lm
                write (ifile, fmt=280)(vins_u(ir,lm), ir=irmin, nr_u)
              end if

            end do

            ! ---  >         write a one to mark the end

            if (lmnr<lmpot) write (ifile, fmt=270) isave
          end if

        end if

      end do
    end do


100 format (' GENERAL POTENTIAL MESH             exc:')
110 format ('#  ', a4, 'POTENTIAL             Z = ', f8.3)
120 format ('#  ', a4, 'POTENTIAL SPIN UP     Z=  ', f8.3)
130 format ('#  ', a4, 'POTENTIAL SPIN DOWN   Z=  ', f8.3)
140 format ('#')
150 format (4f12.8, '   # alat, rmax, rmaxlog, rmt')
160 format (1p, 3i6, 31x, '  # IRWS, IRMT, IRNS ')
170 format (2d15.8, 19x, '  # A , B ')
180 format (3f12.8, 13x, '  # Ef, vbc ')
190 format (1p, 2i5, 39x, '  # NCORE, LMPOT')
200 format (7a4, 6x, '  exc:', a24, 3x, a10)
210 format (3f12.8)
220 format (f10.5, /, f10.5, 2f15.10)
230 format (i3, /, 2d15.8, /, 2i2)
240 format (i5, 1p, d20.11)
250 format (1p, 2d15.6, 1p, d15.8)
260 format (1p, 4d20.12)
270 format (10i5)
280 format (1p, 4d20.13)
  end subroutine generalpot


  !-------------------------------------------------------------------------------
  !> Summary: Remove doubles points in radial mesh
  !> Author: 
  !> Category: KKRhost, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Gets rid of the double-points in the radial mesh, i.e. the points
  !> where RR(I) = RR(I-1). Returns the "new" mesh in the same array,
  !> and rearranges accordingly the WAVEF defined at the same mesh.
  !> IRMAX is also altered to the new value.
  !-------------------------------------------------------------------------------
  subroutine doubleraus1(irmax, irmin, lmpot, rr, drdi, vpot, vins, irmd, irmind, lmpotd)

    implicit none
    integer :: irmd, lmpotd, irmind
    integer :: ncountmax
    parameter (ncountmax=500)
    ! Input and output:
    integer :: irmax
    real (kind=dp) :: rr(irmd), drdi(irmd), vpot(irmd), vins(irmind:irmd, lmpotd)
    ! Inside:
    integer :: ir, icount, nc, ncount, lmpot, irmin
    integer :: lm1
    integer :: idouble(ncountmax)

    ! Find double-points:
    ncount = 0
    do ir = 2, irmax
      if ((rr(ir)-rr(ir-1))<1.e-20_dp) then
        ncount = ncount + 1
        idouble(ncount) = ir
      end if
    end do
    if (ncount+1>ncountmax) stop 'DOUBLERAUS2: Too many double-points.'
    idouble(ncount+1) = irmax + 1  ! To be used below.

    ! Rearrange the arrays.
    do icount = 1, ncount
      do ir = idouble(icount) - icount + 1, idouble(icount+1) - icount
        if ((ir+icount)<=irmax) then
          rr(ir) = rr(ir+icount)
          drdi(ir) = drdi(ir+icount)
          vpot(ir) = vpot(ir+icount)
        end if
      end do
    end do
    irmax = irmax - ncount
    ncount = 0
    do nc = 1, ncountmax
      idouble(nc) = 0
    end do

    do ir = irmin, irmax
      if ((rr(ir)-rr(ir-1))<1.e-20_dp) then
        ncount = ncount + 1
        idouble(ncount) = ir
      end if
    end do
    ! Rearrange the arrays.
    do icount = 1, ncount
      do ir = idouble(icount) - icount + 1, idouble(icount+1) - icount
        do lm1 = 1, lmpot
          vins(ir, lm1) = vins(ir+icount, lm1)
        end do
      end do
    end do
    return
  end subroutine doubleraus1

end module mod_generalpot
