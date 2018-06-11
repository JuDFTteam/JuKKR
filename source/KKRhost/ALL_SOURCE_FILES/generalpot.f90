! **********************************************************************
    Subroutine generalpot(ifile, natps, natyp, nspin, z, alat, rmt, rmtnew, &
      rws, r, drdi, vm2z, irws, a, b, ins, irns, lpot, vins, qbound, irc, &
      kshape, efermi, vbc, ecore, lcore, ncore, lmpotd, irmd, irmind)
      Use mod_datatypes, Only: dp
! **************************************************
! * The subroutine writes out the potential cards
! * in a standard r-mesh that can be read in and
! * interpolated to a different r-mesh from subroutine
! * start No shape function information is needed
! * and all nessecery data are stored in the potential
! * card.
! *                                      ver. 18.5.2000
! ***************************************************
!     ..
      Implicit None
!..
!.. Scalar Arguments ..
      Integer :: lmpotd, irmd, irmind
      Real (Kind=dp) :: alat, qbound
      Integer :: ifile, ins, kshape, lpot, natps, natyp, nspin
!..
!.. Array Arguments ..
      Real (Kind=dp) :: a(*), b(*), drdi(irmd, *), ecore(20, *), efermi, &
        r(irmd, *), rmt(*), rmtnew(*), rws(*), vbc(2), &
        vins(irmind:irmd, lmpotd, *), vm2z(irmd, *), z(*)
      Integer :: irc(*), irns(*), irws(*), lcore(20, *), ncore(*)
!..
!.. Local Scalars ..
      Real (Kind=dp) :: a1, b1, rmax, rmt1, rmtnw1, rv, sum, z1, parsum, &
        parsumderiv, r0, rinter, dr, maxa
      Integer :: i, icore, ih, ip, ir, irmin, irns1, is, isave, j, lm, lmnr, &
        lmpot, ncore1, nr, nz1, nr_u, irmin_u, irns_u, imt1, lm1, irnstot
!..
!.. Local Arrays ..
      Real (Kind=dp) :: dradi(irmd), ecore1(20), ra(irmd), vm2za(irmd), &
        rr_u(irmd), drdi_u(irmd)
      Real (Kind=dp) :: vm2zb(irmd), vm2z_u(irmd), vins_u(irmind:irmd, lmpotd) &
        , vinsa(irmind:irmd, lmpotd), vinsb(irmind:irmd, lmpotd)
      Integer :: lcore1(20)
      Character (Len=4) :: elemname(0:113)
!..
!.. Intrinsic Functions ..
      Intrinsic :: sqrt
!        1      2      3      4      5      6      7      8      9    
      Data elemname/'VAC', 'H   ', 'He  ', 'Li  ', 'Be  ', 'B   ', 'C   ', &
        'N   ', 'O   ', 'F   ', 'Ne  ', 'Na  ', 'Mg  ', 'Al  ', 'Si  ', &
        'P   ', 'S   ', 'Cl  ', 'Ar  ', 'K   ', 'Ca  ', 'Sc  ', 'Ti  ', &
        'V   ', 'Cr  ', 'Mn  ', 'Fe  ', 'Co  ', 'Ni  ', 'Cu  ', 'Zn  ', &
        'Ga  ', 'Ge  ', 'As  ', 'Se  ', 'Br  ', 'Kr  ', 'Rb  ', 'Sr  ', &
        'Y   ', 'Zr  ', 'Nb  ', 'Mo  ', 'Tc  ', 'Ru  ', 'Rh  ', 'Pd  ', &
        'Ag  ', 'Cd  ', 'In  ', 'Sn  ', 'Sb  ', 'Te  ', 'I   ', 'Xe  ', &
        'Cs  ', 'Ba  ', 'La  ', 'Ce  ', 'Pr  ', 'Nd  ', 'Pm  ', 'Sm  ', &
        'Eu  ', 'Gd  ', 'Tb  ', 'Dy  ', 'Ho  ', 'Er  ', 'Tm  ', 'Yb  ', &
        'Lu  ', 'Hf  ', 'Ta  ', 'W   ', 'Re  ', 'Os  ', 'Ir  ', 'Pt  ', &
        'Au  ', 'Hg  ', 'Tl  ', 'Pb  ', 'Bi  ', 'Po  ', 'At  ', 'Rn  ', &
        'Fr  ', 'Ra  ', 'Ac  ', 'Th  ', 'Pa  ', 'U   ', 'Np  ', 'Pu  ', &
        'Am  ', 'Cm  ', 'Bk  ', 'Cf  ', 'Es  ', 'Fm  ', 'Md  ', 'No  ', &
        'Lr  ', 'Rf  ', 'Db  ', 'Sg  ', 'Bh  ', 'Hs  ', 'Mt  ', 'Uun ', &
        'Uuu ', 'Uub ', 'NoE '/

      isave = 1
      lmpot = (lpot+1)*(lpot+1)

      Do ih = 1, natyp
        Do is = 1, nspin
          Do lm = 1, lmpotd
            Do ir = irmind, irmd
              vinsa(ir, lm) = 0.E0_dp
              vinsb(ir, lm) = 0.E0_dp
            End Do
          End Do
          ip = nspin*(ih-1) + is
          rmt1 = rmt(ih)
          rmtnw1 = rmtnew(ih)
          z1 = z(ih)
          rmax = rws(ih)
          If (kshape==0) Then
            nr = irws(ih)
            irns1 = 0
          Else
            nr = irc(ih)
            irns1 = irns(ih)
          End If

          irmin = nr - irns1
          a1 = a(ih)
          b1 = b(ih)
          ncore1 = ncore(ip)

          Do j = 1, nr
            ra(j) = r(j, ih)
            dradi(j) = drdi(j, ih)
            vm2za(j) = vm2z(j, ip)
          End Do
          Do lm1 = 1, lmpot
            Do j = irmind, irmd
              vinsa(j, lm1) = vins(j, lm1, ip)
            End Do
          End Do

          If (ncore1>=1) Then

            Do j = 1, ncore1
              lcore1(j) = lcore(j, ip)
              ecore1(j) = ecore(j, ip)
            End Do
          End If

!  Generate uniform mesh RUNI

          nr_u = nr
          irns_u = irns1
          irmin_u = nr_u
          If (ins>0) irmin_u = nr_u - irns_u

          If (ins==0) Then
            Do i = 1, nr_u
              rr_u(i) = ra(i)
              drdi_u(i) = dradi(i)
            End Do
            imt1 = 0
          Else
            imt1 = anint(log(rmtnw1/b1+1.0E0_dp)/a1) + 1
            Do i = 1, imt1
              rr_u(i) = ra(i)
              drdi_u(i) = dradi(i)
            End Do
            rinter = rmax - rmtnw1
            dr = rinter/real(nr-imt1, kind=dp)
            Do i = 1, nr - imt1
              drdi_u(imt1+i) = dr
              rr_u(imt1+i) = rr_u(imt1) + dr*real(i, kind=dp)
            End Do
            Call doubleraus1(nr, irmin, lmpot, ra, dradi, vm2za, vinsa, irmd, &
              irmind, lmpotd)

!     After this sub the arrays are rearanged and nr is not
!     the same anymore in the case of FP. If ins.eq.0 there is
!     no nead for doubleraus1. IRMIN should remain the same

          End If
! ----------------------------------------------------------------
! Now the new mesh is generated

! test
!  write(6,*) nr_u,imt1,irns_u
!  do i=1,nr_u
!    write(6,*) i,ra(i),rr_u(i)
!  end do
! test

          maxa = 1.E35_dp
          Call spline(irmd, ra, vm2za, nr, maxa, maxa, vm2zb)
          If (ins>0) Then
            Do lm1 = 1, lmpot
              irnstot = nr - irmin ! nr has changed irmin is the same
!write(6,*) ' Testing ',nr,irmin,irnstot,irmind
              Call spline(irmd-irmind, ra(irmind), vinsa(irmind,lm1), irnstot, &
                maxa, maxa, vinsb(irmind,lm1))
            End Do ! LM1
          End If

! OK with spline

          Do ir = 1, nr_u
            r0 = rr_u(ir)
            Call splint(ra, vm2za, vm2zb, nr, r0, parsum, parsumderiv)
            vm2z_u(ir) = parsum
          End Do
          If (ins>0) Then
!IRNSTOT = NR_U - IRMIN_U
            Do lm1 = 1, lmpot
              Do ir = irmin_u, nr_u
                r0 = rr_u(ir)
                Call splint(ra(irmind), vinsa(irmind,lm1), vinsb(irmind,lm1), &
                  irnstot, r0, parsum, parsumderiv)
                vins_u(ir, lm1) = parsum
              End Do
            End Do
          End If
!write(6,*) ' All interpolation ok now write'
!     --------------------------------------------------------------
          Write (ifile, Fmt=100)
          nz1 = z1
          If (nspin==1) Then
            Write (ifile, Fmt=110) elemname(nz1), z1
          Else If (is==1) Then
            Write (ifile, Fmt=130) elemname(nz1), z1
          Else If (is==2) Then
            Write (ifile, Fmt=120) elemname(nz1), z1
          End If
          Write (ifile, Fmt=140)
!          write (ifile,*) ALAT,RMAX,RMTNW1,RMT1
          Write (ifile, Fmt=150) alat, rmax, rmtnw1, rmt1
          Write (ifile, Fmt=160) nr_u, imt1, irns1
          Write (ifile, Fmt=170) a1, b1
          Write (ifile, Fmt=180) efermi, vbc(is)
          Write (ifile, Fmt=190) ncore1, lmpot
          If (ncore1>=1) Write (ifile, Fmt=240)(lcore1(icore), ecore1(icore), &
            icore=1, ncore1)

          If (ins==0 .Or. (ih<natps .And. ins<=2)) Then

!---  >       store only the spherically averaged potential
!     (in mt or as - case)
!     this is done always for the host

            Write (ifile, Fmt=260)(vm2z_u(ir), ir=1, nr_u)
          Else

!---  >     store the full potential , but the non spherical contribution
!     only from irns1 up to irws1 ;
!     remember that the lm = 1 contribution is multiplied
!     by a factor 1/sqrt(4 pi)

            Write (ifile, Fmt=270) nr_u, irns1, lmpot, isave
            Write (ifile, Fmt=280)(vm2z_u(ir), ir=1, nr_u)
            If (lpot>0) Then
              lmnr = 1
              Do lm = 2, lmpot
                sum = 0.0E0_dp
                Do ir = irmin, nr_u
                  rv = vins_u(ir, lm)*rr_u(ir)
                  sum = sum + rv*rv*dradi(ir)
                End Do

                If (sqrt(sum)>qbound) Then
                  lmnr = lmnr + 1
                  Write (ifile, Fmt=270) lm
                  Write (ifile, Fmt=280)(vins_u(ir,lm), ir=irmin, nr_u)
                End If

              End Do

!---  >         write a one to mark the end

              If (lmnr<lmpot) Write (ifile, Fmt=270) isave
            End If

          End If

        End Do
      End Do


100   Format (' GENERAL POTENTIAL MESH             exc:')
110   Format ('#  ', A4, 'POTENTIAL             Z = ', F8.3)
120   Format ('#  ', A4, 'POTENTIAL SPIN UP     Z=  ', F8.3)
130   Format ('#  ', A4, 'POTENTIAL SPIN DOWN   Z=  ', F8.3)
140   Format ('#')
150   Format (4F12.8, '   # alat, rmax, rmaxlog, rmt')
160   Format (1P, 3I6, 31X, '  # IRWS, IRMT, IRNS ')
170   Format (2D15.8, 19X, '  # A , B ')
180   Format (3F12.8, 13X, '  # Ef, vbc ')
190   Format (1P, 2I5, 39X, '  # NCORE, LMPOT')
200   Format (7A4, 6X, '  exc:', A24, 3X, A10)
210   Format (3F12.8)
220   Format (F10.5, /, F10.5, 2F15.10)
230   Format (I3, /, 2D15.8, /, 2I2)
240   Format (I5, 1P, D20.11)
250   Format (1P, 2D15.6, 1P, D15.8)
260   Format (1P, 4D20.12)
270   Format (10I5)
280   Format (1P, 4D20.13)
    End Subroutine
! **********************************************************************

!***********************************************************************

    Subroutine spline(nmax, x, y, n, yp1, ypn, y2)
      Use mod_datatypes, Only: dp
! Given arrays x(1:n) and  y(1:n) containing a tabulated function,
! i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn
! for the rst derivative of the interpolating function at points
! 1 and n, respectively, this routine returns an array y2(1:n) of
! length n which contains the second derivatives of the interpolating
! function at the tabulated points xi.
! If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
! signaled to set the corresponding boundary condition for a natural
! spline, with zero second derivative on that boundary.
! Parameter: NMAX is the largest anticipated value of n.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      Implicit None
      Integer :: n, nmax
      Real (Kind=dp) :: yp1, ypn, x(nmax), y(nmax), y2(nmax)
      Integer :: i, k
      Real (Kind=dp) :: p, qn, sig, un, u(nmax)

      If (n>nmax) Stop 'SPLINE: N > NMAX.'
      If (yp1>0.99E30_dp) Then
! The lower boundary condition is set either to be "natural"
        y2(1) = 0.E0_dp
        u(1) = 0.E0_dp
      Else
! or else to have a specified first derivative.
        y2(1) = -0.5E0_dp
        u(1) = (3.E0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      End If

      Do i = 2, n - 1
! This is the decomposition loop of the tridiagonal algorithm. y2 and u
! are used for temporary storage of the decomposed factors.
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1) + 2.E0_dp
        y2(i) = (sig-1.E0_dp)/p
        u(i) = (6.E0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)- &
          x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      End Do

      If (ypn>.99E30_dp) Then
! The upper boundary condition is set either to be "natural"
        qn = 0.E0_dp
        un = 0.E0_dp
      Else
! or else to have a specified 1rst derivative.
        qn = 0.5E0_dp
        un = (3.E0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      End If
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.E0_dp)
      Do k = n - 1, 1, -1
! This is the backsubstitution loop of the tridiagonal algorithm.
        y2(k) = y2(k)*y2(k+1) + u(k)
      End Do

      Return
    End Subroutine
! **********************************************************************

! **********************************************************************
    Subroutine splint(xa, ya, y2a, n, x, y, yderiv)
      Use mod_datatypes, Only: dp
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the xai's in order), and given the array y2a(1:n), which
! is the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y and the derivative yderiv.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      Implicit None
      Integer :: n
      Real (Kind=dp) :: x, y, yderiv, xa(*), ya(*), y2a(*)
      Integer :: k, khi, klo
      Real (Kind=dp) :: a, b, h
! We will find the right place in the table by means of bisection.
! This is optimal if sequential calls to this routine are at random
! values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of
! klo and khi and test if they remain appropriate on the
! next call.
      klo = 1
      khi = n
100   If (khi-klo>1) Then
        k = (khi+klo)/2
        If (xa(k)>x) Then
          khi = k
        Else
          klo = k
        End If
        Go To 100
      End If
! klo and khi now bracket the input value of x.
      h = xa(khi) - xa(klo)
! The xa's must be distinct.
      If (h==0.E0_dp) Stop 'Bad XA input in SPLINT'
! Cubic spline polynomial is now evaluated.
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2) &
        /6.E0_dp
      yderiv = (ya(khi)-ya(klo))/h - ((3.E0_dp*a*a-1.E0_dp)*y2a(klo)-(3.E0_dp* &
        b*b-1.E0_dp)*y2a(khi))*h/6.E0_dp

      Return
    End Subroutine
! **********************************************************************

! **********************************************************************
!************************************************************************

    Subroutine doubleraus1(irmax, irmin, lmpot, rr, drdi, vpot, vins, irmd, &
      irmind, lmpotd)
      Use mod_datatypes, Only: dp
! Gets rid of the double-points in the radial mesh, i.e. the points
! where RR(I) = RR(I-1). Returns the "new" mesh in the same array,
! and rearranges accordingly the WAVEF defined at the same mesh.
! IRMAX is also altered to the new value.
      Implicit None
      Integer :: irmd, lmpotd, irmind
      Integer :: ncountmax
      Parameter (ncountmax=500)
! Input and output:
      Integer :: irmax
      Real (Kind=dp) :: rr(irmd), drdi(irmd), vpot(irmd), &
        vins(irmind:irmd, lmpotd)
! Inside:
      Integer :: ir, icount, nc, ncount, lmpot, irmin
      Integer :: lm1
      Integer :: idouble(ncountmax)

! Find double-points:
      ncount = 0
      Do ir = 2, irmax
        If ((rr(ir)-rr(ir-1))<1.E-20_dp) Then
          ncount = ncount + 1
          idouble(ncount) = ir
        End If
      End Do
      If (ncount+1>ncountmax) Stop 'DOUBLERAUS2: Too many double-points.'
      idouble(ncount+1) = irmax + 1 ! To be used below.

! Rearrange the arrays.
      Do icount = 1, ncount
        Do ir = idouble(icount) - icount + 1, idouble(icount+1) - icount
          If ((ir+icount)<=irmax) Then
            rr(ir) = rr(ir+icount)
            drdi(ir) = drdi(ir+icount)
            vpot(ir) = vpot(ir+icount)
          End If
        End Do
      End Do
      irmax = irmax - ncount
      ncount = 0
      Do nc = 1, ncountmax
        idouble(nc) = 0
      End Do

      Do ir = irmin, irmax
        If ((rr(ir)-rr(ir-1))<1.E-20_dp) Then
          ncount = ncount + 1
          idouble(ncount) = ir
        End If
      End Do
! Rearrange the arrays.
      Do icount = 1, ncount
        Do ir = idouble(icount) - icount + 1, idouble(icount+1) - icount
          Do lm1 = 1, lmpot
            vins(ir, lm1) = vins(ir+icount, lm1)
          End Do
        End Do
      End Do
      Return
    End Subroutine
