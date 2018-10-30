!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Generate an angular mesh and spherical harmonics for the treatement of the GGA xc-potential
!> Author: R. Zeller, Phivos Mavropoulos
!> Generate an angular mesh and spherical harmonics for the treatement of the GGA 
!> xc-potential. For an angular integration the weights are also generated.
!------------------------------------------------------------------------------------
module mod_sphere_gga
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Generate an angular mesh and spherical harmonics at those mesh points.
  !> Author: R. Zeller
  !> Category: xc-potential, special-functions, KKRhost
  !> Deprecated: False 
  !> Generate an angular mesh and spherical harmonics at those mesh points. For an 
  !> angular integration the weights are also generated. This is needed for the 
  !> calculation of the GGA exchange correlation potentials.
  !-------------------------------------------------------------------------------
  !> @note Phivos Mavropoulos, July 2007: New call to subroutine `ylmderiv` for
  !> accurate derivatives of spherical harmonics.
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine sphere_gga(lmax,yr,wtyr,rij,ijd,lmmaxd,thet,ylm,dylmt1,dylmt2,dylmf1,  &
    dylmf2,dylmtf)

    use :: mod_datatypes, only: dp
    use :: mod_lebedev, only: lebedev
    use :: mod_ymy, only: ymy
    use :: mod_rinit, only: rinit
    use :: mod_constants, only: pi
    implicit none

    ! .. Scalar Arguments ..
    integer :: ijd, lmax, lmmaxd
    ! ..
    ! .. Local Scalars ..
    real (kind=dp) :: r, r1, r2, r3
    integer :: ij, lm1
    ! ..
    ! .. Array Arguments ..
    real (kind=dp) :: dylmf1(ijd, lmmaxd), dylmf2(ijd, lmmaxd), dylmt1(ijd, lmmaxd)
    real (kind=dp) :: dylmt2(ijd, lmmaxd), dylmtf(ijd, lmmaxd), rij(ijd, 3), thet(ijd)
    real (kind=dp) :: wtyr(ijd, *), ylm(ijd, lmmaxd), yr(ijd, *), dydth(lmmaxd)
    real (kind=dp) :: dydfi(lmmaxd), d2ydth2(lmmaxd), d2ydfi2(lmmaxd), d2ydthdfi(lmmaxd)
    ! ..
    ! .. Local Arrays ..
    real (kind=dp) :: wght, y(1000)
    ! ..
    write (1337, *) 'SPHERE for GGA: read LEBEDEV mesh'
    if (ijd>1000) stop 'SPHERE'


    do ij = 1, ijd
      call lebedev(ij, r1, r2, r3, wght)
      rij(ij, 1) = r1
      rij(ij, 2) = r2
      rij(ij, 3) = r3

      ! For the needs of GGA PW91 as implemented here, ylm and derivatives
      ! come with a different sign convention compared to the usual in the
      ! program: sin(fi)**m --> -sin(fi)**m. Thus some signs change
      ! also in array ylm compared to array yr (below).
      call derivylm(r1, r2, r3, lmax, r, y, dydth, dydfi, d2ydth2, d2ydfi2, d2ydthdfi)

      thet(ij) = acos(r3/r)

      do lm1 = 1, (lmax+1)**2
        ylm(ij, lm1) = y(lm1)
        dylmt1(ij, lm1) = dydth(lm1)
        dylmf1(ij, lm1) = dydfi(lm1)
        dylmt2(ij, lm1) = d2ydth2(lm1)
        dylmf2(ij, lm1) = d2ydfi2(lm1)
        dylmtf(ij, lm1) = d2ydthdfi(lm1)
      end do

      ! Call ymy to obtain sher. harmonics with usual convention
      ! ---> multiply the spherical harmonics with the weights
      call ymy(r1, r2, r3, r, y, lmax)
      do lm1 = 1, (lmax+1)**2
        yr(ij, lm1) = y(lm1)
        wtyr(ij, lm1) = yr(ij, lm1)*wght*pi*4.e0_dp
      end do
    end do

  end subroutine sphere_gga

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the 1st and 2nd derivatives of real spherical harmonics with respect to \(\theta\), \(\phi\)
  !> Author: Ph.Mavropoulos
  !> Category: xc-potential, special-functions, KKRhost
  !> Deprecated: False 
  !> Calculate the 1st and 2nd derivatives of real spherical harmonics
  !> with respect to \(\theta\), \(\phi\).
  !> Use recursion relations for the assoc. Legendre functions \(P[l,m]\) to generate
  !> the derivatives. These are (taken from Abramowitz and Stegun, Handbook of 
  !> Mathematical Functions, chapt. 8.):
  !> \begin{equation}
  !> P[l,m+1] = (x^2-1)^{-\frac{1}{2}} ( (l-m)xP[l,m] - (l+m)P[l-1,m] )
  !> \end{equation}
  !> \begin{equation}
  !> (x^2-1)\frac{d}{dx}P[l,m] = (l+m)(l-m+1)(x^2-1)^\frac{1}{2} P[l,m-1] - mxP[l,m]
  !> \end{equation}
  !> \begin{equation}
  !> (x^2-1)\frac{d}{dx}P[l,m] = lxP[l,m] - (l+m)P[l-1,m] 
  !> \end{equation}
  !> where \(x=\cos{\theta}\), \((x^2-1)^\frac{1}{2} = -\sin{\theta}\), \(\frac{d}{d\theta} = -\sin{\theta} \frac{d}{dx}\)
  !> Adding these equations:
  !> \begin{equation}
  !> \frac{d}{d\theta}P[l,m](\cos{\theta}) = \frac{1}{2} ( -(l+m)(l-m+1)P[l,m-1] + P[l,m+1] )
  !> \end{equation}
  !> It is implied that \(P[l,m]=0\) if \(m>l\) or \(m<-l\). Also, the term \((x^2-1)^\frac{1}{2}\)
  !> is ambiguous for real \(x\), \(0<x<1\); here it is interpreted as
  !> \begin{equation}
  !> (x^2-1)^\frac{1}{2}=-\sin{\theta} 
  !> \end{equation}
  !> but
  !> \begin{equation}
  !> (x^2-1)^{-\frac{1}{2}}=\frac{1}{\sin{\theta}}
  !> \end{equation}
  !> otherwise the result from Eq.4 (which is cross-checked and correct) does not follow.
  !> For the 2nd derivative apply Eq.4 twice. Result:
  !> \begin{equation}
  !> \begin{split}
  !> \frac{d^2}{d\theta^2}P[l,m](\cos{\theta}) = \frac{1}{4}((l+m)(l-m+1)(l+m-1)(l-m+2) P[l,m-2]&\\-((l-m)(l+m+1)+(l+m)(l-m+1) )P[l,m]+ P[l,m+2])
  !> \end{split}
  !> \end{equation}
  !> The \(\phi\)-derivatives act on \cos{\phi},\sin{\phi} and are trivial.
  !> For the associated Legendre functions use the recursion formulas:
  !> \begin{equation}
  !> (l-m+1)P[l+1,m] = (2l+1)\cos{\theta}P[l,m] - (l+m)P[l-1,m] 
  !> \end{equation}
  !> \begin{equation}
  !> P[l+1,m] = P[l-1,m] - (2*l+1)\sin{\theta}P[l,m-1] 
  !> \end{equation}
  !> ( with \(x=\cos{\theta} \).
  !> Recursion algorithm for the calculation of \(P[l,m]\) and calculation of \(Y_l^m\)
  !> taken over from subr. `ymy` of KKR program (implemented there by M. Weinert, B. Drittler).
  !> For \(m<0\), use 
  !> \begin{equation}
  !> P[l,-m] = P[l,m] \frac{(l-m)!}{(l+m)!}
  !> \end{equation}
  !> Taking into account the lm-prefactors of the spherical harmonics, we construct 
  !> and use the functions
  !> \begin{equation}
  !> Q[l,m] = \sqrt{\frac{2l+1}{4\pi}} \sqrt{\frac{(l-m)!}{(l+m)!}} P[l,m]
  !> \end{equation}
  !> whence Eq.4 and Eq.7 become
  !> \begin{equation}
  !> \frac{d}{d\theta}Q[l,m]= \frac{1}{2}(-\sqrt{(l+m)(l-m+1)}Q[l,m-1]+ \sqrt{(l+m+1)(l-m)}Q[l,m+1] )
  !> \end{equation}
  !> \begin{equation}
  !> \begin{split}
  !> \frac{d^2}{d\theta^2}Q[l,m] = \frac{1}{4}(\sqrt{(l+m)(l+m-1)(l-m+1)(l-m+2)} Q[l,m-2]&\\+((l-m)(l+m+1)+(l+m)(l-m+1)) Q[l,m]&\\+ \sqrt{(l-m)(l-m-1)(l+m+1)(l+m+2)} Q[l,m+2])
  !> \end{split}
  !> \end{equation}
  !> Note on sign convension:
  !> For the needs of GGA PW91 as implemented here, ylm and derivatives
  !> come with a different sign convention compared to the usual in the
  !> program: \(\sin{\phi}^m \rightarrow (-1)^m \sin{\phi}^m\). Thus some signs change.  
  !-------------------------------------------------------------------------------
  subroutine derivylm(v1, v2, v3, lmax, rabs, ylm, dydth, dydfi, d2ydth2, d2ydfi2, d2ydthdfi)

    use :: mod_datatypes, only: dp
    use :: mod_rinit, only: rinit
    use :: mod_constants, only: pi
    implicit none
    ! Parameters:
    integer :: lmaxd, l4maxd
    parameter (lmaxd=4, l4maxd=4*lmaxd)
    ! Input:
    integer :: lmax                ! up to which l to calculate
    real (kind=dp) :: v1, v2, v3   ! vector where Ylm etc are calculated (not
    ! necessarily normalized)
    ! Output:
    ! Y[l,m], dY/dth, dY/dfi, d(dY/dth)/dth, d(dY/dfi)/dfi, d(dY/dth)/dfi
    real (kind=dp) :: ylm(*), dydth(*), dydfi(*), d2ydth2(*), d2ydfi2(*), d2ydthdfi(*)
    real (kind=dp) :: rabs         ! Norm of input vector (V1,V2,V3)
    ! Inside:
    real (kind=dp) :: cth, sth, cfi, sfi ! cos and sin of th and fi
    real (kind=dp) :: fpi, rtwo ! pi (what else?), 4*pi, sqrt(2)
    real (kind=dp) :: fac          ! factor in construction of polynomials.
    real (kind=dp) :: plm(0:l4maxd, 0:l4maxd) ! Legendre polynomials
    real (kind=dp) :: qlm((l4maxd+1)**2) ! Ylm/cos(m*fi) (m>0) and Ylm/sin(m*fi)
    ! (m<0)
    real (kind=dp) :: cmfi(0:l4maxd), smfi(0:l4maxd) ! cos(m*fi) and sin(m*fi)
    real (kind=dp) :: xy, xyz, sgm, sgmm, fi
    real (kind=dp) :: aux
    real (kind=dp) :: tiny
    parameter (tiny=1.e-20_dp)     ! if th < tiny set th=0
    real (kind=dp) :: tt, aa, cd   ! factors in calcul. of Ylm
    integer :: ll, mm, ii          ! l and m indexes
    integer :: lmmax               ! (lmax+1)**2, total number of spher.
    ! harmonics.
    integer :: imm, ipm, lpm, lmm, lpmp1, lmmp1 ! i-m,i+m,l+m,l-m,l+m+1,l-m-1

    fpi = 4.e0_dp*pi
    rtwo = sqrt(2.e0_dp)
    lmmax = (lmax+1)**2

    if (lmax>l4maxd) stop 'derivylm: lmax out of range.'

    ! --->    calculate sin and cos of theta and phi
    xy = v1**2 + v2**2
    xyz = xy + v3**2

    rabs = sqrt(xyz)
    if (xyz<=0.0e0_dp) stop 'derivylm: v=0.'

    if (xy>tiny*xyz) then
      xy = sqrt(xy)
      xyz = sqrt(xyz)
      cth = v3/xyz
      sth = xy/xyz
      cfi = v1/xy
      sfi = v2/xy

    else

      sth = 0.0e0_dp
      cth = 1.0e0_dp
      if (v3<0) cth = -1.0e0_dp
      cfi = 1.0e0_dp
      sfi = 0.0e0_dp
    end if

    ! First calculate Legendre functions. Use recursion formulas (8.5.3,8.5.5).
    ! Following taken from KKR program (routine ymy, by M.Weinert).
    fac = 1.0e0_dp
    do mm = 0, lmax - 1
      fac = -real(2*mm-1, kind=dp)*fac
      plm(mm, mm) = fac
      plm(mm+1, mm) = real(2*mm+1, kind=dp)*cth*fac

      ! --->    recurse upward in l
      do ll = mm + 2, lmax
        plm(ll, mm) = (real(2*ll-1,kind=dp)*cth*plm(ll-1,mm)-real(ll+mm-1,kind=dp)*plm(ll-2,mm))/real(ll-mm, kind=dp)
      end do
      fac = fac*sth
    end do
    plm(lmax, lmax) = -(2*lmax-1)*fac

    ! Next calculate Ylm and derivatives.
    ! --->    determine powers of sin and cos of phi
    smfi(0) = 0.0e0_dp
    smfi(1) = sfi
    cmfi(0) = 1.0e0_dp
    cmfi(1) = cfi
    do mm = 2, lmax
      smfi(mm) = 2.e0_dp*cfi*smfi(mm-1) - smfi(mm-2)
      cmfi(mm) = 2.e0_dp*cfi*cmfi(mm-1) - cmfi(mm-2)
    end do

    ! For the needs of GGA PW91 as implemented here, ylm and derivatives
    ! come with a different sign convention compared to the usual in the
    ! program: sin(fi)**m --> (-1)**m * sin(fi)**m. Thus some signs change.
    ! This is taken care of here:
    fi = atan2(v2, v1)
    ! THE CHANGE OF SIGN BELOW IS WRONG AND THEREFORE NOT DONE ANYMORE
    ! It was introduced to keep results consistent with older versions which
    ! already yielded wrong results
    ! if (fi.lt.0.d0) then
    ! do mm = 1,lmax
    ! smfi(mm) = -smfi(mm)
    ! enddo
    ! endif
    ! --->    multiply in the normalization factors;
    ! calculate Ylm and derivatives with respect to fi.
    ii = 0
    do ll = 0, lmax
      ii = ii + ll + 1
      aa = sqrt(real(2*ll+1,kind=dp)/fpi)
      cd = 1.e0_dp
      ylm(ii) = aa*plm(ll, 0)
      dydfi(ii) = 0.e0_dp
      d2ydfi2(ii) = 0.e0_dp

      qlm(ii) = rtwo*aa*plm(ll, 0)
      sgm = -rtwo                  ! updated to (-1)**m * rtwo
      sgmm = -1                    ! updated to (-1)**m
      do mm = 1, ll
        ipm = ii + mm
        imm = ii - mm
        tt = real((ll+1-mm)*(ll+mm), kind=dp)
        cd = cd/tt
        tt = aa*sqrt(cd)

        qlm(ipm) = sgm*tt*plm(ll, mm)
        qlm(imm) = sgmm*qlm(ipm)

        ylm(ipm) = qlm(ipm)*cmfi(mm)
        ylm(imm) = sgmm*qlm(imm)*smfi(mm)

        dydfi(ipm) = -real(mm, kind=dp)*ylm(imm)
        dydfi(imm) = real(mm, kind=dp)*ylm(ipm)
        d2ydfi2(ipm) = -real(mm*mm, kind=dp)*ylm(ipm)
        d2ydfi2(imm) = -real(mm*mm, kind=dp)*ylm(imm)

        sgm = -sgm
        sgmm = -sgmm

      end do
      ii = ii + ll
    end do

    ! Derivatives with respect to th
    call rinit(lmmax, dydth)
    call rinit(lmmax, d2ydth2)
    call rinit(lmmax, d2ydthdfi)
    ! The l=0 derivatives are zero (established by initialization above).
    ! Start with l=1.
    do ll = 1, lmax
      ii = ll*(ll+1) + 1           ! (position of m=0 harmonic in array)
      aa = real(ll*(ll+1), kind=dp)
      ! Take special care of m=0 harmonic due to 1/sqrt(2)
      dydth(ii) = -sqrt(aa)*qlm(ii+1)/rtwo
      aux = -2.e0_dp*aa*qlm(ii)
      if (ll>1) aux = aux + (qlm(ii-2)+qlm(ii+2))*sqrt(real((ll-1)*ll*(ll+1)*(ll+2),kind=dp))
      d2ydth2(ii) = 0.25e0_dp*aux/rtwo

      do mm = 1, ll
        ipm = ii + mm
        imm = ii - mm

        lpm = ll + mm
        lmm = ll - mm
        lpmp1 = ll + mm + 1
        lmmp1 = ll - mm + 1
        ! Apply Eq. (A1)
        aux = qlm(ipm-1)*sqrt(real(lpm*lmmp1,kind=dp))
        if (mm<ll) aux = aux - qlm(ipm+1)*sqrt(real(lpmp1*lmm,kind=dp))
        aux = 0.5e0_dp*aux

        dydth(ipm) = aux*cmfi(mm)
        dydth(imm) = aux*smfi(mm)

        d2ydthdfi(ipm) = -real(mm, kind=dp)*aux*smfi(mm)
        d2ydthdfi(imm) = real(mm, kind=dp)*aux*cmfi(mm)

        ! Apply Eq. (B1)
        aux = -qlm(ipm)*real(lmm*lpmp1+lpm*lmmp1, kind=dp)
        if (mm<ll-1) aux = aux + qlm(ipm+2)*sqrt(real((lmm-1)*lmm*lpmp1*(lpm+2),kind=dp))
        aux = aux + qlm(ipm-2)*sqrt(real(lmmp1*(lmm+2)*(lpm-1)*lpm,kind=dp))
        aux = 0.25e0_dp*aux

        d2ydth2(ipm) = aux*cmfi(mm)
        d2ydth2(imm) = aux*smfi(mm)

      end do
    end do

  end subroutine derivylm

end module mod_sphere_gga
