!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculates the charge density inside r(irmin) in case of a non spherical input potential.
!> Author: B. Drittler
!> Calculates the charge density inside r(irmin) in case of a non spherical input potential.
!> Fills the array `cden` for the complex density of states the non spher. 
!> wavefunctions are approximated in that region in the following way:
!> 
!> * The regular one (ir < irmin = irws-irns) :
!> \begin{equation}
!> pns(ir,lm1,lm2) = pz(ir,l1)  ar(lm1,lm2)
!> \end{equation}
!> where \(pz\) is the regular wavefct of the spherically symmetric part of the 
!> potential and ar the alpha matrix (see subroutine `regns`)
!> * The irregular one (ir < irmin) :
!> \begin{equation}
!> qns(ir,lm1,lm2) = pz(ir,l1) cr(lm1,lm2)+ qz(ir,l1) dr(lm1,lm2)
!> \end{equation}
!> where \(pz\) is the regular and \(qz\) is the irregular wavefct of the 
!> spherically symmetric part of the potential and \(cr\) , \(dr\) the matrices calculated
!> at the point `irmin` (see subroutine `irwns()`).
!> The structured part of the greens-function (`gmat`) is symmetric in its lm-indices,
!> therefore only one half of the matrix is calculated in the subroutine for the 
!> back-symmetrisation .
!> Rhe gaunt coeffients are symmetric too (since the are calculated using the real spherical harmonics).
!> Rhat is why the `lm2-` and the `lm02-` loops are only only going up to `lm1` or `lm01`
!> and the summands are multiplied by a factor of 2 in the case of `lm1 .ne. lm2` or `lm01 .ne. lm02`. 
!> (see notes by B. Drittler)
!------------------------------------------------------------------------------------
!> @note
!> 
!> * Remember that the matrices `ar`,`cr`,`dr` are rescaled (see subroutines `irwns()` and `regns()`)
!> * Arrays `rho2ns` and `cden` are initialized in subroutine `rhoout`.
!> @endnote
!> @warning The gaunt coeffients which are used here are ordered in a special way
!> (see subroutine `gaunt()`)
!> @endwarning
!------------------------------------------------------------------------------------
module mod_rhoin

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates the charge density inside r(irmin) in case of a non spherical input potential.
  !> Author: B. Drittler
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculates the charge density inside r(irmin) in case of a non spherical input potential.
  !> Fills the array `cden` for the complex density of states the non spher. 
  !> wavefunctions are approximated in that region in the following way:
  !> 
  !> * The regular one (ir < irmin = irws-irns) :
  !> \begin{equation}
  !> pns(ir,lm1,lm2) = pz(ir,l1)  ar(lm1,lm2)
  !> \end{equation}
  !> where \(pz\) is the regular wavefct of the spherically symmetric part of the 
  !> potential and ar the alpha matrix (see subroutine `regns`)
  !> * The irregular one (ir < irmin) :
  !> \begin{equation}
  !> qns(ir,lm1,lm2) = pz(ir,l1) cr(lm1,lm2)+ qz(ir,l1) dr(lm1,lm2)
  !> \end{equation}
  !> where \(pz\) is the regular and \(qz\) is the irregular wavefct of the 
  !> spherically symmetric part of the potential and \(cr\) , \(dr\) the matrices calculated
  !> at the point `irmin` (see subroutine `irwns()`).
  !> The structured part of the greens-function (`gmat`) is symmetric in its lm-indices,
  !> therefore only one half of the matrix is calculated in the subroutine for the 
  !> back-symmetrisation .
  !> Rhe gaunt coeffients are symmetric too (since the are calculated using the real spherical harmonics).
  !> Rhat is why the `lm2-` and the `lm02-` loops are only only going up to `lm1` or `lm01`
  !> and the summands are multiplied by a factor of 2 in the case of `lm1 .ne. lm2` or `lm01 .ne. lm02`. 
  !> (see notes by B. Drittler)
  !-------------------------------------------------------------------------------
  !> @note 
  !> * Remember that the matrices `ar`,`cr`,`dr` are rescaled (see subroutines `irwns()` and `regns()`)
  !> * Arrays `rho2ns` and `cden` are initialized in subroutine `rhoout`.
  !> @endnote
  !> @warning The gaunt coeffients which are used here are ordered in a special way
  !> (see subroutine `gaunt()`)
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine rhoin(ar,cden,cr,df,gmat,ek,rho2ns,irc1,nsra,efac,pz,fz,qz,sz,cleb,    &
    icleb,jend,iend,ekl,cdenlm)

    use :: mod_datatypes, only: dp
    use :: global_variables, only: lmmaxd, lmaxd, irmd, ncleb, lmpotd
    use :: mod_constants, only: pi, czero
    implicit none
    ! interface
    complex (kind=dp) :: df, ek
    integer :: iend, irc1, nsra
    complex (kind=dp) :: ar(lmmaxd, lmmaxd), cden(irmd, 0:lmaxd), cr(lmmaxd, lmmaxd), efac(lmmaxd), ekl(0:lmaxd), fz(irmd, 0:lmaxd), gmat(lmmaxd, lmmaxd), pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), &
      sz(irmd, 0:lmaxd), cdenlm(irmd, lmmaxd)
    real (kind=dp) :: cleb(ncleb), rho2ns(irmd, lmpotd)
    integer :: icleb(ncleb, 4), jend(lmpotd, 0:lmaxd, 0:lmaxd)
    ! local
    complex (kind=dp) :: efac1, efac2, ffz, gmatl, ppz, v1, v2
    real (kind=dp) :: c0ll
    integer :: i, ir, j, j0, j1, l, l1, l2, lm1, lm2, lm3, lm3max, ln2, ln3, m
    complex (kind=dp) :: vr(lmmaxd, lmmaxd), wf(irmd, 0:lmaxd, 0:lmaxd), wr(lmmaxd, lmmaxd)
#ifdef __GFORTRAN__
    ! for the gfortran compiler zdotu leads to a segfault, then using
    ! dot_product gives the correct result without segfault, this needs these
    ! auxiliary arrays
    complex (kind=dp) :: ar2(lmmaxd**2), vr2(lmmaxd**2)
#endif
    complex (kind=dp), external :: zdotu

    ! C0LL = 1/sqrt(4*pi)
    c0ll = 1.0e0_dp/sqrt(4.0e0_dp*pi)

    lm3max = icleb(iend, 3)

    ! ---> set up array wr(lm1,lm2)
    !      use first vr
    wr = czero
    vr = czero
    do lm2 = 1, lmmaxd
      ln2 = lm2
      v2 = efac(lm2)*efac(lm2)*gmat(ln2, ln2)
      do lm1 = 1, lmmaxd
        vr(lm1, lm2) = ek*cr(lm1, lm2) + v2*ar(lm1, lm2)
      end do
    end do

    ! ---> using symmetry of structural green function
    do lm2 = 2, lmmaxd
      ln2 = lm2
      efac2 = 2.0e0_dp*efac(lm2)
      do lm3 = 1, lm2 - 1
        ln3 = lm3
        v1 = efac2*gmat(ln3, ln2)*efac(lm3)
        do lm1 = 1, lmmaxd
          vr(lm1, lm2) = vr(lm1, lm2) + v1*ar(lm1, lm3)
        end do
      end do
    end do

#ifdef __GFORTRAN__
    ar2 = reshape(ar, [lmmaxd**2])
    vr2 = reshape(vr, [lmmaxd**2])
#endif
    do lm1 = 1, lmmaxd
      efac1 = efac(lm1)
#ifdef __GFORTRAN__
      wr(lm1, lm1) = dot_product(conjg(ar2(lm1::lmmaxd)), vr2(lm1::lmmaxd))/(efac1*efac1) ! this works only with conjugation due to definition of dot_product
#else
      wr(lm1, lm1) = zdotu(lmmaxd, ar(lm1,1), lmmaxd, vr(lm1,1), lmmaxd)/(efac1*efac1) ! this works with the intel compiler
#endif
      do lm2 = 1, lm1 - 1
        ! ---> using symmetry of gaunt coeffients
        efac2 = efac(lm2)
#ifdef __GFORTRAN__
        wr(lm1, lm2) = ( dot_product(conjg(ar2(lm1::lmmaxd)),vr2(lm2::lmmaxd)) + dot_product(conjg(ar2(lm2::lmmaxd)),vr2(lm1::lmmaxd)) )/(efac1*efac2) ! this works with gfortran
#else
        wr(lm1, lm2) = ( zdotu(lmmaxd,ar(lm1,1),lmmaxd,vr(lm2,1),lmmaxd) + zdotu(lmmaxd,ar(lm2,1),lmmaxd,vr(lm1,1),lmmaxd) )/(efac1*efac2) ! this works with the intel compiler
#endif
      end do
    end do

    ! ---> set up array wf(l1,l2) = pz(l1)*pz(l2)
    if (nsra==2) then
      do l1 = 0, lmaxd
        do l2 = 0, l1
          do ir = 2, irc1
            wf(ir, l1, l2) = pz(ir, l1)*pz(ir, l2) + fz(ir, l1)*fz(ir, l2)
          end do
        end do
      end do

    else
      do l1 = 0, lmaxd
        do l2 = 0, l1
          do ir = 2, irc1
            wf(ir, l1, l2) = pz(ir, l1)*pz(ir, l2)
          end do
        end do
      end do
    end if

    ! ---> first calculate only the spherically symmetric contribution
    ! remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
    do l = 0, lmaxd
      gmatl = czero
      do m = -l, l
        lm1 = l*(l+1) + m + 1
        gmatl = gmatl + wr(lm1, lm1)
      end do
      if (nsra==2) then
        do i = 2, irc1
          ppz = pz(i, l)
          ffz = fz(i, l)
          cden(i, l) = ppz*(gmatl*ppz+ekl(l)*qz(i,l)) + ffz*(gmatl*ffz+ekl(l)*sz(i,l))
          rho2ns(i, 1) = rho2ns(i, 1) + c0ll*aimag(df*cden(i,l)) ! Implicit integration over energies
          do m = -l, l
            lm1 = l*(l+1) + m + 1
            cdenlm(i, lm1) = ppz*(wr(lm1,lm1)*ppz+ek*qz(i,l)) + ffz*(wr(lm1,lm1)*ffz+ek*sz(i,l))
          end do
        end do
      else
        do i = 2, irc1
          ppz = pz(i, l)
          cden(i, l) = ppz*(gmatl*ppz+ekl(l)*qz(i,l))
          rho2ns(i, 1) = rho2ns(i, 1) + c0ll*aimag(df*cden(i,l)) ! Implicit integration over energies
          do m = -l, l
            lm1 = l*(l+1) + m + 1
            cdenlm(i, lm1) = ppz*(wr(lm1,lm1)*ppz+ek*qz(i,l))
          end do
        end do
      end if
    end do

    ! ---> calculate the non spherically symmetric contribution
    ! to speed up the pointer jend generated in gaunt is used
    ! remember that the wavefunctions are l and not lm dependent

    j0 = 1

    do lm3 = 2, lm3max

      do l1 = 0, lmaxd
        do l2 = 0, l1
          j1 = jend(lm3, l1, l2)

          if (j1/=0) then

            gmatl = czero

            ! ---> sum over m1,m2 for fixed lm3,l1,l2
            do j = j0, j1
              lm1 = icleb(j, 1)
              lm2 = icleb(j, 2)
              gmatl = gmatl + cleb(j)*wr(lm1, lm2)
            end do

            j0 = j1 + 1

            gmatl = df*gmatl
            do i = 2, irc1
              rho2ns(i, lm3) = rho2ns(i, lm3) + aimag(gmatl*wf(i,l1,l2))
            end do

          end if ! (j1/=0)

        end do ! l2
      end do ! l1
    end do ! lm3

  end subroutine rhoin

end module mod_rhoin
