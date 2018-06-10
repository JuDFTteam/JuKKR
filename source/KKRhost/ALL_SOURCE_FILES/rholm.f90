subroutine rholm(den, df, gmat, nsra, rho2ns, drdi, ipan, ircut, pz, fz, qz, &
  sz, cleb, icleb, iend, jend, ekl)
!-----------------------------------------------------------------------
!     calculate in the paramagnetic case (nspin=1) :
!         the valence charge density times r**2 from the greensfunction
!     calculate in the spin-polarized case (nspin=2) :
!         the valence charge density times r**2 and the valence spin
!         density times r**2 from the greensfunction ,
!         ( convention spin density :=
!                            density(spin up)-density(spin down) )
!     calculate the valence density of states , in the spin-polarized
!      case spin dependent ; splitted into its l-contributions .

!     in this subroutine an implicit energy-spin integration is  done :
!        this subroutine is called for each energy and spin value
!        and n(r,e) times df (the energy weight) is calculated .

!     recognize that the density of states is always complex also in
!      the case of "real-energy-integation" (ief>0) since in that case
!      the energy integration is done parallel to the real energy axis
!      but not on the real energy axis .
!      in the paramagnetic case only rho2ns(irmd,lmxtsq,natypd,1)
!      is used containing  the charge density times r**2 .
!      in the spin-polarized case rho2ns(...,1) contains the charge
!      density times r**2 and rho2ns(...,2) the spin density times
!      r**2 .

!     the charge density is expanded in spherical harmonics :

!             rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)

!          rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
!                                                           unit sphere)
!     in the case of spin-polarization :
!       the spin density is developed in spherical harmonics :

!            sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)

!         sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
!                                                           unit sphere)
!     n(r,e) is developed in

!        n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }

!     therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
!     has to be used to calculate the lm-contribution of the charge
!     density .
!             (see notes by b.drittler)

!     attention : the gaunt coeffients are stored in an index array
!                 (see subroutine gaunt)
!                 the structure part of the greens-function (gmat) is
!                 symmetric in its lm-indices , therefore only one
!                 half of the matrix is calculated in the subroutine
!                 for the back-symmetrisation . the gaunt coeffients
!                 are symmetric too (since the are calculated for
!                 real spherical harmonics) . that is why the lm2-
!                 loop only goes up to lm1 and the summands are
!                 multiplied by a factor of 2 in the case of lm1
!                 not equal to lm2 .

!                               b.drittler   may 1987
!                                   changed  dec 1988
!-----------------------------------------------------------------------
!     .. Parameters ..
  include 'inc.p'
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
!..
!.. Scalar Arguments ..
!..
!.. Array Arguments ..
  integer :: lmmaxd
  parameter (lmmaxd=(krel+1)*(lmaxd+1)**2)
  integer :: lmaxd1
  parameter (lmaxd1=lmaxd+1)
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  double complex :: czero
  parameter (czero=(0.0d0,0.0d0))
!..
!.. Local Scalars ..
  double complex :: df
  integer :: iend, ipan, nsra
!..
!.. Local Arrays ..
  double complex :: den(0:lmaxd1), ekl(0:lmaxd), fz(irmd, 0:lmaxd), &
    gmat(lmmaxd, lmmaxd), pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), &
    sz(irmd, 0:lmaxd)
  double precision :: cleb(*), drdi(irmd), rho2ns(irmd, lmpotd)
  integer :: icleb(ncleb, 4), ircut(0:ipand), jend(lmpotd, 0:lmaxd, 0:lmaxd)
!..
!.. External Subroutines ..
  double complex :: ffz, gmatl, ppz
  double precision :: c0ll, facsym, pi
  integer :: i, j, j0, j1, l, l1, l2, lm3, lm3max, ln1, ln2, lne, lns
!..
!.. Intrinsic Functions ..
  double complex :: denr(irmd), wr(irmd, 0:lmaxd, 0:lmaxd)
!..
!.. Save statement ..
  external :: csimpk
!..

  intrinsic :: atan, dimag, sqrt


  save

  pi = 4.0d0*atan(1.0d0)
  c0ll = 1.0d0/sqrt(4.0d0*pi)
!---> set up of wr(ir,l1,l2) = pz(ir,l1)*pz(ir,l2)

  lm3max = icleb(iend, 3)




  if (nsra==2) then
    do l1 = 0, lmaxd
      do l2 = 0, l1
        do i = 2, ircut(1)
          wr(i, l1, l2) = pz(i, l1)*pz(i, l2) + fz(i, l1)*fz(i, l2)
        end do
      end do
    end do
!---> first calculate only the spherically symmetric contribution
  else

    do l1 = 0, lmaxd
      do l2 = 0, l1
        do i = 2, ircut(1)
          wr(i, l1, l2) = pz(i, l1)*pz(i, l2)
        end do
      end do
    end do

  end if
!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)


  do l = 0, lmaxd
    gmatl = czero
    lns = l*l + 1
    lne = lns + 2*l
    do ln1 = lns, lne
      gmatl = gmatl + gmat(ln1, ln1)
    end do



    denr(1) = czero
    if (nsra==2) then
      do i = 2, ircut(1)
        ppz = pz(i, l)
        ffz = fz(i, l)
        denr(i) = ppz*(gmatl*ppz+ekl(l)*qz(i,l)) + ffz*(gmatl*ffz+ekl(l)*sz(i, &
          l))
        rho2ns(i, 1) = rho2ns(i, 1) + c0ll*dimag(df*denr(i))
      end do
!---> calculate density of states
    else

      do i = 2, ircut(1)
        ppz = pz(i, l)
        denr(i) = ppz*(gmatl*ppz+ekl(l)*qz(i,l))
        rho2ns(i, 1) = rho2ns(i, 1) + c0ll*dimag(df*denr(i))
      end do
    end if

!---> calculate the non spherically symmetric contribution
!        to speed up the pointer jend generated in gaunt is used
!        remember that the wavefunctions are l and not lm dependent
    call csimpk(denr, den(l), ipan, ircut, drdi)
  end do
  den(lmaxd1) = 0.0d0





  j0 = 1

  do i = 1, ircut(1)
    denr(i) = 0.0d0
  end do
  do lm3 = 2, lm3max
    do l1 = 0, lmaxd
      do l2 = 0, l1
!---> sum over m1,m2 for fixed lm3,l1,l2
        j1 = jend(lm3, l1, l2)

        if (j1/=0) then

          gmatl = czero



          do j = j0, j1
            facsym = 2.0d0
            ln1 = icleb(j, 1)
            ln2 = icleb(j, 2)
            if (ln1==ln2) facsym = 1.0d0
            gmatl = gmatl + facsym*cleb(j)*df*gmat(ln2, ln1)
          end do

          j0 = j1 + 1

          do i = 2, ircut(1)
            rho2ns(i, lm3) = rho2ns(i, lm3) + dimag(gmatl*wr(i,l1,l2))
          end do

        end if
!-----------------------------------------------------------------------
      end do
!     calculate in the paramagnetic case (nspin=1) :
    end do
!         the valence charge density times r**2 from the greensfunction
  end do
!     calculate in the spin-polarized case (nspin=2) :
end subroutine
