!------------------------------------------------------------------------------------
!> Summary: Calculate the l-dependent charge density for spherical potential
!> Author: B. Drittler
!> Calculate the l-dependent charge density for spherical potential.
!> * `nspin=1`: The valence charge density times \(r^2\) from the Greens function
!> * `nspin=2`: The valence charge density times \(r^2\) and the valence spin
!> density times \(r^2\) from the greensfunction (convention spin density :=
!> density(spin up)-density(spin down) ) 
!>
!> Calculate the valence density of states, in the spin-polarized case spin dependent;
!> splitted into its l-contributions.
!> 
!> In this subroutine an implicit energy-spin integration is done.
!> This subroutine is called for each energy and spin value and \(n(r,e)\) times 
!> \(df\) (the energy weight) is calculated.
!> Recognize that the density of states is always complex also in the case of 
!> _real-energy-integation_ `(ief>0)` since in that case  the energy integration is 
!> done _parallel to the real energy axis_ but *not on the real energy axis*.
!> In the paramagnetic case only `rho2ns(irmd,lmxtsq,natypd,1)` is used containing 
!> the charge density times \(r^2\) .
!> In the spin-polarized case `rho2ns(...,1)` contains the charge density times \(r^2\)
!> and `rho2ns(...,2)` the spin density times \(r^2\).
!> The charge density is expanded in spherical harmonics:
!> \begin{equation}
!> \rho(r) =   \sum_{l,m}\rho(lm,r) Y(r,lm)
!> \end{equation}
!> \begin{equation}
!> \rho(lm,r) =  \int \rho(r) Y(r,lm)
!> \end{equation}
!> In the case of spin-polarization, the spin density is developed in spherical harmonics:
!> \begin{equation}
!> sden(r) = \sum_{lm} sden(lm,r) Y(r,lm)
!> \end{equation}
!> \begin{equation}
!> sden(lm,r) = \int sden(r) * y(r,lm) 
!> \end{equation}
!> \(n(r,e)\) is developed in
!> \begin{equation}
!> n(r,e) =  Y(r,l'm') n(l'm',lm,r,e) Y(r,lm)
!> \end{equation}
!> therefore a faltung of `n(l'm',lm,r,e)` with the gaunt coeffients
!> has to be used to calculate the lm-contribution of the charg density .
!> (see notes by b.drittler)
!------------------------------------------------------------------------------------
!> @warning The gaunt coeffients are stored in an index array (see subroutine `gaunt`)
!> the structure part of the greens-function (`gmat`) is symmetric in its lm-indices,
!> therefore only one half of the matrix is calculated in the subroutine for the 
!> back-symmetrisation. The gaunt coeffients are symmetric too (since the are calculated for
!> real spherical harmonics). That is why the `lm2`-loop only goes up to `lm1` and the summands are
!> multiplied by a factor of 2 in the case of `lm1` not equal to `lm2`.
!> @endwarning
!------------------------------------------------------------------------------------
module mod_rholm

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate the l-dependent charge density for spherical potential
  !> Author: B. Drittler
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculate the l-dependent charge density for spherical potential.
  !> * `nspin=1`: The valence charge density times \(r^2\) from the Greens function
  !> * `nspin=2`: The valence charge density times \(r^2\) and the valence spin
  !> density times \(r^2\) from the greensfunction (convention spin density :=
  !> density(spin up)-density(spin down) ) 
  !>
  !> Calculate the valence density of states, in the spin-polarized case spin dependent;
  !> splitted into its l-contributions.
  !> 
  !> In this subroutine an implicit energy-spin integration is done.
  !> This subroutine is called for each energy and spin value and \(n(r,e)\) times 
  !> \(df\) (the energy weight) is calculated.
  !> Recognize that the density of states is always complex also in the case of 
  !> _real-energy-integation_ `(ief>0)` since in that case  the energy integration is 
  !> done _parallel to the real energy axis_ but *not on the real energy axis*.
  !> In the paramagnetic case only `rho2ns(irmd,lmxtsq,natypd,1)` is used containing 
  !> the charge density times \(r^2\) .
  !> In the spin-polarized case `rho2ns(...,1)` contains the charge density times \(r^2\)
  !> and `rho2ns(...,2)` the spin density times \(r^2\).
  !> The charge density is expanded in spherical harmonics:
  !> \begin{equation}
  !> \rho(r) =   \sum_{l,m}\rho(lm,r) Y(r,lm)
  !> \end{equation}
  !> \begin{equation}
  !> \rho(lm,r) =  \int \rho(r) Y(r,lm)
  !> \end{equation}
  !> In the case of spin-polarization, the spin density is developed in spherical harmonics:
  !> \begin{equation}
  !> sden(r) = \sum_{lm} sden(lm,r) Y(r,lm)
  !> \end{equation}
  !> \begin{equation}
  !> sden(lm,r) = \int sden(r) * y(r,lm) 
  !> \end{equation}
  !> \(n(r,e)\) is developed in
  !> \begin{equation}
  !> n(r,e) =  Y(r,l'm') n(l'm',lm,r,e) Y(r,lm)
  !> \end{equation}
  !> therefore a faltung of `n(l'm',lm,r,e)` with the gaunt coeffients
  !> has to be used to calculate the lm-contribution of the charg density .
  !> (see notes by b.drittler)
  !-------------------------------------------------------------------------------
  !> @warning The gaunt coeffients are stored in an index array (see subroutine `gaunt`)
  !> the structure part of the greens-function (`gmat`) is symmetric in its lm-indices,
  !> therefore only one half of the matrix is calculated in the subroutine for the 
  !> back-symmetrisation. The gaunt coeffients are symmetric too (since the are calculated for
  !> real spherical harmonics). That is why the `lm2`-loop only goes up to `lm1` and the summands are
  !> multiplied by a factor of 2 in the case of `lm1` not equal to `lm2`.
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine rholm(den,df,gmat,nsra,rho2ns,drdi,ipan,ircut,pz,fz,qz,sz,cleb,icleb,  &
    iend, jend, ekl)

    use :: global_variables
    use :: mod_datatypes, only: dp
    use :: mod_csimpk
    use :: constants, only: czero,pi
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: df
    integer :: iend, ipan, nsra
    ! ..
    ! .. Local Arrays ..
    complex (kind=dp) :: den(0:(lmaxd+1)), ekl(0:lmaxd), fz(irmd, 0:lmaxd), gmat(lmmaxd, lmmaxd), pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd)
    real (kind=dp) :: cleb(*), drdi(irmd), rho2ns(irmd, lmpotd)
    integer :: icleb(ncleb, 4), ircut(0:ipand), jend(lmpotd, 0:lmaxd, 0:lmaxd)
    ! ..
    complex (kind=dp) :: ffz, gmatl, ppz
    real (kind=dp) :: c0ll, facsym
    integer :: i, j, j0, j1, l, l1, l2, lm3, lm3max, ln1, ln2, lne, lns
    ! ..
    complex (kind=dp) :: denr(irmd), wr(irmd, 0:lmaxd, 0:lmaxd)
    ! ..

    c0ll = 1.0e0_dp/sqrt(4.0e0_dp*pi)
    ! ---> set up of wr(ir,l1,l2) = pz(ir,l1)*pz(ir,l2)

    lm3max = icleb(iend, 3)

    if (nsra==2) then
      do l1 = 0, lmaxd
        do l2 = 0, l1
          do i = 2, ircut(1)
            wr(i, l1, l2) = pz(i, l1)*pz(i, l2) + fz(i, l1)*fz(i, l2)
          end do
        end do
      end do
      ! ---> first calculate only the spherically symmetric contribution
    else
      do l1 = 0, lmaxd
        do l2 = 0, l1
          do i = 2, ircut(1)
            wr(i, l1, l2) = pz(i, l1)*pz(i, l2)
          end do
        end do
      end do
    end if
    ! ---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
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
          denr(i) = ppz*(gmatl*ppz+ekl(l)*qz(i,l)) + ffz*(gmatl*ffz+ekl(l)*sz(i,l))
          rho2ns(i, 1) = rho2ns(i, 1) + c0ll*aimag(df*denr(i))
        end do
        ! ---> calculate density of states
      else

        do i = 2, ircut(1)
          ppz = pz(i, l)
          denr(i) = ppz*(gmatl*ppz+ekl(l)*qz(i,l))
          rho2ns(i, 1) = rho2ns(i, 1) + c0ll*aimag(df*denr(i))
        end do
      end if
      ! ---> calculate the non spherically symmetric contribution
      ! to speed up the pointer jend generated in gaunt is used
      ! remember that the wavefunctions are l and not lm dependent
      call csimpk(denr, den(l), ipan, ircut, drdi)
    end do
    den((lmaxd+1)) = 0.0e0_dp

    j0 = 1
    do i = 1, ircut(1)
      denr(i) = 0.0e0_dp
    end do
    do lm3 = 2, lm3max
      do l1 = 0, lmaxd
        do l2 = 0, l1
          ! ---> sum over m1,m2 for fixed lm3,l1,l2
          j1 = jend(lm3, l1, l2)
          if (j1/=0) then
            gmatl = czero
            do j = j0, j1
              facsym = 2.0e0_dp
              ln1 = icleb(j, 1)
              ln2 = icleb(j, 2)
              if (ln1==ln2) facsym = 1.0e0_dp
              gmatl = gmatl + facsym*cleb(j)*df*gmat(ln2, ln1)
            end do
            j0 = j1 + 1
            do i = 2, ircut(1)
              rho2ns(i, lm3) = rho2ns(i, lm3) + aimag(gmatl*wr(i,l1,l2))
            end do
          end if
          ! -----------------------------------------------------------------------
        end do
        ! calculate in the paramagnetic case (nspin=1) :
      end do
      ! the valence charge density times r**2 from the greensfunction
    end do
    ! calculate in the spin-polarized case (nspin=2) :
  end subroutine rholm

end module mod_rholm
