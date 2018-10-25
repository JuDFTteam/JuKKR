!------------------------------------------------------------------------------------
!> Summary: The charge density is developed in spherical harmonics
!> Author: B. Drittler
!> The charge density is developed in spherical harmonics
!> \begin{eqaution}
!> \rho(r) = \sum_{lm} \rho(lm,r) Y(r,lm) }
!> \end{equation}
!> \begin{eqaution}
!> \rho(lm,r) = \int \rho(r) * Y(r,lm) 
!> \end{equation}
!> in the case of spin-polarization :
!> the spin density is developed in spherical harmonics :
!> \begin{eqaution}
!> sden(r) = \sum_{lm}{ sden(lm,r) Y(r,lm) }
!> \end{equation}
!> \begin{eqaution}
!> sden(lm,r) = \int sden(r) T(r,lm)
!> \end{equation}
!> \(n(r,e)\) is developed in
!> \begin{eqaution}
!> n(r,e) = Y(r,l'm') n(l'm',lm,r,e) Y(r,lm)
!> \end{equation}
!> Therefore a faltung of `n(l'm',lm,r,e)` with the gaunt coeffients
!> has to be used to calculate the lm-contribution of the charge density.
!> Calculate the valence density of states , in the spin-polarized case spin dependent.
!> recognize that the density of states is always complex also in the case of 
!> _real-energy-integation_ (`ief>0`) since in that case the energy integration 
!> is done _parallel_ to the real energy axis but *not on the real energy axis*.
!> In the last energy-spin loop the l-contribution of the valence charge is calculated.
!------------------------------------------------------------------------------------
!> @note B. Drittler July 1989: Modified for the use of shape functions
!> @endnote
!> @warning `irmin + 3` has to be less than `imt` if shape functions are used.
!> @endwarning
!------------------------------------------------------------------------------------
module mod_rhons

contains

  !-------------------------------------------------------------------------------
  !> Summary: The charge density is developed in spherical harmonics
  !> Author: B. Drittler
  !> Category: physical-observables, KKRhost
  !> Deprecated: False
  !> The charge density is developed in spherical harmonics
  !> \begin{eqaution}
  !> \rho(r) = \sum_{lm} \rho(lm,r) Y(r,lm) }
  !> \end{equation}
  !> \begin{eqaution}
  !> \rho(lm,r) = \int \rho(r) * Y(r,lm) 
  !> \end{equation}
  !> in the case of spin-polarization :
  !> the spin density is developed in spherical harmonics :
  !> \begin{eqaution}
  !> sden(r) = \sum_{lm}{ sden(lm,r) Y(r,lm) }
  !> \end{equation}
  !> \begin{eqaution}
  !> sden(lm,r) = \int sden(r) T(r,lm)
  !> \end{equation}
  !> \(n(r,e)\) is developed in
  !> \begin{eqaution}
  !> n(r,e) = Y(r,l'm') n(l'm',lm,r,e) Y(r,lm)
  !> \end{equation}
  !> Therefore a faltung of `n(l'm',lm,r,e)` with the gaunt coeffients
  !> has to be used to calculate the lm-contribution of the charge density.
  !> Calculate the valence density of states , in the spin-polarized case spin dependent.
  !> recognize that the density of states is always complex also in the case of 
  !> _real-energy-integation_ (`ief>0`) since in that case the energy integration 
  !> is done _parallel_ to the real energy axis but *not on the real energy axis*.
  !> In the last energy-spin loop the l-contribution of the valence charge is calculated.
  !-------------------------------------------------------------------------------
  !> @note B. Drittler July 1989: Modified for the use of shape functions
  !> @endnote
  !> @warning `irmin + 3` has to be less than `imt` if shape functions are used.
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine rhons(den,df,drdi,gmat,ek,rho2ns,ipan,ircut,irmin,thetas,ifunm,lmsp,   &
    nsra,qns,pns,ar,cr,pz,fz,qz,sz,cleb,icleb,jend,iend,ekl,denlm,gflle_part)

    use :: mod_datatypes
    use :: global_variables
    use :: mod_rhoin
    use :: mod_rhoout
    use :: mod_csimpk
    use :: mod_constants, only: pi
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    complex (kind=dp) :: df, ek
    integer :: iend, ipan, nsra, irmin
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: ar(lmmaxd, lmmaxd), cr(lmmaxd, lmmaxd), den(0:(lmaxd+1)), ekl(0:lmaxd), fz(irmd, 0:lmaxd), gmat(lmmaxd, lmmaxd), pns(lmmaxd, lmmaxd, irmind:irmd, 2), &
      pz(irmd, 0:lmaxd), qns(lmmaxd, lmmaxd, irmind:irmd, 2), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd), denlm(lmmaxd)
#ifndef CPP_MPI
    complex (kind=dp) :: energ     ! lm-dos
#endif
    real (kind=dp) :: cleb(*), drdi(irmd), rho2ns(irmd, lmpotd), thetas(irid, nfund)
    integer :: icleb(ncleb, 4), ifunm(*), ircut(0:ipand), jend(lmpotd, 0:lmaxd, 0:lmaxd), lmsp(*)
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: denns, v1
    integer :: imt1, l, lm, m, irmax, lm1, lm2
    ! ..
    ! .. Local Arrays ..
    complex (kind=dp) :: cden(irmd, 0:lmaxd), cdenns(irmd), efac(lmmaxd), cdenlm(irmd, lmmaxd), cwr(irmd, lmmaxd, lmmaxd) & ! lm-dos
      , gflle_part(lmmaxd, lmmaxd)
    ! ..
    ! .. External Functions ..
    logical :: opt                 ! qdos
    external :: opt                ! qdos

    ! ---> set up efac(lm) = sqrt(e))**l/(2l - 1)!!

    efac(1) = 1.0d0
    v1 = 1.0d0
    do l = 1, lmaxd
      v1 = v1*ek/dble(2*l-1)
      do m = -l, l
        lm = l*(l+1) + m + 1
        efac(lm) = v1
      end do
    end do

    imt1 = ircut(1)
    irmax = ircut(ipan)

    call rhoout(cden, df, gmat, ek, pns, qns, rho2ns, thetas, ifunm, ipan, imt1, irmin, irmax, lmsp, cdenns, nsra, cleb, icleb, & ! Added IRMIN,IRMAX 1.7.2014  &
      iend, cdenlm, cwr)           ! lm-dos

    call rhoin(ar, cden, cr, df, gmat, ek, rho2ns, irmin, nsra, efac, pz, fz, & ! Changed from irmind TO irmin 1.7.2014  &
      qz, sz, cleb, icleb, jend, iend, ekl, cdenlm) ! lm-dos  ! Attention, cwr does not go into rhoin, does lmlm-dos work properly?


    ! ---> calculate complex density of states

    do l = 0, lmaxd

      ! ---> call integration subroutine

      call csimpk(cden(1,l), den(l), ipan, ircut, drdi)
    end do

    do lm1 = 1, lmmaxd             ! lm-dos
      call csimpk(cdenlm(1,lm1), denlm(lm1), ipan, ircut, drdi) ! lm-dos
      if (opt('lmlm-dos') .or. opt('qdos    ') .or. & ! lmlm-dos  &
        opt('LDA+U   ')) then      ! LDAU
        do lm2 = 1, lmmaxd         ! lmlm-dos
          call csimpk(cwr(1,lm1,lm2), gflle_part(lm1,lm2), & ! lmlm-dos  &
            ipan, ircut, drdi)     ! lmlm-dos
        end do
      end if                       ! lmlm-dos
    end do

    ! Energy depends on EK and NSRA:                            ! lm-dos
    ! IF (NSRA.EQ.1) EK = SQRT(E)                           ! lm-dos
    ! IF (NSRA.EQ.2) EK = SQRT(E+E*E/ (CVLIGHT*CVLIGHT))    ! lm-dos
    ! CVLIGHT=274.0720442D0                                 ! lm-dos
    ! Therefore the following is a good approximation           ! lm-dos
    ! for energies of a few Ryd:                                ! lm-dos
#ifndef CPP_MPI
    if (.not. opt('qdos    ')) then
      energ = ek**2                ! lm-dos
      write (30, 100) real(energ, kind=dp), (-aimag(denlm(lm))/pi, lm=1, lmmaxd)
100   format (30e12.4)
    end if                         ! not qdos option
#endif


    if (ipan>1) then
      call csimpk(cdenns, denns, ipan, ircut, drdi)
      den((lmaxd+1)) = denns
    end if

    return
  end subroutine rhons

end module mod_rhons
