!------------------------------------------------------------------------------------
!> Summary: Driving routine to convert the TB-KKR potential from the non-relativistic representation, to the relativistic one
!> Author: B. Popescu
!> A Driving routine to convert the TB-KKR potential from the non-relativistic 
!> representation `VM2Z(IRMD,NPOTD)`, with `IPOTD` the combined index for `ATOM` 
!> and `SPIN` to the relativistic one.
!> \begin{equation}
!> V_{TREL}=\frac{V_{up}+V_{down}}{2}
!> \end{equation}
!> \begin{equation}
!> B_{TREL}=\frac{V_{up}-V_{down}}{2}
!> \end{equation}
!>  Additionally, for compatibility with the relativistic routines included in the 
!> package, `VTREL` includes the Coulomb term, and the auxiliary arrays `ZREL`, 
!> `RMREL`, `JWSREL`, `DRDI`, `R2DRDI` and `IRSHIFT` are created. 
!> `IRSHIFT(NATYPD)` accounts for the shift in the radial mesh, since the first point
!>  (sometimes first two points) of `VM2Z ( = 0D0 )` are skipped.
!> The relativistic routines require an **odd** number of radial points (Simpson integration routine)
!------------------------------------------------------------------------------------
!> @warning 
!> * Because this routine is called only `IF KREL.EQ.0`, the number of spins in `VM2Z` is always 2
!> * So far, only `SPHERICAL` part implemented
!> @endwarning
!------------------------------------------------------------------------------------
module mod_relpotcvt
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Driving routine to convert the TB-KKR potential from the non-relativistic representation, to the relativistic one
  !> Author: B. Popescu
  !> Category: potential, dirac, KKRhost 
  !> Deprecated: False
  !> Driving routine to convert the TB-KKR potential from the non-relativistic 
  !> representation `VM2Z(IRMD,NPOTD)`, with `IPOTD` the combined index for `ATOM` 
  !> and `SPIN` to the relativistic one.
  !> \begin{equation}
  !> V_{TREL}=\frac{V_{up}+V_{down}}{2}
  !> \end{equation}
  !> \begin{equation}
  !> B_{TREL}=\frac{V_{up}-V_{down}}{2}
  !> \end{equation}
  !>  Additionally, for compatibility with the relativistic routines included in the 
  !> package, `VTREL` includes the Coulomb term, and the auxiliary arrays `ZREL`, 
  !> `RMREL`, `JWSREL`, `DRDI`, `R2DRDI` and `IRSHIFT` are created. 
  !> `IRSHIFT(NATYPD)` accounts for the shift in the radial mesh, since the first point
  !>  (sometimes first two points) of `VM2Z ( = 0D0 )` are skipped.
  !> The relativistic routines require an **odd** number of radial points (Simpson integration routine)
  !-------------------------------------------------------------------------------
  !> @warning 
  !> * Because this routine is called only `IF KREL.EQ.0`, the number of spins in `VM2Z` is always 2
  !> * So far, only `SPHERICAL` part implemented
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine relpotcvt(icall,vm2z,zin,rin,drdiin,ircut,vtrel,btrel,zrel,rmrel,      &
    jwsrel,drdirel,r2drdirel,irshift,ipand,irmd,npotd,natypd)

    use :: mod_rinit
    implicit none

    ! PARAMETER definitions
    integer :: nspinpot
    parameter (nspinpot=2)

    ! Scalar arguments
    integer :: icall, ipand, irmd, npotd, natypd

    ! Array arguments
    real (kind=dp) :: vm2z(irmd, npotd)
    real (kind=dp) :: zin(natypd), rin(irmd, natypd)
    real (kind=dp) :: drdiin(irmd, natypd)
    integer :: ircut(0:ipand, natypd)

    real (kind=dp) :: vtrel(irmd, natypd), btrel(irmd, natypd)
    real (kind=dp) :: drdirel(irmd, natypd), r2drdirel(irmd, natypd)
    real (kind=dp) :: rmrel(irmd, natypd)
    integer :: irshift(natypd), jwsrel(natypd), zrel(natypd)

    ! Local scalars
    real (kind=dp) :: vdn, vup
    integer :: it, ir, ip, ipot, ishift, jr

    ! ------------------------------------------------------- INITIALISATION
    if (icall==1) then
      call rinit(irmd*natypd, rmrel)
      call rinit(irmd*natypd, drdirel)
      call rinit(irmd*natypd, r2drdirel)
      do it = 1, natypd
        jwsrel(it) = 0
        irshift(it) = 0
        zrel(it) = 0
      end do
    end if
    call rinit(irmd*natypd, vtrel)
    call rinit(irmd*natypd, btrel)
    ! ------------------------------------------------------- INITIALISATION

    ! *************************************************************** NATYPD
    do it = 1, natypd
      ! ================================================================ ICALL
      ! variables require init only once
      if (icall==1) then

        ! skip first mesh point and also the second if IRCUT(1,IT) = WS-rad odd,
        ! since JWSREL(IT) must be odd

        ishift = 1
        if (mod(ircut(1,it),2)==1) ishift = 2
        ir = 0
        ! ----------------------------------------------------------------------
        do jr = 1 + ishift, ircut(1, it)
          ir = ir + 1
          rmrel(ir, it) = rin(jr, it)
          drdirel(ir, it) = drdiin(jr, it)
          r2drdirel(ir, it) = rmrel(ir, it)*rmrel(ir, it)*drdirel(ir, it)
        end do
        ! ----------------------------------------------------------------------
        jwsrel(it) = ir
        irshift(it) = ishift
        zrel(it) = nint(zin(it))
      end if
      ! ================================================================ ICALL
      ipot = (it-1)*nspinpot + 1
      ishift = irshift(it)
      ! ----------------------------------------------------------------------
      do ir = 1, jwsrel(it)
        ip = ir + ishift
        vdn = -2e0_dp*zin(it)/rin(ip, it) + vm2z(ip, ipot)
        vup = -2e0_dp*zin(it)/rin(ip, it) + vm2z(ip, ipot+1)
        vtrel(ir, it) = (vup+vdn)/2.0e0_dp
        btrel(ir, it) = (vup-vdn)/2.0e0_dp
      end do
      ! ----------------------------------------------------------------------
    end do
    ! *************************************************************** NATYPD
  end subroutine relpotcvt

end module mod_relpotcvt
