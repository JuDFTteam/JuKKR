module mod_rllsllsourceterms
  
  private
  public :: rllsllsourceterms

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of Bessel and Hankel source terms 
  !> Author: 
  !> Category: KKRhost, sigle-site
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Constructs potential matrix (2x2 for SOC) adding SOC potential with proper form 
  !> of small-component in the case of a scalar-relativistic calculation.
  !> Then creates source terms needed to solve Lippmann-Schwinger equations as described 
  !> in the PhD thesis of David Bauer.
  !>
  !> Calculates the source terms J,H and the left solution J2, H2 for:
  !> - non-relativistic
  !> - scalar-relativistic
  !> - full-relativistic
  !> calculations
  !-------------------------------------------------------------------------------
  subroutine rllsllsourceterms(nsra, nvec, eryd, rmesh, nrmax, nrmaxd, lmax, lmsize, use_fullgmat, jlk_index, hlk, jlk, hlk2, jlk2, gmatprefactor)
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: cvlight, ci
    use :: mod_beshank, only: beshank, beshank_smallcomp
    use :: global_variables, only: korbit

    implicit none

    ! inputs
    integer, intent (in) :: nsra, lmax, nrmax, nrmaxd
    integer, intent (in) :: lmsize
    complex (kind=dp), intent (in) :: eryd
    real (kind=dp), dimension (nrmaxd), intent (in) :: rmesh
    integer, intent (in) :: use_fullgmat

    ! outputs
    integer, intent (out) :: nvec
    integer, dimension ((1+korbit)*lmsize), intent (out) :: jlk_index !! index array mapping entries of hlk, jlk (bing/small components one after the other) to L=(l,m,s)
    complex (kind=dp), dimension (1:nsra*(1+korbit)*(lmax+1), nrmax), intent (out) :: hlk, jlk !! right hankel and bessel source functions
    complex (kind=dp), dimension (1:nsra*(1+korbit)*(lmax+1), nrmax), intent (out) :: hlk2, jlk2 !! left hankel and bessel source functions
    complex (kind=dp), intent (out) :: gmatprefactor !! prefactor of the Green function (2M_0\kappa in PhD Bauer, p. 63)

    ! locals
    integer :: l1, lm1, m1, ivec, ispinfullgmat, ir
    complex (kind=dp) :: ek, ek2

    ! just copies value of nsra to nvec
    if (nsra==2) then
      nvec = 2
    else if (nsra==1) then
      nvec = 1
    else
      stop 'error in rllsllsourceterms: nsra is neither 1 nor 2!'
    end if

    lm1 = 1
    do ivec = 1, nvec
      do ispinfullgmat = 0, use_fullgmat
        do l1 = 0, lmax
          do m1 = -l1, l1
            jlk_index(lm1) = l1 + (ivec-1)*(lmax+1) + 1
            lm1 = lm1 + 1
          end do ! m1
        end do ! l1
      end do ! ispinfullgmat
    end do ! nvec

    ! for BdG: here ek, ek2 have to be modified to get source terms for e/h parts that use different kappa=sqrt(E +/- E_F) instead of sqrt(E) (+ some relativistic corrections)
    ! then benhank and beshank_smallcomp should work the same way
    ! add additional loop ofer e/h index and take care of nsra cases below

    if (nsra==1) then
      ek = sqrt(eryd)
      ek2 = sqrt(eryd)
    else if (nsra==2) then
      ek = sqrt(eryd+(eryd/cvlight)**2)
      ek2 = sqrt(eryd+(eryd/cvlight)**2)*(1.0e0_dp+eryd/cvlight**2)
    end if

    do ir = 1, nrmax

      call beshank(hlk(:,ir), jlk(:,ir), ek*rmesh(ir), lmax)
      if (nsra==2) then
        call beshank_smallcomp(hlk(:,ir), jlk(:,ir), ek*rmesh(ir), rmesh(ir), eryd, lmax)
      end if

      ! Attention: here the different definition of Drittler (see Drittler PhD p. 18) for the sperical hankel function is used which gives the additional factor -i
      ! this factor is added here
      do l1 = 1, nvec*(lmax+1)
        hlk(l1, ir) = -ci*hlk(l1, ir)
      end do

      ! use symmetries to get left solutions (minus sign only for NSRA==2 and l1>lmax+1)
      if (nsra==1) then
        do l1 = 1, nvec*(lmax+1)
          jlk2(l1, ir) = jlk(l1, ir)
          hlk2(l1, ir) = hlk(l1, ir)
        end do
      else if (nsra==2) then
        do l1 = 1, lmax + 1
          jlk2(l1, ir) = jlk(l1, ir)
          hlk2(l1, ir) = hlk(l1, ir)
        end do
        do l1 = lmax + 2, (1+korbit)*(lmax+1)
          jlk2(l1, ir) = -jlk(l1, ir)
          hlk2(l1, ir) = -hlk(l1, ir)
        end do
      end if

    end do

    ! store prefactor for Green function, used later on (e.g. equation 3.34 on page 21 of PhD Drittler)
    gmatprefactor = ek2

  end subroutine rllsllsourceterms

end module mod_rllsllsourceterms
