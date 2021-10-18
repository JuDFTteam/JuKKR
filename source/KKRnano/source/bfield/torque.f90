!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculation of the magnetic torque methods
!> 
!> Author: Sascha Brinker, Nicolas Essing
!------------------------------------------------------------------------------------
module mod_torque
  use mod_bfield, only: bfield_data
  use NonCollinearMagnetism_Helpers_mod, only: rotatematrix, intcheb_cell

  implicit none

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the density for the new solver
  !> Author: Sascha Brinker, Nicolas Essing
  !-------------------------------------------------------------------------------
  subroutine calc_torque(bfield, vpot, rho2nsc, theta, phi, lmax, rpan_intervall, ipan_intervall, &
                         npan_tot, ncheb, imt, iend, icleb, cleb, ifunm, thetasnew, &
                         lbfield_constr, lbfield_mt, itscf0, itscf1, iteration, constr_mode)
    type(bfield_data), intent(inout) :: bfield !! Information on the magnetic field
    double precision, dimension(:,:,:), intent(in) :: vpot    !! The potential
    double complex,   dimension(:,:,:), intent(in) :: rho2nsc !! complex density matrix
    double precision, intent(in) :: theta, phi !! Angles of the (old) magnetic moment
    integer,                           intent(in) :: lmax            !! Angular momentum cutoff
    double precision, dimension(:),    intent(in) :: rpan_intervall
    integer,          dimension(:),    intent(in) :: ipan_intervall
    integer,                           intent(in) :: npan_tot
    integer,                           intent(in) :: ncheb
    integer,                           intent(in) :: imt       !! Index of the muffin tin radius
    integer,                           intent(in) :: iend      !! Number of nonzero gaunt coefficients
    integer,          dimension(:,:),  intent(in) :: icleb     !! l and m values for the Gaunt coefficients
    double precision, dimension(:),    intent(in) :: cleb      !! Gaunt coefficients
    integer         , dimension(:),    intent(in) :: ifunm     !! pointer array for shapefun
    double precision, dimension(:, :), intent(in) :: thetasnew !! shapefun on the Cheby mesh
    logical, intent(in) :: lbfield_constr  !! Use constraint fields
    logical, intent(in) :: lbfield_mt      !! Use magnetic fields only inside the muffin tin
    integer, intent(in) :: itscf0, itscf1  !! Apply magnetic fields between these iterations
    integer, intent(in) :: iteration       !! Current iteration
    integer, intent(in) :: constr_mode     !! Mode of the constraining field self-consistency

    double complex, parameter :: cone = (1., 0.)             ! complex one
    double precision, parameter :: rfpi = 3.5449077018110318 ! sqrt(4*pi)

    double precision, parameter :: constraint_bfields_mixing_parameter = 0.03

    integer :: irmd, lmmax, lmpotd     ! radial and angular momentum sizes
    integer :: i, ir, ilm              ! generic, radial and angular momentum loop indices
    double complex :: integrate_result ! to output the result of some integrations
    double precision :: mag_mom_len    ! absolute value of magnetic moment
    double precision, dimension(3) :: dir, mag_mom_dir  ! unit vectors of old and new magnetic moment
    double precision, dimension(3) :: torque, torque_mt ! torque calculated over whole cell or mt only
    double precision, dimension(3) :: old_b_constr      ! temporary storage for the old constraint field
    double complex,   dimension(2,2) :: rho2ns_temp     ! temporary matrix to rotate rho2nsc
    double complex,   dimension(:),     allocatable :: integrand    ! temporary array to integrate stuff
    double precision, dimension(:,:),   allocatable :: bxc          ! xc magnetic field in collinear form
    double precision, dimension(:,:,:), allocatable :: mag_den      ! magnetization density
    double precision, dimension(:,:,:), allocatable :: mag_den_conv ! mag. den. convoluted with shapefun
    double precision, dimension(:,:),   allocatable :: mag          ! integrated lm-resolved magnetic moment
    double precision, dimension(:,:),   allocatable :: mag_mt       ! integrated lm-resolved magnetic moment up to mt

    ! Get some dimensions
    lmmax = (lmax+1)**2
    irmd = size(vpot,1)
    lmpotd = size(vpot,2)

    ! Allocate temporary arrays
    allocate(integrand(irmd))
    allocate(bxc(irmd, lmpotd))
    allocate(mag_den(irmd, lmpotd, 3))
    allocate(mag(lmmax, 3))
    allocate(mag_mt(lmmax, 3))

    ! Get the xc magnetic field. In the rigid spin approximation, the direction
    ! of the field is fixed by the local frame, but it is radial dependent.
    ! Can be calculated as the difference of spin up and spin down part of the
    ! potential.
    bxc(:,:) = (vpot(:,:,1) - vpot(:,:,2)) / 2.

    ! Get magnetization density (in the global frame). For each angular momentum
    ! and radial index, rotate to the global frame.
    do ilm = 1, lmpotd
      do ir = 1, irmd
        rho2ns_temp(1,1) = rho2nsc(ir,ilm,1)
        rho2ns_temp(2,2) = rho2nsc(ir,ilm,2)
        rho2ns_temp(1,2) = rho2nsc(ir,ilm,3)
        rho2ns_temp(2,1) = rho2nsc(ir,ilm,4)

        call rotatematrix(rho2ns_temp, theta, phi, 1, 0)

        mag_den(ir,ilm,1) = aimag( rho2ns_temp(1,2) + rho2ns_temp(2,1) ) 
        mag_den(ir,ilm,2) = real(  rho2ns_temp(2,1) - rho2ns_temp(1,2) )
        mag_den(ir,ilm,3) = aimag( rho2ns_temp(2,2) - rho2ns_temp(1,1) )
      end do
    end do
    
    ! Integrate magnetization density (to mt or end depending on parameter)
    ! to get magnetic moment. First convolute with the shapefun.
    call calc_mag_mom(lmax, imt, iend, icleb, cleb, ifunm, thetasnew, mag_den, mag_den_conv)
    do i = 1, 3
      do ilm = 1, lmmax
          integrand(:) = mag_den_conv(:, ilm, i) * cone
          call intcheb_cell(integrand,integrate_result,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmd)
          mag(ilm,i) = real(integrate_result * rfpi)    
          ! Same for only the muffin tin: Set to zero outside, integrate again
          integrand(imt:) = (0., 0.)
          call intcheb_cell(integrand(:),integrate_result,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmd)
          mag_mt(ilm,i) = real(integrate_result * rfpi)    
      end do
    end do
    
    !TODO Output magnetic moments

    ! Calculate direction and absolute value of the magnetic moment
    if (lbfield_mt) then
      mag_mom_len = sqrt(dot_product(mag_mt(1,:), mag_mt(1,:)))
      mag_mom_dir = mag_mt(1,:) / mag_mom_len
    else
      mag_mom_len = sqrt(dot_product(mag(1,:), mag(1,:)))
      mag_mom_dir = mag(1,:) / mag_mom_len
    end if

    ! Calculate the torque (also mt or end)
    ! Integrate the xc bfield times the magnetization density.
    torque(:) = 0.
    torque_mt(:) = 0.
    do i = 1, 3
      do ilm = 1, lmmax
        integrand(:) = bxc(:, ilm) * mag_den(:, ilm, i) * cone
        call intcheb_cell(integrand(:),integrate_result,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmd)
        torque(i) = torque(i) + real(integrate_result * rfpi) ! rfpi only contained in convoluted quantities
        ! Same for mt
        integrand(imt:) = 0.
        call intcheb_cell(integrand(:),integrate_result,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmd)
        torque_mt(i) = torque_mt(i) + real(integrate_result * rfpi)
      end do
    end do
    
    ! Project torque to perpendicular of magnetic moment
    dir(1) = sin(theta) * cos(phi)
    dir(2) = sin(theta) * sin(phi)
    dir(3) = cos(theta)
    torque(:)    = torque(:)    - dir(:) * dot_product(dir(:), torque(:))
    torque_mt(:) = torque_mt(:) - dir(:) * dot_product(dir(:), torque_mt(:))

    !TODO Output torque

    ! Save torque in bfield_data
    if (lbfield_mt) then
      bfield%mag_torque(:) = torque_mt
    else
      bfield%mag_torque(:) = torque
    end if
    
    ! Scf-cycle for constraint fields, based either on torque or on fields alone
    if (lbfield_constr .and. itscf0 <= iteration .and. iteration <= itscf1) then
      if (constr_mode == 3 .or. constr_mode == 5) then
        bfield%bfield_constr(:) = bfield%bfield_constr(:) - torque(:) / mag_mom_len
      else if (constr_mode == 2 .or. constr_mode == 4) then
        old_b_constr = bfield%bfield_constr(:)
        bfield%bfield_constr(:) = old_b_constr - dot_product(old_b_constr,dir)*dir - &
                (mag_mom_dir - dot_product(mag_mom_dir,dir)*dir)*constraint_bfields_mixing_parameter
      else
        ! There might be other modes that are calculated somewhere else. Do nothing.
      end if
    end if

    ! Deallocate
    deallocate(bxc, integrand, mag_den, mag_den_conv, mag, mag_mt)

  end subroutine calc_torque

  !------------------------------------------------------------------------------------
  !> Summary: Convolute the magnetic moments with the shapefunction.
  !> 
  !> Author: Sascha Brinker, Nicolas Essing
  !------------------------------------------------------------------------------------
  subroutine calc_mag_mom(lmax, imt, iend, icleb, cleb, ifunm, thetasnew, mag_den, mag_den_conv)
    integer,                            intent(in)  :: lmax      !! Angular momentum cut-off
    integer,                            intent(in)  :: imt       !! index muffin-tin radius
    integer,                            intent(in)  :: iend      !! Number of nonzero gaunt coefficients
    integer,          dimension(:,:),   intent(in)  :: icleb     !! l and m values for the Gaunt coefficients
    double precision, dimension(:),     intent(in)  :: cleb      !! Gaunt coefficients
    integer         , dimension(:),     intent(in)  :: ifunm     !! pointer array for shapefun
    double precision, dimension(:, :),  intent(in)  :: thetasnew !! shapefun on the Cheby mesh
    double precision, dimension(:,:,:), intent(in)  :: mag_den   !! uncovoluted magnetization density
    double precision, dimension(:,:,:), allocatable, intent(out) :: mag_den_conv !! result

    double precision, parameter :: rfpi = 3.5449077018110318 ! sqrt(4*pi)
    double precision, parameter :: c0ll = 1.d0 / rfpi
    integer :: lmmax, lmpotd  ! Angular momentum sizes
    integer :: irmdnew        ! Radial size
    integer :: i, j           ! Generic loop indices
    integer :: lm1, lm2, lm3  ! Angular momentum indices
    integer :: ifun           ! Pointer index for thetasnew
    double precision, dimension(:,:), allocatable :: shapefun_mod ! Complete shapefun

    ! Get some dimensions
    lmmax   = (lmax+1)**2
    lmpotd  = size(mag_den, 2) 
    irmdnew = size(thetasnew, 1) ! npan_tot*(ncheb+1)

    ! Allocate temporary shapefun and the result
    allocate(shapefun_mod(irmdnew, lmpotd))
    allocate(mag_den_conv(irmdnew, lmmax, 3))

    ! Build the shapefun array. Start with all zero. Inside muffin tin only l=0
    ! component is /= 0 (and constant), copy l=0 component for r outside of mt
    ! from thetasnew.
    shapefun_mod(:,:)       = 0.d0
    shapefun_mod(1:imt, 1)  = rfpi ! is multipled by C_LL^0
    shapefun_mod(imt+1:, 1) = thetasnew(imt+1:, 1)

    ! Copy other components from thetasnew. Convert from pointer indices to
    ! normal (l,m)-index
    do lm1 = 2, lmpotd
      ifun = ifunm(lm1)
      if(ifun /= 0) then !shapefun%lmused(lm1)==1) then
        shapefun_mod(imt+1:, lm1) = thetasnew(imt+1:, ifun)
      end if
    end do

    ! Initialize convoluted magnetization density.
    mag_den_conv(:,:,:) = 0.

    ! Diagonal part (not contained in gaunt-coeff)
    do i = 1, 3
      do lm1 = 1, lmpotd
        mag_den_conv(:,1,i) = mag_den_conv(:,1,i) + mag_den(:,lm1,i) * shapefun_mod(:,lm1) * c0ll
      end do
    end do

    ! Offdiagonal part. This is effectively a loop over angular momentum indices
    ! lm1,lm2,lm3.  Iterate instead over the flattened array cleb
    ! containing the Gaunt coefficients and extract the angular momentum
    ! indices for each j.
    do i = 1, 3
      do j = 1, iend
        lm1 = icleb(j, 1)
        lm2 = icleb(j, 2)
        lm3 = icleb(j, 3)
        if (lm1 > lmmax .or. lm2 > lmmax .or. lm3 > lmmax) then
          !TODO I think this should not happen, as icleb should only contain
          ! valid values. It might be defined for higher lm3. Check that.
          cycle
        end if
        mag_den_conv(:,lm3,i) = mag_den_conv(:,lm3,i) + cleb(j) * mag_den(:,lm1,i) * shapefun_mod(:,lm2)
        mag_den_conv(:,lm3,i) = mag_den_conv(:,lm3,i) + cleb(j) * mag_den(:,lm2,i) * shapefun_mod(:,lm1)
      end do
    end do

    ! Deallocate temporary shapefun
    deallocate(shapefun_mod)

  end subroutine

end module mod_torque
