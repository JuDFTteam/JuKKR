!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculation of the magnetic torque methods
!> Author: Sascha Brinker
!> 
!> 
!> 
!> 
!> 
!------------------------------------------------------------------------------------
module mod_torque

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the density for the new solver
  !> Author: Sascha Brinker
  !> Category: KKRhost
  !> Deprecated: False 
  !> Calculation of the magnetic torque in the new solver
  !> The routine uses rho2nsc to calculate the magnetic moment in the global
  !> frame from which the torque is calculated using the xc magnetic field.
  !> The XC magnetic field is extracted from the potential (convoluted with the
  !> shape function) thus rho2nsc has to be used instead of cden and cdenns!
  !-------------------------------------------------------------------------------
  subroutine calc_torque(iatom, lmax, irmdnew, nspin, bfield, rpan_intervall, &
                         ipan_intervall, npan_log, npan_eq, npan_tot, ncheb, &
                         theta, phi, rho2nsc, vpot, ifunm, iend, icleb, cleb, &
                         thetasnew, bconstr, fix_direction, icell)
    
    use :: global_variables, only: lmmaxd, ncleb, ntotd, nfund, korbit, lmpotd
    use :: mod_datatypes, only: dp
    use :: mod_bfield, only: type_bfield
    use :: mod_intcheb_cell, only: intcheb_cell
    use :: mod_rotatespinframe, only: rotatematrix!, rotatevector
    use :: mod_types, only: t_inc

    implicit none

    integer, intent (in)                                                        :: iatom  !! Currently treated atom index
    integer, intent (in)                                                        :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in)                                                        :: irmdnew
    integer, intent (in)                                                        :: nspin
    type(type_bfield), intent(inout)                                            :: bfield !! Information about the external and constraint magnetic fields
    real (kind=dp), dimension (0:ntotd), intent (in)                            :: rpan_intervall
    integer, dimension (0:ntotd), intent (in)                                   :: ipan_intervall
    integer, intent (in)                                                        :: npan_log
    integer, intent (in)                                                        :: npan_eq
    integer, intent (in)                                                        :: npan_tot
    integer, intent (in)                                                        :: ncheb  !! Number of Chebychev pannels for the new solver
    real (kind=dp), intent(in)                                                  :: theta ! polar angle of the local frame
    real (kind=dp), intent(in)                                                  :: phi   ! azimuthal angle of the local frame
    complex (kind=dp), dimension (irmdnew,lmpotd        , 2*nspin), intent (in) :: rho2nsc ! complex density matrix (does not contain the shapefun!)
    real (kind=dp), dimension (irmdnew,lmpotd       ,nspin)                     :: vpot  ! vpot in the new mesh!! corresponds to the vinsnew in main1c but only for the two spin channels and vins in rhovalnew!
    integer              , dimension (1:(2*lmax+1)**2)            , intent (in) :: ifunm        ! pointer array for shapefun     ! Check index and dimensions!!!!!
    integer              , intent (in)                                          :: iend         ! Number of nonzero gaunt coefficients
    integer, dimension (ncleb, 4), intent (in)                                  :: icleb !! Pointer array
    real (kind=dp), dimension (ncleb), intent (in)                              :: cleb !! GAUNT coefficients (GAUNT)    ! CHECK THE DIMENSION AND HOW IT IS USED!!!
    real (kind=dp)       , dimension (irmdnew, nfund)   , intent (in)           :: thetasnew    ! shapefun on the Cheby mesh
    real (kind=dp)       , dimension (4)   , intent (out)                       :: bconstr    ! variables to be saved
    logical, intent(in)                                                         :: fix_direction ! Constrain the direction of the magnetic moment of the current atom
    integer, intent(in)                                                         :: icell
    
    !------------------------------------------------------------------------------------
    ! local variables
    !------------------------------------------------------------------------------------
    integer                      :: ispin
    integer                      :: ir
    integer                      :: ilm
    integer                      :: i
    integer                      :: j
    integer                      :: lm1
    integer                      :: lm2
    integer                      :: lm3
    integer                      :: ifun
    integer                      :: imt1  ! muffin-tin radius in the new mesh
    integer                      :: lmmax
    real(kind=dp),    parameter  :: small = 1.d-6
    real(kind=dp)                :: tempreal
    real(kind=dp)                :: tot_mag_moment
    double complex               :: cone
    double complex               :: temp
    double complex               :: temp2
    real(kind=dp),parameter      :: rfpi=3.5449077018110318
    real(kind=dp)                :: totmag  !! total magnetic moment
    real(kind=dp)                :: totmagmoment  !! total magnetic moment
    real(kind=dp)                :: bfac
    double complex,dimension(2,2)                               :: rho2ns_temp
    double complex,dimension(2,2)                               :: vpot_tmp
    real(kind=dp),dimension(3)                                  :: torque    ! magnetic torque
    real(kind=dp),dimension(3)                                  :: torque_mt ! magnetic torque calcualte from mt only
    real(kind=dp),dimension(3)                                  :: magdir !initial magnetization direction
    real(kind=dp),dimension(3)                                  :: magmoment  !! magnetic moment in local frame
    real(kind=dp),dimension(3)                                  :: c_old
    real(kind=dp),dimension(3)                                  :: magdir_it
    real(kind=dp),dimension(3)                                  :: magdir_old
    real(kind=dp),dimension(1:irmdnew,1:lmpotd)          :: bxc
    real(kind=dp),dimension(1:irmdnew,1:lmpotd,1:3)      :: mag_den_glob ! magnetization density in in the global frame (no shapefun)
    real(kind=dp),dimension(1:irmdnew,1:(lmax+1)**2 ,1:3)      :: mag_den_convol ! magnetization density convoluted with the shapefunction
    real(kind=dp),dimension(1:(lmax+1)**2 ,1:3)                :: mag ! integrated lm-resolved magnetic moment
    real(kind=dp),dimension(1:(lmax+1)**2 ,1:3)                :: mag_mt ! integrated lm-resolved magnetic moment in the muffin-tin
    double complex,dimension(1:irmdnew)                         :: integrand
    !------------------------------------------------------------------------------------
    
    !write(*,'(" >>>>>> entering calctorque >>>>>>>>>")')
    if (t_inc%i_write>1) write(1337,'("===============================================================================")')
    if (t_inc%i_write>1) write(1337,'("                      Magnetic torques for atom ",i4)') iatom
    if (t_inc%i_write>1) write(1337,'("===============================================================================")')
    
    lmmax=(lmax+1)**2
    cone  = (1d0,0d0)
    imt1 = ipan_intervall(npan_log + npan_eq) + 1
    
    !write(*,'(" imt1, irmdnew, lmmax",2i4)') imt1, irmdnew, lmmax
    
    ! bxc is constructed from the potential in the new mesh 
    bxc(:,:) = 0.d0
    bxc(:,:) = (vpot(:,:,1)-vpot(:,:,2))/2.d0
    !do ilm= 1,lmpotd
    !  write(*,'("ilm,vpot1=",i4,1000es16.8)') ilm, vpot(:,ilm,1)
    !  write(*,'("ilm,vpot2=",i4,1000es16.8)') ilm, vpot(:,ilm,2)
    !end do
    !do ilm= 1,lmpotd
    !  write(*,'("ilm,bxc=",i4,1000es16.8)') ilm, bxc(:,ilm)
    !end do
    
    !rotate rho2nsc to the global frame and calculate the magnetization density
    mag_den_glob(:,:,:) = 0.d0
    do ilm=1,lmpotd
      do ir=1,irmdnew
        rho2ns_temp(1,1)            = rho2nsc(ir,ilm,1)
        rho2ns_temp(2,2)            = rho2nsc(ir,ilm,2)
        rho2ns_temp(1,2)            = rho2nsc(ir,ilm,3)
        rho2ns_temp(2,1)            = rho2nsc(ir,ilm,4)
        call rotatematrix(rho2ns_temp, theta, phi, 1, 0)
        mag_den_glob(ir,ilm,1)      = aimag( rho2ns_temp(1,2) + rho2ns_temp(2,1) ) 
        mag_den_glob(ir,ilm,2)      = real( rho2ns_temp(2,1) - rho2ns_temp(1,2) )
        mag_den_glob(ir,ilm,3)      = aimag( rho2ns_temp(2,2) - rho2ns_temp(1,1) )
      end do !ir
    end do !ilm

    ! --------------------------------------------------------------------------------------
    ! Mangetic moment in the new mesh (mainly to separate the moment in r<rmt
    ! and get acess to lm components
    mag_den_convol(:,:,:)=0.d0
    call calc_magmom_newmesh(lmax,imt1,irmdnew,icell,iend,ifunm,icleb,cleb,thetasnew,mag_den_glob,mag_den_convol)
    do i= 1,3
      do ilm=1,lmmax 
          integrand(:) = (0.d0,0.d0)
          temp= (0.d0,0.d0)
          integrand(:) = mag_den_convol(:,ilm,i)*cone
          call intcheb_cell(integrand(:),temp,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
          mag(ilm,i) = real(temp*rfpi)    
          ! same for mt only
          integrand(:) = (0.d0,0.d0)
          temp= (0.d0,0.d0)
          integrand(:imt1) = mag_den_convol(:imt1,ilm,i)*cone
          call intcheb_cell(integrand(:),temp,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
          mag_mt(ilm,i) = real(temp*rfpi)    
      end do
    end do

    if (t_inc%i_write>1) then
      write(1337,'("Magnetic moment")') 
      do ilm=1,lmmax  
        write(1337,'("iatom, ilm, mag, mag_mt",2i4,6es16.8)') iatom, ilm, mag(ilm,:), mag_mt(ilm,:)
      end do
    end if

    ! MdSD: muffin-tin or full cell, these will be used below
    if (bfield%lbfield_mt) then
      totmag = sqrt(dot_product(mag_mt(1,:),mag_mt(1,:)))
      magdir_it = mag_mt(1,:)/totmag
    else
      totmag = sqrt(dot_product(mag(1,:),mag(1,:)))
      magdir_it = mag(1,:)/totmag
    end if
    ! MdSD: constraining fields
    bconstr(1:3) = 0.0_dp; bconstr(4) = totmag

    ! --------------------------------------------------------------------------------------
    
    !do ilm= 1,lmpotd
    !  !write(*,'("ilm,mag_den_glob1=",i4,1000es16.8)') ilm, mag_den_glob(:,ilm,1)
    !  !write(*,'("ilm,mag_den_glob2=",i4,1000es16.8)') ilm, mag_den_glob(:,ilm,2)
    !  write(*,'("ilm,mag_den_glob3=",i4,1000es16.8)') ilm, mag_den_glob(:,ilm,3)
    !end do

    ! calculate the torque
    ! sum_L int dr r^2 m_L(r) B_L^xc(r)
    torque(:)=0.d0
    torque_mt(:)=0.d0
    do i=1,3
      if ( .False. ) then
        integrand(:) = (0.d0,0.d0)
        integrand(:) = bxc(:,1)*mag_den_glob(:,1,i)*cone
        call intcheb_cell(integrand(:),temp,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
        torque(i) = Dreal(temp*rfpi) ! this is a spherical quantity =>sqrt(4*pi)
      else
        do ilm=1,lmmaxd
          integrand(:) = (0.d0,0.d0)
          temp= (0.d0,0.d0)
          integrand(:) = bxc(:,ilm)*mag_den_glob(:,ilm,i)*cone
          call intcheb_cell(integrand(:),temp,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
          torque(i) = torque(i) + real(temp*rfpi) ! rfpi only contained in convoluted quantities 
          integrand(:) = (0.d0,0.d0)
          temp= (0.d0,0.d0)
          integrand(:imt1) = bxc(:imt1,ilm)*mag_den_glob(:imt1,ilm,i)*cone
          call intcheb_cell(integrand(:),temp,rpan_intervall,ipan_intervall,npan_tot,ncheb,irmdnew)
          torque_mt(i) = torque_mt(i) + real(temp*rfpi) ! rfpi only contained in convoluted quantities 
        end do !ilm
      end if
    end do !i
    
    if (t_inc%i_write>1) write(1337,'("iatom, torque, torque_mt = ",i4,6es16.8)') iatom, torque(:), torque_mt(:)
    
    
    ! --------------------------------------------------------------------------------------
    ! calculate the perpendicular component of the torque
    ! MdSD: here theta and phi are the input/fixed angles
    magdir(1) = sin(theta)*cos(phi)
    magdir(2) = sin(theta)*sin(phi)
    magdir(3) = cos(theta)
    torque(:) = torque(:) - magdir*dot_product(magdir(:),torque(:))
    torque_mt(:) = torque_mt(:) - magdir*dot_product(magdir(:),torque_mt(:))
    if (t_inc%i_write>1) write(1337,'("iatom, torque_perp, torque_mt_perp = ",i4,6es16.8)') iatom, torque(:), torque_mt(:)

    ! MdSD: muffin-tin or full cell
    if (bfield%lbfield_mt) then
      bfield%mag_torque(iatom,:) = torque_mt(:)
      torque(:) = torque_mt(:)  ! to use below and avoid another if-statement
    else
      bfield%mag_torque(iatom,:) = torque(:)
    end if

    if(bfield%lbfield_constr .and. t_inc%i_iteration>=bfield%itscf0 .and. t_inc%i_iteration<=bfield%itscf1) then !constraining fields
      if(bfield%ibfield_constr == 0 ) then ! constraining fields based on magnetic torques
        ! sum up the torques for all iterations, which yields a scf with constraining fields
        bfac = 1.d0
        if(fix_direction) then
          bfield%bfield_constr(iatom,:) = bfield%bfield_constr(iatom,:) - torque(:)/totmag*bfac
        end if
        !bfield%bfield_strength_constr         = sqrt(dot_product(bfield%bfield_constr(:),bfield%bfield_constr(:))) 
        !bfield%theta_constr                   = acos(bfield%bfield_constr(3)/bfield%bfield_strength_constr)
        !bfield%phi_constr                     = datan2(bfield%bfield_constr(2),bfield%bfield_constr(1))
        if (t_inc%i_write>1) write(1337,'(" itscf, iatom, ibfield_constr, bfield_constr= ",3i4,100f16.8)') t_inc%i_iteration, iatom ,bfield%ibfield_constr , bfield%bfield_constr(iatom,:)
      end if
      if(bfield%ibfield_constr == 1) then ! constraining fields to constrain scf cycle
        c_old                = bfield%bfield_constr(iatom,:)
        bfac                 = 0.030d0
        ! magdir_old(1)         = cos(t_params%phi(iatom))*sin(t_params%theta(iatom))
        ! magdir_old(2)         = sin(t_params%phi(iatom))*sin(t_params%theta(iatom))
        ! magdir_old(3)         = cos(t_params%theta(iatom))
        magdir_old(1)         = cos(phi)*sin(theta)
        magdir_old(2)         = sin(phi)*sin(theta)
        magdir_old(3)         = cos(theta)
        if (t_inc%i_write>1) write(1337,'(" itscf, iatom, magdir, magdir_old = ",2i4,100f16.8)') t_inc%i_iteration, iatom , magdir_it(:), magdir_old(:)
        if(fix_direction) then
          bfield%bfield_constr(iatom,:) = c_old - dot_product(c_old,magdir_old)*magdir_old - (magdir_it - dot_product(magdir_it,magdir_old)*magdir_old)*bfac
        end if
        if (t_inc%i_write>1) write(1337,'(" itscf, iatom, ibfield_constr, bfield_constr= ",3i4,100f16.8)') t_inc%i_iteration, iatom ,bfield%ibfield_constr , bfield%bfield_constr(iatom,:)
        ! MdSD: constraining fields
        bconstr(1:3) = bfield%bfield_constr(iatom,:)
      end if
      !if(density%magmomentfixed == 6) then ! constraining fields
      !  temp2 = 0.d0
      !  do ilm=1,lmmax
      !    integrand(:) = cone*bxc(:,ilm)*(density%rho2ns_complexnew(:,ilm,2)-density%rho2ns_complexnew(:,ilm,1))
      !    call intcheb_cell(integrand,cellnew,temp)
      !    !call simpk(integrand,temp,cell%npan,cell%nrcut(:),cell%drmeshdi(:),cell%npan)
      !    temp2    = temp2 + Dreal(temp*rfpi) ! this is a spherical quantity =>sqrt(4*pi)
      !  end do !ilm
      !end if
    end if

  end subroutine calc_torque
  


  subroutine calc_magmom_newmesh(lmax,imt1,irmdnew,icell,iend,ifunm,icleb,cleb,thetasnew,mag_den,mag_den_convol)
    
    use :: global_variables, only: ncleb, nfund, lmmaxd, lmpotd
    use :: mod_datatypes, only: dp

    implicit none



    integer                                                       , intent (in) :: lmax     
    integer                                                       , intent (in) :: imt1
    integer                                                       , intent (in) :: irmdnew
    integer              , intent (in)                                          :: icell
    integer              , intent (in)                                          :: iend         ! Number of nonzero gaunt coefficients
    integer              , dimension (1:(2*lmax+1)**2)            , intent (in) :: ifunm        ! pointer array for shapefun     ! Check index and dimensions!!!!!
    integer, dimension (ncleb, 4), intent (in)                                  :: icleb        !! Pointer array
    real (kind=dp), dimension (ncleb), intent (in)                              :: cleb         !! GAUNT coefficients (GAUNT)    ! CHECK THE DIMENSION AND HOW IT IS USED!!!
    real (kind=dp)       , dimension (irmdnew, nfund)   , intent (in)           :: thetasnew    ! shapefun on the Cheby mesh
    real(kind=dp),dimension(1:irmdnew,1:lmpotd,1:3)        , intent (in) :: mag_den
    real(kind=dp),dimension(1:irmdnew,1:(lmax+1)**2 ,1:3)        , intent (out):: mag_den_convol

    real(kind=dp),dimension(1:irmdnew,1:lmpotd)          :: shapefun_mod  
    real(kind=dp),parameter                                     :: rfpi=3.5449077018110318
    integer                                                     :: lmmax     
    integer                                                     :: ifun
    integer                                                     :: i
    integer                                                     :: j
    integer                                                     :: ilm
    integer                                                     :: jlm
    integer                                                     :: lm1
    integer                                                     :: lm2
    integer                                                     :: lm3
    integer                                                     :: ir
    real(kind=dp)                                               :: c0ll
    
    lmmax = (lmax+1)**2
    c0ll                                                    = 1.d0/rfpi

    shapefun_mod(:,:)                                       = 0.d0
    shapefun_mod(1:imt1,1)                                  = rfpi ! is multipled by C_LL^0
    shapefun_mod(imt1+1:irmdnew,1)                          = thetasnew(imt1+1:irmdnew,1)
    do ilm=2,lmpotd
      ifun = ifunm(ilm)
      if(.not. ifun == 0) then !shapefun%lmused(ilm)==1) then
        !write(*,'(" converted ifun= ",i4," to ilm= ",i4)') ifun, ilm
        shapefun_mod(imt1+1:irmdnew,ilm)           = thetasnew(imt1+1:irmdnew,ifun)
      end if
    end do

    mag_den_convol(:,:,:)=0.d0

    !write(*,'("c0ll=",es16.8)') c0ll


    !!diagonal part (not contained in gaunt-coeff)
    do i = 1,3
      do ilm = 1,lmpotd
        mag_den_convol(:,1,i) = mag_den_convol(:,1,i) + mag_den(:,ilm,i)*shapefun_mod(:,ilm)*c0ll
      end do
    end do
    do i = 1,3
      do j = 1,iend !gauntcoeff%iend
        lm1 = icleb(j, 1)! gauntcoeff%icleb(j,1) ! lmax
        lm2 = icleb(j, 2)! gauntcoeff%icleb(j,2) ! lmax
        lm3 = icleb(j, 3)! gauntcoeff%icleb(j,3) ! 2*lmax
        !write(*,'("lm1,lm2,lm3,gaunt=",3i4,2es16.8)') lm1,lm2,lm3,cleb(j)
        if(lm1<= lmmax  .and. lm2 <= lmmax  .and. lm3<= lmmax) then
        ! since lm1 and lm2 are only defined up to lmax the convoluted
        ! magnetization density can also be defined only up to there
          do ir = 1,irmdnew!cellnew%nrmaxnew
            mag_den_convol(ir,lm3,i) = mag_den_convol(ir,lm3,i) + cleb(j)*mag_den(ir,lm1,i)*shapefun_mod(ir,lm2)
            mag_den_convol(ir,lm3,i) = mag_den_convol(ir,lm3,i) + cleb(j)*mag_den(ir,lm2,i)*shapefun_mod(ir,lm1)
          end do
        end if
      end do
    end do


  end subroutine calc_magmom_newmesh

end module mod_torque
