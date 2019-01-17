module mod_rotatespin

contains

!-------------------------------------------------------------------------------
!> Summary: Rotate integrated spin density to back to global frame and construct vector of magnetic moment
!> Author: 
!> Category: KKRimp, physical-observables
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!>
!> @note Computes both total and scattering moments (used if noscattergmoment is not given as testflag) @endnote
!-------------------------------------------------------------------------------
  subroutine rotatespin(density,cell,lmax,config,ispin)
    use nrtype, only: pi
    use type_density, only: density_type
    use type_cell, only: cell_type
    use mod_cheb, only: intcheb_complex
    use mod_rotatespinframe, only: rotatematrix, spinmode
    use mod_config, only: config_testflag
    use mod_mathtools, only: rotvector
    use type_config, only: config_type
    implicit none
    
    type(density_type)    :: density
    type(cell_type)       :: cell
    integer               :: lmax,lmpot
    type(config_type)     :: config
    integer               :: ispin
    double complex        :: rho2ns_temp(2,2)
    double precision      :: magmoment(3),magmoment2(3),totmagmoment,totxymagmoment,totmagmoment2,totxymagmoment2
    double precision      :: theta_old,phi_old

    lmpot=(2*lmax+1)**2
    
    theta_old        = density%theta
    phi_old          = density%phi
    density%thetaold = density%theta
    density%phiold   = density%phi
    density%theta2   = density%theta
    density%phi2     = density%phi
    
    
    write(*,*) 'Rotation of the charge density'
    write(*,*) 'The integrated charge density is:'
    ! write(*,*) 'rho_up_up ',density%rho2ns_integrated(1)
    ! write(*,*) 'rho_dn_dn ',density%rho2ns_integrated(2)
    ! write(*,*) 'rho_up_dn ',density%rho2ns_integrated(3)
    ! write(*,*) 'rho_dn_up ',density%rho2ns_integrated(4)
    
    
    rho2ns_temp(1,1)=density%rho2ns_integrated(1)
    rho2ns_temp(2,2)=density%rho2ns_integrated(2)
    rho2ns_temp(1,2)=density%rho2ns_integrated(3)
    rho2ns_temp(2,1)=density%rho2ns_integrated(4)
    
    call rotatematrix(rho2ns_temp,theta_old,phi_old,1, 0) ! 'loc->glob')
    
    density%rho2ns_integrated(1)=rho2ns_temp(1,1)
    density%rho2ns_integrated(2)=rho2ns_temp(2,2)
    density%rho2ns_integrated(3)=rho2ns_temp(1,2)
    density%rho2ns_integrated(4)=rho2ns_temp(2,1)
    
    
    rho2ns_temp(1,1)=density%rho2ns_integrated_scattering(1)
    rho2ns_temp(2,2)=density%rho2ns_integrated_scattering(2)
    rho2ns_temp(1,2)=density%rho2ns_integrated_scattering(3)
    rho2ns_temp(2,1)=density%rho2ns_integrated_scattering(4)
    
    call rotatematrix(rho2ns_temp,theta_old,phi_old,1, 0) !'loc->glob')
    
    density%rho2ns_integrated_scattering(1)=rho2ns_temp(1,1)
    density%rho2ns_integrated_scattering(2)=rho2ns_temp(2,2)
    density%rho2ns_integrated_scattering(3)=rho2ns_temp(1,2)
    density%rho2ns_integrated_scattering(4)=rho2ns_temp(2,1)

    
    density%magmomentold=density%magmoment
    density%magmomentold2=density%magmoment2
    
    
    if (spinmode=='regular') then
      magmoment(1)=-1*DImag( density%rho2ns_integrated(3)+density%rho2ns_integrated(4) )
      magmoment(2)=-1*Dreal(-density%rho2ns_integrated(3)+density%rho2ns_integrated(4) )
      magmoment(3)=-1*DImag( density%rho2ns_integrated(1)-density%rho2ns_integrated(2) )
      magmoment2(1)=-1*DImag( density%rho2ns_integrated_scattering(3)+density%rho2ns_integrated_scattering(4) )
      magmoment2(2)=-1*Dreal(-density%rho2ns_integrated_scattering(3)+density%rho2ns_integrated_scattering(4) )
      magmoment2(3)=-1*DImag( density%rho2ns_integrated_scattering(1)-density%rho2ns_integrated_scattering(2) )
    else if (spinmode=='kkr') then
      magmoment(1)= 1*DImag( density%rho2ns_integrated(3)+density%rho2ns_integrated(4) )
      magmoment(2)=-1*Dreal( density%rho2ns_integrated(3)-density%rho2ns_integrated(4) )
      magmoment(3)= 1*DImag(-density%rho2ns_integrated(1)+density%rho2ns_integrated(2) )
      magmoment2(1)= 1*DImag( density%rho2ns_integrated_scattering(3)+density%rho2ns_integrated_scattering(4) )
      magmoment2(2)=-1*Dreal( density%rho2ns_integrated_scattering(3)-density%rho2ns_integrated_scattering(4) )
      magmoment2(3)= 1*DImag(-density%rho2ns_integrated_scattering(1)+density%rho2ns_integrated_scattering(2) )
    else
      stop '[rotatespin] error in spinmode'
    end if
    
    write(*,*) 'magnetic moment is ',magmoment
    write(1337,*) 'magnetic moment is ',magmoment
    
    write(*,*) 'magnetic moment2 is ',magmoment2
    write(1337,*) 'magnetic moment2 is ',magmoment2
    density%magmoment  = magmoment
    density%magmoment2 = magmoment2

    
    if (config%spinmixfac/=1.0D0 .and. ispin<=config%spinmixbound) then
      density%magmoment =rotvector(theta_old,phi_old,magmoment ,config%spinmixfac)
      density%magmoment2=rotvector(theta_old,phi_old,magmoment2,config%spinmixfac)
      write(*,*) 'magnetic moment with spinmix is ',density%magmoment
      write(1337,*) 'magnetic moment with spinmix is ',density%magmoment
      write(*,*) 'magnetic moment2 with spinmix is ',density%magmoment2
      write(1337,*) 'magnetic moment2 with spinmix is ',density%magmoment2
      magmoment=density%magmoment
      magmoment2=density%magmoment2
    end if
    
    totmagmoment=SQRT(magmoment(1)**2+magmoment(2)**2+magmoment(3)**2)
    totxymagmoment=SQRT(magmoment(1)**2+magmoment(2)**2)
    
    totmagmoment2=SQRT(magmoment2(1)**2+magmoment2(2)**2+magmoment2(3)**2)
    totxymagmoment2=SQRT(magmoment2(1)**2+magmoment2(2)**2)
    

    !if (.false.) then
    !
    !  ! #######################################
    !  ! VERSION UNTIL 2012/04/2012
    !  ! ######################################
    !  if(abs(totxymagmoment) > 1.d-5 .and. density%magmomentfixed/=1)then
    !  
    !    if(abs(magmoment(3)) < 1.d-5)then
    !      density%theta=pi/2.0D0
    !    else
    !      density%theta=acos(magmoment(3)/totmagmoment)
    !    endif 
    !    write(*,*)   'new theta1 [deg]',density%theta*180/pi
    !    write(1337,*)'new theta1 [deg]',density%theta*180/pi
    !    
    !    if(totxymagmoment < 1.d-5)then
    !      density%phi=0.0D0 !phiold
    !      write(*,*) 'implement method here'
    !      write(*,*) 'stop'
    !      stop
    !    else   
    !      density%phi=datan2(magmoment(2),magmoment(1))
    !    endif
    !    write(*,*)   'new phi1 [deg]',density%phi*180/pi
    !    write(1337,*)'new phi1 [deg]',density%phi*180/pi
    !    write(1337,*) 'the charge and spin density is been rotated'
    !  
    !    write(*,*) 'the charge and spin density is been rotated'
    !  !   call rotrho2ns(density%rho2ns_complex,density%rho2ns,cellnew%nrmaxnew,lmpot,4,theta,phi)
    !  !     call rotatevector(density%rho2ns_complex,density%rho2ns,cell%nrmax,lmpot,density%theta,density%phi,theta_old,phi_old)
    !  
    !    write(23452324,'(5000F)') density%theta*180/pi,density%phi*180/pi
    !  
    !  else
    !    if (density%magmomentfixed==1) then
    !      write(*,*) 'no rotation magnetic moment is fixed by config'
    !    else
    !      write(*,*) 'no rotation needed'
    !    end if
    !    density%rho2ns(:,:,1)= dimag ( density%rho2ns_complex(:,:,1) )
    !    density%rho2ns(:,:,2)= dimag ( density%rho2ns_complex(:,:,2) )
    !    write(23452324,'(5000F)') theta_old*180/pi,phi_old*180/pi
    !  
    !  end if
    !  ! #######################################
    !  ! ######################################
    !
    !else
    
    ! #######################################
    ! TEST VERSION
    ! ######################################
    if (density%magmomentfixed/=1) then
      write(1337,*) 'Version that rotates all moments'
        density%theta=acos(magmoment(3)/totmagmoment)
        density%theta2=acos(magmoment2(3)/totmagmoment2)
      write(*,*)   'new theta1 [deg]',density%theta*180/pi
      write(1337,*)'new theta1 [deg]',density%theta*180/pi
      write(*,*)   'new theta2 [deg]',density%theta2*180/pi
      write(1337,*)'new theta2 [deg]',density%theta2*180/pi
    
        density%phi=datan2(magmoment(2),magmoment(1))
        density%phi2=datan2(magmoment2(2),magmoment2(1))
      write(*,*)   'new phi1 [deg]',density%phi*180/pi
      write(1337,*)'new phi1 [deg]',density%phi*180/pi
      write(*,*)   'new phi2 [deg]',density%phi2*180/pi
      write(1337,*)'new phi2 [deg]',density%phi2*180/pi
      write(1337,*) 'the charge and spin density is been rotated'
    
      write(*,*) 'the charge and spin density is been rotated'
    else
      write(*,*) 'no rotation magnetic moment is fixed by config'
      density%rho2ns(:,:,1)= dimag ( density%rho2ns_complex(:,:,1) )
      density%rho2ns(:,:,2)= dimag ( density%rho2ns_complex(:,:,2) )
    ! density%dorotate=1
    end if
    ! call rotatevector(density%rho2ns_complex,density%rho2ns,cell%nrmax,lmpot,density%theta,density%phi,theta_old,phi_old)

    !end if
    
    if (.not. config_testflag('noscatteringmoment')) then
        density%phi=density%phi2
        density%theta=density%theta2
        density%magmoment2=density%magmoment
    end if
    
      write(23452324,'(5000F)') density%theta*180/pi,density%phi*180/pi
      write(23452325,'(5000F)') magmoment

  end subroutine rotatespin

end module mod_rotatespin
