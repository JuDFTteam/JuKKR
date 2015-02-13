module mod_spintools

  implicit none

  private
  public :: spin_expectationvalue, spin_crossterm, Spinrot_AlphaBeta, rotate_wavefunction, &
          & torq_expectationvalue, Spinrot_AlphaBeta_Rashba

contains





  subroutine rotate_wavefunction(lmmaxso, natypd, alpha, beta, rveig_atom_inp, rveig_atom_out)

  implicit none

    ! Arguments
    integer,          intent(in)  :: lmmaxso, natypd
    double precision, intent(in)  :: alpha, beta
    double complex,   intent(in)  :: rveig_atom_inp(lmmaxso, natypd, 2)
    double complex,   intent(out) :: rveig_atom_out(lmmaxso, natypd, 2)

    double complex, parameter :: CI=(0d0,1d0)

    double precision :: alpha2, salph, calph
    double complex   :: efac

    alpha2 = alpha*0.5d0
    salph  = dsin(alpha2)
    calph  = dcos(alpha2)
    efac   = dcos(beta) + CI*sin(beta)

    rveig_atom_out(:,:,1) =  calph * rveig_atom_inp(:,:,1) + salph*efac * rveig_atom_inp(:,:,2)
    rveig_atom_out(:,:,2) = -salph * rveig_atom_inp(:,:,1) + calph*efac * rveig_atom_inp(:,:,2)

  end subroutine rotate_wavefunction





  subroutine spin_expectationvalue(inc, nstates, rhod, rveig_atom, Spin_tot, Spin_atom)
  ! Calculates the spin-expectation value of a state |psi_n> and stores
  !   them into Spin_tot(ixyz)
  !   The corresponding coefficientvector must be already properly normalized.
  !                                put into this subroutine by B.Zimmermann
  !                                                       created: 05.09.12

  use type_inc
  implicit none

    ! Arguments
    type(inc_type), intent(in)  :: inc
    integer,        intent(in)  :: nstates
    double complex, intent(in)  :: rhod(inc%lmmaxso,inc%lmmaxso,inc%natypd,4),&
                                 & rveig_atom(inc%lmmaxso,inc%natypd,nstates)
    double complex, intent(out) :: Spin_tot(3, nstates)
    double complex, intent(out), optional :: Spin_atom(3,inc%natypd,nstates)

    ! Locals
    integer        :: ixyz, iatom, istate
    double complex :: ztmp(inc%natypd), ckhelp(inc%lmmaxso)

    !Parameter
    double complex, parameter :: CONE=(1d0,0d0), CZERO=(0d0, 0d0)

    Spin_tot  = CZERO
    ! calculate Spin expectation value with wavefunctions as stored in rveig_atom
    do istate=1,nstates
      do ixyz=1,3

        ztmp   = CZERO
        do iatom=1,inc%natypd
          ckhelp = CZERO
          call ZGEMM('N', 'N', inc%lmmaxso, 1, inc%lmmaxso, CONE, rhod(:,:,iatom,ixyz+1), inc%lmmaxso, rveig_atom(:,iatom,istate), inc%lmmaxso, CZERO, ckhelp, inc%lmmaxso)
          call ZGEMM('C','N', 1, 1, inc%lmmaxso,CONE, rveig_atom(:,iatom,istate), inc%lmmaxso, ckhelp, inc%lmmaxso, CZERO, ztmp(iatom), 1)
        end do!iatom

        Spin_tot(ixyz,istate) = sum(ztmp)
        if(present(Spin_atom)) Spin_atom(ixyz,:,istate) = ztmp
      end do!ixyz
    end do!istate

  end subroutine spin_expectationvalue




  subroutine torq_expectationvalue(inc, nstates, torq, rveig_atom, Torq_tot, Torq_atom)
  ! Calculates the torque-expectation value of a state |psi_n> and stores
  !   them into Torq_tot(ixyz)
  !   The corresponding coefficientvector must be already properly normalized.
  !                   adapted from spin_expectationvalue by G. Geranton
  !                                                       created: 02.10.14

  use type_inc
  implicit none

    ! Arguments
    type(inc_type), intent(in)  :: inc
    integer,        intent(in)  :: nstates
    double complex, intent(in)  :: torq(inc%lmmaxso,inc%lmmaxso,inc%natypd,3),&
                                 & rveig_atom(inc%lmmaxso,inc%natypd,nstates)
    double complex, intent(out) :: Torq_tot(3,nstates)
    double complex, intent(out), optional :: Torq_atom(3,inc%natypd,nstates)

    ! Locals
    integer        :: ixyz, iatom, istate
    double complex :: ztmp(inc%natypd), ckhelp(inc%lmmaxso)

    !Parameter
    double complex, parameter :: CONE=(1d0,0d0), CZERO=(0d0, 0d0)

    Torq_tot  = CZERO
    ! calculate Torq expectation value with wavefunctions as stored in rveig_atom
    do istate=1,nstates
      do ixyz=1,3

        ztmp   = CZERO
        do iatom=1,inc%natypd
          ckhelp = CZERO
          call ZGEMM('N', 'N', inc%lmmaxso, 1, inc%lmmaxso, CONE, torq(:,:,iatom,ixyz), inc%lmmaxso, rveig_atom(:,iatom,istate), inc%lmmaxso, CZERO, ckhelp, inc%lmmaxso)
          call ZGEMM('C','N', 1, 1, inc%lmmaxso,CONE, rveig_atom(:,iatom,istate), inc%lmmaxso, ckhelp, inc%lmmaxso, CZERO, ztmp(iatom), 1)
        end do!iatom

        Torq_tot(ixyz,istate) = sum(ztmp)
        if(present(Torq_atom)) Torq_atom(ixyz,:,istate) = ztmp
      end do!ixyz
    end do!istate


  end subroutine torq_expectationvalue




  subroutine Spinrot_AlphaBeta(lincombtype, nvec, Spin_ini, Scross, alpha, beta, Spin_estimated, uio)

    use mod_mathtools, only: pi
    implicit none

    ! Arguments
    integer,          intent(in)  :: lincombtype
    double precision, intent(in)  :: nvec(3)
    double complex,   intent(in)  :: Spin_ini(3), Scross(3)
    double precision, intent(out) :: alpha, beta, Spin_estimated
    integer, intent(in), optional :: uio

    ! Locals
    integer          :: icomb, combtake
    double precision :: beta_tmp(4), alphacal_tmp(4), alpha_tmp(4), gvec1(3), gvec2(3), maxv, absv
    double complex   :: Spin_unrot_n, S_12_n, S_unrot_g1, S_unrot_g2, S_12_g1, S_12_g2, xi_help, Spin_rot_test(4)

    alpha= 0d0
    beta = 0d0
    beta_tmp=0d0
    alphacal_tmp=0d0
    alpha_tmp=0d0

    ! project spin onto chosen direction
    Spin_unrot_n = sum(nvec(:)*Spin_ini(:))
    S_12_n       = sum(nvec(:)*Scross(:))

    ! chose condition for linear combination of |psi1> and |psi2>
    if(lincombtype==1) then!MAXIMIZE spin in direction of SQA

      ! calculate combinations between (alpha, alpha+pi) and (beta, beta+pi)
      beta_tmp(1) = -datan2( dimag(S_12_n),  dble(S_12_n) )
      beta_tmp(2) = beta_tmp(1)
      beta_tmp(3) = beta_tmp(1) + pi
      beta_tmp(4) = beta_tmp(1) + pi

      alphacal_tmp = dcos(beta_tmp)*dble(S_12_n) - dsin(beta_tmp)*dimag(S_12_n)

      alpha_tmp(1) = datan2(alphacal_tmp(1), dble(Spin_unrot_n))
      alpha_tmp(2) = alpha_tmp(1) + pi
      alpha_tmp(3) = datan2(alphacal_tmp(2), dble(Spin_unrot_n))
      alpha_tmp(4) = alpha_tmp(3) + pi

    elseif(lincombtype==2) then!VANISHING spin in directions perpendicular to SQA

      !construct perpendicular directions to nvec
      gvec1 = 0d0
      gvec2 = 0d0
      if(abs(nvec(1))>1d-04 .or. abs(nvec(3))>1d-04) then
        gvec1(1) =  nvec(3)
        gvec1(2) =  0d0
        gvec1(3) = -nvec(1)
      else
        gvec1(1) =  0d0
        gvec1(2) = -nvec(3)
        gvec1(3) =  nvec(2)
      end if
      gvec1 = gvec1 / dsqrt( gvec1(1)**2 + gvec1(2)**2 + gvec1(3)**2 )
      gvec2(1) = nvec(2)*gvec1(3) - nvec(3)*gvec1(2)
      gvec2(2) = nvec(3)*gvec1(1) - nvec(1)*gvec1(3)
      gvec2(3) = nvec(1)*gvec1(2) - nvec(2)*gvec1(1)
      !write(*,*) 'check norm gvec1:', dsqrt( gvec1(1)**2 + gvec1(2)**2 + gvec1(3)**2 )
      !write(*,*) 'check norm gvec2:', dsqrt( gvec2(1)**2 + gvec2(2)**2 + gvec2(3)**2 )

      S_12_g1 = gvec1(1)*Scross(1) + gvec1(2)*Scross(2) + gvec1(3)*Scross(3)
      S_12_g2 = gvec2(1)*Scross(1) + gvec2(2)*Scross(2) + gvec2(3)*Scross(3)

      S_unrot_g1 = sum(gvec1(:)*Spin_ini(:))
      S_unrot_g2 = sum(gvec2(:)*Spin_ini(:))

      xi_help = -(S_unrot_g2*conjg(S_12_g1) - S_unrot_g1*conjg(S_12_g2))/(S_unrot_g2*S_12_g1 - S_unrot_g1*S_12_g2)

      ! calculate combinations between (alpha, alpha+pi) and (beta, beta+pi)
      beta_tmp(1) = 0.5d0*datan2( dimag(xi_help),  dble(xi_help) )
      beta_tmp(2) = beta_tmp(1)
      beta_tmp(3) = beta_tmp(1) + pi
      beta_tmp(4) = beta_tmp(1) + pi

      alphacal_tmp = -dcos(beta_tmp)*dble(S_12_g1) + dsin(beta_tmp)*dimag(S_12_g1)
      !write(*,*) 'imag(S_unrot_g1)=0? : ', dimag(S_unrot_g1)

      alpha_tmp(1) = datan2(dble(S_unrot_g1), alphacal_tmp(1))
      alpha_tmp(2) = alpha_tmp(1) + pi
      alpha_tmp(3) = datan2(dble(S_unrot_g1), alphacal_tmp(2))
      alpha_tmp(4) = alpha_tmp(3) + pi

    else
       stop 'no rule to form linear combination'
    end if 

    !Calculate the expected spin-expectation with (alpha_tmp, beta_tmp)
    Spin_rot_test = dcos(alpha_tmp)*Spin_unrot_n + dsin(alpha_tmp)*(dcos(beta_tmp)*dble(S_12_n) - dsin(beta_tmp)*dimag(S_12_n))

    ! pick (alpha, beta) maximizing the absolut value of the spin in direction n
    maxv = 0d0
    do icomb=1,4
      absv = dble(Spin_rot_test(icomb))
      if(absv>maxv) then
        maxv  = absv
        combtake=icomb
      end if!absv>maxv
    end do!icomb

    alpha = alpha_tmp(combtake)
    beta  = beta_tmp(combtake)
    Spin_estimated = Spin_rot_test(combtake)

    if(present(uio)) then
      write(uio,"(6X,A)") 'Combinations of alpha/beta yield:'
      do icomb=1,4
        write(uio,"(8X,A,I0,A,2ES25.16,A)") 'icomb=', icomb, ', Spin_rot_test=(', Spin_rot_test(icomb), ')'
      end do
      write(uio,"(8X,A,I0)") 'Take combination #', combtake
    end if!present(uio)

  end subroutine Spinrot_AlphaBeta


  subroutine Spinrot_AlphaBeta_Rashba(Spin_ini_selected, Scross_selected, alpha, beta)
    
    use mod_mathtools, only: pi
    implicit none

    ! Arguments
    double complex,   intent(in)  :: Spin_ini_selected(3,2), Scross_selected(3)
    double precision, intent(out) :: alpha, beta

    ! Locals
    integer          :: nalpha, nbeta, niter, ialpha, ibeta, iiter, k_candidate, ixyz, k_take, i_kcand, ierr
    double precision :: alpha_tmp, beta_tmp, para_a1, para_b1, para_c1, para_a2, para_b2, para_c2, para_d2, para_e2, para_f2, delta, re_eibS12, im_eibS12, alpha_bounds(2), alpha0, sinalpha, cosalpha, d_alpha, ds_dalpha0, maxspin, maxspin_k, sumS(2), s_ab, s_ab2
    double precision, allocatable :: alpha_get(:), beta_get(:),ds_dalpha(:)
    !Parameter
    double complex, parameter :: CONE=(1d0, 0d0), CZERO=(0d0, 0d0), CI=(0d0, 1d0)



    ! hard coded parameters for loops to find alpha and beta
    nalpha = 1000
    nbeta  = 1000
    niter  = 1

    allocate(ds_dalpha(0:nalpha+1), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating ds_dalpha'
    allocate(alpha_get(0:((nalpha+1)*(nbeta+1))), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating alpha_get'
    allocate(beta_get(0:((nalpha+1)*(nbeta+1))), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating beta_get'
    alpha_get = 0d0 
    beta_get  = 0d0
    k_candidate = 0
    !write(*,*) 'start beta loop'

    do ibeta=1,nbeta+1
       beta_tmp  = 2d0*pi*(ibeta-1)/nbeta
       re_eibS12 = 0d0
       im_eibS12 = 0d0
       para_a1   = 0d0
       para_b1   = 0d0
       para_c1   = 0d0
       para_a2   = 0d0
       para_b2   = 0d0
       para_c2   = 0d0
       para_d2   = 0d0
       para_e2   = 0d0
       para_f2   = 0d0
       do ixyz=1,3
          re_eibS12 = real(exp(CI*beta_tmp)*Scross_selected(ixyz))
          im_eibS12 = aimag(exp(CI*beta_tmp)*Scross_selected(ixyz))
          para_a1   = para_a1+real(Spin_ini_selected(ixyz,2))*im_eibS12
          para_b1   = para_b1+re_eibS12*im_eibS12
          para_c1   = para_c1+real(Spin_ini_selected(ixyz,1))*im_eibS12
          para_a2   = para_a2+real(Spin_ini_selected(ixyz,2))*re_eibS12
          para_b2   = para_b2+real(Spin_ini_selected(ixyz,2))**2
          para_c2   = para_c2+re_eibS12**2
          para_d2   = para_d2+real(Spin_ini_selected(ixyz,1))*real(Spin_ini_selected(ixyz,2))
          para_e2   = para_e2+real(Spin_ini_selected(ixyz,1))**2
          para_f2   = para_f2+real(Spin_ini_selected(ixyz,1))*re_eibS12
       enddo !ixyz=1,3
       ! Now paramters are set, check if alpha can be found
       delta = para_b1**2-para_a1*para_c1
       if (delta.ge.0d0)then
          alpha_tmp = 0d0
          alpha_bounds(1) = 0
          alpha_bounds(2) = 2*pi
          do iiter=1,niter
             ds_dalpha = czero
                !write(*,*) 'start alpha loop',ibeta
             do ialpha=1,nalpha+1
                d_alpha = (alpha_bounds(2)-alpha_bounds(1))/nalpha
                alpha_tmp = alpha_bounds(1)+d_alpha*(ialpha-1)
                sinalpha = sin(alpha_tmp/2d0)
                cosalpha = cos(alpha_tmp/2d0)
                ds_dalpha(ialpha) = -2d0*cosalpha**3*sinalpha*para_e2+2d0*sinalpha**3*cosalpha*para_b2+2d0*cosalpha**3*sinalpha*para_d2-2d0*sinalpha**3*cosalpha*para_d2+4d0*cosalpha**3*sinalpha*para_c2-4d0*sinalpha**3*cosalpha*para_c2+2d0*cosalpha**4*para_f2-6d0*cosalpha**2*sinalpha**2*para_f2+6d0*sinalpha**2*cosalpha**2*para_a2-2d0*sinalpha**4*para_a2
                ! check if ds_dalpha has sign change
                if ((ds_dalpha(ialpha)*ds_dalpha(ialpha-1)).lt.0d0)then
                   alpha0 = alpha_tmp-d_alpha-(ds_dalpha(ialpha-1))/(ds_dalpha(ialpha)-ds_dalpha(ialpha-1))*(d_alpha)
                   sinalpha = sin(alpha0/2d0)
                   cosalpha = cos(alpha0/2d0)
                   ds_dalpha0 =-2d0*cosalpha**3*sinalpha*para_e2+2d0*sinalpha**3*cosalpha*para_b2+2d0*cosalpha**3*sinalpha*para_d2-2d0*sinalpha**3*cosalpha*para_d2+4d0*cosalpha**3*sinalpha*para_c2-4d0*sinalpha**3*cosalpha*para_c2+2d0*cosalpha**4*para_f2-6d0*cosalpha**2*sinalpha**2*para_f2+6d0*sinalpha**2*cosalpha**2*para_a2-2d0*sinalpha**4*para_a2
                   if (ds_dalpha0.le.0d0) then
!                      alphabounds(1) = alpha0
!                      alphabounds(2) = alpha_tmp
!                   else !ds_dalpha0.le.0d0
!                      alphabounds(1) = alpha_tmp-dalpha
!                      alphabounds(2) = alpha0
                      k_candidate = k_candidate+1
                      alpha_get(k_candidate) = alpha0
                      beta_get(k_candidate) = beta_tmp
                   endif !ds_dalpha0.le.0d0
                endif !(ds_dalpha(ialpha)*ds_dalpha(ialpha-1)).lt.0d0
             enddo !ialpha=1,nalpha
                !write(*,*) 'finish alpha loop',ibeta
          enddo !iiter=1,niter
       endif !delta.ge.0d0
    enddo !ibeta=1,nbeta

    ! Now check which of the alpha/beta candidates maximize the spin
    k_take  = 0
    maxspin = 0d0
    if (k_candidate.ge.1)then
       do i_kcand=1,k_candidate
          maxspin_k = 0d0
          do ixyz=1,3
             s_ab = cos(alpha_get(i_kcand)/2d0)**2*real(Spin_ini_selected(ixyz,1))+sin(alpha_get(i_kcand)/2d0)**2*real(Spin_ini_selected(ixyz,2))+2*sin(alpha_get(i_kcand)/2d0)*cos(alpha_get(i_kcand)/2d0)*real(exp(CI*beta_get(i_kcand))*Scross_selected(ixyz))
             maxspin_k = maxspin_k + s_ab**2
          enddo !ixyz=1,3
          if (maxspin_k.gt.maxspin)then
             maxspin = maxspin_k
             k_take = i_kcand
          endif !
                !write(*,*) 'in 2nd loop',i_kcand
       enddo !i_kcand=1,k_candidate
    endif !k_candidate.ge.1
     

    ! Now k_take is the index which maximizes the total spin
    alpha = alpha_get(k_take)
    beta  = beta_get(k_take)

    sumS = 0d0 
    do ixyz=1,3
       s_ab = cos(alpha/2d0)**2*real(Spin_ini_selected(ixyz,1))+sin(alpha/2d0)**2*real(Spin_ini_selected(ixyz,2))+2*sin(alpha/2d0)*cos(alpha/2d0)*real(exp(CI*beta)*Scross_selected(ixyz))
       s_ab2 = cos(alpha/2d0)**2*real(Spin_ini_selected(ixyz,2))+sin(alpha/2d0)**2*real(Spin_ini_selected(ixyz,1))+2*sin(alpha/2d0)*cos(alpha/2d0)*real(exp(CI*beta)*conjg(Scross_selected(ixyz)))
       sumS(1) = sumS(1) + s_ab**2
       sumS(2) = sumS(2) + s_ab2**2
    enddo !ixyz=1,3
    if (sumS(1).lt.sumS(2))then
       alpha = alpha+pi
       beta = beta+pi
    endif

  end subroutine Spinrot_AlphaBeta_Rashba



   subroutine spin_crossterm(inc, rhod, rveig_atom, Spi12, Spi12_atom)
  ! Calculates the spin-crossterm 2 (degenerate) states |psi_i>.
  !   (see PhD thesis of Swantje Heers, chapter 4.6)
  !   The corresponding coefficientvector must be already properly normalized.
  !                                put into this subroutine by B.Zimmermann
  !                                                       created: 05.09.12

  use type_inc
  implicit none

    ! Arguments
    type(inc_type), intent(in)  :: inc
    double complex, intent(in)  :: rhod(inc%lmmaxso,inc%lmmaxso,inc%natypd,4), rveig_atom(inc%lmmaxso,inc%natypd,2)
    double complex, intent(out) :: Spi12(3)
    double complex, intent(out), optional :: Spi12_atom(3,inc%natypd)

    ! Locals
    integer        :: ixyz, iatom
    double complex :: ztmp, ckhelp(inc%lmmaxso), Spi12_tmp(3,inc%natypd)

    !Parameter
    double complex, parameter :: CONE=(1d0,0d0), CZERO=(0d0, 0d0)

    Spi12 = CZERO
    Spi12_atom = CZERO

    ! calculate S^i_{c_1,c_2} = c_1^dagger  \Sigma^i c_2, i=x,y,z     Spi12 = (0d0, 0d0)
    do ixyz = 1,3
      do iatom=1,inc%natypd
        ckhelp = CZERO
        ztmp   = CZERO

        call ZGEMM('N', 'N', inc%lmmaxso, 1, inc%lmmaxso, CONE, rhod(:,:,iatom,ixyz+1), inc%lmmaxso, rveig_atom(:,iatom,2), inc%lmmaxso, CZERO, ckhelp, inc%lmmaxso)
        call ZGEMM('C','N', 1, 1, inc%lmmaxso, CONE, rveig_atom(:,iatom,1), inc%lmmaxso, ckhelp, inc%lmmaxso, CZERO, ztmp, 1)

        Spi12_tmp(ixyz,iatom) = ztmp
      end do
    end do

   do ixyz = 1,3
     Spi12(ixyz) = sum(Spi12_tmp(ixyz,:))
   end do

   if(present(Spi12_atom)) Spi12_atom = Spi12_tmp

  end subroutine spin_crossterm

end module mod_spintools
