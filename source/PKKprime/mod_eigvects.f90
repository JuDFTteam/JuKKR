!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_eigvects

implicit none

  private
  public :: orthgonalize_wavefunc, normalize_wavefunc, calc_norm_wavefunc, calc_proj_wavefunc,&
          & normeigv_new, compare_eigv, rewrite_eigv_atom,compare_eigv_memopt,normeigv_new_memopt



contains

  subroutine rewrite_eigv_atom(inc,nstates,rveig,rveig_atom)

    use type_inc
    implicit none

    type(inc_type), intent(in)  :: inc
    integer,        intent(in)  :: nstates
    double complex, intent(in)  :: rveig(inc%almso, nstates)
    double complex, intent(out) :: rveig_atom(inc%lmmaxso, inc%natypd, nstates)

    integer :: istate, ispin, iatom, lm1, lm1p, lm1s

    rveig_atom = (0d0,0d0)
    do istate=1,nstates
      do ispin=1,2
        do iatom=1,inc%natypd
          do lm1=1,inc%lmmax
            lm1p = (ispin-1)*inc%alm   + (iatom-1)*inc%lmmax + lm1
            lm1s = (ispin-1)*inc%lmmax + lm1
            rveig_atom(lm1s,iatom,istate) = rveig(lm1p,istate)
          end do!lm1
        end do!iatom
      end do!ispin
    end do!ient 

  end subroutine rewrite_eigv_atom



  subroutine compare_eigv(inc, LVref, RV2, LMout, uio)
    ! Compare an reference-eigenvector (LVref) with several other
    !   eigenvectors (contained in RV2) and find the one with the
    !   largest overlap (index LMout). However, only consider
    !   those eigenvectors LV2(:,i2) when take(i2)=.true..


    use type_inc
    use mod_mathtools, only: bubblesort

    implicit none

    type(inc_type), intent(in) :: inc
    double complex, intent(in) :: LVref(inc%almso), &
                                & RV2(inc%almso,inc%almso)
    integer, intent(out)       :: LMout
    integer, optional, intent(in) :: uio

    integer          :: lm2, sorted(inc%almso)
    double complex   :: cprojs(inc%almso)
    double precision :: projs(inc%almso), projmax

    !calculate projections on other eigenvectors
    cprojs = (0d0, 0d0)
    projs  = 0d0
    do lm2=1,inc%almso
        cprojs(lm2) = dot_product(LVref,RV2(:,lm2))
        projs(lm2)  = dble(cprojs(lm2))**2 + dimag(cprojs(lm2))**2
    end do!lm2


    !find corresponding eigenvector
    projmax = 0d0
    LMout =  0
    do lm2=1,inc%almso
      if(projs(lm2)>projmax) then
        projmax = projs(lm2)
        LMout   = lm2
      end if!proj==projmax
    end do!lm2


    if(present(uio)) then
      call bubblesort(inc%almso, projs, sorted)
      write(uio,"(5X,6(I8,16X))") sorted(inc%almso-5:inc%almso)
      write(uio,"(5X,6(8X,E16.6))") projs(sorted(inc%almso-5:inc%almso))
    end if!present(uio)


  end subroutine compare_eigv





  subroutine compare_eigv_memopt(inc, LVref, RV2, nb_ev, LMout, uio)
    ! Compare an reference-eigenvector (LVref) with several other
    !   eigenvectors (contained in RV2) and find the one with the
    !   largest overlap (index LMout). However, only consider
    !   those eigenvectors LV2(:,i2) when take(i2)=.true..


    use type_inc
    use mod_mathtools, only: bubblesort

    implicit none

    type(inc_type), intent(in) :: inc
    double complex, intent(in) :: LVref(inc%almso), &
                                & RV2(inc%almso,inc%neig)
    integer, intent(in)        :: nb_ev
    integer, intent(out)       :: LMout
    integer, optional, intent(in) :: uio

    integer          :: lm2, sorted(inc%neig)
    double complex   :: cprojs(inc%neig)
    double precision :: projs(inc%neig), projmax

    !calculate projections on other eigenvectors
    cprojs = (0d0, 0d0)
    projs  = 0d0
    do lm2=1,nb_ev
        cprojs(lm2) = dot_product(LVref,RV2(:,lm2))
        projs(lm2)  = dble(cprojs(lm2))**2 + dimag(cprojs(lm2))**2
    end do!lm2


    !find corresponding eigenvector
    projmax = 0d0
    LMout =  0
    do lm2=1,nb_ev
      if(projs(lm2)>projmax) then
        projmax = projs(lm2)
        LMout   = lm2
      end if!proj==projmax
    end do!lm2

    if(present(uio)) then
      call bubblesort(inc%neig, projs, sorted)
      write(uio,"(5X,6(I8,16X))") sorted(inc%neig-5:inc%neig)
      write(uio,"(5X,6(8X,E16.6))") projs(sorted(inc%neig-5:inc%neig))
    end if!present(uio)


  end subroutine compare_eigv_memopt





  subroutine orthgonalize_wavefunc(inc, rhod_norm, rveig_atom)
  ! Modifies two coefficient-vectors rveig_atom(:,:,1) and rveig(:,:,2) such that
  !   the corresponding wavefunctions |psi> = \sum_{LL'} R_{LL'} * Y_L * c_L'
  !   are orthogonal to each other.
  !                                put into this subroutine by B.Zimmermann
  !                                                       created: 05.09.12

  use type_inc
  implicit none

    !Arguments
    type(inc_type), intent(in)    :: inc
    double complex, intent(in)    :: rhod_norm(inc%lmmaxso, inc%lmmaxso, inc%natypd)
    double complex, intent(inout) :: rveig_atom(inc%lmmaxso, inc%natypd, 2)

    !Locals
    double complex :: norm(1), proj


      norm = (0d0, 0d0)
      proj = (0d0, 0d0)

      call calc_norm_wavefunc(inc, 1, rhod_norm, rveig_atom(:,:,1), norm)
      call calc_proj_wavefunc(inc, rhod_norm, rveig_atom, proj)

      !orthogonalize
      rveig_atom(:,:,2) = rveig_atom(:,:,2) - proj/norm(1) * rveig_atom(:,:,1)
       
  end subroutine orthgonalize_wavefunc





  subroutine normalize_wavefunc(inc, nstates, rhod_norm, rveig_atom)
  ! Modifies n coefficient-vectors rveig_atom(:,:,i) (i=1..n) such that
  !   the corresponding wavefunctions |psi> = \sum_{LL'} R_{LL'} * Y_L * c_L'
  !   are normalized to 1, i.e. <psi|psi> = <c|rhod_norm|c> 1
  !                                put into this subroutine by B.Zimmermann
  !                                                       created: 05.09.12

  use type_inc
  implicit none

    !Arguments
    type(inc_type), intent(in)    :: inc
    integer,        intent(in)    :: nstates
    double complex, intent(in)    :: rhod_norm(inc%lmmaxso, inc%lmmaxso, inc%natypd)
    double complex, intent(inout) :: rveig_atom(inc%lmmaxso, inc%natypd, nstates)

    !Locals
    integer :: iatom, istate
    double complex :: norm(nstates)

    ! calculate the norm
    call calc_norm_wavefunc(inc, nstates, rhod_norm, rveig_atom, norm)

    ! normalize wavefunctions
    do istate=1,nstates
      rveig_atom(:,:,istate) = rveig_atom(:,:,istate) / CDSQRT(norm(istate))
    end do!ient

  end subroutine normalize_wavefunc





  subroutine calc_norm_wavefunc(inc, nstates, rhod_norm, rveig_atom, norm, norm_atom)
  ! Calculate the norm of n wavefunctions, optionally also atomwise contribution.
  use type_inc
  implicit none

    !Arguments
    type(inc_type), intent(in)  :: inc
    integer,        intent(in)  :: nstates
    double complex, intent(in)  :: rhod_norm(inc%lmmaxso, inc%lmmaxso, inc%natypd)
    double complex, intent(in)  :: rveig_atom(inc%lmmaxso, inc%natypd, nstates)
    double complex, intent(out) :: norm(nstates)
    double complex, intent(out), optional :: norm_atom(inc%natypd,nstates)

    !Locals
    logical :: atomwise
    integer :: iatom, istate
    double complex :: ztmp, ckhelp(inc%lmmaxso)

    !Parameter
    double complex, parameter :: CONE=(1d0,0d0), CZERO=(0d0, 0d0)

    atomwise = present(norm_atom)

    ! calculate the norm
    norm = CZERO
    if(atomwise) norm_atom = CZERO
    do istate=1,nstates
      do iatom=1,inc%natypd
        ckhelp = CZERO
        ztmp   = CZERO
        call ZGEMM('N', 'N', inc%lmmaxso, 1, inc%lmmaxso, CONE, rhod_norm(:,:,iatom), inc%lmmaxso, rveig_atom(:,iatom,istate), inc%lmmaxso, CZERO, ckhelp, inc%lmmaxso)
        call ZGEMM('C','N', 1, 1, inc%lmmaxso,CONE, rveig_atom(:,iatom,istate), inc%lmmaxso, ckhelp, inc%lmmaxso, CZERO, ztmp, 1)
        norm(istate) = norm(istate) + ztmp
        if(atomwise) norm_atom(iatom,istate)= ztmp
      end do!iatom
    end do!istate

  end subroutine calc_norm_wavefunc





  subroutine calc_proj_wavefunc(inc, rhod_norm, rveig_atom, proj)
  ! calculate projection from wavefunction 2 onto wavefunction 1
  use type_inc
  implicit none

    !Arguments
    type(inc_type), intent(in)  :: inc
    double complex, intent(in)  :: rhod_norm(inc%lmmaxso, inc%lmmaxso, inc%natypd)
    double complex, intent(in)  :: rveig_atom(inc%lmmaxso, inc%natypd, 2)
    double complex, intent(out) :: proj

    !Locals
    integer        :: iatom
    double complex :: ckhelp(inc%lmmaxso), ztmp

    !Parameter
    double complex, parameter :: CONE=(1d0,0d0), CZERO=(0d0, 0d0)
 
    proj = CZERO
    do iatom=1,inc%natypd
      ztmp = CZERO
      ckhelp = CZERO

      call ZGEMM('N', 'N', inc%lmmaxso, 1, inc%lmmaxso, CONE, rhod_norm(:,:,iatom), inc%lmmaxso, rveig_atom(:,iatom,2), inc%lmmaxso, CZERO, ckhelp, inc%lmmaxso)
      call ZGEMM('C','N', 1, 1, inc%lmmaxso,CONE, rveig_atom(:,iatom,1), inc%lmmaxso, ckhelp, inc%lmmaxso, CZERO, ztmp, 1)
      proj = proj + ztmp
    end do!iatom

  end subroutine calc_proj_wavefunc




subroutine normeigv_new(almso, LV, RV)

  integer,        intent(in)    :: almso
  double complex, intent(inout) :: LV(almso,almso), RV(almso,almso)

  integer        :: lm1, lm2
  double complex :: snorm, norm, proj

  double complex, parameter :: cone=(1d0,0d0)

  !normalize correctly
  do lm2=1,almso
    snorm=sqrt(dot_product(LV(:,lm2), RV(:,lm2)))
    RV(:,lm2) = RV(:,lm2)/snorm
    LV(:,lm2) = LV(:,lm2)/conjg(snorm)
  end do!lm2

  !perform checks
  do lm1=1,almso
    norm = dot_product(LV(:,lm1),RV(:,lm1))
    if(abs(norm-CONE)>1d-3) write(*,*) "norm not 1, (lm, norm)=", lm1, norm
  end do

  do lm1=1,almso
    do lm2=lm1+1,almso
      proj = dot_product(LV(:,lm2),RV(:,lm1))
      if(abs(proj)>1d-4) write(*,*) "proj not 0, (lm1, lm2, proj)=", lm1, lm2, proj
    end do
  end do

end subroutine normeigv_new




subroutine normeigv_new_memopt(almso, nb_ev, LV, RV)

  integer,        intent(in)    :: almso, nb_ev
  double complex, intent(inout) :: LV(almso,nb_ev), RV(almso,nb_ev)

  integer        :: lm1, lm2
  double complex :: snorm, norm, proj

  double complex, parameter :: cone=(1d0,0d0)

  !normalize correctly
  do lm2=1,nb_ev
    snorm=sqrt(dot_product(LV(:,lm2), RV(:,lm2)))
    RV(:,lm2) = RV(:,lm2)/snorm
    LV(:,lm2) = LV(:,lm2)/conjg(snorm)
  end do!lm2

  !perform checks
  do lm1=1,nb_ev
    norm = dot_product(LV(:,lm1),RV(:,lm1))
    if(abs(norm-CONE)>1d-3) write(*,*) "norm not 1, (lm, norm)=", lm1, norm
  end do

  if(nb_ev>1) then ! with memory optimization number of eigen values can be <= 1
    do lm1=1,nb_ev
      do lm2=lm1+1,nb_ev
        proj = dot_product(LV(:,lm2),RV(:,lm1))
        if(abs(proj)>1d-4) write(*,*) "proj not 0, (lm1, lm2, proj)=", lm1, lm2, proj
      end do
    end do
  end if!almso>1

end subroutine normeigv_new_memopt




SUBROUTINE CHECKDEGENERACY(ALMSO, eps_degenerate, EIGW, ENT, ENTW)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checks, if eigenvalues the eigenvalues contained !
! in 'W' are degenerate. The returned array 'ENT'  !
! contains the order of degeneracy and ENTW(LM1,:) !
! contains the LM-indices of the other degenerate  !
! eigenvalues, corresponding to the one with index !
! LM1.                                             !
!                    Put into the sbroutine by     !
!                       b.zimmermann, jan.2011     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE

  INTEGER,          INTENT (IN)  :: ALMSO
  DOUBLE PRECISION, INTENT (IN)  :: eps_degenerate
  DOUBLE COMPLEX,   INTENT (IN)  :: EIGW(ALMSO)
  INTEGER,          INTENT (OUT) :: ENT(ALMSO), ENTW(ALMSO,ALMSO)

  INTEGER LM1, LM2, DONE(ALMSO), PNTR

  ENT(: )   = 0
  ENTW(:,:) = 0
  DONE(:)   = 0
  DO LM1=1,ALMSO
    IF (DONE(LM1)==0) THEN

      !Find degenerate LM's (and store in row of lowest LM)
      ENT(LM1)    = 1
      ENTW(1,LM1) = LM1
      DONE(LM1)   = 1
      PNTR        = 1
      DO WHILE (ENTW(PNTR,LM1)/=0)
        DO LM2=ENTW(PNTR,LM1)+1,ALMSO
          IF( (ABS(EIGW(ENTW(PNTR,LM1))-EIGW(LM2)) < eps_degenerate) .and. (DONE(LM2)==0) ) THEN
            ENT(LM1)=ENT(LM1)+1
            ENTW(ENT(LM1),LM1)=LM2
            DONE(LM2)=1
          END IF!eps_degenerate
        END DO!LM2
        PNTR = PNTR+1
      END DO!WHILE

      !Copy row from lowest LM to its degenerate partners
      DO LM2=2,ENT(LM1)
        ENTW(:,  ENTW(LM2,LM1)) = ENTW(:,LM1)
        ENTW(1,  ENTW(LM2,LM1)) = ENTW(LM2,LM1)
        ENTW(LM2,ENTW(LM2,LM1)) = LM1
        ENT(ENTW(LM2,LM1)) = ENT(LM1)
      END DO
    END IF!DONE
  END DO!LM1

!!TESTOUTPUT
!  DO LM1=1,ALMSO
!    write(*,*) LM1, ENT(LM1)
!  END DO
!  write(*,*) '---'
!  DO LM1=1,ALMSO
!    write(*,"((8I))") LM1, (ENTW(LM2,LM1), LM2=1,ENT(LM1))
!  END DO
!  STOP

END SUBROUTINE






end module mod_eigvects
