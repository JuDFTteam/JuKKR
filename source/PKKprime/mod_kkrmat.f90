!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


MODULE mod_kkrmat
!----------------------------------------------------------------------------------+
! This module supplies routines to calculate the KKR matrix.                       !
!----------------------------------------------------------------------------------+

IMPLICIT NONE

  PRIVATE
  PUBLIC :: CalcKKRmat, CalcKKRmat2, InvertT, CalcTMAT, &
          & compute_kkr_eigenvectors, compute_kkr_eigenvectors2, &
          & compute_kkr_eigenvectors2_dk

CONTAINS

  SUBROUTINE CalcKKRmat2(inc, GLLKE, TMAT_1, MMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Constructs the KKR-matrix M=[ 1- G_LL'^r(k)*\delta t ],     !
  ! given the greens function and T-matrix.                     !
  !                                                             !
  !                             b.zimmermann, jan.2011          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use type_inc, only: inc_type
    IMPLICIT NONE
    type(inc_type), intent (IN)  :: inc
    DOUBLE COMPLEX, INTENT (IN)  :: GLLKE(inc%ALMSO, inc%ALMSO), TMAT_1(inc%ALMSO,inc%ALMSO)
    DOUBLE COMPLEX, INTENT (OUT) :: MMAT(inc%ALMSO, inc%ALMSO)

    DOUBLE COMPLEX :: CZERO, CONE, CONEM
    PARAMETER (CZERO=(0d0, 0d0), CONE=(1d0, 0d0), CONEM=(-1d0, 0d0))

    INTEGER LM1, LM2

    MMAT=CZERO

    ! calculate (-g*t)
    CALL ZGEMM('N','N',inc%ALMSO,inc%ALMSO,inc%ALMSO,CONEM,GLLKE,inc%ALMSO,TMAT_1,inc%ALMSO,CZERO,MMAT,inc%ALMSO)
    ! calculate (1-g*t)
    DO LM1=1,inc%almso
      MMAT(LM1,LM1)=CONE+MMAT(LM1,LM1)
    ENDDO

!   CALL ZGEMM('N','N',inc%ALMSO,inc%ALMSO,inc%ALMSO,CONE,GLLKE,inc%ALMSO,TMAT_1,inc%ALMSO,CZERO,MMAT,inc%ALMSO)
!   DO LM1=1,inc%almso
!     DO LM2=1,inc%almso
!       IF (LM1.EQ.LM2) THEN
!         MMAT(LM1,LM2)=CONE-MMAT(LM1,LM2)
!       ELSE
!         MMAT(LM1,LM2)=-MMAT(LM1,LM2)
!       ENDIF
!     ENDDO
!   ENDDO

  END SUBROUTINE CalcKKRmat2


  SUBROUTINE CalcTMAT(inc, alat, TMATLL, TMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Constructs T-Matrix T(E)_LL'                                !
  !                             b.zimmermann, mar.2011          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use type_inc, only: inc_type
    IMPLICIT NONE

    type(inc_type),   intent (IN)  :: inc
    DOUBLE PRECISION, INTENT(IN)   :: alat
    DOUBLE COMPLEX,   INTENT (IN)  :: TMATLL(inc%lmmaxso,inc%lmmaxso,inc%naez)
    DOUBLE COMPLEX,   INTENT (OUT) :: TMAT(inc%ALMSO, inc%ALMSO)

    INTEGER I1, IS1, LM1, IS2, LM2, IL1, IL2, IL1T, IL2T
    DOUBLE PRECISION :: RFCTOR

    RFCTOR = alat/(8.D0*ATAN(1.0D0))

    TMAT=(0d0,0d0)
    DO I1=1,inc%NAEZ
      DO IS1=1,inc%NSPD
        DO LM1=1,inc%LMMAX
          DO IS2=1,inc%NSPD
            DO LM2=1,inc%LMMAX
              IL1=(IS1-1)*inc%ALM+inc%LMMAX*(I1-1)+LM1
              IL2=(IS2-1)*inc%ALM+inc%LMMAX*(I1-1)+LM2
              IL1T=(IS1-1)*inc%LMMAX+LM1
              IL2T=(IS2-1)*inc%LMMAX+LM2
              TMAT(IL2,IL1)=TMATLL(IL2T,IL1T,I1)/RFCTOR
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE CalcTMAT

  SUBROUTINE CalcKKRmat(inc, GLLKE, TINVLL, MMAT)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Constructs the KKR-matrix M=[ G_LL'^r(k) - \delta t^{-1} ], !
  ! given the greens function and inverse T-matrix.             !
  !                                                             !
  !                             b.zimmermann, jan.2011          !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use type_inc, only: inc_type
    IMPLICIT NONE

    type(inc_type), intent (IN)  :: inc
    DOUBLE COMPLEX, INTENT (IN)  :: GLLKE(inc%ALMSO, inc%ALMSO), TINVLL(inc%LMMAXSO,inc%LMMAXSO,inc%NAEZ)
    DOUBLE COMPLEX, INTENT (OUT) :: MMAT(inc%ALMSO, inc%ALMSO)

    INTEGER IAT, IS1, LM1, IS2, LM2, IL1, IL2, IL1T, IL2T

    MMAT = GLLKE
    DO IAT=1,inc%NAEZ
     DO IS1=1,inc%NSPD
      DO LM1=1,inc%LMMAX
       DO IS2=1,inc%NSPD
        DO LM2=1,inc%LMMAX
         IL1=(IS1-1)*inc%ALM+inc%LMMAX*(IAT-1)+LM1
         IL2=(IS2-1)*inc%ALM+inc%LMMAX*(IAT-1)+LM2
         IL1T=(IS1-1)*inc%LMMAX+LM1
         IL2T=(IS2-1)*inc%LMMAX+LM2
         MMAT(IL1,IL2) = MMAT(IL1,IL2)-TINVLL(IL1T,IL2T,IAT)
        END DO
       END DO
      END DO
     END DO
    END DO

  END SUBROUTINE CalcKKRmat





  SUBROUTINE InvertT(inc, alat, TMATLL, TINVLL)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculates the inverse of the T-matrix !
  !                 b.zimmermann, Dez.2010 !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use type_inc, only: inc_type
    IMPLICIT NONE

    type(inc_type),  intent (IN)  :: inc
    DOUBLE PRECISION, INTENT(IN)  :: alat
    DOUBLE COMPLEX,   INTENT(IN)  :: TMATLL(inc%lmmaxso,inc%lmmaxso,inc%naez)
    DOUBLE COMPLEX,   INTENT(OUT) :: TINVLL(inc%lmmaxso,inc%lmmaxso,inc%naez)

    INTEGER          :: IPVT1(inc%lmmaxso), INFO, ii, LM1, LM2, IER
    DOUBLE PRECISION :: RFCTOR
    DOUBLE COMPLEX   :: TLL(inc%lmmaxso,inc%lmmaxso)



    RFCTOR = alat/ (8.D0*DATAN(1.0D0))
    TINVLL = (0d0, 0d0)

    DO ii = 1,inc%naez
      IF (inc%ins.EQ.0 .AND. inc%nspo.EQ.1) THEN !for spherical potential
        DO LM1 = 1,inc%lmmaxso
          TINVLL(LM1,LM1,ii) = RFCTOR/TMATLL(LM1,LM1,ii)
        END DO
      ELSE !for full potential (or spin-orbit coupling)
        DO LM1 = 1,inc%lmmaxso
          TINVLL(LM1,LM1,ii) = RFCTOR
        END DO
        CALL ZCOPY(inc%lmmaxso*inc%lmmaxso,TMATLL(1,1,ii),1,TLL,1)
!       invert t-matrix and count negative eigenvalues (real part)
        CALL ZGETRF(inc%lmmaxso,inc%lmmaxso,TLL,inc%lmmaxso,IPVT1,INFO)
        IF (INFO /= 0) STOP 'Error in InvertT, ZGETRF'
        CALL ZGETRS('N',inc%lmmaxso,inc%lmmaxso,TLL,inc%lmmaxso,IPVT1,TINVLL(1,1,ii),inc%lmmaxso,INFO)
        IF (INFO /= 0) STOP 'Error in InvertT, ZGETRS'
      END IF !spherical/nonspherical potential

!     write(985,'("IATOM= ",I0)') ii
!     do lm2=1,inc%lmmaxso
!       do lm1=1,inc%lmmaxso
!         write(985,'(2I8,2ES25.16)') lm1,lm2,TINVLL(lm1,lm2,ii)
!       end do
!     end do

    END DO !ii=1,NAEZD

  END SUBROUTINE InvertT





  subroutine compute_kkr_eigenvectors(inc, lattice, cluster, ginp, tinvll, kpoint, eigw, LVeig, RVeig)

    use type_inc,      only: inc_type
    use type_data,     only: lattice_type, cluster_type, tgmatrx_type
    use mod_dlke,      only: dlke0
    use mod_eigvects,  only: normeigv_new

    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    double precision,   intent(in) :: kpoint(3)
    double complex,     intent(in) :: ginp(inc%naclsd*inc%lmmax,inc%lmmax,inc%nclsd),&
                                    & tinvll(inc%lmmaxso,inc%lmmaxso,inc%naezd)

    double complex, intent(out) :: LVeig(inc%almso,inc%almso),&
                                 & RVeig(inc%almso,inc%almso),&
                                 & eigw(inc%almso)

    double complex   :: gllke(inc%almso,inc%almso),&
                      & gllketmp(inc%alm,inc%alm), &
                      & mmat(inc%almso,inc%almso)

    integer          :: ierr, naux
    double precision, allocatable :: daux(:)
    double complex,   allocatable :: aux(:)


    naux = 2*inc%almso**2+5*inc%almso
    allocate( daux(2*inc%almso), aux(naux), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating aux arrays'

    ! calculate the Greens function
    gllketmp = (0d0, 0d0)
    call dlke0( gllketmp, inc, lattice%alat, cluster%cls, cluster%nacls,&
              & cluster%rr, cluster%ezoa, cluster%atom, kpoint,         &
              & cluster%rcls, ginp                                      )

    ! for spin-orbit coupling, double the array G(k,E)_LL'
    gllke = (0d0,0d0)
    gllke(1:inc%alm, 1:inc%alm) = gllketmp
    if(inc%nspd==2) gllke(inc%alm+1:inc%almso,inc%alm+1:inc%almso) = gllketmp

    ! find the KKR-matrix
    mmat = (0d0,0d0)
    call CalcKKRmat(inc, gllke, tinvll, mmat)

    ! diagonalize the KKR-matrix and store the eigenvalues and eigenvectors
    eigw  = (0d0,0d0)
    LVeig = (0d0,0d0)
    RVeig = (0d0,0d0)
    call ZGEEV( 'V','V', inc%almso, mmat, inc%almso, eigw, LVeig,  &
              & inc%almso, RVeig, inc%almso, aux, naux, daux, ierr )
    if (ierr /= 0 ) stop "Problems with diagonalizing KKR-matrix"

    ! normalize the eigenvectors
    call normeigv_new(inc%almso, LVeig, RVeig)

  end subroutine compute_kkr_eigenvectors





  subroutine compute_kkr_eigenvectors2(inc, lattice, cluster, ginp, tmat, kpoint, eigw, LVeig, RVeig)

    use type_inc,      only: inc_type
    use type_data,     only: lattice_type, cluster_type, tgmatrx_type
    use mod_dlke,      only: dlke0
    use mod_eigvects,  only: normeigv_new

    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    double precision,   intent(in) :: kpoint(3)    
    double complex,     intent(in) :: ginp(inc%naclsd*inc%lmmax,inc%lmmax,inc%nclsd),&
                                    & tmat(inc%almso,inc%almso)

    double complex, intent(out) :: LVeig(inc%almso,inc%almso),&
                                 & RVeig(inc%almso,inc%almso),&
                                 & eigw(inc%almso)

    double complex   :: gllke(inc%almso,inc%almso),&
                      & gllketmp(inc%alm,inc%alm), &
                      & mmat(inc%almso,inc%almso)

    integer          :: ierr, naux
    double precision, allocatable :: daux(:)
    double complex,   allocatable :: aux(:)


    naux = 2*inc%almso**2+5*inc%almso
    allocate( daux(2*inc%almso), aux(naux), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating aux arrays'

    ! calculate the Greens function
    gllketmp = (0d0, 0d0)
    call dlke0( gllketmp, inc, lattice%alat, cluster%cls, cluster%nacls,&
              & cluster%rr, cluster%ezoa, cluster%atom, kpoint,         &
              & cluster%rcls, ginp                                      )

    ! for spin-orbit coupling, double the array G(k,E)_LL'
    gllke = (0d0,0d0)
    gllke(1:inc%alm, 1:inc%alm) = gllketmp
    if(inc%nspd==2) gllke(inc%alm+1:inc%almso,inc%alm+1:inc%almso) = gllketmp

    ! find the KKR-matrix
    mmat = (0d0,0d0)
    call CalcKKRmat2(inc, gllke, tmat, mmat)

    ! diagonalize the KKR-matrix and store the eigenvalues and eigenvectors
    eigw  = (0d0,0d0)
    LVeig = (0d0,0d0)
    RVeig = (0d0,0d0)
    call ZGEEV( 'V','V', inc%almso, mmat, inc%almso, eigw, LVeig,  &
              & inc%almso, RVeig, inc%almso, aux, naux, daux, ierr )
    if (ierr /= 0 ) stop "Problems with diagonalizing KKR-matrix"

    ! normalize the eigenvectors
    call normeigv_new(inc%almso, LVeig, RVeig)

  end subroutine compute_kkr_eigenvectors2



  subroutine compute_kkr_eigenvectors2_dk(inc, lattice, cluster, ginp, tmat, kpoint, LVeig, RVeig, delta_lambda)

    use type_inc,      only: inc_type
    use type_data,     only: lattice_type, cluster_type, tgmatrx_type
    use mod_dlke,      only: dlke0dk

    implicit none

    type(inc_type),     intent(in) :: inc
    type(lattice_type), intent(in) :: lattice
    type(cluster_type), intent(in) :: cluster
    double precision,   intent(in) :: kpoint(3)    
    double complex,     intent(in) :: ginp(inc%naclsd*inc%lmmax,inc%lmmax,inc%nclsd),&
                                    & tmat(inc%almso,inc%almso),&
                                    & LVeig(inc%almso),&
                                    & RVeig(inc%almso)
    double complex,     intent(out) :: delta_lambda(3)

    double complex   :: gllke(inc%almso,inc%almso),  &
                      & gllketmp(inc%alm,inc%alm), &
                      & ctmp1(inc%almso),            &
                      & ctmp2(inc%almso)

    integer          :: ierr, naux, ixyz
    double precision, allocatable :: daux(:)
    double complex,   allocatable :: aux(:)

    double complex, parameter :: CZERO=(0d0,0d0), CONE=(1d0,0d0)

    ! calculate t*RVeig
    ctmp2=(0d0,0d0)
    CALL ZGEMM('N','N',inc%ALMSO,1,inc%ALMSO,CONE,tmat,inc%ALMSO,RVeig,inc%ALMSO,CZERO,ctmp1,inc%ALMSO)

    !init
    delta_lambda = (0d0,0d0)

    do ixyz=1,3
      ! calculate the derivative of the Greens function
      gllketmp = (0d0, 0d0)
      call dlke0dk( gllketmp, inc, lattice%alat, cluster%cls, cluster%nacls,&
                  & cluster%rr, cluster%ezoa, cluster%atom, kpoint,         &
                  & cluster%rcls, ginp, ixyz                                )

      ! for spin-orbit coupling, double the array dG(k,E)_LL'
      gllke = (0d0,0d0)
      gllke(1:inc%alm, 1:inc%alm) = gllketmp
      if(inc%nspd==2) gllke(inc%alm+1:inc%almso,inc%alm+1:inc%almso) = gllketmp

      ! calculate dG*t*RVeig
      ctmp2=(0d0,0d0)
      CALL ZGEMM('N','N',inc%ALMSO,1,inc%ALMSO,CONE,gllke,inc%ALMSO,ctmp1,inc%ALMSO,CZERO,ctmp2,inc%ALMSO)
      ! calculate -LVeig*[dG*t*RVeig]
      delta_lambda(ixyz) = -dot_product(LVeig,ctmp2)
    end do!ixyz


  end subroutine compute_kkr_eigenvectors2_dk



END MODULE mod_kkrmat
