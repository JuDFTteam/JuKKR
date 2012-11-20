#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"

module kkrmat_mod
CONTAINS

! WARNING: Symmetry assumptions might have been used that are
! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)

subroutine KKRMAT01(BZKP,NOFKS,GS,VOLCUB, &
TMATLL,MSSQ, &
ITER, &
ALAT,NSYMAT,NAEZ,CLS,NACLS,RR,EZOA,ATOM, &
GINP,DGINP, &
NUMN0,INDN0,IAT, &
PRSC,EKM,NOITER, &
QMRBOUND,IGUESS,BCP, &
DTDE_LOCAL, &
GSXIJ, &
NXIJ,XCCPL,IXCP,ZKRXIJ, &
BZTR2, &
communicator, comm_size, &
lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
nxijd, nguessd, kpoibz, nrd, ekmd)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD
  implicit none

  !     .. parameters ..
  integer, parameter :: NSYMAXD = 48
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

  ! ************************************************************************
  !   performs k-space integration,
  !   determines scattering path operator (g(k,e)-t**-1)**-1 and
  !   Greens function of the real system -> GS(*,*,*,*),

  !   NEW VERSION 10.99
  !   up -> left , down -> right, for decimation
  ! ------------------------------------------------------------------------

  integer, intent(in) :: communicator
  integer, intent(in) :: comm_size

  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd  ! max. number of atoms in reference cluster
  integer, intent(in) :: nclsd   ! number of reference clusters
  integer, intent(in) :: xdim
  integer, intent(in) :: ydim
  integer, intent(in) :: zdim
  integer, intent(in) :: natbld  ! number of atoms in preconditioning blocks
  integer, intent(in) :: LLY
  integer, intent(in) :: nxijd   ! max. number of atoms in cluster for exchange coupling-calculation
  integer, intent(in) :: nguessd
  integer, intent(in) :: kpoibz
  integer, intent(in) :: nrd
  integer, intent(in) :: ekmd

  !     ..
  !     .. SCALAR ARGUMENTS ..
  double precision:: ALAT

  integer::NAEZ
  integer::NOFKS
  integer::NSYMAT
  integer::IGUESS
  integer::BCP
  integer::ITER
  integer::NXIJ
  integer::EKM
  integer::NOITER
  integer::IAT

  double complex :: TMATLL(lmmaxd,lmmaxd,NAEZ)

  double complex :: DGINP(lmmaxd,lmmaxd,NACLSD,NCLSD)
  double complex :: GINP (lmmaxd,lmmaxd,NACLSD,NCLSD)
  double complex :: GS   (lmmaxd,lmmaxd,NSYMAXD)
  double complex :: GSXIJ(lmmaxd,lmmaxd,NSYMAXD,NXIJD)

  integer        :: IXCP(NXIJD)

  ! .. Lloyd
  double complex :: DTDE_LOCAL(lmmaxd,lmmaxd)
  double complex :: MSSQ(lmmaxd,lmmaxd)

  complex        :: PRSC(NGUESSD*lmmaxd,EKMD) ! array argument

  double precision::BZKP(3,KPOIBZ)
  double precision::VOLCUB(KPOIBZ)
  double precision::RR(3,0:NRD)
  double precision::ZKRXIJ(48,3,NXIJD)

  integer:: NUMN0(NAEZ)
  integer:: INDN0(NAEZ,NACLSD)
  integer:: ATOM(:,:) ! dim naclsd, naez?
  integer:: CLS(:)         ! dim *
  integer:: EZOA(:,:) ! dim naclsd, naez?
  integer:: NACLS(:)

  double complex::BZTR2
  double precision::QMRBOUND

  logical::XCCPL

! ------- local ----------

  double complex :: EIKRP(NACLSD)
  double complex :: EIKRM(NACLSD)
  integer::k_point_index
  double complex::TRACEK

  double complex, allocatable, dimension(:,:) ::GLLKE1
  double complex, allocatable, dimension(:,:,:) ::GLLH
  double complex, allocatable, dimension(:,:) ::GLLHBLCK

!IBM* ALIGN(32, GLLH)
!IBM* ALIGN(32, GLLHBLCK)

  !  .. local arrays for Lloyd's formula
  double complex, allocatable, dimension(:,:) :: DGDE
  double complex, allocatable, dimension(:,:) :: GLLKE_X
  double complex, allocatable, dimension(:,:) :: DPDE_LOCAL

  integer::        site_lm_size
  integer::        NGTBD
  integer::        NBLCKD

  integer :: memory_stat
  logical :: memory_fail

  ! array dimensions
  site_lm_size = NAEZ*LMMAXD
  NGTBD = NACLSD*LMMAXD
  NBLCKD = XDIM*YDIM*ZDIM

  !-----------------------------------------------------------------------
  ! Allocate arrays
  !-----------------------------------------------------------------------
  memory_stat = 0
  memory_fail = .false.

  allocate(GLLKE1(site_lm_size,LMMAXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(GLLH(LMMAXD,NGTBD,NAEZ), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*NBLCKD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (LLY == 1) then
    allocate(DGDE      (site_lm_size,LMMAXD), stat = memory_stat)
    if (memory_stat /= 0) memory_fail = .true.
    allocate(GLLKE_X   (site_lm_size,LMMAXD), stat = memory_stat)
    if (memory_stat /= 0) memory_fail = .true.
    allocate(DPDE_LOCAL(site_lm_size,LMMAXD), stat = memory_stat)
    if (memory_stat /= 0) memory_fail = .true.
  end if

  if (memory_fail .eqv. .true.) then
    write(*,*) "KKRMAT01: FATAL Error, failure to allocate memory."
    write(*,*) "       Probably out of memory."
    stop
  end if


  ! WARNING: Symmetry assumptions might have been used that are
  ! not valid in cases of non-local potential (e.g. for Spin-Orbit coupling)
  ! ---> use sit
  !      G(n,n',L,L')(-k) = G(n',n,L',L)(k)

  BZTR2 = CZERO
  GS = CZERO

  if (XCCPL) then
    GSXIJ = CZERO
  endif

  TESTARRAYLOG(3, GINP)

!==============================================================================
  do k_point_index = 1, NOFKS                       ! K-POINT-LOOP
!==============================================================================

    ! Get the scattering path operator for k-point BZKP(:, k_point_index)
    ! output: GLLKE1, NOITER
    ! inout: PRSC
    ! inout (temporary arrays): GLLH, GLLHBLCK

    call kloopbody( GLLKE1, PRSC(:, EKM + k_point_index), NOITER, &
                   BZKP(:, k_point_index), TMATLL, GINP, ALAT, IGUESS, &
                   BCP, NAEZ, ATOM, EZOA, RR, CLS, INDN0, &
                   NUMN0, EIKRM, EIKRP, GLLH, GLLHBLCK, &
                   IAT, ITER, QMRBOUND, NACLS, lmmaxd, nguessd, naclsd, &
                   natbld, nrd, nclsd, xdim, ydim, zdim)

    ! ----------- Integrate Scattering Path operator over k-points --> GS -----
    ! Note: here k-integration only in irreducible wedge
    call greenKSummation(GLLKE1, GS, VOLCUB(k_point_index), &
                         IAT, NSYMAT, naez, lmmaxd)
    ! -------------------------------------------------------------------------


!    if (LLY == 1) then
!      call lloydTraceK( TRACEK, TMATLL, MSSQ, GLLKE1, GINP, DGINP, &
!                       BZKP(:,k_point_index), DTDE_LOCAL,ALAT, ATOM, &
!                       CLS, EZOA, IAT, INDN0, NACLS, NAEZ, NUMN0, RR, EIKRM, &
!                       EIKRP, GLLH, DPDE_LOCAL, DGDE, GLLKE_X, &
!                       lmmaxd, naclsd, nclsd, nrd)
!    endif
!
!    if (LLY == 1) then
!      BZTR2 = BZTR2 + TRACEK*VOLCUB(k_point_index) ! k-space integration
!    end if
!
!    if (XCCPL) then
!
!       ! ================================================================
!       !       XCCPL communicate off-diagonal elements and multiply with
!       !       exp-factor
!       ! ================================================================
!
!      call KKRJIJ( BZKP(:,k_point_index),VOLCUB(k_point_index), &
!      NSYMAT,NAEZ,IAT, &
!      NXIJ,IXCP,ZKRXIJ, &
!      GLLKE1, &
!      GSXIJ, &
!      communicator, comm_size, &
!      lmmaxd, nxijd)
!
!    endif

!==============================================================================
  end do ! KPT = 1,NOFKS
!==============================================================================

  ! ----------------------------------------------------------------
  ! Deallocate arrays
  ! ----------------------------------------------------------------

  deallocate(GLLKE1)
  deallocate(GLLH)
  deallocate(GLLHBLCK)

  if (LLY == 1) then
    deallocate(DGDE)
    deallocate(GLLKE_X)
    deallocate(DPDE_LOCAL)
  end if

end subroutine KKRMAT01


!------------------------------------------------------------------------------
!> Calculate scattering path operator for 'kpoint'
subroutine kloopbody( GLLKE1, PRSC_k, NOITER, kpoint, TMATLL, GINP, ALAT, IGUESS, &
                     BCP, NAEZ, ATOM, EZOA, RR, CLS, INDN0, &
                     NUMN0, EIKRM, EIKRP, GLLH, GLLHBLCK, &
                     IAT, ITER, QMRBOUND, NACLS, lmmaxd, nguessd,naclsd, natbld, nrd, nclsd, xdim, ydim, zdim)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD
  use TEST_lcutoff_mod
  use initialGuess_store_mod
  implicit none

  integer, intent(in) :: xdim
  integer, intent(in) :: ydim
  integer, intent(in) :: zdim
  integer, intent(in) :: naclsd
  integer :: NAEZ
  integer, intent(in) :: nclsd
  integer, intent(in) :: nguessd
  double precision :: ALAT
  integer :: ATOM(:,:)         ! dim: naclsd, *
  integer :: BCP
  double precision :: kpoint(3)
  integer :: CLS(:)
  double complex :: EIKRM(naclsd)   ! dim: naclsd
  double complex :: EIKRP(naclsd)
  integer :: EZOA(NACLSD,*) ! dim naclsd,*
  doublecomplex :: GINP(lmmaxd,lmmaxd,NACLSD,NCLSD) ! dim: lmmaxd, lmmaxd, naclsd, nclsd
  double complex :: GLLH(LMMAXD,NACLSD*LMMAXD,NAEZ) ! dim: lmmaxd, naclsd*lmmaxd, naez
  double complex :: GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*XDIM*YDIM*ZDIM) ! dim: lmmaxd*natbld, lmmaxd*natbld*xdim*ydim*zdim
  double complex :: GLLKE1(NAEZ*LMMAXD,LMMAXD)
  integer :: IAT
  integer :: IGUESS
  integer :: INDN0(NAEZ,NACLSD)
  integer :: ITER
  integer, intent(in) :: lmmaxd
  integer :: NACLS(:)
  integer, intent(in) :: natbld
  integer :: NOITER
  integer, intent(in) :: nrd
  integer :: NUMN0(NAEZ)
  complex :: PRSC_k (NGUESSD*lmmaxd)
  double precision :: QMRBOUND
  double precision :: RR(3,0:NRD)
  doublecomplex :: TMATLL(lmmaxd,lmmaxd,NAEZ)

  !-------- local ---------
  double complex, parameter :: CONE = ( 1.0D0,0.0D0)
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

  integer :: iteration_counter
  integer :: ref_cluster_index
  integer :: site_index

  !=======================================================================
    
  ! ---> fourier transformation
    
  !     added by h.hoehler 3.7.2002
    
  !                                                     n   0          n
  !     define fourier transform as g mu mu'= ( sum_n g mu mu' exp(-iKR )
  !                                   L  L'             L   L'
    
  !                                             n   0           n
  !                                 +   sum_n g mu'mu exp(-iK(-R ))) *0.5
  !                                             L'  L
    
  !     this operation has to be done to satisfy e.g. the point symmetry!
  !     application of fourier transformation is just an approximation
  !     for the tb system, since the transl. invariance is not satisfied.
    


  ! The same calculation as with lloyds formula is done all over again ???
  ! - NO! EIKRM and EIKRP are SWAPPED in call to DLKE0 !!!!

  GLLH = CZERO

  do site_index = 1,NAEZ
    ref_cluster_index = CLS(site_index)

    call DLKE1(ALAT,NACLS,RR,EZOA(1,site_index), &
               kpoint,ref_cluster_index,EIKRM,EIKRP, &
               nrd, naclsd)

    call DLKE0(site_index,GLLH,EIKRM,EIKRP, &
               ref_cluster_index,NACLS,ATOM(:,site_index),NUMN0,INDN0, &
               GINP(1,1,1,ref_cluster_index), &
               naez, lmmaxd, naclsd)

  end do

  TESTARRAYLOG(3, GLLH)

  !----------------------------------------------------------------------------
  if (cutoffmode == 0 .or. cutoffmode == 1) then
    call generateCoeffMatrix(GLLH, NUMN0, INDN0, TMATLL, NAEZ, lmmaxd, naclsd)
  else
    call generateCoeffMatrixCROPPED(GLLH, NUMN0, INDN0, TMATLL, NAEZ, lmmaxd, naclsd, lmarray)
  end if

  if (cutoffmode == 1) then
    call cropGLLH(GLLH, lmmaxd, naclsd, naez, lmarray, numn0, indn0)
  end if

  TESTARRAYLOG(3, GLLH)
  !----------------------------------------------------------------------------

  ! ==> now GLLH holds (Delta_t * G_ref - 1)


  ! Now solve the linear matrix equation A*X = b (b is also a matrix),
  ! where A = (Delta_t*G_ref - 1) (inverse of scattering path operator)
  ! and b = (-1) * Delta_t

  ! If the initial guess optimisation is used, X and b are modified, but
  ! the form of the matrix equation stays the same

  !===================================================================
  ! 1) if IGUESS is activated, get initial guess from last iteration
    
  if (IGUESS == 1) then
        
    if (ITER > 1) then
      call initialGuess_load(GLLKE1, PRSC_k)
    endif

  endif ! IGUESS == 1

  !===================================================================
  ! 2) if BCP is activated determine preconditioning matrix
  !    GLLHBLCK ..
    
  GLLHBLCK = CZERO
    
  if (BCP == 1) then

    call BCPWUPPER(GLLH,GLLHBLCK,NAEZ,NUMN0,INDN0, &
                   lmmaxd, natbld, xdim, ydim, zdim, naclsd)
  endif

  !===================================================================
  ! 3) solve linear set of equations by iterative TFQMR scheme
  !    solve (\Delta t * G_ref - 1) X = - \Delta t
  !    the solution X is the scattering path operator
    
  call MMINVMOD(GLLH,GLLKE1,TMATLL,NUMN0,INDN0, &
                IAT,ITER,iteration_counter, &
                GLLHBLCK,BCP,IGUESS, &
                QMRBOUND, &
                naez, lmmaxd, naclsd, xdim, ydim, zdim, &
                natbld)
    
  NOITER = NOITER + iteration_counter

  !===================================================================
  ! 4) if IGUESS is activated save solution for next iteration
    
  if (IGUESS == 1) then
    call initialGuess_save(PRSC_k, GLLKE1)
  endif

  !===================================================================
  ! solved. Result in GLLKE1
  TESTARRAYLOG(3, GLLKE1)
  TESTARRAYLOG(3, TMATLL)

end subroutine

!------------------------------------------------------------------------------
!> Generates matrix (\Delta T G_ref - 1).
!> on input: GLLH contains G_ref, on output: GLLH contains coefficient matrix
subroutine generateCoeffMatrix(GLLH, NUMN0, INDN0, TMATLL, NAEZ, lmmaxd, naclsd)
  implicit none

  double complex, parameter :: CONE  = ( 1.0D0,0.0D0)
  double complex, parameter :: CZERO = ( 0.0D0,0.0D0)

  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd
  integer :: NAEZ
  double complex :: GLLH(LMMAXD,NACLSD*LMMAXD,NAEZ)
  integer :: INDN0(NAEZ,NACLSD)
  integer :: NUMN0(NAEZ)
  doublecomplex :: TMATLL(lmmaxd,lmmaxd,NAEZ)

  !---------- local --------------
  double complex :: TGH(lmmaxd)
  integer :: IL1B
  integer :: IL2B
  integer :: LM1
  integer :: LM2
  integer :: LM3
  integer :: site_index
  integer :: site_lm_index
  integer :: cluster_site_index
  integer :: cluster_site_lm_index

  ! -------------- Calculation of (Delta_t * G_ref - 1) ---------------
  !
  !
  ! NUMN0(site_index) is the number of atoms in the reference cluster
  ! of atom/site 'site_index' (inequivalent atoms only!)
  ! INDN0 stores the index of the atom in the basis corresponding to
  ! the reference cluster atom (inequivalent atoms only!)
  ! -------------------------------------------------------------------

  !$omp parallel do private(site_index, site_lm_index, cluster_site_index, &
  !$omp                     cluster_site_lm_index, IL1B, IL2B, &
  !$omp                     LM1, LM2, LM3, TGH)
  do site_index=1,NAEZ
    IL1B=LMMAXD*(site_index-1)
    do cluster_site_index=1,NUMN0(site_index)
      do LM2=1,LMMAXD
        cluster_site_lm_index=LMMAXD*(cluster_site_index-1)+LM2
        IL2B=LMMAXD*(INDN0(site_index,cluster_site_index)-1)+LM2
        do LM1=1,LMMAXD
          TGH(LM1) = CZERO
          do LM3=1,LMMAXD
            TGH(LM1)=TGH(LM1)+TMATLL(LM1,LM3,site_index)*GLLH(LM3,cluster_site_lm_index,site_index)
          enddo
        enddo

        do LM1=1,LMMAXD
          site_lm_index=IL1B+LM1
          GLLH(LM1,cluster_site_lm_index,site_index) = TGH(LM1)

          if (site_lm_index == IL2B) then
            ! substract 1 only at the 'diagonal'
            GLLH(LM1,cluster_site_lm_index,site_index) = GLLH(LM1,cluster_site_lm_index,site_index) - CONE
          endif

        enddo

      enddo
    enddo
  enddo
  !$omp end parallel do
end subroutine


!------------------------------------------------------------------------------
!> Summation of Green's function over k-points. Has to be called for every k-point
!> TODO: it would be better to do the k-space-symmetry treatment separately ???
!> This routine creates NSYMAT copies of the same solution
!> Set GS to 0 before first call
!> in: GLLKE1
!> inout: GS (set to 0 before first call)
subroutine greenKSummation(GLLKE1, GS, k_point_weight, IAT, NSYMAT, naez, lmmaxd)
  implicit none
  integer, parameter :: NSYMAXD = 48

  integer, intent(in) :: naez
  integer, intent(in) :: lmmaxd

  double complex :: GLLKE1(NAEZ*LMMAXD,LMMAXD)
  double complex :: GS(lmmaxd,lmmaxd,NSYMAXD)
  integer :: IAT

  integer :: NSYMAT
  double precision :: k_point_weight

  ! -------- local ------------------
  double complex :: G(lmmaxd,lmmaxd)
  integer :: LM
  integer :: LM1
  integer :: LM2
  integer :: ILM
  integer :: ISYM

  !   combined atom/lm index
  ILM = LMMAXD*(IAT-1) + 1

  !                                      nn
  !         Copy the diagonal elements G_LL' of the Green's-function,
  !         dependent on (k,E) into matrix G
  !         (n = n' = IAT)

  do LM = 1,LMMAXD
    call ZCOPY(LMMAXD,GLLKE1(ILM,LM),1,G(1,LM),1)
  end do

    !         Perform the k-space integration for diagonal element of
    !         Green's function of atom IAT

    do ISYM = 1,NSYMAT
      do LM2=1,LMMAXD
        do LM1=1,LMMAXD
          GS(LM1,LM2,ISYM) = GS(LM1,LM2,ISYM) + k_point_weight * G(LM1,LM2)
        end do
      end do
    end do        ! ISYM = 1,NSYMAT
end subroutine


!------------------------------------------------------------------------------
!> Calculates contribution to the Lloyd's formula trace for k = 'kpoint'
!> in: GLLKE1, MSSQ, kpoint, DGINP, GINP, TMATLL, DTDE_LOCAL
!> out: TRACEK
!> inout (temporary arrays): GLLH, DPDE_LOCAL, EIKRM, EIKRP, GLLKE_X, DGDE
subroutine lloydTraceK( TRACEK, TMATLL, MSSQ, GLLKE1, GINP, DGINP, kpoint, &
                       DTDE_LOCAL, ALAT, ATOM, CLS, EZOA, IAT, INDN0, NACLS, &
                       NAEZ, NUMN0, RR, EIKRM, EIKRP, GLLH, DPDE_LOCAL, DGDE, &
                       GLLKE_X, lmmaxd, naclsd, nclsd, nrd)

  use lloyds_formula_mod
  implicit none

  integer, intent(in) :: nclsd
  double precision :: ALAT
  integer :: ATOM(NACLSD,*)
  double precision :: kpoint(3)
  integer :: CLS(*)

  double complex :: DGDE(NAEZ*LMMAXD, LMMAXD)
  double complex :: DGINP(lmmaxd,lmmaxd,NACLSD,NCLSD)
  double complex :: DPDE_LOCAL(NAEZ*LMMAXD, LMMAXD)

  doublecomplex :: EIKRM(NACLSD)
  doublecomplex :: EIKRP(NACLSD)
  integer :: EZOA(NACLSD,*)
  double complex :: GINP(lmmaxd,lmmaxd,NACLSD,NCLSD)
  double complex :: GLLKE_X(NAEZ*LMMAXD, LMMAXD)
  double complex :: GLLH(LMMAXD,NACLSD*LMMAXD,NAEZ)
  double complex :: GLLKE1(NAEZ*LMMAXD,LMMAXD)
  integer :: IAT
  integer :: INDN0(NAEZ,NACLSD)
  double complex :: DTDE_LOCAL(lmmaxd,lmmaxd)
  integer, intent(in) :: lmmaxd
  double complex :: MSSQ(lmmaxd,lmmaxd)
  integer :: NACLS(*)
  integer, intent(in) :: naclsd
  integer :: NAEZ
  integer, intent(in) :: nrd
  integer :: NUMN0(NAEZ)
  double precision :: RR(3,0:NRD)
  double complex :: TMATLL(lmmaxd,lmmaxd,NAEZ)
  double complex :: TRACEK

  !----- local ------
  double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

  integer :: cluster_site_index
  integer :: cluster_site_lm_index
  integer :: ref_cluster_index
  integer :: site_index
  integer :: site_lm_index
  integer :: site_lm_size
  integer :: LM1
  integer :: LM2

  site_lm_size = NAEZ*lmmaxd

  GLLH = CZERO

  do site_index = 1,NAEZ

    ref_cluster_index = CLS(site_index)

    call DLKE1(ALAT,NACLS,RR,EZOA(1,site_index), &
               kpoint,ref_cluster_index,EIKRM,EIKRP, &
               nrd, naclsd)

    call DLKE0(site_index,GLLH,EIKRP,EIKRM, &
               ref_cluster_index,NACLS,ATOM(1,site_index),NUMN0,INDN0,DGINP(1,1,1,ref_cluster_index), &
               naez, lmmaxd, naclsd)
  end do

  DGDE = CZERO
  do site_index=1,NAEZ
    do cluster_site_index=1,NUMN0(site_index)
      do LM2=1,LMMAXD
        cluster_site_lm_index=LMMAXD*(cluster_site_index-1)+LM2

          if (INDN0(site_index,cluster_site_index) == IAT) then
            do LM1=1,LMMAXD
              site_lm_index=LMMAXD*(site_index-1)+LM1
              DGDE(site_lm_index,LM2)= GLLH(LM1,cluster_site_lm_index,site_index)
            end do
          end if

      enddo
    enddo
  enddo

  !       Fourier transform of reference clusters' Green's function
  !       (from real space to k-space GINP -> GLLH)

  GLLH = CZERO

  do site_index = 1,NAEZ
  ref_cluster_index = CLS(site_index)

  call DLKE1(ALAT,NACLS,RR,EZOA(1,site_index), &
             kpoint,ref_cluster_index,EIKRM,EIKRP, &
             nrd, naclsd)

  call DLKE0(site_index,GLLH,EIKRP,EIKRM, &
             ref_cluster_index,NACLS,ATOM(1,site_index),NUMN0,INDN0, &
             GINP(1,1,1,ref_cluster_index), &
             naez, lmmaxd, naclsd)

  end do

  GLLKE_X = CZERO
  do site_index=1,NAEZ
    do cluster_site_index=1,NUMN0(site_index)
      do LM2=1,LMMAXD
        cluster_site_lm_index=LMMAXD*(cluster_site_index-1)+LM2

        if (INDN0(site_index,cluster_site_index) == IAT) then
          do LM1=1,LMMAXD
            site_lm_index=LMMAXD*(site_index-1)+LM1
            GLLKE_X(site_lm_index,LM2)= GLLH(LM1,cluster_site_lm_index,site_index)
          end do
        endif

      enddo
    enddo
  enddo


  ! dP(E,k)   dG(E,k)                   dT(E)
  ! ------- = ------- * T(E) + G(E,k) * -----
  !   dE        dE                       dE

  call calcDerivativeP(site_lm_size, lmmaxd, alat, &
                       DPDE_LOCAL, GLLKE_X, DGDE, DTDE_LOCAL, TMATLL(1,1,IAT))

  !===================================================================
  !                /  -1    dM  \
  ! calculate  Tr  | M   * ---- |
  !                \        dE  /
  !===================================================================

  call calcLloydTraceXRealSystem(DPDE_LOCAL, GLLKE1, MSSQ, TRACEK, site_lm_size, lmmaxd)

end subroutine

end module
