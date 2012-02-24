subroutine KKRMAT01(BZKP,NOFKS,GS,VOLCUB,VOLBZ, &
TMATLL,MSSQ, &
ITER, &
ALAT,NSYMAT,NAEZ,CLS,NACLS,RR,EZOA,ATOM, &
GINP,DGINP, &
NUMN0,INDN0,IAT, &
SPRS,PRSC,EKM,NOITER, &
QMRBOUND,IGUESS,BCP,CNVFAC, &
DTDE_LOCAL, &
GSXIJ, &
NXIJ,XCCPL,IXCP,ZKRXIJ, &
LLY_GRDT,TR_ALPH, &
LMPIC,LCOMM,LSIZE, &
!new input parameters after inc.p removal
prod_lmpid_smpid_empid, nthrds, &
lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
nxijd, nguessd, kpoibz, nrd, ekmd)

  use lloyds_formula_mod

  implicit none
  include 'mpif.h'
  ! ************************************************************************
  !   performs k-space integration,
  !   determines scattering path operator (g(k,e)-t**-1)**-1 and
  !   Greens function of the real system -> GS(*,*,*,*),

  !   NEW VERSION 10.99
  !   up -> left , down -> right, for decimation
  ! ------------------------------------------------------------------------

  integer, intent(in) :: prod_lmpid_smpid_empid
  integer, intent(in) :: nthrds
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

  !     .. parameters ..
  integer, parameter :: NSYMAXD = 48
  double complex :: CONE
  parameter      (CONE  = ( 1.0D0,0.0D0))
  double complex :: CZERO
  parameter      (CZERO=(0.0D0,0.0D0))

  !     ..
  !     .. SCALAR ARGUMENTS ..
  double precision:: ALAT
  double precision::VOLBZ
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
  !     ..
  !     .. ARRAY ARGUMENTS ..

  double complex :: TMATLL(lmmaxd,lmmaxd,NAEZ)

  double complex :: DGINP(lmmaxd,lmmaxd,NACLSD,NCLSD)
  double complex :: GINP (lmmaxd,lmmaxd,NACLSD,NCLSD)
  double complex :: GS   (lmmaxd,lmmaxd,NSYMAXD)
  double complex :: GSXIJ(lmmaxd,lmmaxd,NSYMAXD,NXIJD)

  integer        :: IXCP(NXIJD)
  double complex :: EIKRP(NACLSD)
  double complex :: EIKRM(NACLSD)

  ! .. Lloyd
  double complex :: DTDE_LOCAL(lmmaxd,lmmaxd)
  double complex :: MSSQ(lmmaxd,lmmaxd)

  double complex :: TR_ALPH
  double complex :: LLY_GRDT

  complex        :: PRSC(NGUESSD*lmmaxd,EKMD) ! array argument
  integer        :: SPRS(NGUESSD*lmmaxd+1,EKMD+1)! array argument

  double precision::BZKP(3,KPOIBZ)
  double precision::VOLCUB(KPOIBZ)
  double precision::RR(3,0:NRD)
  double precision::ZKRXIJ(48,3,NXIJD)
  double precision::CNVFAC(EKMD)

  integer:: NUMN0(NAEZ)
  integer:: INDN0(NAEZ,NACLSD)
  integer:: ATOM(NACLSD,*)
  integer:: CLS(*)
  integer:: EZOA(NACLSD,*)
  integer:: NACLS(*)


  !     .. LOCAL SCALARS ..
  double complex::TRACE
  double complex::TRACEK
  double complex::BZTR2
  double complex::CFCTORINV
  double precision::TWO_PI
  double precision::QMRBOUND
  double precision::DZNRM2   ! external function for 2-norm
  integer::ref_cluster_index
  integer::site_index
  integer::cluster_site_index
  integer::ILM
  integer::ISYM
  integer::symmetry_index
  integer::site_lm_index
  integer::IL1B
  integer::cluster_site_lm_index
  integer::IL2B
  integer::k_point_index
  integer::LM
  integer::LM1
  integer::LM2
  integer::LM3
  integer::XIJ
  integer::iteration_counter
  logical::XCCPL

  !     .. LOCAL ARRAYS ..
  double complex  ::G  (lmmaxd,lmmaxd) ! small
  double complex  ::TGH(lmmaxd) ! small
  double precision::N2B(lmmaxd) ! small

  !    double complex::GLLKE1(ALM,LMMAXD) !large
  !    double complex::DUMMY(ALM,LMMAXD)  !large
  !    double complex::GLLH(LMMAXD,NGTBD,NAEZD) ! large!
  !    double complex::GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*NBLCKD) ! large
  !    double complex::X0(ALM,LMMAXD) ! large

  double complex, allocatable, dimension(:,:) ::GLLKE1
  double complex, allocatable, dimension(:,:) ::DUMMY
  double complex, allocatable, dimension(:,:,:) ::GLLH
  double complex, allocatable, dimension(:,:) ::GLLHBLCK
  double complex, allocatable, dimension(:,:) ::X0

  !  .. local arrays for Lloyd's formula
  !double complex :: DGDE(LLYALM,LMMAXD)
  !double complex :: GLLKE_X(LLYALM,LMMAXD)
  !double complex :: DPDE_LOCAL(LLYALM,LMMAXD)
  double complex, allocatable, dimension(:,:) :: DGDE
  double complex, allocatable, dimension(:,:) :: GLLKE_X
  double complex, allocatable, dimension(:,:) :: DPDE_LOCAL

  !     ..
  !     .. EXTERNAL SUBROUTINES ..
  external CINIT,DLKE0,OUTTIME,DZNRM2
  !     ..
  !     .. INTRINSIC FUNCTIONS ..
  intrinsic ATAN,EXP
  !     ..
  !     .. MPI ..
    
  !     .. L-MPI
  integer:: LCOMM(prod_lmpid_smpid_empid)
  integer:: LSIZE(prod_lmpid_smpid_empid)
  integer:: LMPIC

  !     .. N-MPI
  integer:: IERR

  integer::        LMGF0D
  integer::        ALM
  integer::        NGTBD
  integer::        NBLCKD

  integer :: memory_stat
  logical :: memory_fail

  ! array dimensions
  LMGF0D= lmmaxd
  ALM = NAEZ*LMMAXD
  NGTBD = NACLSD*LMMAXD
  NBLCKD = XDIM*YDIM*ZDIM

  !-----------------------------------------------------------------------
  ! Allocate arrays
  !-----------------------------------------------------------------------
  memory_stat = 0
  memory_fail = .false.

  allocate(GLLKE1(ALM,LMMAXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(DUMMY (ALM,LMMAXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(GLLH(LMMAXD,NGTBD,NAEZ), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*NBLCKD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(X0(ALM,LMMAXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (LLY == 1) then
    allocate(DGDE      (ALM,LMMAXD), stat = memory_stat)
    if (memory_stat /= 0) memory_fail = .true.
    allocate(GLLKE_X   (ALM,LMMAXD), stat = memory_stat)
    if (memory_stat /= 0) memory_fail = .true.
    allocate(DPDE_LOCAL(ALM,LMMAXD), stat = memory_stat)
    if (memory_stat /= 0) memory_fail = .true.
  end if

  if (memory_fail .eqv. .true.) then
    write(*,*) "KKRMAT01: FATAL Error, failure to allocate memory."
    write(*,*) "       Probably out of memory."
    stop
  end if

  TWO_PI = 8.D0*ATAN(1.D0)    ! = 2*PI
  CFCTORINV = (CONE*TWO_PI)/ALAT

  BZTR2 = CZERO

  do 20 symmetry_index = 1,NSYMAXD
    call CINIT(LMMAXD**2,GS(1,1,symmetry_index))
    20 end do

    if (XCCPL) then
      do XIJ = 1, NXIJ
        do ISYM = 1,NSYMAT
          call CINIT(LMMAXD*LMMAXD,GSXIJ(1,1,ISYM,XIJ))
        enddo             ! ISYM = 1,NSYMAT
      enddo
    endif


    ! ---> use sit
    !      G(n,n',L,L')(-k) = G(n',n,L',L)(k)


!=======================================================================
    do 300 k_point_index = 1, NOFKS                       ! K-POINT-LOOP
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
    

      !================= Lloyd's Formula ====================================

      if (LLY == 1) then
        call CINIT(ALM*NGTBD,GLLH)

        do site_index = 1,NAEZ

          ref_cluster_index = CLS(site_index)

          call DLKE1(ALAT,NACLS,RR,EZOA(1,site_index), &
                     BZKP(1,k_point_index),ref_cluster_index,EIKRM,EIKRP, &
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

      endif


      !============ END Lloyd's Formula =====================================
    

      !       Fourier transform of reference clusters' Green's function
      !       (from real space to k-space GINP -> GLLH)

      call CINIT(ALM*NGTBD,GLLH)

      do site_index = 1,NAEZ
        ref_cluster_index = CLS(site_index)

        call DLKE1(ALAT,NACLS,RR,EZOA(1,site_index), &
                   BZKP(1,k_point_index),ref_cluster_index,EIKRM,EIKRP, &
                   nrd, naclsd)

        call DLKE0(site_index,GLLH,EIKRP,EIKRM, &
                   ref_cluster_index,NACLS,ATOM(1,site_index),NUMN0,INDN0, &
                   GINP(1,1,1,ref_cluster_index), &
                   naez, lmmaxd, naclsd)

      end do

      !=========== Lloyd's Formula ==========================================

      if (LLY == 1) then

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

      end if

      !===========END Lloyd's Formula=========================================


      ! The same calculation as some lines above is done all over again ???
      ! - NO! EIKRM and EIKRP are SWAPPED in call to DLKE0 !!!!

      call CINIT(ALM*NGTBD,GLLH)

      do site_index = 1,NAEZ
        ref_cluster_index = CLS(site_index)

        call DLKE1(ALAT,NACLS,RR,EZOA(1,site_index), &
                   BZKP(1,k_point_index),ref_cluster_index,EIKRM,EIKRP, &
                   nrd, naclsd)

        call DLKE0(site_index,GLLH,EIKRM,EIKRP, &
                   ref_cluster_index,NACLS,ATOM(1,site_index),NUMN0,INDN0, &
                   GINP(1,1,1,ref_cluster_index), &
                   naez, lmmaxd, naclsd)

      end do

      ! -------------- Calculation of (Delta_t * G_ref - 1) ---------------
      !                = inverse of scattering path operator
      !
      ! NUMN0(site_index) is the number of atoms in the reference cluster
      ! of atom/site 'site_index'
      ! INDN0 stores the index of the atom in the basis corresponding to
      ! the reference cluster atom
      ! -------------------------------------------------------------------

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

! ==> now GLLH holds the inverse of the scattering path operator (Delta_t * G_ref - 1)


! ####################################################################
! Start Lloyd's formula here
! ####################################################################
    
      ! dP(E,k)   dG(E,k)                   dT(E)
      ! ------- = ------- * T(E) + G(E,k) * -----
      !   dE        dE                       dE

      if (LLY == 1) then

       !calcDerivativeP(site_lm_size, lmmaxd, alat, &
       !                       DPDE_LOCAL, GLLKE_X, DGDE, DTmatDE_LOCAL, Tmat_local)
        call calcDerivativeP(ALM, lmmaxd, alat, &
                             DPDE_LOCAL, GLLKE_X, DGDE, DTDE_LOCAL, TMATLL(1,1,IAT))

      endif
!#######################################################################
! LLOYD
!#######################################################################


! Now solve the linear matrix equation A*X = b (b is also a matrix),
! where A = (Delta_t*G_ref - 1) (inverse of scattering path operator)
! and b = (-1) * Delta_t

! If the initial guess optimisation is used, X and b are modified, but
! the form of the matrix equation stays the same

!===================================================================
!===================================================================

! solve linear matrix equation in 5 subsequent steps

!===================================================================

! 1) find true residual tolerance by calculation of |b|
!    meaning of true residual tolerance:
!    the norm of |b| is used as convergence criterion
!    note that if the initial guess optimisation is used, the matrix
!    equation is solved using a modified b' - nevertheless the norm
!    |b| of the original right hand side of the equation is used
!    to test for convergence
    
      call CINIT(ALM*LMMAXD,DUMMY)
    
      do LM1=1,LMMAXD
        site_lm_index=LMMAXD*(IAT-1)+LM1
        do LM2=1,LMMAXD
          DUMMY(site_lm_index,LM2)=-TMATLL(LM1,LM2,IAT)
        enddo
      enddo

      !     store the norms of the columns of b in an array

      do LM2=1,LMMAXD
        N2B(LM2) = DZNRM2(NAEZ*LMMAXD,DUMMY(1,LM2),1)
      enddo
    
! ..
!===================================================================

!===================================================================
! 2) if IGUESS is activated perform intitial step of
!    intial guess - store original b in DUMMY and set up modified b'
    
      if (IGUESS == 1) then
        
        call CINIT(ALM*LMMAXD,X0)
        call CINIT(ALM*LMMAXD,DUMMY)
        
        do site_index=1,NAEZ
          do LM1=1,LMMAXD
            site_lm_index=LMMAXD*(site_index-1)+LM1
            do LM2=1,LMMAXD
              DUMMY(site_lm_index,LM2)=TMATLL(LM1,LM2,site_index)
            enddo
          enddo
        enddo
        
        if (ITER > 1) then
          call initialGuess_start( IAT, NUMN0, INDN0, &
                        TMATLL,GLLH,X0, &
                        PRSC(1,EKM + k_point_index),SPRS(1,EKM + k_point_index), &
                        naez, lmmaxd, naclsd, nguessd, nthrds)
        endif

      endif ! IGUESS == 1

! ..
!===================================================================

!===================================================================
! 3) if BCP is activated determine preconditioning matrix
!    GLLHBLCK ..
    
      call CINIT(LMMAXD*NATBLD*LMMAXD*NATBLD*NBLCKD,GLLHBLCK)
    
      if (BCP == 1) then

        call BCPWUPPER(GLLH,GLLHBLCK,NAEZ,NUMN0,INDN0, &
                       lmmaxd, nthrds, natbld, xdim, ydim, zdim, naclsd)
      endif
! ..
!===================================================================

!===================================================================
! 4) solve linear set of equations by iterative TFQMR scheme
    
      call CINIT(ALM*LMMAXD,GLLKE1)
    
      call MMINVMOD(GLLH,GLLKE1,TMATLL,NUMN0,INDN0,N2B, &
                    IAT,ITER,iteration_counter, &
                    GLLHBLCK,BCP,IGUESS,CNVFAC(EKM+k_point_index), &
                    QMRBOUND, &
                    naez, lmmaxd, naclsd, xdim, ydim, zdim, &
                    natbld, nthrds)
    
      NOITER = NOITER + iteration_counter
    
! ..
!===================================================================

!===================================================================
! 5) if IGUESS is activated perform final step of
!    intial guess           ~
!                  X = X  + X
!                       0          ..
    
      if (IGUESS == 1) then
        call initialGuess_finish(X0, PRSC(1,EKM+k_point_index),SPRS(1,EKM+k_point_index), &
                      GLLKE1, naez, lmmaxd, nguessd)
        
        do site_index=1,NAEZ
          do LM1=1,LMMAXD
            site_lm_index=LMMAXD*(site_index-1)+LM1
            do LM2=1,LMMAXD
              TMATLL(LM1,LM2,site_index)=DUMMY(site_lm_index,LM2)
            enddo
          enddo
        enddo
        
      endif
! ..
!===================================================================

! solved. Result in GLLKE1

!===================================================================
!===================================================================



!#######################################################################
! LLOYD calculate Trace of matrix ...
!#######################################################################
      if (LLY == 1) then
      !===================================================================
      !                /  -1    dM  \
      ! calculate  Tr  | M   * ---- |
      !                \        dE  /
      !===================================================================
        
        !call calcLloydTraceXRealSystem(DPDE_LOCAL, GLLKE1, inv_Tmat, TRACEK, site_lm_size, lmmaxd)
        call calcLloydTraceXRealSystem(DPDE_LOCAL, GLLKE1, MSSQ, TRACEK, ALM, lmmaxd)
        
        BZTR2 = BZTR2 + TRACEK*VOLCUB(k_point_index)  ! k-space integration
        
      endif
!#######################################################################
! LLOYD .
!#######################################################################
    
      !===================================================================

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
          do LM1=1,LMMAXD
            do LM2=1,LMMAXD
              GS(LM1,LM2,ISYM) = GS(LM1,LM2,ISYM) + &
              VOLCUB(k_point_index) * G(LM1,LM2)
            end do
          end do
        end do        ! ISYM = 1,NSYMAT
    
    
    
! ================================================================
          if (XCCPL) then

            ! ================================================================
            !       XCCPL communicate off-diagonal elements and multiply with
            !       exp-factor
            ! ================================================================
        
            call KKRJIJ( BZKP,VOLCUB,k_point_index, &
            NSYMAT,NAEZ,IAT, &
            NXIJ,IXCP,ZKRXIJ, &
            GLLKE1, &
            GSXIJ, &
            LMPIC, LCOMM, LSIZE, &
            lmmaxd, nxijd, prod_lmpid_smpid_empid)
        
! ================================================================
          endif
! ================================================================
    
!=======================================================================
    300 end do ! KPT = 1,NOFKS
!=======================================================================

  !=========== Lloyd's Formula =====================================

  if(LLY == 1)  then
    BZTR2 = BZTR2*NSYMAT/VOLBZ + TR_ALPH
    TRACE=CZERO

    call MPI_ALLREDUCE(BZTR2,TRACE,1,MPI_DOUBLE_COMPLEX,MPI_SUM, &
                       LCOMM(LMPIC),IERR)

    LLY_GRDT = TRACE
  endif
  !========== END Lloyd's Formula ==================================

  ! ----------------------------------------------------------------
  ! Deallocate arrays
  ! ----------------------------------------------------------------

  deallocate(GLLKE1)
  deallocate(DUMMY)
  deallocate(GLLH)
  deallocate(GLLHBLCK)
  deallocate(X0)

  if (LLY == 1) then
    deallocate(DGDE)
    deallocate(GLLKE_X)
    deallocate(DPDE_LOCAL)
  end if

end subroutine KKRMAT01
