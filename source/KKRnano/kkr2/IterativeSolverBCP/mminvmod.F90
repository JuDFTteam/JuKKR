! comment E.R.:
! use the following strategy:
! IGUESS = 0: always use 0 as starting solution
! IGUESS = 1 and SCITER=1: use 0
! IGUESS = 1 and SCITER /= 1: use X2 as initial guess
! on output: solution in X2

#define TESTNAN(X) if (any(X /= X)) then; write(*,*) #X, __LINE__; STOP; endif

subroutine MMINVMOD(GLLH1,X2, B ,NUMN0,INDN0, &
                    SCITER,ITCOUNT, &
                    GLLHBLCK,BCP,IGUESS, &
                    TOL, &
                    naezd, lmmaxd, naclsd, xdim, ydim, zdim, &
                    natbld)


  implicit none

  integer, intent(in) :: naezd
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd
  integer, intent(in) :: xdim
  integer, intent(in) :: ydim
  integer, intent(in) :: zdim
  integer, intent(in) :: natbld

  double complex :: GLLH1(LMMAXD,NACLSD*LMMAXD,NAEZD) ! in?
  double complex :: X2(NAEZD*LMMAXD,LMMAXD)           ! inout
  double complex :: B(NAEZD*LMMAXD,LMMAXD)            ! in, right-hand side
  integer:: NUMN0(NAEZD)                              ! in
  integer:: INDN0(NAEZD,NACLSD)                       ! in
  integer::SCITER                                     ! in
  integer::ITCOUNT                                    ! out
  double complex :: GLLHBLCK(NATBLD*LMMAXD, NATBLD*XDIM*YDIM*ZDIM*LMMAXD) ! inout, work-array
  integer::BCP                                        ! in
  integer::IGUESS                                     ! in
  double precision::CNVFAC                            ! inout - useless
  double precision::TOL                               ! in

!----------------- local variables --------------------------------------------

  INTEGER :: NLEN
  INTEGER :: NDIM

  double complex :: CONE
  double complex :: CZERO
  parameter         (CONE = (1.0D0,0.0D0),CZERO = (0.0D0,0.0D0))
  ! ..

  logical::QMRABS

  ! local scalars ..
  double complex :: ZTMP
  double precision::DTMP
  double precision::RESNAV
  double precision::TOLAV
  integer::NLIM
  integer::LM2
  integer::I
  integer::DIMVEC
  integer::IT
  integer::PROBE
  integer::ITPROBE
  logical::EXITIT
  ! ..

  ! local arrays ..

  !     small, local arrays with dimension LMMAXD
  double precision::N2B(LMMAXD)
  double complex :: RHO(LMMAXD)
  double complex :: ETA(LMMAXD)
  double complex :: BETA(LMMAXD)
  double complex :: ALPHA(LMMAXD)

  double precision::R0(LMMAXD)
  double precision::RESN(LMMAXD)
  double precision::VAR(LMMAXD)
  double precision::TAU(LMMAXD)
  double precision::COS(LMMAXD)
  double precision::TOLB(LMMAXD)

  integer::ITER(LMMAXD)
  logical::DONE(LMMAXD)

  double precision::HSTR(10)
  integer::HSTI(3)

  ! some large local arrays
  double complex, allocatable, dimension(:,:,:) :: VECS
  double complex, allocatable, dimension(:,:)   :: DUMMY
!IBM* ALIGN(32, VECS)

  ! external ..
  external           DZNRM2,ZDOTU,ZRANDN
  double complex :: ZDOTU
  double precision::   DZNRM2

  !------------------------------------------------------------------
  ! convergence parameters

  double precision :: rhs_norm
  double precision :: target_upper_bound
  double precision :: max_upper_bound
  double precision, parameter :: TEST_FACTOR = 100d0

  !=======================================================================
  ! INITIALIZATION I
  !=======================================================================
  rhs_norm = 0.0d0

  NDIM = NAEZD*LMMAXD
  NLEN = NAEZD*LMMAXD

  !     Allocate arrays
  allocate(VECS(NDIM,LMMAXD,9))
  allocate(DUMMY(LMMAXD*NAEZD,LMMAXD))

  TOLAV = 0.0
  QMRABS = .true.           ! QMRABS tolerance for residual norm is defined globally

  target_upper_bound = TOL * TEST_FACTOR

  CNVFAC  = 1000

  do I=1,3
    HSTI(I) = 1             ! HSTI and HSTR store the history of QMR-convergency,
  enddo                     ! which is the basis to identify stagnation
  do I=1,10
    HSTR(I) = 99+I
  enddo

  ITCOUNT = 0
  ITPROBE = 0

  !=======================================================================
  ! INITIALIZATION II
  !=======================================================================

  do I = 1,NDIM
    do LM2=1,LMMAXD
      VECS(I,LM2,2) = - B(I,LM2)
      !VECS(I,LM2,1) = CZERO
      do DIMVEC=3,9
        VECS(I,LM2,DIMVEC) = CZERO
      enddo
    enddo
  enddo

  NLIM = 2000  ! limit of max. 2000 iterations

  !--------------
  do LM2=1,LMMAXD
  !--------------
    
    N2B(LM2)   = DZNRM2(NDIM,VECS(1,LM2,2),1)           ! norm of right-hand-sight

    if (QMRABS) then
      TOLB(LM2)  = TOL
    else
      TOLB(LM2)  = MAX(TOL*N2B(LM2),1.0D-10)          ! adapted relative tolerance
    endif
    
    TOLAV      = TOLAV  + TOLB(LM2)/LMMAXD
    
    DONE(LM2) = .false.
    ITER(LM2) = 0
    PROBE     = 1
    
    if (IGUESS == 0 .or. SCITER == 1) then   ! start with zero initial guess
      VECS(:,LM2,5) = VECS(:,LM2,2)
      VECS(:,LM2,1) = CZERO
    endif

  !--------------
  enddo
  !--------------

  if (IGUESS == 1 .and. SCITER > 1) then     ! start with initial guess passed as X2
    ! VECS(:,:,9) = CZERO
    VECS(:,:,1) = X2

    do LM2 = 1,LMMAXD
      call ZCOPY(NDIM,VECS(1,LM2,1),1,DUMMY(1,LM2),1)
    enddo

    if (BCP == 1) then
      call APPBLCKCIRC(DUMMY,GLLHBLCK, &
      naezd,lmmaxd,natbld,xdim,ydim,zdim)
    endif

    call SPRSZMM( GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE, &
                  CONE,CZERO, &
                  VECS(1,1,9), &
                  naezd, lmmaxd, naclsd)

    ! r0 = b - Ax0 = v2 - v9
    VECS(:,:,5) = VECS(:,:,2) - VECS(:,:,9)

  end if

  !--------------
  do LM2=1,LMMAXD
  !--------------
    
    R0(LM2) = DZNRM2(NLEN,VECS(1,LM2,5),1)
    
    ! Check whether the auxiliary vector must be supplied.
    
    call ZRANDN (NLEN,VECS(1,LM2,3),1)
    
    !     Initialize the variables.
    
    RESN(LM2) = 1.0D0
    RHO(LM2)  = CONE
    VAR(LM2)  = 0.0D0
    ETA(LM2)  = CZERO
    TAU(LM2)  = R0(LM2) * R0(LM2)

    VECS(:,LM2,8) = CZERO
    VECS(:,LM2,4) = CZERO
    VECS(:,LM2,6) = CZERO
    
  !--------------
  enddo
  !--------------

  !============================================================================
  !============================================================================
  ! ITERATION

  do IT=1, NLIM
    
    !============================================================================
    !============================================================================
    
    !--------------
    do LM2=1,LMMAXD
      if ( .not. DONE(LM2)) then
    !--------------
            
        ZTMP      = ZDOTU(NLEN,VECS(1,LM2,3),1,VECS(1,LM2,5),1)
        BETA(LM2) = ZTMP / RHO(LM2)
        RHO(LM2)  = ZTMP

        call ZAXPBY(NLEN,VECS(1,LM2,4), &
                    BETA(LM2),VECS(1,LM2,4),CONE,VECS(1,LM2,8))

        call ZAXPBY(NLEN,VECS(1,LM2,6), &
                    CONE,VECS(1,LM2,5),BETA(LM2),VECS(1,LM2,6))
            
      !--------------
      endif
    enddo
    !--------------

    do LM2 = 1,LMMAXD
      call ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
    enddo

    if (BCP == 1) then
      call APPBLCKCIRC(DUMMY,GLLHBLCK, &
      naezd,lmmaxd,natbld,xdim,ydim,zdim)
    endif
    
    call SPRSZMM( GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE, &
                  CONE,CZERO, &
                  VECS(1,1,9), &
                  naezd, lmmaxd, naclsd)
    
    !     VECS(:,:,6) input vector to be multiplied by A = GLLH1
    !     VECS(:,:,9) result
    
    !--------------
    do LM2=1,LMMAXD
      if ( .not. DONE(LM2)) then
        !--------------
            
        call ZAXPBY(NLEN,VECS(1,LM2,4), &
        BETA(LM2),VECS(1,LM2,4),CONE,VECS(1,LM2,9))
            
        ZTMP = ZDOTU(NLEN,VECS(1,LM2,3),1,VECS(1,LM2,4),1)
            
        ALPHA(LM2) = RHO(LM2) / ZTMP

        ZTMP       = VAR(LM2) * ETA(LM2) / ALPHA(LM2)
        call ZAXPBY (NLEN,VECS(1,LM2,7), &
                     CONE,VECS(1,LM2,6),ZTMP,VECS(1,LM2,7))
        call ZAXPBY (NLEN,VECS(1,LM2,5), &
                     CONE,VECS(1,LM2,5),-ALPHA(LM2),VECS(1,LM2,9))
            
        DTMP = DZNRM2(NLEN,VECS(1,LM2,5),1)
        DTMP = DTMP * DTMP
        VAR(LM2)  = DTMP / TAU(LM2)
        COS(LM2)  = 1.0D0 / ( 1.0D0 + VAR(LM2) )
        TAU(LM2)  = DTMP * COS(LM2)
        ETA(LM2)  = ALPHA(LM2) * COS(LM2)
            
        call ZAXPBY(NLEN,VECS(1,LM2,1), &
                    CONE,VECS(1,LM2,1),ETA(LM2),VECS(1,LM2,7))
            
      !--------------
      endif
    enddo
    !--------------

    do LM2=1,LMMAXD
      if ( .not. DONE(LM2)) then
            
        call ZAXPBY(NLEN,VECS(1,LM2,6), &
                    CONE,VECS(1,LM2,6),-ALPHA(LM2),VECS(1,LM2,4))

        ZTMP = VAR(LM2) * COS(LM2)

        call ZAXPBY(NLEN,VECS(1,LM2,7), &
                    CONE,VECS(1,LM2,6),ZTMP,VECS(1,LM2,7))
            
      endif
    enddo

    
    EXITIT = .true.
    do LM2 = 1, LMMAXD
      if ( .not. DONE(LM2)) EXITIT = .false.
    enddo
    if (EXITIT) go to 67

    
    do LM2 = 1,LMMAXD
      call ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
    enddo

    if (BCP == 1) then
      call APPBLCKCIRC(DUMMY,GLLHBLCK, &
      naezd,lmmaxd,natbld,xdim,ydim,zdim)
    endif
    
    call SPRSZMM( GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE, &
                  CONE,CZERO, VECS(1,1,8), &
                  naezd, lmmaxd, naclsd)
    
    !     VECS(:,:,6) input vector to be multiplied by A = GLLH1
    !     VECS(:,:,8) result
    

    max_upper_bound = 0.0d0
    !--------------
    do LM2=1,LMMAXD
      if ( .not. DONE(LM2)) then
        !--------------
            
        call ZAXPBY(NLEN,VECS(1,LM2,5), &
                    CONE,VECS(1,LM2,5),-ALPHA(LM2),VECS(1,LM2,8))
            
        DTMP = DZNRM2(NLEN,VECS(1,LM2,5),1)
        DTMP = DTMP * DTMP
        VAR(LM2)  = DTMP / TAU(LM2)
        COS(LM2)  = 1.0D0 / ( 1.0D0 + VAR(LM2) )
        TAU(LM2)  = DTMP * COS(LM2)
        ETA(LM2)  = ALPHA(LM2) * COS(LM2)
            
        call ZAXPBY(NLEN,VECS(1,LM2,1), &
        CONE,VECS(1,LM2,1),ETA(LM2),VECS(1,LM2,7))
        
        rhs_norm = N2B(LM2)
        if (rhs_norm <= epsilon(1.0d0)) then
          rhs_norm = 1.0d0
        end if

        max_upper_bound = max(max_upper_bound, &
                          sqrt( (2*IT+1)*TAU(LM2)) / rhs_norm)
    !--------------
      endif
    enddo
    !--------------
    
    if (max_upper_bound <= target_upper_bound) then
      PROBE = IT ! probe residual
    else
      PROBE = IT+1 ! don't probe residual
    end if
    
    ! >>>>>>>>>>>
    if (MOD(IT,PROBE) == 0) then
        
      ! in case of right-preconditioning
      !                  -1
      !         r = A * M  * y - b      otherwise   r = A * y - b
      !                  2
      ! >>>>>>>>>>>
      ! has to be performed ..
        
      if (BCP == 1) then
            
        do LM2 = 1,LMMAXD
          call ZCOPY(NDIM,VECS(1,LM2,1),1,DUMMY(1,LM2),1)
        enddo
            
        call APPBLCKCIRC(DUMMY,GLLHBLCK,naezd,lmmaxd,natbld,xdim,ydim,zdim)
            
      else
            
        do LM2 = 1,LMMAXD
          call ZCOPY(NDIM,VECS(1,LM2,1),1,DUMMY(1,LM2),1)
        enddo
            
      endif
        
      call SPRSZMM( GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE, &
                    CONE,CZERO, &
                    VECS(1,1,9), &
                    naezd, lmmaxd, naclsd)
        
      !     VECS(:,:,1) input vector to be multiplied by A = GLLH1
      !     VECS(:,:,9) result
        
      !--------------
      ITPROBE = ITPROBE + 1
      HSTI(MOD(ITPROBE,3)+1) = IT
        
      RESNAV = 0.0D0

      do LM2=1,LMMAXD
        if ( .not. DONE(LM2)) then
          !--------------
                
          call ZAXPBY(NLEN,VECS(1,LM2,9), &
                      CONE,VECS(1,LM2,2),-CONE,VECS(1,LM2,9))

          RESN(LM2) = DZNRM2(NLEN,VECS(1,LM2,9),1) / R0(LM2)
                
          if (RESN(LM2) <= TOLB(LM2)) then
            DONE(LM2) = .true.
            ITER(LM2) = 2*IT
          endif
          RESNAV = RESNAV + RESN(LM2)/LMMAXD
                
        !--------------
        endif
      enddo
      !--------------

      HSTR(MOD(ITPROBE,10)+1) = RESNAV
        
      ! check if QMR stagnated .. if yes leave cycle ..
      if (MAX(HSTR(1),HSTR(2),HSTR(3),HSTR(4),HSTR(5), &
      HSTR(6),HSTR(7),HSTR(8),HSTR(9),HSTR(10)) &
      == MIN(HSTR(1),HSTR(2),HSTR(3),HSTR(4),HSTR(5), &
      HSTR(6),HSTR(7),HSTR(8),HSTR(9),HSTR(10))) then
        do LM2=1,LMMAXD
          if ( .not. DONE(LM2)) then
            DONE(LM2) = .true.
            ITER(LM2) = 2*IT
            write(6,*) 'stgn.atm. ',0 ,LM2,' w. ',RESN(LM2),TOLB(LM2)
          endif
        enddo
      endif
      ! ..

      ! PROBE  = MAX(1,INT( LOG( RESNAV/(TOLAV*CNVFAC) ) ) )
        
      EXITIT = .true.
      do LM2 = 1, LMMAXD
        if ( .not. DONE(LM2)) EXITIT = .false.
      enddo
      if (EXITIT) goto 67
        
    ! <<<<<<<<<<<<
    endif
  ! <<<<<<<<<<<<

!============================================================================
!============================================================================
  enddo

! ITERATION
!============================================================================
!============================================================================

67 continue

   write(*,*) RESN
   write (*,*) "Converged: ", IT
   !     Done.

   ! >>>>>>>>>>>

   ! in case of right-preconditioning
   !              -1
   !         x = M  * y
   !              2
   ! has to be performed ..

   if (BCP == 1) then
    
     do LM2 = 1,LMMAXD
       call ZCOPY(NDIM,VECS(1,LM2,1),1,DUMMY(1,LM2),1)
     enddo
    
     call APPBLCKCIRC(DUMMY,GLLHBLCK,naezd,lmmaxd,natbld,xdim,ydim,zdim)
    
     do LM2 = 1,LMMAXD
       call ZCOPY(NDIM,DUMMY(1,LM2),1,VECS(1,LM2,1),1)
     enddo
    
   endif

   ! <<<<<<<<<<

   X2 = CZERO

   do LM2 = 1,LMMAXD
     call ZCOPY(NDIM,VECS(1,LM2,1),1,X2(1,LM2),1)
     ITCOUNT = ITCOUNT + ITER(LM2)
   enddo


   ITPROBE = MAX(HSTI(1),HSTI(2),HSTI(3)) &
   - MIN(HSTI(1),HSTI(2),HSTI(3))
   if (CNVFAC < 1D+8 .and. CNVFAC > 1D-8) then
     CNVFAC = CNVFAC * 10.0**(MIN(ITPROBE,8)-5)
   elseif (CNVFAC >= 1D+8 .and. (MIN(ITPROBE,8)-5) < 0) then
     CNVFAC = CNVFAC * 10.0**(MIN(ITPROBE,8)-5)
   elseif (CNVFAC <= 1D-8 .and. (MIN(ITPROBE,8)-5) > 0) then
     CNVFAC = CNVFAC * 10.0**(MIN(ITPROBE,8)-5)
   endif

   deallocate(VECS)
   deallocate(DUMMY)

 end subroutine MMINVMOD
