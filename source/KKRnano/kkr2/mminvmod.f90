!***********************************************************************

subroutine MMINVMOD(GLLH1,X2,TMATLL,NUMN0,INDN0,N2B, &
                    IAT,SCITER,ITCOUNT, &
                    GLLHBLCK,BCP,IGUESS,CNVFAC, &
                    TOL, &
                    !   new parameters after inc.p removal
                    naezd, lmaxd, naclsd, xdim, ydim, zdim, &
                    natbld, nthrds)


  implicit none


  integer:: naezd
  integer:: lmaxd
  integer:: naclsd
  integer:: xdim
  integer:: ydim
  integer:: zdim
  integer:: natbld
  integer:: nthrds

  INTEGER :: LMMAXD
  !     PARAMETER         (LMMAXD= (LMAXD+1)**2)
  !     INTEGER            NDIM,NAEZ
  !     PARAMETER         (NAEZ=NAEZD,NDIM=NAEZD*LMMAXD)
  !     INTEGER            NGTBD
  !     PARAMETER         (NGTBD = NACLSD*LMMAXD)
  !     INTEGER            NBLCKD
  !     PARAMETER         (NBLCKD= XDIM * YDIM * ZDIM)
  INTEGER :: NLEN
  INTEGER :: NDIM
  !     PARAMETER         (NLEN=NAEZD*LMMAXD)

  double complex :: CONE
  double complex :: CZERO
  parameter         (CONE = (1.0D0,0.0D0),CZERO = (0.0D0,0.0D0))
  ! ..

  ! global scalars ..
  double precision::CNVFAC
  double precision::TOL
  integer::IAT
  integer::SCITER
  integer::ITCOUNT
  integer::BCP
  integer::IGUESS
  logical::QMRABS
  ! ..

  ! global arrays ..
  !     DOUBLE COMPLEX     TMATLL(LMMAXD,LMMAXD,NAEZD),
  !    +                   X2(LMMAXD*NAEZD,LMMAXD),
  !    +                   GLLH1(LMMAXD,NGTBD,NAEZD),
  !    +                   GLLHBLCK(LMMAXD*NATBLD,LMMAXD*NATBLD*NBLCKD)

  double complex :: TMATLL((LMAXD+1)**2,(LMAXD+1)**2,NAEZD)
  double complex :: X2(NAEZD*(LMAXD+1)**2,(LMAXD+1)**2)
  double complex :: GLLH1((LMAXD+1)**2,NACLSD*(LMAXD+1)**2,NAEZD)
  double complex :: GLLHBLCK(NATBLD*(LMAXD+1)**2, NATBLD*XDIM*YDIM*ZDIM*(LMAXD+1)**2)

  integer:: NUMN0(NAEZD)
  integer:: INDN0(NAEZD,NACLSD)
  ! ..

  ! local scalars ..
  double complex :: ZTMP
  double precision::DTMP
  double precision::RESNAV
  double precision::TOLAV
  integer::NLIM
  integer::INIT
  integer::LM2
  integer::LM1
  integer::IL1
  integer::I
  integer::DIMVEC
  integer::IT
  integer::PROBE
  integer::ITPROBE
  logical::EXITIT
  logical::CNVCONST
  ! ..

  ! local arrays ..

  !     small, local arrays with dimension LMMAXD = (LMAXD+1)**2
  double complex :: RHO((LMAXD+1)**2)
  double complex :: ETA((LMAXD+1)**2)
  double complex :: BETA((LMAXD+1)**2)
  double complex :: ALPHA((LMAXD+1)**2)

  double precision::R0((LMAXD+1)**2)
  double precision::RESN((LMAXD+1)**2)
  double precision::VAR((LMAXD+1)**2)
  double precision::TAU((LMAXD+1)**2)
  double precision::COS((LMAXD+1)**2)
  double precision::N2B((LMAXD+1)**2)
  double precision::TOLB((LMAXD+1)**2)

  integer::ITER((LMAXD+1)**2)
  logical::DONE((LMAXD+1)**2)

  double precision::HSTR(10)
  integer::HSTI(3)

  ! some large local arrays
  !     DOUBLE COMPLEX     B(NDIM,LMMAXD)
  !     DOUBLE COMPLEX     VECS(NDIM,LMMAXD,9)
  !     DOUBLE COMPLEX     DUMMY(LMMAXD*NAEZD,LMMAXD)
  double complex, allocatable, dimension(:,:)   :: B
  double complex, allocatable, dimension(:,:,:) :: VECS
  double complex, allocatable, dimension(:,:)   :: DUMMY

  ! external ..
  external           DZNRM2,ZDOTU,ZRANDN
  double complex :: ZDOTU
  double precision::   DZNRM2

  integer :: memory_stat
  logical :: memory_fail

  !=======================================================================
  ! INITIALIZATION I
  !=======================================================================

  LMMAXD = (LMAXD+1)**2
  NDIM = NAEZD*LMMAXD
  NLEN = NAEZD*LMMAXD

  !     Allocate arrays
  memory_stat = 0
  memory_fail = .false.

  allocate(B(NDIM,LMMAXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(VECS(NDIM,LMMAXD,9), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(DUMMY(LMMAXD*NAEZD,LMMAXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (memory_fail == .true.) then
    write(*,*) "MMINVMOD: FATAL Error, failure to allocate memory."
    write(*,*) "       Probably out of memory."
    stop
  end if

  TOLAV = 0.0
  QMRABS = .true.           ! QMRABS tolerance for residual norm is defined globally
  CNVCONST = .true.         ! CNVFAC is constant for all sc-steps

  if (CNVCONST) then
    CNVFAC  = 1000
  endif

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

  call CINIT(NDIM*LMMAXD,B)

  do LM1=1,LMMAXD
    IL1=LMMAXD*(IAT-1)+LM1
    do LM2=1,LMMAXD
      B(IL1,LM2)=-TMATLL(LM1,LM2,IAT)
    enddo
  enddo

  if (IGUESS == 1 .and. SCITER > 1) then
    do I=1,NAEZD
      do LM1=1,LMMAXD
        IL1=LMMAXD*(I-1)+LM1
        do LM2=1,LMMAXD
          B(IL1,LM2)=-TMATLL(LM1,LM2,I)
        enddo
      enddo
    enddo
  endif

  do I = 1,NDIM
    do LM2=1,LMMAXD
      VECS(I,LM2,2) = - B(I,LM2)
      VECS(I,LM2,1) = CZERO
      do DIMVEC=3,9
        VECS(I,LM2,DIMVEC) = CZERO
      enddo
    enddo
  enddo

  NLIM = 2000
  INIT = 0

  !--------------
  do LM2=1,LMMAXD
  !--------------
    
    ! N2B(LM2)   = DZNRM2(NDIM,B(1,LM2),1)           ! norm of right-hand-sight
    if (QMRABS) then
      TOLB(LM2)  = TOL
    else
      TOLB(LM2)  = MAX(TOL*N2B(LM2),1.0D-10)          ! adapted relative tolerance
    endif
    
    TOLAV      = TOLAV  + TOLB(LM2)/LMMAXD
    
    DONE(LM2) = .false.
    ITER(LM2) = 0
    PROBE     = 1
    
    call ZAXPBY(NLEN,VECS(1,LM2,5), &
                CONE,VECS(1,LM2,2),CZERO,VECS(1,LM2,5))
    call ZAXPBY(NLEN,VECS(1,LM2,1), &
                CZERO,VECS(1,LM2,1),CZERO,VECS(1,LM2,1))
    
  !--------------
  enddo
  !--------------


  !--------------
  do LM2=1,LMMAXD
  !--------------
    
    ! R0(LM2) = DZNRM2(NLEN,VECS(1,LM2,5),1)
    R0(LM2) = N2B(LM2)
    
    ! Check whether the auxiliary vector must be supplied.
    
    if (INIT == 0) call ZRANDN (NLEN,VECS(1,LM2,3),1)
    
    !     Initialize the variables.
    
    RESN(LM2) = 1.0D0
    RHO(LM2)  = CONE
    VAR(LM2)  = 0.0D0
    ETA(LM2)  = CZERO
    TAU(LM2)  = R0(LM2) * R0(LM2)

    call ZAXPBY(NLEN,VECS(1,LM2,8), &
                CZERO,VECS(1,LM2,8),CZERO,VECS(1,LM2,8))

    call ZAXPBY(NLEN,VECS(1,LM2,4), &
                CZERO,VECS(1,LM2,4),CZERO,VECS(1,LM2,4))

    call ZAXPBY(NLEN,VECS(1,LM2,6), &
                CZERO,VECS(1,LM2,6),CZERO,VECS(1,LM2,6))
    
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
    
    if (BCP == 1) then
        
      do LM2 = 1,LMMAXD
        call ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
      !            WRITE(*,*) IT,LM2,DZNRM2(NDIM,DUMMY(1,LM2),1),
      !     +                 BETA(LM2)
      enddo
        
      call APPBLCKCIRC(DUMMY,GLLHBLCK, &
      naezd,lmaxd,nthrds,natbld,xdim,ydim,zdim)
        
    else
        
      do LM2 = 1,LMMAXD
        call ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
      enddo
        
    endif
    
    call SPRSZMM( &
    IAT,GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE, &
    CONE,CZERO, &
    VECS(1,1,9))
    
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
    
    !--------------
    do LM2=1,LMMAXD
      if ( .not. DONE(LM2)) then
        !--------------
            
        call ZAXPBY(NLEN,VECS(1,LM2,6), &
                    CONE,VECS(1,LM2,6),-ALPHA(LM2),VECS(1,LM2,4))

        ZTMP = VAR(LM2) * COS(LM2)

        call ZAXPBY(NLEN,VECS(1,LM2,7), &
                    CONE,VECS(1,LM2,6),ZTMP,VECS(1,LM2,7))
            
      !--------------
      endif
    enddo
    !--------------
    
    EXITIT = .true.
    do LM2 = 1, LMMAXD
      if ( .not. DONE(LM2)) EXITIT = .false.
    enddo
    if (EXITIT) go to 67
    
    
    if (BCP == 1) then
        
      do LM2 = 1,LMMAXD
        call ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
      enddo
        
      call APPBLCKCIRC(DUMMY,GLLHBLCK, &
      naezd,lmaxd,nthrds,natbld,xdim,ydim,zdim)
        
    else
        
      do LM2 = 1,LMMAXD
        call ZCOPY(NDIM,VECS(1,LM2,6),1,DUMMY(1,LM2),1)
      enddo
        
    endif
    
    call SPRSZMM( &
    IAT,GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE, &
    CONE,CZERO, &
    VECS(1,1,8))
    
    !     VECS(:,:,6) input vector to be multiplied by A = GLLH1
    !     VECS(:,:,8) result
    
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
            
    !--------------
      endif
    enddo
    !--------------
    
    
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
            
        call APPBLCKCIRC(DUMMY,GLLHBLCK,naezd,lmaxd,nthrds,natbld,xdim,ydim,zdim)
            
      else
            
        do LM2 = 1,LMMAXD
          call ZCOPY(NDIM,VECS(1,LM2,1),1,DUMMY(1,LM2),1)
        enddo
            
      endif
        
      call SPRSZMM( &
      IAT,GLLH1,NUMN0,INDN0,DUMMY(1,1),DONE, &
      CONE,CZERO, &
      VECS(1,1,9))
        
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
            write(6,*) 'stgn.atm. ',IAT,LM2,' w. ',RESN(LM2),TOLB(LM2)
          endif
        enddo
      endif
      ! ..
        
      PROBE  = MAX(1,INT( LOG( RESNAV/(TOLAV*CNVFAC) ) ) )
        
      EXITIT = .true.
      do LM2 = 1, LMMAXD
        if ( .not. DONE(LM2)) EXITIT = .false.
      enddo
      if (EXITIT) go to 67
        
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
    
     call APPBLCKCIRC(DUMMY,GLLHBLCK,naezd,lmaxd,nthrds,natbld,xdim,ydim,zdim)
    
     do LM2 = 1,LMMAXD
       call ZCOPY(NDIM,DUMMY(1,LM2),1,VECS(1,LM2,1),1)
     enddo
    
   endif

   ! <<<<<<<<<<

   call CINIT(NDIM*LMMAXD,X2)

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

   deallocate(B)
   deallocate(VECS)
   deallocate(DUMMY)

 end subroutine MMINVMOD
