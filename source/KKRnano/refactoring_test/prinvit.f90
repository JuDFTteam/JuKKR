    subroutine PRINVIT( &
    INFO,IAT, &
    NUMN0,INDN0, &
    TMATLL,GLLH1,X0,PRSC,SPRS, &
    GLLKE1)

    use inc_p_wrapper_module

    implicit none
    ! ------------------------------------------------------------------------
    ! Initial guess optimization: Using the result of previous
    ! self-consistency step as starting vector for iterative solver
    ! Options:
    ! a) sc prec
    !    store solution obtained at the last self-consistency iteration
    !                                                         A. Thiess Nov'09
    ! ------------------------------------------------------------------------
    !     .. parameters ..


    integer::          LMMAXD
    parameter        (LMMAXD= (LMAXD+1)**2)
    integer::          ALM
    parameter        (ALM   = NAEZD*LMMAXD)
    integer::          NGTBD
    parameter        (NGTBD = NACLSD*LMMAXD)
    real::             CUT
    parameter        (CUT   = 1.0D-5)
    double complex :: CZERO
    parameter        (CZERO = ( 0.0D0,0.0D0 ))
    double complex :: CONE
    parameter        (CONE  = ( 1.0D0,0.0D0 ))
    !     ..
    !     .. SCALAR ARGUMENTS ..
    character :: INFO
    integer::    IAT
    !     ..
    !     .. ARRAY ARGUMENTS ..
    double complex :: GLLH1(LMMAXD,NGTBD,NAEZD)
    double complex :: TMATLL(LMMAXD,LMMAXD,NAEZD)
    double complex :: TMATP(ALM,LMMAXD)
    double complex :: GLLKE1(ALM,LMMAXD)
    complex::         PRSC(NGUESSD*LMMAXD)
    integer::         NUMN0(NAEZD)
    integer::         INDN0(NAEZD,NACLSD)
    integer::         SPRS(NGUESSD*LMMAXD+1)

    !     .. LOCAL SCALARS ..
    integer::I1
    integer::LM1
    integer::LM2
    integer::JLM
    integer::IL1
    integer::JSP
    integer::ISP
    logical::LSAME
    !     ..
    !     .. LOCAL ARRAYS ..
    double complex :: X0(ALM,LMMAXD)
    logical        :: DONE(LMMAXD)
!     ..
!     .. EXTERNAL FUNCTIONS ..
!   LSAME is a routine from LAPACK
    external         LSAME
!     ..
!     ..
!-----------------------------------------------------------------------

! ================================================================
    if (LSAME(INFO,'I')) then
    ! ================================================================
    ! Ia)  translate sparse format of PRSC back to X0 ..
    
        call CINIT(ALM*LMMAXD,X0)
    
        do JSP = 1, NGUESSD*LMMAXD
        
            if (SPRS(JSP) == (NAEZD*LMMAXD*LMMAXD+9999)) &
            goto 99
        
            JLM = INT((SPRS(JSP)-1)/LMMAXD) + 1
            LM2 = MOD((SPRS(JSP)-1),LMMAXD) + 1
        
            X0(JLM,LM2) = &
            DCMPLX(REAL(PRSC(JSP)),AIMAG(PRSC(JSP)))
        
        enddo
    
        99 continue
    ! ..
    ! ================================================================

    ! ================================================================
    ! Ib)  obtain b' = b - Ax
    !                        0            ..
    

    !       The DONE-array keeps track of which lm-component has
    !       already converged - initialised here
         
        do LM1=1,LMMAXD
            DONE(LM1) = .false.
        enddo
    
    !--?------  initialize TMATLL for the next IE by GLLKE0 ??? ---------
    ! sparse matrix multiplication: (-1) * GLLH1 * X0 = TMATP

        call SPRSZMM( &
        IAT,GLLH1,NUMN0,INDN0,X0,DONE, &
        -CONE,CZERO, &
        TMATP)

    !----------   \Delta t' = \Delta t - X0 + \Delta t * G_ref * X0------
    
        do I1=1,NAEZD
            do LM1=1,LMMAXD
                IL1=LMMAXD*(I1-1)+LM1
                do LM2=1,LMMAXD

                    if (I1 == IAT) then
                        TMATLL(LM1,LM2,I1) = TMATP(IL1,LM2) + TMATLL(LM1,LM2,I1)
                    else
                        TMATLL(LM1,LM2,I1) = TMATP(IL1,LM2)
                    endif

                enddo
            enddo
        enddo

    endif
! ..
! ================================================================



    if (LSAME(INFO,'F')) then
    ! ================================================================
    ! Fa) determine true solution by adding initial guess ..
    
        do I1=1,NAEZD
            do LM1=1,LMMAXD
                IL1=LMMAXD*(I1-1)+LM1
                do LM2=1,LMMAXD
                    GLLKE1(IL1,LM2) = GLLKE1(IL1,LM2) + X0(IL1,LM2)
                enddo
            enddo
        enddo
    ! ..
    ! ================================================================

    ! ================================================================
    ! Fb) store new result as initial guess for the next self-consistency
    !     iteration in sparse format ..
    
        ISP = 1
    
        do I1=1,NAEZD
            do LM1=1,LMMAXD
                JLM=LMMAXD*(I1-1)+LM1
                do LM2=1,LMMAXD

                ! sparse >>
                    if ((ABS(DREAL(GLLKE1(JLM,LM2))) > CUT) .or. &
                    (ABS(DIMAG(GLLKE1(JLM,LM2))) > CUT)) then
                    
                        SPRS(ISP)=  LMMAXD*(JLM-1) + LM2
                    
                    !             Convert to single precision
                        PRSC(ISP) = &
                        CMPLX(DREAL(GLLKE1(JLM,LM2)),DIMAG(GLLKE1(JLM,LM2)))
                    
                        ISP  =  ISP + 1
                    
                    endif
                ! sparse <<

                enddo
            enddo
        enddo
    
        SPRS(ISP) = NAEZD*LMMAXD*LMMAXD + 9999 !STOP signature
    ! ..
    ! ================================================================
    
    endif

    return

    end subroutine PRINVIT
