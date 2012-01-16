    subroutine PRINVIT(INFO,IAT, &
                       NUMN0,INDN0, &
                       TMATLL,GLLH1,X0,PRSC,SPRS, &
                       GLLKE1, &
                       ! new parameters after inc.p removal
                       naez, lmax, naclsd, nguessd, nthrds)

    !use inc_p_wrapper_module

    implicit none
    ! ------------------------------------------------------------------------
    ! Initial guess optimization: Using the result of previous
    ! self-consistency step as starting vector for iterative solver
    ! Options:
    ! a) sc prec
    !    store solution obtained at the last self-consistency iteration
    !                                                         A. Thiess Nov'09
    ! ------------------------------------------------------------------------

    integer, intent(in) :: naez
    integer, intent(in) :: lmax
    integer, intent(in) :: naclsd
    integer, intent(in) :: nguessd
    integer, intent(in) :: nthrds

    !integer::          LMMAXD
    !parameter        (LMMAXD= (LMAXD+1)**2)
    !integer::          ALM
    !parameter        (ALM   = NAEZD*LMMAXD) ! NAEZ*(LMAX+1)**2
    !integer::          NGTBD
    !parameter        (NGTBD = NACLSD*LMMAXD) ! NACLSD*(LMAX+1)**2

    real::            CUT
    parameter        (CUT   = 1.0D-5)
    double complex :: CZERO
    parameter        (CZERO = ( 0.0D0,0.0D0 ))
    double complex :: CONE
    parameter        (CONE  = ( 1.0D0,0.0D0 ))

    !     .. SCALAR ARGUMENTS ..
    character :: INFO
    integer::    IAT
    !     ..
    !     .. ARRAY ARGUMENTS ..
    !double complex :: GLLH1(LMMAXD,NGTBD,NAEZD)
    double complex :: GLLH1((LMAX+1)**2,NACLSD*(LMAX+1)**2,NAEZ)
    double complex :: TMATLL((LMAX+1)**2,(LMAX+1)**2,NAEZ)
    double complex :: TMATP(NAEZ*(LMAX+1)**2,(LMAX+1)**2)
    double complex :: GLLKE1(NAEZ*(LMAX+1)**2,(LMAX+1)**2)
    complex::         PRSC(NGUESSD*(LMAX+1)**2)
    integer::         NUMN0(NAEZ)
    integer::         INDN0(NAEZ,NACLSD)
    integer::         SPRS(NGUESSD*(LMAX+1)**2+1)
    double complex :: X0(NAEZ*(LMAX+1)**2,(LMAX+1)**2)

    !     .. LOCAL SCALARS ..
    integer::site_index
    integer::LM1
    integer::LM2
    integer::site_lm_index
    integer::site_lm_index2
    integer::sparse_index
    logical::LSAME
    !     ..
    !     .. LOCAL ARRAYS ..
    logical        :: DONE((LMAX+1)**2)
!     ..
!     .. EXTERNAL FUNCTIONS ..
!   LSAME is a routine from LAPACK
    external         LSAME

    integer::          LMMAXD
    integer::          ALM
    integer::          NAEZD

    LMMAXD= (LMAX+1)**2
    ALM   = NAEZ*LMMAXD
    NAEZD = NAEZ

!-----------------------------------------------------------------------

! ================================================================
    if (LSAME(INFO,'I')) then
    ! ================================================================
    ! Ia)  translate sparse format of PRSC back to X0 ..
    
        call CINIT(ALM*LMMAXD,X0)
    
        do sparse_index = 1, NGUESSD*LMMAXD
        
            if (SPRS(sparse_index) == (NAEZD*LMMAXD*LMMAXD+9999)) &
            goto 99
        
            site_lm_index = INT((SPRS(sparse_index)-1)/LMMAXD) + 1
            LM2 = MOD((SPRS(sparse_index)-1),LMMAXD) + 1
        
            X0(site_lm_index,LM2) = &
            DCMPLX(REAL(PRSC(sparse_index)),AIMAG(PRSC(sparse_index)))
        
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

        call SPRSZMM( IAT,GLLH1,NUMN0,INDN0,X0,DONE, &
                      -CONE,CZERO, &
                      TMATP, &
                      naez, lmmaxd, naclsd, nthrds)

    !----------   \Delta t' = \Delta t - X0 + \Delta t * G_ref * X0------
    
        do site_index=1,NAEZD
            do LM1=1,LMMAXD
                site_lm_index2=LMMAXD*(site_index-1)+LM1
                do LM2=1,LMMAXD

                    if (site_index == IAT) then
                        TMATLL(LM1,LM2,site_index) = TMATP(site_lm_index2,LM2) + TMATLL(LM1,LM2,site_index)
                    else
                        TMATLL(LM1,LM2,site_index) = TMATP(site_lm_index2,LM2)
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
    
        do site_index=1,NAEZD
            do LM1=1,LMMAXD
                site_lm_index2=LMMAXD*(site_index-1)+LM1
                do LM2=1,LMMAXD
                    GLLKE1(site_lm_index2,LM2) = GLLKE1(site_lm_index2,LM2) + X0(site_lm_index2,LM2)
                enddo
            enddo
        enddo
    ! ..
    ! ================================================================

    ! ================================================================
    ! Fb) store new result as initial guess for the next self-consistency
    !     iteration in sparse format ..
    
        sparse_index = 1
    
        do site_index=1,NAEZD
            do LM1=1,LMMAXD
                site_lm_index=LMMAXD*(site_index-1)+LM1  ! use a combined site and lm-index
                do LM2=1,LMMAXD

                ! sparse >>
                    if ((ABS(DREAL(GLLKE1(site_lm_index,LM2))) > CUT) .or. &
                    (ABS(DIMAG(GLLKE1(site_lm_index,LM2))) > CUT)) then
                    
                        SPRS(sparse_index)=  LMMAXD*(site_lm_index-1) + LM2
                    
                    !             Convert to single precision
                        PRSC(sparse_index) = &
                        CMPLX(DREAL(GLLKE1(site_lm_index,LM2)),DIMAG(GLLKE1(site_lm_index,LM2)))
                    
                        sparse_index  =  sparse_index + 1
                    
                    endif
                ! sparse <<

                enddo
            enddo
        enddo
    
        SPRS(sparse_index) = NAEZD*LMMAXD*LMMAXD + 9999 !STOP signature
    ! ..
    ! ================================================================
    
    endif

    return

    end subroutine PRINVIT
