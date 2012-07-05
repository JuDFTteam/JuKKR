    subroutine initialGuess_start(IAT,NUMN0,INDN0, &
                       TMATLL,GLLH1,X0,PRSC, &
                       !GLLKE1, &
                       ! new parameters after inc.p removal
                       naez, lmmaxd, naclsd, nguessd, nthrds)

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
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: naclsd
    integer, intent(in) :: nguessd
    integer, intent(in) :: nthrds

    double complex :: CZERO
    parameter        (CZERO = ( 0.0D0,0.0D0 ))
    double complex :: CONE
    parameter        (CONE  = ( 1.0D0,0.0D0 ))

    !     .. SCALAR ARGUMENTS ..
    integer::    IAT
    !     ..
    !     .. ARRAY ARGUMENTS ..
    !double complex :: GLLH1(LMMAXD,NGTBD,NAEZD)
    double complex :: GLLH1(LMMAXD,NACLSD*LMMAXD,NAEZ)
    double complex :: TMATLL(LMMAXD,LMMAXD,NAEZ)

    complex::         PRSC(NGUESSD*LMMAXD)
    integer::         NUMN0(NAEZ)
    integer::         INDN0(NAEZ,NACLSD)

    double complex :: X0(NAEZ*LMMAXD,LMMAXD)

    !     .. LOCAL SCALARS ..
    integer::site_index
    integer::LM1
    integer::LM2
    integer::site_lm_index
    integer::site_lm_index2
    integer::ind
    !     .. LOCAL ARRAYS ..
    logical        :: DONE(LMMAXD)
    double complex, allocatable, dimension(:,:) :: TMATP
    !double complex :: TMATP(NAEZ*LMMAXD,LMMAXD)     ! LARGE, local !!!

    integer::          ALM
    integer::          NAEZD

    integer :: memory_stat

    ALM   = NAEZ*LMMAXD
    NAEZD = NAEZ

!-----------------------------------------------------------------------

    ! Ia)  translate sparse format of PRSC back to X0 ..

        ! allocate array
        allocate(TMATP(NAEZ*LMMAXD,LMMAXD), stat=memory_stat)

        if (memory_stat /= 0) then
          write(*,*) "initialGuess: FATAL Error, failure to allocate memory."
          write(*,*) "       Probably out of memory."
          stop
        end if

        X0 = CZERO

        do ind = 1, NGUESSD*LMMAXD

            site_lm_index = INT((ind-1)/LMMAXD) + 1
            LM2 = MOD((ind-1),LMMAXD) + 1

            X0(site_lm_index,LM2) = &
            DCMPLX(REAL(PRSC(ind)),AIMAG(PRSC(ind)))

        enddo
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

    ! deallocate work array
    deallocate(TMATP)

    end subroutine initialGuess_start

! =====================================================================================================================
    subroutine initialGuess_finish(X0,PRSC, GLLKE1, &
                       naez, lmmaxd, nguessd)

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
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: nguessd

    double complex :: CZERO
    parameter        (CZERO = ( 0.0D0,0.0D0 ))
    double complex :: CONE
    parameter        (CONE  = ( 1.0D0,0.0D0 ))

    !     ..
    !     .. ARRAY ARGUMENTS ..

    double complex :: GLLKE1(NAEZ*LMMAXD,LMMAXD)
    complex::         PRSC(NGUESSD*LMMAXD)

    double complex :: X0(NAEZ*LMMAXD,LMMAXD)

    !     .. LOCAL SCALARS ..
    integer::ind
    integer::LM1
    integer::LM2
    integer::site_lm_index
    integer::site_lm_index2
    integer::site_index

    integer::          ALM
    integer::          NAEZD

    ALM   = NAEZ*LMMAXD
    NAEZD = NAEZ

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
    !     iteration in

        ind = 1

        do site_index=1,NAEZD
            do LM1=1,LMMAXD
                site_lm_index=LMMAXD*(site_index-1)+LM1  ! use a combined site and lm-index
                do LM2=1,LMMAXD

                    !             Convert to single precision
                        PRSC(ind) = &
                        CMPLX(DREAL(GLLKE1(site_lm_index,LM2)),DIMAG(GLLKE1(site_lm_index,LM2)))

                        ind =  ind + 1

                enddo
            enddo
        enddo
    ! ..
    ! ================================================================

    end subroutine initialGuess_finish
