    subroutine initialGuess_start(X0,PRSC, &
                                  naez, lmmaxd, nguessd)

    implicit none
    ! ------------------------------------------------------------------------
    ! Initial guess optimization: Using the result of previous
    ! self-consistency step as starting vector for iterative solver
    ! Options:
    ! a) sc prec
    !    store solution obtained at the last self-consistency iteration
    !                                                         A. Thiess Nov'09
    ! simplified: Elias Rabel, 2012
    ! ------------------------------------------------------------------------

    integer, intent(in) :: naez
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: nguessd

    double complex :: CZERO
    parameter        (CZERO = ( 0.0D0,0.0D0 ))
    double complex :: CONE
    parameter        (CONE  = ( 1.0D0,0.0D0 ))

    !     .. ARRAY ARGUMENTS ..
    complex::         PRSC(NGUESSD*LMMAXD)

    double complex :: X0(NAEZ*LMMAXD,LMMAXD)

    !     .. LOCAL SCALARS ..
    integer::LM2
    integer::site_lm_index

    integer::ind

!-----------------------------------------------------------------------

    ! translate sparse format of PRSC back to X0 ..

        X0 = CZERO

        do ind = 1, NGUESSD*LMMAXD

            site_lm_index = INT((ind-1)/LMMAXD) + 1
            LM2 = MOD((ind-1),LMMAXD) + 1

            X0(site_lm_index,LM2) = &
            DCMPLX(REAL(PRSC(ind)),AIMAG(PRSC(ind)))

        enddo

    end subroutine initialGuess_start

! =====================================================================================================================
    subroutine initialGuess_finish(PRSC, GLLKE1, &
                       naez, lmmaxd, nguessd)

    implicit none
    ! ------------------------------------------------------------------------
    ! Initial guess optimization: Using the result of previous
    ! self-consistency step as starting vector for iterative solver
    ! Options:
    ! a) sc prec
    !    store solution obtained at the last self-consistency iteration
    !                                                         A. Thiess Nov'09
    ! simplified: Elias Rabel, 2012
    ! ------------------------------------------------------------------------

    integer, intent(in) :: naez
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: nguessd

    !     .. ARRAY ARGUMENTS ..

    double complex :: GLLKE1(NAEZ*LMMAXD,LMMAXD)
    complex::         PRSC(NGUESSD*LMMAXD)

    !     .. LOCAL SCALARS ..
    integer::ind
    integer::LM1
    integer::LM2
    integer::site_lm_index
    integer::site_index

    ! ================================================================
    ! Fb) store new result as initial guess for the next self-consistency
    !     iteration in

        ind = 1

        do site_index=1,NAEZ
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
