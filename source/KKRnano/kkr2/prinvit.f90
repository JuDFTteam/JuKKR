    subroutine initialGuess_start(IAT,NUMN0,INDN0, &
                       TMATLL,GLLH1,X0,PRSC,SPRS, &
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
    integer::    IAT
    !     ..
    !     .. ARRAY ARGUMENTS ..
    !double complex :: GLLH1(LMMAXD,NGTBD,NAEZD)
    double complex :: GLLH1(LMMAXD,NACLSD*LMMAXD,NAEZ)
    double complex :: TMATLL(LMMAXD,LMMAXD,NAEZ)

    complex::         PRSC(NGUESSD*LMMAXD)
    integer::         NUMN0(NAEZ)
    integer::         INDN0(NAEZ,NACLSD)
    integer::         SPRS(NGUESSD*LMMAXD+1)
    double complex :: X0(NAEZ*LMMAXD,LMMAXD)

    !     .. LOCAL SCALARS ..
    integer::site_index
    integer::LM1
    integer::LM2
    integer::site_lm_index
    integer::site_lm_index2
    integer::sparse_index
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
          write(*,*) "KKRMAT01: FATAL Error, failure to allocate memory."
          write(*,*) "       Probably out of memory."
          stop
        end if

        call CINIT(ALM*LMMAXD,X0)

        do sparse_index = 1, NGUESSD*LMMAXD

            if (SPRS(sparse_index) == (NAEZD*LMMAXD*LMMAXD+9999)) &
            goto 99

            site_lm_index = INT((SPRS(sparse_index)-1)/LMMAXD) + 1
            LM2 = MOD((SPRS(sparse_index)-1),LMMAXD) + 1

            ! TODO: debug - remove
            if (SPRS(sparse_index) < 1 .or. SPRS(sparse_index) > NAEZD*LMMAXD*LMMAXD) then
              write(*,*) "illegal value for SPRS: index value", sparse_index, SPRS(sparse_index)
              stop
            end if

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

    ! deallocate work array
    deallocate(TMATP)

    end subroutine initialGuess_start

! =====================================================================================================================
    subroutine initialGuess_finish(X0,PRSC,SPRS, GLLKE1, &
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

    !integer::          LMMAXD
    !parameter        (LMMAXD= (LMAXD+1)**2)
    !integer::          ALM
    !parameter        (ALM   = NAEZD*LMMAXD) ! NAEZ*(LMAX+1)**2

    real::            CUT
    parameter        (CUT   = 1.0D-5)
    double complex :: CZERO
    parameter        (CZERO = ( 0.0D0,0.0D0 ))
    double complex :: CONE
    parameter        (CONE  = ( 1.0D0,0.0D0 ))

    !     ..
    !     .. ARRAY ARGUMENTS ..

    double complex :: GLLKE1(NAEZ*LMMAXD,LMMAXD)
    complex::         PRSC(NGUESSD*LMMAXD)
    integer::         SPRS(NGUESSD*LMMAXD+1)
    double complex :: X0(NAEZ*LMMAXD,LMMAXD)

    !     .. LOCAL SCALARS ..
    integer::site_index
    integer::LM1
    integer::LM2
    integer::site_lm_index
    integer::site_lm_index2
    integer::sparse_index

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

    end subroutine initialGuess_finish
