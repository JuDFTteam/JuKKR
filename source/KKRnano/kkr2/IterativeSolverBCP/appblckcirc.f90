subroutine APPBLCKCIRC(VECS,GLLHBLCK, &
                       naez,lmmaxd, &
                       natbld, xdim, ydim, zdim, num_columns)

  implicit none

  INCLUDE 'fftw3.f'

  integer naez
  integer lmmaxd
  !     number of atoms per preconditioning block
  integer natbld
  !     number of preconditioning blocks in each direction
  integer xdim
  integer ydim
  integer zdim
  integer num_columns


  double complex CONE
  parameter      (CONE  = ( 1.0D0,0.0D0))
  double complex CZERO
  parameter      (CZERO = ( 0.0D0,0.0D0))

  double complex VECS(NAEZ*LMMAXD,num_columns)
  double complex GLLHBLCK(NATBLD*LMMAXD, &
  XDIM*YDIM*ZDIM*NATBLD*LMMAXD)

  !     local arrays - large, stack based!!!
  double complex :: TBLCK(NATBLD*LMMAXD,NATBLD*LMMAXD)
  double complex :: TXK(NATBLD*LMMAXD)
  double complex :: TYK(NATBLD*LMMAXD)

  !     local arrays - large, stack based!!!
  double complex :: X(XDIM,YDIM,ZDIM)
  double complex :: XK(NATBLD*LMMAXD,XDIM,YDIM,ZDIM)
  double complex :: Y(XDIM,YDIM,ZDIM)
  double complex :: YK(NATBLD*LMMAXD,XDIM,YDIM,ZDIM)

!IBM* ALIGN(32, XK, YK)

  ! ..
  ! local scalars ..
  double precision FAC
  integer        LM1,LMATBL,IX,IY,IZ,I,J
  integer(8)     FFTWPLAN_fwd
  integer(8)     FFTWPLAN_bwd
  integer        num

  num = lmmaxd*natbld

  FAC = 1/FLOAT(XDIM*YDIM*ZDIM)

  ! all threads use the same fftw3-plans
  ! - possible according to fftw3 documentation
  call DFFTW_PLAN_DFT_3D(FFTWPLAN_bwd,XDIM,YDIM,ZDIM,X,X, &
  FFTW_BACKWARD,FFTW_ESTIMATE)
  call DFFTW_PLAN_DFT_3D(FFTWPLAN_fwd,XDIM,YDIM,ZDIM,Y,Y, &
  FFTW_FORWARD,FFTW_ESTIMATE)


  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  !$OMP PARALLEL PRIVATE(X,Y,TBLCK,TXK,TYK,XK,YK)
  !$omp do
  do LM1=1,num_columns

    !=======================================================================

    !-----------------------------------------------------------------------
    ! perform Fourier-transform backward              begin
    !-----------------------------------------------------------------------

    do LMATBL = 1,num

        do IZ =1,ZDIM
          do IY =1,YDIM
            do IX =1,XDIM

              X(IX,IY,IZ) = VECS(((IX-1)+(IY-1)*XDIM+ &
              (IZ-1)*XDIM*YDIM)*num+LMATBL,LM1)

            enddo
          enddo
        enddo

        call DFFTW_EXECUTE_DFT(FFTWPLAN_bwd,X,X)

        do IZ =1,ZDIM
          do IY =1,YDIM
            do IX =1,XDIM
              XK(LMATBL,IX,IY,IZ) = X(IX,IY,IZ)*FAC
            enddo
          enddo
        enddo

    enddo

    !-----------------------------------------------------------------------
    ! perform Fourier-transform backward              end
    !-----------------------------------------------------------------------

    do IZ =1,ZDIM
        do IY =1,YDIM
          do IX =1,XDIM

            do I=1,num
              TXK(I) = XK(I,IX,IY,IZ)
            end do

            do J=1,num
              do I=1,num
                TBLCK(I,J) = GLLHBLCK(I,num* &
                ((IX-1)+(IY-1)*XDIM+(IZ-1)*XDIM*YDIM)+J)
              enddo
            enddo

            call ZGEMV('N',num,num, &
            CONE,TBLCK, &
            num,TXK, &
            1,CZERO,TYK,1)

            do I=1,num
              YK(I,IX,IY,IZ) = TYK(I)
            enddo

          enddo
        enddo
    enddo

    !-----------------------------------------------------------------------
    ! perform Fourier-transform forward               begin
    !-----------------------------------------------------------------------

    do LMATBL = 1,num

        do IZ =1,ZDIM
          do IY =1,YDIM
            do IX =1,XDIM
              Y(IX,IY,IZ) = YK(LMATBL,IX,IY,IZ)
            enddo
          enddo
        enddo

        call DFFTW_EXECUTE_DFT(FFTWPLAN_fwd,Y,Y)

        do IZ =1,ZDIM
          do IY =1,YDIM
            do IX =1,XDIM

              VECS(((IX-1)+(IY-1)*XDIM+ &
              (IZ-1)*XDIM*YDIM)*num+LMATBL,LM1) = Y(IX,IY,IZ)

            enddo
          enddo
        enddo

    enddo ! LM1

  !-----------------------------------------------------------------------
  ! perform Fourier-transform forward               end
  !-----------------------------------------------------------------------

  enddo
  !$omp end do
  !$OMP END PARALLEL
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call DFFTW_DESTROY_PLAN(FFTWPLAN_fwd)
  call DFFTW_DESTROY_PLAN(FFTWPLAN_bwd)

end
