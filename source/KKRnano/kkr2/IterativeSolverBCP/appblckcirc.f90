!> apply block circulant preconditioner on vecs_in and put result into vecs_out
subroutine appblckcirc(vecs_in, vecs_out, gllhblck, naez,lmmaxd, natbld, xdim, ydim, zdim, num_columns)
  implicit none

  include 'fftw3.f'

  integer naez
  integer lmmaxd
  !     number of atoms per preconditioning block
  integer natbld
  !     number of preconditioning blocks in each direction
  integer xdim
  integer ydim
  integer zdim
  integer num_columns


  double complex, parameter :: cone  = (1.d0, 0.d0)
  double complex, parameter :: czero = (0.d0, 0.d0)

  double complex vecs_in(naez*lmmaxd,num_columns)
  double complex vecs_out(naez*lmmaxd,num_columns)
  double complex gllhblck(natbld*lmmaxd,xdim*ydim*zdim*natbld*lmmaxd)

  !     local arrays - large, stack based!!!
  double complex :: tblck(natbld*lmmaxd,natbld*lmmaxd)
  double complex :: txk(natbld*lmmaxd)
  double complex :: tyk(natbld*lmmaxd)

  !     local arrays - large, stack based!!!
  double complex :: x(xdim,ydim,zdim)
  double complex :: xk(natbld*lmmaxd,xdim,ydim,zdim)
  double complex :: y(xdim,ydim,zdim)
  double complex :: yk(natbld*lmmaxd,xdim,ydim,zdim)

!ibm* align(32, xk, yk)

  ! ..
  ! local scalars ..
  double precision :: fac
  integer :: lm1, lmatbl, ix, iy, iz, j, num
  integer(kind=8) :: FFTwplan_fwd, FFTwplan_bwd

  num = lmmaxd*natbld

  fac = 1.d0/dble(xdim*ydim*zdim)

  ! all threads use the same FFTw3-plans
  ! - possible according to FFTw3 documentation
  call dFFTw_plan_dft_3d(FFTwplan_bwd, xdim, ydim, zdim, x, x,  FFTw_backward, FFTw_estimate)
  call dFFTw_plan_dft_3d(FFTwplan_fwd, xdim, ydim, zdim, y, y,  FFTw_forward,  FFTw_estimate)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !$omp parallel private(x,y,tblck,txk,tyk,xk,yk)
  !$omp do
  do lm1 = 1, num_columns

    !=======================================================================

    !-----------------------------------------------------------------------
    ! perform fourier-transform backward              begin
    !-----------------------------------------------------------------------

    do lmatbl = 1, num

      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 1, xdim

            x(ix,iy,iz) = vecs_in(((ix-1)+(iy-1)*xdim+(iz-1)*xdim*ydim)*num+lmatbl,lm1)

          enddo ! ix
        enddo ! iy
      enddo ! iz

      call dFFTw_execute_dft(FFTwplan_bwd, x, x)

      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 1, xdim
            xk(lmatbl,ix,iy,iz) = x(ix,iy,iz)*fac
          enddo ! ix
        enddo ! iy
      enddo ! iz

    enddo ! lmatbl

    !-----------------------------------------------------------------------
    ! perform fourier-transform backward              end
    !-----------------------------------------------------------------------

    do iz = 1, zdim
      do iy = 1, ydim
        do ix = 1, xdim

          txk(1:num) = xk(1:num,ix,iy,iz)

          do j = 1, num
            tblck(1:num,j) = gllhblck(1:num,num*((ix-1)+(iy-1)*xdim+(iz-1)*xdim*ydim)+j)
          enddo ! j

          call zgemv('n',num, num, cone, tblck, num, txk, 1, czero, tyk, 1)

          yk(1:num,ix,iy,iz) = tyk(1:num)

        enddo ! ix
      enddo ! iy
    enddo ! iz

    !-----------------------------------------------------------------------
    ! perform fourier-transform forward               begin
    !-----------------------------------------------------------------------

    do lmatbl = 1, num

      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 1, xdim
            y(ix,iy,iz) = yk(lmatbl,ix,iy,iz)
          enddo ! ix
        enddo ! iy
      enddo ! iz

      call dFFTw_execute_dft(FFTwplan_fwd, y, y)

      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 1, xdim

            vecs_out(((ix-1)+(iy-1)*xdim+(iz-1)*xdim*ydim)*num+lmatbl,lm1) = y(ix,iy,iz)

          enddo ! ix
        enddo ! iy
      enddo ! iz

    enddo ! lmatbl

  !-----------------------------------------------------------------------
  ! perform fourier-transform forward               end
  !-----------------------------------------------------------------------

  enddo ! lm1
  !$omp enddo
  !$omp endparallel
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  call dFFTw_destroy_plan(FFTwplan_fwd)
  call dFFTw_destroy_plan(FFTwplan_bwd)

endsubroutine ! appblckcirc
