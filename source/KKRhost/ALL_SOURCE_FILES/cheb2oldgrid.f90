subroutine cheb2oldgrid(nrmax, nrmaxnew, lmmaxpot, rmesh, ncheb, npan_tot, &
  rpan_intervall, ipan_intervall, arrayin, arrayout, irmd)

!use mod_cheb, only: getCCmatrix, getCinvmatrix
  implicit none

! Interpolate from Chebyshev mesh to the old radial mesh
! Programmed by David Bauer, 2011-2013
  integer :: ncheb, npan_tot, nrmax, nrmaxnew, lmmaxpot, irmd
  double precision :: rmesh(irmd)
  double precision :: rpan_intervall(0:npan_tot)
  integer :: ipan_intervall(0:npan_tot)
  double complex :: arrayin(nrmaxnew, lmmaxpot)
  double complex :: arrayout(nrmax, lmmaxpot)
!local
  integer :: in, ir, ir2, ir3, it, ilm
  integer :: intsub(2, nrmax)
  double precision :: cinvmatrix(0:ncheb, 0:ncheb)
  double precision, allocatable :: ccmatrix(:, :)
  double complex :: alphaparams(0:ncheb, lmmaxpot)
  double precision :: rmeshnorm(nrmax)
  double precision :: tol
  parameter (tol=1.d-09)

! divide the mesh into subintervals
  intsub(1, :) = 0
  intsub(2, :) = -1
  in = 1

  ir = 1
  it = 1
  do while (ir<=nrmax .and. in<=npan_tot)
    if (abs(rmesh(ir)-rpan_intervall(in))<tol) then
      intsub(it, in) = ir
      intsub(2, in) = ir
      it = 1
      ir = ir + 1
      in = in + 1
    else if (((rmesh(ir)-rpan_intervall(in)<-tol)) .and. ((rmesh(ir)- &
        rpan_intervall(in-1))>-tol)) then
      intsub(it, in) = ir
      intsub(2, in) = ir
      it = 2
      ir = ir + 1
    else if ((rmesh(ir)-rpan_intervall(in-1))<-tol) then
      intsub(1, in) = 0
      intsub(2, in) = -1
      it = 1
      ir = ir + 1
    else
      it = 1
      in = in + 1
    end if
  end do
  if (abs(rmesh(nrmax)-rpan_intervall(npan_tot))>tol) then
    write (*, *) rmesh(nrmax), rpan_intervall(npan_tot), &
      rmesh(nrmax) - rpan_intervall(npan_tot)
    stop 'error with the new and old mesh'
  end if


  call getcinvmatrix(ncheb, cinvmatrix)

  in = 0
  do while (in<=npan_tot)
    in = in + 1
    if (intsub(2,in)<intsub(1,in)) cycle
    alphaparams = 0.0d0
    do ilm = 1, lmmaxpot
      do ir2 = 0, ncheb
        do ir = 0, ncheb
          alphaparams(ir2, ilm) = alphaparams(ir2, ilm) + &
            cinvmatrix(ir2, ir)*arrayin(ipan_intervall(in-1)+1+ir, ilm)
        end do
      end do
    end do

! Transform to normalized coordinates between [-1,1]
! Shift upper & lower end to be centered around zero
    do ir = intsub(1, in), intsub(2, in)
      ir2 = ir + 1 - intsub(1, in)
      rmeshnorm(ir2) = (2d0*rmesh(ir)-(rpan_intervall(in)+rpan_intervall(in- &
        1)))/(rpan_intervall(in)-rpan_intervall(in-1))
      if (abs(rmeshnorm(ir2))>1.0d0) then
        if (abs(abs(rmeshnorm(ir2))-1.0d0)<tol) then
          rmeshnorm(ir2) = sign(1.0d0, rmeshnorm(ir2))
        else
          write (*, *) 'ir, rmeshnorm(ir2)', ir, rmeshnorm(ir2)
          write (*, *) 'rmeshnorm not in interval [-1,1] in cheb2oldgrid'
          write (*, *) 'Check radial mesh (panels) or increase tolerance.'
          stop 'rmeshnorm not in interval [-1,1] in cheb2oldgrid'
        end if
      end if
    end do

    allocate (ccmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

    call getccmatrix(ncheb, rmeshnorm(1:intsub(2,in)-intsub(1, &
      in)+1), intsub(2,in)-intsub(1,in)+1, ccmatrix)
    do ilm = 1, lmmaxpot
      do ir2 = intsub(1, in), intsub(2, in)
        do ir = 0, ncheb
          ir3 = ir2 - intsub(1, in) + 1
          arrayout(ir2, ilm) = arrayout(ir2, ilm) + ccmatrix(ir3, ir)* &
            alphaparams(ir, ilm)
        end do
      end do
    end do !ilm=1,lmmaxpot
    deallocate (ccmatrix) !(CCmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

  end do !in

end subroutine

