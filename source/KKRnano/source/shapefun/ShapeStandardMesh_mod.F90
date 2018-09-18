!> Module for interstitial mesh generation given critical points
!> (= panel positions).

module ShapeStandardMesh_mod
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
#include "macros.h"
  implicit none
  private
  public :: mesh

  contains

  !----------------------------------------------------------------------------
  !> radial mesh generation for non-spherical part.
  !> @brief
  !>
  !-----------------------------------------------------------------------
  !>     this routine defines a unique suitable radial mesh 'xrn,drn' of
  !>     'meshn' points,distributed into 'npan' pannels  defined  by the
  !>     critical points 'crt'
  !-----------------------------------------------------------------------
  subroutine mesh(crt, npan, nm, xrn, drn, meshn, npoi, keypan, nmin, meshnd, npand, o)
    double precision, intent(inout) :: crt(npand) ! reordered on exit
    integer, intent(in) :: npan !! number of critical points
    integer, intent(inout) :: nm(npand) ! number of points in panels/critical intervals, if keypan /= 0, there could be an intent(in) here
    double precision, intent(out) :: xrn(meshnd) !! radial mesh points
    double precision, intent(out) :: drn(meshnd) !! weights for the radial mesh points
    integer, intent(out) :: meshn !! number of mesh points generated
    integer, intent(in) :: npoi !! number of mesh points requested
    integer, intent(in) :: keypan !! set to 0 to generate mesh
    integer, intent(in) :: nmin !! minimal number of points for panels/critical intervals
    integer, intent(in) :: meshnd !! dimension of crt, xrn and drn arrays
    integer, intent(in) :: npand !! dimension of nm array
    integer, intent(in) :: o !! output unit

    integer :: iord, ipan, ir_start, ir, ir_end, ir_off
    double precision :: crt_temp, off, step

    !-----------------------------------------------------------------------

    ! sort critical points with a stupid algorithm
    do iord = 1, npan
      crt_temp = crt(iord)
      do ipan = npan, iord, -1
        if (crt(ipan) <= crt_temp) then
          crt_temp  = crt(ipan)
          crt(ipan) = crt(iord)
          crt(iord) = crt_temp
        endif ! unordered
      enddo ! ipan
    enddo ! iord

    if (keypan == 0) call mesh0(crt, npan, npoi, nmin, nm)

    if (o>0) write(o, fmt="(/50('-')/'i',13x,'suitable radial mesh',15x,'i'/'i',13x,20('*'),15x,'i'/'i',3x,'ipan',7x,'from',7x,'to',13x,'points  i'/'i',48x,'i')")

    ir_off = 0
    do ipan = 1, npan-1

      if (o>0) write(o, fmt="('i',2x,i5,2e14.7,i10,'   i')") ipan, crt(ipan:ipan+1), nm(ipan)

      ir_start = ir_off + 1 ! ir start index
      ir_end   = ir_off + nm(ipan) ! prelim. end index
      if (meshnd < ir_end) die_here("nxr ="+meshnd+"is too small!")
      
      step = (crt(ipan+1) - crt(ipan))/dble(ir_end-ir_start)
      off = crt(ipan) - step*dble(ir_start)
      do ir = ir_start, ir_end
        xrn(ir) = step*ir + off
        drn(ir) = step
      enddo ! ir
      ir_off = ir_off + nm(ipan)
    enddo ! ipan
    meshn = ir_off

    if (o>0) write(o,fmt="(50('-'))")

  endsubroutine ! mesh

  subroutine mesh0(crt, npan, naprox, nmin, nm)
    ! ***********************************************************
    ! *  this subroutine calculates an appropriate mesh for
    ! *  the shape functions. more than nmin points between two
    ! *  critical points
    ! *  in case of more dense mesh increase nmin
    ! ***********************************************************
    double precision, intent(in) :: crt(:)
    integer, intent(in) :: npan
    integer, intent(in) :: naprox
    integer, intent(in) :: nmin
    integer, intent(out) :: nm(:)
    
    double precision :: dist, d1
    integer :: n, ntot, i, na

    if ((npan-1)*nmin > naprox) die_here("npan ="+npan+"  nmin ="+nmin+"  naprox ="+naprox-", increase number of points!")

    dist = abs(crt(1) - crt(npan))
    do i = 1, npan-1
      d1 = abs(crt(i) - crt(i+1))
      nm(i) = floor((naprox*d1)/dist)
      if (nm(i) < nmin) nm(i) = nmin
    enddo ! i (panels)
    n = nmin*(npan-1)
    n = naprox-n
    if (n <= 0) die_here("increase number of mesh points! n ="+n)
    d1 = n/dble(naprox)
    ntot = n
    do i = 1, npan-1
      na = nint(d1*nm(i))
      if (nm(i) > nmin .and. (ntot-na) > 0) then
        nm(i) = nmin + na
        ntot = ntot - na
      endif
    enddo ! i
    nm(1) = nm(1) + ntot
    
  endsubroutine ! mesh0

endmodule ! ShapeStandardMesh_mod
