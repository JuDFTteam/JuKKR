  subroutine spline_panels(r,y,nr,y2,numpan,numrcut)
! Generalized version of spline, which treats the panels
!!!!!!!!!!
! copied from 'Numerical recipes in Fortran 77', page 109
! subroutine calculates second derivative of input data
! input: data points x(n),y(n)
!       first derivatives y1
! output: second derivative y2
  use global, only: i4b,r8b,c8b,npanat,ircutat, nrpts0
  use derivative_panels

  implicit none
  
  integer(kind=i4b), intent(in)       :: nr
  real(kind=r8b),    intent(in)       :: r(1:nr)
  complex(kind=c8b), intent(in)       :: y(1:nr)
  complex(kind=c8b), intent(out)      :: y2(1:nr)
! --> Number of panels > 1
  integer(kind=i4b), intent(in)       :: numpan, numrcut(numpan+1) 
! #############################################################
  integer(kind=i4b)                   :: ip,ist,ien
  complex(kind=c8b)                   :: y1(1:nr)
  
  call calc_derivative_panels(y,y1,r,nr,numpan,numrcut(:))
! Loop over panels
  do ip = 1, numpan
! Begin and end of panels
    if (ip == 1) then
      ist = numrcut(ip) + 1            ! Start from 1  
    else
      ist = numrcut(ip) - nrpts0(1) + 2                                
    endif
    ien = numrcut(ip+1) - nrpts0(1) + 1                               
  
    call spline_panels2(r,y,ist,ien,y1(ist),y1(ien),y2(ist:ien))
  end do ! panels

  end subroutine spline_panels
