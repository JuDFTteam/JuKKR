  module mod_derivative_panels 

  interface calc_derivative_panels
    module procedure calc_derivative_panels_real, calc_derivative_panels_complex
  end interface calc_derivative_panels

  contains

    subroutine calc_derivative_panels_complex(f,dfdr,r,nr,numpan,numrcut)
!   calculates the radial derivative of a function f 
!   takes care of the panels, if SOC from host
    use global, only: i4b, r8b, c8b, nrpts0

!   radial function 
    complex(kind=c8b), intent(in)    :: f(1:nr)
!   derivative of function
    complex(kind=c8b), intent(inout)   :: dfdr(1:nr)
!   --> Number of panels > 1
    integer(kind=i4b), intent(in)    :: numpan, numrcut(numpan+1)
!   rmesh
    integer(kind=i4b), intent(in)    :: nr
    real(kind=r8b), intent(in)       :: r(1:nr) 
!   ------------------------------------------------------------------------
    integer(kind=i4b)                :: ip, ist, ien
    real(kind=r8b)                   :: dri,dri1


    dfdr = 0.d0
!   Loop over panels
    do ip = 1, numpan
!   Begin and end of panels
      if (ip == 1) then
        ist = numrcut(ip) + 1            ! Start from 1  
      else
        ist = numrcut(ip) - nrpts0(1) + 2                                
      endif
      ien = numrcut(ip+1) - nrpts0(1) + 1                               
!     forward difference
      dfdr(ist)= (f(ist+1)-f(ist))/(r(ist+1)-r(ist))
!     central differences with a second order approach (for non uniform r grid)
      do ir=ist+1, ien-1
        dri=r(ir+1)-r(ir)
        dri1=r(ir)-r(ir-1)
        dfdr(ir)=-dri/dri1/(dri+dri1)*f(ir-1)
        dfdr(ir)=dfdr(ir)+(dri-dri1)/dri/dri1*f(ir)
        dfdr(ir)=dfdr(ir)+dri1/dri/(dri+dri1)*f(ir+1)
      end do
!     backward difference
      dfdr(ien)= (f(ien)-f(ien-1))/(r(ien)-r(ien-1))
    end do ! panels

    end subroutine calc_derivative_panels_complex

    subroutine calc_derivative_panels_real(f,dfdr,r,nr,numpan,numrcut)
!   calculates the radial derivative of a function f 
!   takes care of the panels, if SOC from host
    use global, only: i4b, r8b, c8b, nrpts0

!   radial function 
    real(kind=c8b), intent(in)       :: f(1:nr)
!   derivative of function
    real(kind=c8b), intent(inout)    :: dfdr(1:nr)
!   --> Number of panels > 1
    integer(kind=i4b), intent(in)    :: numpan, numrcut(numpan+1)
!   rmesh
    integer(kind=i4b), intent(in)    :: nr
    real(kind=r8b), intent(in)       :: r(1:nr) 
!   ------------------------------------------------------------------------
    integer(kind=i4b)                :: ip, ist, ien
    real(kind=r8b)                   :: dri,dri1


    dfdr = 0.d0
!   Loop over panels
    do ip = 1, numpan
!   Begin and end of panels
      if (ip == 1) then
        ist = numrcut(ip) + 1            ! Start from 1  
      else
        ist = numrcut(ip) - nrpts0(1) + 2                                
      endif
      ien = numrcut(ip+1) - nrpts0(1) + 1                               
!     forward difference
      dfdr(ist)= (f(ist+1)-f(ist))/(r(ist+1)-r(ist))
!     central differences with a second order approach (for non uniform r grid)
      do ir=ist+1, ien-1
        dri=r(ir+1)-r(ir)
        dri1=r(ir)-r(ir-1)
        dfdr(ir)=-dri/dri1/(dri+dri1)*f(ir-1)
        dfdr(ir)=dfdr(ir)+(dri-dri1)/dri/dri1*f(ir)
        dfdr(ir)=dfdr(ir)+dri1/dri/(dri+dri1)*f(ir+1)
      end do
!     backward difference
      dfdr(ien)= (f(ien)-f(ien-1))/(r(ien)-r(ien-1))
    end do ! panels

    end subroutine calc_derivative_panels_real
end module mod_derivative_panels
