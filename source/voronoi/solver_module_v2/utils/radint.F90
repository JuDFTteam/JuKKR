! ************** Modified Function to include NPAN > 1 *************
  pure function radint(n,f,dx,numpan,numrcut)
! Integrates array x with weights dx coming from radial mesh
! If x(n1:n2) and n = n2 - n1 + 1, can be used to integrate over an interval
! dx has the effect of the change of variable to an equidistant mesh
! Same radial shift for all atoms  
  use global

  implicit none

  integer(kind=i4b), intent(in) :: n
  complex(kind=c8b), intent(in) :: f(n)
  real(kind=r8b),    intent(in) :: dx(n)
  complex(kind=c8b) :: radint
! Number of panels > 1
  integer(kind=i4b), intent(in) :: numpan
  integer(kind=i4b), intent(in) :: numrcut(numpan+1)
! -----------------------------------------------------------------
  complex(kind=c8b) :: work(n)
  integer(kind=i4b) :: ip, ist, ien
  complex(kind=c8b) :: radint_panel(numpan)

  radint = 0.d0

! dr/di weights   
  work = f*dx

! Looping over panels
  do ip = 1, numpan 

!   start and end of panels 
    if (ip == 1) then
      ist = numrcut(ip) + 1            ! Start from 1  
    else
      ist = numrcut(ip) - nrpts0(1) + 2                                
    endif
      ien = numrcut(ip+1) - nrpts0(1) + 1                                

    if ((ien-ist+1) < 2) then
!     nothing
      radint_panel(ip) = 0.d0
    else if ((ien-ist+1) == 2) then
!     trapezoidal rule
      radint_panel(ip) = 0.5d0*(work(ist) + work(ien))
    else if (mod((ien-ist+1),2) == 1) then
!   normal extended simpson rule
      radint_panel(ip) = 2.d0*sum(work(ist:ien:2))/3.d0
      radint_panel(ip) = radint_panel(ip) + 4.d0*sum(work(ist+1:ien:2))/3.d0
      radint_panel(ip) = radint_panel(ip) - (work(ist) + work(ien))/3.d0
    else
!   trapezoidal rule for first two points, might be inaccurate
      radint_panel(ip) = 0.5d0*(work(ist) + work(ist+1))
!   extended simpson rule on n-1 points 
      radint_panel(ip) = radint_panel(ip) + 2.d0*sum(work(ist+1:ien:2))/3.d0
      radint_panel(ip) = radint_panel(ip) + 4.d0*sum(work(ist+1:ien:2))/3.d0
      radint_panel(ip) = radint_panel(ip) - (work(ist+1) + work(ien))/3.d0
    end if
    radint = radint + radint_panel(ip)

  enddo ! panels     
      
! All done!
  end function radint


