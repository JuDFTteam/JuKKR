module mod_chebyshev
!-------------------------------------------------------------------------------
!> Summary: Routines which are needed to use the Chebyshev mesh 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
double complex,allocatable :: intweight(:)

contains 
!-------------------------------------------------------------------------------
!> Summary:  calculates the C matrix according to equation 5.36 in Bauer, PhD
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine getCmatrix(Ncheb,Cmatrix)
! calculates the C matrix according to
! equation 5.36 in Bauer, PhD
! function arguments are Ncheb the Chebyshev nodes x_k
! more information seeGonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
use nrtype, only: pi
implicit none
integer, intent(in)           :: ncheb
double precision, intent(out) :: Cmatrix(0:Ncheb,0:Ncheb)
!local
integer                       :: icheb1,icheb2

do icheb1=0,ncheb
  do icheb2=0,ncheb
    Cmatrix(icheb2,icheb1)=dcos(icheb1*pi*((Ncheb-icheb2)+0.5D0)/(Ncheb+1))
  end do
end do
end subroutine getCmatrix
!-------------------------------------------------------------------------------
!> Summary:  calculates the C matrix according to equation 5.36 in Bauer, PhD
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine getCCmatrix(Ncheb,rmesh,nrmesh,Cmatrix)
! calculates the C matrix according to:
! equation 5.36 in Bauer, PhD
!
! Note: function arguments are the values stored in rmesh
!       This routine is mostly used to interpolate from Chebyshev expansion
!       to an arbitrary mesh
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
implicit none
integer  :: ncheb
double precision :: rmesh(nrmesh)
double precision :: Cmatrix(1:nrmesh,0:Ncheb)
integer  :: icheb,nrmesh,ir

do ir=1,nrmesh
  do icheb=0,ncheb
    Cmatrix(ir,icheb)=cos(icheb*acos(rmesh(ir)))
  end do
end do
end subroutine getCCmatrix

!-------------------------------------------------------------------------------
!> Summary: set up the Lambda matrix which differentiates the coefficients of a 
!> Chebyshev expansion (equation 5.58 of Bauer, PhD) 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine getLambda(Ncheb,Lambda)
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion (equation 5.58 of Bauer, PhD)
implicit none
integer          :: Ncheb
double precision :: Lambda(0:Ncheb,0:Ncheb)
!local
integer icheb,icheb2
do icheb2=1,Ncheb,2
  Lambda(0,icheb2)=icheb2
end do
do icheb=1,Ncheb
  do icheb2=icheb+1,Ncheb,2
    Lambda(icheb,icheb2)=icheb2*2
  end do
end do
end subroutine

!-------------------------------------------------------------------------------
!> Summary: calculates a matrix 5.59 of Bauer, PhD 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine getCLambdaCinv(Ncheb,CLambdaCinv)
implicit none
! calculates a matrix 5.59 of Bauer, PhD
! This matrix can be used to calculate the derivative of a function fn, which
! is given in terms of the function values at the Chebyshev roots
!
integer          :: Ncheb
double precision :: CLambdaCinv(0:Ncheb,0:Ncheb)
!local
double precision :: Lambda(0:Ncheb,0:Ncheb)
double precision :: Cmatrix(0:Ncheb,0:Ncheb)
double precision :: Cinvmatrix(0:Ncheb,0:Ncheb)
double precision :: temp1(0:Ncheb,0:Ncheb)
 Lambda=(0.0D0,0.0D0)
 Cmatrix=(0.0D0,0.0D0)
 Cinvmatrix=(0.0D0,0.0D0)
 Lambda=(0.0D0,0.0D0)
 temp1=(0.0D0,0.0D0)

call getLambda(Ncheb,Lambda)
call getCinvmatrix(Ncheb,Cinvmatrix)
call getCmatrix(Ncheb,Cmatrix)

 temp1=matmat_dmdm(Lambda,Cinvmatrix)
 CLambdaCinv=matmat_dmdm(Cmatrix,temp1)

end subroutine
!-------------------------------------------------------------------------------
!> Summary: calculates the C**-1 matrix according to equation 5.39 in Bauer, PhD
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine getCinvmatrix(Ncheb,Cinvmatrix)
! calculates the C**-1 matrix according to equation 5.39 in Bauer, PhD
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
use nrtype, only: pi
implicit none
integer, intent(in)           :: ncheb
double precision, intent(out) :: Cinvmatrix(0:Ncheb,0:Ncheb)
!local
integer                       :: icheb1,icheb2
double precision              :: fac

fac=1.0D0/(Ncheb+1)
do icheb1=0,ncheb
  do icheb2=0,ncheb
    Cinvmatrix(icheb1,icheb2)=fac*dcos(icheb1*pi*((Ncheb-icheb2)+0.5D0)/(Ncheb+1))
  end do
  fac=2.0D0/(Ncheb+1)
end do

end subroutine getCinvmatrix

!-------------------------------------------------------------------------------
!> Summary: integrates an array arr1 containing function values of the Ncheb
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine intcheb_complex(Ncheb,arr1,result1)
! integrates an array arr1 containing function values of the Ncheb 
! Chebyshev roots and stores the results in result1
use nrtype, only: pi
implicit none
integer, intent(in)         :: ncheb
double complex, intent(in)  :: arr1(0:Ncheb)
double complex, intent(out) :: result1
integer :: icheb1,icheb2

if (.not. allocated(intweight)) then
  allocate(intweight(0:ncheb))
  intweight=1.0D0
  do icheb1=0,ncheb
    do icheb2=2,ncheb,2
      intweight(icheb1)=intweight(icheb1)+(-2.0D0/(icheb2**2-1.0D0))*dcos(icheb2*pi*(icheb1+0.5D0)/(Ncheb+1))
    end do
    intweight(icheb1)=intweight(icheb1)*2.0D0/(Ncheb+1)
  end do
end if

if (ubound(intweight,1)/=Ncheb) stop 'error1234'
result1=(0.0D0,0.0D0)
do icheb1=0,ncheb
  result1=result1+intweight(icheb1)*arr1(icheb1)
end do

end subroutine

!-------------------------------------------------------------------------------
!> Summary: differentiates a funtion stored in an array fn containing the function values 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine diff_cell(fn,cellnew,dfndr)
! differentiates a funtion stored in an array fn containing the function values.
! The array fn is divided into panels. Each panel is expanded into Chebyshev
! polynomials. 
! Each panel is differentiated separately by making use of equation 5.59 of Bauer, PhD
! of the Ncheb 
! Chebyshev roots and stores the resulting function in an array dfndr
 
use type_cellnew

implicit none
double precision        :: fn(:),dfndr(:)
type(cell_typenew)    :: cellnew
!local
double precision          :: CLambdaCinv(0:cellnew%ncheb,0:cellnew%ncheb)
double precision          :: widthfac
integer                    :: irstart,irstop,ipan
call getCLambdaCinv(cellnew%Ncheb,CLambdaCinv)
do ipan=1,cellnew%npan_tot
  irstart=cellnew%ipan_intervall(ipan-1)+1
  irstop = cellnew%ipan_intervall(ipan)
  widthfac = 2.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
  dfndr(irstart:irstop) = matvec_dmdm(CLambdaCinv,fn(irstart:irstop))
  dfndr(irstart:irstop) = dfndr(irstart:irstop)*widthfac
end do
end subroutine diff_cell


! subroutine diff2_cell(fn,cellnew,dfndr)
! use type_cellnew
! implicit none
! double precision        :: fn(:),dfndr(:)
! type(cell_typenew)    :: cellnew
! double precision        :: fntemp(0:cellnew%ncheb)
! !local
! double precision          :: CLambda2Cinv(0:cellnew%ncheb,0:cellnew%ncheb)
! double precision          :: widthfac
! integer                    :: irstart,irstop,ipan
! call getCLambda2Cinv(cellnew%Ncheb,CLambda2Cinv)
! do ipan=1,cellnew%npan_tot
!   irstart=cellnew%ipan_intervall(ipan-1)+1
!   irstop = cellnew%ipan_intervall(ipan)
!   widthfac = 4.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))**2
!   fntemp=fn(irstart:irstop)*widthfac
!   dfndr(irstart:irstop) = matvec_dmdm(CLambda2Cinv,fntemp)
! end do
! end subroutine diff2_cell

!-------------------------------------------------------------------------------
!> Summary: differentiates a complex funtion stored in an array fn containing the function values 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine diff_cell_complex(fn,cellnew,dfndr)
use type_cellnew

implicit none
double complex        :: fn(:),dfndr(:)
type(cell_typenew)    :: cellnew
!local
double precision          :: CLambdaCinv(0:cellnew%ncheb,0:cellnew%ncheb)
double complex            :: CLambdaCinv2(0:cellnew%ncheb,0:cellnew%ncheb)
double precision          :: widthfac
integer                    :: irstart,irstop,ipan

 call getCLambdaCinv(cellnew%Ncheb,CLambdaCinv)
 CLambdaCinv2=CLambdaCinv

do ipan=1,cellnew%npan_tot
  irstart=cellnew%ipan_intervall(ipan-1)+1
  irstop = cellnew%ipan_intervall(ipan)
  widthfac = 2.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
  dfndr(irstart:irstop) = matvec_zmzm(CLambdaCinv2,fn(irstart:irstop))
  dfndr(irstart:irstop) = dfndr(irstart:irstop)*widthfac
end do

end subroutine diff_cell_complex





!-------------------------------------------------------------------------------
!> Summary: integrate a function and returns its integral value
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine intcheb_cell(fn,cellnew,rho2ns_integrated)
!***********************************************************************
! integrate a function and returns its integral value
! The function is assumed to be devided into panels
! Each panel is expanded in Chebyshev polynomials and
! the array fn contains the function values at the roots
! transformed to the interval  cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1)
! see eq 5.40 in Bauer, PhD
! For each panel, the integration routine rho2ns_integrated is called 
! which calculates the integral over a single panel
!***********************************************************************
use type_cellnew
use mod_checknan
implicit none

double complex        :: fn(:)
type(cell_typenew)    :: cellnew
integer               :: irstart,irstop,ipan
double precision      :: widthfac
double complex        :: rho2ns_integrated,sum1,int1
rho2ns_integrated=(0.0D0,0.0D0)
sum1=(0.0D0,0.0D0)
do ipan=1,cellnew%npan_tot
  irstart=cellnew%ipan_intervall(ipan-1)+1
  irstop = cellnew%ipan_intervall(ipan)
  widthfac = 0.5D0*(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
  call intcheb_complex(cellnew%Ncheb,fn(irstart:irstop),int1)
    rho2ns_integrated=rho2ns_integrated+int1*widthfac
  end do

end subroutine intcheb_cell




!-------------------------------------------------------------------------------
!> Summary: double precision matrix-vector multiplication 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
      function matvec_dmdm(mat1,vec1)
      implicit none
      real(8), intent(in) :: mat1(:,:),vec1(:)
      real(8)             :: matvec_dmdm(size(mat1,1))
      integer             :: n,m
      m = size(mat1,1)
      n = size(mat1,2)
      if(size(vec1,1).ne.n) stop 'matmat_dmdm: dimensions of first input array differ.'
!       if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
!       if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
      call DGEMV('N',M,N,1.0D0,mat1,M,vec1,1,0.0D0,matvec_dmdm,1)
      end function matvec_dmdm

!-------------------------------------------------------------------------------
!> Summary: complex matrix-vector multiplication 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
      function matvec_zmzm(mat1,vec1)
      implicit none
      double complex, intent(in) :: mat1(:,:),vec1(:)
      double complex             :: matvec_zmzm(size(mat1,1))
      integer             :: n,m
      m = size(mat1,1)
      n = size(mat1,2)
      if(size(vec1,1).ne.n) stop 'matmat_dmdm: dimensions of first input array differ.'
!       if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
!       if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
      call ZGEMV('N',M,N,(1.0D0,0.0D0),mat1,M,vec1,1,(0.0D0,0.0D0),matvec_zmzm,1)
      end function matvec_zmzm


!-------------------------------------------------------------------------------
!> Summary: double precision matrix-matrix multiplication 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
      function matmat_dmdm(mat1,mat2)
      implicit none
      real(8), intent(in) :: mat1(:,:),mat2(:,:)
      real(8)             :: matmat_dmdm(size(mat1,1),size(mat2,2))
      integer             :: n
      n = size(mat1,1)
      if(size(mat1,2).ne.n) stop 'matmat_dmdm: dimensions of first input array differ.'
      if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
      if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
      call dgemm('N','N',n,n,n,1d0,mat1,n,mat2,n,0d0,matmat_dmdm,n)
      end function matmat_dmdm





end module mod_chebyshev
