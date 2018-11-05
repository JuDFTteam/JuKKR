module mod_cheb2oldgridc
!-------------------------------------------------------------------------------
!> Summary: Interpolates from the Chebychev mesh to the 'old' radial mesh 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
use mod_chebyshev
integer :: first=1
contains
!-------------------------------------------------------------------------------
!> Summary: Interpolates from the Chebychev mesh to the 'old' radial mesh 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
subroutine cheb2oldgridc(cell,cellnew,ncheb,lmmaxpot,arrayin,arrayout)
use type_cell
use type_cellnew
implicit none
type(cell_type) :: cell
type(cell_typenew) :: cellnew
integer  :: ncheb,lmmaxpot
double complex :: arrayin(cellnew%nrmaxnew,lmmaxpot)
double complex :: arrayout(cell%nrmax,lmmaxpot)
!local
integer :: in,ir,ir2,ir3,it,ilm
integer :: intsub(2,cell%nrmax)
double precision :: Cinvmatrix(0:cellnew%ncheb,0:cellnew%ncheb)
double precision,allocatable :: CCmatrix(:,:)
double complex :: alphaparams(0:ncheb,lmmaxpot)
double precision  :: rmeshnorm(cell%nrmax)
double precision  :: halfsum,halfdiffinv,rshift,tol
parameter(tol=1.d-10)
!parameter(tol=1.d-13)

! divide the mesh into subintervals
intsub(1,:)=0
intsub(2,:)=-1
in=1




ir=1
it=1
do while (ir<=cell%nrmax .and. in<=cellnew%npan_tot)
  if ( abs(cell%rmesh(ir)-cellnew%rpan_intervall(in))<10e-14) then
    intsub(it,in)=ir
    intsub(2,in)=ir
    it=1
    ir=ir+1
    in=in+1
  elseif ( (cell%rmesh(ir)<= cellnew%rpan_intervall(in) ) .and. &
           (cell%rmesh(ir)>= cellnew%rpan_intervall(in-1)) ) then
  intsub(it,in)=ir
  intsub(2,in)=ir
  it=2
  ir=ir+1
  elseif (cell%rmesh(ir)< cellnew%rpan_intervall(in-1))  then
    intsub(1,in)=0
    intsub(2,in)=-1
    it=1
    ir=ir+1
  else
    it=1
    in=in+1
  end if
end do
 
if (abs(cell%rmesh(cell%nrmax)- cellnew%rpan_intervall(cellnew%npan_tot))>1d-15) then
   write(*,*) cell%rmesh(cell%nrmax),cellnew%rpan_intervall(cellnew%npan_tot),cell%rmesh(cell%nrmax)-cellnew%rpan_intervall(cellnew%npan_tot)
  stop 'error with the new and old mesh'
else
  
end if


do in=1,cellnew%npan_tot
  do ir=intsub(1,in),intsub(2,in)
    ! added checking against tolerance level since sometimes a 10**-15 differece gave a wrong positive in this error checking
    if (     cell%rmesh(ir)-cellnew%rpan_intervall(in-1)<-tol &
        .or. cell%rmesh(ir)-cellnew%rpan_intervall(in)>tol   ) then
      write(*,*) '------------------------------------------------------------',in
      write(*,*) 'panel ',in
      write(*,*) 'borders r=',cellnew%rpan_intervall(in-1),cellnew%rpan_intervall(in)
      write(*,*) 'meshpoints in rmesh n=',intsub(1,in),intsub(2,in)
      write(*,*) 'rval :'
      !write(*,'(50000E)') cell%rmesh(intsub(1,in):intsub(2,in))
      write(*,'(50000E)') cell%rmesh(intsub(1,in))
      write(*,'(50000E)') cell%rmesh(intsub(2,in))
      write(*,*) '------------------------------------------------------------',in
      write(*,*) '[cheb2oldgridc] something went wrong in this panel'
      stop
    end if
  end do
end do


if (first==1) then
  do in=1,cellnew%npan_tot
    write(1337,*) '------------------------------------------------------------',in
    write(1337,*) 'panel ',in
    write(1337,*) 'borders r=',cellnew%rpan_intervall(in-1),cellnew%rpan_intervall(in)
    write(1337,*) 'meshpoints in rmesh n=',intsub(1,in),intsub(2,in)
    write(1337,*) 'rval :'
    write(1337,'(50000E)') cell%rmesh(intsub(1,in):intsub(2,in)) 
    write(1337,*) '------------------------------------------------------------',in
    do ir=intsub(1,in),intsub(2,in)
      write(1337,*) cell%rmesh(ir)
    end do
  end do
end if

call getCinvmatrix(Ncheb,Cinvmatrix)

in=0
do while (in<=cellnew%npan_tot)
  in=in+1
  if (intsub(2,in)<intsub(1,in)) cycle
  alphaparams=0.0D0
  do ilm=1,lmmaxpot
  do ir2=0,ncheb
    do ir=0,ncheb
      alphaparams(ir2,ilm)=alphaparams(ir2,ilm)+ Cinvmatrix(ir2,ir)* arrayin(cellnew%ipan_intervall(in-1)+1+ir,ilm)
    end do
  end do
  end do


  ! Transform to normalized coordinates between [-1,1]
  ! Shift upper & lower end to be centered around zero
  halfsum = 0.5d0 * (cellnew%rpan_intervall(in) + cellnew%rpan_intervall(in-1)) ! 0.5 * (upper + lower)
  halfdiffinv = 1.d0/(cellnew%rpan_intervall(in) - halfsum)             ! 1. / (0.5*(upper - lower))

  do ir=intsub(1,in),intsub(2,in)
    ir2=ir+1-intsub(1,in)
    !    The following was numerically inaccurate, changed to halfsum & halfdiffinv. Phivos 22.07.2014
    !    rmeshnorm(ir2)=(  2.d0*rmesh(ir)-( rpan_intervall(in) + rpan_intervall(in-1) )  ) &
    !                                / ( rpan_intervall(in) - rpan_intervall(in-1) )
    !    The following was accurate but recalculating sum and difference too many times.
    !    rmeshnorm(ir2)=(  rmesh(ir)- rpan_intervall(in) + rmesh(ir) - rpan_intervall(in-1)   ) &
    !                                / ( rpan_intervall(in) - rpan_intervall(in-1) )
    rshift = cell%rmesh(ir) - halfsum
    rmeshnorm(ir2) = rshift * halfdiffinv
    ! Value should not be outside interval [-1,+1]:
    if (abs(rmeshnorm(ir2))>1.0D0) then
       if (abs(rmeshnorm(ir2))-1.0D0<tol) then
          rmeshnorm(ir2)=sign(1.0D0,rmeshnorm(ir2)) ! set to +/- 1
       else
          write(*,*) 'cheb2oldgrid: argument of chebyshev plynomial exceeds interval [-1,+1]:',rmeshnorm(ir2)
          stop 'cheb2oldgrid'
       end if
    end if
  end do

  allocate(CCmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

  call getCCmatrix(Ncheb,rmeshnorm(1:intsub(2,in)-intsub(1,in)+1),intsub(2,in)-intsub(1,in)+1,CCmatrix)
  do ilm=1,lmmaxpot
  do ir2=intsub(1,in),intsub(2,in)
    do ir=0,ncheb
      ir3=ir2-intsub(1,in)+1
      arrayout(ir2,ilm)=arrayout(ir2,ilm)+ CCmatrix(ir3,ir)* alphaparams(ir,ilm)
    end do
  end do
  end do !ilm=1,lmmaxpot
  deallocate(CCmatrix) !(CCmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

end do !in
!   stop
first=0
end subroutine cheb2oldgridc

!-------------------------------------------------------------------------------
!> Summary: Complex matrix multiplication 
!> Author:
!> Category: KKRimp, radial-grid
!>           
!-------------------------------------------------------------------------------
      function matmat_zmzm(mat1,mat2)
      implicit none
      complex(16), intent(in) :: mat1(:,:),mat2(:,:)
      complex(16)             :: matmat_zmzm(size(mat1,1),size(mat2,2))
      integer                :: n1,n,n2
      n1 = size(mat1,1)
      n  = size(mat1,2)
      n2 = size(mat2,2)
      if(size(mat2,1).ne.n) stop 'matmat_zmzm: dimensions of matrices are inconsistent.'
      call zgemm('N','N',n1,n2,n,(1d0,0d0),mat1,n1,mat2,n,(0d0,0d0),matmat_zmzm,n1)
      end function matmat_zmzm

end module mod_cheb2oldgridc
