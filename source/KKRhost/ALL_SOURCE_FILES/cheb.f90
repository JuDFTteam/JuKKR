!module mod_chebyshev
!double complex,allocatable :: intweight(:)


!contains 


subroutine getCmatrix(Ncheb,Cmatrix)
! calculates the C matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
implicit none
integer, intent(in)           :: ncheb
double precision, intent(out) :: Cmatrix(0:Ncheb,0:Ncheb)
double precision              :: pi
!local
integer                       :: icheb1,icheb2

pi=4d0*datan(1d0)
do icheb1=0,ncheb
  do icheb2=0,ncheb
    ! maybe incorrect
    Cmatrix(icheb2,icheb1)=dcos(icheb1*pi*((Ncheb-icheb2)+0.5D0)/(Ncheb+1))
  end do
end do
end subroutine getCmatrix


subroutine getCinvmatrix(Ncheb,Cinvmatrix)
! calculates the C**-1 matrix according to:
! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
implicit none
integer, intent(in)           :: ncheb
double precision, intent(out) :: Cinvmatrix(0:Ncheb,0:Ncheb)
!local
double precision              :: pi
integer                       :: icheb1,icheb2
double precision              :: fac

pi=4d0*datan(1d0)
fac=1.0D0/(Ncheb+1)
do icheb1=0,ncheb
  do icheb2=0,ncheb
    Cinvmatrix(icheb1,icheb2)=fac*dcos(icheb1*pi*((Ncheb-icheb2)+0.5D0)/(Ncheb+1))
  end do
  fac=2.0D0/(Ncheb+1)
end do

end subroutine getCinvmatrix


!subroutine diff_cell(fn,cellnew,dfndr)
!use type_cellnew
!
!implicit none
!double precision        :: fn(:),dfndr(:)
!type(cell_typenew)    :: cellnew
!!local
!double precision          :: CLambdaCinv(0:cellnew%ncheb,0:cellnew%ncheb)
!double precision          :: widthfac
!integer                    :: irstart,irstop,ipan
!call getCLambdaCinv(cellnew%Ncheb,CLambdaCinv)
!do ipan=1,cellnew%npan_tot
!  irstart=cellnew%ipan_intervall(ipan-1)+1
!  irstop = cellnew%ipan_intervall(ipan)
!  widthfac = 2.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
!  dfndr(irstart:irstop) = matvec_dmdm(CLambdaCinv,fn(irstart:irstop))
!  dfndr(irstart:irstop) = dfndr(irstart:irstop)*widthfac
!end do
!end subroutine diff_cell


! subroutine diff_cell_alternative(fn,cellnew,dfndr)
! use type_cellnew
! 
! implicit none
! double precision        :: fn(:),dfndr(:)
! type(cell_typenew)    :: cellnew
! !local
! double precision          :: Cinvmatrix(0:cellnew%ncheb,0:cellnew%ncheb)
! double precision          :: Cmatrix(0:cellnew%ncheb,0:cellnew%ncheb)
! double precision          :: ccoeff(0:cellnew%ncheb)
! 
! double precision          :: widthfac
! integer                    :: irstart,irstop,ipan
! 
! 
! call getCinvmatrix(Ncheb,Cinvmatrix)
! call getCmatrix(Ncheb,Cmatrix)
! 
! 
! ! do ipan=0,cellnew%ncheb
! !   write(*,*)
! ! end do
! do ipan=1,cellnew%npan_tot
!   irstart=cellnew%ipan_intervall(ipan-1)+1
!   irstop = cellnew%ipan_intervall(ipan)
! !   write(*,*) irstart, irstop
!   widthfac = 2.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
! 
! 
!   ccoeff = matvec_zmzm(Cinvmatrix,fn(irstart:irstop))
!   ccoeffd(cellnew%ncheb)=0.0D0
!   ccoeffd(cellnew%ncheb-1)=2*(n-1)*ccoeff(cellnew%ncheb)
!   do ival=cellnew%ncheb,1,-1
!     ccoeffd(ival-1)=ccoeffd(cellnew%ncheb+1)+2.0D0*ival*ccoeff(ival)
!   end do  
! 
!   dfndr(irstart:irstop) = matvec_zmzm(Cmatrix,ccoeff)
! 
! 
! 
! !   dfndr(irstart:irstop) = matvec_dmdm(CLambdaCinv,fn(irstart:irstop))
!   dfndr(irstart:irstop) = dfndr(irstart:irstop)*widthfac
! end do
! ! stop
! end subroutine diff_cell_alternative







!subroutine diff2_cell(fn,cellnew,dfndr)
!use type_cellnew
!implicit none
!double precision        :: fn(:),dfndr(:)
!type(cell_typenew)    :: cellnew
!double precision        :: fntemp(0:cellnew%ncheb)
!!local
!double precision          :: CLambda2Cinv(0:cellnew%ncheb,0:cellnew%ncheb)
!double precision          :: widthfac
!integer                    :: irstart,irstop,ipan
!call getCLambda2Cinv(cellnew%Ncheb,CLambda2Cinv)
!do ipan=1,cellnew%npan_tot
!  irstart=cellnew%ipan_intervall(ipan-1)+1
!  irstop = cellnew%ipan_intervall(ipan)
!  widthfac = 4.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))**2
!  fntemp=fn(irstart:irstop)*widthfac
!  dfndr(irstart:irstop) = matvec_dmdm(CLambda2Cinv,fntemp)
!!   dfndr(irstart:irstop) = dfndr(irstart:irstop)*widthfac
!end do
!end subroutine diff2_cell










!subroutine diff_cell_complex(fn,cellnew,dfndr)
!use type_cellnew
!
!implicit none
!double complex        :: fn(:),dfndr(:)
!type(cell_typenew)    :: cellnew
!!local
!double precision          :: CLambdaCinv(0:cellnew%ncheb,0:cellnew%ncheb)
!double complex            :: CLambdaCinv2(0:cellnew%ncheb,0:cellnew%ncheb)
!double precision          :: widthfac
!integer                    :: irstart,irstop,ipan
!
! call getCLambdaCinv(cellnew%Ncheb,CLambdaCinv)
! CLambdaCinv2=CLambdaCinv
!
!do ipan=1,cellnew%npan_tot
!  irstart=cellnew%ipan_intervall(ipan-1)+1
!  irstop = cellnew%ipan_intervall(ipan)
!  widthfac = 2.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
!  dfndr(irstart:irstop) = matvec_zmzm(CLambdaCinv2,fn(irstart:irstop))
!  dfndr(irstart:irstop) = dfndr(irstart:irstop)*widthfac
!end do

!end subroutine diff_cell_complex

subroutine getCCmatrix(Ncheb,rmesh,nrmesh,Cmatrix)
! calculates the C matrix according to:
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


subroutine getLambda(Ncheb,Lambda)
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion 
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

subroutine getCLambdaCinv(Ncheb,CLambdaCinv)
implicit none
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion
integer          :: Ncheb
double precision :: CLambdaCinv(0:Ncheb,0:Ncheb)
!local
double precision :: Lambda(0:Ncheb,0:Ncheb)
double precision :: Cmatrix(0:Ncheb,0:Ncheb)
double precision :: Cinvmatrix(0:Ncheb,0:Ncheb)
double precision :: temp1(0:Ncheb,0:Ncheb)
integer n
 Lambda=(0.0D0,0.0D0)
 Cmatrix=(0.0D0,0.0D0)
 Cinvmatrix=(0.0D0,0.0D0)
 Lambda=(0.0D0,0.0D0)
 temp1=(0.0D0,0.0D0)

call getLambda(Ncheb,Lambda)
call getCinvmatrix(Ncheb,Cinvmatrix)
call getCmatrix(Ncheb,Cmatrix)
n=Ncheb+1
 call dgemm('N','N',n,n,n,1d0,Lambda,n,Cinvmatrix,n,0d0,temp1,n)
 call dgemm('N','N',n,n,n,1d0,Cmatrix,n,temp1,n,0d0,CLambdaCinv,n)
! temp1=matmat_dmdm(Lambda,Cinvmatrix,Ncheb)
! CLambdaCinv=matmat_dmdm(Cmatrix,temp1,Ncheb)

end subroutine


subroutine getCLambda2Cinv(Ncheb,CLambda2Cinv)
implicit none
! set up the Lambda matrix which differentiates the coefficients of an
! Chebyshev expansion
integer          :: Ncheb
double precision :: CLambda2Cinv(0:Ncheb,0:Ncheb)
!local
double precision :: Lambda(0:Ncheb,0:Ncheb)
double precision :: Cmatrix(0:Ncheb,0:Ncheb)
double precision :: Cinvmatrix(0:Ncheb,0:Ncheb)
double precision :: temp1(0:Ncheb,0:Ncheb)
double precision :: temp2(0:Ncheb,0:Ncheb)
double precision :: matmat_dmdm
 Lambda=(0.0D0,0.0D0)
 Cmatrix=(0.0D0,0.0D0)
 Cinvmatrix=(0.0D0,0.0D0)
 Lambda=(0.0D0,0.0D0)
 temp1=(0.0D0,0.0D0)

call getLambda(Ncheb,Lambda)
call getCinvmatrix(Ncheb,Cinvmatrix)
call getCmatrix(Ncheb,Cmatrix)
! write(*,'(50000E)') Lambda
! write(*,'(50000E)') Cinvmatrix
! write(*,'(50000E)') Cmatrix

 temp1=matmat_dmdm(Lambda,Lambda,Ncheb)

 temp2=matmat_dmdm(temp1,Cinvmatrix,Ncheb)
!  temp2=matmat_dmdm(Lambda,temp1)
 CLambda2Cinv=matmat_dmdm(Cmatrix,temp2,Ncheb)

end subroutine


subroutine diffCheb(fn,ncheb,dfndr)
implicit none
integer :: ncheb
double precision    :: fn(0:ncheb)
double precision    :: dfndr(0:ncheb)
double precision :: CLambdaCinv(0:Ncheb,0:Ncheb)
double precision, external :: matvec_dmdm
! real(8), external :: matvec_dmdm

!needs to be checked!!!!!!1
call getCLambdaCinv(Ncheb,CLambdaCinv(0:ncheb,0:ncheb))
dfndr(0:ncheb)=matvec_dmdm(CLambdaCinv(0:ncheb,0:ncheb),fn(0:ncheb))
end subroutine



      double precision function matvec_dmdm(mat1,vec1)
      implicit none
      double precision, intent(in) :: mat1(:,:),vec1(:)
!       double precision             :: matvec_dmdm(size(mat1,1))
!       real(8), intent(in) :: mat1(:,:),vec1(:)
!       real(8)             :: matvec_dmdm(size(mat1,1))
      integer             :: n,m
      m = size(mat1,1)
      n = size(mat1,2)
      if(size(vec1,1).ne.n) stop 'matmat_dmdm: dimensions of first input array differ.'
!       if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
!       if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
      call DGEMV('N',M,N,1.0D0,mat1,M,vec1,1,0.0D0,matvec_dmdm,1)
      end function matvec_dmdm


      double complex function matvec_zmzm(mat1,vec1)
      implicit none
      double complex, intent(in) :: mat1(:,:),vec1(:)
!       double complex             :: matvec_zmzm(size(mat1,1))
      integer             :: n,m
      m = size(mat1,1)
      n = size(mat1,2)
      if(size(vec1,1).ne.n) stop 'matvec_zmzm: dimensions of first input array differ.'
      call ZGEMV('N',M,N,(1.0D0,0.0D0),mat1,M,vec1,1,(0.0D0,0.0D0),matvec_zmzm,1)
      end function matvec_zmzm



      double precision function matmat_dmdm(mat1,mat2,Ncheb)
      implicit none
      integer             :: Ncheb,n
      double precision, intent(in) :: mat1(0:Ncheb,0:Ncheb),mat2(0:Ncheb,0:Ncheb)
!       double precision             :: matmat_dmdm(Ncheb+1,Ncheb+1)
!      n = size(mat1,1)
!      if(size(mat1,2).ne.n) stop 'matmat_dmdm: dimensions of first input array differ.'
!      if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
!      if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
      n=Ncheb+1
      call dgemm('N','N',n,n,n,1d0,mat1,n,mat2,n,0d0,matmat_dmdm,n)
      end function matmat_dmdm





!end module mod_chebyshev
