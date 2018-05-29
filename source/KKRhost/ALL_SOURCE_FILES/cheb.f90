!module mod_cheb

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


subroutine getCCmatrix(Ncheb,rmesh,nrmesh,Cmatrix)
   ! calculates the C matrix according to:
   ! Gonzalez et al, Journal of Computational Physics 134, 134-149 (1997)
   implicit none
   integer, intent(in) :: ncheb,nrmesh
   double precision, intent(in)  :: rmesh(nrmesh)
   double precision, intent(out) :: Cmatrix(1:nrmesh,0:Ncheb)
   integer  :: icheb,ir
   
   do ir=1,nrmesh
     do icheb=0,ncheb
       Cmatrix(ir,icheb)=cos(dfloat(icheb)*dacos(rmesh(ir)))
     end do
   end do
end subroutine getCCmatrix


subroutine getLambda(Ncheb,Lambda)
   ! set up the Lambda matrix which differentiates the coefficients of an
   ! Chebyshev expansion 
   implicit none
   integer, intent(in)           :: Ncheb
   double precision, intent(out) :: Lambda(0:Ncheb,0:Ncheb)
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
   !function
   double precision :: matmat_dmdm

   Lambda=(0.0D0,0.0D0)
   Cmatrix=(0.0D0,0.0D0)
   Cinvmatrix=(0.0D0,0.0D0)
   Lambda=(0.0D0,0.0D0)
   temp1=(0.0D0,0.0D0)
  
   call getLambda(Ncheb,Lambda)
   call getCinvmatrix(Ncheb,Cinvmatrix)
   call getCmatrix(Ncheb,Cmatrix)
  
   temp1=matmat_dmdm(Lambda,Lambda,Ncheb)
   temp2=matmat_dmdm(temp1,Cinvmatrix,Ncheb)
   CLambda2Cinv=matmat_dmdm(Cmatrix,temp2,Ncheb)
end subroutine


subroutine diffCheb(fn,ncheb,dfndr)
   implicit none
   integer :: ncheb
   double precision    :: fn(0:ncheb)
   double precision    :: dfndr(0:ncheb)
   double precision :: CLambdaCinv(0:Ncheb,0:Ncheb)
   !function
   double precision :: matvec_dmdm
  
   !needs to be checked!!!!!!1
   call getCLambdaCinv(Ncheb,CLambdaCinv(0:ncheb,0:ncheb))
   dfndr(0:ncheb)=matvec_dmdm(ncheb,CLambdaCinv(0:ncheb,0:ncheb),fn(0:ncheb))
end subroutine


! helper functions

double precision function matvec_dmdm(ncheb,mat1,vec1)
   implicit none
   integer, intent(in) :: ncheb
   double precision, intent(in) :: mat1(0:ncheb,0:ncheb),vec1(0:ncheb)
   integer             :: n,m
   m = size(mat1,1)
   n = size(mat1,2)
   if(size(vec1,1).ne.n) stop 'matmat_dmdm: dimensions of first input array differ.'
   call DGEMV('N',M,N,1.0D0,mat1,M,vec1,1,0.0D0,matvec_dmdm,1)
end function matvec_dmdm

double precision function matmat_dmdm(mat1,mat2,Ncheb)
   implicit none
   integer             :: Ncheb,n
   double precision, intent(in) :: mat1(0:Ncheb,0:Ncheb),mat2(0:Ncheb,0:Ncheb)
   n=Ncheb+1
   call dgemm('N','N',n,n,n,1d0,mat1,n,mat2,n,0d0,matmat_dmdm,n)
end function matmat_dmdm


!end module mod_cheb
