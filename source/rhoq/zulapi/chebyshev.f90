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

subroutine intcheb_cell(cden,den,rpan_intervall,ipan_intervall, &
                        npan_tot,ncheb,irmdnew)
!***********************************************************************
! integrate the complex density of states for LM=1 
! gives the total complex charge which is then
! transformed to the xyz component of the magnetic 
! moment
!***********************************************************************
implicit none

integer           :: ncheb,npan_tot,irmdnew
integer           :: ipan_intervall(0:npan_tot)
double precision  :: rpan_intervall(0:npan_tot)
double complex    :: cden(irmdnew),den
integer           :: ir,irstart,irstop,ipan
double precision  :: widthfac
double complex    :: int1

den=(0.0D0,0.0D0)

  do ipan=1,npan_tot
    irstart=ipan_intervall(ipan-1)+1
    irstop = ipan_intervall(ipan)
    widthfac = 0.5D0*(rpan_intervall(ipan)-rpan_intervall(ipan-1))
    call intcheb_complex(ncheb,cden(irstart:irstop),int1)
    den=den+int1*widthfac
    end do

end subroutine intcheb_cell 

subroutine intcheb_complex(Ncheb,arr1,result1)
implicit none
integer, intent(in)         :: ncheb
double complex, intent(in)  :: arr1(0:Ncheb)
double complex, intent(out) :: result1
double precision            :: pi
double precision,allocatable  :: intweight(:)
integer :: icheb1,icheb2

pi=4d0*datan(1d0)
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

end subroutine intcheb_complex

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
double precision    :: fn(0:ncheb)
double precision    :: dfndr(0:ncheb)
double precision :: CLambdaCinv(0:Ncheb,0:Ncheb)

!needs to be checked!!!!!!1
call getCLambdaCinv(Ncheb,CLambdaCinv)
dfndr=matvec_dmdm(CLambdaCinv,fn)
end subroutine



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


      function matvec_zmzm(mat1,vec1)
      implicit none
      double complex, intent(in) :: mat1(:,:),vec1(:)
      double complex             :: matvec_zmzm(size(mat1,1))
      integer             :: n,m
      m = size(mat1,1)
      n = size(mat1,2)
      if(size(vec1,1).ne.n) stop 'matvec_zmzm: dimensions of first input array differ.'
      call ZGEMV('N',M,N,(1.0D0,0.0D0),mat1,M,vec1,1,(0.0D0,0.0D0),matvec_zmzm,1)
      end function matvec_zmzm



      function matmat_dmdm(mat1,mat2,Ncheb)
      implicit none
      integer             :: Ncheb,n
      double precision, intent(in) :: mat1(0:Ncheb,0:Ncheb),mat2(0:Ncheb,0:Ncheb)
      double precision             :: matmat_dmdm(Ncheb+1,Ncheb+1)
!      n = size(mat1,1)
!      if(size(mat1,2).ne.n) stop 'matmat_dmdm: dimensions of first input array differ.'
!      if(size(mat2,1).ne.n) stop 'matmat_dmdm: second input array has wrong dimensions.'
!      if(size(mat2,2).ne.n) stop 'matmat_dmdm: dimensions of second input array differ.'
      n=Ncheb+1
      call dgemm('N','N',n,n,n,1d0,mat1,n,mat2,n,0d0,matmat_dmdm,n)
      end function matmat_dmdm





