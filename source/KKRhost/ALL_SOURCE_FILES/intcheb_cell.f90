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

subroutine intcheb_complex(ncheb,arr1,result1)
implicit none
integer, intent(in)         :: ncheb
double complex, intent(in)  :: arr1(0:ncheb)
double complex, intent(out) :: result1
double precision            :: pi
double precision  :: intweight(0:ncheb)
integer :: icheb1,icheb2

pi=4d0*datan(1d0)
intweight=1.0D0
  do icheb1=0,ncheb
    do icheb2=2,ncheb,2
      intweight(icheb1)=intweight(icheb1)+(-2.0D0/(icheb2**2-1.0D0))*dcos(icheb2*pi*(icheb1+0.5D0)/(Ncheb+1))
    end do
    intweight(icheb1)=intweight(icheb1)*2.0D0/(Ncheb+1)
  end do

result1=(0.0D0,0.0D0)
do icheb1=0,ncheb
  result1=result1+intweight(icheb1)*arr1(icheb1)
end do

end subroutine
