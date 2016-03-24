subroutine calc_orbitalmoment(lmax,lmsize,Loperator)
implicit none
integer                         :: lmax,lmsize
double complex                  :: Loperator(lmsize,lmsize,3)
integer,save                    :: first=1
integer                         :: lval,ilm
double complex,allocatable      :: lorbit_onel(:,:,:)
integer                         :: lmmax,lstart,lstop
lmmax=(lmax+1)**2
Loperator=(0.0D0,0.0D0)

Loperator(1,1,1)=(0.0D0,0.0D0)
Loperator(1,1,2)=(0.0D0,0.0D0)
Loperator(1,1,3)=(0.0D0,0.0D0)

do lval=1,lmax
  allocate(lorbit_onel(  2*lval+1,2*lval+1,3 )  )  
  lorbit_onel=(0.0D0,0.0D0)
  call calc_orbit_onel(lval,Lorbit_onel)
  lstart=((lval-1)+1)**2+1
  lstop=((lval  )+1)**2
  Loperator(lstart:lstop,lstart:lstop,1) = Lorbit_onel(:,:,1)
  Loperator(lstart:lstop,lstart:lstop,2) = Lorbit_onel(:,:,2)
  Loperator(lstart:lstop,lstart:lstop,3) = Lorbit_onel(:,:,3)
  deallocate(lorbit_onel)
end do
if (lmsize/=lmmax) then
  Loperator(lmmax+1:lmsize,lmmax+1:lmsize,:) = Loperator(:lmmax,:lmmax,:)
end if

! if (first==1) then
!  open(unit=423492157,file='out_Lx')
!  open(unit=423492158,file='out_Ly')
!  open(unit=423492159,file='out_Lz')
!  do ilm=1,lmsize
!    write(423492157,'(5000F)'),Loperator(ilm,:,1)
!    write(423492158,'(5000F)'),Loperator(ilm,:,2)
!    write(423492159,'(5000F)'),Loperator(ilm,:,3)
!  end do
!  close(423492157);close(423492158);close(423492159)
! end if

first=0
end subroutine calc_orbitalmoment





subroutine calc_orbit_onel(lval,Lorbit_onel)
implicit none
integer        :: lval
double complex Lorbit_onel(2*lval+1,2*lval+1,3)
double complex L_z (-lval:lval,-lval:lval),L_x (-lval:lval,-lval:lval),L_y(-lval:lval,-lval:lval)
double complex L_dn(-lval:lval,-lval:lval),L_up(-lval:lval,-lval:lval)
integer :: i1
double precision :: lfac
double complex,parameter :: icompl=(0.0D0,1.0D0)


L_z=(0.0D0,0.0D0)
L_x=(0.0D0,0.0D0)
L_y=(0.0D0,0.0D0)
L_dn=(0.0D0,0.0D0)
L_up=(0.0D0,0.0D0)

Lorbit_onel=(0.0D0,0.0D0)
       do i1=-lval,lval
         L_z(-i1,i1)=-icompl*i1 
       end do 





       IF (lval>0) then

         lfac=sqrt(lval*(lval+1d0))/sqrt(2d0)
         l_dn(0,-1)=-icompl*lfac
         l_dn(0,1)=lfac
         l_dn(-1,0)=icompl*lfac
         l_dn(1,0)=-lfac
       
         IF (lval > 1) then

            do i1=2,lval

              lfac=0.5d0*SQRT(lval*(lval+1d0)-i1*(i1-1d0))
              l_dn(-i1,-i1+1)=-lfac
              l_dn(-i1,i1-1)=icompl*lfac
              l_dn(i1,-i1+1)=-icompl*lfac
              l_dn(i1,i1-1)=-lfac

              lfac=0.5d0*SQRT(lval*(lval+1d0)-(i1-1d0)*i1)
              l_dn(-i1+1,-i1)=lfac
              l_dn(-i1+1,i1)=icompl*lfac
              l_dn(i1-1,-i1)=-icompl*lfac
              l_dn(i1-1,i1)=lfac

            end do

         END IF
       END IF




       IF (lval>0) then

         lfac=sqrt(lval*(lval+1d0))/sqrt(2d0)
         l_up(0,-1)=-icompl*lfac
         l_up(0,1)=-lfac
         l_up(-1,0)=icompl*lfac
         l_up(1,0)=lfac
       
         IF (lval > 1) then

            do i1=2,lval

              lfac=0.5d0*SQRT(lval*(lval+1d0)-i1*(i1-1d0))
              l_up(-i1,-i1+1)=lfac
              l_up(-i1,i1-1)=icompl*lfac
              l_up(i1,-i1+1)=-icompl*lfac
              l_up(i1,i1-1)=lfac

              lfac=0.5d0*SQRT(lval*(lval+1d0)-(i1-1d0)*i1)
              l_up(-i1+1,-i1)=-lfac
              l_up(-i1+1,i1)=icompl*lfac
              l_up(i1-1,-i1)=-icompl*lfac
              l_up(i1-1,i1)=-lfac

            end do

         END IF
       END IF


L_x =  0.5D0          * (L_up+L_dn)
L_y = -0.5D0 * icompl * (L_up-L_dn)


Lorbit_onel(:,:,1)=L_x
Lorbit_onel(:,:,2)=L_y
Lorbit_onel(:,:,3)=L_z

end subroutine calc_orbit_onel

