subroutine calc_orbitalmoment(lmax,lmsize,Loperator)

   use Constants
   use Profiling

   implicit none
   integer, intent(in)                                      :: lmax,lmsize
   double complex, dimension(lmsize,lmsize,3), intent(out)  :: Loperator
   integer,save                                             :: first=1
   integer                                                  :: lval
   double complex, dimension(:,:,:), allocatable            :: lorbit_onel
   integer                                                  :: lmmax,lstart,lstop,i_stat,i_all
   lmmax=(lmax+1)**2
   Loperator=CZERO

   Loperator(1,1,1)=CZERO
   Loperator(1,1,2)=CZERO
   Loperator(1,1,3)=CZERO

   do lval=1,lmax
      allocate(lorbit_onel(  2*lval+1,2*lval+1,3 ),stat=i_stat)
      call memocc(i_stat,product(shape(lorbit_onel))*kind(lorbit_onel),'lorbit_onel','calc_orbitalmoment')
      lorbit_onel=CZERO
      call calc_orbit_onel(lval,Lorbit_onel)
      lstart=((lval-1)+1)**2+1
      lstop=((lval  )+1)**2
      Loperator(lstart:lstop,lstart:lstop,1) = Lorbit_onel(:,:,1)
      Loperator(lstart:lstop,lstart:lstop,2) = Lorbit_onel(:,:,2)
      Loperator(lstart:lstop,lstart:lstop,3) = Lorbit_onel(:,:,3)
      !
      i_all=-product(shape(lorbit_onel))*kind(lorbit_onel)
      deallocate(lorbit_onel,stat=i_stat)
      call memocc(i_stat,i_all,'lorbit_onel','calc_orbitalmoment')
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

   use Constants

   implicit none
   integer, intent(in) :: lval
   double complex, dimension(2*lval+1,2*lval+1,3), intent(out) :: Lorbit_onel
   double complex, dimension(-lval:lval,-lval:lval)   :: L_z
   double complex, dimension(-lval:lval,-lval:lval)   :: L_x
   double complex, dimension(-lval:lval,-lval:lval)   :: L_y
   double complex, dimension(-lval:lval,-lval:lval)   :: L_dn
   double complex, dimension(-lval:lval,-lval:lval)   :: L_up
   integer :: i1
   double precision :: lfac
   !
   L_z   =CZERO
   L_x   =CZERO
   L_y   =CZERO
   L_dn  =CZERO
   L_up  =CZERO
   !
   Lorbit_onel=CZERO
   do i1=-lval,lval
      L_z(-i1,i1)=-CI*i1
   end do
   !
   if (lval>0) then
      !
      lfac=sqrt(lval*(lval+1d0))/sqrt(2d0)
      l_dn(0,-1)=-CI*lfac
      l_dn(0,1)=lfac
      l_dn(-1,0)=CI*lfac
      l_dn(1,0)=-lfac
      !
      if (lval > 1) then
         do i1=2,lval
            lfac=0.5d0*SQRT(lval*(lval+1d0)-i1*(i1-1d0))
            l_dn(-i1,-i1+1)=-lfac
            l_dn(-i1,i1-1)=CI*lfac
            l_dn(i1,-i1+1)=-CI*lfac
            l_dn(i1,i1-1)=-lfac
            !
            lfac=0.5d0*SQRT(lval*(lval+1d0)-(i1-1d0)*i1)
            l_dn(-i1+1,-i1)=lfac
            l_dn(-i1+1,i1)=CI*lfac
            l_dn(i1-1,-i1)=-CI*lfac
            l_dn(i1-1,i1)=lfac
         end do
      end if
   end if
   !
   if (lval>0) then
      !
      lfac=sqrt(lval*(lval+1d0))/sqrt(2d0)
      l_up(0,-1)=-CI*lfac
      l_up(0,1)=-lfac
      l_up(-1,0)=CI*lfac
      l_up(1,0)=lfac
      !
      if (lval > 1) then
         do i1=2,lval
            lfac=0.5d0*SQRT(lval*(lval+1d0)-i1*(i1-1d0))
            l_up(-i1,-i1+1)=lfac
            l_up(-i1,i1-1)=CI*lfac
            l_up(i1,-i1+1)=-CI*lfac
            l_up(i1,i1-1)=lfac
            !
            lfac=0.5d0*SQRT(lval*(lval+1d0)-(i1-1d0)*i1)
            l_up(-i1+1,-i1)=-lfac
            l_up(-i1+1,i1)=CI*lfac
            l_up(i1-1,-i1)=-CI*lfac
            l_up(i1-1,i1)=-lfac
         end do
      end if
   end if
   !
   L_x =  0.5D0      * (L_up+L_dn)
   L_y = -0.5D0 * CI * (L_up-L_dn)
   !
   Lorbit_onel(:,:,1)=L_x
   Lorbit_onel(:,:,2)=L_y
   Lorbit_onel(:,:,3)=L_z

end subroutine calc_orbit_onel
