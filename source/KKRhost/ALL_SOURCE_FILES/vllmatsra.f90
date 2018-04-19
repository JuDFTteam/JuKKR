subroutine vllmatsra(vll0,vll,rmesh,lmsize,nrmax,eryd,lmax,lval_in,cmode)

   use Constants
   !************************************************************************************
   ! The perubation matrix for the SRA-equations are set up
   !************************************************************************************
   implicit none

   integer, intent(in)           :: lmax !< Maximum l component in wave function expansion
   integer, intent(in)           :: nrmax !< NTOTD*(NCHEBD+1)
   integer, intent(in)           :: lmsize
   integer, intent(in)           :: lval_in
   double complex, intent(in)    :: eryd
   character(len=*), intent(in)  :: cmode

   double precision, dimension(nrmax), intent(in)              :: rmesh
   double complex, dimension(lmsize,lmsize,nrmax), intent(in)  :: VLL0
   ! .. Output variables
   double complex, dimension(2*lmsize,2*lmsize,nrmax), intent(out) :: VLL
   ! .. Local variables
   integer                     :: ilm,lval,mval,ival,ir
   integer, dimension(lmsize)  :: loflm
   double complex              :: Mass,Mass0

   !************************************************************************************
   ! determine the bounds of the matricies to get the lm-expansion and the max. number
   ! of radial points
   !************************************************************************************

   !************************************************************************************
   ! calculate the index array to determine the L value of an LM index
   ! in case of spin-orbit coupling 2*(LMAX+1)**2 are used instead of (LMAX+1)**2
   ! the second half refers to the second spin and has the the same L value
   !************************************************************************************
   ilm=0

   if (lmsize==1) then
      loflm(1)=lval_in
   elseif ((lmax+1)**2 == lmsize) then
      do lval=0,lmax
         do mval = -lval,lval
            ilm=ilm+1
            loflm(ilm)=lval
         end do
      end do
   elseif (2* (lmax+1)**2 ==lmsize ) then
      do ival=1,2
         do lval=0,lmax
            do mval = -lval,lval
               ilm=ilm+1
               loflm(ilm)=lval
            end do
         end do
      end do
   else
      stop '[vllmatsra] error'
   end if

   vll=(0.0D0,0d0)

   if (cmode=='Ref=0') then
      vll(1:lmsize,1:lmsize,:)= vll0 !/cvlight

      do ir=1,nrmax
         do ival=1,lmsize
            lval=loflm(ival)
            Mass =cone+(eryd-vll0(ival,ival,ir))/cvlight**2
            Mass0=cone+eryd/cvlight**2

            !************************************************************************************
            ! Conventional potential matrix
            !************************************************************************************

            vll(lmsize+ival,lmsize+ival,ir)= -vll0(ival,ival,ir)/cvlight**2 ! TEST 9/22/2011
            vll(ival,ival,ir)=vll(ival,ival,ir)+ (1.0D0/Mass-1.0D0/Mass0)*lval*(lval+1)/rmesh(ir)**2

            !************************************************************************************
            ! The pertubation matrix is changed in the following way
            !
            !     from  / V11  V12 \   to    / V21  V22 \
            !           \ V21  V22 /         \-V11 -V12 /
            ! because of the convention used for the left solution
            !************************************************************************************
         end do !ival

      end do !ir
   elseif     (cmode=='Ref=Vsph') then
      vll(lmsize+1:2*lmsize,1:lmsize,:)=vll0
   endif

end subroutine vllmatsra
