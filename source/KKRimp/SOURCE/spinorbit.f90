module mod_spinorbit
  
  contains

subroutine spinorbit(lmax,zatom,eryd,cellnew,nrmaxnew,nspin,vpotll,theta,phi,ncoll,mode)
use mod_rotatespinframe, only: rotatematrix
use type_cellnew
use mod_mathtools
use mod_chebyshev, only: getCLambdaCinv
use mod_physic_params, only: cvlight
use mod_config, only: config_testflag
use mod_basistransform
implicit none
!interface
integer            :: lmax
double precision   :: zatom
double complex     :: eryd
type(cell_typenew) :: cellnew
integer            :: nrmaxnew
integer            :: nspin
double complex     :: vpotll(2*(lmax+1)**2,2*(lmax+1)**2,nrmaxnew)
double precision   :: theta
double precision   :: phi
integer            :: ncoll
character(len=*)   :: mode
!local
double complex            :: LSHAM(2*(lmax+1)**2,2*(lmax+1)**2)
integer                   :: lmmax
integer                   :: irstart,irstop,ipan,ir
double precision          :: widthfac
double precision          :: vpot(nrmaxnew),vpot2(nrmaxnew)
double precision          :: dvpotdr(nrmaxnew)
double complex            :: RMASS,temp !(nrmaxnew)
double complex            :: HSOFAC(nrmaxnew)
double precision          :: CLambdaCinv(0:cellnew%ncheb,0:cellnew%ncheb)
integer                    :: ilm1,ilm2
double precision           :: RNUCL,ATN

lmmax=(lmax+1)**2
vpot=0.0D0

if (nspin==1) then
  stop 'function spinorbit used but NSPIN=1'
else
  do ir=1,nrmaxnew
    vpot(ir)=0.5D0*(cellnew%vpotnew(ir,1,1)+cellnew%vpotnew(ir,1,2))
  end do 
end if


!************************************************************************************
! derivative of VPOT (without the core potential)
!************************************************************************************

! integration is done by a matrix-matrix multiplication using
! equation 5.53-5.54 of Bauer, PhD thesis
call getCLambdaCinv(cellnew%Ncheb,CLambdaCinv)

do ipan=1,cellnew%npan_tot
  irstart=cellnew%ipan_intervall(ipan-1)+1
  irstop = cellnew%ipan_intervall(ipan)
  widthfac = 2.0D0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
  dvpotdr(irstart:irstop) = matvec_dmdm(CLambdaCinv,vpot(irstart:irstop))
  dvpotdr(irstart:irstop) = dvpotdr(irstart:irstop)*widthfac
end do


!************************************************************************************
! add the derivate of the core potential
!************************************************************************************

if ( config_testflag('finitecore') ) then

  if (zatom>24) then
    atn=-16.1532921+2.70335346*zatom
  else 
    atn=0.03467714+2.04820786*zatom
  end if

  rnucl=1.2d0/0.529177d0*atn**(1./3d0)*1.d-5  

  write(1337,*) 'size of the core : ',rnucl
  do ir=1,nrmaxnew
    if (cellnew%rmeshnew(ir) .le. rnucl) then
      dvpotdr(ir)=dvpotdr(ir)+2d0*zatom*cellnew%rmeshnew(ir)/rnucl**3d0
    else
      dvpotdr(ir)=dvpotdr(ir)+2d0*zatom/cellnew%rmeshnew(ir)**2d0
    end if
  end do

else

  do ir=1,nrmaxnew
      dvpotdr(ir)=dvpotdr(ir)+2d0*zatom/cellnew%rmeshnew(ir)**2d0
  end do

end if

!************************************************************************************
! calucate the L*S operator in the real spherical harmonics basis
!************************************************************************************
CALL SPIN_ORBIT_COMPL(LMAX,LMMAX,LSHAM)
! write(*,*) 'test so',eryd
! LSHAM=(1.0D0,1.0D0)

if ( config_testflag('nospinflip') ) then
  lsham((lmax+1)**2+1:,:(lmax+1)**2)=(0.0D0,0.0D0)
  lsham(:(lmax+1)**2,(lmax+1)**2+1:)=(0.0D0,0.0D0)
end if

if (ncoll==1) then                                        ! added 2.2.2012
  call rotatematrix(LSHAM,theta, phi,LMMAX,'glob->loc')   ! added 2.2.2012
end if                                                    ! added 2.2.2012

if (mode=='conjg') then                                    ! not needed, might be deleted
  do ilm1=1,2*(lmax+1)**2
    do ilm2=1,2*(lmax+1)**2
      LSHAM(ilm2,ilm1)=conjg(LSHAM(ilm2,ilm1))
    end do
  end do
elseif (mode=='transpose') then                            ! used for the left solution
  do ilm1=1,2*(lmax+1)**2
    do ilm2=1,ilm1-1
      temp            =LSHAM(ilm2,ilm1)
      LSHAM(ilm2,ilm1)=LSHAM(ilm1,ilm2)
      LSHAM(ilm1,ilm2)=temp
    end do
  end do
elseif (mode=='transpose2x2') then                         ! not needed, might be deleted
  do ilm1=1,(lmax+1)**2
    do ilm2=(lmax+1)**2+1,2*(lmax+1)**2
      temp            =LSHAM(ilm2,ilm1)
      LSHAM(ilm2,ilm1)=LSHAM(ilm1,ilm2)
      LSHAM(ilm1,ilm2)=temp
    end do
  end do
elseif (mode=='transpose2x2conjg') then                    ! not needed, might be deleted
  do ilm1=1,(lmax+1)**2
    do ilm2=(lmax+1)**2+1,2*(lmax+1)**2
      temp            =LSHAM(ilm2,ilm1)
      LSHAM(ilm2,ilm1)=LSHAM(ilm1,ilm2)
      LSHAM(ilm1,ilm2)=temp
    end do
  end do
  do ilm1=1,2*(lmax+1)**2
    do ilm2=1,2*(lmax+1)**2
      LSHAM(ilm2,ilm1)=conjg(LSHAM(ilm2,ilm1))
    end do
  end do
elseif (mode=='1') then
else
  stop'[spinorbit] mode not known'
end if

!************************************************************************************
! calculate the prefacor of the spin-orbit coupling hamiltonian
! and adding the resulting term to the potential
! check Heers, PhD thesis  or Bauer, PhD thesis for details on the prefactor
!************************************************************************************
DO ir=1,nrmaxnew

  RMASS=0.5d0 - 0.5d0/cvlight**2*((VPOT(IR)-ERYD) - 2.0D0*ZATOM/cellnew%RMESHNEW(IR))

  HSOFAC(ir)=1d0/(2d0*RMASS**2*cvlight**2*cellnew%rmeshnew(ir))*DVPOTDR(ir)

  VPOTLL(:,:,ir)=VPOTLL(:,:,ir)+HSOFAC(ir)*LSHAM

END DO

end subroutine spinorbit

!************************************************************************************
! the subroutine REL_MASS calculates the relativistic mass,
! depending on the radial position.   
!************************************************************************************
SUBROUTINE REL_MASS(RMASS,RM,VPOT,ZATOM,ERYD,cvlight,IRMD)
  implicit none 
!interface
  integer          ::   IRMD
  DOUBLE PRECISION ::   RM(IRMD),VPOT(IRMD)
  DOUBLE COMPLEX   ::   ERYD
  DOUBLE PRECISION ::   cvlight,ZATOM
!output 
  DOUBLE PRECISION ::     RMASS(IRMD), &! radial derivative of the spherical input potential
                         DIFFV
!local 
  integer     ::     nr

do nr=1,IRMD
  DIFFV = (VPOT(NR)-ERYD) - 2.0D0*ZATOM/RM(NR)
  RMASS(nr)=0.5d0 - 0.5d0/cvlight**2*DIFFV
end do

END SUBROUTINE REL_MASS


! C ************************************************************************
      SUBROUTINE SPIN_ORBIT_COMPL(LMAX,LMMAXD,L_S)

      implicit none
! C ************************************************************************
! c      in this subroutine the matrix L*S is calculated for the basis of
! c      real spherical harmonics 

! c

       integer, intent(in)     ::     lmax,LMMAXD

       double complex, intent(out)    ::     L_S(LMMAXD*2,LMMAXD*2)
! c       complex                  ::     L_S(LMMAXD*2,LMMAXD*2)

! c  local variables 
       integer                 ::     i1,i2,i1l,RL,lm1,lm2
       double complex                 ::     icompl
       double complex,allocatable     ::     LS_L(:,:)


       icompl=(0d0,1d0) 

! c       write(6,*) "lmax:",lmax
! c       write(6,*) "LMMAXD",LMMAXD

! c       allocate(LS_L((2*LMAX+1)*2,(2*LMAX+1)*2))

!        CALL CINIT((2*LMMAXD)**2,L_S)
         L_S=(0.0D0,0.0D0)
! c       write(6,*) "after CINIT"

       do RL=0,lmax
! c         write(6,*) "rl:",rl

         ALLOCATE(LS_L((2*RL+1)*2,(2*RL+1)*2))
!          CALL CINIT(((2*RL+1)*2)**2,LS_L)
         LS_L=(0.0D0,0.0D0)
         CALL SPIN_ORBIT_ONE_L(RL,LS_L)
! c        write(6,*) "SPIN_ORBIT_ONE_L"
! c        do lm1=1,(2*RL+1)*2
! c          do lm2=1,2*(2*RL+1)
! c            write(6,*) lm1,lm2,LS_L(lm1,lm2)
! c          end do
! c        end do

         do lm1=1,(2*RL+1)*2

           IF (lm1 <= 2*RL+1 ) THEN
             do lm2=1,(2*RL+1)
               L_S(RL**2+Lm1,RL**2+lm2)=0.5d0*LS_L(lm1,lm2)
             end do 
             do lm2=(2*RL+1)+1,(2*RL+1)*2
               L_S(RL**2+Lm1,LMMAXD+RL**2-(2*RL+1)+lm2)= &
                                           0.5d0*LS_L(lm1,lm2)
             end do 
           ELSE
             do lm2=1,(2*RL+1)
               L_S(LMMAXD+RL**2-(2*RL+1)+Lm1,RL**2+lm2)= &
                                           0.5d0*LS_L(lm1,lm2)
             end do 
             do lm2=(2*RL+1)+1,(2*RL+1)*2
               L_S(LMMAXD+RL**2-(2*RL+1)+Lm1,LMMAXD+RL**2-(2*RL+1)+lm2)= &
                                           0.5d0*LS_L(lm1,lm2)
             end do
           END IF

         end do    !lm1

         DEALLOCATE(LS_L)
! c        write(6,*) "SPIN_ORBIT_GESAMT"
! c        do lm1=1,LMMAXD*2
! c          do lm2=1,LMMAXD*2
! c            write(16+rl,*) lm1,lm2,L_S(lm1,lm2)
! c          end do
! c        end do


       end do     !rl=0,lmax

! c        write(16,*) "SPIN_ORBIT_GESAMT"
! c        do lm1=1,LMMAXD*2
! c          do lm2=1,LMMAXD*2
! c            write(16,*) lm1,lm2,L_S(lm1,lm2)
! c          end do
! c        end do

! c       deallocate(LS_L)
! c       write(6,*) "end of spin_orbit alle l"

       END SUBROUTINE





! C ************************************************************************
      SUBROUTINE SPIN_ORBIT_ONE_L(lmax,L_S)

      implicit none
! C ************************************************************************
! c      in this subroutine the matrix L*S is calculated for the basis of
! c      real spherical harmonics 
! 
! c      schematically it has the form
! c      (  -L_z    L_+  )
! c      (  L_-     L_z  )
! c

       integer, intent(in)         ::    lmax
       double complex, intent(out)        ::    L_S((2*lmax+1)*2,(2*lmax+1)*2)

! c  local variables 
       integer                     ::    i1,i2,i1l 
       double complex                     ::    icompl
       double complex,allocatable         ::    l_min(:,:)
       double complex,allocatable         ::    l_up(:,:)
       double precision            ::    lfac



       icompl=(0d0,1d0) 


       allocate(l_min(-lmax:lmax,-lmax:lmax))
       allocate(l_up(-lmax:lmax,-lmax:lmax))

! c  initialize the matrix 

       do i1=1,(2*lmax+1)*2
         do i2=1,(2*lmax+1)*2
           L_S(i2,i1)=0d0
         end do
       end do

       do i1=-lmax,lmax
         do i2=-lmax,lmax
           L_min(i2,i1)=0d0
           L_up(i2,i1)=0d0
         end do
       end do

! c  fill the second and the forth quadrant with L_z
! c (-L_z,respectively) 


       do i1=1,2*lmax+1
         i1l=i1-lmax-1       ! the value of m (varies from -l to +l)  
         i2=2*lmax+1-(i1-1)  

         L_S(i2,i1)=-icompl*i1l 

       end do 

       do i1=2*lmax+2,(2*lmax+1)*2
         i1l=i1-lmax-1-(2*lmax+1)       ! the value of m (varies from -l to +l)  
         i2=(2*lmax+1)*2-(i1-(2*lmax+2))  

         L_S(i2,i1)=icompl*i1l 

       end do 


! c      write(6,*) "after lz"

! c  implement now L_- in the third quadrant

       IF (lmax>0) then

         lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
         l_min(0,-1)=-icompl*lfac
! c         l_min(0,-1)=icompl*lfac
         l_min(0,1)=lfac
         l_min(-1,0)=icompl*lfac
         l_min(1,0)=-lfac
       
         IF (lmax > 1) then

            do i1=2,lmax

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-i1*(i1-1d0))
              l_min(-i1,-i1+1)=-lfac
              l_min(-i1,i1-1)=icompl*lfac
              l_min(i1,-i1+1)=-icompl*lfac
              l_min(i1,i1-1)=-lfac

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-(i1-1d0)*i1)
              l_min(-i1+1,-i1)=lfac
              l_min(-i1+1,i1)=icompl*lfac
              l_min(i1-1,-i1)=-icompl*lfac
              l_min(i1-1,i1)=lfac

            end do

         END IF
       END IF

! c      write(6,*) "after l down"

! c       open(unit=13,file="l_down", form="formatted")
 
        do i1=-lmax,lmax
         do i2=-lmax,lmax
!           L_S(i2+lmax+1,i1+3*lmax+2)=l_min(i2,i1)
!     transpose l_min          
!           L_S(i2+lmax+1,i1+3*lmax+2)=l_min(i1,i2)
           L_S(i2+3*lmax+2,i1+lmax+1)=l_min(i1,i2)
! c           write(13,"((2I5),(4e17.9))") i2,i1, l_min(i2,i1)
         end do
       end do

! c       close(13)

! c  implement now L_+ in the   quadrant

       IF (lmax>0) then

         lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
         l_up(0,-1)=-icompl*lfac
         l_up(0,1)=-lfac
         l_up(-1,0)=icompl*lfac
         l_up(1,0)=lfac
       
         IF (lmax > 1) then

            do i1=2,lmax

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-i1*(i1-1d0))
              l_up(-i1,-i1+1)=lfac
              l_up(-i1,i1-1)=icompl*lfac
              l_up(i1,-i1+1)=-icompl*lfac
              l_up(i1,i1-1)=lfac

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-(i1-1d0)*i1)
              l_up(-i1+1,-i1)=-lfac
              l_up(-i1+1,i1)=icompl*lfac
              l_up(i1-1,-i1)=-icompl*lfac
              l_up(i1-1,i1)=-lfac

            end do

         END IF
       END IF

! c      write(6,*) "after l up"

! c      open(unit=13,file="l_up", form="formatted")
! c
        do i1=-lmax,lmax
         do i2=-lmax,lmax
!           L_S(i2+3*lmax+2,i1+lmax+1)=l_up(i2,i1)
!     transpose l_up          
           L_S(i2+lmax+1,i1+3*lmax+2)=l_up(i1,i2)
!           L_S(i2+3*lmax+2,i1+lmax+1)=l_up(i1,i2)
! c           write(13,"((2I5),(4e17.9))") i2,i1, l_up(i2,i1)
         end do
       end do

! c      close(13)


       deallocate(l_min)
       deallocate(l_up)


       END SUBROUTINE














SUBROUTINE SPIN_ORBIT_ONE_L_OLD(lmax,L_S)
! here the 1x1 block is still spin up
! and the 2x2 block is spin down ( swantje's convention )
! not used anymore
  implicit none
! C ************************************************************************
! c      in this subroutine the matrix L*S is calculated for the basis of
! c      real spherical harmonics 

! c      schematically it has the form
! c      (  L_z    L_-  )
! c      (  L_+   -L_z  )
! c

       integer, intent(in)         ::    lmax
       double complex, intent(out)        ::    L_S((2*lmax+1)*2,(2*lmax+1)*2)

! c  local variables 
       integer                     ::    i1,i2,i1l 
       double complex                     ::    icompl
       double complex,allocatable         ::    l_min(:,:)
       double complex,allocatable         ::    l_up(:,:)
       double precision            ::    lfac



       icompl=(0d0,1d0) 


       allocate(l_min(-lmax:lmax,-lmax:lmax))
       allocate(l_up(-lmax:lmax,-lmax:lmax))

! c  initialize the matrix 

       do i1=1,(2*lmax+1)*2
         do i2=1,(2*lmax+1)*2
           L_S(i2,i1)=0d0
         end do
       end do

       do i1=-lmax,lmax
         do i2=-lmax,lmax
           L_min(i2,i1)=0d0
           L_up(i2,i1)=0d0
         end do
       end do

!  fill the second and the forth quadrant with L_z
! (-L_z,respectively) 


       do i1=1,2*lmax+1
         i1l=i1-lmax-1       ! the value of m (varies from -l to +l)  
         i2=2*lmax+1-(i1-1)  

         L_S(i2,i1)=icompl*i1l 

       end do 

       do i1=2*lmax+2,(2*lmax+1)*2
         i1l=i1-lmax-1-(2*lmax+1)       ! the value of m (varies from -l to +l)  
         i2=(2*lmax+1)*2-(i1-(2*lmax+2))  

         L_S(i2,i1)=-icompl*i1l 

       end do 

! c       open(unit=13,file="LS_lz", form="formatted")
! c
! c       do i1=1,(2*lmax+1)*2
! c        do i2=1,(2*lmax+1)*2
! c          write(13,"((2I5),(4e17.9))") i2,i1, L_S(i2,i1)
! c        end do
! c      end do

!      close(13)

!      write(6,*) "after lz"

!  implement now L_- in the third quadrant

       IF (lmax>0) then

         lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
         l_min(0,-1)=-icompl*lfac
!         l_min(0,-1)=icompl*lfac
         l_min(0,1)=lfac
         l_min(-1,0)=icompl*lfac
         l_min(1,0)=-lfac
       
         IF (lmax > 1) then

            do i1=2,lmax

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-i1*(i1-1d0))
              l_min(-i1,-i1+1)=-lfac
              l_min(-i1,i1-1)=icompl*lfac
              l_min(i1,-i1+1)=-icompl*lfac
              l_min(i1,i1-1)=-lfac

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-(i1-1)*(i1))
              l_min(-i1+1,-i1)=lfac
              l_min(-i1+1,i1)=icompl*lfac
              l_min(i1-1,-i1)=-icompl*lfac
              l_min(i1-1,i1)=lfac

            end do

         END IF
       END IF

! c      write(6,*) "after l down"

! c       open(unit=13,file="l_down", form="formatted")
 
        do i1=-lmax,lmax
         do i2=-lmax,lmax
!            L_S(i2+lmax+1,i1+3*lmax+2)=l_min(i2,i1)
!     transpose l_min          
           L_S(i2+lmax+1,i1+3*lmax+2)=l_min(i1,i2)
! c           write(13,"((2I5),(4e17.9))") i2,i1, l_min(i2,i1)
         end do
       end do

! c       close(13)

! c  implement now L_+ in the   quadrant

       IF (lmax>0) then

         lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
         l_up(0,-1)=-icompl*lfac
         l_up(0,1)=-lfac
         l_up(-1,0)=icompl*lfac
         l_up(1,0)=lfac
       
         IF (lmax > 1) then

            do i1=2,lmax

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-i1*(i1-1d0))
              l_up(-i1,-i1+1)=lfac
              l_up(-i1,i1-1)=icompl*lfac
              l_up(i1,-i1+1)=-icompl*lfac
              l_up(i1,i1-1)=lfac

              lfac=0.5d0*SQRT(lmax*(lmax+1d0)-(i1-1)*(i1))
              l_up(-i1+1,-i1)=-lfac
              l_up(-i1+1,i1)=icompl*lfac
              l_up(i1-1,-i1)=-icompl*lfac
              l_up(i1-1,i1)=-lfac

            end do

         END IF
       END IF

!      write(6,*) "after l up"

!      open(unit=13,file="l_up", form="formatted")
!
        do i1=-lmax,lmax
         do i2=-lmax,lmax
!            L_S(i2+3*lmax+2,i1+lmax+1)=l_up(i2,i1)
!     transpose l_up          
           L_S(i2+3*lmax+2,i1+lmax+1)=l_up(i1,i2)
!           write(13,"((2I5),(4e17.9))") i2,i1, l_up(i2,i1)
         end do
       end do

!      close(13)


       deallocate(l_min)
       deallocate(l_up)


       END SUBROUTINE

subroutine transposematrix(matrix1)
implicit none
double complex matrix1(:,:),temp1
integer :: dim1,ival1,ival2
dim1=ubound(matrix1,1)
if (ubound(matrix1,2)/=dim1) stop'error in transposematrix'
do ival1=1,dim1
  do ival2=ival1+1,dim1
  temp1=matrix1(ival1,ival2)
  matrix1(ival1,ival2)=matrix1(ival2,ival1)
  matrix1(ival2,ival1)=temp1
  end do
end do
end subroutine transposematrix

end module mod_spinorbit

