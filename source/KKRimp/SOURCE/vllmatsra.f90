module mod_vllmatsra

contains

subroutine vllmatsra(vll,rmesh,eryd,lmax,lval_in,cmode)
!************************************************************************************
! The potential for the SRA-equations are set up
! according to formula 4.107 of Bauer, PhD thesis
!
!V = ( 1/2M = 1/2M_0 l(l+1)/r**2 + V_LL;     0     )
!    (             0                       2M-2M_0 )
!
!************************************************************************************
use mod_physic_params, only: cvlight
use mod_mathtools, only: inverse
implicit none
!interface
  DOUBLE COMPLEX,allocatable :: VLL(:,:,:)
  double precision,allocatable :: rmesh(:)
  double complex              :: eryd
  integer                      :: lmax
  character(len=*)             :: cmode
!local
  double complex,allocatable :: vlltemp(:,:,:)
  integer                     :: nrmax
  integer                      :: lmsize
  integer                     :: ilm,lval,mval,ival,ir
  integer,allocatable       :: loflm(:)
  double complex              :: Mass,Mass0
  double complex,parameter    :: cone=(1.0D0,0.0D0)
  double complex,parameter    :: czero=(0.0D0,0.0D0)
!  double precision,parameter :: cvlight = 274.0720442d0
  integer                     :: lval_in
  double complex,allocatable   :: rmass_matrix(:,:),rmass0_matrix(:,:),rmass0_inv(:,:)
  character(len=*),parameter  :: rmass_mode='spherical'


!************************************************************************************
! determine the bounds of the matricies to get the lm-expansion and the max. number
! of radial points
!************************************************************************************
nrmax=ubound(vll,3)
lmsize=ubound(vll,1)

if (rmass_mode=='non-spherical') then
  allocate(rmass_matrix(lmsize,lmsize), rmass0_matrix(lmsize,lmsize),rmass0_inv(lmsize,lmsize) )
end if

!************************************************************************************
! calculate the index array to determine the L value of an LM index
! in case of spin-orbit coupling 2*(LMAX+1)**2 are used instead of (LMAX+1)**2
! the second half refers to the second spin and has the the same L value
!************************************************************************************
allocate(loflm(lmsize))
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

allocate(vlltemp(lmsize,lmsize,nrmax))
vlltemp=vll
deallocate(vll)
allocate(vll(2*lmsize,2*lmsize,nrmax))
vll=(0.0D0,0.0D0)

if     (cmode=='Ref=0') then
!************************************************************************************
! Reference system is free space full V_LL is used
!************************************************************************************

  vll(1:lmsize,1:lmsize,:)= vlltemp !/cvlight

  if (rmass_mode=='non-spherical') then
    rmass0_matrix = czero
    rmass0_inv    = czero
    do ival=1, lmsize  
      rmass0_matrix(ival,ival)=cone+eryd/cvlight**2
      rmass0_inv(ival,ival)   =cone / rmass0_matrix(ival,ival)
    end do

    vll(lmsize+1:2*lmsize,lmsize+1:2*lmsize,:)= -vlltemp/cvlight**2

  end if

  do ir = 1, nrmax

    if (rmass_mode=='non-spherical') then ! option not used in the code
                                          ! might be deleted in the future
                                          ! uses the off-diag. elements of B_LL
                                          ! in eq. 4.67 section 4.3 of Bauer, PhD

      rmass_matrix = rmass0_matrix - vlltemp(:,:,ir)/cvlight**2

      call inverse(lmsize,rmass_matrix) 

      do ival=1, lmsize  
        rmass_matrix(ival,ival) = rmass_matrix(ival,ival) - rmass0_inv(ival,ival) 
      end do

      do ival=1, lmsize  
        lval=loflm(ival)
        rmass_matrix(:,ival) =  rmass_matrix(:,ival)*lval*(lval+1)/rmesh(ir)**2
      end do

      vll(1:lmsize,1:lmsize,ir)= vll(1:lmsize,1:lmsize,ir) + rmass_matrix

    elseif (rmass_mode=='spherical') then  ! standart mode which is always used

      do ival=1, lmsize  
        lval=loflm(ival)
        Mass =cone+(eryd-vlltemp(ival,ival,ir))/cvlight**2
        Mass0=cone+eryd/cvlight**2

        vll(lmsize+ival,lmsize+ival,ir)= -vlltemp(ival,ival,ir)/cvlight**2 ! TEST 9/22/2011

        vll(ival,ival,ir)=vll(ival,ival,ir)+ (1.0D0/Mass-1.0D0/Mass0)*lval*(lval+1)/rmesh(ir)**2

      end do !ival

   else
      stop '[vllmatsra] rmass mode not known'
    end if

  end do

elseif (cmode=='Ref=Vsph') then
!************************************************************************************
! Reference system are radial solutions
! obmit the diagonal part of the potential V_LL=0  
!   => M = M0 => M-M_0 =0
!************************************************************************************

  vll(lmsize+1:2*lmsize,1:lmsize,:)= vlltemp

else
  stop '[vllmatsra] mode not known'
end if




end subroutine vllmatsra

end module mod_vllmatsra
