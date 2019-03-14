module mod_interpolatecell
contains

!-------------------------------------------------------------------------------
!> Summary: ! This routine 
!> 1) distributes panels over the radial mesh of the single site problem. 
!> 2) Interpolates the potential or shape function of the old mesh points
!>    to the new Chebyshev panels
!> 1) The mesh is devided into 3 different parts:
!> | |  logarithmic   |     equidistant  |      shapefn       |     regions
!> | |  ------------- | ---------------  | -----------------  |     regions
!> | |||| |  |    |       |    |    |    |    |    |     |    |     panels
!> 0 rmin          rlogpan        rmesh(cell%NRMIN_NS)      R_max   limits
!>
!>1st region: logarithmic distribution of the panels (like in mesh points in ASA)
!>2nd region: equistitant distribution of the panels
!>3rd region: distribution of the panels according to the kinks for the shape function
!>            (like it is done for the old mesh)
!>
!> 2) Interpolation is done using a spline method
!> 
!> Author: Who wrote this subroutine
!> Category: new-mesh, old-mesh, radial-grid, shape-functions, potential, kkrimp
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------

subroutine interpolatecell(nstart,vpot,lmvpot,cell,lmaxatom,cellnew,ispin,nspin,config,cmode,testpot)
use mod_gauntharmonics, only: gauntcoeff
use type_cell
use type_cellnew
use type_config
use mod_interpolpot
use mod_vllmat
use mod_gauntharmonics, only: gauntcoeff
use mod_rllsll
use mod_config, only: config_testflag
use mod_interpolspline

implicit none
integer                                   :: nstart
type(cell_type)                           ::  cell
type(cell_typenew)                        ::  cellnew
integer                                   :: lmaxatom
integer                                   :: lmvpot
real(kind=dp),intent(in)                  ::  vpot(nstart:cell%nrmaxd,lmvpot)
type(config_type)                         :: config

real(kind=dp),allocatable                 ::  vpotnew(:,:)
real(kind=dp),allocatable                 ::  testpot(:,:)
integer                                   :: ispin,nspin,ishift
integer                                   :: ir,lm1,lmpot,lmmax,ipan,ipan2
integer                                   :: irmin,irmax,irminnew,irmaxnew
integer                                   :: npan_logtemp
double precision                          :: rmin,rmax,rval,fac !,rpan_logmesh
integer                                   :: first=1
integer                                   :: ilinmesh_panshift,ilogmesh_panshift,itemp
character(len=*)                          :: cmode

! This routine 
! 1) distributes panels over the radial mesh of the single site problem. 
! 2) Interpolates the potential or shape function of the old mesh points
!    to the new Chebyshev panels

! 1) The mesh is devided into 3 different parts:

! | |  logarithmic   |     equidistant  |      shapefn       |     regions
! | |  ------------- | ---------------  | -----------------  |     regions
! | |||| |  |    |       |    |    |    |    |    |     |    |     panels
! 0 rmin          rlogpan        rmesh(cell%NRMIN_NS)      R_max   limits

!1st region: logarithmic distribution of the panels (like in mesh points in ASA)
!2nd region: equistitant distribution of the panels
!3rd region: distribution of the panels according to the kinks for the shape function
!            (like it is done for the old mesh)

! 2) Interolation is done using a spline method


!  ############################################################################
!                    CREATE RADIAL PANELS
!  ############################################################################

lmmax=(lmaxatom+1)**2
lmpot=(2*lmaxatom+1)**2

 cellnew%npan_log  = config%npan_log ! panels in the log. part
 cellnew%npan_eq   = config%npan_eq  ! panels in the eq. part
 cellnew%npan_inst = cell%npan-1     ! panels in the shapefn part
 cellnew%npan_tot  = cellnew%npan_log+cellnew%npan_eq+cellnew%npan_inst

 cellnew%ncheb=config%ncheb !64 !16 !16 !16
 cellnew%nrmaxnew=cellnew%npan_tot*(cellnew%ncheb+1)

if (.not. allocated(cellnew%rpan_intervall)) then
  allocate(cellnew%rpan_intervall(0:cellnew%npan_tot))
  allocate(cellnew%ipan_intervall(0:cellnew%npan_tot))
  cellnew%rpan_intervall=0.0D0
  cellnew%ipan_intervall=0
end if


! set the region for the region where panels are distributed 
! logarithmically
rmin= config%rmin     ! 
rmax= config%rlogpan  ! set the

! set most inner panel to 0 if a negative value is set in the config file
if ( config%rmin< -1.0D-10 ) then
  rmin=cell%rmesh(2) ! because cell%rmesh(1) is always 0!
end if

rval=0
fac=config%npan_logfac !4.0D0/3.0D0
ishift=0


! the old mesh has a discontinuity where the non-spherical mesh is set to
! 0 for the inner region. To avoid an inaccurate interpolation, we set here
! an additional separation in panels

! if discontinuity is within the logarithmic part then set the variables 
! accordingly
if (config%rlogpan>cell%rmesh(cell%NRMIN_NS)) then 
  ilogmesh_panshift=1 ! additional panel inside the log. region
  ilinmesh_panshift=0 ! additional panel not inside the eq. region
else
  ilogmesh_panshift=0 ! additional panel not inside the log. region
  ilinmesh_panshift=1! additional panel inside the eq. region
end if

! check if the discontinuity is within the logarithmic mesh
! if its not, then calculation stops, since the inclusion of the discontinuity
! in the equidistant part is not implemented yet
if (ilinmesh_panshift==1) then 
    stop '[interpolatecell] non-spherical part of the potential needs to be in the inner log panel'
!  ############################################################################
!   NOT IMPLEMENTED JET: Code fraction leads to problem using -03
!  lines below need to be put into the DO LOOP for the EQ mesh
!  ############################################################################
!     if (ilinmesh_panshift==1 .and. ishift==0 .and. cellnew%rpan_intervall(cellnew%npan_log+ipan)>cell%rmesh(cell%NRMIN_NS)) then
!       ishift=1
!       npan_logtemp=cellnew%npan_log+ipan
!   
!       cellnew%rpan_intervall(cellnew%npan_log+ipan+1) = cellnew%rpan_intervall(cellnew%npan_log+ipan)
!       cellnew%ipan_intervall(cellnew%npan_log+ipan+1) = (cellnew%npan_log+ipan+ishift)*(cellnew%ncheb+1)
!       cellnew%rpan_intervall(cellnew%npan_log+ipan) = cell%rmesh(cell%NRMIN_NS)
!       cellnew%ipan_intervall(cellnew%npan_log+ipan) = (cellnew%npan_log+ipan)*(cellnew%ncheb+1)
!     end if
!  ############################################################################

else  !( additional panel in the log. part)

!  ############################################################################
  ! creates the panels for the logarithmic region
!  ############################################################################
  do ipan=0,cellnew%npan_log-ilogmesh_panshift
    rval=(fac**(ipan)-1.0D0) / (fac**(cellnew%npan_log-ilogmesh_panshift)-1.0D0)
    cellnew%rpan_intervall(ipan+ishift) = rmin+rval*(rmax-rmin)                       ! set radial value
    cellnew%ipan_intervall(ipan+ishift) = (ipan+ishift)*(cellnew%ncheb+1)             ! and the index
    if (ishift==0 .and. cellnew%rpan_intervall(ipan)>cell%rmesh(cell%NRMIN_NS)) then
      ishift=1                                                                        ! insert the additional panel
      npan_logtemp=ipan                                                               ! 
  
      cellnew%rpan_intervall(ipan+1) = cellnew%rpan_intervall(ipan)
      cellnew%ipan_intervall(ipan+1) = (ipan+ishift)*(cellnew%ncheb+1)
      cellnew%rpan_intervall(ipan) = cell%rmesh(cell%NRMIN_NS)
      cellnew%ipan_intervall(ipan) = ipan*(cellnew%ncheb+1)
    end if
  end do !npan_log
  if (ishift==0 .and. ilogmesh_panshift==1) stop 'error in interpolatecell'

  ishift=0
  rmin= config%rlogpan 
  rmax= cell%rcore

!  ############################################################################
  ! creates the panels for the equidistant region
!  ############################################################################
  do ipan=0,cellnew%npan_eq-ilinmesh_panshift
    cellnew%rpan_intervall(cellnew%npan_log+ipan+ishift) = rmin+ipan*(rmax-rmin)/(cellnew%npan_eq-ilinmesh_panshift)
    cellnew%ipan_intervall(cellnew%npan_log+ipan+ishift) = (cellnew%npan_log+ipan+ishift)*(cellnew%ncheb+1)
  end do !npan_log

  if (ishift==0 .and. ilinmesh_panshift==1) stop 'error in interpolatecell'

!  ############################################################################
  ! creates the panels for the shape function region
!  ############################################################################
  do ipan=1, cellnew%npan_inst
    cellnew%rpan_intervall(cellnew%npan_log+cellnew%npan_eq+ipan)=cell%rmesh(cell%nrcut(ipan+1))
    cellnew%ipan_intervall(cellnew%npan_log+cellnew%npan_eq+ipan) = (cellnew%npan_log+cellnew%npan_eq+ipan)*(cellnew%ncheb+1)
  end do
  cellnew%npan_eq= cellnew%npan_eq+ cellnew%npan_log-npan_logtemp
  cellnew%npan_log=npan_logtemp

end if ! ilinmesh_panshift==1

! do only the first time the routine is called
if (first==1) then
  write(1337,*) '#########################################'
  write(1337,*) '          Panel division'
  write(1337,*) '#########################################'
    write(1337,*) '        panel # |  start n  |     end r   |     panel width '
  ipan=0
  write(1337,*) ipan,cellnew%ipan_intervall(ipan),cellnew%rpan_intervall(ipan)
  do ipan=1,cellnew%npan_tot
    write(1337,*) ipan,cellnew%ipan_intervall(ipan),cellnew%rpan_intervall(ipan),cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1)
  end do
  write(1337,*) '#########################################'
end if

! check for errors
ipan=1
itemp = cellnew%ipan_intervall(ipan)-cellnew%ipan_intervall(ipan-1)
do ipan=1,cellnew%npan_tot
  if (itemp/=cellnew%ipan_intervall(ipan)-cellnew%ipan_intervall(ipan-1) ) then
    print *,'[interpolatecell] error in distributing the panel width'
    print *,'                  panel',1,'points',cellnew%ipan_intervall(1)-cellnew%ipan_intervall(0)
    print *,'                  panel',ipan,'points',cellnew%ipan_intervall(ipan)-cellnew%ipan_intervall(ipan-1)
    stop
  end if
end do


!  ############################################################################
!                    INTERPOLATE POTENTIAL/SHAPEFUNCTION
!  ############################################################################


if (.not. allocated(cellnew%rmeshnew)) then
  allocate(cellnew%rmeshnew(cellnew%nrmaxnew))
end if
allocate(vpotnew(1:cellnew%nrmaxnew,lmvpot))

! creates the radial nodes for npan_tot panels with the boundaries defined in 
! rpan_intervall and saves them into cellnew%rmeshnew
call create_newmesh(cellnew%npan_tot,cellnew%ncheb,cellnew%rpan_intervall,cellnew%rmeshnew)


do lm1=1,lmvpot
  if (cmode/='shapefun') then ! if mode is not shapefun
    ! interpolate for log. region
    do ipan=1,cellnew%npan_log
      irmin=1
      irmax=cell%NRMIN_NS
      irminnew=cellnew%ipan_intervall(ipan-1)+1
      irmaxnew=cellnew%ipan_intervall(ipan)
      if (config_testflag('lsqinterpol')) then
        call interpolpot(cell%rmesh(irmin:irmax),cellnew%rmeshnew(irminnew:irmaxnew), &
                        vpot(irmin:irmax,lm1),vpotnew(irminnew:irmaxnew,lm1),&
                        irmax-irmin+1,irmaxnew-irminnew+1,8)
      else
        call interpolspline(cell%rmesh(irmin:irmax),cellnew%rmeshnew(irminnew:irmaxnew), &
                        vpot(irmin:irmax,lm1),vpotnew(irminnew:irmaxnew,lm1),&
                        irmax-irmin+1,irmaxnew-irminnew+1)
      end if
    end do
    ! interpolate for eq. region
    do ipan=cellnew%npan_log+1,cellnew%npan_log+cellnew%npan_eq
      irmin=cell%NRMIN_NS
      irmax=cell%NRCORE
      irminnew=cellnew%ipan_intervall(ipan-1)+1
      irmaxnew=cellnew%ipan_intervall(ipan)
      if (config_testflag('lsqinterpol')) then
        call interpolpot(cell%rmesh(irmin:irmax),cellnew%rmeshnew(irminnew:irmaxnew), &
                        vpot(irmin:irmax,lm1),vpotnew(irminnew:irmaxnew,lm1),&
                        irmax-irmin+1,irmaxnew-irminnew+1,8)
      else
        call interpolspline(cell%rmesh(irmin:irmax),cellnew%rmeshnew(irminnew:irmaxnew), &
                        vpot(irmin:irmax,lm1),vpotnew(irminnew:irmaxnew,lm1),&
                        irmax-irmin+1,irmaxnew-irminnew+1)
      end if
  
    end do

  end if

  ! interpolate for shapefun region
  ipan2=0
  do ipan=cellnew%npan_log+cellnew%npan_eq+1,cellnew%npan_log+cellnew%npan_eq+cellnew%npan_inst
    ipan2=ipan2+1
    irmin=cell%nrcut(ipan2)+1
    irmax=cell%nrcut(ipan2+1)
    irminnew=cellnew%ipan_intervall(ipan-1)+1
    irmaxnew=cellnew%ipan_intervall(ipan)
    if (config_testflag('lsqinterpol')) then
      call interpolpot(cell%rmesh(irmin:irmax),cellnew%rmeshnew(irminnew:irmaxnew), &
                       vpot(irmin:irmax,lm1),vpotnew(irminnew:irmaxnew,lm1),&
                       irmax-irmin+1,irmaxnew-irminnew+1,8)
    else
      call interpolspline(cell%rmesh(irmin:irmax),cellnew%rmeshnew(irminnew:irmaxnew), &
                       vpot(irmin:irmax,lm1),vpotnew(irminnew:irmaxnew,lm1),&
                       irmax-irmin+1,irmaxnew-irminnew+1)
    end if
  end do
end do

if (cmode=='potential') then
  if (.not. allocated(cellnew%vpotnew) ) then
    allocate(cellnew%vpotnew(1:cellnew%nrmaxnew,lmvpot,nspin))
  end if
  cellnew%vpotnew(:,:,ispin)= vpotnew
else if  (cmode=='shapefun') then
  if (.not. allocated(cellnew%shapefun) ) then
    allocate(cellnew%shapefun(cellnew%nrmaxnew,lmvpot))
  end if  
  cellnew%shapefun(:,:)= vpotnew(:,:)
else if  (cmode=='test') then
  if (.not. allocated(testpot) ) then
    allocate(testpot(cellnew%nrmaxnew,lmvpot))
  end if  
  testpot(:,:)= vpotnew(:,:)
else
  stop '[interpolatecell] cmode not known'
end if

if (config_testflag('rmeshdensity')) then
  open(unit=234918173,file='out_rmeshdensity')
  write(234918173,'(50000E25.14)') (cellnew%rmeshnew(ir+1)-cellnew%rmeshnew(ir),ir=1,cellnew%nrmaxnew-1)
end if

write(1337,*) 'exiting interpolecell'

first=0
end subroutine interpolatecell

!-------------------------------------------------------------------------------
!> Summary: creates the new mesh
!> Author: Who wrote this subroutine
!> Category: new-mesh, radial-grid, kkrimp
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------


subroutine create_newmesh(nintv,ncheb,ri,r)
use nrtype
implicit none
integer :: nintv
integer :: ncheb
double precision :: ri(0:nintv)
double precision :: r(nintv*(ncheb+1))
integer :: i,k,ik
double precision :: tau
      DO I = 1,NINTV
         DO K = 0,NCHEB
             IK = I*NCHEB + I - K
             TAU = DCOS(((2*K+1)*PI)/ (2*(NCHEB+1)))
             TAU = 0.5D0*((RI(I) - RI(I-1))*TAU + RI(I) + RI(I-1))
             R(IK) = TAU 
         END DO
      END DO

end subroutine create_newmesh!(nintv,ncheb,ri,r)


end module mod_interpolatecell
