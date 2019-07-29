#ifdef TASKLOCAL_FILES
#define FILEWRITE(X,Y) write(X)
#define FILEREAD(X,Y) read(X)
#else
#define FILEWRITE(X,Y) write(X,Y)
#define FILEREAD(X,Y) read(X,Y)
#endif

module ChebMeshData_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: ChebMeshData, create, destroy, load, construct
  public :: getMinReclenChebMesh
 ! public :: initChebMesh, initMuffinTinMesh, initInterstitialMesh,
  public :: writeChebMeshDataDA           
  public :: readChebMeshDataDA, openChebMeshDataDAFile, closeChebMeshDataDAFile, writeChebMeshDataIndexDA      
  public :: readChebMeshDataIndexDA, openChebMeshDataIndexDAFile, closeChebMeshDataIndexDAFile 
  public :: interpolate_poten  
 ! public :: readChebMeshDataHeader    
  
  type ChebMeshData
  
  integer                       :: ncheb               ! Number of points per panel in new mesh
  integer                       :: irmd_new            ! Number of radial mesh points in new mesh
  integer                       :: npan_lognew         ! Number of intervals from nucleus to [R_LOG] 
  integer                       :: npan_eqnew          ! Number of intervals from [R_LOG] to muffin-tin radius
  integer                       :: npan_tot            ! Number of overall intervals 
  integer, allocatable          :: ipan_intervall(:)   ! Global index of starting point of a certain panel 
  double precision, allocatable :: rpan_intervall(:)   ! R value of starting point of a certain panel
  double precision, allocatable :: rnew(:)             ! New radial mesh points 
  double precision, allocatable :: thetasnew(:,:)      ! New shapefunctions in spherical harmonics expansion
  
  endtype
  
  integer, parameter, private :: MAGIC_NUMBER = -889271554
  
  interface create
    module procedure createChebMeshData
  endinterface

  interface load
    module procedure createChebMeshDataFromFile
  endinterface
  
  interface construct
    module procedure ConstructChebMesh
  endinterface
  
  interface destroy
    module procedure destroyChebMeshData
  endinterface

contains

!----------------------------------------------------------------------------
subroutine createChebMeshData(self, irmd_new, npan_tot, nfun, params)
    use InputParams_mod

    type(ChebMeshData), intent(inout)      :: self
    integer, intent(in)                :: irmd_new ! Number of radial mesh points in new mesh
    integer, intent(in)                :: npan_tot ! Number of overall intervals 
    integer, intent(in)                :: nfun     ! Number of non-zero shape function components
    type(InputParams), intent(in)      :: params
    
    self%ncheb    = params%ncheb
    self%npan_tot = npan_tot
    self%irmd_new = irmd_new

    allocate(self%RNEW(self%irmd_new))
    allocate(self%RPAN_INTERVALL(0:self%npan_tot))
    allocate(self%IPAN_INTERVALL(0:self%npan_tot))
    allocate(self%THETASNEW(self%irmd_new,nfun))

    self%NPAN_LOGNEW     = 0
    self%NPAN_EQNEW      = 0
    self%IPAN_INTERVALL  = 0

    self%RNEW            = 0.0d0
    self%RPAN_INTERVALL  = 0.0d0
    self%THETASNEW       = 0.0d0
    
endsubroutine ! create

subroutine destroyChebMeshData(self)
    type(ChebMeshData), intent(inout) :: self
    
    deallocate(self%RNEW)
    deallocate(self%RPAN_INTERVALL)
    deallocate(self%IPAN_INTERVALL)
    deallocate(self%THETASNEW)

endsubroutine ! destroy

!----------------------------------------------------------------------------
!> Write mesh data to direct access file 'fileunit' at record 'recnr'
subroutine writeChebMeshDataDA(meshdata, fileunit, recnr)
    type(ChebMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr


#ifdef TASKLOCAL_FILES
    character(len=16) :: num 
    write(unit=num, fmt='(a,i7.7)') "cheb_mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='write')

    call writeChebMeshDataIndexDA(meshdata, fileunit, recnr, 0)
#endif

    FILEWRITE (fileunit, rec=recnr) MAGIC_NUMBER + recnr, &
                                    meshdata%NPAN_LOGNEW, &
                                    meshdata%NPAN_EQNEW, &
                                    meshdata%RNEW(:), &
                                    meshdata%RPAN_INTERVALL(0:), &
                                    meshdata%IPAN_INTERVALL(0:), &
                                    meshdata%THETASNEW(:,:), &
                                    MAGIC_NUMBER + recnr
#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif
endsubroutine ! write

!===========  Index file ======================================================

  !----------------------------------------------------------------------------
  !> Write mesh dimension data to direct access file 'fileunit' at record
  !'recnr'
subroutine writeChebMeshDataIndexDA(meshdata, fileunit, recnr, max_reclen)
    type(ChebMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(in) :: max_reclen

    FILEWRITE (fileunit, rec=recnr) meshdata%npan_tot*(meshdata%ncheb+1), & ! #points in new radial mesh
                                    meshdata%npan_tot, &                  ! #panels in new radial mesh
                                    max_reclen, &
                                    MAGIC_NUMBER + recnr
endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Opens ChebMeshData index file.
subroutine openChebMeshDataIndexDAFile(meshdata, fileunit, filename)
    type(ChebMeshData), intent(in) :: meshdata
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

#ifndef TASKLOCAL_FILES
    integer :: reclen, max_reclen

    max_reclen = 0
    inquire (iolength=reclen) meshdata%npan_tot*(meshdata%ncheb+1), &
                              meshdata%npan_tot, &
                              max_reclen, &
                              MAGIC_NUMBER

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')
#endif
endsubroutine ! open

!> Opens ChebMeshData direct access file.
subroutine openChebMeshDataDAFile(fileunit, filename, max_reclen)
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: max_reclen

#ifndef TASKLOCAL_FILES
    integer :: reclen

      reclen = max_reclen

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')
#endif
endsubroutine ! open

  !---------------------------------------------------------------------------
  !> Return MINIMUM record length needed to store mesh of this atom
  !>  
  !> Note: for a file containing meshes of several atoms, the maximum
  !> of all their record lengths has to be determined
integer function getMinReclenChebMesh(meshdata) result(reclen)
    type(ChebMeshData), intent(in) :: meshdata

    inquire (iolength=reclen) MAGIC_NUMBER, &
                              meshdata%NPAN_LOGNEW, &
                              meshdata%NPAN_EQNEW, &
                              meshdata%RNEW(:), &
                              meshdata%RPAN_INTERVALL(0:), &
                              meshdata%IPAN_INTERVALL(0:), &
                              meshdata%THETASNEW(:,:), &
                              MAGIC_NUMBER

  endfunction ! get

  !----------------------------------------------------------------------------
  !> Read and create radial mesh data from files <filename> and <filename>.idx
  !>  
  !> The index file with extension .idx is necessary because different
  !> atoms can have a different number of radial mesh points
subroutine createChebMeshDataFromFile(meshdata, filename, recnr, nfun, params)
    use InputParams_mod
    
    type(ChebMeshData), intent(inout) :: meshdata
    character(len=*), intent(in) :: filename
    integer, intent(in) :: recnr
    integer, intent(in) :: nfun
    type(InputParams), intent(in) :: params
    integer, parameter :: fu = 37
    integer :: irmd_new, ntotd, max_reclen

    ! index file has extension .idx but the filename is ignored for task-local
    ! files
    call openChebMeshDataIndexDAFile(meshdata, fu, filename-".idx")
    call readChebMeshDataIndexDA(fu, recnr, irmd_new, ntotd, max_reclen)
    call closeChebMeshDataIndexDAFile(fu)

    call createChebMeshData(meshdata, irmd_new, ntotd, nfun, params)

    call openChebMeshDataDAFile(fu, filename, max_reclen)
    call readChebMeshDataDA(meshdata, fu, recnr)
    call closeChebMeshDataDAFile(fu)

endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Read mesh dimension data from direct access file 'fileunit' at record
  !'recnr'
  !>  
  !> Returns dimensions irmd and ipand
subroutine readChebMeshDataIndexDA(fileunit, recnr, irmd_new, ntotd, max_reclen)
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr
    integer, intent(out) :: irmd_new
    integer, intent(out) :: ntotd
    integer, intent(out) :: max_reclen

#ifdef TASKLOCAL_FILES
    character(len=16) :: num 

    write(unit=num, fmt='(a,i7.7)') "cheb_mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='read', status='old')
#endif

    call readChebMeshDataHeader(fileunit, recnr, irmd_new, ntotd, max_reclen)

#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif
endsubroutine ! read

subroutine readChebMeshDataHeader(fileunit, recnr, irmd_new, ntotd, max_reclen)
    integer, intent(in) :: fileunit, recnr
    integer, intent(out) :: irmd_new, ntotd, max_reclen

    integer :: magic, checkmagic

    checkmagic = MAGIC_NUMBER + recnr

    FILEREAD  (fileunit, rec=recnr) irmd_new, &
                                    ntotd, &
                                    max_reclen, &
                                    magic

    if (magic /= checkmagic) die_here("Invalid mesh index data read.")
endsubroutine ! read

!----------------------------------------------------------------------------
!> Closes ChebMeshData direct access file.
subroutine closeChebMeshDataDAFile(fileunit)
    integer, intent(in) :: fileunit
#ifndef TASKLOCAL_FILES
    close(fileunit)
#endif
endsubroutine ! close

!----------------------------------------------------------------------------
!> Closes ChebMeshData index file.
subroutine closeChebMeshDataIndexDAFile(fileunit)
    integer, intent(in) :: fileunit

#ifndef TASKLOCAL_FILES
    close(fileunit)
#endif

endsubroutine ! close

  !----------------------------------------------------------------------------
  !> Read mesh data from direct access file 'fileunit' at record 'recnr'
subroutine readChebMeshDataDA(meshdata, fileunit, recnr)
    type(ChebMeshData), intent(inout) :: meshdata
    integer, intent(in) :: fileunit, recnr

    integer :: magic, magic2, checkmagic

#ifdef TASKLOCAL_FILES
    character(len=16) :: num 
    integer :: ipand, irmd, max_reclen    

    write(unit=num, fmt='(a,i7.7)') "cheb_mesh.",recnr
    open(fileunit, file=num, form='unformatted', action='read', status='old')

    ! read header at beginning of file
    call readChebMeshDataHeader(meshdata, fileunit, recnr, irmd, ipand, max_reclen)
#endif

    checkmagic = MAGIC_NUMBER + recnr

    FILEREAD (fileunit, rec=recnr)  magic, &
                                    meshdata%NPAN_LOGNEW, &
                                    meshdata%NPAN_EQNEW, &
                                    meshdata%RNEW(:), &
                                    meshdata%RPAN_INTERVALL(0:), &
                                    meshdata%IPAN_INTERVALL(0:), &
                                    meshdata%THETASNEW(:,:), &
                                    magic2
#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif

    if (magic /= checkmagic .or. magic2 /= checkmagic) die_here("Invalid mesh data read.")

    
    !call test_mesh_consistency(self)
endsubroutine ! read


  !----------------------------------------------------------------------------
  !> Construct Chebyshev mesh
subroutine ConstructChebMesh(r_log,npan_log,npan_eq,ncheb,  &
                              npan_lognew,npan_eqnew,  &
                              npan_tot,rnew,rpan_intervall,ipan_intervall,  &
                              thetasnew,thetas,nfu,radial_mesh) ! new parameters

!use read_formatted_shapefun_mod, only: shapefunfile
use RadialMeshData_mod, only: RadialMeshData

double precision, intent(in)             :: r_log 
integer, intent(in)                      :: npan_log
integer, intent(in)                      :: npan_eq
integer, intent(in)                      :: ncheb
integer, intent(out)                     :: npan_lognew
integer, intent(out)                     :: npan_eqnew
integer, intent(out)                     :: npan_tot
double precision, intent(out)            :: rnew(:)
double precision, intent(out)            :: rpan_intervall(0:)
integer, intent(out)                     :: ipan_intervall(0:)
double precision, intent(out)            :: thetasnew(:,:)
double precision, intent(in)             :: thetas(:,:)
integer, intent(in)                      :: nfu
type(RadialMeshData), intent(in)         :: radial_mesh

double precision, parameter :: fac=2d0
integer :: ipotm,ir2,ip,  &
    ishift,ilogpanshift,ilinpanshift,npan_logtemp,npan_inst,imin,imax,iminnew,imaxnew,lm1
double precision :: rmin,rmax,rval

thetasnew=0d0
ipotm=0

  npan_inst = radial_mesh%ipan-1
  npan_tot  = npan_log+npan_eq+npan_inst
  
  
! log panel
  rmin=radial_mesh%r(2)
  rmax=r_log
  rval=0d0
  ishift=0
  if (r_log > radial_mesh%r(radial_mesh%irmin)) then
    ilogpanshift=1
    ilinpanshift=0
  else
    ilogpanshift=0
    ilinpanshift=1
  end if
  
  if (ilinpanshift == 1) then
    stop 'non-spherical part of the potential needs to be inside the log panel'
  end if
  
  do ip=0,npan_log-ilogpanshift
    rval=(fac**ip-1d0)/(fac**(npan_log-ilogpanshift)-1d0)
    rpan_intervall(ip+ishift)= rmin+rval*(rmax-rmin)
    ipan_intervall(ip+ishift)= (ip+ishift)*(ncheb+1)
    if (ishift == 0.and. rpan_intervall(ip) > radial_mesh%r(radial_mesh%irmin)) then
      ishift=1
      npan_logtemp=ip
      rpan_intervall(ip+1)=rpan_intervall(ip)
      ipan_intervall(ip+1)=(ip+ishift)*(ncheb+1)
      rpan_intervall(ip)=radial_mesh%r(radial_mesh%irmin)
      ipan_intervall(ip)=ip*(ncheb+1)
    end if
  end do ! npan_log
  
! equivalent panel
  ishift=0
  rmin=r_log
  rmax=radial_mesh%r(radial_mesh%ircut(1))
  do ip=0,npan_eq-ilinpanshift
    rpan_intervall(ip+ishift+npan_log)=rmin+ip*(rmax-rmin)/  &
        (npan_eq-ilinpanshift)
    ipan_intervall(ip+ishift+npan_log)=(npan_log+ip+ishift)* (ncheb+1)
  end do ! npan_eq
  
! intersection zone
  do ip=1,npan_inst
    rpan_intervall(npan_log+npan_eq+ip)=radial_mesh%r(radial_mesh%ircut(ip+1))
    ipan_intervall(npan_log+npan_eq+ip)=(npan_log+npan_eq+ip)* (ncheb+1)
  end do ! npan_inst
  
  npan_eqnew=npan_eq+npan_log-npan_logtemp
  npan_lognew=npan_logtemp
  
  call chebmesh(npan_tot,ncheb,rpan_intervall(0:), rnew(1:))
  
! interpolate shape function thetas to new shape function thetasnew
! save thetas to thetasin
!  icell = i1 ! #shapefunctions = #atoms in kkrnano, jm: icell = ntcell(i1)
  do lm1=1,nfu
!    thetasin(:,lm1,icell)=sfile%shapes(icell)%thetas(:,lm1) ! get thetas for specific atom and lm1
    ir2=0
    do ip=npan_lognew+npan_eqnew+1,npan_tot
      ir2=ir2+1
      imin=radial_mesh%ircut(ir2)+1
      imax=radial_mesh%ircut(ir2+1)
      iminnew=ipan_intervall(ip-1)+1
      imaxnew=ipan_intervall(ip)
      call interpolspline(radial_mesh%r(imin:imax),rnew(iminnew:imaxnew),  &
           thetas(imin-radial_mesh%ircut(1):imax-radial_mesh%ircut(1),lm1),  &
           thetasnew(iminnew:imaxnew,lm1), imax-imin+1,imaxnew-iminnew+1)
    end do
  end do

end subroutine ConstructChebMesh


  !----------------------------------------------------------------------------
  !> Interpolate potential to Chebyshev mesh
subroutine interpolate_poten(nspin,r,irmin,irws,ircut,vins,  &
        visp,npan_log,npan_eq,  &
        npan_tot,rnew,  &
        ipan_intervall,vinsnew, &
        irmd, lmpotd) ! new parameters
 
implicit none

integer, intent(in)                      :: nspin
double precision, intent(in)             :: r(:)
integer, intent(in)                      :: irmin
integer, intent(in)                      :: irws
integer, intent(in)                      :: ircut(0:)
double precision, intent(in)             :: vins(irmin:,:,:)
double precision, intent(in)             :: visp(:,:)
integer, intent(in)                      :: npan_log
integer, intent(in)                      :: npan_eq
integer, intent(in)                      :: npan_tot
double precision, intent(in)             :: rnew(:)
integer, intent(in)                      :: ipan_intervall(0:)
double precision, intent(out)            :: vinsnew(:,:,:)
integer, intent(in)                      :: irmd
integer, intent(in)                      :: lmpotd



double precision, allocatable :: vinsin(:,:,:)

integer :: ipot,ipotm,imin,imax,ip,ir,lm1,ispin,iminnew,imaxnew,ir2

vinsnew=0D0
ipotm=0
allocate(vinsin(irmd,lmpotd,nspin))! irmd=#points in old mesh

ipot=1
  
! save input potential to VINSIN
  vinsin=0D0
  do ir=1,irws
    if (ir < irmin) then
      vinsin(ir,1,1)=visp(ir,ipot)
      vinsin(ir,1,nspin)=visp(ir,ipot+nspin-1)
    else
      do lm1=1,lmpotd
        if (lm1 == 1) then
          vinsin(ir,lm1,1)=visp(ir,ipot)
          vinsin(ir,lm1,nspin)=visp(ir,ipot+nspin-1)
        else
          vinsin(ir,lm1,1)=vins(ir,lm1,ipot)
          vinsin(ir,lm1,nspin)=vins(ir,lm1,ipot+nspin-1)
        end if
      end do
    end if
  end do
  
  do ispin=1,nspin
    
    ipotm=ipotm+1
    
    do lm1=1,lmpotd
      
      imin=1
      imax=irmin
      do ip=1,npan_log
        iminnew=ipan_intervall(ip-1)+1
        imaxnew=ipan_intervall(ip)
        call interpolspline(r(imin:imax),rnew(iminnew:imaxnew),  &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm),  &
            imax-imin+1,imaxnew-iminnew+1)
      end do
      
      imin=irmin
      imax=ircut(1)
      do ip=npan_log+1,npan_log+npan_eq
        iminnew=ipan_intervall(ip-1)+1
        imaxnew=ipan_intervall(ip)
        call interpolspline(r(imin:imax),rnew(iminnew:imaxnew),  &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm),  &
            imax-imin+1,imaxnew-iminnew+1)
      end do
      
      ir2=0
      do ip=npan_log+npan_eq+1,npan_tot
        ir2=ir2+1
        imin=ircut(ir2)+1
        imax=ircut(ir2+1)
        iminnew=ipan_intervall(ip-1)+1
        imaxnew=ipan_intervall(ip)
        call interpolspline(r(imin:imax),rnew(iminnew:imaxnew),  &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm),  &
            imax-imin+1,imaxnew-iminnew+1)
      end do
    end do ! lm1
  end do ! ispin
end subroutine interpolate_poten

  !----------------------------------------------------------------------------
  !> Helping routine for ConstructChebMesh 
subroutine chebmesh(npan,ncheb,ri,ro)

integer, intent(in)                      :: npan
integer, intent(in)                      :: ncheb
double precision, intent(in)             :: ri(0:npan)
double precision, intent(out)            :: ro(npan*(ncheb+1))
!implicit none



integer :: i,k,ik
double precision :: tau,pi

pi=4d0*datan(1d0)
do i=1,npan
  do k=0,ncheb
    ik=i*ncheb+i-k
    tau=dcos(((2*k+1)*pi)/(2*(ncheb+1)))
    tau=0.5d0*((ri(i)-ri(i-1))*tau+ri(i)+ri(i-1))
    ro(ik)=tau
  end do
end do
end subroutine chebmesh

  !----------------------------------------------------------------------------
  !> Helping routine for ConstructChebMesh

subroutine spline_real(nmax,x,y,n,yp1,ypn,y2)

implicit none

integer, intent(in)                      :: nmax
real*8, intent(in)                       :: x(nmax)
real*8, intent(in)                       :: y(nmax)
integer, intent(in)                      :: n
real*8, intent(in out)                   :: yp1
real*8, intent(in out)                   :: ypn
real*8, intent(out)                      :: y2(nmax)


!      REAL*8           x(NMAX)
!      COMPLEX          yp1,ypn,y(NMAX),y2(NMAX)
! Given arrays x(1:n) and  y(1:n) containing a tabulated function,
! i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn
! for the 1rst derivative of the interpolating function at points
! 1 and n, respectively, this routine returns an array y2(1:n) of
! length n which contains the second derivatives of the interpolating
! function at the tabulated points xi.
! If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
! signaled to set the corresponding boundary condition for a natural
! spline, with zero second derivative on that boundary.
! Parameter: NMAX is the largest anticipated value of n.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
integer :: i,k
real*8          p,qn,sig,un,u(nmax)
!      complex          p,qn,sig,un,u(nmax)

if (n > nmax) stop 'spline: n > nmax.'
if (abs(yp1) > 0.99d30) then
! the lower boundary condition is set either to be "natural"
  y2(1) = 0.d0
  u(1) = 0.d0
else
! or else to have a specified first derivative.
  y2(1) = -0.5d0
  u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
end if

do i = 2,n-1
! this is the decomposition loop of the tridiagonal algorithm. y2 and u
! are used for temporary storage of the decomposed factors.
  sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
  p = sig * y2(i-1) + 2.d0
  y2(i) = (sig-1.d0)/p
  u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))  &
      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1)) / p
end do

if (abs(ypn) > 0.99d30) then
! the upper boundary condition is set either to be "natural"
  qn = 0.d0
  un = 0.d0
else
! or else to have a specified 1rst derivative.
  qn = 0.5d0
  un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
end if
y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.d0)
do k = n-1,1,-1
! this is the backsubstitution loop of the tridiagonal algorithm.
  y2(k)=y2(k)*y2(k+1)+u(k)
end do

return
end subroutine spline_real
 

  !----------------------------------------------------------------------------
  !> Helping routine for ConstructChebMesh
subroutine splint_real(xa,ya,y2a,n,x,y,yderiv)

implicit none

real*8, intent(in)                       :: xa(*)
real*8, intent(in)                       :: ya(*)
real*8, intent(in out)                   :: y2a(*)
integer, intent(in)                      :: n
real*8, intent(in)                       :: x
real*8, intent(out)                      :: y
real*8, intent(out)                      :: yderiv


! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the xai's in order), and given the array y2a(1:n), which
! is the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y and the derivative yderiv.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
integer :: k,khi,klo
real*8         a,b,h
! We will  nd the right place in the table by means of bisection.
! This is optimal if sequential calls to this routine are at random
! values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of
! klo and khi and test if they remain appropriate on the
! next call.
klo=1
khi=n
1    if (khi-klo > 1) then
  k=(khi+klo)/2
  if(xa(k) > x)then
    khi=k
  else
    klo=k
  end if
  go to 1
end if
! klo and khi now bracket the input value of x.
h=xa(khi)-xa(klo)
! the xa's must be distinct.
if (h == 0.d0) STOP 'bad xa input in splint' !! used to be PAUSE
! cubic spline polynomial is now evaluated.
a = (xa(khi)-x)/h
b = (x-xa(klo))/h
y = a*ya(klo) + b*ya(khi) +  &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * (h**2)/6.d0
yderiv = (ya(khi)-ya(klo))/h -  &
    ((3.d0*a*a-1.d0)*y2a(klo) - (3.d0*b*b-1.d0)*y2a(khi))*h/6.d0

return
end subroutine splint_real
  
  !----------------------------------------------------------------------------
  !> Helping routine for ConstructChebMesh
subroutine interpolspline(rmesh,rmeshnew,vpot,vpotnew,  &
        nrmax,nrmaxnew)
 
! code converted using to_f90 by alan miller
! date: 2016-04-01  time: 12:05:24
 
implicit none
!interface
integer :: nrmax
integer :: nrmaxnew
double precision :: rmesh(nrmax)
double precision :: rmeshnew(nrmaxnew)
double precision :: vpot(nrmax)
double precision :: vpotnew(nrmaxnew)
!local
double precision :: maxa
double precision :: spline(nrmax)
double precision :: parsum, parsumderiv,r0
integer :: ir
maxa = 1.d35
call spline_real(nrmax,rmesh,vpot,nrmax,maxa,maxa,spline)
!           call spline(irmdjj,r,vm2z,nr,maxa,maxa,vm2zb)

do ir = 1,nrmaxnew
  r0 = rmeshnew(ir)
  call splint_real(rmesh,vpot,spline,nrmax,r0,parsum,parsumderiv)
  vpotnew(ir) = parsum
end do
end subroutine  interpolspline

end module ! ChebMeshData_mod
