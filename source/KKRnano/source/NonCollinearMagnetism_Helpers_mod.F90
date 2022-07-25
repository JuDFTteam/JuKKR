module NonCollinearMagnetism_Helpers_mod

  private
  public :: rotatematrix, create_Umatrix, intcheb_cell
  
  contains
  
  subroutine rotatematrix(mat,theta,phi,lmmax,mode)
  ! rotates a matrix in the local frame pointing in
  ! the direction of phi and theta to the global frame
  implicit none
  !interface
  double complex,intent(inout)    ::  mat(2*lmmax,2*lmmax)
  double precision,intent(in)     :: phi
  double precision,intent(in)     :: theta
  integer                         :: lmmax
  integer                         :: mode
  !local
  double complex   :: Umat(2*lmmax,2*lmmax)
  double complex   :: Udeggamat(2*lmmax,2*lmmax)
  double complex   :: mattemp(2*lmmax,2*lmmax)
  !double precision :: matmat_zmzm

  !***********************************************************************
  ! create the rotation matrix:
  !     | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
  !  U= |                                                          |
  !     | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |
  !
  !  Udegga = transpose(complex conjug ( U ) )
  !***********************************************************************


  call create_Umatrix(theta,phi,lmmax,Umat,Udeggamat)
  !***********************************************************************
  ! calculate matrix in the global frame:
  !
  !  t_glob = U * t_loc * Udegga
  !***********************************************************************


  if (mode==0) then ! 'loc->glob'
    call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),mat,2*lmmax,Udeggamat,2*lmmax,(0d0,0d0),mattemp,2*lmmax)
    call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Umat,2*lmmax,mattemp,2*lmmax,(0d0,0d0),mat,2*lmmax)
  elseif (mode==1) then !'glob->loc'
    call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),mat,2*lmmax,Umat,2*lmmax,(0d0,0d0),mattemp,2*lmmax)
    call zgemm('N','N',2*lmmax,2*lmmax,2*lmmax,(1d0,0d0),Udeggamat,2*lmmax,mattemp,2*lmmax,(0d0,0d0),mat,2*lmmax)
  else
    stop '[rotatematrix] mode not known'
  end if
  !  writE(324,'(5000F)') tmat
  ! stop

  end subroutine rotatematrix

  subroutine create_Umatrix(theta,phi,lmmax,Umat,Udeggamat)
  implicit none
  !***********************************************************************
  ! create the rotation matrix:
  !     | cos(theta/2) exp(-i/2 phi)   -sin(theta/2) exp(-i/2 phi) |
  !  U= |                                                          |
  !     | sin(theta/2) exp( i/2 phi)    cos(theta/2) exp( i/2 phi) |
  !
  !  Udegga = transpose(complex conjug ( U ) )
  !***********************************************************************double
  !precision :: phi
  !interface
  double precision,intent(in)     :: phi 
  double precision,intent(in)     :: theta
  integer,intent(in)              :: lmmax
  double complex,intent(out)      :: Umat(2*lmmax,2*lmmax)
  double complex,intent(out)      :: Udeggamat(2*lmmax,2*lmmax)
  !local
  double complex                  :: Umat11,Umat12,Umat21,Umat22
  double complex                  :: Udeggamat11,Udeggamat12,Udeggamat21,Udeggamat22
  integer                         :: ival
  double complex,parameter        :: ci=(0.0D0,1.0D0)
  character*25               :: spinmode

  spinmode='kkr'
  if (spinmode=='regular') then
    Umat11      =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
    Umat12      = -sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
    Umat21      =  sin(theta/2.0D0)*exp( ci/2.0D0*phi)
    Umat22      =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
  else if (spinmode=='kkr') then
    Umat11      =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
    Umat12      =  sin(theta/2.0D0)*exp( ci/2.0D0*phi)
    Umat21      = -sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
    Umat22      =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
  else 
    stop '[create_Umatrix] mode not known'
  end if

  Umat=(0.0D0,0.0D0)
  do ival=1,lmmax
    Umat(      ival,      ival) = Umat11
    Umat(      ival,lmmax+ival) = Umat12
    Umat(lmmax+ival,ival)       = Umat21
    Umat(lmmax+ival,lmmax+ival) = Umat22
  end do

  if (spinmode=='regular') then
  Udeggamat11 =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
  Udeggamat12 =  sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
  Udeggamat21 = -sin(theta/2.0D0)*exp( ci/2.0D0*phi)
  Udeggamat22 =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
  else if (spinmode=='kkr') then
  Udeggamat11 =  cos(theta/2.0D0)*exp(-ci/2.0D0*phi)
  Udeggamat12 = -sin(theta/2.0D0)*exp( ci/2.0D0*phi)
  Udeggamat21 =  sin(theta/2.0D0)*exp(-ci/2.0D0*phi)
  Udeggamat22 =  cos(theta/2.0D0)*exp( ci/2.0D0*phi)
  else 
    stop '[create_Umatrix] mode not known'
  end if



  Udeggamat=(0.0D0,0.0D0)
  do ival=1,lmmax
    Udeggamat(      ival,      ival) = Udeggamat11
    Udeggamat(      ival,lmmax+ival) = Udeggamat12
    Udeggamat(lmmax+ival,ival)       = Udeggamat21
    Udeggamat(lmmax+ival,lmmax+ival) = Udeggamat22
  end do

  end subroutine create_Umatrix

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
    integer           :: irstart,irstop,ipan
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

  end subroutine

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

end module
