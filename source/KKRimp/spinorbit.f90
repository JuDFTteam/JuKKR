!-------------------------------------------------------------------------------
!> Summary: Contains subroutines needed to construct SOC Hamiltonian
!> Author: 
!> Deprecated: false ! this needs to be set to true for deprecated subroutines
!>
!-------------------------------------------------------------------------------
module mod_spinorbit

  private
  public :: spinorbit
  
  contains

  !-------------------------------------------------------------------------------
  !> Summary: Computes SOC Hamiltonian
  !> Author: 
  !> Category: KKRimp, potential, spin-orbit-coupling
  !> Deprecated: false ! this needs to be set to true for deprecated subroutines
  !> 
  !> @note Chebychev mesh is used for the integration using a matrix-matrix multiplication (equation 5.53-5.54 of Bauer, Phd thesis) @endnote
  !-------------------------------------------------------------------------------
  subroutine spinorbit(lmax,zatom,eryd,cellnew,nrmaxnew,nspin,vpotll,theta,phi,ncoll,mode)
    use mod_rotatespinframe, only: rotatematrix
    use type_cellnew, only: cell_typenew
    use mod_mathtools, only: matvec_dmdm
    use mod_chebyshev, only: getclambdacinv
    use mod_physic_params, only: cvlight
    use mod_config, only: config_testflag
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
    double complex   :: lsham(2*(lmax+1)**2,2*(lmax+1)**2)
    integer          :: lmmax
    integer          :: irstart,irstop,ipan,ir
    double precision :: widthfac
    double precision :: vpot(nrmaxnew)
    double precision :: dvpotdr(nrmaxnew)
    double complex   :: rmass,temp !(nrmaxnew)
    double complex   :: hsofac(nrmaxnew)
    double precision :: clambdacinv(0:cellnew%ncheb,0:cellnew%ncheb)
    integer          :: ilm1,ilm2
    double precision :: rnucl,atn

    lmmax=(lmax+1)**2
    vpot=0.0d0

    if (nspin==1) then
      stop 'function spinorbit used but nspin=1'
    else
      do ir=1,nrmaxnew
        vpot(ir)=0.5d0*(cellnew%vpotnew(ir,1,1)+cellnew%vpotnew(ir,1,2))
      end do 
    end if


    !************************************************************************************
    ! derivative of vpot (without the core potential)
    !************************************************************************************

    ! integration is done by a matrix-matrix multiplication using
    ! equation 5.53-5.54 of bauer, phd thesis
    call getclambdacinv(cellnew%ncheb,clambdacinv)

    do ipan=1,cellnew%npan_tot
      irstart=cellnew%ipan_intervall(ipan-1)+1
      irstop = cellnew%ipan_intervall(ipan)
      widthfac = 2.0d0/(cellnew%rpan_intervall(ipan)-cellnew%rpan_intervall(ipan-1))
      dvpotdr(irstart:irstop) = matvec_dmdm(clambdacinv,vpot(irstart:irstop))
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
    ! calucate the l*s operator in the real spherical harmonics basis
    !************************************************************************************

    call spin_orbit_compl(lmax,lmmax,lsham)

    if ( config_testflag('nospinflip') ) then
      lsham((lmax+1)**2+1:,:(lmax+1)**2)=(0.0d0,0.0d0)
      lsham(:(lmax+1)**2,(lmax+1)**2+1:)=(0.0d0,0.0d0)
    end if

    if (ncoll==1) then                                        ! added 2.2.2012
      call rotatematrix(lsham,theta, phi,lmmax,1) !'glob->loc')   ! added 2.2.2012
    end if                                                    ! added 2.2.2012

    if (mode=='conjg') then ! not needed, might be deleted
      do ilm1=1,2*(lmax+1)**2
        do ilm2=1,2*(lmax+1)**2
          lsham(ilm2,ilm1)=conjg(lsham(ilm2,ilm1))
        end do
      end do
    elseif (mode=='transpose') then ! used for the left solution
      do ilm1=1,2*(lmax+1)**2
        do ilm2=1,ilm1-1
          temp            =lsham(ilm2,ilm1)
          lsham(ilm2,ilm1)=lsham(ilm1,ilm2)
          lsham(ilm1,ilm2)=temp
        end do
      end do
    elseif (mode=='transpose2x2') then ! not needed, might be deleted
      do ilm1=1,(lmax+1)**2
        do ilm2=(lmax+1)**2+1,2*(lmax+1)**2
          temp            =lsham(ilm2,ilm1)
          lsham(ilm2,ilm1)=lsham(ilm1,ilm2)
          lsham(ilm1,ilm2)=temp
        end do
      end do
    elseif (mode=='transpose2x2conjg') then ! not needed, might be deleted
      do ilm1=1,(lmax+1)**2
        do ilm2=(lmax+1)**2+1,2*(lmax+1)**2
          temp            =lsham(ilm2,ilm1)
          lsham(ilm2,ilm1)=lsham(ilm1,ilm2)
          lsham(ilm1,ilm2)=temp
        end do
      end do
      do ilm1=1,2*(lmax+1)**2
        do ilm2=1,2*(lmax+1)**2
          lsham(ilm2,ilm1)=conjg(lsham(ilm2,ilm1))
        end do
      end do
    elseif (mode=='1') then
      ! do nothing
    else
      stop'[spinorbit] mode not known'
    end if

    !************************************************************************************
    ! calculate the prefacor of the spin-orbit coupling hamiltonian
    ! and adding the resulting term to the potential
    ! check heers, phd thesis  or bauer, phd thesis for details on the prefactor
    !************************************************************************************
    do ir=1,nrmaxnew
      rmass=0.5d0 - 0.5d0/cvlight**2*((vpot(ir)-eryd) - 2.0d0*zatom/cellnew%rmeshnew(ir))
      hsofac(ir)=1d0/(2d0*rmass**2*cvlight**2*cellnew%rmeshnew(ir))*dvpotdr(ir)
      vpotll(:,:,ir)=vpotll(:,:,ir)+hsofac(ir)*lsham
    end do

  end subroutine spinorbit


  !-------------------------------------------------------------------------------
  !> Summary: Calculates the relativistic mass depending on the radial position
  !> Author: 
  !> Category: KKRimp, potential, physical-observables, spin-orbit-coupling
  !> Deprecated: false ! this needs to be set to true for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine rel_mass(rmass,rm,vpot,zatom,eryd,cvlight,irmd)
    implicit none 
    !interface
    integer          :: irmd
    double precision :: rm(irmd),vpot(irmd)
    double complex   :: eryd
    double precision :: cvlight,zatom
    !output 
    double precision :: rmass(irmd), diffv ! radial derivative of the spherical input potential
    !local 
    integer :: nr

    do nr=1,irmd
      diffv = (vpot(nr)-eryd) - 2.0d0*zatom/rm(nr)
      rmass(nr)=0.5d0 - 0.5d0/cvlight**2*diffv
    end do

  end subroutine rel_mass


  !-------------------------------------------------------------------------------
  !> Summary: Claculate the matrix L*s in the basis of real spherical harmonics
  !> Author: 
  !> Category: KKRimp, special-functions, physical-observables
  !> Deprecated: false ! this needs to be set to true for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine spin_orbit_compl(lmax,lmmaxd,l_s)

    implicit none
    integer, intent(in) :: lmax,lmmaxd

    double complex, intent(out) :: l_s(lmmaxd*2,lmmaxd*2)

    ! local variables 
    integer :: rl,lm1,lm2
    double complex :: icompl
    double complex,allocatable :: ls_l(:,:)


    icompl=(0d0,1d0) 

    l_s=(0.0d0,0.0d0)

    do rl=0,lmax

      allocate(ls_l((2*rl+1)*2,(2*rl+1)*2))
      ls_l=(0.0d0,0.0d0)
      call spin_orbit_one_l(rl,ls_l)

      do lm1=1,(2*rl+1)*2

        if (lm1 <= 2*rl+1 ) then
          do lm2=1,(2*rl+1)
            l_s(rl**2+lm1,rl**2+lm2)=0.5d0*ls_l(lm1,lm2)
          end do 
          do lm2=(2*rl+1)+1,(2*rl+1)*2
            l_s(rl**2+lm1,lmmaxd+rl**2-(2*rl+1)+lm2)= 0.5d0*ls_l(lm1,lm2)
          end do 
        else
          do lm2=1,(2*rl+1)
            l_s(lmmaxd+rl**2-(2*rl+1)+lm1,rl**2+lm2)= 0.5d0*ls_l(lm1,lm2)
          end do 
          do lm2=(2*rl+1)+1,(2*rl+1)*2
            l_s(lmmaxd+rl**2-(2*rl+1)+lm1,lmmaxd+rl**2-(2*rl+1)+lm2)= 0.5d0*ls_l(lm1,lm2)
          end do
        end if

      end do !lm1

      deallocate(ls_l)

    end do !rl=0,lmax

  end subroutine


  !-------------------------------------------------------------------------------
  !> Summary: Calculate the matrix l*s is calculated for the basis of real spherical harmonics for a single l channel
  !> Author: 
  !> Category: KKRimp, spin-orbit-coupling
  !> Deprecated: false ! this needs to be set to true for deprecated subroutines
  !>
  !> schematically the matrix has the form
  !> (  -l_z    l_+  )
  !> (  l_-     l_z  )
  !-------------------------------------------------------------------------------
  subroutine spin_orbit_one_l(lmax,l_s)

    implicit none
    ! i

    integer, intent(in) :: lmax
    double complex, intent(out) :: l_s((2*lmax+1)*2,(2*lmax+1)*2)

    ! c  local variables 
    integer :: i1,i2,i1l 
    double complex :: icompl
    double complex,allocatable :: l_min(:,:)
    double complex,allocatable :: l_up(:,:)
    double precision :: lfac

    icompl=(0d0,1d0) 

    allocate(l_min(-lmax:lmax,-lmax:lmax))
    allocate(l_up(-lmax:lmax,-lmax:lmax))

    !  initialize the matrix 
    do i1=1,(2*lmax+1)*2
      do i2=1,(2*lmax+1)*2
        l_s(i2,i1)=0d0
      end do
    end do

    do i1=-lmax,lmax
      do i2=-lmax,lmax
        l_min(i2,i1)=0d0
        l_up(i2,i1)=0d0
      end do
    end do

    ! fill the second and the forth quadrant with l_z
    ! (-l_z,respectively) 
    do i1=1,2*lmax+1
      i1l=i1-lmax-1       ! the value of m (varies from -l to +l)  
      i2=2*lmax+1-(i1-1)  
      l_s(i2,i1)=-icompl*i1l 
    end do 

    do i1=2*lmax+2,(2*lmax+1)*2
      i1l=i1-lmax-1-(2*lmax+1)       ! the value of m (varies from -l to +l)  
      i2=(2*lmax+1)*2-(i1-(2*lmax+2))  
      l_s(i2,i1)=icompl*i1l 
    end do 

    ! implement now l_- in the third quadrant
    if (lmax>0) then

      lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
      l_min(0,-1)=-icompl*lfac
      ! l_min(0,-1)=icompl*lfac
      l_min(0,1)=lfac
      l_min(-1,0)=icompl*lfac
      l_min(1,0)=-lfac

      if (lmax > 1) then

        do i1=2,lmax

          lfac=0.5d0*sqrt(lmax*(lmax+1d0)-i1*(i1-1d0))
          l_min(-i1,-i1+1)=-lfac
          l_min(-i1,i1-1)=icompl*lfac
          l_min(i1,-i1+1)=-icompl*lfac
          l_min(i1,i1-1)=-lfac

          lfac=0.5d0*sqrt(lmax*(lmax+1d0)-(i1-1d0)*i1)
          l_min(-i1+1,-i1)=lfac
          l_min(-i1+1,i1)=icompl*lfac
          l_min(i1-1,-i1)=-icompl*lfac
          l_min(i1-1,i1)=lfac

        end do

      end if
    end if

    do i1=-lmax,lmax
      do i2=-lmax,lmax
        ! l_s(i2+lmax+1,i1+3*lmax+2)=l_min(i2,i1)
        ! transpose l_min          
        ! l_s(i2+lmax+1,i1+3*lmax+2)=l_min(i1,i2)
        l_s(i2+3*lmax+2,i1+lmax+1)=l_min(i1,i2)
      end do
    end do

    ! implement now l_+ in the   quadrant
    if (lmax>0) then

      lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
      l_up(0,-1)=-icompl*lfac
      l_up(0,1)=-lfac
      l_up(-1,0)=icompl*lfac
      l_up(1,0)=lfac

      if (lmax > 1) then

        do i1=2,lmax

          lfac=0.5d0*sqrt(lmax*(lmax+1d0)-i1*(i1-1d0))
          l_up(-i1,-i1+1)=lfac
          l_up(-i1,i1-1)=icompl*lfac
          l_up(i1,-i1+1)=-icompl*lfac
          l_up(i1,i1-1)=lfac

          lfac=0.5d0*sqrt(lmax*(lmax+1d0)-(i1-1d0)*i1)
          l_up(-i1+1,-i1)=-lfac
          l_up(-i1+1,i1)=icompl*lfac
          l_up(i1-1,-i1)=-icompl*lfac
          l_up(i1-1,i1)=-lfac

        end do

      end if
    end if

    do i1=-lmax,lmax
      do i2=-lmax,lmax
      ! l_s(i2+3*lmax+2,i1+lmax+1)=l_up(i2,i1)
      ! transpose l_up          
        l_s(i2+lmax+1,i1+3*lmax+2)=l_up(i1,i2)
      ! l_s(i2+3*lmax+2,i1+lmax+1)=l_up(i1,i2)
      end do
    end do

    deallocate(l_min)
    deallocate(l_up)

  end subroutine


  !-------------------------------------------------------------------------------
  !> Summary: Old form of spinorbit matrix (unused)
  !> Author: 
  !> Category: KKRimp, spin-orbit-coupling
  !> Deprecated: True ! this needs to be set to true for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine spin_orbit_one_l_old(lmax,l_s)
    ! here the 1x1 block is still spin up
    ! and the 2x2 block is spin down ( swantje's convention )
    ! not used anymore
    implicit none
    ! ************************************************************************
    !      in this subroutine the matrix l*s is calculated for the basis of
    !      real spherical harmonics 
    !
    !      schematically it has the form
    !      (  l_z    l_-  )
    !      (  l_+   -l_z  )
    !

    integer, intent(in) :: lmax
    double complex, intent(out) :: l_s((2*lmax+1)*2,(2*lmax+1)*2)

    ! c  local variables 
    integer :: i1,i2,i1l 
    double complex :: icompl
    double complex,allocatable :: l_min(:,:)
    double complex,allocatable :: l_up(:,:)
    double precision :: lfac


    icompl=(0d0,1d0) 

    allocate(l_min(-lmax:lmax,-lmax:lmax))
    allocate(l_up(-lmax:lmax,-lmax:lmax))

    ! initialize the matrix 
    do i1=1,(2*lmax+1)*2
      do i2=1,(2*lmax+1)*2
        l_s(i2,i1)=0d0
      end do
    end do

    do i1=-lmax,lmax
      do i2=-lmax,lmax
        l_min(i2,i1)=0d0
        l_up(i2,i1)=0d0
      end do
    end do

    ! fill the second and the forth quadrant with l_z
    ! (-l_z,respectively) 
    do i1=1,2*lmax+1
      i1l=i1-lmax-1       ! the value of m (varies from -l to +l)  
      i2=2*lmax+1-(i1-1)  
      l_s(i2,i1)=icompl*i1l 
    end do 

    do i1=2*lmax+2,(2*lmax+1)*2
      i1l=i1-lmax-1-(2*lmax+1)       ! the value of m (varies from -l to +l)  
      i2=(2*lmax+1)*2-(i1-(2*lmax+2))  
      l_s(i2,i1)=-icompl*i1l 
    end do 

    !  implement now l_- in the third quadrant
    if (lmax>0) then

      lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
      l_min(0,-1)=-icompl*lfac
      !         l_min(0,-1)=icompl*lfac
      l_min(0,1)=lfac
      l_min(-1,0)=icompl*lfac
      l_min(1,0)=-lfac

      if (lmax > 1) then

        do i1=2,lmax

          lfac=0.5d0*sqrt(lmax*(lmax+1d0)-i1*(i1-1d0))
          l_min(-i1,-i1+1)=-lfac
          l_min(-i1,i1-1)=icompl*lfac
          l_min(i1,-i1+1)=-icompl*lfac
          l_min(i1,i1-1)=-lfac

          lfac=0.5d0*sqrt(lmax*(lmax+1d0)-(i1-1)*(i1))
          l_min(-i1+1,-i1)=lfac
          l_min(-i1+1,i1)=icompl*lfac
          l_min(i1-1,-i1)=-icompl*lfac
          l_min(i1-1,i1)=lfac

        end do

      end if
    end if

    do i1=-lmax,lmax
      do i2=-lmax,lmax
        !  l_s(i2+lmax+1,i1+3*lmax+2)=l_min(i2,i1)
        ! transpose l_min          
        l_s(i2+lmax+1,i1+3*lmax+2)=l_min(i1,i2)
      end do
    end do


    ! implement now l_+ in the   quadrant
    if (lmax>0) then

      lfac=sqrt(lmax*(lmax+1d0))/sqrt(2d0)
      l_up(0,-1)=-icompl*lfac
      l_up(0,1)=-lfac
      l_up(-1,0)=icompl*lfac
      l_up(1,0)=lfac

      if (lmax > 1) then

        do i1=2,lmax

          lfac=0.5d0*sqrt(lmax*(lmax+1d0)-i1*(i1-1d0))
          l_up(-i1,-i1+1)=lfac
          l_up(-i1,i1-1)=icompl*lfac
          l_up(i1,-i1+1)=-icompl*lfac
          l_up(i1,i1-1)=lfac

          lfac=0.5d0*sqrt(lmax*(lmax+1d0)-(i1-1)*(i1))
          l_up(-i1+1,-i1)=-lfac
          l_up(-i1+1,i1)=icompl*lfac
          l_up(i1-1,-i1)=-icompl*lfac
          l_up(i1-1,i1)=-lfac

        end do

      end if
    end if
    
    do i1=-lmax,lmax
      do i2=-lmax,lmax        
        l_s(i2+3*lmax+2,i1+lmax+1)=l_up(i1,i2)
      end do
    end do

    deallocate(l_min)
    deallocate(l_up)

  end subroutine


  !-------------------------------------------------------------------------------
  !> Summary: Transpose double complex matrix
  !> Author: 
  !> Category: KKRimp, numerical-tools
  !> Deprecated: True ! this needs to be set to true for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
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
