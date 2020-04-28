!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Module storing the run options and the paramters for bfields and constraining fields
!> Author: Sascha Brinker
!> 
!> 
!> 
!> 
!> 
!------------------------------------------------------------------------------------
module mod_bfield
  !-------------------------------------------------------------------------------
  !> Summary: Type used in t_params to store all relevant information for bfields and constraining fields
  !> Author: Sascha Brinker
  !> Category: communication, KKRhost, bfield
  !> Deprecated: False
  !> 
  !-------------------------------------------------------------------------------
  type :: type_bfield
    logical :: lbfield = .False. ! external magnetic field (turned on via runoption <noncobfield>) non-collinear magnetic field
    logical :: lbfield_constr = .False. ! constraining fields (turned on via runoption <noncobfield>) non-collinear magnetic field
    logical :: lbfield_all = .False. ! apply same field to all atoms (True) or individual fields to each atom
    integer :: ibfield = 0 ! spin (0), orbital (1), spin+orbial (2) fields
    integer :: ibfield_constr = 0 ! type of contraint (0 = torque, 1 = magnetic moment)
    integer :: itscf0 = 0    ! start magnetic field at iteration itscf0
    integer :: itscf1 = 10000 ! stop applying magnetic field after iteration itscf1
    real (kind=dp), dimension (:,:), allocatable :: bfield ! external magnetic field in cartesian coordinates, dimensions (natom,3)
    real (kind=dp), dimension (:), allocatable :: bfield_strength ! absolute value of the external magnetic field, dimensions (natom)
    real (kind=dp), dimension (:,:), allocatable :: bfield_constr ! constraining field in cartesian coordinates, dimensions (natom,3)
    real (kind=dp), dimension (:), allocatable :: theta ! polar angle of the magnetic field
    real (kind=dp), dimension (:), allocatable :: phi   ! azimuthal angle of the magnetic field
  end type type_bfield

  type (type_bfield), save :: bfield

contains
  !-------------------------------------------------------------------------------
  !> Summary: Allocate initial magnetic field parameters to be broadcasted via mpi
  !> Author: Sascha Brinker
  !> Category: memory-management, profiling, KKRhost, bfield
  !> Deprecated: False
  !> Allocate initial parameters to be broadcasted via mpi. allocate arrays, has to
  !> be done after `bcast t_params_scalars` for myrank<>master otherwise are the parameters not set
  !-------------------------------------------------------------------------------
  subroutine init_bfield(natyp,lbfield,lbfield_constr,lbfield_all,ibfield,ibfield_constr,itscf0, &
    itscf1)

      implicit none

      integer, intent(in) :: natyp ! external magnetic field (turned on via runoption <noncobfield>) non-collinear magnetic field
      integer :: i_stat
      
      logical, intent(in) :: lbfield 
      logical, intent(in) :: lbfield_constr 
      logical, intent(in) :: lbfield_all 
      integer, intent(in) :: ibfield 
      integer, intent(in) :: ibfield_constr 
      integer, intent(in) :: itscf0 
      integer, intent(in) :: itscf1 
      real (kind=dp), dimension (natyp,3), intent(in) :: bfield_in ! external magnetic field in cartesian coordinates, dimensions (natom,3)
      real (kind=dp), dimension (natyp),   intent(in) :: bfield_strength ! absolute value of the external magnetic field, dimensions (natom)
      real (kind=dp), dimension (natyp,3), intent(in) :: bfield_constr ! constraining field in cartesian coordinates, dimensions (natom,3)
      real (kind=dp), dimension (natyp),   intent(in) :: theta ! polar angle of the magnetic field
      real (kind=dp), dimension (natyp),   intent(in) :: phi   ! azimuthal angle of the magnetic field

      ! init basic parameters
      bfield%lbfield = lbfield
      bfield%lbfield_constr = lbfield_constr
      bfield%lbfield_all = lbfield_all
      bfield%ibfield = ibfield
      bfield%ibfield_constr = ibfield_constr
      bfield%itscf0 = itscf0
      bfield%itscf1 = itscf1
      ! allocate arrays and add to memory screening routine
      allocate (bfield%theta(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%theta))*kind(bfield%theta), 'bfield%theta', 'init_bfield')
      allocate (bfield%phi(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%phi))*kind(bfield%phi), 'bfield%phi', 'init_bfield')
      allocate (bfield%bfield(natyp,3), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%bfield))*kind(bfield%bfield), 'bfield%bfield', 'init_bfield')
      allocate (bfield%bfield_strength(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%bfield_strength))*kind(bfield%bfield_strength), 'bfield%bfield_strength', 'init_bfield')
      allocate (bfield%bfield_constr(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(bfield%bfield_constr))*kind(bfield%bfield_constr), 'bfield%bfield_constr', 'init_bfield')
      ! init allocated arrays
      bfield%bfield_constr(:,:) = 0.d0
      call read_bfield(bfield,natyp)
  end subroutine init_bfield


  !-------------------------------------------------------------------------------
  !> Summary: Reads the atom-wise magnetic field from bfield.dat
  !> Author: Sascha Brinker
  !> Category: KKRhost, bfield
  !> Deprecated: False
  !> the file has the format:
  !> theta  phi  bfield_strength (in Tesla)
  !-------------------------------------------------------------------------------
  subroutine read_bfield(bfield,natyp)
    implicit none

    integer            , intent(in)     :: natyp
    type(type_bfield)  , intent(inout)  :: bfield
    !local
    integer                 :: iatom,ios,ierror
    character(len=200)      :: string1
    integer,save            :: first=1
    
    if (first==1) then
    
      open(unit=57493215, file='bfield.dat', status='old', iostat=ierror)
      if (ierror/=0) then
        write(*,*) '[read_bfield] bfield file does not exist'
        write(*,*) '              setting all bfields to zero'
        write(*,*) '              disabling magnetic field lbfield = F'
        bfield%lbfield        = .false.
        do iatom=1,natyp
          bfield%theta(iatom)            = 0.0D0
          bfield%phi(iatom)              = 0.0D0
          bfield%bfield_strength(iatom)  = 0.0D0
          bfield%bfield(iatom,:)         = 0.0D0
        end do
        return
      end if
      call read_numbofbfields(natyp)
    
    end if
    
    
    
    write(*,*) '  ###############################################'
    write(*,*) '  non-collinear magnetic fields'
    write(*,*) '  ###############################################'
    write(*,*) '        iatom      theta               phi                  bfield (in Ry   )'
    write(1337,*) 'Angles for bfield'
    write(1337,*) 'iatom  theta                    phi              bfield (in Ry   )'
    bfield%lbfield                = .true.
    do iatom=1,natyp
       string1=this_readline(57493215,ios)
       if (ios/=0) stop '[read_bfield] Error reading atom info2'
       if (ios==-1) stop '[read_bfield] EOF'
       read(string1,*) iatom,bfield%theta(iatom),bfield%phi(iatom),bfield%bfield_strength(iatom)
       write(1337,*) iatom,bfield%theta(iatom),bfield%phi(iatom),bfield%bfield_strength(iatom)
       write(*,*) iatom,bfield%theta(iatom),bfield%phi(iatom),bfield%bfield_strength(iatom)
       bfield%theta(iatom)                  = bfield%theta(iatom)/360.0D0*8.0D0*datan(1.0D0)
       bfield%phi(iatom)                    = bfield%phi(iatom)  /360.0D0*8.0D0*datan(1.0D0)
       bfield%bfield_strength(iatom)        = bfield%bfield_strength(iatom)!/235051.787_dp !conversion from Tesla to Ry
       bfield%bfield(iatom,1)               = bfield%bfield_strength(iatom)*cos(bfield%phi(iatom))*sin(bfield%theta(iatom))
       bfield%bfield(iatom,2)               = bfield%bfield_strength(iatom)*sin(bfield%phi(iatom))*sin(bfield%theta(iatom))
       bfield%bfield(iatom,3)               = bfield%bfield_strength(iatom)*cos(bfield%theta(iatom))
       write(*,*) iatom,bfield%theta(iatom),bfield%phi(iatom),bfield%bfield_strength(iatom),bfield%bfield(iatom,1),bfield%bfield(iatom,2),bfield%bfield(iatom,3)
    end do
    ! close(33952084)
    first=0
  end subroutine read_bfield
  
  subroutine read_numbofbfields(natom)
    implicit none
    integer                 :: ios,linecount
    character(len=200)      :: string1
    
    open(unit=55493215, file='bfield.dat', status='old', iostat=ierror)
    ios=0
    linecount = -1
    do while (ios/=-1)
       string1=this_readline(55493215,ios)
       linecount=linecount+1
    end do
    
    if (natom/=linecount) then
      print *,'[read_bfield] number of angles given in file bfield.dat'
      print *,'             is not correct'
      print *,'linecount',linecount
      print *,'natyp',natyp
      stop
    end if
    close(55493215)
  end subroutine read_numbofbfields
  
  function this_readline(ifile,ios)
    !--------------------------------------------------------
    !--  reads the next line in file unit ifile            --
    !--------------------------------------------------------
    !--  files starting with a dash (#) are treated as     --
    !--  a comment !!!                                     --
    !--  OUTPUT: next line which is not commented out      --
    !--          error variable IOS (should be zero)       --
    !--------------------------------------------------------
    ! input variables
      implicit none
    integer,intent(in)               :: ifile
    integer,intent(out)              :: ios
    ! local variables
    character(len=200)  ::this_readline
    do
      read(unit=ifile,fmt='(A)', iostat=ios) this_readline
      if (ios/=0 .or. this_readline(1:1)/='#') exit
    end do
  end function this_readline

end module mod_bfield





