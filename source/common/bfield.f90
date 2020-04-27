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
    itscf1,bfield_in,bfield_strength,bfield_constr,theta,phi)
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
      bfield%theta(:) = theta(:)
      bfield%phi(:) = phi(:)
      bfield%bfield_strength(:) = bfield_strength(:)
      ! define bfield out of theta, phi and strength
      bfield%bfield(:,:) = bfield_in(:,:)
      bfield%bfield_constr(:,:) = 0.d0
  end subroutine init_bfield
end module mod_bfield





