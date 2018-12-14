!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_symmetries
! This module provides routines regardinf the symmetry of the lattice and the irreducible brillouin zone.


  use mod_mympi, only: nranks, myrank, master
  use mod_pointgrp, only: pointgrp
  use mod_findgroup, only: findgroup
  implicit none

  private
  public :: symmetries_type, set_symmetries, get_IBZwedge_faces, points_in_wedge, singlepoint_in_wedge, rotate_kpoints, expand_visarrays, expand_areas, expand_spinvalues, expand_torqvalues, unfold_visarrays, get_2DIBZwedge_lines

    type :: symmetries_TYPE

      integer :: N = 5

      integer              :: nsym_used
      integer              :: nsym_found
      double precision     :: rotmat(64,3,3)
      character(len=10)    :: rotname(64)
      integer, allocatable :: isym_used(:)
      integer, allocatable :: isym_found(:)
      integer, allocatable :: isym_notused(:)

    end type symmetries_TYPE 


contains


  subroutine unfold_visarrays(nkpts_all, nkpts_irr, kpt2irr_in, kpoints_in, kpt2irr_out, irr2kpt_out, kpoints_out)
    implicit none
    integer,                       intent(in)  :: nkpts_irr, nkpts_all
    integer,                       intent(in)  :: kpt2irr_in(nkpts_all)
    double precision,              intent(in)  :: kpoints_in(3,nkpts_irr)
    integer,          allocatable, intent(out) :: irr2kpt_out(:), kpt2irr_out(:)
    double precision, allocatable, intent(out) :: kpoints_out(:,:)

    integer :: ierr, ikp

    allocate(kpoints_out(3,nkpts_all), irr2kpt_out(nkpts_all), kpt2irr_out(nkpts_all), STAT=ierr)
    if(ierr/=0) stop 'problem allocating kpoints, irr2kpt etc.'

    do ikp=1,nkpts_all
      kpoints_out(:,ikp) = kpoints_in(:,kpt2irr_in(ikp))
    end do

    do ikp=1,nkpts_all
      irr2kpt_out(ikp) = ikp
      kpt2irr_out(ikp) = ikp
    end do

  end subroutine unfold_visarrays




  subroutine expand_visarrays(nsym, nkpts1, nkpts_irr1, kpt2irr1, irr2kpt1, kpt2irr, irr2kpt)
    implicit none
    integer, intent(in) :: nsym, nkpts1, nkpts_irr1
    integer, intent(in) :: kpt2irr1(nkpts1), irr2kpt1(nkpts_irr1)
    integer, allocatable, intent(out) :: kpt2irr(:), irr2kpt(:)

    integer :: ierr, isy, ub, lb

    allocate(kpt2irr(nsym*nkpts1), irr2kpt(nsym*nkpts_irr1), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating kpt2irr in expand_visarrays'

    do isy=1,nsym
      lb = (isy-1)*nkpts1+1
      ub = isy*nkpts1
      kpt2irr(lb:ub) = kpt2irr1 + (isy-1)*nkpts_irr1
    end do!isy

    do isy=1,nsym
      lb = (isy-1)*nkpts_irr1+1
      ub = isy*nkpts_irr1
      irr2kpt(lb:ub) = irr2kpt1 + (isy-1)*nkpts1
    end do!isy

  end subroutine expand_visarrays





  subroutine expand_areas(nsym,nkpts_in,areas_in,areas_out)
    implicit none
    integer, intent(in) :: nsym, nkpts_in
    double precision, intent(in) :: areas_in(nkpts_in)
    double precision, allocatable, intent(out) :: areas_out(:)
    integer :: ierr, isy, ub, lb

    allocate(areas_out(nsym*nkpts_in), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating areas_out'

    do isy=1,nsym
      lb = (isy-1)*nkpts_in+1
      ub = isy*nkpts_in
      areas_out(lb:ub) = areas_in
    end do!isy

  end subroutine expand_areas



  subroutine expand_spinvalues(nsym,ndegen,nsqa,nkpts_in,spinval_in,spinval_out)
    implicit none
    integer, intent(in) :: nsym, ndegen, nsqa, nkpts_in
    double precision, intent(in) :: spinval_in(ndegen,nsqa,nkpts_in)
    double precision, allocatable, intent(out) :: spinval_out(:,:,:)
    integer :: ierr, isy, ub, lb

    allocate(spinval_out(ndegen,nsqa,nsym*nkpts_in), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating spinval_out'

    do isy=1,nsym
      lb = (isy-1)*nkpts_in+1
      ub = isy*nkpts_in
      spinval_out(:,:,lb:ub) = spinval_in
    end do!isy

  end subroutine expand_spinvalues




  subroutine expand_torqvalues(nsym,ndegen,nkpts_in,torqval_in,torqval_out)
    implicit none
    integer, intent(in) :: nsym, ndegen, nkpts_in
    double precision, intent(in) :: torqval_in(3,ndegen,nkpts_in)
    double precision, allocatable, intent(out) :: torqval_out(:,:,:)
    integer :: ierr, isy, ub, lb

    allocate(torqval_out(3,ndegen,nsym*nkpts_in), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating torqval_out'

    do isy=1,nsym
      lb = (isy-1)*nkpts_in+1
      ub = isy*nkpts_in
      torqval_out(:,:,lb:ub) = torqval_in
    end do!isy

  end subroutine expand_torqvalues




  subroutine rotate_kpoints(rotmat, nkpts_in, kpoints_in, nsym, isym, nkpts_out, kpoints_out)

    implicit none
    integer,          intent(in)      :: nkpts_in, nsym, isym(nsym)
    double precision, intent(in)      :: kpoints_in(3,nkpts_in), rotmat(64,3,3)
    integer,                       intent(out) :: nkpts_out
    double precision, allocatable, intent(out) :: kpoints_out(:,:)

    integer :: ierr, ipoint, isy, iout

    allocate(kpoints_out(3,nkpts_in*nsym), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating kpoints_out'

    do isy=1,nsym
      do ipoint=1,nkpts_in
        iout = (isy-1)*nkpts_in+ipoint
        kpoints_out(:,iout) =  rotmat(isym(isy),1,:)*kpoints_in(1,ipoint) &
                           & + rotmat(isym(isy),2,:)*kpoints_in(2,ipoint) &
                           & + rotmat(isym(isy),3,:)*kpoints_in(3,ipoint)
      end do!ipoint
    end do!isy

    nkpts_out = nkpts_in*nsym

  end subroutine rotate_kpoints



  function points_in_wedge(nfaces, nvec, dscal, npts, points, mode)

    logical :: points_in_wedge
    integer, intent(in) :: nfaces, npts
    double precision, intent(in) :: nvec(3,nfaces), dscal(nfaces), points(3,npts)
    character(len=*), intent(in) :: mode

    integer :: ipt
    logical :: lpoints(npts)

    do ipt=1,npts
      lpoints(ipt) = singlepoint_in_wedge(nfaces, nvec, dscal, points(:,ipt))
    end do`!ipt

    select case( mode )
    case( 'any' ) ; points_in_wedge = any(lpoints)
    case( 'all' ) ; points_in_wedge = all(lpoints)
    case default  ; stop "Mode for 'points_in_wedge' not known."
    end select!icase


  end function points_in_wedge


  function singlepoint_in_wedge(nfaces, nvec, dscal, point)

    logical :: singlepoint_in_wedge
    integer, intent(in) :: nfaces
    double precision, intent(in) :: nvec(3,nfaces), dscal(nfaces), point(3)

    integer :: iface
    double precision :: dtmp

    double precision, parameter :: eps=1d-16

    singlepoint_in_wedge=.true.
    do iface=1,nfaces
      dtmp = sum(nvec(:,iface)*point)-dscal(iface)+eps
      if(dtmp<0d0)then
        singlepoint_in_wedge=.false.
        return
      end if!eps
    end do!iface

  end function singlepoint_in_wedge



  subroutine get_2DIBZwedge_lines(recbv, nsym, rotmat, isym, nBZlines, nvec, dscal, bounds )
    ! Generates the equations for the faces of the IBZ-wedge (2D-Version),
    !  such that for all points in the IBZ the inequality
    !      nvec(:,i) * k(:) >= dvec(i)
    !  holds.
    !
    ! NOTE: nvec, dvec etc. are 3D-vectors due to compatibility with the bulk-version of the program
    !
    ! In the end, the IBZ and its equvalents by symmetry
    !  are plotted as matplotlib-file.
    !
    !                                 B.Zimmermann, Mai 2014
    use mod_ioinput,   only: IoInput
    use mod_mathtools, only: crossprod
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer,                       intent(in)  :: nsym, isym(nsym)
    double precision,              intent(in)  :: rotmat(64,3,3), recbv(3,3)
    integer,                       intent(out) :: nBZlines
    double precision, allocatable, intent(out) :: nvec(:,:), dscal(:)
    double precision,              intent(out) :: bounds(3,2)

    integer                       :: npoints!, isymtmp(1)
    double precision              :: dtest, diff(3), vecinp(3), testpoint(3), diffvec(3,2)
    integer, allocatable          :: verts2d(:)
    double precision, allocatable :: points(:,:)
    logical :: cartesian=.false.

    integer :: ipoint, iline, ierr,isy
    character(len=80)    :: uio

    !=====
    ! Read in high symmetry points and number of BZ face lines
    !=====
    if(myrank==master)then
      call IoInput('NSYMPTS   ',uio,1,7,ierr)
      read (unit=uio,fmt=*) npoints

      call IoInput('CARTESIAN ',uio,1,7,ierr)
      read (unit=uio,fmt=*) cartesian
!     write(*,*) 'cartesian=', cartesian

      allocate(points(3,npoints), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating high symmetry points etc.'

      do ipoint=1,npoints
        call IoInput('SYMCOORD  ',uio,ipoint,7,ierr)
        read (unit=uio,fmt=*) vecinp(:)
        if(cartesian)then
          points(:,ipoint) = vecinp(:)
        else!cartesian
          points(:,ipoint) = recbv(:,1)*vecinp(1) + recbv(:,2)*vecinp(2) + recbv(:,3)*vecinp(3)
        end if!cartesian

        if(abs(points(3,ipoint))>1d-16) stop 'z-coordinate of SYMCOORD not 0 for 2D system'
      end do!ipoint


      call IoInput('NIBZ2D    ',uio,1,7,ierr)
      read (unit=uio,fmt=*) nBZlines

      allocate(verts2D(nBZlines+1), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating nfaceverts'

      call IoInput('VERTS2D   ',uio,1,7,ierr)
      read (unit=uio,fmt=*) verts2D

      call IoInput('POINTIBZ  ',uio,1,7,ierr)
      read (unit=uio,fmt=*) vecinp(:)
      if(cartesian)then
        testpoint(:) = vecinp(:)
      else!cartesian
        testpoint(:) = recbv(:,1)*vecinp(1) + recbv(:,2)*vecinp(2) + recbv(:,3)*vecinp(3)
      end if!cartesian
!     write(*,*) testpoint

      bounds(:,1) = minval(points, dim=2)
      bounds(:,2) = maxval(points, dim=2)

    end if!myrank==master


    !=====
    ! Allocate result arrays
    !=====
#ifdef CPP_MPI
    call MPI_Bcast(nBZlines,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nfaces'
#endif
    allocate(nvec(3,nBZlines), dscal(nBZlines), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating nvec etc.'


    !=====
    ! Compute the equations for the BZ lines
    !=====
    if(myrank==master)then
      do iline=1,nBZlines
        !get difference vectors between three points of a face
        diffvec(:,1) = points(:,verts2D(iline+1)) - points(:,verts2D(iline))
        diffvec(:,2) = (/ 0d0, 0d0, 1d0 /)

        !calculate the normal of the face
        call crossprod(diffvec(:,1),diffvec(:,2),nvec(:,iline))

        !calculate the distance of plane to origin
        dscal(iline) = sum( nvec(:,iline)*points(:,verts2D(iline)) )

        !calculate the distance of the testpoint
        dtest = sum( nvec(:,iline)*testpoint )

        if(dtest<dscal(iline))then
          nvec(:,iline) = -nvec(:,iline)
          dscal(iline)  = -dscal(iline)
        end if

      end do!iface
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(nvec,3*nBZlines,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nvec'
    call MPI_Bcast(dscal,nBZlines,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nvec'
    call MPI_Bcast(bounds,6,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nvec'
#endif


    !=====
    ! Write info to the screen
    !=====
    if(myrank==master)then

      write(*,*)
      write(*,'("Reciprocal lattice:")')
      write(*,'(3ES25.16)') recbv
      write(*,*)

      write(*,'("Coordinates of the corner points of the IBZ:")')
      do ipoint=1,npoints
        write(*,'(2X,"Point ",I0,": (",3F10.6,")")') ipoint, points(:,ipoint)
      end do!ipoint
      write(*,*)

      write(*,'("Equations for the lines of the IBZ:")')
      do iline=1,nBZlines
        write(*,'(2X,"n=(",3F10.6,"), d=",F10.6)') nvec(:,iline), dscal(iline)
      end do!iface
      write(*,*)

      write(*,'("Bounds for the IBZ:")')
      write(*,'(" min:",3F12.6)') bounds(:,1)
      write(*,'(" max:",3F12.6)') bounds(:,2)

      open(unit=1359,file='bounds',form='formatted',action='write')
      write(1359,'(3ES25.16)') bounds
      close(1359)

      diff = bounds(:,2)-bounds(:,1)
      diff = diff/minval(diff(1:2))
      write(*,'(" fac:",3F12.1)') diff
      write(*,*)

    end if!myrank==master

    !=====
    ! Create TXT-File containing the IBZ
    !=====
    if(myrank==master)then
      open(unit=19535, file='ibz_wedge.txt', form='formatted', action='write')
      write(19535,'(A)') '# 2D IBZ. Each line represents one IBZ point'
      do ipoint=1,nBZlines+1
        write(19535,'(3ES25.16)') points(:,verts2D(ipoint))
      end do
!      write(19535,'(A)') '# 2D IBZ. Each line represents one IBZ straigt line, containing the start and end point'
!      do iline=1,nBZlines
!        write(19535,'(6ES25.16)') points(:,verts2D(iline)), points(:,verts2D(iline+1))
!      end do
      close(19535)


      open(unit=19536, file='ibz_fullbz.txt', form='formatted', action='write')
      write(19536,'(A)') '# 2D IBZ. Each line represents one BZ point, first number is number of symmetry operation'
      write(19536,'(2I8)') nsym, nBZlines+1
      do isy=1,nsym
        do ipoint=1,nBZlines+1
          write(19536,'(I8,3ES25.16)') isy, rotmat(isym(isy),1,:)*points(1,verts2D(ipoint)) &
                                     &    + rotmat(isym(isy),2,:)*points(2,verts2D(ipoint)) &
                                     &    + rotmat(isym(isy),3,:)*points(3,verts2D(ipoint))
        end do!ipoint
      end do!isy
      close(19536)

    end if!myrank==master


  end subroutine get_2DIBZwedge_lines





  subroutine get_IBZwedge_faces(recbv, nsym, rotmat, isym, nfaces, nvec, dscal, bounds)
    ! Generates the equations for the faces of the IBZ-wedge,
    !  such that for all points in the IBZ the inequality
    !      nvec(:,i) * k(:) >= dvec(i)
    !  holds.
    ! In the end, the IBZ and its equvalents by symmetry
    !  are plotted as vtk-file.
    !
    !                                 B.Zimmermann, Mar.2013

    use mod_ioinput,   only: IoInput
    use mod_mathtools, only: crossprod
    use mod_vtkxml,    only: write_IBZ_rot
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    integer,                       intent(in)  :: nsym, isym(nsym)
    double precision,              intent(in)  :: rotmat(64,3,3), recbv(3,3)
    integer,                       intent(out) :: nfaces
    double precision, allocatable, intent(out) :: nvec(:,:), dscal(:)
    double precision,              intent(out) :: bounds(3,2)

    integer                       :: npoints, isymtmp(1)
    double precision              :: dtest, diff(3), vecinp(3), testpoint(3), diffvec(3,2)
    integer, allocatable          :: nfaceverts(:), ifaceverts(:,:)
    double precision, allocatable :: points(:,:)
    logical :: cartesian=.false.

    integer :: ipoint, iface, ierr
    character(len=80)    :: uio



    !=====
    ! Read in high symmetry points and number of faces
    !=====
    if(myrank==master)then
      call IoInput('NSYMPTS   ',uio,1,7,ierr)
      read (unit=uio,fmt=*) npoints

      call IoInput('CARTESIAN ',uio,1,7,ierr)
      read (unit=uio,fmt=*) cartesian
!     write(*,*) 'cartesian=', cartesian

      allocate(points(3,npoints), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating high symmetry points etc.'

      do ipoint=1,npoints
        call IoInput('SYMCOORD  ',uio,ipoint,7,ierr)
        read (unit=uio,fmt=*) vecinp(:)
        if(cartesian)then
          points(:,ipoint) = vecinp(:)
        else!cartesian
          points(:,ipoint) = recbv(:,1)*vecinp(1) + recbv(:,2)*vecinp(2) + recbv(:,3)*vecinp(3)
        end if!cartesian
!       write(*,*) points(:,ipoint)
      end do!ipoint

      call IoInput('NIBZFAC   ',uio,1,7,ierr)
      read (unit=uio,fmt=*) nfaces

      allocate(nfaceverts(nfaces), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating nfaceverts'

      do iface=1,nfaces
        call IoInput('NVERT     ',uio,iface,7,ierr)
        read (unit=uio,fmt=*) nfaceverts(iface)
!       write(*,*) nfaceverts(iface)
      end do!iface

      allocate(ifaceverts(maxval(nfaceverts),nfaces), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating ifaceverts'
      ifaceverts = 0

      do iface=1,nfaces
        call IoInput('VERTIDS   ',uio,iface,7,ierr)
        read (unit=uio,fmt=*) ifaceverts(1:nfaceverts(iface),iface)
!       write(*,'(10I3)') ifaceverts(:,iface)
      end do!iface

      call IoInput('POINTIBZ  ',uio,1,7,ierr)
      read (unit=uio,fmt=*) vecinp(:)
      if(cartesian)then
        testpoint(:) = vecinp(:)
      else!cartesian
        testpoint(:) = recbv(:,1)*vecinp(1) + recbv(:,2)*vecinp(2) + recbv(:,3)*vecinp(3)
      end if!cartesian
!     write(*,*) testpoint

      bounds(:,1) = minval(points, dim=2)
      bounds(:,2) = maxval(points, dim=2)

    end if!myrank==master



    !=====
    ! Allocate result arrays
    !=====
#ifdef CPP_MPI
    call MPI_Bcast(nfaces,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nfaces'
#endif
    allocate(nvec(3,nfaces), dscal(nfaces), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating nvec etc.'



    !=====
    ! Compute the equations for the planes of the BZ faces
    !=====
    if(myrank==master)then
      do iface=1,nfaces
        !get difference vectors between three points of a face
        diffvec(:,1) = points(:,ifaceverts(2,iface)) - points(:,ifaceverts(1,iface))
        diffvec(:,2) = points(:,ifaceverts(3,iface)) - points(:,ifaceverts(1,iface))

        !calculate the normal of the face
        call crossprod(diffvec(:,1),diffvec(:,2),nvec(:,iface))

        !calculate the distance of plane to origin
        dscal(iface) = sum( nvec(:,iface)*points(:,ifaceverts(1,iface)) )

        !calculate the distance of the testpoint
        dtest = sum( nvec(:,iface)*testpoint )

        if(dtest<dscal(iface))then
          nvec(:,iface) = -nvec(:,iface)
          dscal(iface)  = -dscal(iface)
        end if

      end do!iface
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(nvec,3*nfaces,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nvec'
    call MPI_Bcast(dscal,nfaces,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nvec'
    call MPI_Bcast(bounds,6,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nvec'
#endif



    !=====
    ! Write info to the screen
    !=====
    if(myrank==master)then

      write(*,'("Coordinates of the corner points of the IBZ:")')
      do ipoint=1,npoints
        write(*,'(2X,"Point ",I0,": (",3F10.6,")")') ipoint, points(:,ipoint)
      end do!ipoint
      write(*,*)

      write(*,'("Equations for the faces of the IBZ:")')
      do iface=1,nfaces
        write(*,'(2X,"n=(",3F10.6,"), d=",F10.6)') nvec(:,iface), dscal(iface)
      end do!iface
      write(*,*)

      write(*,'("Bounds for the IBZ:")')
      write(*,'(" min:",3F12.6)') bounds(:,1)
      write(*,'(" max:",3F12.6)') bounds(:,2)

      open(unit=1359,file='bounds',form='formatted',action='write')
      write(1359,'(3ES25.16)') bounds
      close(1359)

      diff = bounds(:,2)-bounds(:,1)
      diff = diff/minval(diff)
      write(*,'(" fac:",3F12.1)') diff
      write(*,*)

    end if!myrank==master


    !=====
    ! Create VTK-File containing the IBZ
    !=====
    if(myrank==master) call write_IBZ_rot('ibz_fullbz.vtp',npoints,points,nfaces,nfaceverts,ifaceverts,nsym,rotmat,isym)
    isymtmp = 1
    if(myrank==master) call write_IBZ_rot('ibz_wedge.vtp',npoints,points,nfaces,nfaceverts,ifaceverts,1,rotmat,isymtmp)

  end subroutine get_IBZwedge_faces



  subroutine set_symmetries(inc, lattice, symmetries)
  ! This subroutine sets the symmetry infortmation
  !  and stores them in the derived datatype symmetries.
  ! The symmetry information as found from 'findgroup',
  !  as well as the user-defined symmetries from the
  !  input-file are concerned.

    use type_inc,  only: inc_type
    use type_data, only: lattice_type
#ifdef CPP_MPI
    use mpi
#endif
    implicit none

    type(inc_type),     intent(in)     :: inc
    type(lattice_type), intent(in)     :: lattice
    type(symmetries_type), intent(out) :: symmetries

    integer :: ierr, isy1, isy2, icount

    !=====
    ! Find symmetry information from the lattice
    !=====
    call findgroup( inc%naez,inc%naezd,inc%nembd,lattice%bravais,         &
                  & lattice%rbasis,lattice%alat,lattice%recbv,inc%nBZdim, &
                  & symmetries%rotmat,symmetries%rotname,                 &
                  & symmetries%nsym_found,symmetries%isym_found           )
    !=====
    ! Find symmetry information from the input-file
    !=====
    if(myrank==master) then
      call read_sym_inp( symmetries%rotname,symmetries%nsym_used,symmetries%isym_used )
    end if!myrank==master

#ifdef CPP_MPI
    call MPI_Bcast(symmetries%nsym_used,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting nsym'
    if(myrank/=master)then
      allocate( symmetries%isym_used(symmetries%nsym_used), STAT=ierr )
      if(ierr/=0) stop 'Problem allocating symmetries%isym_used on slaves'
    end if
    call MPI_Bcast(symmetries%isym_used,symmetries%nsym_used,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem broadcasting symmetries%isym_used'
#endif

    !=====
    ! Find symmetries that are exluded (i.e. valid ones for the lattice but not specified by the user)
    !=====
    allocate(symmetries%isym_notused(symmetries%nsym_found - symmetries%nsym_used), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating symmetries%isym_notused'
    icount=0
    outer: do isy1=1,symmetries%nsym_found
      inner: do isy2=1,symmetries%nsym_used
        if(symmetries%isym_found(isy1) == symmetries%isym_used(isy2)) cycle outer
      end do inner!isy2
      icount = icount+1
      symmetries%isym_notused(icount) = symmetries%isym_found(isy1)
    end do outer!isy1

    if(myrank==master)then
      write(6,*)
      write(6,'(" Symmetries NOT used: ",I0)') icount
      write(6,1050) (symmetries%rotname(symmetries%isym_notused(isy1)), isy1=1,icount )
      write(6,*)
    end if!myrank==master

  1050 FORMAT(5(A10,2X))
  end subroutine set_symmetries



  subroutine read_sym_inp(rotname, nsym, isymindex)
  !This subroutine reads symmetry information from the inputfile.
  !  In the inputfile the names of symmetry transformations
  !  according to the subroutine 'pointgrp' must be given.
  !  INPUT: the array 'rotname' containing the names of
  !         all possible symmetry operations.
  !  OUTPUT: nsym contains the number of symmetry operations used
  !          isymindex contains the id's of the used symmetries

    use mod_ioinput, only: IoInput
    implicit none

    character(len=10),    intent(in)  :: rotname(64)
    integer,              intent(out) :: nsym
    integer, allocatable, intent(out) :: isymindex(:)

    integer              :: irow, icol, isym, istore, rest, ierr
    integer, allocatable :: linespercol(:)
    character(len=10)    :: colname
    character(len=80)    :: uio
    character(len=10)    :: charstr(64)

    integer, parameter :: ncols = 5

    !read in number of symmetries to use
    call IoInput('NSYM      ',uio,1,7,ierr)
    read (unit=uio,fmt=*) nsym

    !create storage array
    allocate( isymindex(nsym), STAT=ierr )
    if(ierr/=0) stop 'Problem allocating isymindex'

    !create temporary variables to read in symmetries
    allocate(linespercol(ncols))
    linespercol = nsym/ncols
    rest = nsym - linespercol(1)*ncols
    linespercol(1:rest) = linespercol(1:rest)+1

    do icol=1,ncols
      write(colname,'("SYMFIELD",I0)') icol
      do irow = 1,linespercol(icol)
        call IoInput(colname,uio,1+irow,7,ierr)
        symloop: do isym=1,64
          if(trim(uio(1:10))==trim(rotname(isym)))then
            istore = (irow-1)*ncols+icol
            isymindex(istore) = isym
            exit symloop
          end if
        end do symloop!isym
      end do!irow
    end do!icol

    !write info to the screen
    write(6,1040) nsym
    do isym=1,nsym
      istore        = isymindex(isym)
      charstr(isym) = rotname(istore)
    end do
    write(6,1030) (charstr(isym),isym=1,nsym)
    write(6,*)

  1030 FORMAT(5(A10,2X))
  1040 FORMAT(' Symmetries set by hand: ',I5)
  end subroutine read_sym_inp

end module mod_symmetries
