!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


module mod_symmetries
! This module provides routines regardinf the symmetry of the lattice and the irreducible brillouin zone.


  use mod_mympi, only: nranks, myrank, master
  implicit none

  private
  public :: symmetries_type, set_symmetries, pointgrp, get_IBZwedge_faces, points_in_wedge, singlepoint_in_wedge, rotate_kpoints, expand_visarrays, expand_areas, expand_spinvalues, expand_torqvalues, unfold_visarrays, get_2DIBZwedge_lines

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
    end do!ipt

    selectcase( mode )
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



  !-------------------------------------------------------------------------------
  !> Summary: Find symmetry operations that leave crystal lattice invariant
  !> Author: 
  !> Category: PKKprime, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> @note copied from host code @endnote
  !>
  !> This subroutine finds the rotation matrices that leave the
  !> real lattice unchanged. 
  !> input:  bravais(i,j)    true bravais lattice vectors
  !>                         i = x,y,z ; j = A, B, C (a.u.)
  !>         recbv(i,j)      reciprocal basis vectors
  !>         rbasis          coordinates of basis atoms
  !>         nbasis          number of basis atoms
  !>         rfctor          alat/4/pi
  !> output: rotmat          all 64 rotation matrices.
  !>         rotname         names for the rotation matrices
  !>         nsymat          number of rotations that restore the lattice.
  !>         isymindex       index for the symmeties found
  !>
  !> This sub makes all 64 rotations in the basis vectors and bravais
  !> vectors and checks if the new rotated vectror belongs in the 
  !> lattice. The proper rotation must bring all vectors to a lattice
  !> vector. Information about the rotations found is printed in the end.
  !> The array isymindex holds the numbers of the symmetry operations
  !> that are stored in array RSYMAT
  !-------------------------------------------------------------------------------
  subroutine findgroup( nbasis,naezd,nembd,bravais,rbasis,  &
                      & rfctor,recbv,nbzdim,                &
                      & rotmat,rotname,nsymat,isymindex_out )

      implicit none

      integer,          intent(in) :: nbasis, naezd, nembd, nbzdim
      double precision, intent(in) :: bravais(3,3),rbasis(3,naezd+nembd)
      double precision, intent(in) :: rfctor,recbv(3,3)

      double precision,     intent(out) :: rotmat(64,3,3)
      character(len=10),    intent(out) :: rotname(64)
      integer,              intent(out) :: nsymat
      integer, allocatable, intent(out) :: isymindex_out(:)

      integer, parameter :: NSYMAXD=48
      integer :: isymindex(NSYMAXD)

      ! Local variables
      double precision :: r(3,4),rotrbas(3,naezd+nembd)
      double precision :: alat,bravais1(3,3)
      integer :: i,j,isym,nsym,i0,ia
      character(len=10) :: charstr(64)
      logical :: llatbas,lbulk
!     -------------------------------------------------------------
      nsym = 0
      call pointgrp(rotmat,rotname)
!     alat = rfctor*8.d0*datan(1.0d0)
!     - ---------------------------------
      do i=1,3
         do j=1,3
            bravais1(j,i) = bravais(j,i) !/alat
         end do
      end do
      !Check for surface mode. If so, set bravais1(3,3) very large, so
      !that only the in-plane symmetries are found. Not checked, be careful of z--> -z!
      if(nbzdim==2)then
        lbulk=.false.
      else!nbzdim==2
        lbulk=.true.
      end if!nbzdim==2
!     !Now check the bravais vectors if they have a z component
!     if ((bravais(1,3).eq.0.d0).and.(bravais(2,3).eq.0.d0).and.&
!        & (bravais(3,3).eq.0.d0)) then
!        lbulk=.false.
!     end if

!     write(100,*) 'bravais:'
!     write(100,'(3ES25.16)') bravais
!     write(100,*) 'lbulk=', lbulk

      do isym=1,64
         !rotate bravais lattice vectors

         !In the case of slab/interface geometry look only for
         !symmetry opperations that preserve the z axis..
         if (lbulk .or. (rotmat(isym,3,3).eq.1) ) then
            !do rotation only in case bulk or if slab and z axis is restored..
            do i=1,3            ! Loop on bravais vectors
               do j=1,3         ! Loop on coordinates
                  r(j,i) = rotmat(isym,j,1)*bravais1(1,i) + &
                         & rotmat(isym,j,2)*bravais1(2,i) + &
                         & rotmat(isym,j,3)*bravais1(3,i)
               enddo
            enddo

            !rotate the basis atoms p and take RSYMAT.p - p then
            !find if R = (RSYMAT.bravais + RSYMAT.p - p) belongs to the
            !lattice. This is done by function latvec by checking
            !if R.q = integer (q reciprocal lattice vector)

            llatbas = .true.
            do ia=1,nbasis      ! Loop on basis atoms
               do j=1,3         ! Loop on coordinates
                  rotrbas(j,ia) = rotmat(isym,j,1)*rbasis(1,ia) + &
                                & rotmat(isym,j,2)*rbasis(2,ia) + &
                                & rotmat(isym,j,3)*rbasis(3,ia)

                  rotrbas(j,ia) = rotrbas(j,ia) - rbasis(j,ia)
                  r(j,4) = rotrbas(j,ia) 
               enddo
               if (.not.latvec(4,recbv,r)) llatbas=.false.
            enddo               ! ia=1,nbasis

            !if llatbas=.true. the rotation does not change the lattice 
            if (llatbas) then
               nsym = nsym + 1
               isymindex(nsym) = isym
            end if
         end if                 ! (lbulk .OR. (rotmat(isym,3,3).EQ.1) )
      end do                    ! isym=1,nmatd

      !nsym symmetries were found
      !the isymindex array has the numbers of the symmetries found
      nsymat = nsym
      allocate(isymindex_out(nsymat))
      isymindex_out(:) = isymindex(1:nsymat)

      !write info to the screen
      if(myrank==master)then
        write(6,*) 'Information from FindGroup'
        if (.not.lbulk) write(6,*) 'Surface Symmetries '
        write(6,1020) nsymat
        do i=1,nsymat
           I0 = isymindex(i)
           charstr(i) =  rotname(I0)
        end do
        write(6,1010) (charstr(i),i=1,nsymat)
        write(6,*) '----------- * findgroup ends here * ---------------'
        write(6,*)
      end if!myrank==master
 1010 FORMAT(5(A10,2X))
 1020 FORMAT(' Symmetries found for this lattice: ',I5)

  end subroutine findgroup



  !-------------------------------------------------------------------------------
  !> Summary: Rotation matrices of 32 point groups
  !> Author: 
  !> Category: PKKprime, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> @note copied from host code @endnote
  !> This subroutine defines the rotation matrices for
  !> all the 32 point groups and names them after 
  !> J.F. Cornwell (Group Theory?) second edition 
  !> Appendix D, p 324-325
  !-------------------------------------------------------------------------------
  subroutine pointgrp(rotmat,rotname)

      implicit none

      double precision,  intent(out) :: ROTMAT(64,3,3)
      character(len=10), intent(out) :: ROTNAME(64)

      !Locals
      integer i,j,i1,is
      double precision RTHREE,HALF 
      integer, parameter :: iou=3514
      logical :: writesymfile = .true.

      RTHREE = sqrt(3.d0)/2.d0
      HALF = 0.5d0
! set to zero
      do i1=1,64  
          do i=1,3
             do j=1,3
                ROTMAT(i1,i,j) = 0.d0
             end do
          end do
      end do

      ROTMAT(1,1,1) =  1.d0
      ROTMAT(1,2,2) =  1.d0
      ROTMAT(1,3,3) =  1.d0
      ROTNAME(1) = 'E'

      ROTMAT(2,1,2) =  1.d0
      ROTMAT(2,2,3) = -1.d0
      ROTMAT(2,3,1) = -1.d0
      ROTNAME(2) = 'C3alfa'          

      ROTMAT(3,1,2) = -1.d0
      ROTMAT(3,2,3) = -1.d0
      ROTMAT(3,3,1) =  1.d0
      ROTNAME(3) = 'C3beta '

      ROTMAT(4,1,2) = -1.d0
      ROTMAT(4,2,3) =  1.d0
      ROTMAT(4,3,1) = -1.d0
      ROTNAME(4) = 'C3gamma'

      ROTMAT(5,1,2) = 1.d0
      ROTMAT(5,2,3) = 1.d0
      ROTMAT(5,3,1) = 1.d0
      ROTNAME(5) = 'C3delta '

      ROTMAT(6,1,3) = -1.d0
      ROTMAT(6,2,1) =  1.d0
      ROTMAT(6,3,2) = -1.d0
      ROTNAME(6) = 'C3alfa-1'

      ROTMAT(7,1,3) =  1.d0
      ROTMAT(7,2,1) = -1.d0
      ROTMAT(7,3,2) = -1.d0
      ROTNAME(7) = 'C3beta-1 '

      ROTMAT(8,1,3) = -1.d0
      ROTMAT(8,2,1) = -1.d0
      ROTMAT(8,3,2) =  1.d0
      ROTNAME(8) = 'C3gamma-1'

      ROTMAT(9,1,3) =  1.d0
      ROTMAT(9,2,1) =  1.d0
      ROTMAT(9,3,2) =  1.d0
      ROTNAME(9) = 'C3delta-1'

      ROTMAT(10,1,1) =  1.d0
      ROTMAT(10,2,2) = -1.d0
      ROTMAT(10,3,3) = -1.d0
      ROTNAME(10) = 'C2x'
           
      ROTMAT(11,1,1) = -1.d0
      ROTMAT(11,2,2) =  1.d0
      ROTMAT(11,3,3) = -1.d0
      ROTNAME(11) = 'C2y'

      ROTMAT(12,1,1) = -1.d0
      ROTMAT(12,2,2) = -1.d0
      ROTMAT(12,3,3) =  1.d0
      ROTNAME(12) = 'C2z'

      ROTMAT(13,1,1) =  1.d0
      ROTMAT(13,2,3) =  1.d0
      ROTMAT(13,3,2) = -1.d0
      ROTNAME(13) = 'C4x'
           
      ROTMAT(14,1,3) = -1.d0
      ROTMAT(14,2,2) =  1.d0
      ROTMAT(14,3,1) =  1.d0
      ROTNAME(14) = 'C4y '

      ROTMAT(15,1,2) =  1.d0
      ROTMAT(15,2,1) = -1.d0
      ROTMAT(15,3,3) =  1.d0
      ROTNAME(15) = 'C4z'
           
      ROTMAT(16,1,1) =  1.d0
      ROTMAT(16,2,3) = -1.d0
      ROTMAT(16,3,2) =  1.d0
      ROTNAME(16) = 'C4x-1 '

      ROTMAT(17,1,3) =  1.d0
      ROTMAT(17,2,2) =  1.d0
      ROTMAT(17,3,1) = -1.d0
      ROTNAME(17) = 'C4y-1'

      ROTMAT(18,1,2) = -1.d0
      ROTMAT(18,2,1) =  1.d0
      ROTMAT(18,3,3) =  1.d0
      ROTNAME(18) = 'C4z-1'
           
      ROTMAT(19,1,2) =  1.d0
      ROTMAT(19,2,1) =  1.d0
      ROTMAT(19,3,3) = -1.d0
      ROTNAME(19) = 'C2a'

      ROTMAT(20,1,2) = -1.d0
      ROTMAT(20,2,1) = -1.d0
      ROTMAT(20,3,3) = -1.d0
      ROTNAME(20) = 'C2b'

      ROTMAT(21,1,3) =  1.d0
      ROTMAT(21,2,2) = -1.d0
      ROTMAT(21,3,1) =  1.d0
      ROTNAME(21) = 'C2c'

      ROTMAT(22,1,3) = -1.d0
      ROTMAT(22,2,2) = -1.d0
      ROTMAT(22,3,1) = -1.d0
      ROTNAME(22) = 'C2d'

      ROTMAT(23,1,1) = -1.d0
      ROTMAT(23,2,3) =  1.d0
      ROTMAT(23,3,2) =  1.d0
      ROTNAME(23) = 'C2e'

      ROTMAT(24,1,1) = -1.d0
      ROTMAT(24,2,3) = -1.d0
      ROTMAT(24,3,2) = -1.d0
      ROTNAME(24) = 'C2f'
      do i1=1,24
         do i=1,3
            do j=1,3 
              ROTMAT(i1+24,i,j) = -ROTMAT(i1,i,j)
            end do
         end do
      ROTNAME(i1+24) = 'I'//ROTNAME(i1)
      end do

!*********************************************
! Trigonal and hexagonal groups
!*********************************************

      ROTMAT(49,1,1) = -HALF
      ROTMAT(49,1,2) =  RTHREE
      ROTMAT(49,2,1) = -RTHREE
      ROTMAT(49,2,2) = -HALF
      ROTMAT(49,3,3) =  1.d0
      ROTNAME(49) = 'C3z'  

      ROTMAT(50,1,1) = -HALF
      ROTMAT(50,1,2) = -RTHREE
      ROTMAT(50,2,1) =  RTHREE
      ROTMAT(50,2,2) = -HALF
      ROTMAT(50,3,3) =  1.d0
      ROTNAME(50) = 'C3z-1'

      ROTMAT(51,1,1) =  HALF
      ROTMAT(51,1,2) =  RTHREE
      ROTMAT(51,2,1) = -RTHREE
      ROTMAT(51,2,2) =  HALF
      ROTMAT(51,3,3) =  1.d0
      ROTNAME(51) = 'C6z'

      ROTMAT(52,1,1) =  HALF
      ROTMAT(52,1,2) = -RTHREE
      ROTMAT(52,2,1) =  RTHREE
      ROTMAT(52,2,2) =  HALF
      ROTMAT(52,3,3) =  1.d0
      ROTNAME(52) = 'C6z-1'

      ROTMAT(53,1,1) = -HALF
      ROTMAT(53,1,2) =  RTHREE
      ROTMAT(53,2,1) =  RTHREE
      ROTMAT(53,2,2) =  HALF
      ROTMAT(53,3,3) = -1.d0
      ROTNAME(53) = 'C2A'    

      ROTMAT(54,1,1) = -HALF
      ROTMAT(54,1,2) = -RTHREE
      ROTMAT(54,2,1) = -RTHREE
      ROTMAT(54,2,2) =  HALF
      ROTMAT(54,3,3) = -1.d0
      ROTNAME(54) = 'C2B'

      ROTMAT(55,1,1) =  HALF
      ROTMAT(55,1,2) = -RTHREE
      ROTMAT(55,2,1) = -RTHREE
      ROTMAT(55,2,2) = -HALF
      ROTMAT(55,3,3) = -1.d0
      ROTNAME(55) = 'C2C'

      ROTMAT(56,1,1) =  HALF
      ROTMAT(56,1,2) =  RTHREE
      ROTMAT(56,2,1) =  RTHREE
      ROTMAT(56,2,2) = -HALF
      ROTMAT(56,3,3) = -1.d0
      ROTNAME(56) = 'C2D'
      do is=1,8
          do i=1,3
             do j=1,3     
                ROTMAT(56+is,i,j) = -ROTMAT(48+is,i,j)
             end do
          end do
          ROTNAME(56+is) = 'I'//ROTNAME(48+is) 
      end do
      

      if(myrank==master .and. writesymfile)then
        open(unit=iou,file='sym.out',form='formatted',action='write')
        write(iou,'(I0)') 64
        do is=1,64
          write(iou,'(A)') ROTNAME(is)
          write(iou,'(3ES25.16)') ROTMAT(is,:,:)
        end do
        close(iou)

        writesymfile = .false.
      end if!myrank==master .and. writesymfile
     

  end subroutine pointgrp



  !-------------------------------------------------------------------------------
  !> Summary: Checks if a set of vectors are lattice vectors
  !> Author: 
  !> Category: PKKprime, geometry
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> @note copied from host code @endnote
  !> Inputs:                                                              
  !>   n     :number of vectors                                           
  !>   qlat  :primitive translation vectors in reciprocal space           
  !>   vec   :double-precision vector                                     
  !> Outputs:                                                             
  !>   latvec:.true. if all vectors are lattice vectors
  !-------------------------------------------------------------------------------
  logical function latvec(n,qlat,vec) 
  
      implicit none 
      ! Passed parameters:                                                    
      integer,          intent(in) :: n 
      double precision, intent(in) :: qlat(3,3),vec(3,n) 
      ! Local parameters:                                                     
      integer :: i,m 
      double precision :: vdiff 
      double precision, parameter :: tol=1.d-6
      ! Intrinsic functions:                                                  
      intrinsic  dabs,dnint

      latvec=.false. 
      do i=1,n
       do m=1,3
         vdiff=vec(1,i)*qlat(1,m)+vec(2,i)*qlat(2,m)+vec(3,i)*qlat(3,m) 
         vdiff=dabs(vdiff-dnint(vdiff)) 
         if (vdiff.gt.tol) return 
       enddo 
      enddo 
      latvec=.true. 

  end function latvec

end module mod_symmetries
