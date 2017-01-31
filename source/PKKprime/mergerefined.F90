!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


program mergerefined

  use mod_read,       only: read_kpointsfile_vis, read_kpointsfile_int
  use mod_ioformat,   only: filename_cubesinfo, filename_fsdata, filename_vtktest,&
                          & fmt_fn_ext, fmt_fn_sub_ext,                           &
                          & filemode_ref, filemode_vis, filemode_int,             &
                          & ext_formatted, ext_vtkxml, ext_new, ext_refined, ext_orig
  use mod_mympi,      only: mympi_init, myrank, nranks, master
  use mod_symmetries, only: pointgrp
  use mod_fermisurf_basic,  only: find_kpoints_irredset, save_kpointsfile_vis, save_kpointsfile_int
  use mod_vtkxml,     only: write_pointdata_rot
#ifdef CPP_MPI
  use mpi
#endif

  implicit none

  integer :: nBZdim, nCub3(3)
  double precision :: bounds(3,2), areatot

  integer :: nkpts_refined, nkpts_irr_refined, nsym_refined, &
           & nkpts_orig,    nkpts_irr_orig,    nsym_orig,    &
           & nkpts_new,     nkpts_irr_new,     nkpts_orig_filtered

  integer, allocatable :: isym_refined(:), kpt2irr_refined(:), irr2kpt_refined(:), cubeids_refined(:), &
                        & isym_orig(:),    kpt2irr_orig(:),    irr2kpt_orig(:),    cubeids_orig(:),    &
                        & isym_new(:),     kpt2irr_new(:),     irr2kpt_new(:),     cubeids_new(:),     &
                        & cubeids_orig_save(:)

  double precision, allocatable :: kpoints_irr_refined(:,:), kpoints_refined(:,:), areas_refined(:), &
                                 & kpoints_irr_orig(:,:),    kpoints_orig(:,:),    areas_orig(:),    &
                                 & kpoints_irr_new(:,:),     kpoints_new(:,:),     areas_new(:)


  integer :: ierr, ikp, itri, itmparr(1)
  character(len=256) :: filename
  integer, parameter :: ifile=1656

  double precision   :: rotmat(64,3,3)
  character(len=10)  :: rotname(64)
  character(len=256) :: scalarstring(3), vectorstring(1)
  double precision, allocatable :: vectordata(:,:,:), &
                                 & scalardata(:,:)


  !initialize MPI
#ifdef CPP_MPI
  call MPI_Init ( ierr )
#endif
  call mympi_init()

  !initializations
  call pointgrp(rotmat,rotname)
  areatot = 0d0


  write(*,*) '****************************************************'
  write(*,*) '* What is the dimensionality of your BZ            *'
  write(*,*) '* - - - - - - - - - - - - - - - - - - - - - - - - -*'
  write(*,*) '* possible choices implemented: 2 or 3             *'
  read(*,*) nBZdim

  !read in the cubes header
  write(filename,fmt_fn_ext) filename_cubesinfo, ext_formatted
  open(unit=ifile,file=trim(filename),form='formatted',action='read')
  read(ifile,'(3I8)') nCub3
  close(ifile)

  !read the bounds
  open(unit=1359,file='bounds',form='formatted',action='read')
  read(1359,'(3ES25.16)') bounds
  close(1359)




  !#################################
  !# work on the visualization set #
  !#################################

  !read in the k-points data
  write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_vis, ext_orig
  call read_kpointsfile_vis(nkpts_orig, nkpts_irr_orig, kpoints_irr_orig, nsym_orig, isym_orig, kpt2irr_orig, irr2kpt_orig, filenamein=filename)
  write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_vis, ext_refined
  call read_kpointsfile_vis(nkpts_refined, nkpts_irr_refined, kpoints_irr_refined, nsym_refined, isym_refined, kpt2irr_refined, irr2kpt_refined, filenamein=filename)

  if(nsym_refined/=nsym_orig)then
    stop 'number of symmetries not consistent for vis'
  else
   if(any(isym_refined/=isym_orig)) stop 'symmetries not consistent for vis'
  end if

  !unfold the compressed data to a big array
  allocate(kpoints_refined(3,nkpts_refined), kpoints_orig(3,nkpts_orig), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating big kpoint-arrays'
  do ikp=1,nkpts_orig
    kpoints_orig(:,ikp) = kpoints_irr_orig(:,kpt2irr_orig(ikp))
  end do!ikp
  do ikp=1,nkpts_refined
    kpoints_refined(:,ikp) = kpoints_irr_refined(:,kpt2irr_refined(ikp))
  end do!ikp

  !compute the cube-ids for each triangle
  call find_cubeids_vis(nBZdim, nkpts_orig,    kpoints_orig,    bounds, nCub3, cubeids_orig,    areas_orig   )
  call find_cubeids_vis(nBZdim, nkpts_refined, kpoints_refined, bounds, nCub3, cubeids_refined, areas_refined)

  !go through the old set and delete triangles that are in the new set
! allocate(cubeids_orig_save(nkpts_orig/3))
! cubeids_orig_save = cubeids_orig
  do itri=1,nkpts_orig/nBZdim
    if(any(cubeids_refined==cubeids_orig(itri))) cubeids_orig(itri) = -1
  end do!itri

! if(myrank==master)then
!   write(*,*) 'deleted triangles:'
!   do itri=1,nkpts_orig/3
!     if(cubeids_orig(itri)==-1) then 
!       write(*,*) cubeids_orig_save(itri)
!     end if!cubeids_orig(itri)==-1
!   end do!itri
!   write(1002,*) cubeids_refined
! end if!myrank==master

  !allocate the new arrays
  nkpts_orig_filtered = nBZdim*count(cubeids_orig/=-1)
  nkpts_new = nkpts_orig_filtered + nkpts_refined
  if(myrank==master) write(*,'("#original= ", I0, " #additional= ",I0, " #deleted= ",I0, " #new= ",I0)') &
                     &          nkpts_orig/nBZdim,   nkpts_refined/nBZdim, count(cubeids_orig==-1),   nkpts_new/nBZdim
  allocate(kpoints_new(3,nkpts_new), cubeids_new(nkpts_new/nBZdim), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating kpoints_new etc.'

  !merge the two arrays
  call merge_arrays_vis()

  !save the visualization set
  call find_kpoints_irredset( bounds, nkpts_new, kpoints_new, nkpts_irr_new, kpt2irr_new, irr2kpt_new )
  if(myrank==master)then
    
    write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_vis, ext_new
    call save_kpointsfile_vis(nkpts_new, nkpts_irr_new, kpoints_new, nsym_orig, isym_orig, kpt2irr_new, irr2kpt_new, filenamein=filename )

    itmparr(1) = 1
    write(filename,fmt_fn_sub_ext) filename_vtktest, filemode_ref, ext_vtkxml
    call write_pointdata_rot( trim(filename),nkpts_new,kpoints_new, &
                            & 0,scalardata,scalarstring,            &
                            & 0,vectordata,vectorstring,            &
                            & 1,rotmat, itmparr                     )
    write(*,*) 'Visualization data written!'
    write(*,'(A,ES18.9)') 'Total area is ', areatot
  end if!myrank==master


  deallocate( isym_refined, kpt2irr_refined, irr2kpt_refined, cubeids_refined, &
            & isym_orig,    kpt2irr_orig,    irr2kpt_orig,    cubeids_orig,    &
            &               kpt2irr_new,     irr2kpt_new,     cubeids_new,     &
            & kpoints_irr_refined, kpoints_refined, areas_refined,             &
            & kpoints_irr_orig,    kpoints_orig,    areas_orig,                &
            &                      kpoints_new                                 )





  !###############################
  !# work on the integration set #
  !###############################

  !read in the k-points data
  write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_int, ext_orig
  call read_kpointsfile_int(nkpts_orig, kpoints_orig, areas_orig, nsym_orig, isym_orig, filename)
  write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_int, ext_refined
  call read_kpointsfile_int(nkpts_refined, kpoints_refined, areas_refined, nsym_refined, isym_refined, filename)

  if(nsym_refined/=nsym_orig)then
    stop 'number of symmetries not consistent for int'
  else
   if(any(isym_refined/=isym_orig)) stop 'symmetries not consistent for int'
  end if

  !compute the cubeids for each kpoint 
  call find_cubeids_int(nkpts_orig,    kpoints_orig,    bounds, nCub3, cubeids_orig   )
  call find_cubeids_int(nkpts_refined, kpoints_refined, bounds, nCub3, cubeids_refined)

  !go through the old set and delete points that are in the new set
! allocate(cubeids_orig_save(nkpts_orig))
! cubeids_orig_save = cubeids_orig
  do ikp=1,nkpts_orig
    if(any(cubeids_refined==cubeids_orig(ikp))) cubeids_orig(ikp) = -1
  end do!itri

! if(myrank==master)then
!   write(*,*) 'deleted points:'
!   do itri=1,nkpts_orig
!     if(cubeids_orig(itri)==-1) then 
!       write(*,*) cubeids_orig_save(itri)
!       write(*,'(8X,4ES25.16)') kpoints_orig(:,itri), areas_orig(itri)
!     end if!cubeids_orig(itri)==-1
!   end do!itri
!   write(1002,*) cubeids_refined
! end if!myrank==master

  !allocate the new arrays
  nkpts_orig_filtered = count(cubeids_orig/=-1)
  nkpts_new = nkpts_orig_filtered + nkpts_refined
  if(myrank==master) write(*,'("#original= ", I0, " #additional= ",I0, " #deleted= ",I0, " #new= ",I0)') &
                     &          nkpts_orig,     nkpts_refined,   count(cubeids_orig==-1),   nkpts_new
  allocate(kpoints_new(3,nkpts_new), cubeids_new(nkpts_new), areas_new(nkpts_new), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating kpoints_new etc.'

  !merge the two arrays
  call merge_arrays_int()

  if(myrank==master)then
    
    write(filename,fmt_fn_sub_ext) filename_fsdata, filemode_int, ext_new
    call save_kpointsfile_int(nkpts_new, kpoints_new, areas_new, nsym_orig, isym_orig, filename)

    write(*,*) 'Integration data written!'
    write(*,'(A,ES18.9)') 'Total area is ', sum(areas_new)
  end if!myrank==master

#ifdef CPP_MPI
  call MPI_Finalize(ierr)
#endif


contains


  subroutine merge_arrays_int()
    implicit none

    integer :: ikp, ipointer

    ipointer = 0
    do ikp=1,nkpts_orig
      if(cubeids_orig(ikp)==-1)then
        write(*,'(A,3ES16.9)') 'skipping cube: k=', kpoints_orig(:,ikp)
        cycle
      end if
      ipointer = ipointer+1
      kpoints_new(:,ipointer) = kpoints_orig(:,ikp)
      cubeids_new(ipointer)   = cubeids_orig(ikp)
      areas_new(ipointer)     = areas_orig(ikp)
    end do!itri

    if(ipointer/=nkpts_orig_filtered) stop 'inconsistency in number of kpoints when filtering old array, int'

    kpoints_new(:,ipointer+1:nkpts_new) = kpoints_refined
    cubeids_new(ipointer+1:nkpts_new)   = cubeids_refined
    areas_new(ipointer+1:nkpts_new)     = areas_refined

  end subroutine merge_arrays_int


  subroutine merge_arrays_vis()
    implicit none

    integer :: itri, ipointer, istart1, istart2, istop1, istop2

    areatot = 0d0

    ipointer = 0
    do itri=1,nkpts_orig/nBZdim
      if(cubeids_orig(itri)==-1) cycle
      ipointer = ipointer+1
      istart1 = (ipointer-1)*nBZdim+1
      istop1  = ipointer*nBZdim
      istart2 = (itri-1)*nBZdim+1
      istop2  = itri*nBZdim
      kpoints_new(:,istart1:istop1) = kpoints_orig(:,istart2:istop2)
      cubeids_new(ipointer) = cubeids_orig(itri)
      areatot = areatot + areas_orig(itri)
    end do!itri

    if(ipointer/=nkpts_orig_filtered/nBZdim) stop 'inconsistency in number of triangles when filtering old array, vis'

    kpoints_new(:,nBZdim*ipointer+1:nkpts_new) = kpoints_refined
    cubeids_new(ipointer+1:nkpts_new/nBZdim)   = cubeids_refined
    areatot = areatot + sum(areas_refined)

  end subroutine merge_arrays_vis


  subroutine find_cubeids_int(nkpts_in, kpoints_in, bounds, nCub3, cubeids)
    implicit none
    integer,                       intent(in)  :: nkpts_in, nCub3(3)
    double precision,              intent(in)  :: kpoints_in(3,nkpts_in), bounds(3,2)
    integer,          allocatable, intent(out) :: cubeids(:)

    integer :: ierr, ikp, ii3(3)

    allocate(cubeids(nkpts_in), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating cubeids in int'

    do ikp=1,nkpts_in
      call kpoint_to_indices(nCub3,bounds,kpoints_in(:,ikp),ii3,cubeids(ikp))
    end do!itri

  end subroutine find_cubeids_int


  subroutine find_cubeids_vis(nBZdim, nkpts_in, kpoints_in, bounds, nCub3, cubeids, areas)
    implicit none
    integer,                       intent(in)  :: nBZdim, nkpts_in, nCub3(3)
    double precision,              intent(in)  :: kpoints_in(3,nkpts_in), bounds(3,2)
    integer,          allocatable, intent(out) :: cubeids(:)
    double precision, allocatable, intent(out) :: areas(:)

    integer :: ierr, itri, ik1, ik2, ikp, ii3(3)
    double precision :: ktriangle(3,nBZdim), kmedium(3)

    allocate(cubeids(nkpts_in/nBZdim), areas(nkpts_in/nBZdim), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating cubeids in vis'

    do itri=1,nkpts_in/nBZdim
      ik1=(itri-1)*nBZdim+1
      ik2=itri*nBZdim
      ktriangle = kpoints_in(:,ik1:ik2)

      if(nBZdim==3)then
        areas(itri) = area_triangle(ktriangle)
      elseif(nBZdim==2)then
        areas(itri) = sqrt( sum((ktriangle(:,2)-ktriangle(:,1))**2) )
      else
        stop 'nBZdim shall be 2 or 3'
      endif

      kmedium = 0
      do ikp=1,nBZdim
        kmedium = kmedium + ktriangle(:,ikp)
      end do!ikp
      kmedium = kmedium/nBZdim

      call kpoint_to_indices(nCub3,bounds,kmedium,ii3,cubeids(itri))

    end do!itri

  end subroutine find_cubeids_vis



  subroutine kpoint_to_indices(nCub3,bounds,kpoint,ii3,ixyz)
    implicit none
    integer,          intent(in) :: nCub3(3)
    double precision, intent(in) :: bounds(3,2), kpoint(3)
    integer,          intent(out) :: ii3(3), ixyz

    double precision :: ktemp(3)

    ktemp = (kpoint - bounds(:,1))/(bounds(:,2) - bounds(:,1))
    ii3 = int(ktemp*nCub3+1d-14)

    ixyz = 1 + ii3(1) + nCub3(1)*ii3(2) + nCub3(1)*nCub3(2)*ii3(3)

  end subroutine kpoint_to_indices



  double precision function area_triangle(kpoints)

    use mod_mathtools, only: crossprod

    double precision, intent(in) :: kpoints(3,3)
    double precision :: k21(3), k31(3), kcross(3)

    k21 = kpoints(:,3) - kpoints(:,1)
    k31 = kpoints(:,3) - kpoints(:,2)
    call crossprod(k21, k31, kcross)

    area_triangle = 0.5d0*sqrt(sum(kcross**2))

  end function area_triangle

end program mergerefined
