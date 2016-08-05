!-----------------------------------------------------------------------------------------!
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of kk-prime@juKKR and available as free software under the conditions !
! of the MIT license as expressed in the LICENSE file in more detail.                     !
!-----------------------------------------------------------------------------------------!


program translate_cubesfile

  use mod_mathtools, only: bubblesort_int
  use mod_ioformat,  only: filename_cubesinfo, filename_cubesrefine, ext_formatted, fmt_fn_ext
  implicit none

  integer :: nCub3(3), ii, iptr, curval, nlines, nmarked, ierr
  integer, allocatable :: cubesin(:), sortind(:), cubes_filtered(:), imarked(:)
  integer, parameter :: ifile=14
  double precision :: bounds(3,2)

  !read in the cubes header
  write(filename,fmt_fn_ext) filename_cubesinfo, ext_formatted
  open(unit=ifile,file=trim(filename),form='formatted',action='read')
  read(ifile,'(3I8)') nCub3
  close(ifile)

  !read in the cubesoutput file
  open(unit=ifile,file='cubes_inp',form='formatted',action='read')
  nlines = get_nLines(ifile)
  allocate(cubesin(nlines), sortind(nlines), cubes_filtered(nlines), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating cubesin'
  rewind(ifile)
  do ii=1,nlines
    read(ifile,*) cubesin(ii)
  end do!ii
  close(ifile)

  !filter the cubes
  call bubblesort_int(nlines, cubesin, sortind)
  iptr=1
  curval = cubesin(sortind(1))
  cubes_filtered(1) = cubesin(sortind(1))
  do ii=2,nlines
    if(cubesin(sortind(ii))>curval)then
      iptr=iptr+1
      cubes_filtered(iptr) = cubesin(sortind(ii))
      curval = cubesin(sortind(ii))
    end if
  end do
  nmarked = iptr
  allocate(imarked(nmarked), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating imarked'
  imarked = cubes_filtered(1:iptr)
  deallocate(cubes_filtered, cubesin, sortind)

  !save the filtered cubesfile
  open(unit=ifile,file='cubesinfo.new',form='formatted',action='write',)
  write(ifile,'(3I0)') nCub3
  write(ifile,'(I0)') nmarked
  write(ifile,'(10(I0,1X))') imarked
  close(ifile)

  open(unit=1359,file='bounds',form='formatted',action='read')
  read(1359,'(3ES25.16)') bounds
  close(1359)

  call cubes2VTK('lost_cubes.vtp', nCub3, nmarked, imarked, bounds)


contains


  integer function get_nFiles()
    use mod_ioformat, only: fmt_fn_rank_ext, filename_outinfo, ext_formatted
    implicit none

    integer :: nFiles
    logical :: fileexists
    character(len=80) :: filename


    nFiles=0
    fileexists = .TRUE.
    do while (fileexists)
      write(filename,fmt_fn_rank_ext) filename_outinfo, nFiles, ext_formatted
      inquire(file=filename, exist=fileexists)
      nFiles = nFiles+1
    end do
    get_nFiles = nFiles-1
    write(*,'(A,I0)') 'Number of output-files found = ', get_nFiles
  end function get_nFiles


  function get_nLines( UIO ) result( nLines )
  implicit none
  integer :: nLines, IOERR, UIO
  character(len=50) :: tmpstring
    nLines = 0
    Read_Loop: DO
      READ( UIO, *, IOSTAT = IOERR ) tmpstring
      IF (IOERR < 0) THEN
        EXIT Read_Loop
      END IF
      nLines = nLines + 1
    END DO Read_Loop
  end function get_nLines

  subroutine cubes2VTK(filename,nCub3,nmarked,imarked,bounds)
    use mod_vtkxml
    implicit none

    character(len=*), intent(in) :: filename
    integer,          intent(in) :: nCub3(3), nmarked, imarked(nmarked)
    double precision, intent(in) :: bounds(3,2)

    integer          :: faces(4,6), offsets(6)
    double precision :: kverts(3,8)
    integer :: icub, ivert, iface
    character(len=80) :: fmtstr

    open(unit=123894, file=trim(filename), action='write', form='formatted')

    !write header
    write(123894, FMT=vtkfmt_ivtkfile)
    write(123894, FMT=vtkfmt_ipolydata)
    write(123894, FMT=vtkfmt_ipiece) 8*nmarked, 6*nmarked

    !write points
    write(123894, FMT=vtkfmt_ipoints)
    write(123894, FMT=vtkfmt_idata_points)
    write(fmtstr,'("(",I0,"X,3ES14.5)")') vtkfmxXdata
    do icub=1,nmarked
      call generate_cubevertices(nCub3,imarked(icub),bounds,kverts)
      do ivert=1,8
        write(123894, FMT=fmtstr) kverts(:,ivert)
      end do!ivert
    end do!icub
    write(123894, FMT=vtkfmt_fdata)
    write(123894, FMT=vtkfmt_fpoints)

    faces(:,1) = (/ 0,1,3,2 /)
    faces(:,2) = (/ 0,1,5,4 /)
    faces(:,3) = (/ 0,2,6,4 /)
    faces(:,4) = (/ 1,3,7,5 /)
    faces(:,5) = (/ 2,3,7,6 /)
    faces(:,6) = (/ 4,5,7,6 /)

    offsets = (/ 4, 8, 12, 16, 20, 24 /)

    !write polys
    write(123894, FMT=vtkfmt_ipolys)
    write(123894, FMT=vtkfmt_idata_connectivity)
    do icub=1,nmarked
      do iface=1,6
        write(123894, FMT='(10X,4I8)' ) faces(:,iface)+(icub-1)*8
      end do!iface
    end do!icub
    write(123894, FMT=vtkfmt_fdata)
    write(123894, FMT=vtkfmt_idata_offsets)
    write(fmtstr, '("(",I0,"X,6(I0,X))")') vtkfmxXdata
    do icub=1,nmarked
      write(123894, FMT=fmtstr) offsets + (icub-1)*24
    end do!icub
    write(123894, FMT=vtkfmt_fdata)
    write(123894, FMT=vtkfmt_fpolys)

    !write tail
    write(123894, FMT=vtkfmt_fpiece)
    write(123894, FMT=vtkfmt_fpolydata)
    write(123894, FMT=vtkfmt_fvtkfile)

    close(123894)

  end subroutine cubes2VTK


  subroutine generate_cubevertices(nCub3,cubeid,bounds,kverts)
    !generates the vertices of a cube given by the cubeid

    implicit none

    integer, intent(in) :: nCub3(3), cubeid
    double precision, intent(in)  :: bounds(3,2)
    double precision, intent(out) :: kverts(3,8)

    integer :: ii3(3), ivert

    call unroll_ixyz(cubeid, nCub3, ii3)

    !generate vertices, limits are from 0 to 1
    kverts(:,1) = dble( (/ ii3(1)  , ii3(2)  , ii3(3)   /) )/nCub3
    kverts(:,2) = dble( (/ ii3(1)+1, ii3(2)  , ii3(3)   /) )/nCub3
    kverts(:,3) = dble( (/ ii3(1)  , ii3(2)+1, ii3(3)   /) )/nCub3
    kverts(:,4) = dble( (/ ii3(1)+1, ii3(2)+1, ii3(3)   /) )/nCub3
    kverts(:,5) = dble( (/ ii3(1)  , ii3(2)  , ii3(3)+1 /) )/nCub3
    kverts(:,6) = dble( (/ ii3(1)+1, ii3(2)  , ii3(3)+1 /) )/nCub3
    kverts(:,7) = dble( (/ ii3(1)  , ii3(2)+1, ii3(3)+1 /) )/nCub3
    kverts(:,8) = dble( (/ ii3(1)+1, ii3(2)+1, ii3(3)+1 /) )/nCub3

    !shift the vertices to the correct bounds
    do ivert=1,8
      kverts(:,ivert) = kverts(:,ivert)*(bounds(:,2)-bounds(:,1))+bounds(:,1)
    end do!ivert

  end subroutine generate_cubevertices


  subroutine unroll_ixyz(ixyz, nCub3, ii3)
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !+ helper function to divide the UID to an integer triple +
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    implicit none

    integer, intent(in)  :: ixyz, nCub3(3)
    integer, intent(out) :: ii3(3)

    integer :: irest

    ii3(3) = int((ixyz-1)/(nCub3(1)*nCub3(2)))
    irest= ixyz-1-ii3(3)*(nCub3(1)*nCub3(2))
    ii3(2) = int(irest/nCub3(1))
    irest= ixyz-1-ii3(3)*(nCub3(1)*nCub3(2)) - ii3(2)*nCub3(1)
    ii3(1) = irest

  end subroutine unroll_ixyz

end program translate_cubesfile
