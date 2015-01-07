program refineBZparts

  use mod_mathtools, only: bubblesort_int
  use mod_ioinput,   only: IoInput
  use mod_ioformat,  only: filename_cubesinfo, filename_cubesrefine, ext_formatted, fmt_fn_ext

  implicit none

  integer :: nCub3(3), nCub3_inp(3), nCub3_file(3), nFSiter, ii3(3), imarked(1), nmarked, nverts, nBZdimen=3

  integer, allocatable :: nCut_iter(:), nCub3_steps(:,:)
  double precision :: bounds(3,2)

  integer :: idomain, ierr, istep, iselect
  character(len=80)    :: uio
  character(len=256) :: filename
  integer, parameter :: ifile=14, iofile=16

  write(*,*) '****************************************************'
  write(*,*) '* What is the dimensionality of your BZ            *'
  write(*,*) '* - - - - - - - - - - - - - - - - - - - - - - - - -*'
  write(*,*) '* possible choices implemented: 2 or 3             *'
  read(*,*) nBZdimen
  select case (nBZdimen)
    case (3); nverts=8
    case (2); nverts=4
    case default; stop 'nBZdimen =/= 2 or 3'
  end select

  write(filename,fmt_fn_ext) filename_cubesrefine, ext_formatted
  open(unit=iofile,file=filename,form='formatted',action='write',access='append')


  !read in the cubes header
  write(filename,fmt_fn_ext) filename_cubesinfo, ext_formatted
  open(unit=ifile,file=trim(filename),form='formatted',action='read')
  read(ifile,'(3I8)') nCub3_file
  close(ifile)

  write(*,'(A)') 'BZdivisions from cubesinfo.txt:'
  write(*,'(3I8)') nCub3_file
  nCub3 = nCub3_file

  !read the bounds
  open(unit=1359,file='bounds',form='formatted',action='read')
  read(1359,'(3ES25.16)') bounds
  close(1359)

  !*************************************
  !*** get infos from the input card ***
  !*************************************
  call IoInput('NKPTCUBES ',uio,1,7,ierr)
  read(unit=uio,fmt=*) nCub3_inp(:)

  call IoInput('NFSITER   ',uio,1,7,ierr)
  read(unit=uio,fmt=*) nFSiter

  allocate(nCut_iter(nFSiter), STAT=ierr)
  if(ierr/=0) stop 'Problem allocating nCut_iter etc.'
  call IoInput('NREFINE   ',uio,1,7,ierr)
  read(unit=uio,fmt=*) nCut_iter


  !***********************************
  !*** compute the divisions of BZ ***
  !***********************************
  allocate(nCub3_steps(3,nFSiter+1),STAT=ierr)
  if(ierr/=0) stop 'Problem allocating nCub3_steps'
  nCub3_steps(:,1) = nCub3_inp
  do istep=1,nFSiter
    nCub3_steps(:,istep+1) = nCub3_steps(:,istep)*nCut_iter(istep)
    if(nBZdimen==2) nCub3_steps(3,istep+1)=1
  end do


  !loop to input multiple choices
  idomain=0
  do while(idomain==0)

    write(*,*) '****************************************************'
    write(*,*) '* What do you want to do?                          *'
    write(*,*) '* - - - - - - - - - - - - - - - - - - - - - - - - -*'
    write(*,*) '* 0) exit program                                  *'
    write(*,*) '* 1) change the cube size                          *'
    write(*,*) '* 2) add a single cube containing a single point   *'
    write(*,*) '* 3) add multiple cubes in a sphere around a point *'
    write(*,*) '* 4) include the incomplete cubes from outfiles    *'
    read(*,*) iselect

    select case(iselect)
      case(0); idomain=1
      case(1)
        nCub3 = select_cubesize(nFSiter, nCub3_steps, bounds)
      case(2); call add_point ( nCub3, nCub3_file, nFSiter, nCub3_steps, bounds )
      case(3); call add_sphere( nverts, nCub3, nCub3_file, nFSiter, nCub3_steps, bounds )
      case(4); call add_from_outfile(nCub3_file, bounds)
      case default; write(*,*) '--> illegal input option. Try again.'
    end select

  end do!while

  close(iofile)

contains

  subroutine add_from_outfile(nCub3,bounds)
    use mod_ioformat, only: fmt_fn_rank_ext, filename_outinfo, ext_formatted
    implicit none

    integer, intent(in) :: nCub3(3)
    double precision, intent(in) :: bounds(3,2)

    integer :: nFiles, nmarked, iselect, ipt, ii, ifile, iline, icube
    integer, allocatable :: nLines(:), imarked(:)
    character(len=256) :: filename

    integer, parameter :: iounit=6592

    nFiles = get_nFiles()
    allocate(nLines(nFiles))

    do ifile=1,nFiles
      write(filename,fmt_fn_rank_ext) filename_outinfo, ifile-1, ext_formatted
      open(iounit, file=trim(filename), form='formatted', action='read')
      nLines(ifile) = get_nLines( iounit )
      write(*,'(A,A,I0)') trim(filename), ' --- ', nLines(ifile)
      close(iounit)
    end do

    nmarked=sum(nLines)
    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked'

    ipt=0
    do ifile=1,nFiles
      write(filename,fmt_fn_rank_ext) filename_outinfo, ifile-1, ext_formatted
      open(iounit, file=trim(filename), form='formatted', action='read')
      do iline=1,nLines(ifile)
        ipt=ipt+1
        read(iounit,*) imarked(ipt)
      end do!iline
      close(iounit)
    end do!ifile

    if(ipt/=nmarked) then
      write(*,'(I0,A,I0)') ipt, ' = ipt does not equal nmarked = ', nmarked
      stop 'Data inconsitency with outfiles'
    end if

    select case (nBZdimen)
      case (3); call cubes2VTK('test_cube.vtp',     nCub3, nmarked, imarked, bounds)
      case (2); call squares2TXT('test_square.txt', nCub3, nmarked, imarked, bounds)
      case default; stop 'nBZdimen =/= 2 or 3'
    end select

    write(*,'(A,I0,A,I0)') ' .. number of cubes found in the ', nFiles, ' files =', nmarked
    write(*,*) '***************************'
    write(*,*) '* Add these points ?      *'
    write(*,*) '* 0: no                   *'
    write(*,*) '* 1: yes                  *'
    read(*,*) iselect

    select case(iselect)
      case(0)
        write(*,*) ' ... Cubes ignored!'
      case(1)
        do ii=1,nmarked
          write(iofile,'(5I8)') nCub3, 1, imarked(ii)
        end do!ii
        write(*,'(A,I0,A)') '  ... ', ii-1, ' cubes added to refinement-file'
      case default
        write(*,*) ' ... This is was illegal input option.'
    end select

  end subroutine add_from_outfile



  subroutine add_point(nCub3, nCub3_file, nFSiter, nCub3_steps, bounds )

    implicit none

    integer,          intent(in) :: nCub3(3), nCub3_file(3), nFSiter, nCub3_steps(3,nFSiter+1)
    double precision, intent(in) :: bounds(3,2)

    integer :: ido, nmarked, imarked(1), ii3(3), nCub3tmp(3), iselect
    double precision :: kpoint(3)

    nCub3tmp = nCub3

    !**************************
    !*** enter the k-point ***!
    !**************************
    write(*,*) '***************************'
    write(*,*) '* Enter the k-point:      *'
    read(*,*) kpoint
    nmarked = 1

    ido=0
    do while(ido==0) 
      call kpoint_to_indices(nCub3tmp,bounds,kpoint,ii3,imarked(1))

      !**************************
      !*** visualize the cube ***
      !**************************
      select case (nBZdimen)
        case (3); call cubes2VTK('test_cube.vtp',     nCub3tmp, nmarked, imarked, bounds)
        case (2); call squares2TXT('test_square.txt', nCub3tmp, nmarked, imarked, bounds)
        case default; stop 'nBZdimen =/= 2 or 3'
      end select

      !*********************************
      !*** Ask user about what to do ***
      !*********************************
      write(*,*) '.. number of this cube: imarked=', imarked
      write(*,*) '.. cube written to file test_cube.vtp'
      write(*,*) ''
      write(*,*) '*************************************************'
      write(*,*) '* Shall I add this cube to the refinement-file? *'
      write(*,*) '* 0: no                                         *'
      write(*,*) '* 1: yes                                        *'
      write(*,*) '* 2: change cubesize and try again              *'
      read(*,*) iselect

      select case(iselect)
        case(0)
          write(*,*) ' ... Cube ignored!'
          ido=1
        case(1)
          write(iofile,'(5I8)') nCub3tmp, nCub3_file(1)/nCub3tmp(1), imarked(1)
          write(*,*) ' ... Cube added to refinement-file'
          ido=1
        case(2)
          !change size of cubes
          nCub3tmp = select_cubesize(nFSiter, nCub3_steps, bounds)
          ido=0
        case default
          write(*,*) ' ... This is an illegal input option.'
          write(*,*) ' ..... Your input was: ', iselect
          write(*,*) ' ..... Try again!'
          ido=0
      end select


    end do!while ido==0

  end subroutine add_point



  subroutine add_sphere(nverts,nCub3,nCub3_file,nFSiter,nCub3_steps,bounds)
    implicit none

    integer,          intent(in) :: nverts, nCub3(3), nCub3_file(3), nFSiter, nCub3_steps(3,nFSiter+1)
    double precision, intent(in) :: bounds(3,2)

    integer :: nmarked, ii, ido, nCub3tmp(3)
    double precision :: kcenter(3), radius
    integer, allocatable :: imarked(:)

    nCub3tmp = nCub3

    write(*,*) nCub3tmp
    !**************************
    !*** enter the sphere ***!
    !**************************
    write(*,*) '***************************************'
    write(*,*) '* Enter the center of the sphere      *'
    read(*,*) kcenter
    write(*,*) '* Enter the radius of the sphere      *'
    read(*,*) radius

    ido=0
    do while(ido==0) 
      call mark_cubes_in_sphere(nverts, nCub3tmp,bounds,kcenter,radius,nmarked,imarked)

      !**************************
      !*** visualize the cube ***
      !**************************
      select case (nBZdimen)
        case (3); call cubes2VTK('test_cube.vtp',     nCub3tmp, nmarked, imarked, bounds)
        case (2); call squares2TXT('test_square.txt', nCub3tmp, nmarked, imarked, bounds)
        case default; stop 'nBZdimen =/= 2 or 3'
      end select


      !*********************************
      !*** Ask user about what to do ***
      !*********************************
      write(*,*) '.. number of cubes found = ', nmarked
      write(*,*) '.. cubes written to file test_cube.vtp'
      write(*,*) ''
      write(*,*) '***************************************************'
      write(*,*) '* Shall I add these cubes to the refinement-file? *'
      write(*,*) '* 0: no                                           *'
      write(*,*) '* 1: yes                                          *'
      write(*,*) '* 2: change cubesize and try again                *'
      read(*,*) iselect

      select case(iselect)
        case(0)
          write(*,*) ' ... Cube ignored!'
          ido=1
        case(1)
          do ii=1,nmarked
            write(iofile,'(5I8)') nCub3tmp, nCub3_file(1)/nCub3tmp(1), imarked(ii)
          end do!ii
          write(*,'(A,I0,A)') '  ... ', ii-1, ' cubes added to refinement-file'
          ido=1
        case(2)
          !change size of cubes
          nCub3tmp = select_cubesize(nFSiter, nCub3_steps, bounds)
          ido=0
        case default
          write(*,*) ' ... This is an illegal input option.'
          write(*,*) ' ..... Your input was: ', iselect
          write(*,*) ' ..... Try again!'
          ido=0
      end select


    end do!while ido==0

  end subroutine add_sphere



  function select_cubesize(nFSiter, nCub3_steps, bounds) result(nCub3)
    implicit none

    integer,          intent(in) :: nFSiter, nCub3_steps(3,nFSiter+1)
    double precision, intent(in) :: bounds(3,2)
    integer :: nCub3(3)
    integer :: ido, iselect, istep

    write(*,*) '**********************************************'
    write(*,*) '* The BZ-Divisios for the steps are:         *'
    write(*,*) '*                                            *'
    do istep=1,nFSiter+1
      write(*,'(" * step= ",I2," -> nCub= ",3I8," *")') istep, nCub3_steps(:,istep)
    end do
    write(*,*) '* - - - - - - - - - - - - - - - - - - - - - -*'
    write(*,*) '* Enter the cubesize you want to use.        *'
    write(*,*) '*   [Press 0 to exit program]                *'

    ido=0
    do while (ido==0)
      read(*,*) iselect
      if(iselect>nFSiter+1)then
        write(*,'(A,I0)') ' ... Only integers between 1 and ', nFSiter+1, ' are allowed! Try again:'
      elseif(iselect==0)then
        stop '!!! Aborted by user !!!'
      else
        nCub3 = nCub3_steps(:,iselect)
        write(*,'(A,3I8)') '... your selected cubesize: ', nCub3
        write(*,'(A,3F8.4)') '... this corresponds to cubelengths of: ', (bounds(:,2)-bounds(:,1))/nCub3
        ido=1
      end if
    end do!while

  end function select_cubesize



  subroutine kpoint_to_indices(nCub3,bounds,kpoint,ii3,ixyz)
    implicit none

    integer,          intent(in) :: nCub3(3)
    double precision, intent(in) :: bounds(3,2), kpoint(3)
    integer,          intent(out) :: ii3(3), ixyz

    double precision :: ktemp(3)

    ktemp = (kpoint - bounds(:,1))/(bounds(:,2) - bounds(:,1))
    ii3 = int(ktemp*nCub3)

    ixyz = 1 + ii3(1) + nCub3(1)*ii3(2) + nCub3(1)*nCub3(2)*ii3(3)

  end subroutine kpoint_to_indices



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


  subroutine squares2TXT(filename,nCub3,nmarked,imarked,bounds)
!    use mod_fermisurf_basic, only: generate_squarevertices
    implicit none

    character(len=*), intent(in) :: filename
    integer,          intent(in) :: nCub3(3), nmarked, imarked(nmarked)
    double precision, intent(in) :: bounds(3,2)

    double precision :: kverts(3,4)
    integer :: isqu, ivert
    character(len=80) :: fmtstr

    open(unit=123895, file=trim(filename), action='write', form='formatted')

    !write header
    write(123895,'(A)') '# four lines define one square'

    !write points
    do isqu=1,nmarked
      call generate_squarevertices(nCub3,imarked(isqu),bounds,kverts)
      do ivert=1,4
        write(123895,'(3ES25.16)') kverts(:,ivert)
      end do!ivert
    end do!icub

    close(123895)

  end subroutine squares2TXT

  subroutine mark_cubes_in_sphere(nverts,nCub3,bounds,kcenter,radius,nmarked,imarked)
    implicit none
    integer,          intent(in) :: nverts, nCub3(3)
    double precision, intent(in) :: bounds(3,2), kcenter(3), radius
    integer,              intent(out) :: nmarked
    integer, allocatable, intent(out) :: imarked(:)

    integer :: ierr, icube, ntotmax
    double precision :: kverts(3,nverts)
    integer, allocatable :: imarkedtmp(:)

    ntotmax = product(nCub3)
    allocate(imarkedtmp(ntotmax), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarkedtmp(ntotmax) in mark_cubes_in_sphere'

    nmarked=0
    do icube=1,ntotmax
      select case (nBZdimen)
        case (3); call generate_cubevertices(nCub3,icube,bounds,kverts)
        case (2); call generate_squarevertices(nCub3,icube,bounds,kverts)
        case default; stop 'nBZdimen =/= 2 or 3'
      end select

      if(cube_in_sphere(nverts,kverts,kcenter,radius)) then
        nmarked = nmarked+1
        imarkedtmp(nmarked) = icube
      end if
    end do!icube

    if(nmarked==0) write(*,*) 'WARNING: no cubes in sphere found'

    allocate(imarked(nmarked), STAT=ierr)
    if(ierr/=0) stop 'Problem allocating imarked(nmarked) in mark_cubes_in_sphere'
    imarked = imarkedtmp(1:nmarked)

    deallocate(imarkedtmp)

  end subroutine mark_cubes_in_sphere



  logical function cube_in_sphere(nverts,kverts,kcenter,radius)
    implicit none
    integer,          intent(in) :: nverts
    double precision, intent(in) :: kverts(3,nverts), kcenter(3), radius

    integer :: ipoint
    double precision :: dtmp, rsquare
    logical :: points_in_sphere(nverts)

    rsquare = radius**2

    do ipoint=1,nverts
      dtmp = sum((kverts(:,ipoint) - kcenter(:))**2)
      points_in_sphere(ipoint) = (dtmp<rsquare)
    end do!ipoint
    cube_in_sphere = all(points_in_sphere)

  end function cube_in_sphere



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



  subroutine generate_squarevertices(nCub3,cubeid,bounds,kverts)
    !generates the vertices of a square given by the id
    implicit none

    integer, intent(in) :: nCub3(3), cubeid
    double precision, intent(in)  :: bounds(3,2)
    double precision, intent(out) :: kverts(3,4)

    integer :: ii3(3), ivert

    call unroll_ixyz(cubeid, nCub3, ii3)

    !generate vertices, limits are from 0 to 1
    kverts(:,1) = dble( (/ ii3(1)  , ii3(2)  , 0 /) )/nCub3
    kverts(:,2) = dble( (/ ii3(1)+1, ii3(2)  , 0 /) )/nCub3
    kverts(:,3) = dble( (/ ii3(1)  , ii3(2)+1, 0 /) )/nCub3
    kverts(:,4) = dble( (/ ii3(1)+1, ii3(2)+1, 0 /) )/nCub3

    !shift the vertices to the correct bounds
    do ivert=1,4
      kverts(:,ivert) = kverts(:,ivert)*(bounds(:,2)-bounds(:,1))+bounds(:,1)
    end do!ivert

  end subroutine generate_squarevertices



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

end program refineBZparts
