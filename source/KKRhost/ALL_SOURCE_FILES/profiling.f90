! -------------------------------------------------------------------------------
! MODULE: Profiling

! DESCRIPTION:
!> @brief Subroutine to handle memory profiling
!> @details Allows one to track the memory consumption by making use of the
! memocc()
!> subroutine

!> @author
!> Jonathan Chico
!> @date 10.12.2017
! -------------------------------------------------------------------------------
module mod_profiling

  use :: mod_datatypes
  implicit none

contains

  ! control the memory occupation by calculating the overall size in bytes of
  ! the allocated arrays
  ! usage:
  ! when allocating allocating an array "stuff" of dimension n in the routine
  ! "dosome"
  ! allocate(stuff(n),stat=i_stat)
  ! call memocc(i_stat,product(shape(stuff))*kind(stuff),'stuff','dosome')
  ! when deallocating
  ! i_all=-product(shape(stuff))*kind(stuff)
  ! deallocate(stuff,stat=i_stat)
  ! call memocc(i_stat,i_all,'stuff','dosome')
  ! the counters are initialized with
  ! call memocc(0,iproc,'count','start') (iproc = mpi rank, nproc=mpi size)
  ! and stopped with
  ! call memocc(0,0,'count','stop')
  ! at the end of the calculation a short report is printed on the screen
  ! some information can be also written on disk following the needs
  ! This file is distributed under the terms of the
  ! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
  ! Copyright (C) Luigi Genovese, CEA Grenoble, France, 2007
  !> Memory profiling routine
  subroutine memocc(istat, isize, array, routine)
    use :: mod_types, only: t_inc
    implicit none
    character (len=*), intent (in) :: array
    character (len=*), intent (in) :: routine
    integer, intent (in) :: istat
    integer, intent (in) :: isize

    ! Local variables
    character (len=20) :: filename
    character (len=36) :: maxroutine, locroutine
    character (len=36) :: maxarray                             ! , locarray

    integer :: nalloc, ndealloc, locpeak, locmemory, iproc
    integer :: dblsize, mfileno
    integer (kind=di) :: memory, maxmemory
    character (len=1) :: allocationflag

    save :: memory, nalloc, ndealloc, maxroutine, maxarray, maxmemory
    save :: locroutine, locpeak, locmemory, iproc ! , locarray

    mfileno = 77
    dblsize = 1

    if (t_inc%i_write>0) then

      write (filename, '(A)') 'meminfo.txt'

      select case (array)
      case ('count')
        if (routine=='start') then
          memory = 0
          maxmemory = 0
          nalloc = 0
          ndealloc = 0
          locroutine = 'routine'
          ! locarray = 'array'
          locmemory = 0
          locpeak = 0
          iproc = isize
          ! open the writing file for the root process
          if (iproc==0) then
            open (unit=mfileno, file=trim(filename), status='unknown')
            write (mfileno, '(a32,1x,a20,3(1x,a12))') '(Data in kB)             Routine', '    Peak Array', 'Routine Mem', 'Total Mem', 'Action'
          end if
        else if (routine=='stop' .and. iproc==0) then
          ! write(mfileno,'(a32,a16,4(1x,i12))')&
          ! trim(locroutine),trim(locarray),locmemory/1024,locpeak/1024,memory/1024,&
          ! (memory+locpeak-locmemory)/1024
          ! close(mfileno)
          write (*, '(1x,a)') '-------------------MEMORY CONSUMPTION REPORT-----------------------'
          write (*, '(1x,2(i0,a),i0)') nalloc, ' allocations and ', ndealloc, ' deallocations, remaining memory(B):', memory
          write (*, '(1x,a,i0,a)') 'Memory occupation peak: ', maxmemory/1024, ' kB'
          write (*, '(1x,5(a))') 'For the array: "', trim(maxarray), '" in routine "', trim(maxroutine), '"'
          write (*, '(1x,a)') '-----------------END MEMORY CONSUMPTION REPORT---------------------'
        end if

      case default
        ! control of the allocation/deallocation status
        if (istat/=0) then
          if (isize>=0) then
            write (*, *) ' subroutine ', routine, ': problem of allocation of array ', array
            stop
          else if (isize<0) then
            write (*, *) ' subroutine ', routine, ': problem of deallocation of array ', array
            stop
          end if
        end if
        select case (iproc)
        case (0)
          ! To be used for inspecting an array which is not deallocated
          if (isize>0) then
            allocationflag = 'A'
          else if (isize<0) then
            allocationflag = 'D'
          else
            allocationflag = '?'
          end if

          write (mfileno, '(a32,1x,a20,3(1x,i12),1x,a12)') trim(routine), trim(array), isize*dblsize, memory, maxmemory, allocationflag
          if (trim(locroutine)/=routine) then
            ! write(mfileno,'(a32,a14,4(1x,i12))')&
            ! trim(locroutine),trim(locarray),locmemory/1024,locpeak/1024,memory/1024,&
            ! (memory+locpeak-locmemory)/1024
            locroutine = routine
            locmemory = isize*dblsize
            locpeak = isize*dblsize
          else
            locmemory = locmemory + isize*dblsize
            if (locmemory>locpeak) then
              locpeak = locmemory
            end if
          end if
          ! locarray = array
          memory = memory + isize*dblsize
          if (memory>maxmemory) then
            maxmemory = memory
            maxroutine = routine
            maxarray = array
          end if
          if (isize>0) then
            nalloc = nalloc + 1
          else if (isize<0) then
            ndealloc = ndealloc + 1
          end if
        case default
          return
        end select

      end select

    end if                         ! (t_inc%i_write>0)


  end subroutine memocc


end module mod_profiling
