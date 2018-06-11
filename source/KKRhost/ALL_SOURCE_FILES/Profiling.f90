!-------------------------------------------------------------------------------
! MODULE: Profiling
!
! DESCRIPTION:
!> @brief Subroutine to handle memory profiling
!> @details Allows one to track the memory consumption by making use of the memocc()
!> subroutine
!
!> @author
!> Jonathan Chico
!> @date 10.12.2017
!-------------------------------------------------------------------------------
    Module profiling

      Implicit None

    Contains

!control the memory occupation by calculating the overall size in bytes of the allocated arrays
!usage:
! when allocating allocating an array "stuff" of dimension n in the routine "dosome"
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
      Subroutine memocc(istat, isize, array, routine)
        Implicit None
        Character (Len=*), Intent (In) :: array
        Character (Len=*), Intent (In) :: routine
        Integer, Intent (In) :: istat
        Integer, Intent (In) :: isize

! Local variables
        Character (Len=36) :: maxroutine, locroutine
        Character (Len=36) :: maxarray, locarray
        Integer :: nalloc, ndealloc, locpeak, locmemory, iproc
        Integer :: dblsize, mfileno
        Integer *8 :: memory, maxmemory
        Character (Len=1) :: allocationflag

        Save :: memory, nalloc, ndealloc, maxroutine, maxarray, maxmemory
        Save :: locroutine, locarray, locpeak, locmemory, iproc

        mfileno = 77
        dblsize = 1
!print *,'MEMOCC',isize,array
!
        Select Case (array)
        Case ('count')
          If (routine=='start') Then
            memory = 0
            maxmemory = 0
            nalloc = 0
            ndealloc = 0
            locroutine = 'routine'
            locarray = 'array'
            locmemory = 0
            locpeak = 0
            iproc = isize
!open the writing file for the root process
            If (iproc==0) Then
              Open (Unit=mfileno, File='meminfo', Status='unknown')
              Write (mfileno, '(a32,1x,a20,3(1x,a12))') &
                '(Data in kB)             Routine', '    Peak Array', &
                'Routine Mem', 'Total Mem', 'Action'
            End If
          Else If (routine=='stop' .And. iproc==0) Then
!write(mfileno,'(a32,a16,4(1x,i12))')&
!     trim(locroutine),trim(locarray),locmemory/1024,locpeak/1024,memory/1024,&
!     (memory+locpeak-locmemory)/1024
!close(mfileno)
            Write (*, '(1x,a)') '-------------------MEMORY CONSUMPTION &
              &REPORT-----------------------'
            Write (*, '(1x,2(i0,a),i0)') nalloc, ' allocations and ', &
              ndealloc, ' deallocations, remaining memory(B):', memory
            Write (*, '(1x,a,i0,a)') 'Memory occupation peak: ', &
              maxmemory/1024, ' kB'
            Write (*, '(1x,5(a))') 'For the array: "', trim(maxarray), &
              '" in routine "', trim(maxroutine), '"'
            Write (*, '(1x,a)') '-----------------END MEMORY CONSUMPTION &
              &REPORT---------------------'
          End If

        Case Default
!control of the allocation/deallocation status
          If (istat/=0) Then
            If (isize>=0) Then
              Write (*, *) ' subroutine ', routine, &
                ': problem of allocation of array ', array
              Stop
            Else If (isize<0) Then
              Write (*, *) ' subroutine ', routine, &
                ': problem of deallocation of array ', array
              Stop
            End If
          End If
          Select Case (iproc)
          Case (0)
! To be used for inspecting an array which is not deallocated
            If (isize>0) Then
              allocationflag = 'A'
            Else If (isize<0) Then
              allocationflag = 'D'
            Else
              allocationflag = '?'
            End If

            Write (mfileno, '(a32,1x,a20,3(1x,i12),1x,a12)') trim(routine), &
              trim(array), isize*dblsize, memory, maxmemory, allocationflag
            If (trim(locroutine)/=routine) Then
!          write(mfileno,'(a32,a14,4(1x,i12))')&
!               trim(locroutine),trim(locarray),locmemory/1024,locpeak/1024,memory/1024,&
!               (memory+locpeak-locmemory)/1024
              locroutine = routine
              locmemory = isize*dblsize
              locpeak = isize*dblsize
            Else
              locmemory = locmemory + isize*dblsize
              If (locmemory>locpeak) Then
                locpeak = locmemory
              End If
            End If
            locarray = array
            memory = memory + isize*dblsize
            If (memory>maxmemory) Then
              maxmemory = memory
              maxroutine = routine
              maxarray = array
            End If
            If (isize>0) Then
              nalloc = nalloc + 1
            Else If (isize<0) Then
              ndealloc = ndealloc + 1
            End If
          Case Default
            Return
          End Select

        End Select

      End Subroutine


    End Module
