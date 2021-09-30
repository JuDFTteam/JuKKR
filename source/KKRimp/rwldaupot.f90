module mod_rwldaupot
contains

  !-------------------------------------------------------------------------------
  !> Summary: Reads or writes the LDA+U arrays to the disk.
  !> Author: 
  !> Category: KKRimp, input-output
  !> Deprecated: False ! this needs to be set to True for deprecated subroutines
  !>
  !> lwrite = .true.  : write out
  !> lwrite = .false. : read in
  !>
  !> ipos = 1,2,3,4 : up to which position to read or write, in the order:
  !>
  !> irunldau,lopt,ueff,jeff,erefldau,natyp
  !> ----- 1 -----
  !> wldau
  !> ----- 2 -----
  !> uldau
  !> ----- 3 -----
  !> phi
  !> ----- 4 -----
  !-------------------------------------------------------------------------------
  subroutine rwldaupot(lwrite,ipos,natyp,nspin,lmaxd,irmd,irunldau,ldau)
    use type_ldau, only: ldau_type
    use mod_version_info, only: version_print_header
    implicit none
    ! input:
    integer                        :: ipos,natyp,nspin,lmaxd,irmd
    logical                        :: lwrite
    ! input/output:
    type(ldau_type), allocatable   :: ldau(:)  ! lda+u variables, intended dimension: (natyp) ! lda+u
    integer                        :: irunldau
    ! inside
    character*20                   :: text,text1
    integer                        :: ir,m1,m2,m3,m4,is,i1,mmaxd,lmmaxd,nspind

    mmaxd = 2*lmaxd + 1
    lmmaxd = (lmaxd + 1)**2
    nspind = 2

    open (67,file='ldaupot',form='formatted')
    rewind (67)

    if (lwrite) then

      write (67,*) irunldau, (ldau(i1)%lopt,i1=1,natyp),(ldau(i1)%ueff,i1=1,natyp), (ldau(i1)%jeff,i1=1,natyp),(ldau(i1)%erefldau,i1=1,natyp)   

      if (ipos.eq.1) return

      write (67,9010) 'wldau'
      do i1 = 1,natyp
        if (ldau(i1)%lopt.ge.0) then
          write (67,9020) 'atom ', i1
          write (67,9000) (((ldau(i1)%wldau(m1,m2,is),m1=1,mmaxd),m2=1,mmaxd),is=1,nspin)
        endif
      enddo

      if (ipos.eq.2) return

      write (67,9010) 'uldau'
      do i1 = 1,natyp
        if (ldau(i1)%lopt.ge.0) then
          write (67,9020) 'atom ', i1
          write (67,9000) ((((ldau(i1)%uldau(m1,m2,m3,m4),m1=1,mmaxd),m2=1,mmaxd),m3=1,mmaxd),m4=1,mmaxd)
        endif
      enddo

      if (ipos.eq.3) return

      write (67,9010) 'phi  '
      do i1 = 1,natyp
        if (ldau(i1)%lopt.ge.0) then
          write (67,9020) 'atom ', i1
          write (67,9005) (ldau(i1)%phi(ir),ir=1,irmd)
        endif
      enddo

    else  ! read from file

      read (67,*) irunldau, (ldau(i1)%lopt,i1=1,natyp),(ldau(i1)%ueff,i1=1,natyp), (ldau(i1)%jeff,i1=1,natyp),(ldau(i1)%erefldau,i1=1,natyp) 

      ! allocate everything that is not allocated but should be:
      do i1 = 1,natyp
        if (ldau(i1)%lopt.ge.0) then
          if (.not.allocated(ldau(i1)%wldau)) allocate( ldau(i1)%wldau(mmaxd,mmaxd,nspin) )
          if (.not.allocated(ldau(i1)%uldau)) allocate( ldau(i1)%uldau(mmaxd,mmaxd,mmaxd,mmaxd) )
          if (.not.allocated(ldau(i1)%phi))   allocate( ldau(i1)%phi(irmd) )
          if (.not.allocated(ldau(i1)%cutoff)) allocate( ldau(i1)%cutoff(irmd) )
          if (.not.allocated(ldau(i1)%denmatc)) allocate( ldau(i1)%denmatc(mmaxd,mmaxd,nspin) )
        endif
      enddo

      if (ipos.eq.1) return

      read (67,*) text
      do i1 = 1,natyp
        if (ldau(i1)%lopt.ge.0) then
          read (67,*) text,text1
          read (67,*) (((ldau(i1)%wldau(m1,m2,is),m1=1,mmaxd),m2=1,mmaxd),is=1,nspin)
        endif
      enddo

      if (ipos.eq.2) return

      read (67,*) text
      do i1 = 1,natyp
        if (ldau(i1)%lopt.ge.0) then
          read (67,*) text,text1
          read (67,*) ((((ldau(i1)%uldau(m1,m2,m3,m4),m1=1,mmaxd),m2=1,mmaxd),m3=1,mmaxd),m4=1,mmaxd)
        endif
      enddo

      if (ipos.eq.3) return

      read (67,*) text
      do i1 = 1,natyp
        if (ldau(i1)%lopt.ge.0) then
          read (67,*) text,text1
          read (67,9005) (ldau(i1)%phi(ir),ir=1,irmd)
        endif
      enddo

    endif

    close (67)

    9000 format (4d22.14)
    9005 format (4e20.12)
    9010 format (a5)
    9020 format (a5,i5)

    return
  end subroutine rwldaupot


end module mod_rwldaupot
