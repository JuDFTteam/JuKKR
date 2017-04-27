  subroutine spin_directions(magdir0,magdir1)
! Sets up the spin quantization axes
  use global

  implicit none

! initial and final spin quantization axes
  real(kind=r8b), intent(out) :: magdir0(3,nasusc), magdir1(3,nasusc)
! ----------------------------------------------------------------------
  real(kind=r8b), parameter :: invr2 = 1.d0/sqrt(2.d0)
  integer(kind=i4b) :: ia, is, i, ia2
  real(kind=r8b)    :: bfield(3)
  logical           :: exists

! ****************************************
! Defaults
! ****************************************
! New non collinear magnetism from kkrflex
! ****************************************
  if (lsoc_new) then 
    if (ldos) write(*,'("spin_directions: DOS flag => magdir set to kkrflex_angles")')
!   Passing the magnetization   
    do ia = 1, nasusc
      magdir0(:,ia) = magdir(:,ia)  
      magdir1(:,ia) = magdir(:,ia)
    end do      
  else  
    if (ldos) write(*,'("spin_directions: DOS flag => magdir set to +z")')
    do ia=1,nasusc
      magdir0(:,ia) = (/0.d0,0.d0,1.d0/)
      magdir1(:,ia) = (/0.d0,0.d0,1.d0/)
      if (ldos) magdir(:,ia) = (/0.d0,0.d0,1.d0/)
    end do
  end if 
! ****************************************
! Non-defaults
! read_rho2ns sets magdir according to sign of spin moment
  if (lrot) then
!   *********************
!   global spin rotations
    if (ispinrot == 0) then
!   *********************
      do ia=1,nasusc
!        magdir0(:,ia) = (/0.d0,0.d0,1.d0/)
        magdir0(:,ia) = magdir(:,ia)
!        magdir1(:,ia) = urot
        magdir1(:,ia) = urot*magdir(3,ia)/abs(magdir(3,ia))
        magdir(:,ia) = magdir1(:,ia)
!       no constraining fields
        bconlen(ia) = 0.d0
        bcondir(:,ia) = (/0.d0,0.d0,1.d0/)
      end do
!   ********************
!   local spin rotations
    else
!   ********************
      inquire(file='magdir.dat',exist=exists)
      if (exists) then
!     from file
        open(file='magdir.dat',unit=iofile,status='old')
        do ia=1,nasusc
!          magdir0(:,ia) = (/0.d0,0.d0,1.d0/)
          magdir0(:,ia) = magdir(:,ia)
          read(iofile,*) i, magdir1(:,ia), iarot(ia)
          magdir1(:,ia) = magdir1(:,ia)/sqrt(dot_product(magdir1(:,ia),magdir1(:,ia)))
          magdir(:,ia) = magdir1(:,ia)
        end do
!       Constraining fields
        if (lbconst) then
          do ia=1,nasusc
            read(iofile,*) ia2, bcondir(:,ia), bconlen(ia)
          end do
        end if
        close(iofile)
      else
!     from previous call to groundstate
        do ia=1,nasusc
!          magdir0(:,ia) = (/0.d0,0.d0,1.d0/)
          magdir0(:,ia) = magdir(:,ia)
          magdir1(:,ia) = magdir(:,ia)
          iarot(ia) = 3
!         no constraining fields
          bconlen(ia) = 0.d0
          bcondir(:,ia) = (/0.d0,0.d0,1.d0/)
        end do
      end if
!   ******
    end if
!   ******
! ******************************
! Spin orientations from kkrflex
! ******************************
  else 
    do ia = 1, nasusc
      magdir0(:,ia) = magdir(:,ia)
      magdir1(:,ia) = magdir(:,ia)       
    end do ! ia
  end if ! lrot 
! Up and down spinors according to local spin quantization axes
  do ia=1,nasusc
    write(iodb,'("ia=",i4,"  magdir=",3f8.4,"  magdir0=",3f8.4,"  magdir1=",3f8.4)') ia, magdir(:,ia), magdir0(:,ia), magdir1(:,ia)
!   up
    spinproj(2,2,ia) = invr2*cmplx(magdir(1,ia)+magdir(2,ia),0.d0) + magdir(3,ia)
    spinproj(1,2,ia) = invr2*cmplx(magdir(1,ia),magdir(2,ia))
!   dn
    spinproj(2,1,ia) = invr2*cmplx(magdir(1,ia)-magdir(2,ia),0.d0)
    spinproj(1,1,ia) = invr2*cmplx(magdir(1,ia),-magdir(2,ia))     + magdir(3,ia)
  end do
! All done!
  end subroutine spin_directions
