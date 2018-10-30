  subroutine new_directions(magdir0,magdir1)
! Sets up new spin quantization axes
  use global

  implicit none

! initial and final spin quantization axes
  real(kind=r8b), intent(in) :: magdir0(3,nasusc), magdir1(3,nasusc)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: small = 1.d-6, torqmix = 1.d0, dmax = 0.10d0
  integer(kind=i4b) :: ia, nr
  real(kind=r8b)    :: tvec(3), dmagdir(3), dmag, tfac, tlen

  open(file='newdirections.dat',unit=iofile,status='replace')
  do ia=1,nasusc
!   fixed global direction
    if (ispinrot == 0) then
      dmagdir = 0.d0
      magdir(:,ia) = magdir1(:,ia)
!     no constraining field
      bconlen(ia) = 0.d0; bcondir(:,ia) = (/0.d0,0.d0,1.d0/)
!   spin moment too small: don't rotate
!    else if (abs(gs_mlm(1,ia)) < small) then
!      tvec = 0.d0
!      dmagdir = 0.d0
!      magdir(:,ia) = (/0.d0,0.d0,1.d0/)
!   rotate spin moment using gradient of band energy
    else
!     initial assumption about energy scale of torque (1 mRy)
      tfac = 1.d3
!     torque <=> gradient of band energy
      tvec = sum(etorque(:,0:nlmax,ia),dim=2)
      tlen = sqrt(dot_product(tvec,tvec))
      if (tlen*tfac > dmax) tfac = dmax/tlen
      write(iodb,'("new_directions: ia=",i4,"  tfac=",es16.8)') ia, tfac
      dmagdir = torqmix*tfac*tvec
      dmag = sqrt(dot_product(dmagdir,dmagdir))
      if (dmag > dmax) dmagdir = dmagdir*dmax/dmag
!     fixed direction
      if (iarot(ia) == 0) then
        newdir(:,ia) = magdir(:,ia)
!     use both newdir and torque
      else if (iarot(ia) == 1) then
        magdir(:,ia) = magdir(:,ia) + dirmix2*dmagdir(:)
        magdir(:,ia) = magdir(:,ia)/sqrt(dot_product(magdir(:,ia),magdir(:,ia)))
!        magdir(:,ia) = (1.d0 - dirmix3)*magdir(:,ia) + dirmix3*newdir(:,ia)
        magdir(:,ia) = magdir(:,ia) + dirmix3*newdir(:,ia)
        magdir(:,ia) = magdir(:,ia)/sqrt(dot_product(magdir(:,ia),magdir(:,ia)))
!     use only torque
      else if (iarot(ia) == 2) then
        magdir(:,ia) = magdir(:,ia) + dirmix2*dmagdir(:)
        magdir(:,ia) = magdir(:,ia)/sqrt(dot_product(magdir(:,ia),magdir(:,ia)))
!     use only newdir
      else if (iarot(ia) == 3) then
!        magdir(:,ia) = (1.d0 - dirmix3)*magdir(:,ia) + dirmix3*newdir(:,ia)
!        magdir(:,ia) = magdir(:,ia) + dirmix3*newdir(:,ia)
        magdir(:,ia) = magdir(:,ia) + dirmix3*(newdir(:,ia) - magdir(:,ia))
        magdir(:,ia) = magdir(:,ia)/sqrt(dot_product(magdir(:,ia),magdir(:,ia)))
      end if
!     no constraining fields
      if (iarot(ia) /= 0) then
        bconlen(ia) = 0.d0; bcondir(:,ia) = (/0.d0,0.d0,1.d0/)
      end if
    end if
    write(iofile,'("  ia,iarot=",2i4,"  mspin=",f8.5," dmagdir=",3f8.5," newdir=",3f8.5," magdir=",3f8.5)')  &
        ia, iarot(ia), gs_mlm(1,ia), dmagdir, newdir(:,ia), magdir(:,ia)
  end do
  close(iofile)
! spin density multiplied by sign of spin moment from kkrimp
  do ia=1,nasusc
    nr = nrpts(ia)
    new_rho2ns(1:nr,1:lmmax2,2,ia) = new_rho2ns(1:nr,1:lmmax2,2,ia)*magdir0(3,ia)/abs(magdir0(3,ia))
  end do
! write new spin directions to file
  open(file='magdir.dat',unit=iofile,status='replace')
  do ia=1,nasusc
    write(iofile,'(i4,3f12.8,i4)') ia, magdir(:,ia), iarot(ia)
  end do
  if (lbconst) then
    do ia=1,nasusc
      write(iofile,'(i4,4f12.8)') ia, bcondir(:,ia), bconlen(ia)
    end do
  end if
  close(iofile)
! All done!
  end subroutine new_directions
