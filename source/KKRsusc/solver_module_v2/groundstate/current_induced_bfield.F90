  subroutine current_induced_bfield(curr_lm,nr,ia,lmmaxJ,numpan,numrcut)

  use global

! atom number ia
  integer(kind=i4b), intent(in)       :: ia, nr, lmmaxJ
! --> Number of panels > 1
  integer(kind=i4b), intent(in)       :: numpan, numrcut(numpan+1) 
! charge current lm decomposed
  complex(kind=c8b), intent(in)       :: curr_lm(1:3,1:nr,1:lmmaxJ)
! -----------
  integer(kind=i4b)                   :: i,j, ilm, jlm, ilmxyz(1:3), a, b, c
  real(kind=r8b)                      :: tmprgaunt
  complex(kind=c8b)                   :: curr(1:3),curr_r(1:3,1:3), curr_r_r(1:3,1:3,1:3)
  complex(kind=c8b), external         :: radint


  ilmxyz(1)=lm2i(1,1);ilmxyz(2)=lm2i(-1,1);ilmxyz(3)=lm2i(0,1)

  curr(:)=0.d0;curr_r(:,:)=0.d0;curr_r_r(:,:,:)=0.d0
  do a=1,3 !current component
    curr(a) = radint(nr,curr_lm(a,1:nr,1)*r(1:nr)**2,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
    do ilm = 1,lmmaxJ
      do b=1,3 !r component
        tmprgaunt=rgaunt(ilm,ilmxyz(b),1)
	if(abs(tmprgaunt) > ylmtol) then
          curr_r(a,b)=curr_r(a,b)+radint(nr,sqrt(1.d0/3.d0)*tmprgaunt*curr_lm(a,1:nr,ilm)*r(1:nr)**3,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
        end if
        do c=1,3
          tmprgaunt=rgaunt(ilm,ilmxyz(b),ilmxyz(c))
          if(abs(tmprgaunt) > ylmtol) then
            curr_r_r(a,b,c)=curr_r_r(a,b,c)+radint(nr,1.d0/3.d0*tmprgaunt*curr_lm(a,1:nr,ilm)*r(1:nr)**4,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
          end if
        end do
      end do
    end do
  end do

  write(777,'(100f16.8)') (curr(a),a=1,3), ((curr_r(a,b),a=1,3),b=1,3), (((curr_r_r(a,b,c),a=1,3),b=1,3),c=1,3)! IO file for current induced magnetic fields


  end subroutine current_induced_bfield
