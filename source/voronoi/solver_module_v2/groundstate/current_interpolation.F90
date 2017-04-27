  subroutine current_interpolation(curr_lm,lmmaxJ,r,nr,curr_int,ia,spin,numpan,numrcut)
  use global, only: i4b,r8b,c8b,i2lm,ri,lpositions,n_int
  implicit none
  
  integer(kind=i4b), intent(in)  :: nr,lmmaxJ,ia,spin
  real(kind=r8b),    intent(in)  :: r(nr)
  complex(kind=c8b), intent(in)  :: curr_lm(1:3,1:nr,1:lmmaxJ)
  complex(kind=c8b), intent(out) :: curr_int(1:3,1:n_int,1:n_int,1:n_int)
! --> Number of panels > 1
  integer(kind=i4b), intent(in) :: numpan, numrcut(numpan+1) 
! ###################################################
  
  complex(kind=c8b)           :: dcurr_lm_dr2(1:3,1:nr,1:lmmaxJ),curr_lm_tmp
! complex(kind=c8b)           :: dcurr_lm_dr(1:3,1:2,1:lmmaxJ)
  real(kind=r8b)              :: x(n_int),y(n_int),z(n_int),xmin,xmax,step,uvec(1:3)
  integer(kind=i4b)           :: i,j,k,ilm,jlm,i2(2),lmax,n
  real(kind=r8b)              :: rtmp,ylmtmp(1:lmmaxJ),rx,ry,rz
  character(len=1024) :: filename
  
  
  
! calculate lmax
  i2 = i2lm(:,lmmaxJ)
  lmax =i2(2)
! lmax=sqrt(float(lmmaxJ))-1
  
! set up the grid
  xmax=r(nr)
  xmin=-xmax
  step=(xmax-xmin)/float(n_int-1)
  do i=1,n_int
    x(i)=xmin+step*float(i-1)
  end do
  y(:)=x(:)
  z(:)=x(:)
  
! setting up the cubic spline interpolation
!  dcurr_lm_dr(:,1,:)=(curr_lm(:,2,:)-curr_lm(:,1,:))/(r(2)-r(1))         !derivative at left bound
 ! dcurr_lm_dr(:,2,:)=(curr_lm(:,nr,:)-curr_lm(:,nr-1,:))/(r(nr)-r(nr-1)) !derivative at right bound
  do ilm=1,lmmaxJ
    do n=1,3
!     calculate second derivate of j_lm, which is used for the spline interpolation
!      call spline(r,curr_lm(n,:,ilm),nr,dcurr_lm_dr(n,1,ilm),dcurr_lm_dr(n,2,ilm),dcurr_lm_dr2(n,:,ilm))
      call spline_panels(r,curr_lm(n,:,ilm),nr,dcurr_lm_dr2(n,:,ilm),numpan,numrcut)
    end do
  end do
  
  curr_int=0.d0
  do k=1,n_int
    do j=1,n_int
      do i=1,n_int
        rtmp=sqrt(x(i)**2+y(j)**2+z(k)**2)
        if (rtmp < r(nr)) then
          uvec(1)=x(i)
          uvec(2)=y(j)
          uvec(3)=z(k)
          call ymy(uvec,lmax,lmmaxJ,ylmtmp)
          !write(*,*) i,j,k
          do ilm=1,lmmaxJ
            do n=1,3
              !interpolate j_lm(r) at position rtmp
              !consider only components unequal (-2,2) (2,2)!!!!!!!!!
  !            if (ilm /= 2 .and. ilm/=4) then 
                call splint(r,curr_lm(n,:,ilm),dcurr_lm_dr2(n,:,ilm),nr,rtmp,curr_lm_tmp)
  !              write(*,*) n,i,j,k 
                curr_int(n,i,j,k)=curr_int(n,i,j,k)+curr_lm_tmp*ylmtmp(ilm)
  !            end if
            end do
          end do
        end if
      end do
    end do
  end do
  
  
! #############################################################
! write interpolated current to file ##########################
  
  write(filename,"(A12,I0.3,A1,I0.1,A4)") "current_int_",ia,"_",spin,".dat"
  open(unit=10,file=filename)
  5000 format(9e18.9)
  write(10,"('#',I5)") n_int
  !atom positions
  rx=ri(1,ia)
  ry=ri(2,ia)
  rz=ri(3,ia)
  do k=1,n_int
    do j=1,n_int
      do i=1,n_int
        write(10,5000) x(i)+rx,y(j)+ry,z(k)+rz, curr_int(:,i,j,k)
      end do
    end do 
  end do
  
  close(10)
  
  if(lpositions/=.true.) then 
    write(*,*) "####################"
    write(*,*) "lpositions= F"
    write(*,*) "current interpolation failed!"
  end if
  
  if(spin==0) then 
    do k=1,n_int
      do j=1,n_int
        do i=1,n_int
          write(220,5000) x(i)+rx,y(j)+ry,z(k)+rz, curr_int(:,i,j,k)
        end do
      end do 
    end do
  end if
  
  if(spin==1) then 
    do k=1,n_int
      do j=1,n_int
        do i=1,n_int
          write(221,5000) x(i)+rx,y(j)+ry,z(k)+rz, curr_int(:,i,j,k)
        end do
      end do 
    end do
  end if
  
  if(spin==2) then 
    do k=1,n_int
      do j=1,n_int
        do i=1,n_int
          write(222,5000) x(i)+rx,y(j)+ry,z(k)+rz, curr_int(:,i,j,k)
        end do
      end do 
    end do
  end if
  
  if(spin==3) then 
    do k=1,n_int
      do j=1,n_int
        do i=1,n_int
          write(223,5000) x(i)+rx,y(j)+ry,z(k)+rz, curr_int(:,i,j,k)
        end do
      end do 
    end do
  end if
  
  
  
  end subroutine current_interpolation
