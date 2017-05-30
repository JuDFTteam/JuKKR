  subroutine rho_interpolation(numpan,numrcut)
  use global

  implicit none

! --> Number of panels > 1
  integer(kind=i4b), intent(in)       :: numpan, numrcut(numpan+1) 
! ---------------------------------------------------------------------------
  complex(kind=c8b)           :: work(1:nrmax,1:lmmax2,0:3,1:nasusc2) ! ,work2(1:nrmax,1:lmmax2,0:3,1:nasusc2)
  complex(kind=c8b)           :: work3(1:nrmax,1:lmmax2,0:3,1:nasusc2)
  integer(kind=i4b)           :: i4(4),i3(3),ilm,jlm,ib,jb,is,js,ia,ja,ia2,ja2,lmax4,jq,iq,i,j,k,nr,ir
  complex(kind=c8b)           :: block(1:4,1:4),czero=0.d0,cone=(1.d0,0.d0),ci=(0.d0,1.d0)
  real(kind=r8b)              :: r(1:nrmax),x(1:n_int),y(1:n_int),z(1:n_int),xmax,xmin,step,rmax,uvec(1:3),ylmtmp(1:lmmax4),rtmp,rx,ry,rz
  complex(kind=c8b)           :: rho_int(0:3,1:n_int,1:n_int,1:n_int),ctmp
  character(len=1024) :: filename

  
! ----------------------------------------------------- 
! only for ia = 1 
! Generalization easily possible if needed
  ia = 1
  r(:) = rmesh(:,ia)
  nr = nrpts(ia)
  rmax = rmesh(nr,1) 
  !set up the grid
  xmax=rmax  
  xmin=-xmax
  step=(xmax-xmin)/float(n_int-1)
  do i=1,n_int
    x(i)=xmin+step*float(i-1)
  end do
  y(:)=x(:)
  z(:)=x(:)

! calculate lmax for ymy routine
  lmax4 = i2lm(2,lmmax4)

! setting up the cubic spline interpolation
  do is = 0,3
    do ilm = 1,lmmax2
!     calculate second derivate of m_lm, which is used for the spline interpolation
!     work3 is array with second derivative
!     call spline(r,work(1:nr,ilm,is,ia),nr,work2(1,ilm,is,ia),work2(2,ilm,is,ia),work3(1:nr,ilm,is,ia))
      write(*,*) numpan, numrcut(:)
      call spline_panels(r,rho_lm(1:nr,ilm,is,ia),nr,work3(1:nr,ilm,is,ia),numpan, numrcut(:))
    end do
  end do

! interpolation on regular grid
  rho_int=0.d0
  do k=1,n_int
    do j=1,n_int
      do i=1,n_int
        rtmp=sqrt(x(i)**2+y(j)**2+z(k)**2)
        if (rtmp < rmax) then
          uvec(1)=x(i)
          uvec(2)=y(j)
          uvec(3)=z(k)
          call ymy(uvec,lmax4,lmmax4,ylmtmp)
          !write(*,*) i,j,k
            do is = 0,3
              do ilm=1,lmmax2
                call splint(r,rho_lm(:,ilm,is,ia),work3(:,ilm,is,ia),nr,rtmp,ctmp)
                rho_int(is,i,j,k)=rho_int(is,i,j,k)+ctmp*ylmtmp(ilm)
              end do
            end do
        end if
      end do
    end do
  end do


! #############################################################
! write interpolated rho to file
! atom positions
  rx=ri(1,ia)
  ry=ri(2,ia)
  rz=ri(3,ia)
  write(filename,"(A8,I0.3,A4)") "rho_int_",ia,".dat"
  open(unit=123,file=filename)
  write(123,"('# x  y  z  rho charge,x,y,z')")
  write(123,"('#',I5)") n_int
  do k=1,n_int
    do j=1,n_int
      do i=1,n_int
        write(123,"(100e18.9)") x(i)+rx,y(j)+ry,z(k)+rz, (rho_int(is,i,j,k),is=0,3)
      end do
    end do 
  end do
  close(123)




  end subroutine rho_interpolation
