  subroutine susc_interpolation(enhanced_susc,iw)
! interpolates the susceptibility
! also calclates m x dm/dt which is used for the spin pumped current in the equation of Tserkovnyak
  use global

  implicit none
 
  integer(kind=c8b), intent(in)  :: iw
  complex(kind=c8b), intent(in)  :: enhanced_susc(ndensum,ndensum)
! ---------------------------------------------------------------------------
  complex(kind=c8b)           :: work(1:nrmax,1:lmmax4,1:4,1:4,1:nasusc2) ! ,work2(1:2,1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)           :: work3(1:nrmax,1:lmmax4,1:4,1:4,1:nasusc2),m_cross_susc(1:nrmax,1:lmmax4,1:3,1:nasusc2)
  integer(kind=i4b)           :: i4(4),i3(3),ilm,jlm,ib,jb,is,js,ia,ja,ia2,ja2,lmax4,jq,iq,i,j,k,nr,ir,klm,a,b,c
  complex(kind=c8b)           :: block(1:4,1:4),czero=0.d0,cone=(1.d0,0.d0),ci=(0.d0,1.d0)
  real(kind=r8b)              :: r(1:nrmax),x(1:n_int),y(1:n_int),z(1:n_int),xmax,xmin,step,rmax,uvec(1:3),ylmtmp(1:lmmax4),rtmp,rx,ry,rz,tmprgaunt
  real(kind=r8b)              :: eps(1:3,1:3,1:3)
  complex(kind=c8b)           :: susc_int(1:4,1:4,1:n_int,1:n_int,1:n_int),ctmp
  character(len=1024) :: filename

  if(iw/=1) return


  eps(:,:,:)=0.d0
  eps(1,2,3)=1.d0;eps(2,3,1)=1.d0;eps(3,1,2)=1.d0
  eps(3,2,1)=-1.d0;eps(1,3,2)=-1.d0;eps(2,1,3)=-1.d0


  work = 0.d0
  do jq = 1, ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4); ja = iasusc2(ja2)
    if(jlm == 1 .AND. ja == 1) then !uniform magnetic field
      do iq = 1, ndensum
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4); ia = iasusc2(ja2)
        nr = nrpts(ia)
        work(1:nr,ilm,is,js,ia2)=work(1:nr,ilm,is,js,ia2) + suscbasis(1:nr,ib,ilm,is,ia2)*enhanced_susc(iq,jq)*suscnorm(jq)
      end do
    end if
  end do

  if (lcartesian) then
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      nr = nrpts(ia)
      do ilm=1,lmmax4
        do ir = 1,nr
          block(:,:) = czero
          do j=1,4
          do i=1,4
            do js=1,4
            do is=1,4
                block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*work(ir,ilm,is,js,ia2)*pc2s(i2is(1,js),i2is(2,js),j)
            end do
            end do
          end do
          end do
          work(ir,ilm,:,:,ia2) = block(:,:)
        end do
      end do
    end do
    
! Tserkovnyak spin current
! first part: m x dm/dt = -i omega (m x susc)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      nr = nrpts(ia)
      do ilm=1,lmmax2
      do jlm=1,lmmax2
      do klm=1,lmmax2
        tmprgaunt=rgaunt(klm,jlm,ilm)
        if(abs(tmprgaunt)>ylmtol) then
          do a = 1,3
          do b = 1,3
          do c = 1,3
            m_cross_susc(1:nr,ilm,a,ia2)=eps(a,b,c)*tmprgaunt*rho_lm(1:nr,jlm,b,ia)*work(1:nr,klm,c,1,ia2)
          end do
          end do
          end do
        end if
      end do
      end do
      end do
    end do


  end if

! ----------------------------------------------------- only for ia = 1
  ia = 1
  r(:) = rmesh(:,ia)
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
! work2(1,:,:,:,ia)=(work(2,:,:,:,ia)-work(1,:,:,:,ia))/(r(2)-r(1))         !derivative at left bound
! work2(2,:,:,:,ia)=(work(nr,:,:,:,ia)-work(nr-1,:,:,:,ia))/(r(nr)-r(nr-1)) !derivative at right bound
  do js = 1,4
    do is = 1,4
      do ilm = 1,lmmax4
!       calculate second derivate of m_lm, which is used for the spline interpolation
!       work3 is array with second derivative
!        call spline(r,work(1:nr,ilm,is,js,ia),nr,work2(1,ilm,is,js,ia),work2(2,ilm,is,js,ia),work3(1:nr,ilm,is,js,ia))
        call spline_panels(r,work(1:nr,ilm,is,js,ia),nr,work3(1:nr,ilm,is,js,ia),npanat(ia),ircutat(:,ia))
      end do
    end do
  end do

! interpolation on regular grid
  susc_int=0.d0
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
          do js = 1,4
            do is = 1,4
              do ilm=1,lmmax4
                call splint(r,work(:,ilm,is,js,ia),work3(:,ilm,is,js,ia),nr,rtmp,ctmp)
                susc_int(is,js,i,j,k)=susc_int(is,js,i,j,k)+ctmp*ylmtmp(ilm)
              end do
            end do
          end do
        end if
      end do
    end do
  end do


! #############################################################
! write interpolated susc to file
! atom positions
  rx=ri(1,ia)
  ry=ri(2,ia)
  rz=ri(3,ia)
  do js = 1,4
    do is = 1,4 
      write(filename,"(A9,I0.3,A1,I0.1,A1,I0.1,A4)") "susc_int_",ia,"_",is,"_",js,".dat"
      open(unit=11,file=filename)
      write(11,"('# x  y  z  susc')")
      write(11,"('#',I5)") n_int
      do k=1,n_int
        do j=1,n_int
          do i=1,n_int
            write(11,"(100e18.9)") x(i)+rx,y(j)+ry,z(k)+rz, susc_int(is,js,i,j,k)
          end do
        end do 
      end do
      close(11)
    end do 
  end do
  end subroutine susc_interpolation
