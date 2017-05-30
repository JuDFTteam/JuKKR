  subroutine current_correlation_interpolation(enhanced_curr_corr,iw)
! Interpolation routine for the spin current at a specific frequency (iw)
! Uniform magnetic field acting on the first atom is assumed => (0,0) component of r' dependence
! routine is only working for ia ==1
  use global
  use mod_derivative_panels
  use mod_spline_panels
  implicit none 
 
  integer(kind=i4b), intent(in)  :: iw
  complex(kind=c8b), intent(in)  :: enhanced_curr_corr(ngradsum,ndensum)
! ----------------------------------------------------------------------
  complex(kind=c8b)     :: work(1:nrmax,1:3,1:lmmax4,1:4,1:4,1:nasusc2),work2(1:2,1:3,1:lmmax4,1:4,1:4),work3(1:nrmax,1:3,1:lmmax4,1:4,1:4) 
  integer(kind=i4b)     :: jq,iq,i4(4),i3(3),i2(2),i,ib,ilm,is,ia,ia2,j,jb,jlm,js,ja,ja2,q1,q2,a,k
  integer(kind=i4b)     :: ilm1,ilm2,b1,b2,s1,s2,nr,lmax4,n,ir
  real(kind=r8b)        :: rmax,xmax,xmin,step,rtmp,uvec(1:3),ylmtmp(1:lmmax4),rx,ry,rz,r(1:nrmax)
  complex(kind=c8b)     :: curr_int(1:3,1:4,1:4,1:n_int,1:n_int,1:n_int),ctmp,block(1:4,1:4),czero=0.d0
  real(kind=r8b)        :: x(1:n_int),y(1:n_int),z(1:n_int)
  character(len=1024) :: filename

  write(*,"('current interpolation started')")
  if(iw /=1) return

  work = 0.d0
  do jq = 1, ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4); ja = iasusc2(ja2)
    if(jlm == 1 .AND. ja == 1) then !uniform magnetic field
      do iq = 1, ngradsum
        i3 = i2almsbgrad(:,iq)
        q1 = i3(1); q2 = i3(2); ia2 = i3(3); ia = iasusc2(ia2); nr = nrpts(ia)
        i3 = i2lmsb(:,q1,ia)
        b1 = i3(1); ilm1 = i3(2); s1 = i3(3)
        i3 = i2lmsb(:,q2,ia)
        b2 = i3(1); ilm2 = i3(2); s2 = i3(3)
        is = is2i(s1,s2)
        work(1:nr,1:3,1:lmmax4,is,js,ia2)=work(1:nr,1:3,1:lmmax4,is,js,ia2) + gradbasis_lm(1:nr,1:3,1:lmmax4,iq)*enhanced_curr_corr(iq,jq)*suscnorm(jq)
      end do
    end if
  end do

  if (lcartesian) then
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      nr = nrpts(ia)
      do ilm=1,lmmax4
        do a = 1,3
          do ir = 1,nr
            block(:,:) = czero
            do j=1,4
            do i=1,4
              do js=1,4
              do is=1,4
                  block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*work(ir,a,ilm,is,js,ia2)*pc2s(i2is(1,js),i2is(2,js),j)
              end do
              end do
            end do
            end do
            work(ir,a,ilm,:,:,ia2) = block(:,:)
          end do
        end do
      end do
    end do
  end if
! work array is j_alpha_lm^ab of eq. (34) in notes

  if(lcurrcorrdiv) call current_correlation_divergence_interpolation(work)

! #############################################################
! from here only for ia==1
! #############################################################

  do ia2 = 1,nasusc2
  ia = iasusc2(ia2)
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
! work2(1,:,:,:,:)=(work(2,:,:,:,:,ia2)-work(1,:,:,:,:,ia2))/(r(2)-r(1))         !derivative at left bound
! work2(2,:,:,:,:)=(work(nr,:,:,:,:,ia2)-work(nr-1,:,:,:,:,ia2))/(r(nr)-r(nr-1)) !derivative at right bound
  do js = 1,4
    do is = 1,4
      do ilm = 1,lmmax4
        do n = 1,3
!         calculate second derivate of j_lm, which is used for the spline interpolation
!         work3 is array with second derivative
          call spline_panels(r,work(1:nr,n,ilm,is,js,ia2),nr,work3(1:nr,n,ilm,is,js),npanat(ia),ircutat(:,ia))
        end do
      end do
    end do
  end do

! interpolation on regular grid
  curr_int=0.d0
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
                do n=1,3
                    call splint(r,work(:,n,ilm,is,js,ia2),work3(:,n,ilm,is,js),nr,rtmp,ctmp)
                    curr_int(n,is,js,i,j,k)=curr_int(n,is,js,i,j,k)+ctmp*ylmtmp(ilm)
                end do
              end do
            end do
          end do
        end if
      end do
    end do
  end do
  write(*,"('current interpolation done')")


! #############################################################
! write interpolated current to file
! atom positions
  rx=ri(1,ia)
  ry=ri(2,ia)
  rz=ri(3,ia)
  do js = 1,4
    do is = 1,4 
      write(filename,"(A14,I0.3,A1,I0.1,A1,I0.1,A4)") "curr_corr_int_",ia,"_",is,"_",js,".dat"
      open(unit=10,file=filename)
      write(10,"('# x  y  z  curr_x  curr_y  curr_z')")
      write(10,"('#',I5)") n_int
      do k=1,n_int
        do j=1,n_int
          do i=1,n_int
            write(10,"(100e18.9)") x(i)+rx,y(j)+ry,z(k)+rz, curr_int(:,is,js,i,j,k)
          end do
        end do 
      end do
      close(10)
    end do 
  end do

  end do !atom loop

  end subroutine current_correlation_interpolation
  
