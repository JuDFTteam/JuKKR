  subroutine current_correlation_divergence_interpolation(curr_lm)
! calculates the divergence of the induced current
! lm decompositon and interpolation on regular grid

  use global
  use derivative_panels
  implicit none

  complex(kind=c8b), intent(in)   :: curr_lm(1:nrmax,1:3,1:lmmax4,1:4,1:4,1:nasusc2)
! -------------------------------------------------------------------------------------------------
  complex(kind=c8b)               :: dcurr_lm_dr(1:nrmax,1:3,1:lmmax4,1:4,1:4,1:nasusc2),curr_lm_div(1:nrmax,1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: curr_div_int(1:4,1:4,1:n_int,1:n_int,1:n_int),curr_lm_div_dec(1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: curr_curl_int(1:3,1:4,1:4,1:n_int,1:n_int,1:n_int),curr_lm_curl(1:nrmax,1:3,1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: curr_lm_curl_dec(1:3,1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: work(1:nrmax), work2(1:nrmax,1:lmmax4,1:4,1:4), work3(1:nrmax,1:lmmax4,1:4,1:4), ctmp
  complex(kind=c8b)               :: work4(1:nrmax,1:3),work5(1:2,1:3,1:lmmax4,1:4,1:4),work6(1:nrmax,1:3,1:lmmax4,1:4,1:4)
  complex(kind=c8b)               :: tmp1(1:nrmax), tmp2(1:nrmax)
  integer(kind=i4b)               :: ia2, ia, nr, ir, is, js, a, ilm, jlm, ilmxyz(1:3), il, lm3, m3, lmax4, b,c,d,e,i,j,k,n
  real(kind=r8b)                  :: eps(1:3,1:3,1:3)
  real(kind=r8b)                  :: r(1:nrmax), tmprgaunt, tmprgaunt2, maxelem
  real(kind=r8b)                  :: x(1:n_int), y(1:n_int), z(1:n_int), step, rmax, uvec(1:3), rtmp, ylmtmp(1:lmmax4),rx,ry,rz, dri, dri1,xmin,xmax
  complex(kind=c8b)               :: ci = (0.d0,1.d0)
  character(len=1024)             :: filename
  complex(kind=c8b), external :: radint

 
  write(*,'("divergence calculation started")') 
! -----------------------------------------------------------------------
  write(filename,'("currentspdf.dat")')
  open(file=filename,unit=iofile,status='old',access='append')
! -----------------------------------------------------------------------

! initialize Levi-Civita symbol
  eps(:,:,:)=0.d0
  eps(1,2,3)=1.d0;eps(2,3,1)=1.d0;eps(3,1,2)=1.d0
  eps(3,2,1)=-1.d0;eps(1,3,2)=-1.d0;eps(2,1,3)=-1.d0
  ilmxyz(1)=lm2i(1,1);ilmxyz(2)=lm2i(-1,1);ilmxyz(3)=lm2i(0,1)

  do ia2 = 1,nasusc2
    ia = iasusc2(ia2)
    nr=nrpts(ia)
    r=rmesh(1:nr,ia) !r-points
    do js = 1,4
    do is = 1,4
      do ilm = 1,lmmax4 
        do a = 1,3
          call calc_derivative_panels(curr_lm(1:nr,a,ilm,is,js,ia2),dcurr_lm_dr(1:nr,a,ilm,is,js,ia2),r,nr,npanat(ia),ircutat(:,ia))
!          work(1:nr) = curr_lm(1:nr,a,ilm,is,js,ia2)
!          dcurr_lm_dr(1,a,ilm,is,js,ia2) = (work(2)-work(1))/(r(2)-r(1))
!!         central differences with a second order approach (non uniform r grid)
!          do ir=2,nr-1
!            dri=r(ir+1)-r(ir)
!            dri1=r(ir)-r(ir-1)
!      close(10)
!            dcurr_lm_dr(ir,a,ilm,is,js,ia2)=-dri/dri1/(dri+dri1)*work(ir-1)
!            dcurr_lm_dr(ir,a,ilm,is,js,ia2)=dcurr_lm_dr(ir,a,ilm,is,js,ia2)+(dri-dri1)/dri/dri1*work(ir)
!            dcurr_lm_dr(ir,a,ilm,is,js,ia2)=dcurr_lm_dr(ir,a,ilm,is,js,ia2)+dri1/dri/(dri+dri1)*work(ir+1)
!          end do
!!         backward difference
!          dcurr_lm_dr(nr,a,ilm,is,js,ia2) = (work(nr)-work(nr-1))/(r(nr)-r(nr-1))
        end do
      end do
    end do
    end do
  end do


! calculate the divergence using radial/spherical decomposition of nabla operator
  

  curr_lm_div = 0.d0
  do ia2 = 1,nasusc2
    ia = iasusc2(ia2)
    nr=nrpts(ia)
    r=rmesh(1:nr,ia) !r-points
    do js = 1,4
    do is = 1,4
    do jlm = 1,lmmax4     ! l'm' sum in notes
      work = 0.d0
      do ilm = 1,lmmax2   ! lm sum in notes
        il = i2lm(2,ilm)
        do a = 1,3
!          tmp1 = 0.d0
!          tmp2 = 0.d0
          tmprgaunt = rgaunt(ilm,ilmxyz(a),jlm)
          if(abs(tmprgaunt) > ylmtol) then
            work(1:nr)=work(1:nr) + 1.d0/sqrt(3.d0)*tmprgaunt*dcurr_lm_dr(1:nr,a,ilm,is,js,ia2)
!            tmp1(1:nr) = 1.d0/sqrt(3.d0)*tmprgaunt*dcurr_lm_dr(1:nr,a,ilm,is,js,ia2)
          end if
          do c = 1,3
            do m3 = -il,il !m''' in notes
              lm3 = lm2i(m3,il) ! l m'''
              tmprgaunt2 = rgaunt(lm3,ilmxyz(c),jlm)  
              if(abs(tmprgaunt2) > ylmtol) then 
                do b =1,3
                  work(1:nr) = work(1:nr) + ci/sqrt(3.d0)*tmprgaunt2*eps(a,b,c)*lorb(lm3,ilm,b)*curr_lm(1:nr,a,ilm,is,js,ia2)/r(1:nr)
!                  tmp2(1:nr) = tmp2(1:nr) + ci/sqrt(3.d0)*tmprgaunt2*eps(a,b,c)*lorb(lm3,ilm,b)*curr_lm(1:nr,a,ilm,is,js,ia2)/r(1:nr)
                end do
              end if
            end do
          end do
!          tmp2(1:nr) = tmp1(1:nr) + tmp2(1:nr)
!          do ir = 1,nr
!            if( abs(tmp2(ir)) > 0.01*abs(tmp1(ir)) ) work(ir) = work(ir) + tmp2(ir)
!          end do
        end do      
      end do
      curr_lm_div(1:nr,jlm,is,js,ia2) = work(1:nr)
      curr_lm_div_dec(jlm,is,js,ia2) = radint(nr,r(1:nr)**2*work(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
    end do
    end do
    end do
  end do
! ---------------------------------------------------------------------------
  maxelem = maxval(abs(curr_lm_div_dec))
  where (abs(curr_lm_div_dec) < susctol*maxelem) curr_lm_div_dec = 0.d0
  write(iofile,'("lm decompostion of the divergence")')
  write(iofile,'("ia, im, il")')
  do ia2 = 1,nasusc2
    ia = iasusc2(ia2)
    do jlm = 1,lmmax0  
      do is = 1,4
        write(iofile,'(3i4,10es16.8)') ia, i2lm(:,jlm), (curr_lm_div_dec(jlm,is,js,ia2),js=1,4)
      end do
    end do
  end do
! ---------------------------------------------------------------------------
! calculate the curl using radial/spherical decomposition of nabla operator
  
  curr_lm_curl = 0.d0
  do ia2 = 1,nasusc2
    ia = iasusc2(ia2)
    nr=nrpts(ia)
    r=rmesh(1:nr,ia) !r-points
    do js = 1,4
    do is = 1,4
    do jlm = 1,lmmax4     ! l'm' sum in notes
      work4 = 0.d0
      do a = 1,3
        do ilm = 1,lmmax2   ! lm sum in notes
          il = i2lm(2,ilm)
          do b = 1,3
            tmprgaunt = rgaunt(ilm,ilmxyz(b),jlm)
            do c = 1,3
              if(abs(tmprgaunt) > ylmtol) then
                work4(1:nr,a)=work4(1:nr,a) + eps(a,b,c)/sqrt(3.d0)*tmprgaunt*dcurr_lm_dr(1:nr,c,ilm,is,js,ia2)
              end if
              do e = 1,3
                do m3 = -il,il !m''' in notes
                  lm3 = lm2i(m3,il) ! l m'''
                  tmprgaunt2 = rgaunt(lm3,ilmxyz(e),jlm)  
                  if(abs(tmprgaunt2) > ylmtol) then 
                    do d = 1,3
                      work4(1:nr,a) = work4(1:nr,a) + ci/sqrt(3.d0)*tmprgaunt2*eps(a,b,c)*eps(b,d,e)*lorb(lm3,ilm,d)*curr_lm(1:nr,c,ilm,is,js,ia2)/r(1:nr)
                    end do
                  end if
                end do
              end do
            end do
          end do      
        end do
        curr_lm_curl(1:nr,a,jlm,is,js,ia2) = work4(1:nr,a)
        curr_lm_curl_dec(a,jlm,is,js,ia2) = radint(nr,r(1:nr)**2*work4(1:nr,a),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
      end do
    end do
    end do
    end do
  end do

  maxelem = maxval(abs(curr_lm_curl_dec))
  where (abs(curr_lm_curl_dec) < susctol*maxelem) curr_lm_curl_dec = 0.d0
  write(iofile,'("lm decompostion of the curl")')
  write(iofile,'("ia, im, il, direction")')
  do ia2 = 1,nasusc2
    ia = iasusc2(ia2)
    do jlm = 1,lmmax0  
    do a = 1,3
      do is = 1,4
        write(iofile,'(4i4,10es16.8)') ia, i2lm(:,jlm),a, (curr_lm_curl_dec(a,jlm,is,js,ia2),js=1,4)
      end do
    end do
    end do
  end do

  close(iofile)
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

! setting up the cubic spline interpolation for divergence
!  work2(1,:,:,:)=(curr_lm_div(2,:,:,:,ia2)-curr_lm_div(1,:,:,:,ia2))/(r(2)-r(1))         !derivative at left bound
!  work2(2,:,:,:)=(curr_lm_div(nr,:,:,:,ia2)-curr_lm_div(nr-1,:,:,:,ia2))/(r(nr)-r(nr-1)) !derivative at right bound
  do js = 1,4
    do is = 1,4
      do ilm = 1,lmmax4
!       calculate second derivate of j_lm^div, which is used for the spline interpolation
!       work3 is array with second derivative
!        call spline(r,curr_lm_div(1:nr,ilm,is,js,ia2),nr,work2(1,ilm,is,js),work2(2,ilm,is,js),work3(1:nr,ilm,is,js))
        call spline_panels(r,curr_lm_div(1:nr,ilm,is,js,ia2),nr,work3(1:nr,ilm,is,js),npanat(ia),ircutat(:,ia))
      end do
    end do
  end do


! ---------------------------------------------------------------------------------------------------
! setting up the cubic spline interpolation for curl
!  work5(1,:,:,:,:)=(curr_lm_curl(2,:,:,:,:,ia2)-curr_lm_curl(1,:,:,:,:,ia2))/(r(2)-r(1))         !derivative at left bound
!  work5(2,:,:,:,:)=(curr_lm_curl(nr,:,:,:,:,ia2)-curr_lm_curl(nr-1,:,:,:,:,ia2))/(r(nr)-r(nr-1)) !derivative at right bound
  do js = 1,4
    do is = 1,4
      do ilm = 1,lmmax4
        do n = 1,3
!         calculate second derivate of j_lm, which is used for the spline interpolation
!         work6 is array with second derivative
!          call spline(r,curr_lm_curl(1:nr,n,ilm,is,js,ia2),nr,work5(1,n,ilm,is,js),work5(2,n,ilm,is,js),work6(1:nr,n,ilm,is,js))
          call spline_panels(r,curr_lm_curl(1:nr,n,ilm,is,js,ia2),nr,work6(1:nr,n,ilm,is,js),npanat(ia),ircutat(:,ia))
        end do
      end do
    end do
  end do

! ---------------------------------------------------------------------------------------------------
! interpolation on regular grid
  curr_div_int=0.d0
  curr_curl_int=0.d0
  do k=1,n_int
    do j=1,n_int
      do i=1,n_int
        rtmp=sqrt(x(i)**2+y(j)**2+z(k)**2)
        if (rtmp < rmax) then
          uvec(1)=x(i)
          uvec(2)=y(j)
          uvec(3)=z(k)
          call ymy(uvec,lmax4,lmmax4,ylmtmp)
          do js = 1,4
            do is = 1,4
              do ilm=1,lmmax4
                call splint(r,curr_lm_div(1:nr,ilm,is,js,ia2),work3(1:nr,ilm,is,js),nr,rtmp,ctmp)
                curr_div_int(is,js,i,j,k)=curr_div_int(is,js,i,j,k)+ctmp*ylmtmp(ilm)
                do n = 1,3
                  call splint(r,curr_lm_curl(1:nr,n,ilm,is,js,ia2),work6(1:nr,n,ilm,is,js),nr,rtmp,ctmp)
                  curr_curl_int(n,is,js,i,j,k)=curr_curl_int(n,is,js,i,j,k)+ctmp*ylmtmp(ilm)
                end do
              end do
            end do
          end do
        end if
      end do
    end do
  end do
  write(*,"('current divergence/curl interpolation done')")

! #############################################################
! write interpolated current to file
! atom positions
  rx=ri(1,ia)
  ry=ri(2,ia)
  rz=ri(3,ia)
  do js = 1,4
    do is = 1,4 
      write(filename,"(A18,I0.3,A1,I0.1,A1,I0.1,A4)") "curr_corr_div_int_",ia,"_",is,"_",js,".dat"
      open(unit=10,file=filename)
      write(10,"('# x  y  z  div(curr)')")
      write(10,"('#',I5)") n_int
      do k=1,n_int
        do j=1,n_int
          do i=1,n_int
            write(10,"(100e18.9)") x(i)+rx,y(j)+ry,z(k)+rz, curr_div_int(is,js,i,j,k)
          end do
        end do 
      end do
      close(10)
    end do 
  end do

  do js = 1,4
    do is = 1,4 
      write(filename,"(A19,I0.3,A1,I0.1,A1,I0.1,A4)") "curr_corr_curl_int_",ia,"_",is,"_",js,".dat"
      open(unit=10,file=filename)
      write(10,"('# x  y  z  curl(curr)')")
      write(10,"('#',I5)") n_int
      do k=1,n_int
        do j=1,n_int
          do i=1,n_int
            write(10,"(100e18.9)") x(i)+rx,y(j)+ry,z(k)+rz, curr_curl_int(:,is,js,i,j,k)
          end do
        end do 
      end do
      close(10)
    end do 
  end do

  end do !atom loop
  end subroutine current_correlation_divergence_interpolation
