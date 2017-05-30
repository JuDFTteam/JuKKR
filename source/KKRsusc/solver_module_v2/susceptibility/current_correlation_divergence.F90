  subroutine current_correlation_divergence(enhanced_curr_corr,curr_corr_div_lm,iw)
! calculates the divergence of the induced current
! lm decompositon and interpolation on regular grid

  use global
  use mod_derivative_panels
  implicit none

  integer(kind=i4b), intent(in)   :: iw
  complex(kind=c8b), intent(in)   :: enhanced_curr_corr(ngradsum,ndensum)
  complex(kind=c8b),intent(inout) :: curr_corr_div_lm(1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega)
! -------------------------------------------------------------------------------------------------
  complex(kind=c8b)               :: dcurr_lm_dr(1:nrmax,1:3,1:lmmax4,1:4,1:4,1:nasusc2),curr_lm_div(1:nrmax,1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: curr_div_int(1:4,1:4,1:n_int,1:n_int,1:n_int),curr_lm_div_dec(1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: curr_curl_int(1:3,1:4,1:4,1:n_int,1:n_int,1:n_int),curr_lm_curl(1:nrmax,1:3,1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: curr_lm_curl_dec(1:3,1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: work(1:nrmax), work2(1:nrmax)
  complex(kind=c8b)               :: tmp1(1:nrmax), tmp2(1:nrmax)
  integer(kind=i4b)               :: ilmxyz(1:3), il, lm3, m3, lmax4, b,c,d,e,n,b1,ilm1,s1,b2,ilm2,s2
  integer(kind=i4b)               :: jq,iq,i4(4),i3(3),i2(2),i,ib,ilm,is,ia,ia2,j,jb,jlm,js,ja,ja2,q1,q2,a,k,nr,ir
  real(kind=r8b)                  :: eps(1:3,1:3,1:3)
  real(kind=r8b)                  :: r(1:nrmax), tmprgaunt, tmprgaunt2
  real(kind=r8b)                  :: dri, dri1
  complex(kind=c8b)               :: ci = (0.d0,1.d0), czero= (0.d0,0.d0)
  complex(kind=c8b), external     :: radint
  complex(kind=c8b)               :: curr_lm(1:nrmax,1:3,1:lmmax4,1:4,1:4,1:nasusc2)
  complex(kind=c8b)               :: block(1:4,1:4)
  integer(kind=i4b)               :: numpan
  integer(kind=i4b), allocatable  :: numrcut(:) 



  numpan = npanat(ia)
  allocate(numrcut(1:numpan+1))
  numrcut(:) = ircutat(:,ia)

  curr_lm = 0.d0
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
        curr_lm(1:nr,1:3,1:lmmax4,is,js,ia2)=curr_lm(1:nr,1:3,1:lmmax4,is,js,ia2) + gradbasis_lm(1:nr,1:3,1:lmmax4,iq)*enhanced_curr_corr(iq,jq)*suscnorm(jq)
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
                  block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*curr_lm(ir,a,ilm,is,js,ia2)*pc2s(i2is(1,js),i2is(2,js),j)
              end do
              end do
            end do
            end do
            curr_lm(ir,a,ilm,:,:,ia2) = block(:,:)
          end do
        end do
      end do
    end do
  end if

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
          call calc_derivative_panels(curr_lm(1:nr,a,ilm,is,js,ia2),dcurr_lm_dr(1:nr,a,ilm,is,js,ia2),r,nr,numpan,numrcut)
!          dcurr_lm_dr(1,a,ilm,is,js,ia2) = (work(2)-work(1))/(r(2)-r(1))
!!         central differences with a second order approach (non uniform r grid)
!          do ir=2,nr-1
!            dri=r(ir+1)-r(ir)
!            dri1=r(ir)-r(ir-1)
!!      close(10)
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
    do jlm = 1,lmmax2     ! l'm' sum in notes
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
      curr_corr_div_lm(is,js,jlm,1,ia2,1,iw) = curr_lm_div_dec(jlm,is,js,ia2)
    end do
    end do
    end do
  end do

  deallocate(numrcut)
  end subroutine current_correlation_divergence
