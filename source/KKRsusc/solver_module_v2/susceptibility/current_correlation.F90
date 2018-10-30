  subroutine current_correlation(enhanced_susc,zw,iw,curr_corr_lm,curr_corr_0_lm,curr_corr_div_lm,curr_corr_div_0_lm)
! Calculates the spin-current spin correlation function
  use global
  
  implicit none 

  complex(kind=c8b),intent(inout) :: curr_corr_lm(1:3,1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega)
  complex(kind=c8b),intent(inout) :: curr_corr_0_lm(1:3,1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega)
  complex(kind=c8b),intent(inout) :: curr_corr_div_lm(1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega)
  complex(kind=c8b),intent(inout) :: curr_corr_div_0_lm(1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2,1:nomega)
  integer(kind=i4b),intent(in) :: iw
  complex(kind=c8b),intent(in) :: zw(nomega)
  complex(kind=c8b),intent(in) :: enhanced_susc(ndensum,ndensum)
! --------------------------------------------------------------------
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
! work arrays
  complex(kind=c8b) :: work(ngradsum,ndensum)
  complex(kind=c8b) :: work2(ngradsum,ndensum)
  complex(kind=c8b) :: block(1:4,1:4), block2(1:4,1:4)
  complex(kind=c8b) :: enhanced_curr_corr(ngradsum,ndensum)
  real(kind=r8b)    :: start, finish, tmprgaunt
  integer(kind=i4b) :: i4(1:4),jb,jlm,js,ja2,ib,ilm,is,ia2,a,q1,q2,s1,s2,iq,jq,ia,i3(1:3),iq1,iq0,i,j,klm,ja
  complex(kind=c8b), external         :: radint
  real(kind=r8b)    :: r(1:nrmax)
  integer(kind=i4b) :: nr
  complex(kind=c8b) :: mass_correction_r(1:nrmax,1:4,1:4,1:nasusc2,1:nasusc2),mass_correction(1:4,1:4)
  
  kscurrcorr=0.d0

  call cpu_time(start)
! KS correlation function from taylor expansion
  call zcopy(ngradsum*ndensum,kscurrcorr0,1,kscurrcorr,1)
  call zaxpy(ngradsum*ndensum,zw(iw),kscurrcorr1,1,kscurrcorr,1)
  call zaxpy(ngradsum*ndensum,0.5d0*zw(iw)**2,kscurrcorr2,1,kscurrcorr,1)
  call cpu_time(finish)
  write(*,'("Vector multiplication time=",es16.8)') finish-start

  call cpu_time(start)
! Matrix multiplication of chi_0 * xc kernel 
  call zgemm('N','N',ngradsum,ndensum,ndensum,cone,kscurrcorr,ngradsum,kernel,ndensum,czero,work,ngradsum)
  call cpu_time(finish)
  write(*,'("chi_0 kernel multiplication time=",es16.8)') finish-start

  call cpu_time(start)
! Matrix multiplication of current_correlation * xc kernel * chi 
  call zgemm('N','N',ngradsum,ndensum,ndensum,cone,work,ngradsum,enhanced_susc,ndensum,czero,work2,ngradsum)
  call cpu_time(finish)
  write(*,'("chi_0 kernel chi multiplication time=",es16.8)') finish-start

  enhanced_curr_corr(:,:)= kscurrcorr(:,:) + work2(:,:)
 
! output interpolated current distribution on regular grid 
!  if(lcurrcorrint) call current_correlation_interpolation(kscurrcorr,iw)
  if(lcurrcorrint) call current_correlation_interpolation(enhanced_curr_corr,iw)
  call current_correlation_divergence(enhanced_curr_corr,curr_corr_div_lm,iw)
  call current_correlation_divergence(kscurrcorr,curr_corr_div_0_lm,iw)

  curr_corr_lm(:,:,:,:,:,:,:,iw)=0.d0
  curr_corr_0_lm(:,:,:,:,:,:,:,iw)=0.d0

  do jq = 1, ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
    if(jlm < lmmax2+1) then
    do iq = 1, ngradsum
      i3 = i2almsbgrad(:,iq)
!    do ia2 = 1, nasusc2
!      ia = iasusc2(ia2)
!      iq0 = sum(nalmsbgf(1:ia-1))
!      iq1 = sum(nalmsbgf(1:ia))
!      do i=1,iq1-iq0 
!        iq = iq0 + i
!        i3 = i2almsbgf(:,iq)     
      q1 = i3(1); q2 = i3(2); ia2 = i3(3); ia = iasusc2(ia2)
      i3 = i2lmsb(:,q1,ia)
      s1 = i3(3)  !check spin indices once again!!!
      i3 = i2lmsb(:,q2,ia)
      s2 = i3(3)  !check spin indices once again!!!
      is = is2i(s1,s2)
      do ilm=1,lmmax2
        do a = 1,3
          curr_corr_lm(a,is,js,ilm,jlm,ia2,ja2,iw)=curr_corr_lm(a,is,js,ilm,jlm,ia2,ja2,iw)+gradnorm(a,ilm,iq)*enhanced_curr_corr(iq,jq)*suscnorm(jq)
          curr_corr_0_lm(a,is,js,ilm,jlm,ia2,ja2,iw)=curr_corr_0_lm(a,is,js,ilm,jlm,ia2,ja2,iw)+gradnorm(a,ilm,iq)*kscurrcorr(iq,jq)*suscnorm(jq)
        end do
      end do
    end do    
    end if
  end do

!  write(*,'("curr corr lm calculated")')
  if (lcartesian) then
    do ja2=1,nasusc2
    do ia2=1,nasusc2
      do jlm=1,lmmax2
      do ilm=1,lmmax2
        do a = 1,3
          block(:,:) = czero
          block2(:,:) = czero
          do j=1,4
          do i=1,4
            do js=1,4
            do is=1,4
                block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*curr_corr_lm(a,is,js,ilm,jlm,ia2,ja2,iw)*pc2s(i2is(1,js),i2is(2,js),j)
                block2(i,j) = block2(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*curr_corr_0_lm(a,is,js,ilm,jlm,ia2,ja2,iw)*pc2s(i2is(1,js),i2is(2,js),j)
            end do
            end do
          end do
          end do
          curr_corr_lm(a,:,:,ilm,jlm,ia2,ja2,iw) = block
          curr_corr_0_lm(a,:,:,ilm,jlm,ia2,ja2,iw) = block2
        end do
      end do
      end do
    end do
    end do
  end if



! Scalar relativistic mass correction to the continuity equation
! only the sphere average ( ilm=jlm=0 ) is calculated
! chi_jm * (d/dk ln(M(r)))
  if(iw == 1) then
  if(lscalarcorr) then
    mass_correction_r(:,:,:,:,:)=0.d0
    do jq = 1, ndensum
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
      if(jlm == 1) then
      do iq = 1, ngradsum
        i3 = i2almsbgrad(:,iq)
!      do ia2 = 1, nasusc2
!        ia = iasusc2(ia2)
!        iq0 = sum(nalmsbgf(1:ia-1))
!        iq1 = sum(nalmsbgf(1:ia))
!        do i=1,iq1-iq0 
!          iq = iq0 + i
!          i3 = i2almsbgf(:,iq)     
        q1 = i3(1); q2 = i3(2); ia2 = i3(3); ia = iasusc2(ia2)
        i3 = i2lmsb(:,q1,ia)
        s1 = i3(3)  !check spin indices once again!!!
        i3 = i2lmsb(:,q2,ia)
        s2 = i3(3)  !check spin indices once again!!!
        is = is2i(s1,s2)
        nr = nrpts(ia)
        r(1:nr) = rmesh(1:nr,ia)
        do ilm=1,lmmax
          do a = 1,3
!           scalar relativistic effects to the continuity equation (d/dk ln(M(r)))
            do klm=1,lmmax
              tmprgaunt = rgaunt(ilm,klm,1)
              if(abs(tmprgaunt) > ylmtol) then 
                mass_correction_r(1:nr,is,js,ia2,ja2) = mass_correction_r(1:nr,is,js,ia2,ja2) - gradbasis_lm(1:nr,a,ilm,iq)*enhanced_curr_corr(iq,jq)*suscnorm(jq)*tmprgaunt*grad_mass(a,1:nr,klm,ia,nesusc)
              end if
            end do
          end do
        end do
      end do    
      end if
    end do

    ia = 1 
    do js=1,4
      do is=1,4
        mass_correction(is,js) = radint(nr,mass_correction_r(1:nr,is,js,ia,ia)*r(1:nr)**2,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
      end do
    end do

    if (lcartesian) then
      block(:,:) = czero
      do j=1,4
      do i=1,4
        do js=1,4
        do is=1,4
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*mass_correction(is,js)*pc2s(i2is(1,js),i2is(2,js),j)
        end do
        end do
      end do
      end do
      mass_correction(:,:) = block
    end if

!   output results
    write(*,'(/,"Mass correction")')
    ia = 1
    ja = 1
    write(iofile,'("Atom i=",i4," j=",i4)') ia,ja
    do i=1,4
      write(*,'(6i4,8es16.8)') ia, ja, i2lm(:,1), i2lm(:,1), (mass_correction(i,j),j=1,4)
    end do
  end if
  end if

  end subroutine current_correlation
