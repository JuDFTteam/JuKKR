  subroutine current_zeeman(curr_lm,nr,lmmaxJ,r,phi,eps,ia)
  
  use global
  
  implicit none
  
  integer(kind=i4b), intent(in) :: nr, lmmaxJ, eps(1:3,1:3,1:3), ia
  complex(kind=c8b), intent(inout) :: curr_lm(1:3,1:nr,1:lmmaxJ)
  real(kind=r8b), intent(in) :: r(1:nr), phi(1:nr,1:nbmax,0:nlmax,1:nsmax)
! ------------------------------------------------------------------------------------
  complex(kind=c8b) :: curr_zeeman_lm(1:3,1:nr,1:lmmaxJ),curr_zeeman_lm_sum(1:3,1:lmmaxJ)
  complex(kind=c8b) :: prefac,tmpfac(1:3)
  integer(kind=i4b) :: ir,i2(2),i3(3),n,ill,ilmxyz(1:3)
  integer(kind=i4b) :: i, ilm, il, im, ib, is, ilm2
  integer(kind=i4b) :: j, jlm, jl, jm, jb, js, jlm2
  integer(kind=i4b) :: k, klm, kl, km, mlm, a, b, c, klm2, mlm2, m2, j2, k2
  complex(kind=c8b) :: cone=(1.d0,0.d0),ci=(0.d0,1.d0)
  real(kind=r8b)    :: tmpgaunt, tmpgaunt2,tmpdphidr,dri,dri1,norm
  integer(kind=i4b) :: tmpeps, tmpeps2
  character(len=1024) :: filename
  complex(kind=c8b), external :: radint
! ------------------------------------------------------------------------------------
  
  
  !prefactor i*e/m*sqrt(4pi/3)
  prefac = ci/sqrt(3.0)
  
  !ilm(x/y/z)
  ilmxyz(1)=lm2i(1,1);ilmxyz(2)=lm2i(-1,1);ilmxyz(3)=lm2i(0,1)
  
  !test of the implementation with a simple input for the magnetization
  !rho_lm(:,:,:,ia)=0.d0
  !do ir = 1,nr
  !  rho_lm(ir,lm2i(0,0),3,ia)=exp(-1.d0/2.d0*r(ir)**2)*sqrt(3.d0)
  !  rho_lm(ir,lm2i(0,1),3,ia)=exp(-1.d0/2.d0*r(ir)**2)
  !end do
  
  
  
  !enjoy the loops...
  curr_zeeman_lm=0.d0
  
  do j=1,3
    do i= 1,3
      do k= 1,3
        tmpeps = eps(k,i,j)
        if(abs(tmpeps) > ylmtol) then
         do ilm= 1,lmmaxJ !l'm' loop
           do jlm = 1, lmmax2 !lm loop
             !longitudinal part first
             tmpgaunt= rgaunt(jlm,ilmxyz(i),ilm)
             if(abs(tmpgaunt) > ylmtol) then
               !forward differences
               curr_zeeman_lm(k,1,ilm)=curr_zeeman_lm(k,1,ilm)-ci*tmpeps*tmpgaunt*(rho_lm(2,jlm,j,ia)-rho_lm(1,jlm,j,ia))/(r(2)-r(1))
               !central differneces
               do ir=2,nr-1
                 dri=r(ir+1)-r(ir)
                 dri1=r(ir)-r(ir-1)
                 curr_zeeman_lm(k,ir,ilm)=curr_zeeman_lm(k,ir,ilm)+ci*tmpeps*tmpgaunt*dri/dri1/(dri+dri1)*(rho_lm(ir-1,jlm,j,ia))
                 curr_zeeman_lm(k,ir,ilm)=curr_zeeman_lm(k,ir,ilm)-ci*tmpeps*tmpgaunt*(dri-dri1)/dri/dri1*(rho_lm(ir,jlm,j,ia))
                 curr_zeeman_lm(k,ir,ilm)=curr_zeeman_lm(k,ir,ilm)-ci*tmpeps*tmpgaunt*dri1/dri/(dri+dri1)*(rho_lm(ir+1,jlm,j,ia))
               end do
               !forward differences
               curr_zeeman_lm(k,nr,ilm)=curr_zeeman_lm(k,nr,ilm)-ci*tmpeps*tmpgaunt*(rho_lm(nr,jlm,j,ia)-rho_lm(nr-1,jlm,j,ia))/(r(nr)-r(nr-1))
             end if
             !second term in the braket (transversal)
             i2 = i2lm(:,jlm)
             jm = i2(1); jl = i2(2)
             do j2 = 1,3
               do k2 = 1,3
                 tmpeps2 = eps(i,j2,k2)
                 if(abs(tmpeps2) > ylmtol) then
                   do km = -jl,jl !sum over m''
                     klm = lm2i(km,jl) !index for lm''
                     tmpgaunt2 = rgaunt(klm,ilmxyz(j2),ilm)
                     if(abs(tmpgaunt2) > ylmtol) then
                       do ir=1,nr
                         curr_zeeman_lm(k,ir,ilm) = curr_zeeman_lm(k,ir,ilm) - tmpeps* rho_lm(ir,jlm,j,ia)/r(ir)*tmpeps2*tmpgaunt2*lorb(klm,jlm,k2)
                       end do
                     end if
                   end do
                 end if
               end do
             end do
           end do
         end do
        end if
      end do
    end do
  end do
  
  curr_zeeman_lm(:,:,:) = prefac * curr_zeeman_lm(:,:,:)
  
  
  !analytic result for the simple input magnetization from above
  !do ir=1,nr
  !  curr_zeeman_lm(1,ir,lm2i(-1,1)) = curr_zeeman_lm(1,ir,lm2i(-1,1)) + exp(-1.d0/2.d0*r(ir)**2)*r(ir)
  !  curr_zeeman_lm(2,ir,lm2i(1,1)) = curr_zeeman_lm(2,ir,lm2i(1,1)) - exp(-1.d0/2.d0*r(ir)**2)*r(ir)
  !  curr_zeeman_lm(1,ir,lm2i(-1,2)) = curr_zeeman_lm(1,ir,lm2i(-1,2)) + exp(-1.d0/2.d0*r(ir)**2)*r(ir)**2/sqrt(5.d0)*(1.d0/r(ir)+1.d0/r(ir)**3)
  !  curr_zeeman_lm(2,ir,lm2i(1,2)) = curr_zeeman_lm(2,ir,lm2i(1,2)) - exp(-1.d0/2.d0*r(ir)**2)*r(ir)**2/sqrt(5.d0)*(1.d0/r(ir)+1.d0/r(ir)**3)
  !end do
  
  
  curr_lm(:,:,:)=curr_lm(:,:,:)+curr_zeeman_lm(:,:,:)
  
  !#############################################################
  !integrate the radial part of j_soc_lm
  curr_zeeman_lm_sum=0.d0
  do n=1,3
    do mlm=1,lmmaxJ
      !integrate r^2*j_lm
      curr_zeeman_lm_sum(n,mlm)=radint(nr,curr_zeeman_lm(n,1:nr,mlm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
    end do
  end do
  
  !#############################################################
  !write integrated current to file
  
  write(filename,"(A22,I0.3,A4)") "current_zeeman_lm_sum_",ia,".dat"
  open(unit=10,file=filename)
  4000 format(2i4,6e18.9)
  
  write(iofile,'("  curr_zeeman_lm_sum=")')
  do mlm=1,lmmaxJ
    write(10,4000) i2lm(:,mlm), curr_zeeman_lm_sum(1,mlm),curr_zeeman_lm_sum(2,mlm),curr_zeeman_lm_sum(3,mlm)
    norm = sqrt(dot_product(real(curr_zeeman_lm_sum(1:3,mlm)),real(curr_zeeman_lm_sum(1:3,mlm))))
    if (norm > atol) then
      write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_zeeman_lm_sum(1:3,mlm))
    end if
  end do
  
  close(10)
  
  end subroutine current_zeeman
