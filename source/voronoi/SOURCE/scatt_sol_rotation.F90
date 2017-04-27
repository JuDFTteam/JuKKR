   subroutine scatt_sol_rotation(ie,natom,lmmaxd,theta,phi)

!  This rotates the onsite green function and regular solutions to the global frame

   use mod_rotatespinframe
   use global

   implicit none 

!  ---------------------------------------------------------------------------------------------------------
   integer(kind=i4b), intent(in) :: ie, natom, lmmaxd
   real(kind=r8b),    intent(in) :: theta(natom), phi(natom)
!  ---------------------------------------------------------------------------------------------------------
   complex(kind=c8b)   :: pzl_temp(nasusc,nbmax,2*lmmaxd,2*lmmaxd), pzr_temp(nasusc,nbmax,2*lmmaxd,2*lmmaxd)
   complex(kind=c8b)   :: gf_onsite(nasusc,nbmax,nbmax,2*lmmaxd,2*lmmaxd)
   real(kind=r8b)      :: theta_susc(nasusc), phi_susc(nasusc)
   integer(kind=i4b)   :: i, j, ia, ih, i2(2), i3(3), ib, jb, ni, nj
   integer(kind=i4b)   :: is, js, il, jl,ilm, jlm, ilms, jlms, ilmsn, jlmsn
!  --------------------------------------------------------------------------------------------------------

   ! Initialize all to zero
   pzl_temp  = 0.d0
   pzr_temp  = 0.d0   
   gf_onsite = 0.d0

   ! Loop over all atoms
   do ia = 1, nasusc
     do i = 1, nlmsba(ia) 
       i3 = i2lmsb(:,i,ia) 
       ib = i3(1); ilm = i3(2); is = i3(3)
       il = i2lm(2,ilm)

       ! Onsite Green function 
       do j = 1, nlmsba(ia)
         i3 = i2lmsb(:,j,ia)    
         jb = i3(1); jlm = i3(2); js = i3(3)
         jl = i2lm(2,jlm)
         if (il == jl) then
           ! Back to David's representation
           ni = lms2i_new(ilm,is)
           nj = lms2i_new(jlm,js)
           gf_onsite(ia,ib,jb,ni,nj) = gfpq(i,j,ia,ie) 
         end if  
       end do ! j
 
      ! Regular solutions
       do j = 1, nlms
         i2 = i2lms(:,j)       
         jlm = i2(1); js = i2(2) 
         jl  = i2lm(2,jlm)
         if (il == jl) then
           ! Back to David's representation
           ni = lms2i_new(ilm,is) 
           nj = lms2i_new(jlm,js)
           pzl_temp(ia,ib,ni,nj)= pzl(i,j,ia,ie)
           pzr_temp(ia,ib,ni,nj)= pzr(i,j,ia,ie)   
         end if
       end do ! j

     end do ! i

     ! Angles in kkrflex --> kkrsusc
     ih = iasusc(ia)
     theta_susc(ia) = theta(ih) 
     phi_susc(ia)   = phi(ih)
 
     ! Rotation to global
     do ib = 1, nbmax
       ! Rotate regular solutions
       call rotatematrix(pzr_temp(ia,ib,:,:),theta_susc(ia),phi_susc(ia),lmmaxd,'loc->glob')
       pzl_temp(ia,ib,:,:) = transpose(pzl_temp(ia,ib,:,:))
       call rotatematrix(pzl_temp(ia,ib,:,:),theta_susc(ia),phi_susc(ia),lmmaxd,'loc->glob')
       pzl_temp(ia,ib,:,:) = transpose(pzl_temp(ia,ib,:,:))
       ! Rotate onsite Green function
       do jb = 1, nbmax 
!         if (ib == jb) then
           call rotatematrix(gf_onsite(ia,ib,jb,:,:),theta_susc(ia),phi_susc(ia),lmmaxd,'loc->glob') 
!         end if
       end do ! jb
     end do ! ib
      
     ! Put back the solution to kkrsusc in global
     do i = 1, nlmsba(ia) 
       i3 = i2lmsb(:,i,ia) 
       ib = i3(1); ilm = i3(2); is = i3(3)
       il = i2lm(2,ilm)

       ! Onsite Green function
       do j = 1, nlmsba(ia)
         i3 = i2lmsb(:,j,ia)    
         jb = i3(1); jlm = i3(2); js = i3(3)
         jl = i2lm(2,jlm)
         if (il == jl) then
           ! Back to David's representation
           ni = lms2i_new(ilm,is)
           nj = lms2i_new(jlm,js)
           gfpq(i,j,ia,ie) = gf_onsite(ia,ib,jb,ni,nj) 
         end if  
       end do ! j

       ! Regular solutions
       do j = 1, nlms
         i2 = i2lms(:,j)       
         jlm = i2(1); js = i2(2) 
         jl  = i2lm(2,jlm)
         if (il == jl) then
           ! Back to David's representation
           ni = lms2i_new(ilm,is) 
           nj = lms2i_new(jlm,js)
           pzl(i,j,ia,ie) = pzl_temp(ia,ib,ni,nj) 
           pzr(i,j,ia,ie) = pzr_temp(ia,ib,ni,nj)    
         end if
       end do ! j

     end do ! i
   end do ! ia
   ! All done                 
   end subroutine scatt_sol_rotation
