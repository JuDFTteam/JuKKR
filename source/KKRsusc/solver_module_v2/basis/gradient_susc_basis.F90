  subroutine gradient_susc_basis()
! Analysis of the gradient basis used in the spin-current spin correlation function
! gradnorm is the used array with a basis grad(phi*Y)*phi*Y
! gradnorm2 is a testarray for the basis  grad(phi*Y)*phi*Y - phi*Y*grad(phi*Y)
  use global
  use mod_derivative_panels

  implicit none


  complex(kind=c8b) :: work(1:nrmax,1:3)
  integer(kind=i4b) :: ia,ia2,iq,i4(4),i3(3),i2(2),i,q1,q2,b1,b2,lm1,lm2,lm3,s1,s2,l1,l2,m1,m2,m3,nr,ir
  integer(kind=i4b) :: a,b,c, ilm,jlm
  real(kind=r8b)    :: dri,dri1
! Levi-Civita symbol
  real(kind=r8b) :: eps(1:3,1:3,1:3)
  complex(kind=c8b), external :: radint
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0), ci = (0.d0,1.d0)
  real(kind=r8b)    :: dphidr(1:nrmax,1:nbmax,0:nlmax,1:nsmax,1:nasusc),r(1:nrmax),phi(1:nrmax)
  integer(kind=i4b) :: ilmxyz(1:3)
  real(kind=r8b)    :: tmprgaunt, tmprgaunt2, tmprgaunt3


  if(allocated(gradnorm)) deallocate(gradnorm)
  allocate(gradnorm(1:3,1:lmmax4,1:ngfsum))

  if(allocated(gradbasis_lm)) deallocate(gradbasis_lm)
  allocate(gradbasis_lm(1:nrmax,1:3,1:lmmax4,1:ngradsum)) !gradbasis(radial mesh,direction,lm decomp, l1 m1 b1, l2 m2 b2)

!  if(allocated(gfnorm)) deallocate(gfnorm)
!  allocate(gfnorm(1:lmmax4,1:ngfsum))

! initialize Levi-Civita symbol
  eps(:,:,:)=0.d0
  eps(1,2,3)=1.d0;eps(2,3,1)=1.d0;eps(3,1,2)=1.d0
  eps(3,2,1)=-1.d0;eps(1,3,2)=-1.d0;eps(2,1,3)=-1.d0
  ilmxyz(1)=lm2i(1,1);ilmxyz(2)=lm2i(-1,1);ilmxyz(3)=lm2i(0,1)

  do ia2 = 1,nasusc2
    ia = iasusc2(ia2)
    nr=nrpts(ia)
    r=rmesh(1:nr,ia) !r-points
    do i = 1,nlmsba(ia)
      i3 = i2lmsb(:,i,ia)
      b1 = i3(1); lm1 = i3(2); s1 = i3(3)
      i2 = i2lm(:,lm1)
      m1 = i2(1); l1 = i2(2)
      phi(1:nr) = phiref(1:nr,b1,l1,s1,ia)/r(1:nr)
      call calc_derivative_panels(phi,dphidr(1:nr,b1,l1,s1,ia),r,nr,npanat(ia),ircutat(:,ia))
!      dphidr(1,b1,l1,s1,ia) = (phi(2)-phi(1))/(r(2)-r(1))
!!     central differences with a second order approach (non uniform r grid)
!      do ir=2,nr-1
!        dri=r(ir+1)-r(ir)
!        dri1=r(ir)-r(ir-1)
!        dphidr(ir,b1,l1,s1,ia)=-dri/dri1/(dri+dri1)*phi(ir-1)
!        dphidr(ir,b1,l1,s1,ia)=dphidr(ir,b1,l1,s1,ia)+(dri-dri1)/dri/dri1*phi(ir)
!        dphidr(ir,b1,l1,s1,ia)=dphidr(ir,b1,l1,s1,ia)+dri1/dri/(dri+dri1)*phi(ir+1)
!      end do
!!     backward difference
!      dphidr(nr,b1,l1,s1,ia) = (phi(nr)-phi(nr-1))/(r(nr)-r(nr-1))
!      dphidr(1:nr,b1,l1,s1,ia) = dphidr(1:nr,b1,l1,s1,ia)*r(1:nr)
    end do
  end do

!  write(*,'("gradnorm calculated")')
! In the following an array to calculate the spherical average is constructed
  do iq=1,ngradsum
    i3 = i2almsbgrad(:,iq)     
    q1 = i3(1); q2 = i3(2); ia2 = i3(3); ia = iasusc2(ia2)
    i3 = i2lmsb(:,q1,ia)
    b1 = i3(1); lm1 = i3(2); s1 = i3(3)
    i2 = i2lm(:,lm1)
    m1 = i2(1); l1 = i2(2)
    i3 = i2lmsb(:,q2,ia)
    b2 = i3(1); lm2 = i3(2); s2 = i3(3)
    i2 = i2lm(:,lm2)
    m2 = i2(1); l2 = i2(2)

    nr=nrpts(ia)
    r=rmesh(1:nr,ia) !r-points
    do jlm = 1,lmmax4     ! l''m'' sum in notes
      work = 0.d0
      do ilm = 1,lmmax2   ! l'm' sum in notes
        tmprgaunt = rgaunt(ilm,lm2,jlm)
        if(abs(tmprgaunt) > ylmtol) then
        do a = 1,3
          tmprgaunt2 = rgaunt(lm1,ilmxyz(a),ilm)
          if(abs(tmprgaunt2) > ylmtol) then
            work(1:nr,a)=work(1:nr,a) + 1.d0/sqrt(3.d0)*tmprgaunt*tmprgaunt2*dphidr(1:nr,b1,l1,s1,ia)*phiref(1:nr,b2,l2,s2,ia)
          end if
          do c = 1,3
            do m3 = -l1,l1 !m''' in notes
              lm3 = lm2i(m3,l1) ! l_1 m'''
              tmprgaunt3 = rgaunt(lm3,ilmxyz(c),ilm)  
              if(abs(tmprgaunt3) > ylmtol) then 
                do b =1,3
                  work(1:nr,a) = work(1:nr,a) + ci/sqrt(3.d0)*tmprgaunt*tmprgaunt3*eps(a,b,c)*lorb(lm3,lm1,b)*phiref(1:nr,b1,l1,s1,ia)*phiref(1:nr,b2,l2,s2,ia)/r(1:nr)
                end do
              end if
            end do
          end do
        end do      
        end if
      end do

!     integrate work and put to gradnorm
      do a = 1,3
        gradnorm(a,jlm,iq) = radint(nr,work(1:nr,a),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
        gradbasis_lm(1:nr,a,jlm,iq)=work(1:nr,a)/(r(1:nr)**2)
      end do
!      if(abs(sum(gradnorm(:,jlm,iq))) > 0.001) write(*,'(2i4,6es16.8)') i2lm(:,jlm), (gradnorm(a,jlm,iq),a=1,3)
!      if(jlm == 1) then
!         write(*,'("gradnorm  =")')
!         write(*,'(2i4,6es16.8)') i2lm(:,jlm), (gradnorm(a,jlm,iq),a=1,3)
!      end if

!      work = 0.d0
!      work(1:nr,1)=phiref(1:nr,b1,l1,s1,ia)*phiref(1:nr,b2,l2,s2,ia)*rgaunt(lm1,lm2,jlm)
!      gfnorm(jlm,iq)=radint(nr,work(1:nr,1),drmesh(1:nr,ia))
    end do 
  end do
   
!  gradnorm(:,:,:)=real(gradnorm(:,:,:))
!  gradnorm2(:,:,:)=real(gradnorm2(:,:,:))

  end subroutine gradient_susc_basis
