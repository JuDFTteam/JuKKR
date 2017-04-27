  subroutine current_soc(curr_lm,nr,lmmaxJ,r,phi,eps,ia,numpan,numrcut)
  
  use global
  
  implicit none
  
  integer(kind=i4b), intent(in) :: nr, lmmaxJ, eps(1:3,1:3,1:3),ia
  complex(kind=c8b), intent(inout) :: curr_lm(0:3,1:3,1:nr,1:lmmaxJ)
  real(kind=r8b), intent(in) :: r(1:nr), phi(1:nr,1:nbmax,0:nlmax,1:nsmax)
! --> Number of panels > 1
    integer(kind=i4b), intent(in) :: numpan, numrcut(numpan+1) 
! ------------------------------------------------------------------------------------
  complex(kind=c8b) :: curr_soc_lm(0:3,1:3,1:nr,1:lmmaxJ), paulimat(1:3,1:2,1:2),curr_soc_lm_sum(0:3,1:3,1:lmmaxJ)
  real(kind=r8b)    :: prefac, E_soc(1:3,1:nr), tmpgaunt, tmpgaunt2,norm
  integer(kind=i4b) :: ir,i2(2),i3(3),n,ill
  integer(kind=i4b) :: i, ilm, il, im, ib, is
  integer(kind=i4b) :: j, jlm, jl, jm, jb, js
  integer(kind=i4b) :: k, klm, kl, km, mlm, a, b, c
  integer(kind=i4b) :: ilmxyz(1:3), tmpeps
  complex(kind=c8b) :: cone=(1.d0,0.d0),ci=(0.d0,1.d0)
  character(len=1024) :: filename
  complex(kind=c8b), external :: radint
! Complex SOC potential
  complex(kind=c8b) :: vsoc(nr)
! Speed of light
  real(kind=r8b), parameter :: clight = 274.0720442d0
! ------------------------------------------------------------------------------------
  
  ilmxyz(1)=lm2i(1,1);ilmxyz(2)=lm2i(-1,1);ilmxyz(3)=lm2i(0,1)
  
  
! prefactor e^2*hbar*sqrt(4pi/3)  at this point e=sqrt(2),hbar=1,m=1/2
! e maybe in scalar potential
  prefac = 1.d0/sqrt(3.d0)
  
! Spin-orbit electric field
  call build_vsoc(clight,esusc(nescf),zat(ia),nrpts(ia),rmesh(:,ia),vr(:,ia),socscaling(ia),vsoc(:),numpan,numrcut)
  
! Charge current correction
  
  curr_soc_lm=0.d0
! outer loop l'm'
  do ilm=1,lmmaxJ
    do jlm=1,lmmax2  !innerloop lm
      do j=1,3
        tmpgaunt=rgaunt(jlm,ilmxyz(j),ilm)
        if(abs(tmpgaunt) > ylmtol) then 
          do a=1,3
            do ir=1,nr
!              tmpeps=eps(a,j,k)
!              if(abs(tmpeps) > ylmtol) then
                do k=1,3
                  tmpeps=eps(a,j,k)
                  curr_soc_lm(0,k,ir,ilm)=curr_soc_lm(0,k,ir,ilm)+prefac*r(ir)*vsoc(ir)*rho_lm(ir,jlm,a,ia)*tmpeps*tmpgaunt
                end do
!              end if
            end do
          end do
        end if
      end do
    end do
  end do
  
  
! outer loop l'm'
  do ilm=1,lmmaxJ
    do jlm=1,lmmax2  !innerloop lm
      do j=1,3
        tmpgaunt=rgaunt(jlm,ilmxyz(j),ilm)
        if(abs(tmpgaunt) > ylmtol) then 
          do ir=1,nr
!            tmpeps=eps(a,j,k)
!            if(abs(tmpeps) > ylmtol) then
              do k=1,3
                do a=1,3
                  !Different signs in front of prefac have to be checked again!!!!
                  tmpeps=eps(a,j,k)
                  curr_soc_lm(a,k,ir,ilm)=curr_soc_lm(a,k,ir,ilm)+prefac*r(ir)*vsoc(ir)*rho_lm(ir,jlm,0,ia)*tmpeps*tmpgaunt
                end do
              end do
  !          end if
          end do
        end if
      end do
    end do
  end do
  
  
  curr_lm(:,:,:,:)=curr_lm(:,:,:,:)+curr_soc_lm(:,:,:,:)
  
  
! #############################################################
! integrate the radial part of j_soc_lm
  curr_soc_lm_sum=0.d0
  do n=1,3
    do mlm=1,lmmaxJ
      !integrate r^2*j_lm
      do a=0,3
        curr_soc_lm_sum(a,n,mlm)=radint(nr,curr_soc_lm(a,n,1:nr,mlm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
      end do
    end do
  end do
  
! #############################################################
! write integrated current to file
  
! write(filename,"(A19,I0.3,A4)") "current_soc_lm_sum_",ia,".dat"
! open(unit=10,file=filename)
! 4000 format(2i4,6e18.9)
  
  do a=0,3
    write(iofile,'("  curr_soc_lm_sum_",I1,"=")') a
    do mlm=1,lmmaxJ
!     write(10,4000) i2lm(:,mlm), curr_soc_lm_sum(a,1,mlm),curr_soc_lm_sum(a,2,mlm),curr_soc_lm_sum(a,3,mlm)
      norm = sqrt(dot_product(real(curr_soc_lm_sum(a,1:3,mlm)),real(curr_soc_lm_sum(a,1:3,mlm))))
      if (norm > atol) then
        write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_soc_lm_sum(a,1:3,mlm))
      end if
    end do
  end do
  
! close(10)
  
  end subroutine current_soc
