  subroutine current_density(ia,nr,gfsum,curr_sum,morb_clust_sum,numpan,numrcut)
! Calculates the groundstate charge and spin currents and
! the net-current contribution to the orbital magnetic moment for atom ia

  use global
  use mod_derivative_panels

  implicit none

! atom number ia
  integer(kind=i4b), intent(in)       :: ia, nr
! density matrix
  complex(kind=c8b), intent(in)       :: gfsum(nlmsb,nlmsb)
! summing net current and orbital moment from net current
  real(kind=r8b),    intent(inout)    :: curr_sum(0:3,1:3),morb_clust_sum(1:3)
! --> Number of panels > 1
  integer(kind=i4b), intent(in)       :: numpan, numrcut(numpan+1) 
! -----------
  integer(kind=i4b)                   :: ir,i2(2),i3(3),k,n,ill
  integer(kind=i4b)                   :: i, ilm, il, im, ib, is
  integer(kind=i4b)                   :: j, jlm, jl, jm, jb, js
  integer(kind=i4b)                   :: klm, kl, km, mlm
  integer(kind=i4b)                   :: lmmaxP,lmmaxJ,lmmaxPdiv
  integer(kind=i4b)                   :: ilmxyz(1:3)
! Levi-Civita symbol
  integer(kind=i4b)                   :: eps(1:3,1:3,1:3),tmpeps
! dphidr = derivative of phi wrt r | r = radial points | radbasis= tmp array of phi | phir = phi/r
  real(kind=r8b)                      :: dphidr(1:nr,1:nbmax,0:nlmax,1:nsmax),phir(1:nr,1:nbmax,0:nlmax,1:nsmax),phi(1:nr,1:nbmax,0:nlmax,1:nsmax)
  real(kind=r8b)                      :: r(1:nr),radbasis(1:nr),dbasisdr(1:nr)
! curr_lm = decomposition of current into spherical harmonics | curr_lm_sum radial part is integrated | curr pointwise current
  complex(kind=c8b), allocatable      :: curr_lm(:,:,:,:),curr_lm_sum(:,:,:),curr(:,:,:,:)
! longditudinal and transverse components of the current
  complex(kind=c8b), allocatable      :: curr_lm_long(:,:,:,:),curr_lm_trans(:,:,:,:),curr_lm_sum_long(:,:,:),curr_lm_sum_trans(:,:,:)
! divergence of the current in spherical harmonics+integrated divergence | derivative of current_lm
  complex(kind=c8b), allocatable      :: curr_div_lm(:,:,:),curr_div_lm_sum(:,:),dcurr_lm_dr(:,:,:,:) 
  complex(kind=c8b), allocatable      :: curr_div_long_lm(:,:,:),curr_div_long_lm_sum(:,:),curr_div_trans_lm(:,:,:),curr_div_trans_lm_sum(:,:) 
! curl of the current
  complex(kind=c8b), allocatable      :: curr_curl_lm(:,:,:,:),curr_curl_lm_sum(:,:,:)
! array for r cross j (maybe magnetization density)
  complex(kind=c8b), allocatable      :: r_cross_j_lm(:,:,:,:)
! array for the matrix elements P(alpha,b,lm,lm1,r) lm1 has maximum l1=l+1
  complex(kind=c8b), allocatable      :: P(:,:,:,:,:),Plong(:,:,:,:,:),Ptrans(:,:,:,:,:)
! array for matrix elements for the divergence and curl
  complex(kind=c8b), allocatable      :: Pdiv(:,:,:,:,:),Pcurl(:,:,:,:),Pdiv_long(:,:,:,:,:),Pdiv_trans(:,:,:,:,:)
! orbital magnetization
  complex(kind=c8b)                   :: morb(1:3),morb_cluster(1:3)
  real(kind=r8b)                      :: morb_abs, morb_cluster_abs
! Helmholtz decomposition
  complex(kind=c8b),allocatable       :: phi_lm(:,:),A_lm(:,:,:)
  complex(kind=c8b)                   :: tmpVint(1:3), tmpSint(1:3), M_helm(1:3)
! constant pi
  real(kind=r8b), parameter           :: pi=4.d0*atan(1.d0)
! interpolation 
  complex(kind=c8b)                   :: curr_int(0:3,1:3,1:n_int,1:n_int,1:n_int)
! logicals
  logical                             :: lcurrentpts,lcurrentdiv,lcurrenthelmholtzdec,lrcrossj,loutputfiles
  real(kind=r8b)                      :: start,finish,norm
  character(len=1024)                 :: filename
  complex(kind=c8b), external         :: radint
  complex(kind=c8b)                   :: cone=(1.d0,0.d0),ci=(0.d0,1.d0)
  complex(kind=c8b)                   :: tmppsum(3),tmpfac,delta_gfsum,tmp_gij,tmp_gji,work(1:nr)
  real(kind=r8b)                      :: tmprgaunt,dri,dri1,rx,ry,rz,rxyz(1:3)
  integer(kind=i4b)                   :: ip, ist, ien

! #############################################################
! Options:

  lcurrentpts=.false.          !Calculation of the  current in each point of the grid (old option, better use interpolation)
  lcurrentdiv=.false.          !Calculation of the divergence
  lcurrenthelmholtzdec=.false. !Helmholtz decomposition
  lrcrossj=.false.             !r cross j is calculated
  loutputfiles=.false.         !output of current lm decomposition to files (general output in chargespdf.dat is not changed)

! #############################################################
! open iofile chargespdf.dat

  write(filename,'("chargespdf.dat")')
  open(file=filename,unit=iofile,status='old',access='append')
  write(iofile,'("GS current calculation results")')
  
! #############################################################
! write(*,*) "Current density calculation started"
  call cpu_time(start)

! initialize Levi-Civita symbol
  eps(:,:,:)=0
  eps(1,2,3)=1;eps(2,3,1)=1;eps(3,1,2)=1;
  eps(3,2,1)=-1;eps(1,3,2)=-1;eps(2,1,3)=-1;

!radial mesh
  r=rmesh(1:nr,ia) !r-points

! allocate P array (matrix elements of momentum operator acting on basis)
! lmmax of the second component is (l+2)**2 (combination of l and l'=1)
  i2=i2lm(:,lmmax)
  il=i2(2)
  lmmaxP=(il+2)**2 
  allocate(P(1:3,1:nr,1:nbmax,1:lmmax,1:lmmaxP))
  allocate(Ptrans(1:3,1:nr,1:nbmax,1:lmmax,1:lmmaxP),Plong(1:3,1:nr,1:nbmax,1:lmmax,1:lmmaxP))
! current_lm has maximum index l+1+l
  lmmaxJ=(il+il+2)**2 !! combined index l+1 and l
  allocate(curr_lm(0:3,1:3,1:nr,1:lmmaxJ),curr_lm_sum(0:3,1:3,1:lmmaxJ))
  allocate(curr(0:3,1:3,1:nr,1:nll),dcurr_lm_dr(0:3,1:3,1:nr,1:lmmaxJ))
  allocate(curr_lm_long(0:3,1:3,1:nr,1:lmmaxJ),curr_lm_sum_long(0:3,1:3,1:lmmaxJ),curr_lm_trans(0:3,1:3,1:nr,1:lmmaxJ),curr_lm_sum_trans(0:3,1:3,1:lmmaxJ))

! maximum of lm in P matrix of the divergence calculation
  lmmaxPdiv=(il+il+3)**2 !lmmaxJ + 1
  if (lcurrentdiv) then
    allocate(Pdiv(0:3,1:3,1:nr,1:lmmaxJ,1:lmmaxPdiv),curr_div_lm(0:3,1:nr,1:lmmaxPdiv),curr_div_lm_sum(0:3,1:lmmaxPdiv))
    allocate(curr_div_long_lm(0:3,1:nr,1:lmmaxPdiv),curr_div_long_lm_sum(0:3,1:lmmaxPdiv))
    allocate(curr_div_trans_lm(0:3,1:nr,1:lmmaxPdiv),curr_div_trans_lm_sum(0:3,1:lmmaxPdiv))
    allocate(Pdiv_trans(0:3,1:3,1:nr,1:lmmaxJ,1:lmmaxPdiv),Pdiv_long(0:3,1:3,1:nr,1:lmmaxJ,1:lmmaxPdiv))
    allocate(Pcurl(1:3,1:nr,1:lmmaxJ,1:lmmaxPdiv),curr_curl_lm(0:3,1:3,1:nr,1:lmmaxPdiv),curr_curl_lm_sum(0:3,1:3,1:lmmaxPdiv))
  end if
  if (lrcrossj) allocate(r_cross_j_lm(0:3,1:3,1:nr,1:lmmaxPdiv))

  if (lcurrenthelmholtzdec) then
!   arrays for Helmholtz decomposition
    allocate(phi_lm(1:nr,1:lmmax2),A_lm(1:3,1:nr,1:lmmax2)) 
  end if

! #############################################################
! calculate dphi/dr (and phi/r) for all basis functions and all radial points
! loop over panels if SOC from host
  radbasis(:)=0.d0
  dbasisdr(:)=0.d0
  do j=1,nlmsba(ia)
    !copied from get_xc routine
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    i2 = i2lm(:,jlm)
    jm = i2(1); jl = i2(2) 
    radbasis(:)=phiref(:,jb,jl,js,ia)/r(:)

! MODEL BASIS
!  radbasis(:)=0.d0
!  if(jb==1) then
!    do ir=1,nr
!      radbasis(ir)=r(ir)**2*exp(-r(ir))
!      radbasis(ir)=4.d0*3.0d0**(7/2)/sqrt(15.d0)/sqrt(sqrt(pi))*r(ir)**2*exp(-9.d0*r(ir)**2) !alpha=3**4
!    end do
!  end if

    call calc_derivative_panels(radbasis,dbasisdr,r,nr,numpan,numrcut)
    dphidr(1:nr,jb,jl,js)=dbasisdr(1:nr)

    do ir=1,nr
      phir(ir,jb,jl,js)=radbasis(ir)/r(ir)
      phi(ir,jb,jl,js)=radbasis(ir)
    end do
  end do

! #############################################################
! ilmx=4;ilmy=2;ilmz=3    ! (l,m)=(1,1)=x   (1,-1)=y   (1,0)=z
! ilmx=lm2i(1,1);ilmy=lm2i(-1,1);ilmz=lm2i(0,1)
  ilmxyz(1)=lm2i(1,1);ilmxyz(2)=lm2i(-1,1);ilmxyz(3)=lm2i(0,1)

  P=0.d0;Plong=0.d0;Ptrans=0.d0
! calculate P matrix for x,y,z
  do jlm=1,lmmaxP !sum over l_1,m_1 
    i2 = i2lm(:,jlm)
    jm = i2(1); jl = i2(2)
    do i=1,nlmsba(ia) !sum over basis l,m,b,s
      i3 = i2lmsb(:,i,ia)
      ib = i3(1); ilm = i3(2); is = i3(3)
      i2 = i2lm(:,ilm)
      im = i2(1); il = i2(2)
      if (il+jl == 1 .OR. abs(il-jl) == 1) then    
        !sum over m'
        tmppsum(:)=0.d0 
        do km=-il,il !sum over m' in range of l
          klm=lm2i(km,il)
          tmppsum(1)=tmppsum(1)+lorb(klm,ilm,3)*rgaunt(klm,ilmxyz(2),jlm)-lorb(klm,ilm,2)*rgaunt(klm,ilmxyz(3),jlm)
          tmppsum(2)=tmppsum(2)+lorb(klm,ilm,1)*rgaunt(klm,ilmxyz(3),jlm)-lorb(klm,ilm,3)*rgaunt(klm,ilmxyz(1),jlm)
          tmppsum(3)=tmppsum(3)+lorb(klm,ilm,2)*rgaunt(klm,ilmxyz(1),jlm)-lorb(klm,ilm,1)*rgaunt(klm,ilmxyz(2),jlm)
        end do
       
        
        do ir=1,nr
          Plong(:,ir,ib,ilm,jlm)=-ci/sqrt(3.d0)*rgaunt(ilm,ilmxyz(:),jlm)*dphidr(ir,ib,il,is)
          Ptrans(:,ir,ib,ilm,jlm)=-1.d0/sqrt(3.d0)*phir(ir,ib,il,is)*tmppsum(:) !Minus sign according to Manuels notes
          P(:,ir,ib,ilm,jlm)=Plong(:,ir,ib,ilm,jlm)+Ptrans(:,ir,ib,ilm,jlm)
  
        end do
      end if
    end do
  end do

! #############################################################
  curr_lm=0.d0;curr_lm_long=0.d0;curr_lm_trans=0.d0
! calculate the lm component of the current for all radial points
  do j=1,nlmsba(ia)                 !L_2 summation in notes
    i3 = i2lmsb(:,j,ia)
    jb = i3(1); jlm = i3(2); js = i3(3)
    i2 = i2lm(:,jlm)
    jm = i2(1); jl = i2(2)
    do klm=1,lmmaxP                !L_3 summation in notes
      i2 = i2lm(:,klm)
      km = i2(1); kl = i2(2)
      do i=1,nlmsba(ia)            !corresponds to L_1 summation in notes
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
        tmp_gij=gfsum(i,j)
        tmp_gji=gfsum(j,i)
        do mlm=1,lmmaxJ !! lmmaxJ is (lmax+lmaxP+1)^2 see notes
          tmprgaunt=rgaunt(jlm,klm,mlm)
          if(abs(tmprgaunt) > ylmtol) then
            do ir=1,nr
              do k=1,4
                delta_gfsum=tmp_gij*ds2c(k,js,is)-tmp_gji*ds2c(k,is,js)  
                if(abs(delta_gfsum) > ylmtol) then
                  tmpfac=phi(ir,jb,jl,js)*delta_gfsum*tmprgaunt
                  curr_lm_long(mod(k,4),:,ir,mlm)=curr_lm_long(mod(k,4),:,ir,mlm)+Plong(:,ir,ib,ilm,klm)*tmpfac
                  curr_lm_trans(mod(k,4),:,ir,mlm)=curr_lm_trans(mod(k,4),:,ir,mlm)+Ptrans(:,ir,ib,ilm,klm)*tmpfac
                  curr_lm(mod(k,4),:,ir,mlm)=curr_lm(mod(k,4),:,ir,mlm)+P(:,ir,ib,ilm,klm)*tmpfac
                end if
!  model calculation       2*xi/gamma=1 => atan(..)=pi/4.d0
!              if (ib==1 .and. ilm==5 .and. jlm==9 .and. js==1) then
!                curr_lm(:,ir,mlm)=curr_lm(:,ir,mlm)+1.d0/2.d0*P(:,ir,ib,ilm,klm)*phi(ir,jb,jl,js)*ci*pi/4.d0*rgaunt(jlm,klm,mlm)
!              else if(ib==1 .and. ilm==9 .and. jlm==5 .and. js==1) then
!                curr_lm(:,ir,mlm)=curr_lm(:,ir,mlm)-1.d0/2.d0*P(:,ir,ib,ilm,klm)*phi(ir,jb,jl,js)*ci*pi/4.d0*rgaunt(jlm,klm,mlm)
!              end if
              end do
            end do
          end if
        end do
      end do    
    end do
  end do
  

! #############################################################
! integrate the radial part of j_lm
  curr_lm_sum=0.d0;curr_lm_sum_long=0.d0;curr_lm_sum_trans=0.d0
  do mlm=1,lmmaxJ
    do n=1,3
      !integrate r^2*j_lm
      do k=0,3
        curr_lm_sum(k,n,mlm)=radint(nr,curr_lm(k,n,1:nr,mlm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
        curr_lm_sum_long(k,n,mlm)=radint(nr,curr_lm_long(k,n,1:nr,mlm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
        curr_lm_sum_trans(k,n,mlm)=radint(nr,curr_lm_trans(k,n,1:nr,mlm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
      end do
    end do
  end do
  
! #############################################################
! write integrated paramagnetic current to file

  write(iofile,'("  curr_para_lm_sum_0=")')
  do mlm=1,lmmaxJ
    norm = sqrt(dot_product(real(curr_lm_sum(0,1:3,mlm)),real(curr_lm_sum(0,1:3,mlm))))
    if (norm > atol) then
      write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_lm_sum(0,1:3,mlm))
    end if
  end do
  
  write(iofile,'("  curr_para_lm_sum_1=")')
  do mlm=1,lmmaxJ
    norm = sqrt(dot_product(real(curr_lm_sum(1,1:3,mlm)),real(curr_lm_sum(1,1:3,mlm))))
    if (norm > atol) then
      write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_lm_sum(1,1:3,mlm))
    end if
  end do
  
  write(iofile,'("  curr_para_lm_sum_2=")')
  do mlm=1,lmmaxJ
    norm = sqrt(dot_product(real(curr_lm_sum(2,1:3,mlm)),real(curr_lm_sum(2,1:3,mlm))))
    if (norm > atol) then
      write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_lm_sum(2,1:3,mlm))
    end if
  end do
  
  write(iofile,'("  curr_para_lm_sum_3=")')
  do mlm=1,lmmaxJ
    norm = sqrt(dot_product(real(curr_lm_sum(3,1:3,mlm)),real(curr_lm_sum(3,1:3,mlm))))
    if (norm > atol) then
      write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_lm_sum(3,1:3,mlm))
    end if
  end do
  
  4000 format(2i4,6e18.9)
  if(lcurrentoutput) then 
    write(filename,"(A20,I0.3,A4)") "current_para_lm_sum_",ia,".dat"
    open(unit=10,file=filename)
    do mlm=1,lmmaxJ
      write(10,4000) i2lm(:,mlm), curr_lm_sum(0,1,mlm),curr_lm_sum(0,2,mlm),curr_lm_sum(0,3,mlm)
    end do
    close(10)
  end if
  
! #############################################################
! orbital magnetization
  morb(:)=0.d0
  do k=1,3
    do j=1,3
      do i=1,3
        morb(i)=morb(i)+0.5d0*eps(i,j,k)*sqrt(1.d0/3.d0)*rgaunt(ilmxyz(j),ilmxyz(j),1)*radint(nr,curr_lm(0,k,1:nr,ilmxyz(j))*r(1:nr)**3,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
      end do
    end do
  end do
  
! write(*,'("morb from current = ",6e18.9)') morb
  
  write(iofile,'("  mox         =",100f16.8)') real(morb(1))
  write(iofile,'("  moy         =",100f16.8)') real(morb(2))
  write(iofile,'("  moz         =",100f16.8)') real(morb(3))
  
! contribution from the net current flowing out of the sphere 
! Positions have to be read in from positions.dat file into array ri(1:3,ia)
  morb_cluster(:)=0.d0
  if(lpositions) then
    rxyz(1)=ri(1,ia)!-(ri(1,1)+ri(1,2)+ri(1,3))/3.d0
    rxyz(2)=ri(2,ia)!-(ri(2,1)+ri(2,2)+ri(2,3))/3.d0
    rxyz(3)=ri(3,ia)
    do k=1,3
      do j=1,3
        do i=1,3
          morb_cluster(i)=morb_cluster(i)+0.5d0*eps(i,j,k)*rxyz(j)*curr_lm_sum(0,k,1)
        end do
      end do
    end do
!   write(*,'("morb net current contribution = ",6e18.9)') morb_cluster
  
    write(iofile,'("  mo_non_loc_x=",100f16.8)') real(morb_cluster(1))
    write(iofile,'("  mo_non_loc_y=",100f16.8)') real(morb_cluster(2))
    write(iofile,'("  mo_non_loc_z=",100f16.8)') real(morb_cluster(3))
  
    morb_abs = sqrt(dot_product(real(morb(1:3)),real(morb(1:3))))
    morb_cluster_abs = sqrt(dot_product(real(morb_cluster(1:3)),real(morb_cluster(1:3))))
    write(440,'(100e18.9)') ri(1,ia),ri(2,ia),ri(3,ia),morb_abs,morb_cluster_abs
  end if
  
! Summation over the whole cluster
  morb_clust_sum(:)=morb_clust_sum(:)+real(morb_cluster(:))
  
! calculation of r cross j
  if(lrcrossj) then 
    r_cross_j_lm=0.d0
    do ilm=1,lmmaxPdiv
      do jlm=1,lmmaxJ
        do k=1,3
          do j=1,3
            tmprgaunt=1.d0/sqrt(3.d0)*rgaunt(ilmxyz(j),jlm,ilm)
            if(abs(tmprgaunt) > ylmtol) then
              do ir=1,nr    
                do i=1,3
                  r_cross_j_lm(0,i,ir,ilm)=r_cross_j_lm(0,i,ir,ilm)+tmprgaunt*eps(i,j,k)*r(ir)*curr_lm(0,k,ir,jlm)
                end do
              end do
            end if
          end do
        end do
      end do
    end do
!  Integrated magnetization is the same as above 
!  do i=1,3
!    morb(i)=radint(nr,r_cross_j_lm(i,1:nr,1)*r(1:nr)**2,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
!  end do
!  write(*,'("morb from current = ",6e18.9)') morb
  end if
  
! #############################################################
! SOC and Zeeman contribution
  
  if(lcurrentsoc) then
    curr_lm=0.d0
    call current_soc(curr_lm(:,:,:,:),nr,lmmaxJ,r,phi,eps,ia,numpan,numrcut)
  end if
  
  if(lcurrentzeeman) then 
    call current_zeeman(curr_lm(0,:,:,:),nr,lmmaxJ,r,phi,eps,ia)
  end if
  
  
! #############################################################
! integrate the radial part of j_lm (with zeeman and soc contribution now)
  curr_lm_sum=0.d0
  do mlm=1,lmmaxJ
    do n=1,3
      !integrate r^2*j_lm
      do i=0,3
        curr_lm_sum(i,n,mlm)=radint(nr,curr_lm(i,n,1:nr,mlm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
      end do
    end do
  end do
  
! Summation of the net current in the whole cluster
  curr_sum(:,:)=curr_sum(:,:)+real(curr_lm_sum(:,:,1))
  
! #############################################################
! calculate current in each point
  if(lcurrentpts) then
    do mlm=1,lmmaxJ
      do ill=1,nll
        do ir=1,nr
          curr(0,:,ir,ill)=curr(0,:,ir,ill)+curr_lm(0,:,ir,mlm)*ylm(ill,mlm)
        end do
      end do
    end do
  end if
  
! #############################################################
! calculate the divergence of the current
! first the radial derivative of j_lm is calculated, then the P matrix for the divergence
  
  if(lcurrentdiv) then 
!    !calculate first the derivative of all j_lm's
!    !forward differences
!    dcurr_lm_dr(:,:,1,:)= (curr_lm(:,:,2,:)-curr_lm(:,:,1,:))/(r(2)-r(1))
!    !central differences
!    do ir=2,nr-1
!      dri=r(ir+1)-r(ir)
!      dri1=r(ir)-r(ir-1)
!  !    dcurr_lm_dr(:,ir,:)=-dri/dri1/(dri+dri1)*curr_lm(:,ir-1,:)
!  !    dcurr_lm_dr(:,ir,:)=dcurr_lm_dr(:,ir,:)+(dri-dri1)/dri/dri1*curr_lm(:,ir,:)
!  !    dcurr_lm_dr(:,ir,:)=dcurr_lm_dr(:,ir,:)+dri1/dri/(dri+dri1)*curr_lm(:,ir+1,:)
!      dcurr_lm_dr(:,:,ir,:)=0.5d0*(curr_lm(:,:,ir+1,:)-curr_lm(:,:,ir,:))/(r(ir+1)-r(ir))
!      dcurr_lm_dr(:,:,ir,:)= dcurr_lm_dr(:,:,ir,:) + 0.5d0*(curr_lm(:,:,ir,:)-curr_lm(:,:,ir-1,:))/(r(ir)-r(ir-1))
!    end do
!    !backward difference
!    dcurr_lm_dr(:,:,nr,:)= (curr_lm(:,:,nr,:)-curr_lm(:,:,nr-1,:))/(r(nr)-r(nr-1))
    do ilm=1,lmmaxJ
      do i=1,3
        do k =0,3
          call calc_derivative_panels(curr_lm(k,i,1:nr,ilm),work,r,nr,numpan,numrcut)
          dcurr_lm_dr(k,i,1:nr,ilm)=work(1:nr)
        end do
      end do
    end do
    
    Pdiv=0.d0;Pdiv_long=0.d0;Pdiv_trans=0.d0
    !calculate Pdiv matrix for x,y,z
    do jlm=1,lmmax2 !lmmaxPdiv !sum over l_1,m_1 
      i2 = i2lm(:,jlm)
      jm = i2(1); jl = i2(2)
      do ilm=1,lmmax2  !lmmaxJ !sum over basis l,m,b,s
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
       !sum over m'
        tmppsum(:)=0.d0 
        do km=-il,il !sum over m' in range of l
          !klm=lm2i(km,il)
          klm = (il+1)**2-il+km
          do j=1,3
            tmprgaunt=rgaunt(klm,ilmxyz(j),jlm)
            if(abs(tmprgaunt) > ylmtol) then
              do k=1,3
                do i=1,3
                  tmpeps=eps(i,j,k)
                  if(abs(tmpeps) > ylmtol) then
                    tmppsum(i)=tmppsum(i)+tmpeps*tmprgaunt*lorb(klm,ilm,k)
                  end if
                end do
              end do
            end if
          end do
        end do
        do ir=1,nr
          do i=1,3
!            Pdiv_long(i,ir,ilm,jlm) = -ci/sqrt(3.d0)*rgaunt(ilm,ilmxyz(i),jlm)*dcurr_lm_dr(i,ir,ilm)
!     Pdiv_long is not the same as in the notes !!!!
!     This to test the stability of the algorithm
!     instead of calculating the divergence in each radial point we look at the integrated divergence and make use of partial integration
!     so only the radial integrated part is meaningful from here on      
            Pdiv_long(:,i,ir,ilm,jlm) = -ci/sqrt(3.d0)*rgaunt(ilm,ilmxyz(i),jlm)*curr_lm(:,i,ir,ilm)
            Pdiv_trans(:,i,ir,ilm,jlm) =-1.d0/sqrt(3.d0)*curr_lm(:,i,ir,ilm)/r(ir)*tmppsum(i)
            Pdiv(:,i,ir,ilm,jlm) = Pdiv_long(:,i,ir,ilm,jlm)+Pdiv_trans(:,i,ir,ilm,jlm)
          end do
        end do
      end do
    end do
  
    !Calculate div j in basis of spherical harmonics
    curr_div_lm=0.d0;curr_div_long_lm=0.d0;curr_div_trans_lm=0.d0
    do jlm=1,lmmax2!lmmaxPdiv
      i2 = i2lm(:,jlm)
      jm = i2(1); jl = i2(2)
      do ilm=1,lmmax2!lmmaxJ !sum over basis l,m,b,s
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
        do ir=1,nr
          do n=1,3
            curr_div_long_lm(:,ir,jlm)=curr_div_long_lm(:,ir,jlm)+ci*Pdiv_long(:,n,ir,ilm,jlm)
            curr_div_trans_lm(:,ir,jlm)=curr_div_trans_lm(:,ir,jlm)+ci*Pdiv_trans(:,n,ir,ilm,jlm)
            curr_div_lm(:,ir,jlm)=curr_div_lm(:,ir,jlm)+ci*Pdiv(:,n,ir,ilm,jlm)
          end do
        end do
      end do
    end do
    
    !integrate the lm elements of the divergence over the radial points
    curr_div_lm_sum=0.d0;curr_div_long_lm_sum=0.d0;curr_div_trans_lm_sum=0.d0
    do ilm=1,lmmax2 !lmmaxPdiv
      do k=0,3
        curr_div_long_lm_sum(k,ilm)=radint(nr,curr_div_long_lm(k,1:nr,ilm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
  !      curr_div_long_lm_sum(:,ilm)=r(nr)**2*curr_div_long_lm(:,nr,ilm)-r(1)**2*curr_div_long_lm(:,1,ilm)-radint(nr,2*curr_div_long_lm(:,1:nr,ilm)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
        curr_div_trans_lm_sum(k,ilm)=radint(nr,curr_div_trans_lm(k,1:nr,ilm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
  !      curr_div_lm_sum(ilm)=radint(nr,curr_div_lm(1:nr,ilm)*r(1:nr)*r(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
        curr_div_lm_sum(k,ilm)=curr_div_long_lm_sum(k,ilm)+curr_div_trans_lm_sum(k,ilm)
      end do
    end do
    
    write(iofile,'("  curr_para_div_lm_sum_0=")')
    do mlm=1,lmmax2
      norm = abs(real(curr_div_lm_sum(0,mlm)))
      if (norm > atol) then
        write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_div_lm_sum(0,mlm))
      end if
    end do
      
    write(iofile,'("  curr_para_div_lm_sum_1=")')
    do mlm=1,lmmax2
      norm = abs(real(curr_div_lm_sum(1,mlm)))
      if (norm > atol) then
        write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_div_lm_sum(1,mlm))
      end if
    end do
      
    write(iofile,'("  curr_para_div_lm_sum_2=")')
    do mlm=1,lmmax2
      norm = abs(real(curr_div_lm_sum(2,mlm)))
      if (norm > atol) then
        write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_div_lm_sum(2,mlm))
      end if
    end do
      
    write(iofile,'("  curr_para_div_lm_sum_3=")')
    do mlm=1,lmmax2
      norm = abs(real(curr_div_lm_sum(3,mlm)))
      if (norm > atol) then
        write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_div_lm_sum(3,mlm))
      end if
    end do
    
  end if
    
! #############################################################
! interpolation on a regular grid
  if (lcurrentint) then
    do k=0,3
      call current_interpolation(curr_lm(k,:,:,:),lmmaxJ,r,nr,curr_int(k,:,:,:,:),ia,k,numpan,numrcut)
    end do
  end if
  
! #############################################################
! write current to file (lm decomposition)
  if(loutputfiles) then
    
    write(filename,"(A11,I0.3,A4)") "current_lm_",ia,".dat"
    open(unit=9,file=filename)
    1000 format(2i4,7e18.9)
    2000 format(2i4)
    
    write(9,*) "# lmmax nr"
    write(9,2000) lmmaxJ, nr
    write(9,*) "# m, l , r ,j_x,j_y,j_z"
    do mlm=1,lmmaxJ
      do ir=1,nr
        write(9,1000) i2lm(:,mlm) , r(ir), curr_lm(0,1,ir,mlm),curr_lm(0,2,ir,mlm),curr_lm(0,3,ir,mlm)
      end do
    end do
    
    close(9)
    
  
!   #############################################################
!   write integrated current to file
    
    write(filename,"(A15,I0.3,A4)") "current_lm_sum_",ia,".dat"
    open(unit=10,file=filename)
    
    write(iofile,'("  curr_lm_sum=")')
    do mlm=1,lmmaxJ
      write(10,4000) i2lm(:,mlm), curr_lm_sum(0,1,mlm),curr_lm_sum(0,2,mlm),curr_lm_sum(0,3,mlm)
      norm = sqrt(dot_product(real(curr_lm_sum(0,1:3,mlm)),real(curr_lm_sum(0,1:3,mlm))))
      if (norm > atol) then
        write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_lm_sum(0,1:3,mlm))
      end if
    end do
    
    close(10)
    
  end if
  
! #############################################################
! write net current file
  
  rx=ri(1,ia)
  ry=ri(2,ia)
  rz=ri(3,ia)
  write(330,'(100e18.9)') rx,ry,rz,curr_lm_sum(0,1,1),curr_lm_sum(0,2,1),curr_lm_sum(0,3,1)
  write(331,'(100e18.9)') rx,ry,rz,curr_lm_sum(1,1,1),curr_lm_sum(1,2,1),curr_lm_sum(1,3,1)
  write(332,'(100e18.9)') rx,ry,rz,curr_lm_sum(2,1,1),curr_lm_sum(2,2,1),curr_lm_sum(2,3,1)
  write(333,'(100e18.9)') rx,ry,rz,curr_lm_sum(3,1,1),curr_lm_sum(3,2,1),curr_lm_sum(3,3,1)
    
! #############################################################
! write integrated transversal current to file
 
! write(filename,"(A21,I0.3,A4)") "current_lm_trans_sum_",ia,".dat"
! open(unit=10,file=filename)
!
! do mlm=1,lmmaxJ
!   write(10,4000) i2lm(:,mlm), curr_lm_sum_trans(1,mlm),curr_lm_sum_trans(2,mlm),curr_lm_sum_trans(3,mlm)
! end do
!
! close(10)

! #############################################################
! write integrated longditudinal current to file

! write(filename,"(A20,I0.3,A4)") "current_lm_long_sum_",ia,".dat"
! open(unit=10,file=filename)
!
! do mlm=1,lmmaxJ
!    write(10,4000) i2lm(:,mlm), curr_lm_sum_long(1,mlm),curr_lm_sum_long(2,mlm),curr_lm_sum_long(3,mlm)
! end do

! close(10)


! #############################################################
! write integrated divergence of the current to file
  if(lcurrentdiv) then 
    write(filename,"(A19,I0.3,A4)") "current_div_lm_sum_",ia,".dat"
    open(unit=10,file=filename)
  
    write(iofile,'("  curr_div_lm_sum=")')
    do mlm=1,lmmax2!lmmaxPdiv
      write(10,4000) i2lm(:,mlm), curr_div_lm_sum(0,mlm)
      norm = abs(real(curr_div_lm_sum(0,mlm)))
      if (norm > atol) then
        write(iofile,'(2i4,4f16.8)') i2lm(:,mlm), real(curr_div_lm_sum(0,mlm))
      end if
    end do
  
    close(10)
    
!   #############################################################
!   write integrated longditudinal divergence of the current to file
    
    write(filename,"(A24,I0.3,A4)") "current_div_long_lm_sum_",ia,".dat"
    open(unit=10,file=filename)
    
    do mlm=1,lmmax2!lmmaxPdiv
      write(10,4000) i2lm(:,mlm), curr_div_long_lm_sum(0,mlm)
    end do
    
    close(10)
    
!   #############################################################
!   write integrated transversal divergence of the current to file
  
    write(filename,"(A25,I0.3,A4)") "current_div_trans_lm_sum_",ia,".dat"
    open(unit=10,file=filename)
    
    do mlm=1,lmmax2!lmmaxPdiv
      write(10,4000) i2lm(:,mlm), curr_div_trans_lm_sum(0,mlm)
    end do
    
    close(10)
    !#############################################################
    !write divergence of a specific lm current component to file
    !
    !write(filename,"(A15,I0.3,A4)") "current_div_lm_",ia,".dat"
    !open(unit=10,file=filename)
    !
    !do ir=1,nr!lmmaxPdiv
    !  write(10,'(i3,2e18.9)') ir, curr_div_lm(ir,lm2i(-4,4))!component m=-4,l=4
    !end do
    !
    !close(10)
  end if
  !#############################################################
  
  deallocate(P,Ptrans,Plong)
  deallocate(curr_lm,curr_lm_sum)
  deallocate(curr,dcurr_lm_dr)
  deallocate(curr_lm_long,curr_lm_sum_long,curr_lm_trans,curr_lm_sum_trans)
  if (lcurrentdiv) then 
    deallocate(Pdiv,curr_div_lm,curr_div_lm_sum)
    deallocate(curr_div_long_lm,curr_div_long_lm_sum)
    deallocate(curr_div_trans_lm,curr_div_trans_lm_sum)
    deallocate(Pdiv_trans,Pdiv_long)
    deallocate(Pcurl,curr_curl_lm,curr_curl_lm_sum)
  end if
  if (lrcrossj)  deallocate(r_cross_j_lm)
  if (lcurrenthelmholtzdec) deallocate(phi_lm,A_lm) 
! write(*,*) "Current density calculation finished."
  call cpu_time(finish)
! write(*,'(/," Current calculation   time=",f10.3," s",/)') finish-start 
  close(iofile)
  end subroutine current_density
