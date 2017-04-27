  subroutine overlaps_susc2()
! Analysis of the product basis for the susceptibility
! Orthogonalization of radial basis functions for non-zero Gaunts and spinors
  use global

  implicit none

  logical,        parameter :: printbasis = .false.
  real(kind=r8b), parameter :: tol = 1.d-3
  real(kind=r8b), parameter :: fourpi = 16.d0*atan(1.d0)
  integer(kind=i4b) :: i3(3), i4(4), nb, nr, i, j, lm, ib, is, ievals, ir, ilm, iq, jq, ia, ia2
  integer(kind=i4b) :: nbnew, ngaunt(nasusc2), nden0, nden1, ngf0, ngf1
  integer(kind=i4b) :: nlmpairs(lmmax2,4), nonzero(3,lmmax2,lmmax2,4)
  integer(kind=i4b) :: q1, b1, lm1, s1
  integer(kind=i4b) :: q2, b2, lm2, s2
  integer(kind=i4b) :: q3, b3, lm3, s3
  integer(kind=i4b) :: q4, b4, lm4, s4
  integer(kind=i4b) :: info, lrwork, liwork
  integer(kind=i4b), allocatable :: iwork(:), lmsb2i2(:,:)
  real(kind=r8b),    allocatable :: basis(:,:)
  complex(kind=c8b) :: work(nrmax), norm, rho2, rho2proj, rho, rhoc
  real(kind=r8b)    :: dr(nrmax), doublegaunt, doublegaunt0, tmp, start, finish, ram
  complex(kind=c8b), external :: radint
! Mpi 
  integer(kind=i4b) :: my_rank
  my_rank = 0

! Susceptibility?
  if (.not.lsusc) return

  call cpu_time(start)
! Storage for number of basis functions per atom
! max # of basis functions for density
! *****************************************
  nbmax2 = maxval(sum(iwsusc(:,:,iasusc2(1:nasusc2))**2,dim=1))
  write(*,'("overlaps_susc2: nbmax2 should be ",i8)') nbmax2
  if (allocated(nalmsbden)) deallocate(nalmsbden)
  allocate(nalmsbden(nasusc2))
! how many radial basis functions per lms for density
  if (allocated(iwsusc2)) deallocate(iwsusc2)
  allocate(iwsusc2(lmmax2,nsmax2,nasusc2))
! the actual radial basis functions
  if (allocated(suscbasis)) deallocate(suscbasis)
  allocate(suscbasis(nrmax,nbmax2,lmmax2,4,nasusc2))
  iwsusc2 = 0; nalmsbden = 0
! **************
  do ia2=1,nasusc2
! **************
    ia = iasusc2(ia2)
!   First determine which angular momentum couplings are allowed
    write(iodb,'("overlaps_susc2: ia=",i4,"  nonzero Gaunts")') ia
    nonzero = 0; nlmpairs = 0
    do is=1,nsmax2
      s1 = i2is(1,is); s2 = i2is(2,is)
      do lm=1,lmmax0
        ib = 0
        do lm2=1,lmmax
          do lm1=1,lmmax
!            write(*,'("lm1,lm2,lm=",4i4)') lm1, lm2, lm
            if (abs(rgaunt(lm1,lm2,lm)) > ylmtol) then
              ib = ib + 1
              nb = iwsusc(i2lm(2,lm1),s1,ia)*iwsusc(i2lm(2,lm2),s2,ia)
              nonzero(:,ib,lm,is) = (/lm1,lm2,nb/)
!              write(iodb,'(10i4)') ia, i2is(:,is), i2lm(:,lm), i2lm(:,lm1), i2lm(:,lm2), nb
            end if
          end do
        end do
        nlmpairs(lm,is) = ib
        write(iodb,'("overlaps_susc2: ia,is,lm=",4i4," nblm=",2i4)') ia, is, i2lm(:,lm), ib, sum(nonzero(3,1:ib,lm,is))
      end do
    end do
!   Now construct the basis
    nb = maxval(sum(nonzero(3,:,:,:),dim=1))
    nbnew = maxval(sum(iwsusc(:,:,ia)*iwsusc(:,:,ia),dim=1))
    write(iodb,'("overlaps_susc2: nbmax,nbnew=",2i4)') nb, nbnew
    allocate(basis(nrmax,nb))
    nr = nrpts(ia); dr(1:nr) = drmesh(1:nr,ia)/rmesh(1:nr,ia)**2
    do is=1,nsmax2
      s1 = i2is(1,is); s2 = i2is(2,is)
      do lm=1,lmmax0
        ib = 0
!       products of regular wfns
        do i=1,nlmpairs(lm,is)
          lm1 = nonzero(1,i,lm,is); lm2 = nonzero(2,i,lm,is)
          do b2=1,iwsusc(i2lm(2,lm2),s2,ia)
            do b1=1,iwsusc(i2lm(2,lm1),s1,ia)
              ib = ib + 1
              basis(1:nr,ib) = phiref(1:nr,b1,i2lm(2,lm1),s1,ia)*phiref(1:nr,b2,i2lm(2,lm2),s2,ia)
            end do
          end do
        end do
!       basis
        if (ib > 0) call find_basis2(ib,nb,basis,nr,nrmax,dr,tol,npanat(ia),ircutat(:,ia),printbasis,my_rank)
        write(iodb,'("overlaps_susc2: ia,is,lm=",4i4," nblm=",2i4)') ia, is, i2lm(:,lm), nlmpairs(lm,is), ib
!       save basis
        iwsusc2(lm,is,ia2) = min(ib,nbnew)
        do i=1,iwsusc2(lm,is,ia2)
!          write(*,'(4i4)') i, lm, is, ia2
          suscbasis(1:nr,i,lm,is,ia2) = basis(1:nr,i)
        end do
!       test output
        if (lm == 1 .and. is == 4) then
          do ir=1,nr
            write(999,'(100es16.8)') rmesh(ir,ia), basis(ir,1:ib)
          end do
        end if
!       compare with GS density
        do i=1,nlmpot(ia)
          if (i2lmpot(i,ia) == lm) then
!         ********  up  ******************
          if (s1 == 2 .and. s2 == s1) then
            work(1:nr) = 0.5d0*(old_rho2ns(1:nr,lm,2,ia) + old_rho2ns(1:nr,lm,1,ia))
            tmp = radint(nr,work,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
!            write(iodb,'("nup=",f16.8)') real(tmp)
            norm = 0.d0
            do j=1,ib
              work(1:nr) = basis(1:nr,j)
              tmp = radint(nr,work,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
              work(1:nr) = 0.5d0*(old_rho2ns(1:nr,lm,2,ia) + old_rho2ns(1:nr,lm,1,ia))*basis(1:nr,j)
              norm = norm + radint(nr,work,dr,npanat(ia),ircutat(:,ia))*tmp
            end do
            write(iodb,'("GS up lm vs projected: (l,m)=",2i4,2f16.8)') i2lm(:,lm), 0.5d0*(gs_qlm(lm,ia) + gs_mlm(lm,ia)), real(norm)
          end if
!         ********  down  ****************
          if (s1 == 1 .and. s2 == s1) then
            work(1:nr) = 0.5d0*(old_rho2ns(1:nr,lm,1,ia) - old_rho2ns(1:nr,lm,2,ia))
            tmp = radint(nr,work,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
!            write(iodb,'("ndn=",f16.8)') real(tmp)
            norm = 0.d0
            do j=1,ib
              work(1:nr) = basis(1:nr,j)
              tmp = radint(nr,work,drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
              work(1:nr) = 0.5d0*(old_rho2ns(1:nr,lm,1,ia) - old_rho2ns(1:nr,lm,2,ia))*basis(1:nr,j)
              norm = norm + radint(nr,work,dr,npanat(ia),ircutat(:,ia))*tmp
            end do
            write(iodb,'("GS dn lm vs projected: (l,m)=",2i4,2f16.8)') i2lm(:,lm), 0.5d0*(gs_qlm(lm,ia) - gs_mlm(lm,ia)), real(norm)
          end if
!         ******
          end if
        end do
      end do
    end do
    deallocate(basis)
!   Total number of basis functions for density
    nalmsbden(ia2) = sum(iwsusc2(:,:,ia2))
    write(iodb,'("overlaps_susc2: ia=",i4," nalmsbden=",i4)') ia2, nalmsbden(ia2)
! ******
  end do
! ******
! Total number of basis functions
  ndensum = sum(nalmsbden); ndenmax = maxval(nalmsbden)
  write(*,'("overlaps_susc2: ndensum=",i8,"; ndenmax=",i8)') ndensum, ndenmax
! +++++++++++++++++ POINTERS +++++++++++++++++++++
! Tight upper bound for number of basis functions
  allocate(i2almsbden(4,ndensum))
  i = 0
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ilm=1,lmmax0
      do is=1,nsmax2
        do ib=1,iwsusc2(ilm,is,ia2)
          i = i + 1
          i2almsbden(:,i) = (/ib,ilm,is,ia2/)
        end do
      end do
    end do
  end do
  write(*,'("overlaps_susc2: ndensum should be ",i8)') i
! ++++++++++++++++++++++++++++++++++++++++++++++++
  tmp = 8.d0*ndenmax*ngfmax*nasusc2
  tmp = tmp/(1024.d0**2)
  write(*,'("overlaps_susc2: den Gaunt RAM=",f16.3," MB")') tmp
! Gaunt-like coefficients for density
  if (allocated(dengaunt)) deallocate(dengaunt)
  allocate(dengaunt(ngfmax,ndenmax,nasusc2))
  dengaunt = 0.d0
! potential in the density basis
  if (allocated(vlmsbden)) deallocate(vlmsbden)
  allocate(vlmsbden(ndensum))
  vlmsbden = 0.d0
! radial integral of each basis function
  if (allocated(suscnorm)) deallocate(suscnorm)
  allocate(suscnorm(ndensum))
  suscnorm = 0.d0
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    nr = nrpts(ia)
    dr(1:nr) = drmesh(1:nr,ia)
    nden0 = sum(nalmsbden(1:ia2-1))
    nden1 = sum(nalmsbden(1:ia2))
    do iq=1+nden0,nden1
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
      work(1:nr) = suscbasis(1:nr,ib,ilm,is,ia2)
      norm = real(radint(nr,work,dr,npanat(ia),ircutat(:,ia)))
      write(iodb,'("overlaps_susc2: ia,is,ilm,ib=",4i4," suscnorm=",es16.8)') ia, is, ilm, ib, real(norm)
      suscnorm(iq) = norm
    end do
  end do
  ngaunt = 0
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    nr = nrpts(ia)
    dr(1:nr) = drmesh(1:nr,ia)/rmesh(1:nr,ia)**2
    nden0 = sum(nalmsbden(1:ia2-1))
    nden1 = sum(nalmsbden(1:ia2))
    do jq=1+nden0,nden1
      i4 = i2almsbden(:,jq)
      ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
      ngf0 = sum(nalmsbgf(1:ia-1))
      ngf1 = sum(nalmsbgf(1:ia))
      do iq=1+ngf0,ngf1
        i3 = i2almsbgf(:,iq)
        q1 = i3(1); q2 = i3(2)!; ia = i3(3)
        i3 = i2lmsb(:,q2,ia)
        b2 = i3(1); lm2 = i3(2); s2 = i3(3)
        i3 = i2lmsb(:,q1,ia)
        b1 = i3(1); lm1 = i3(2); s1 = i3(3)
!       same spin combination
        if (is == is2i(s1,s2)) then
!       non-zero gaunt
        if (abs(rgaunt(lm1,lm2,ilm)) > ylmtol) then
!       ----------------------------------------------------------------
!         overlaps * Gaunt numbers
          ngaunt(ia2) = ngaunt(ia2) + 1
          work(1:nr) = phiref(1:nr,b1,i2lm(2,lm1),s1,ia)*phiref(1:nr,b2,i2lm(2,lm2),s2,ia)*suscbasis(1:nr,ib,ilm,is,ia2)
          dengaunt(iq-ngf0,jq-nden0,ia2) = rgaunt(lm1,lm2,ilm)*radint(nr,work,dr,npanat(ia),ircutat(:,ia))!/sqrt(fourpi)
!          write(*,'(2i4,2es16.8)') iq - ngf0, jq - nden0, dengaunt(iq-ngf0,jq-nden0,ia2)
!       ----------------------------------------------------------------
        end if
        end if
      end do
    end do
    write(*,'("overlaps_susc2: ia=",i4," non-zero den Gaunts=",i8)') ia, ngaunt(ia2)
  end do
! ************
! KS susc RAM
  ram = 16.d0*3*ndensum**2
  ram = ram/(1024.d0**3)
  write(*,'("overlaps_susc2: KS susc   RAM=",f16.3," GB")') ram
  if (ram > 4.d0) stop 'KS susc RAM > 4 GB'
! storage for KS susc in density basis
  if (allocated(kssusc0)) deallocate(kssusc0)
  if (allocated(kssusc1)) deallocate(kssusc1)
  if (allocated(kssusc2)) deallocate(kssusc2)
  allocate(kssusc0(ndensum,ndensum),kssusc1(ndensum,ndensum),kssusc2(ndensum,ndensum))
  kssusc0 = 0.d0; kssusc1 = 0.d0; kssusc2 = 0.d0
! Added by Sascha:
! storage for KS spin-current spin correlation function in mixed basis
  if (lcurrcorr) then
!   KS susc RAM
    ram = 16.d0*3*ndensum*ngradsum
    ram = ram/(1024.d0**3)
    write(*,'("overlaps_susc2: KS curr corr,   RAM=",f16.3," GB")') ram
    if (ram > 4.d0) stop 'KS curr corr RAM > 4 GB'
    if (allocated(kscurrcorr)) deallocate(kscurrcorr)
    if (allocated(kscurrcorr0)) deallocate(kscurrcorr0)
    if (allocated(kscurrcorr1)) deallocate(kscurrcorr1)
    if (allocated(kscurrcorr2)) deallocate(kscurrcorr2)
    allocate(kscurrcorr(ngradsum,ndensum),kscurrcorr0(ngradsum,ndensum),kscurrcorr1(ngradsum,ndensum),kscurrcorr2(ngradsum,ndensum))
    kscurrcorr = 0.d0; kscurrcorr0 = 0.d0; kscurrcorr1 = 0.d0; kscurrcorr2 = 0.d0
  end if
! kernel
  if (lkxc .or. lkha) then
    ram = 16.d0*ndensum**2
    ram = ram/(1024.d0**3)
    write(*,'("overlaps_susc2: kernel    RAM=",f16.3," GB")') ram
    if (ram > 4.d0) stop 'kernel RAM > 4 GB'
!   storage for kronecker KS
    if (allocated(kernel)) deallocate(kernel)
    allocate(kernel(ndensum,ndensum))
    kernel = 0.d0
  end if
! denominator
  if (lsusc .and. lenhanced) then
    ram = 16.d0*2.d0*ndensum**2
    ram = ram/(1024.d0**3)
    write(*,'("overlaps_susc2: susc den  RAM=",f16.3," GB")') ram
    if (ram > 4.d0) stop 'susc den RAM > 4 GB'
!   storage for kronecker KS
    if (allocated(denominator)) deallocate(denominator)
    if (allocated(kssusc)) deallocate(kssusc)
    allocate(denominator(ndensum,ndensum),kssusc(ndensum,ndensum))
    denominator = 0.d0; kssusc = 0.d0
  end if
  call cpu_time(finish)
  write(*,'("overlaps_susc2: time=",f10.3," s")') finish - start
! All done!
  end subroutine overlaps_susc2
