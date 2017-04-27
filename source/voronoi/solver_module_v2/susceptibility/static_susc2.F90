  subroutine static_susc2(onsite,struct)
! Static KS susceptibility
! Also first and second frequency derivatives
! Using density basis
  use global

  implicit none

  logical,           intent(in)  :: onsite, struct
! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0), ci= (0.d0,1.d0)
! -----------------------------------------------------------------
  complex(kind=c8b), allocatable :: gfijw(:,:), gfjiw(:,:), gfij0(:,:), gfji0(:,:), gf(:,:)
  complex(kind=c8b), allocatable :: gfijdef(:,:), gfjidef(:,:), gfijef(:,:), gfjief(:,:)
  complex(kind=c8b), allocatable :: gfijdeb(:,:), gfjideb(:,:), gfijeb(:,:), gfjieb(:,:)
  complex(kind=c8b), allocatable :: suscblock0(:,:), suscblock1(:,:), suscblock2(:,:), suscwork(:,:)
  complex(kind=c8b), allocatable :: currcorrblock0(:,:), currcorrblock1(:,:), currcorrblock2(:,:), currcorrwork(:,:)
  complex(kind=c8b), allocatable :: suscylm0(:,:,:,:,:,:), suscylm1(:,:,:,:,:,:), suscylm2(:,:,:,:,:,:), jijsusc(:,:,:,:)
  complex(kind=c8b), allocatable :: suscy00(:,:,:,:,:,:), dosylm(:,:,:), suscden(:,:,:), rhoden(:), rhodenylm(:,:,:)
  complex(kind=c8b), allocatable :: curr_corr_lm0(:,:,:,:,:,:,:),curr_corr_lm1(:,:,:,:,:,:,:),curr_corr_lm2(:,:,:,:,:,:,:)
  complex(kind=c8b) :: norm, de1, de2, e, norm2
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1, b1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2, b2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3, b3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4, b4
  integer(kind=i4b) :: ib, jb, is, js
  integer(kind=i4b) :: i2(2), i3(3), i4(4), ie, je, ia, ja, ia2, ja2, jlm, ilm, j, i, klm, ip, ie0, k, a
  integer(kind=i4b) :: iq, jq, iq0, iq1, jq0, jq1
  integer(kind=i4b) :: ipiv(ngfsum), info
  real(kind=r8b)    :: gaunti, gauntj, doublegaunt, maxelem0, maxelem1, maxelem2, start, finish, gftime, basistime, re, im
  real(kind=r8b)    :: maxelemcurr0, maxelemcurr1, maxelemcurr2
  real(kind=r8b)    :: dr(nrmax), rea, ima, reb, imb, ram
  complex(kind=c8b) :: work(nrmax), block(4,4)
  integer(kind=i4b) :: ne, nr
  complex(kind=c8b), external :: radint
  integer(kind=i4b) :: nepan, netot
  integer(kind=i4b) :: ipan(nescf)
  complex(kind=c8b) :: ea(nescf), eb(nescf), z(nescf), w(nescf), w0(nescf), w1(nescf)
  complex(kind=c8b) :: w00(nescf,nescf), w10(nescf,nescf), w01(nescf,nescf), w11(nescf,nescf) 
!  real(kind=r8b)    :: r(nalmsb), c(nalmsb), rwork(2*nalmsb), rcond, ferr(nalmsb), berr(nalmsb)
!  complex(kind=c8b) :: work(2*nalmsb), temp(nalmsb,nalmsb)
!  complex(kind=c8b) :: x(nalmsb,nalmsb)
!  character*1       :: equed
  character(len=1024) :: filename
  integer(kind=i4b)  :: orb22(1:2)
  complex(kind=c8b)  :: gfprint(1:4,1:4)
  gfprint = czero
  orb22(1)=5
  orb22(2)=9

  ram = 16*(4*ngfmax**2 + 13*nlmsb**2 + 16*(3*lmmax0**2+lmmax**2)*nasusc2**2 + (4*lmmax+6*lmmax0)*nasusc2 + ndensum)
  ram = ram/1024.d0**3
  write(*,'("RAM for static_susc2: ",f8.3," GB")') ram
  allocate(suscblock0(ngfmax,ngfmax),suscblock1(ngfmax,ngfmax),suscblock2(ngfmax,ngfmax),suscwork(ngfmax,ngfmax))
  allocate(gf(nlmsb,nlmsb),gfijw(nlmsb,nlmsb),gfjiw(nlmsb,nlmsb),gfij0(nlmsb,nlmsb),gfji0(nlmsb,nlmsb))
  allocate(gfijdef(nlmsb,nlmsb),gfjidef(nlmsb,nlmsb),gfijef(nlmsb,nlmsb),gfjief(nlmsb,nlmsb))
  allocate(gfijdeb(nlmsb,nlmsb),gfjideb(nlmsb,nlmsb),gfijeb(nlmsb,nlmsb),gfjieb(nlmsb,nlmsb))
  allocate(suscylm0(4,4,lmmax0,lmmax0,nasusc2,nasusc2),suscylm1(4,4,lmmax0,lmmax0,nasusc2,nasusc2))
  allocate(suscylm2(4,4,lmmax0,lmmax0,nasusc2,nasusc2),suscy00(4,4,lmmax,lmmax,nasusc2,nasusc2))
  allocate(jijsusc(4,4,nasusc2,nasusc2))
  allocate(dosylm(2,lmmax0,nasusc2),suscden(4,lmmax,nasusc2),rhoden(ndensum),rhodenylm(4,lmmax0,nasusc2))
  suscylm0 = 0.d0; suscylm1 = 0.d0; suscylm2 = 0.d0; suscy00 = 0.d0; suscden = 0.d0
  kssusc0 = 0.d0; kssusc1 = 0.d0; kssusc2 = 0.d0; dosylm = 0.d0
  gftime = 0.d0; basistime = 0.d0
  if(lcurrcorr) then
    allocate(currcorrblock0(ngfmax,ngfmax),currcorrblock1(ngfmax,ngfmax),currcorrblock2(ngfmax,ngfmax),currcorrwork(ngfmax,ngfmax))
    allocate(curr_corr_lm0(1:3,1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2))
    allocate(curr_corr_lm1(1:3,1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2))
    allocate(curr_corr_lm2(1:3,1:4,1:4,1:lmmax2,1:lmmax2,1:nasusc2,1:nasusc2))
    currcorrblock0 = 0.d0; currcorrblock1 = 0.d0; currcorrblock2 = 0.d0; currcorrwork = 0.d0; curr_corr_lm0 = 0.d0; curr_corr_lm1 = 0.d0; curr_corr_lm2 = 0.d0
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  open(file='meshpanels.dat',unit=iofile,status='old')
  read(iofile,*) nepan, netot
  if (netot /= nescf) stop 'static_susc2: netot /= nescf'
  do ie=1,nepan
    read(iofile,*) rea, ima, reb, imb, ipan(ie)
    write(*,'(4f10.6,i8)') rea, ima, reb, imb, ipan(ie)
    ea(ie) = cmplx(rea,ima); eb(ie) = cmplx(reb,imb)
  end do
  if (netot /= sum(ipan(1:nepan))) stop 'static_susc2: netot /= sum(ipan)'
  close(iofile)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Added by Sascha:
  if(lcurrcorr) then 
    write(filename,'("currentspdf.dat")')
    open(file=filename,unit=iofile)
    write(iofile,'("#Current correlation function calculation results")')
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ja2=1,nasusc2      ! atom j
    ja = iasusc2(ja2)
    jq0 = sum(nalmsbgf(1:ja-1))
    jq1 = sum(nalmsbgf(1:ja))
    do ia2=1,nasusc2    ! atom i
      ia = iasusc2(ia2)
      iq0 = sum(nalmsbgf(1:ia-1))
      iq1 = sum(nalmsbgf(1:ia))
      suscblock0 = czero; suscblock1 = czero; suscblock2 = czero
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call cpu_time(start)
!     Integration on panels
      ie0 = 0
!     Extrapolation to EB
      gfijeb = czero; gfjieb = czero; gfijdeb = czero; gfjideb = czero
!     Extrapolation to EF
      gfijef = czero; gfjief = czero; gfijdef = czero; gfjidef = czero
!     *************
      do ip=1,nepan
!     *************
        call chebyweights(ea(ip),eb(ip),ipan(ip),nescf,z,w,w0,w1,w00,w10,w01,w11)
!       Go over energies in each panel
!       ~~~~~~~~~~~~~~~~
        do je=1,ipan(ip)     ! energy
!       ~~~~~~~~~~~~~~~~
!         --------------------------------------------------------------
!         Get the projected GFs
!         Gji(E + i0)
          call projected_gf(ie0+je,ja,ia,gfji0,lsusconsite,lsuscstruct)
          if (.true.) call local_frame(ja,ia,magdir,gfji0,gf)   !!!!!!!!!!!!!!!!!! this only needs to be called if lrot or soc from host!!!!!!!!!!!!!!!!!!!!!
!         Gji(E - i0) = Gij(E + i0)^\dagger
          call projected_gf(ie0+je,ia,ja,gfjiw,lsusconsite,lsuscstruct)
          if (.true.) call local_frame(ia,ja,magdir,gfjiw,gf)  !!!! same here !!!!!!!!!!!!!!!!!
          gfjiw = conjg(transpose(gfjiw))
!         --------------------------------------------------------------
!         ~~~~~~~~~~~~~~~~
          do ie=1,ipan(ip)     ! energy
!         ~~~~~~~~~~~~~~~~
!          ie = je
!            write(*,*) "ie, je", ie, je
!           Get the projected GFs
!           ------------------------------------------------------------
!           Gij(E + i0)
            call projected_gf(ie0+ie,ia,ja,gfijw,lsusconsite,lsuscstruct)
            if (.true.) call local_frame(ia,ja,magdir,gfijw,gf)
!           Gij(E - i0) = Gji(E + i0)^\dagger
            call projected_gf(ie0+ie,ja,ia,gfij0,lsusconsite,lsuscstruct)
            if (.true.) call local_frame(ja,ia,magdir,gfij0,gf)
            gfij0 = conjg(transpose(gfij0))
!           ------------------------------------------------------------
            if (ie == je) then
!             Extrapolation to EB + i0
              if (ip == 1) then
!               Gij(EB + i0)
                gfijeb  = gfijeb  + w0(ipan(1)-ie+1)*gfijw
!               Gji(EB + i0)
                gfjieb  = gfjieb  + w0(ipan(1)-ie+1)*gfji0
!               d/dE Gij(EB + i0)
                gfijdeb = gfijdeb - w1(ipan(1)-ie+1)*gfijw
!               d/dE Gji(EB + i0)
                gfjideb = gfjideb - w1(ipan(1)-ie+1)*gfji0
              end if
!             Extrapolation to EF + i0
              if (ip == nepan) then
!               Gij(EF + i0)
                gfijef  = gfijef  + w0(ie)*gfijw
!               Gji(EF + i0)
                gfjief  = gfjief  + w0(ie)*gfji0
!               d/dE Gij(EF + i0)
                gfijdef = gfijdef + w1(ie)*gfijw
!               d/dE Gji(EF + i0)
                gfjidef = gfjidef + w1(ie)*gfji0
              end if
            end if
!           Static susc + frequency derivatives
            do j=1,jq1-jq0
              jq = jq0 + j
              i3 = i2almsbgf(:,jq)
              q3 = i3(1); q4 = i3(2)!; ja = i3(3)
              i3 = i2lmsb(:,q3,ia)
              b3 = i3(1); lm3 = i3(2)
              i3 = i2lmsb(:,q4,ia)
              b4 = i3(1); lm4 = i3(2)
              do i=1,iq1-iq0
                iq = iq0 + i
                i3 = i2almsbgf(:,iq)
                q1 = i3(1); q2 = i3(2)!; ia = i3(3)
                i3 = i2lmsb(:,q1,ia)
                b1 = i3(1); lm1 = i3(2)
                i3 = i2lmsb(:,q2,ia)
                b2 = i3(1); lm2 = i3(2)
!               +1/i2pi \int dz^* Gji(z^*)^\dagger Gij(z^*)^\dagger - 1/i2pi \int dz Gij(z) Gji(z)
                norm = (conjg(w00(ie,je))*gfij0(q2,q3)*gfjiw(q4,q1) - w00(ie,je)*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
                suscblock0(i,j) = suscblock0(i,j) + norm
!               -1/i2pi \int dz^* Gji(z^*)^\dagger d/dz^* Gij(z^*)^\dagger - 1/i2pi \int dz d/dz Gij(z) Gji(z)
                norm = -(conjg(w01(ie,je))*gfij0(q2,q3)*gfjiw(q4,q1) + w10(ie,je)*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
                suscblock1(i,j) = suscblock1(i,j) + norm
!               -1/i2pi \int dz^* d/dz^* Gji(z^*)^\dagger d/dz^* Gij(z^*)^\dagger + 1/i2pi \int dz d/dz Gij(z) d/dz Gji(z)
                norm = -(conjg(w11(ie,je))*gfij0(q2,q3)*gfjiw(q4,q1) - w11(ie,je)*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
                suscblock2(i,j) = suscblock2(i,j) + norm
              end do
            end do
!         ~~~~~~
          end do            ! energy
!         ~~~~~~
!       ~~~~~~
        end do            ! energy
!       ~~~~~~
        ie0 = ie0 + ipan(ip)
!     ******
      end do
!     ******
!     Contributions from EF and EB
      do j=1,jq1-jq0
        jq = jq0 + j
        i3 = i2almsbgf(:,jq)
        q3 = i3(1); q4 = i3(2)!; ja = i3(3)
        i3 = i2lmsb(:,q3,ia)
        b3 = i3(1); lm3 = i3(2)
        i3 = i2lmsb(:,q4,ia)
        b4 = i3(1); lm4 = i3(2)
        do i=1,iq1-iq0
          iq = iq0 + i
          i3 = i2almsbgf(:,iq)
          q1 = i3(1); q2 = i3(2)!; ia = i3(3)
          i3 = i2lmsb(:,q1,ia)
          b1 = i3(1); lm1 = i3(2)
          i3 = i2lmsb(:,q2,ia)
          b2 = i3(1); lm2 = i3(2)
!         Nothing for static susc
          norm = czero
          suscblock0(i,j) = suscblock0(i,j) + norm
!         1/i2pi Gij(EF+i0) Gij(EF+i0)^\dagger - 1/i2pi Gij(EB+i0) Gij(EB+i0)^\dagger
          norm = czero
          norm = norm + (gfijef(q2,q3)*conjg(gfijef(q1,q4)) - gfijeb(q2,q3)*conjg(gfijeb(q1,q4)))/i2pi
          suscblock1(i,j) = suscblock1(i,j) + norm
!         -1/i2pi d/dE Gij(EF+i0) * (Gji(EF+i0) - Gij(EF+i0)^\dagger)
!         +1/i2pi d/dE Gij(EB+i0) * (Gji(EB+i0) - Gij(EB+i0)^\dagger)
!         -1/i2pi (Gij(EF+i0) - Gji(EF+i0)^\dagger) * d/dE Gij(EF+i0)^\dagger
!         +1/i2pi (Gij(EB+i0) - Gji(EB+i0)^\dagger) * d/dE Gij(EB+i0)^\dagger
          norm = czero
          norm = norm - gfijdef(q2,q3)*(gfjief(q4,q1) - conjg(gfijef(q1,q4)))/i2pi
          norm = norm + gfijdeb(q2,q3)*(gfjieb(q4,q1) - conjg(gfijeb(q1,q4)))/i2pi
          norm = norm - (gfijef(q2,q3) - conjg(gfjief(q3,q2)))*conjg(gfijdef(q1,q4))/i2pi
          norm = norm + (gfijeb(q2,q3) - conjg(gfjieb(q3,q2)))*conjg(gfijdeb(q1,q4))/i2pi
          suscblock2(i,j) = suscblock2(i,j) + norm
        end do
      end do
!     ******
!  Print Gf elements
!      do b1 = 1,2
!        do b2 = 1,2
!          do i = 1,2
!            do j =1,2
!              lm1 = orb22(i)
!              lm2 = orb22(j)
!              q1 = lmsb2i(b1,lm1,1,ia)
!              q2 = lmsb2i(b2,lm2,1,ia)
!              gfprint(i+2*b1-2,j+2*b2-2)=gfijef(q1,q2)             
!            end do
!          end do
!        end do
!      end do
!      write(iofile,'("Gfef for lm=5,9 and b=1,2")')
!      write(iofile,'("lm=5 b=1 , lm=9 b=1, lm=5 b=2, lm=9 b=2")')
!      do i=1,4
!        write(iofile,'(8es16.8)') (gfprint(i,j),j=1,4)
!      end do
!     ******
!     Added by Sascha:
!     KS spin-current spin correlation function
!     Difference of susceptibility with swapped basis in the first argument ( L_1 b1 <-> L_4 b_4 )
      if(lcurrcorr) then
!      write(*,'("KS curr corr in gf basis")')
      do j=1,jq1-jq0
        jq = jq0 + j
        i3 = i2almsbgf(:,jq)
        q3 = i3(1); q4 = i3(2)!; ja = i3(3)
        i3 = i2lmsb(:,q3,ia)
        lm3 = i3(2)
        i3 = i2lmsb(:,q4,ia)
        lm4 = i3(2)
        do i=1,iq1-iq0
          iq = iq0 + i
          i3 = i2almsbgf(:,iq)
          q1 = i3(1); q2 = i3(2)!; ia = i3(3)
          i3 = i2lmsb(:,q1,ia)
          b1 = i3(1); lm1 = i3(2); s1 = i3(3)
          i3 = i2lmsb(:,q2,ia)
          b2 = i3(1); lm2 = i3(2); s2 = i3(3)
          q1 = lmsb2i(b1,lm1,s2,ia)
          q2 = lmsb2i(b2,lm2,s1,ia)
!         Swapped indices in the first argument of the susceptibility in the gfbasis
          k = almsbgf2i(q2,q1,ia) - iq0
          currcorrblock0(i,j) =  ci*( suscblock0(i,j) - suscblock0(k,j) )
          currcorrblock1(i,j) =  ci*( suscblock1(i,j) - suscblock1(k,j) )
          currcorrblock2(i,j) =  ci*( suscblock2(i,j) - suscblock2(k,j) )
!          if(abs(currcorrblock0(i,j)) > ylmtol) write(*,'(3i6,2e16.8)') 0,i,j, currcorrblock0(i,j)
!          if(abs(currcorrblock1(i,j)) > ylmtol) write(*,'(3i6,2e16.8)') 1,i,j, currcorrblock1(i,j)
!          if(abs(currcorrblock2(i,j)) > ylmtol) write(*,'(3i6,2e16.8)') 2,i,j, currcorrblock2(i,j)
        end do
      end do
!      if (abs(sum(currcorrblock0)) < ylmtol) write(*,'("KS curr corr0 in GF basis wrong!",2e16.8)') abs(sum(currcorrblock0))
!      if (abs(sum(currcorrblock1)) < ylmtol) write(*,'("KS curr corr1 in GF basis wrong!",2e16.8)') abs(sum(currcorrblock1))
!      if (abs(sum(currcorrblock2)) < ylmtol) write(*,'("KS curr corr2 in GF basis wrong!",2e16.8)') abs(sum(currcorrblock2))
      end if
!     ------------------------------------------------------------------
!     DOS at Fermi energy
      if (ia2 == ja2) then
        do klm=1,lmmax0
          do jq=1,nlmsba(ia)
            i3 = i2lmsb(:,jq,ia)
            jb = i3(1); jlm = i3(2); js = i3(3)
            do iq=1,nlmsba(ia)
              i3 = i2lmsb(:,iq,ia)
              ib = i3(1); ilm = i3(2); is = i3(3)
              if (abs(rgaunt(ilm,jlm,klm)) > ylmtol) then
                if (is == 1 .and. js == 1) dosylm(2,klm,ia2) = dosylm(2,klm,ia2) &
                                                               + rgaunt(ilm,jlm,klm)*overlap(iq,jq,ia)*(gfijef(iq,jq) - conjg(gfjief(jq,iq)))/i2pi
                if (is == 2 .and. js == 2) dosylm(1,klm,ia2) = dosylm(1,klm,ia2) &
                                                               + rgaunt(ilm,jlm,klm)*overlap(iq,jq,ia)*(gfijef(iq,jq) - conjg(gfjief(jq,iq)))/i2pi
              end if
            end do
          end do
        end do
      end if
!     ------------------------------------------------------------------
      call cpu_time(finish)
      gftime = gftime + finish - start
      call cpu_time(start)
      call add_susc_block(ia2,ja2,suscblock0,suscwork,kssusc0)
      call add_susc_block(ia2,ja2,suscblock1,suscwork,kssusc1)
      call add_susc_block(ia2,ja2,suscblock2,suscwork,kssusc2)
!     Added by Sascha:
!     KS spin-current spin correlation function in density basis only in the second argument r' 
      if(lcurrcorr) then
        call add_susc_block_one_sided(ia2,ja2,currcorrblock0,currcorrwork,kscurrcorr0)
        call add_susc_block_one_sided(ia2,ja2,currcorrblock1,currcorrwork,kscurrcorr1)
        call add_susc_block_one_sided(ia2,ja2,currcorrblock2,currcorrwork,kscurrcorr2)
!        if (abs(sum(kscurrcorr0)) < ylmtol) write(*,'("KS curr corr0 in den basis wrong!",2e16.8)') abs(sum(kscurrcorr0))
!        if (abs(sum(kscurrcorr1)) < ylmtol) write(*,'("KS curr corr1 in den basis wrong!",2e16.8)') abs(sum(kscurrcorr1))
!        if (abs(sum(kscurrcorr2)) < ylmtol) write(*,'("KS curr corr2 in den basis wrong!",2e16.8)') abs(sum(kscurrcorr2))
      end if
      call cpu_time(finish)
      basistime = basistime + finish - start
!      maxelem = maxval(abs(kssusc0))
!      where (abs(kssusc0) < susctol*maxelem) kssusc0 = 0.d0
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end do            ! atom i
  end do              ! atom j
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  do jq=1,ndensum
!    do iq=1,ndensum
!      re = real(kssusc0(iq,jq)); im = aimag(kssusc0(iq,jq))
!      if (abs(re) < susctol) re = 0.d0
!      if (abs(im) < susctol) im = 0.d0
!      kssusc0(iq,jq) = cmplx(re,im)
!    end do
!  end do
!  where (abs(kssusc0) < susctol) kssusc0 = 0.d0
! ------------------------------------------------------------------
!  if (lenhanced) then
!    call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc0,ngfsum,kernel,ngfsum,czero,denominator,ngfsum)
!    do iq=1,ndensum
!      denominator(iq,iq) = denominator(iq,iq) + 1.d0
!    end do
!    call zgesv(ndensum,ndensum,denominator,ngfsum,ipiv,kssusc0,ngfsum,info)
!    if (info /= 0) stop 'dyn_susc_real2: failure in zgesv'
!    call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc0,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
!    write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
!    write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
!    if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
!    kssusc0 = x
!  end if
! ------------------------------------------------------------------
! Symmetrize
!  call symmetrize(nalmsb,kssusc0,susctol)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! multipoles
  suscylm0 = 0.d0; suscylm1 = 0.d0; suscylm2 = 0.d0
  do jq=1,ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
    do iq=1,ndensum
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
      suscylm0(is,js,ilm,jlm,ia2,ja2) = suscylm0(is,js,ilm,jlm,ia2,ja2) + suscnorm(iq)*kssusc0(iq,jq)*suscnorm(jq)
      suscylm1(is,js,ilm,jlm,ia2,ja2) = suscylm1(is,js,ilm,jlm,ia2,ja2) + suscnorm(iq)*kssusc1(iq,jq)*suscnorm(jq)
      suscylm2(is,js,ilm,jlm,ia2,ja2) = suscylm2(is,js,ilm,jlm,ia2,ja2) + suscnorm(iq)*kssusc2(iq,jq)*suscnorm(jq)
      if (ilm == 1 .and. jlm == 1) then
        jijsusc(is,js,ia2,ja2) = jijsusc(is,js,ia2,ja2) + suscnorm(iq)*kssusc0(iq,jq)*suscnorm(jq)
      end if
    end do
  end do
! Added by Sascha:
! multipoles
  if(lcurrcorr) then
  curr_corr_lm0 = 0.d0; curr_corr_lm1 = 0.d0; curr_corr_lm2 = 0.d0
  do jq=1,ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
    if(jlm < lmmax2 +1) then
    do iq = 1,ngradsum
      i3 = i2almsbgrad(:,iq)
      q1 = i3(1); q2 = i3(2); ia2 = i3(3); ia = iasusc2(ia2)
      i3 = i2lmsb(:,q1,ia)
      s1 = i3(3)  !check spin indices once again!!!
      i3 = i2lmsb(:,q2,ia)
      s2 = i3(3)  !check spin indices once again!!!
      is = is2i(s1,s2)
      do ilm = 1,lmmax2
        do a= 1,3
          curr_corr_lm0(a,is,js,ilm,jlm,ia2,ja2) = curr_corr_lm0(a,is,js,ilm,jlm,ia2,ja2) + gradnorm(a,ilm,iq)*kscurrcorr0(iq,jq)*suscnorm(jq)
          curr_corr_lm1(a,is,js,ilm,jlm,ia2,ja2) = curr_corr_lm1(a,is,js,ilm,jlm,ia2,ja2) + gradnorm(a,ilm,iq)*kscurrcorr1(iq,jq)*suscnorm(jq)
          curr_corr_lm2(a,is,js,ilm,jlm,ia2,ja2) = curr_corr_lm2(a,is,js,ilm,jlm,ia2,ja2) + gradnorm(a,ilm,iq)*kscurrcorr2(iq,jq)*suscnorm(jq)
        end do
      end do
    end do
    end if
  end do
!  if (abs(sum(curr_corr_lm0)) < ylmtol) write(*,'("KS curr_corr_lm0 wrong!",2e16.8)') abs(sum(curr_corr_lm0))
!  if (abs(sum(curr_corr_lm1)) < ylmtol) write(*,'("KS curr_corr_lm1 wrong!",2e16.8)') abs(sum(curr_corr_lm1))
!  if (abs(sum(curr_corr_lm2)) < ylmtol) write(*,'("KS curr_corr_lm2 wrong!",2e16.8)') abs(sum(curr_corr_lm2))
  end if
! ----------------------------------------------------------------------
! Final filtering (global)
  maxelem0 = maxval(abs(suscylm0))
  where (abs(suscylm0) < susctol*maxelem0) suscylm0 = 0.d0
  maxelem1 = maxval(abs(suscylm1))
  where (abs(suscylm1) < susctol*maxelem1) suscylm1 = 0.d0
  maxelem2 = maxval(abs(suscylm2))
  where (abs(suscylm2) < susctol*maxelem2) suscylm2 = 0.d0
! Added by Sascha
! Final filtering (global)
  if(lcurrcorr) then
    maxelemcurr0 = maxval(abs(curr_corr_lm0))
    where (abs(curr_corr_lm0) < susctol*maxelemcurr0) curr_corr_lm0 = 0.d0
    maxelemcurr1 = maxval(abs(curr_corr_lm1))
    where (abs(curr_corr_lm1) < susctol*maxelemcurr1) curr_corr_lm1 = 0.d0
    maxelemcurr2 = maxval(abs(curr_corr_lm2))
    where (abs(curr_corr_lm2) < susctol*maxelemcurr2) curr_corr_lm2 = 0.d0
  end if
!  maxelem = maxval(abs(suscy00))
!  where (abs(suscy00) < susctol*maxelem) suscy00 = 0.d0
! ----------------------------------------------------------------------
! DOS at EF test
  where (abs(dosylm) < susctol) dosylm = 0.d0
  write(*,'(/,"DOS at EF sum rule: ia, dos, susc (up then dn)")')
  do ia2=1,nasusc2
    do ilm=1,lmmax0
      write(*,'(2i4,8es16.8)') ia2, ilm, real(dosylm(1,ilm,ia2)), real(sum(suscylm0(3,3:4,ilm,1,ia2,:))),  &
          real(dosylm(2,ilm,ia2)), real(sum(suscylm0(4,3:4,ilm,1,ia2,:)))
!      write(*,'(2i4,8es16.8)') ia2, ilm, dosylm(1,ilm,ia2), sum(suscylm(3,3:4,ilm,1,ia2,:)),  &
!          dosylm(2,ilm,ia2), sum(suscylm(4,3:4,ilm,1,ia2,:))
    end do
  end do
! ----------------------------------------------------------------------
! change from spin to cartesian
  if (lcartesian) then
    do ja2=1,nasusc2
    do ia2=1,nasusc2
      do jlm=1,lmmax0
      do ilm=1,lmmax0
        block(:,:) = czero
        do j=1,4
        do i=1,4
          do js=1,4
          do is=1,4
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*suscylm0(is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
          end do
          end do
        end do
        end do
        suscylm0(:,:,ilm,jlm,ia2,ja2) = block
        block(:,:) = czero
        do j=1,4
        do i=1,4
          do js=1,4
          do is=1,4
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*suscylm1(is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
          end do
          end do
        end do
        end do
        suscylm1(:,:,ilm,jlm,ia2,ja2) = block
        block(:,:) = czero
        do j=1,4
        do i=1,4
          do js=1,4
          do is=1,4
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*suscylm2(is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
          end do
          end do
        end do
        end do
        suscylm2(:,:,ilm,jlm,ia2,ja2) = block
      end do
      end do
    end do
    end do
  end if

! Added by Sascha:
! change from spin to cartesian
  if (lcartesian .and. lcurrcorr) then
    do ja2=1,nasusc2
    do ia2=1,nasusc2
      do jlm=1,lmmax2
      do ilm=1,lmmax2
      do a = 1,3
        block(:,:) = czero
        do j=1,4
        do i=1,4
          do js=1,4
          do is=1,4
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*curr_corr_lm0(a,is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
          end do
          end do
        end do
        end do
        curr_corr_lm0(a,:,:,ilm,jlm,ia2,ja2) = block
        block(:,:) = czero
        do j=1,4
        do i=1,4
          do js=1,4
          do is=1,4
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*curr_corr_lm1(a,is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
          end do
          end do
        end do
        end do
        curr_corr_lm1(a,:,:,ilm,jlm,ia2,ja2) = block
        block(:,:) = czero
        do j=1,4
        do i=1,4
          do js=1,4
          do is=1,4
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*curr_corr_lm2(a,is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
          end do
          end do
        end do
        end do
        curr_corr_lm2(a,:,:,ilm,jlm,ia2,ja2) = block
      end do
      end do
      end do
    end do
    end do
!  if (abs(sum(curr_corr_lm0)) < ylmtol) write(*,'("KS curr_corr_lm0 in cartesian coodinates wrong!")')
!  if (abs(sum(curr_corr_lm1)) < ylmtol) write(*,'("KS curr_corr_lm1 in cartesian coodinates wrong!")')
!  if (abs(sum(curr_corr_lm2)) < ylmtol) write(*,'("KS curr_corr_lm2 in cartesian coodinates wrong!")')
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*,'(/,"KS susc 0 ia, ja, im, il, jm, jl, susc ylm")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm=1,1!lmmax0
        do ilm=1,1!lmmax0
          norm = 0.d0
          do j=1,4
            do i=1,4
              re = real(suscylm0(i,j,ilm,jlm,ia2,ja2))
              if (abs(re) < maxelem0*susctol) re = 0.d0
              im = aimag(suscylm0(i,j,ilm,jlm,ia2,ja2))
              if (abs(im) < maxelem0*susctol) im = 0.d0
              suscylm0(i,j,ilm,jlm,ia2,ja2) = cmplx(re,im)
              norm = max(abs(norm),abs(re)+abs(im))
            end do
          end do
!          norm = sum(abs(suscylm0(:,:,ilm,jlm,ia2,ja2)))
          if (abs(norm) > maxelem0*susctol) then
            do i=1,4
              write(*,'(6i4,8es16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscylm0(i,j,ilm,jlm,ia2,ja2),j=1,4)
            end do
          end if
        end do
      end do
    end do
  end do
! ----------------------------------------------------------------------
  write(*,'(/,"KS susc 1 ia, ja, im, il, jm, jl, susc ylm")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm=1,1!lmmax0
        do ilm=1,1!lmmax0
          norm = 0.d0
          do j=1,4
            do i=1,4
              re = real(suscylm1(i,j,ilm,jlm,ia2,ja2))
              if (abs(re) < maxelem1*susctol) re = 0.d0
              im = aimag(suscylm1(i,j,ilm,jlm,ia2,ja2))
              if (abs(im) < maxelem1*susctol) im = 0.d0
              suscylm1(i,j,ilm,jlm,ia2,ja2) = cmplx(re,im)
              norm = max(abs(norm),abs(re)+abs(im))
            end do
          end do
!          norm = sum(abs(suscylm1(:,:,ilm,jlm,ia2,ja2)))
          if (abs(norm) > maxelem1*susctol) then
            do i=1,4
              write(*,'(6i4,8es16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscylm1(i,j,ilm,jlm,ia2,ja2),j=1,4)
            end do
          end if
        end do
      end do
    end do
  end do
! ----------------------------------------------------------------------
  write(*,'(/,"KS susc 2 ia, ja, im, il, jm, jl, susc ylm")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm=1,1!lmmax0
        do ilm=1,1!lmmax0
          norm = 0.d0
          do j=1,4
            do i=1,4
              re = real(suscylm0(i,j,ilm,jlm,ia2,ja2))
              if (abs(re) < maxelem0*susctol) re = 0.d0
              im = aimag(suscylm0(i,j,ilm,jlm,ia2,ja2))
              if (abs(im) < maxelem0*susctol) im = 0.d0
              suscylm0(i,j,ilm,jlm,ia2,ja2) = cmplx(re,im)
              norm = max(abs(norm),abs(re)+abs(im))
            end do
          end do
!          norm = sum(abs(suscylm0(:,:,ilm,jlm,ia2,ja2)))
          if (abs(norm) > maxelem2*susctol) then
            do i=1,4
              write(*,'(6i4,8es16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscylm2(i,j,ilm,jlm,ia2,ja2),j=1,4)
            end do
          end if
        end do
      end do
    end do
  end do
! Added by Sascha:
! write ks curr corr
  if(lcurrcorr) then
  write(*,'(/,"KS curr corr 0 ia, ja, im, il, jm, jl, direction, curr corr lm")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      write(iofile,'("Atom i=",i4," j=",i4)') ia,ja
      write(iofile,'("static KScurr im, il, direction, current")')
      do jlm=1,1!lmmax0
        do ilm=1,lmmax0!lmmax0
          do a = 1,3
            norm = 0.d0
            do j=1,4
              do i=1,4
                re = real(curr_corr_lm0(a,i,j,ilm,jlm,ia2,ja2))
                if (abs(re) < maxelemcurr0*susctol) re = 0.d0
                im = aimag(curr_corr_lm0(a,i,j,ilm,jlm,ia2,ja2))
                if (abs(im) < maxelemcurr0*susctol) im = 0.d0
                curr_corr_lm0(a,i,j,ilm,jlm,ia2,ja2) = cmplx(re,im)
                norm = max(abs(norm),abs(re)+abs(im))
              end do
            end do
            if (abs(norm) > maxelemcurr0*susctol) then
              do i=1,4
                write(*,'(7i4,8es16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), a, (curr_corr_lm0(a,i,j,ilm,jlm,ia2,ja2),j=1,4)
                write(iofile,'(3i4,8es16.8)') i2lm(:,ilm), a, (curr_corr_lm0(a,i,j,ilm,jlm,ia2,ja2),j=1,4)
              end do
            end if
          end do
        end do
      end do
    end do
  end do
  write(*,'(/,"KS curr corr 1 ia, ja, im, il, jm, jl, direction, curr corr lm")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm=1,1!lmmax0
        do ilm=1,lmmax0
          do a = 1,3
            norm = 0.d0
            do j=1,4
              do i=1,4
                re = real(curr_corr_lm1(a,i,j,ilm,jlm,ia2,ja2))
                if (abs(re) < maxelemcurr1*susctol) re = 0.d0
                im = aimag(curr_corr_lm1(a,i,j,ilm,jlm,ia2,ja2))
                if (abs(im) < maxelemcurr1*susctol) im = 0.d0
                curr_corr_lm1(a,i,j,ilm,jlm,ia2,ja2) = cmplx(re,im)
                norm = max(abs(norm),abs(re)+abs(im))
              end do
            end do
            if (abs(norm) > maxelemcurr1*susctol) then
              do i=1,4
                write(*,'(7i4,8es16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), a, (curr_corr_lm1(a,i,j,ilm,jlm,ia2,ja2),j=1,4)
              end do
            end if
          end do
        end do
      end do
    end do
  end do
  write(*,'(/,"KS curr corr 2 ia, ja, im, il, jm, jl, direction, curr corr lm")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm=1,1!lmmax0
        do ilm=1,lmmax0
          do a = 1,3
            norm = 0.d0
            do j=1,4
              do i=1,4
                re = real(curr_corr_lm2(a,i,j,ilm,jlm,ia2,ja2))
                if (abs(re) < maxelemcurr2*susctol) re = 0.d0
                im = aimag(curr_corr_lm2(a,i,j,ilm,jlm,ia2,ja2))
                if (abs(im) < maxelemcurr2*susctol) im = 0.d0
                curr_corr_lm2(a,i,j,ilm,jlm,ia2,ja2) = cmplx(re,im)
                norm = max(abs(norm),abs(re)+abs(im))
              end do
            end do
            if (abs(norm) > maxelemcurr2*susctol) then
              do i=1,4
                write(*,'(7i4,8es16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), a, (curr_corr_lm2(a,i,j,ilm,jlm,ia2,ja2),j=1,4)
              end do
            end if
          end do
        end do
      end do
    end do
  end do
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    norm2 = czero
!   Now construct the potential in the susceptibility basis
    nr = nrpts(ia)
    dr(1:nr) = drmesh(1:nr,ia)!/rmesh(1:nr,ia)**2
    do iq=1+sum(nalmsbden(1:ia2-1)),sum(nalmsbden(1:ia2))
      vlmsbden(iq) = 0.d0
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
      s1 = i2is(1,is); s2 = i2is(2,is)
!     spherical input potential
!      if (ilm == 1 .and. is == 1) then
      if (ilm == 1) then
        work(1:nr) = suscbasis(1:nr,ib,ilm,is,ia2)*br(1:nr,ia)
        norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
        vlmsbden(iq) = norm*2.d0
        norm2 = norm2 + vlmsbden(iq)*suscnorm(iq)
      end if
    end do
    write(*,'("KS susc ia=",i4," potnorm=",2f16.8)') ia, norm2
  end do
! ----------------------------------------------------------------------
! Jij from Bxc * KS susc * Bxc
  jijsusc = 0.d0
  do jq=1,ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
    do iq=1,ndensum
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
      jijsusc(is,js,ia2,ja2) = jijsusc(is,js,ia2,ja2) - vlmsbden(iq)*kssusc0(iq,jq)*vlmsbden(jq)
    end do
  end do
  write(*,'(/,"Jijsusc 0 ia, ja, Jij")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      if (lcartesian) then
        block(:,:) = czero
        do j=1,4
        do i=1,4
          do js=1,4
          do is=1,4
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*jijsusc(is,js,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
          end do
          end do
        end do
        end do
        jijsusc(:,:,ia2,ja2) = block
      end if
      do j=1,4
        do i=1,4
          re = real(jijsusc(i,j,ia2,ja2))
          if (abs(re) < 1.d-8) re = 0.d0
          im = aimag(jijsusc(i,j,ia2,ja2))
          if (abs(im) < 1.d-8) im = 0.d0
          jijsusc(i,j,ia2,ja2) = cmplx(re,im)
        end do
      end do
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (jijsusc(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
! ----------------------------------------------------------------------
! KS susc * Bxc in density basis
  call zgemv('N',ndensum,ndensum,cone,kssusc0,ndensum,vlmsbden,1,czero,rhoden,1)
! multipoles
  rhodenylm = 0.d0
  do iq=1,ndensum
    i4 = i2almsbden(:,iq)
    ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
    rhodenylm(is,ilm,ia2) = rhodenylm(is,ilm,ia2) + suscnorm(iq)*rhoden(iq)
  end do
  write(*,'(/,"KS susc ia, im, il, rhoden ylm")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ilm=1,lmmax0
      norm = sum(abs(rhodenylm(:,ilm,ia2)))
      if (abs(norm) > susctol) then
        write(*,'(3i4,8f16.8)') ia, i2lm(:,ilm), rhodenylm(:,ilm,ia2)
      end if
    end do
  end do
!  if (lenhanced) denominator = kssusc0
!  call zgesv(ndensum,1,kssusc0,ndensum,ipiv,rhoden,ndensum,info)
!  if (lenhanced) kssusc0 = denominator
!  rhodenylm = 0.d0
!  do iq=1,ndensum
!    i4 = i2almsbden(:,iq)
!    ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
!    rhodenylm(is,ilm,ia2) = rhodenylm(is,ilm,ia2) + suscnorm(iq)*rhoden(iq)
!  end do
!  write(*,'(/,"KS susc ia, im, il, potden ylm")')
!  do ia2=1,nasusc2
!    ia = iasusc2(ia2)
!    do ilm=1,lmmax0
!      norm = sum(abs(rhodenylm(:,ilm,ia2)))
!      if (abs(norm) > susctol) then
!        write(*,'(3i4,8f16.8)') ia, i2lm(:,ilm), rhodenylm(:,ilm,ia2)
!      end if
!    end do
!  end do
! ----------------------------------------------------------------------
  deallocate(suscblock0,suscblock1,suscblock2,suscwork)
  deallocate(gf,gfijw,gfjiw,gfij0,gfji0,gfijdef,gfjidef,gfijef,gfjief,gfijdeb,gfjideb,gfijeb,gfjieb)
  deallocate(suscylm0,suscylm1,suscylm2,jijsusc,suscy00,dosylm,suscden,rhoden,rhodenylm)
  if(lcurrcorr) deallocate(currcorrblock0,currcorrblock1,currcorrblock2,currcorrwork,curr_corr_lm0,curr_corr_lm1,curr_corr_lm2)
  write(*,'(/," Static  KS susc chiGF time=",f10.3," s")') gftime
  write(*,'(" Static  KS susc basis time=",f10.3," s")') basistime
  close(iofile)
! All done!
  end subroutine static_susc2
