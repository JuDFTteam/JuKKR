  subroutine static_susc_rsa(onsite,struct,chiks0,chiks1,chiks2)
! Static KS susceptibility
! Also first and second frequency derivatives
! Using density basis
  use global

  implicit none

  logical,           intent(in)  :: onsite, struct
! static KS susceptibility
  complex(kind=c8b), intent(out) :: chiks0(4,4,nasusc2,nasusc2)
! first frequency derivative
  complex(kind=c8b), intent(out) :: chiks1(4,4,nasusc2,nasusc2)
! second frequency derivative
  complex(kind=c8b), intent(out) :: chiks2(4,4,nasusc2,nasusc2)
! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
! -----------------------------------------------------------------
  complex(kind=c8b), allocatable :: gfijw(:,:), gfjiw(:,:), gfij0(:,:), gfji0(:,:), gf(:,:), dti(:,:,:,:)
  complex(kind=c8b), allocatable :: gfijdef(:,:), gfjidef(:,:), gfijef(:,:), gfjief(:,:)
  complex(kind=c8b), allocatable :: gfijdeb(:,:), gfjideb(:,:), gfijeb(:,:), gfjieb(:,:)
  complex(kind=c8b), allocatable :: work1(:,:), work2(:,:)
  complex(kind=c8b), allocatable :: suscylm0(:,:,:,:,:,:), suscylm1(:,:,:,:,:,:), suscylm2(:,:,:,:,:,:), jijsusc(:,:,:,:)
  complex(kind=c8b), allocatable :: suscy00(:,:,:,:,:,:), dosylm(:,:,:), suscden(:,:,:), rhoden(:), rhodenylm(:,:,:)
  real(kind=r8b),    allocatable :: rotmat(:,:,:)
  complex(kind=c8b) :: norm, de1, de2, e, norm2, trace
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4
  integer(kind=i4b) :: ib, jb, is, js
  integer(kind=i4b) :: i2(2), i3(3), i4(4), ie, je, ia, ja, ia2, ja2, jlm, ilm, j, i, klm, ip, ie0
  integer(kind=i4b) :: iq, jq, iq0, iq1, jq0, jq1
  integer(kind=i4b) :: ipiv(ngfsum), info
  real(kind=r8b)    :: gaunti, gauntj, doublegaunt, maxelem0, maxelem1, maxelem2, start, finish, gftime, basistime, re, im
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

  ram = 16*(4*ngfmax**2 + 13*nlmsb**2 + 16*(3*lmmax0**2+lmmax**2)*nasusc2**2 + (4*lmmax+6*lmmax0)*nasusc2 + ndensum)
  ram = ram/1024.d0**3
  write(*,'("RAM for static_susc_rsa: ",f8.3," GB")') ram
  allocate(gf(nlmsb,nlmsb),gfijw(nlmsb,nlmsb),gfjiw(nlmsb,nlmsb),gfij0(nlmsb,nlmsb),gfji0(nlmsb,nlmsb))
  allocate(gfijdef(nlmsb,nlmsb),gfjidef(nlmsb,nlmsb),gfijef(nlmsb,nlmsb),gfjief(nlmsb,nlmsb))
  allocate(gfijdeb(nlmsb,nlmsb),gfjideb(nlmsb,nlmsb),gfijeb(nlmsb,nlmsb),gfjieb(nlmsb,nlmsb))
  allocate(dti(nlmsb,nlmsb,7,nasusc2),work1(nlmsb,nlmsb),work2(nlmsb,nlmsb),rotmat(3,3,nasusc2))
  gftime = 0.d0; basistime = 0.d0
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
! Bxc and Pauli matrices in GF basis, rotation matrices to local frame
  dti = czero
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
!   dti = Ri Bi Ri \vec\sigma
    call build_dti(1,ia,.false.,dti(:,:,:,ia2),jiout(ia2))
    write(*,'(" ia=",i4,"  Jiiout=",es16.8)') ia, jiout(ia2)
!   rotations to local frame
    if (lrotsusc) then
      call rotvec(uz,magdir(:,ia),rotmat(:,:,ia2))
    else
      call rotvec(uz,uz,rotmat(:,:,ia2))
    end if
!    do i=1,3
!      write(*,'("rotmat, ia2=",i4,3es16.8)') ia2, rotmat(i,:,ia2)
!    end do
  end do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  chiks0 = czero; chiks1 = czero; chiks2 = czero
  do ja2=1,nasusc2      ! atom j
    ja = iasusc2(ja2)
    do ia2=1,nasusc2    ! atom i
      ia = iasusc2(ia2)
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
          call projected_gf(ie0+je,ja,ia,gfji0,onsite,struct)
!          if (lrot) call local_frame(ja,ia,magdir,gfji0,gf)
!         Gji(E - i0) = Gij(E + i0)^\dagger
          call projected_gf(ie0+je,ia,ja,gfjiw,onsite,struct)
!          if (lrot) call local_frame(ia,ja,magdir,gfjiw,gf)
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
            call projected_gf(ie0+ie,ia,ja,gfijw,onsite,struct)
!            if (lrot) call local_frame(ia,ja,magdir,gfijw,gf)
!           Gij(E - i0) = Gji(E + i0)^\dagger
            call projected_gf(ie0+ie,ja,ia,gfij0,onsite,struct)
!            if (lrot) call local_frame(ja,ia,magdir,gfij0,gf)
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
!           ------------------------------------------------------------
            do j=1,4
              do i=1,4
!               \sigma_i G_ij(Ei+i0)  -->  work1
                call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ia),cone,dti(:,:,i+3,ia2),nlmsb,gijw,nlmsb,czero,work1,nlmsb)
!               (\sigma_i * Gij(Ei+i0)) * \sigma_j --> work2
                call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,work1,nlmsb,dti(:,:,j+3,ja2),nlmsb,czero,work2,nlmsb)
!               (\sigma_i * Gij(Ei+i0) * \sigma_j) * Gji(Ej+i0) --> work1
                call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ja),cone,work2,nlmsb,gji0,nlmsb,czero,work1,nlmsb)
                trace = czero
                do ilms=1,nlmsba(ia)
                  trace = trace + work1(ilms,ilms)
                end do
!               Static susceptibility: - 1/i2pi \int dz Gij(z) Gji(z)
                kssusc0(i,j,ia2,ja2) = kssusc0(i,j,ia2,ja2) - w00(ie,je)*trace/i2pi
!               First frequency derivative: - 1/i2pi \int dz d/dz Gij(z) Gji(z)
                kssusc1(i,j,ia2,ja2) = kssusc1(i,j,ia2,ja2) - w10(ie,je)*trace/i2pi
!               Second frequency derivative: -1/i2pi \int dz^* d/dz^* Gji(z^*)^\dagger d/dz^* Gij(z^*)^\dagger
                kssusc2(i,j,ia2,ja2) = kssusc2(i,j,ia2,ja2) + w11(ie,je)*trace/i2pi
!               \sigma_i G_ij(Ei-i0)  -->  work1
                call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ia),cone,dti(:,:,i+3,ia2),nlmsb,gij0,nlmsb,czero,work1,nlmsb)
!               (\sigma_i * Gij(Ei-i0)) * \sigma_j --> work2
                call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,work1,nlmsb,dti(:,:,j+3,ja2),nlmsb,czero,work2,nlmsb)
!               (\sigma_i * Gij(Ei-i0) * \sigma_j) * Gji(Ej-i0) --> work1
                call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ja),cone,work2,nlmsb,gjiw,nlmsb,czero,work1,nlmsb)
                trace = czero
                do ilms=1,nlmsba(ia)
                  trace = trace + work1(ilms,ilms)
                end do
!               Static susceptibility: +1/i2pi \int dz^* Gji(z^*)^\dagger Gij(z^*)^\dagger
                kssusc0(i,j,ia2,ja2) = ksusc0(i,j,ia2,ja2) + conjg(w00(ie,je))*trace/i2pi
!               First frequency derivative: -1/i2pi \int dz^* Gji(z^*)^\dagger d/dz^* Gij(z^*)^\dagger
                kssusc1(i,j,ia2,ja2) = ksusc1(i,j,ia2,ja2) - conjg(w01(ie,je))*trace/i2pi
!               Second frequency derivative: -1/i2pi \int dz^* d/dz^* Gji(z^*)^\dagger d/dz^* Gij(z^*)^\dagger
                kssusc2(i,j,ia2,ja2) = ksusc2(i,j,ia2,ja2) - conjg(w11(ie,je))*trace/i2pi
              end do
            end do
!           ------------------------------------------------------------
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
      do j=1,4
        do i=1,4
!         \sigma_i G_ij(EF+i0)  -->  work1
          call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ia),cone,dti(:,:,i+3,ia2),nlmsb,gij0,nlmsb,czero,work1,nlmsb)
!         (\sigma_i * Gij(Ei-i0)) * \sigma_j --> work2
          call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,work1,nlmsb,dti(:,:,j+3,ja2),nlmsb,czero,work2,nlmsb)
!         (\sigma_i * Gij(Ei-i0) * \sigma_j) * Gji(Ej-i0) --> work1
          call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ja),cone,work2,nlmsb,gjiw,nlmsb,czero,work1,nlmsb)
          trace = czero
          do ilms=1,nlmsba(ia)
            trace = trace + work1(ilms,ilms)
          end do
!         Nothing for static susc
!         First derivative: 1/i2pi Gij(EF+i0) Gij(EF+i0)^\dagger - 1/i2pi Gij(EB+i0) Gij(EB+i0)^\dagger
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
! ----------------------------------------------------------------------
! Final filtering (global)
  maxelem0 = maxval(abs(suscylm0))
  where (abs(suscylm0) < susctol*maxelem0) suscylm0 = 0.d0
  maxelem1 = maxval(abs(suscylm1))
  where (abs(suscylm1) < susctol*maxelem1) suscylm1 = 0.d0
  maxelem2 = maxval(abs(suscylm2))
  where (abs(suscylm2) < susctol*maxelem2) suscylm2 = 0.d0
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
  write(*,'(/," Static  KS susc chiGF time=",f10.3," s")') gftime
  write(*,'(" Static  KS susc basis time=",f10.3," s")') basistime
! All done!
  end subroutine static_susc_rsa
