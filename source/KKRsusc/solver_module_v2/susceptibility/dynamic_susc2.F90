  subroutine dynamic_susc2(onsite,struct)
! Dynamic KS susceptibility
! Positive and negative frequencies
! Using density basis
! No fit of GF
  use global

  implicit none

  logical,           intent(in)  :: onsite, struct
! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
  complex(kind=c8b), parameter :: e0up = (0.481995d0,-0.05d0), e0dn = (0.711995,-0.007d0)
! -----------------------------------------------------------------
  complex(kind=c8b) :: suscylm(4,4,lmmax0,lmmax0,nasusc2,nasusc2)
  complex(kind=c8b) :: suscy00(4,4,lmmax,lmmax,nasusc2,nasusc2)
  complex(kind=c8b), allocatable :: gfijw(:,:), gfjiw(:,:)
  complex(kind=c8b), allocatable :: gfij0(:,:), gfji0(:,:)
  complex(kind=c8b), allocatable :: suscblock(:,:), suscwork(:,:), gf(:,:)
  complex(kind=c8b), allocatable :: zw(:), chi0escf(:,:), chi0eb(:,:), chi0ef(:,:)
  complex(kind=c8b) :: norm, de, e, norm2
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4
  integer(kind=i4b) :: ib, jb, is, js
  integer(kind=i4b) :: i2(2), i3(3), i4(4), ie, ia, ja, ia2, ja2, jlm, ilm, j, i
  integer(kind=i4b) :: iw, jw, iew, iew0, jew, jew0
  integer(kind=i4b) :: iq, jq, iq0, iq1, jq0, jq1
  integer(kind=i4b) :: ipiv(ndensum), info
  real(kind=r8b)    :: gaunti, gauntj, doublegaunt, maxelem, start, finish, gftime, basistime
  real(kind=r8b)    :: dr(nrmax), uparam
  complex(kind=c8b) :: work(nrmax), block(4,4)
  complex(kind=c8b) :: gfij0up, gfij0dn, gfji0up, gfji0dn, gfijwup, gfijwdn, gfjiwup, gfjiwdn, desum
  complex(kind=c8b) :: chi0_dnup, chi0escf_dnup, chi0eb_dnup, chi0ef_dnup
  complex(kind=c8b) :: chi0_updn, chi0escf_updn, chi0eb_updn, chi0ef_updn
  integer(kind=i4b) :: ne, nr, nescf2, nw, nweb, nwef
  real(kind=r8b)    :: wmax
  complex(kind=c8b), external :: radint
!  real(kind=r8b)    :: r(nalmsb), c(nalmsb), rwork(2*nalmsb), rcond, ferr(nalmsb), berr(nalmsb)
!  complex(kind=c8b) :: work(2*nalmsb), temp(nalmsb,nalmsb)
!  complex(kind=c8b) :: x(nalmsb,nalmsb)
!  character*1       :: equed

! First get info about the frequency mesh
  open(file='wmesh.dat',unit=iofile,status='old')
  read(iofile,*) ! comment line
  read(iofile,*) wmax, nescf2, nw, nweb, nwef, uparam
  uparam = 1.d0/uparam
  write(*,*) "uparam=", uparam
  close(iofile)
! Equidistant frequencies
  allocate(zw(0:nw))
  do iw=0,nw
    zw(iw) = 0.5d0*wmax*(1.d0 + (2.d0*iw-nw)/nw)
  end do
! Check options and allocate memory
  write(*,'(" Dynamic KS max RAM=",f10.3," GB")') 16.d0*(3*ndensum**2 + 2*ngfmax**2)/1024.d0**3
  allocate(suscblock(ngfmax,ngfmax),suscwork(ngfmax,ngfmax),gf(nlmsb,nlmsb))
  allocate(gfij0(nlmsb,nlmsb),gfji0(nlmsb,nlmsb),gfijw(nlmsb,nlmsb),gfjiw(nlmsb,nlmsb))
!  if (lanalytic) allocate(chi0escf(ndensum,ndensum))
!  if (nonanalytic .and. nweb > 0) allocate(chi0eb(ndensum,ndensum))
!  if (nonanalytic .and. nwef > 0) allocate(chi0ef(ndensum,ndensum))
!  gftime = 0.d0; basistime = 0.d0
  open(file='dynsusc.dat',unit=iofile,status='replace')
! ********************
! Negative frequencies
  do iw=nw,1,1
! ********************
    write(*,'(i4,2es16.8)') iw, -zw(iw)
!    write(*,'("test element=",2i4)') lmsb2i(1,5,2,1),lmsb2i(1,5,1,1)
    suscylm = czero; suscy00 = czero; maxelem = 0.d0
    chi0_dnup = czero; chi0escf_dnup = czero; chi0eb_dnup = czero; chi0ef_dnup = czero
    chi0_updn = czero; chi0escf_updn = czero; chi0eb_updn = czero; chi0ef_updn = czero
    kssusc0 = czero
!    if (lanalytic) chi0escf = czero
!    if (nonanalytic .and. nweb > 0) chi0eb = czero
!    if (nonanalytic .and. nwef > 0) chi0ef = czero
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    do ja2=1,nasusc2      ! atom j
      ja = iasusc2(ja2)
      jq0 = sum(nalmsbgf(1:ja-1))
      jq1 = sum(nalmsbgf(1:ja))
      do ia2=1,nasusc2    ! atom i
        ia = iasusc2(ia2)
        iq0 = sum(nalmsbgf(1:ia-1))
        iq1 = sum(nalmsbgf(1:ia))
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       --------------------------------------------------------------
        if (lanalytic) then
!       --------------------------------------------------------------
          call cpu_time(start)
          suscblock = czero; desum = czero
!         Integration on usual box contour
          iew0 = nescf + 2*(iw-1)*(nescf+nweb+nwef) 
          do ie=1,nescf     ! energy
!            write(*,*) "ie", ie
            iew = iew0 + ie
            jew = iew0 + nescf + ie
            de = desusc(ie)
            desum = desum + de
!            write(*,'(2i8,8es16.8)') iew, jew, esusc(iew), esusc(jew), esusc(iew) - esusc(jew), de
!           get the projected GFs
!           E - w + i0
            call projected_gf(iew,ia,ja,gfijw,onsite,struct)
            if (lrot) call local_frame(ia,ja,magdir,gfijw,gf)
!           E + i0
            call projected_gf(ie,ja,ia,gfji0,onsite,struct)
            if (lrot) call local_frame(ia,ja,magdir,gfij0,gf)
!           E - i0
            gfij0 = conjg(transpose(gfji0))
!           E + w - i0
            call projected_gf(jew,ia,ja,gfjiw,onsite,struct)
            if (lrot) call local_frame(ja,ia,magdir,gfjiw,gf)
            gfjiw = conjg(transpose(gfjiw))
            do j=1,jq1-jq0
              jq = jq0 + j
              i3 = i2almsbgf(:,jq)
              q3 = i3(1); q4 = i3(2)!; ja = i3(3)
              do i=1,iq1-iq0
                iq = iq0 + i
                i3 = i2almsbgf(:,iq)
                q1 = i3(1); q2 = i3(2)!; ia = i3(3)
                norm = (conjg(de)*gfij0(q2,q3)*gfjiw(q4,q1) - de*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
                suscblock(i,j) = suscblock(i,j) + norm
              end do
            end do
            gfijwup = 1.d0/(esusc(iew) - e0up); gfijwdn = 1.d0/(esusc(iew) - e0dn)
            gfjiwup = 1.d0/(esusc(jew) - e0up); gfjiwdn = 1.d0/(esusc(jew) - e0dn)
            gfij0up = 1.d0/(esusc(ie) - e0up);  gfij0dn = 1.d0/(esusc(ie) - e0dn)
            gfji0up = 1.d0/(esusc(ie) - e0up);  gfji0dn = 1.d0/(esusc(ie) - e0dn)
            gfij0up = conjg(gfij0up); gfij0dn = conjg(gfij0dn)
            gfjiwup = conjg(gfjiwup); gfjiwdn = conjg(gfjiwdn)
!            gfij0up = gfij0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfij0dn = gfij0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!            gfji0up = gfji0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfji0dn = gfji0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!            gfijwup = gfijw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfijwdn = gfijw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!            gfjiwup = gfjiw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfjiwdn = gfjiw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
            chi0escf_updn = chi0escf_updn + (conjg(de)*gfij0up*gfjiwdn - de*gfijwup*gfji0dn)/i2pi
            chi0escf_dnup = chi0escf_dnup + (conjg(de)*gfij0dn*gfjiwup - de*gfijwdn*gfji0up)/i2pi
          end do            ! energy
!          write(*,'("Analytic desum=",2es16.8)') desum
          call cpu_time(finish)
          gftime = gftime + finish - start
          call cpu_time(start)
!          call add_susc_block(ia2,ja2,suscblock,suscwork,chi0escf)
          call add_susc_block(ia2,ja2,suscblock,suscwork,kssusc0)
          call cpu_time(finish)
          basistime = basistime + finish - start
        end if
!       ----------------------------------------------------------------
        if (lnonanalytic .and. nweb > 0) then
!       ----------------------------------------------------------------
          call cpu_time(start)
!         Integration on real axis around bottom of box contour
          suscblock = czero; desum = czero
!         ++++++++++++++++
!         loop over panels
          do jw=1,iw
!         ++++++++++++++++
!           panel for Gij(E)
            iew0 = nescf + 2*(jw-1)*(nescf+nweb+nwef) + 2*nescf + nweb
!           panel for Gji(E-w)
            jew0 = nescf + 2*(iw-jw)*(nescf+nweb+nwef) + 2*nescf + nweb
            do ie=0,nweb  ! energy
!              write(*,*) "ie", ie
!             find energies in storage
              iew = iew0 + ie
              jew = jew0 - ie
              de = desusc(iew)
!             handle ebottom (not repeated in energy mesh)
              if (ie == 0) then
                if (jw == 1) then
                  iew = 1
                  de = desusc(iew0+nweb)
                else
                  iew = iew0 - 2*(nescf+nwef) - nweb
                  de = desusc(iew)
                end if
              end if
              if (ie == nweb) then
                jew = jew0 - 2*(nescf+nweb+nwef) 
                if (jw == iw) jew = 1
              end if
              desum = desum + de
!              write(*,'(4i8,8es16.8)') jew0, iew0, jew, iew, esusc(jew), esusc(iew), esusc(jew) - esusc(iew), de
!             get the projected GFs
!             E - w + i0
              call projected_gf(jew,ia,ja,gfijw,onsite,struct)
              if (lrot) call local_frame(ia,ja,magdir,gfijw,gf)
!             E - i0
              call projected_gf(iew,ia,ja,gfji0,onsite,struct)
              if (lrot) call local_frame(ia,ja,magdir,gfji0,gf)
              gfji0 = conjg(transpose(gfji0))
              do j=1,jq1-jq0
                jq = jq0 + j
                i3 = i2almsbgf(:,jq)
                q3 = i3(1); q4 = i3(2)!; ja = i3(3)
                do i=1,iq1-iq0
                  iq = iq0 + i
                  i3 = i2almsbgf(:,iq)
                  q1 = i3(1); q2 = i3(2)!; ia = i3(3)
                  norm = de*gfijw(q2,q3)*gfji0(q4,q1)/i2pi
                  suscblock(i,j) = suscblock(i,j) + norm
                end do
              end do
              gfijwup = 1.d0/(esusc(jew) - e0up); gfijwdn = 1.d0/(esusc(jew) - e0dn)
              gfji0up = 1.d0/(esusc(iew) - e0up); gfji0dn = 1.d0/(esusc(iew) - e0dn)
              gfji0up = conjg(gfji0up); gfji0dn = conjg(gfji0dn)
!              gfijwup = gfijw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfijwdn = gfijw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!              gfji0up = gfji0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfji0dn = gfji0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
              chi0eb_updn = chi0eb_updn + de*gfijwup*gfji0dn/i2pi
              chi0eb_dnup = chi0eb_dnup + de*gfijwdn*gfji0up/i2pi
            end do  ! energy
!         ++++++
          end do
!         ++++++
          write(*,'("Ebottom desum=",2es16.8)') desum
          call cpu_time(finish)
          gftime = gftime + finish - start
          call cpu_time(start)
!          call add_susc_block(ia2,ja2,suscblock,suscwork,chi0eb)
          call add_susc_block(ia2,ja2,suscblock,suscwork,kssusc0)
          call cpu_time(finish)
          basistime = basistime + finish - start
!       --------------------------------------------------------------
        end if
!       ----------------------------------------------------------------
        if (lnonanalytic .and. nwef > 0) then
!       ----------------------------------------------------------------
          call cpu_time(start)
!         Integration on real axis around Efermi
          suscblock = czero; desum = czero
!         ++++++++++++++++
!         loop over panels
          do jw=1,iw
!         ++++++++++++++++
!           panel for Gij(E)
            iew0 = nescf + 2*(jw-1)*(nescf+nweb+nwef) + 2*(nescf+nweb) + nwef
!           panel for Gji(E-w)
            jew0 = nescf + 2*(iw-jw)*(nescf+nweb+nwef) + 2*(nescf+nweb) + nwef
            do ie=0,nwef  ! energy
!              write(*,*) "ie", ie
!             find energies in storage
              iew = iew0 + ie
              jew = jew0 - ie
              de = desusc(iew)
!             handle efermi (not repeated in energy mesh)
              if (ie == 0) then
                if (jw == 1) then
                  iew = nescf
                  de = desusc(iew0+nwef)
                else
                  iew = iew0 - 2*(nescf+nweb) - nwef
                  de = desusc(iew)
                end if
              end if
              if (ie == nwef) then
                jew = jew0 - 2*(nescf+nweb+nwef)
                if (jw == iw) jew = nescf
              end if
              desum = desum + de
!              write(*,'(4i8,8es16.8)') jew0, iew0, jew, iew, esusc(jew), esusc(iew), esusc(jew) - esusc(iew), de
!             get the projected GFs
!             E - w + i0
              call projected_gf(jew,ia,ja,gfijw,onsite,struct)
              if (lrot) call local_frame(ia,ja,magdir,gfijw,gf)
!             E - i0
              call projected_gf(iew,ia,ja,gfji0,onsite,struct)
              if (lrot) call local_frame(ia,ja,magdir,gfji0,gf)
              gfji0 = conjg(transpose(gfji0))
              do j=1,jq1-jq0
                jq = jq0 + j
                i3 = i2almsbgf(:,jq)
                q3 = i3(1); q4 = i3(2)!; ja = i3(3)
                do i=1,iq1-iq0
                  iq = iq0 + i
                  i3 = i2almsbgf(:,iq)
                  q1 = i3(1); q2 = i3(2)!; ia = i3(3)
                  norm = -de*gfijw(q2,q3)*gfji0(q4,q1)/i2pi
                  suscblock(i,j) = suscblock(i,j) + norm
                end do
              end do
              gfijwup = 1.d0/(esusc(jew) - e0up); gfijwdn = 1.d0/(esusc(jew) - e0dn)
              gfji0up = 1.d0/(esusc(iew) - e0up); gfji0dn = 1.d0/(esusc(iew) - e0dn)
              gfji0up = conjg(gfji0up); gfji0dn = conjg(gfji0dn)
!              gfijwup = gfijw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfijwdn = gfijw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!              gfji0up = gfji0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfji0dn = gfji0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
              chi0ef_updn = chi0ef_updn - de*gfijwup*gfji0dn/i2pi
              chi0ef_dnup = chi0ef_dnup - de*gfijwdn*gfji0up/i2pi
            end do  ! energy
!         ++++++
          end do
!         ++++++
!          write(*,'("EFermi  desum=",2es16.8)') desum
          call cpu_time(finish)
          gftime = gftime + finish - start
          call cpu_time(start)
!          call add_susc_block(ia2,ja2,suscblock,suscwork,chi0ef)
          call add_susc_block(ia2,ja2,suscblock,suscwork,kssusc0)
          call cpu_time(finish)
          basistime = basistime + finish - start
!       ----------------------------------------------------------------
        end if
!       ----------------------------------------------------------------
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      end do            ! atom i
    end do              ! atom j
!   Add contributions to KS susc
!    if (lanalytic) kssusc0 = kssusc0 + chi0escf
!    if (nonanalytic .and. nweb > 0) kssusc0 = kssusc0 + chi0eb
!    if (nonanalytic .and. nwef > 0) kssusc0 = kssusc0 + chi0ef
    chi0_updn = chi0escf_updn + chi0eb_updn + chi0ef_updn
    chi0_dnup = chi0escf_dnup + chi0eb_dnup + chi0ef_dnup
    write(555,'(6es16.8)') -zw(iw), chi0_updn, chi0_dnup
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    where (abs(kssusc0) < susctol) kssusc0 = 0.d0
!   ------------------------------------------------------------------
    if (lenhanced) then
      call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc0,ndensum,kernel,ndensum,czero,denominator,ndensum)
      do iq=1,ndensum
        denominator(iq,iq) = denominator(iq,iq) + 1.d0
      end do
      call zgesv(ndensum,ndensum,denominator,ndensum,ipiv,kssusc0,ndensum,info)
      if (info /= 0) stop 'dyn_susc_real2: failure in zgesv'
!      call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc0,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
!      write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
!      write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
!      if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
!      kssusc0 = x
    end if
!   ------------------------------------------------------------------
!   Symmetrize
!    call symmetrize(nalmsb,kssusc0,susctol)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   multipoles
    suscylm = czero
    do jq=1,ndensum
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
      do iq=1,ndensum
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
        suscylm(is,js,ilm,jlm,ia2,ja2) = suscylm(is,js,ilm,jlm,ia2,ja2) + suscnorm(iq)*kssusc0(iq,jq)*suscnorm(jq)
      end do
    end do
!   Temporary xc kernel
!    suscylm(1,2,1,1,1,1) = suscylm(1,2,1,1,1,1)/(1.d0 - uparam*suscylm(1,2,1,1,1,1))
!    suscylm(2,1,1,1,1,1) = suscylm(2,1,1,1,1,1)/(1.d0 - uparam*suscylm(2,1,1,1,1,1))
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
              block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*suscylm(is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
            end do
            end do
          end do
          end do
          suscylm(:,:,ilm,jlm,ia2,ja2) = block
        end do
        end do
      end do
      end do
    end if
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Final filtering (global)
!    maxelem = maxval(abs(suscylm))
!    where (abs(suscylm) < susctol*maxelem) suscylm = 0.d0
!    maxelem = maxval(abs(suscy00))
!    where (abs(suscy00) < susctol*maxelem) suscy00 = 0.d0
!    write(*,'(/,"KS susc ia, ja, im, il, jm, jl, susc ylm")')
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do ia2=1,nasusc2
        ia = iasusc2(ia2)
        do jlm=1,lmmax0
          do ilm=1,lmmax0
            norm = sum(abs(suscylm(:,:,ilm,jlm,ia2,ja2)))
            if (abs(norm) > susctol) then
              do j=1,4
!                write(*,'(6i4,8f16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscylm(i,j,ilm,jlm,ia2,ja2),i=1,4)
              end do
            end if
          end do
        end do
      end do
    end do
    write(iofile,'(1000es16.8)') -real(zw(iw)), ((((suscylm(i,j,1,1,ia2,ja2),i=1,4),j=1,4),ia2=1,nasusc2),ja2=1,nasusc2)
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ******
  end do
! ******
! ********************
! Zero frequency
! ********************
  write(*,'(i4,2es16.8)') 0, czero
!  write(*,'("test element=",2i4)') lmsb2i(1,5,2,1),lmsb2i(1,5,1,1)
  suscylm = czero; suscy00 = czero; maxelem = 0.d0
  chi0_dnup = czero; chi0escf_dnup = czero; chi0eb_dnup = czero; chi0ef_dnup = czero
  chi0_updn = czero; chi0escf_updn = czero; chi0eb_updn = czero; chi0ef_updn = czero
  kssusc0 = czero!; chi0escf = czero
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ja2=1,nasusc2      ! atom j
    ja = iasusc2(ja2)
    jq0 = sum(nalmsbgf(1:ja-1))
    jq1 = sum(nalmsbgf(1:ja))
    do ia2=1,nasusc2    ! atom i
      ia = iasusc2(ia2)
      iq0 = sum(nalmsbgf(1:ia-1))
      iq1 = sum(nalmsbgf(1:ia))
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ----------------------------------------------------------------
      if (lanalytic) then
!     ----------------------------------------------------------------
        call cpu_time(start)
        suscblock = czero; desum = czero
!       Integration on usual box contour
        do ie=1,nescf     ! energy
!          write(*,*) "ie", ie
          de = desusc(ie)
          desum = desum + de
!          write(*,'(2i8,8es16.8)') iew, jew, esusc(iew), esusc(jew), esusc(iew) - esusc(jew), de
!         get the projected GFs
!         E + w + i0
          call projected_gf(ie,ia,ja,gfijw,onsite,struct)
          if (lrot) call local_frame(ia,ja,magdir,gfijw,gf)
!         E + i0
          call projected_gf(ie,ja,ia,gfji0,onsite,struct)
          if (lrot) call local_frame(ja,ia,magdir,gfji0,gf)
!         E - i0
          gfij0 = conjg(transpose(gfji0))
!         E - w - i0
          gfjiw = conjg(transpose(gfijw))
          do j=1,jq1-jq0
            jq = jq0 + j
            i3 = i2almsbgf(:,jq)
            q3 = i3(1); q4 = i3(2)!; ja = i3(3)
            do i=1,iq1-iq0
              iq = iq0 + i
              i3 = i2almsbgf(:,iq)
              q1 = i3(1); q2 = i3(2)!; ia = i3(3)
              norm = (conjg(de)*gfij0(q2,q3)*gfjiw(q4,q1) - de*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
              suscblock(i,j) = suscblock(i,j) + norm
            end do
          end do
!         E + w + i0
          gfijwup = 1.d0/(esusc(ie) - e0up); gfijwdn = 1.d0/(esusc(ie) - e0dn)
!         E + i0
          gfji0up = 1.d0/(esusc(ie) - e0up); gfji0dn = 1.d0/(esusc(ie) - e0dn)
!         E - i0
          gfij0up = 1.d0/(esusc(ie) - e0up); gfij0dn = 1.d0/(esusc(ie) - e0dn)
          gfij0up = conjg(gfij0up); gfij0dn = conjg(gfij0dn)
!         E - w - i0
          gfjiwup = 1.d0/(esusc(jew) - e0up); gfjiwdn = 1.d0/(esusc(jew) - e0dn)
          gfjiwup = conjg(gfjiwup); gfjiwdn = conjg(gfjiwdn)
!          gfij0up = gfij0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfij0dn = gfij0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!          gfji0up = gfji0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfji0dn = gfji0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!          gfijwup = gfijw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfijwdn = gfijw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!          gfjiwup = gfjiw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfjiwdn = gfjiw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
          chi0escf_updn = chi0escf_updn + (conjg(de)*gfij0up*gfjiwdn - de*gfijwup*gfji0dn)/i2pi
          chi0escf_dnup = chi0escf_dnup + (conjg(de)*gfij0dn*gfjiwup - de*gfijwdn*gfji0up)/i2pi
        end do            ! energy
!        write(*,'("Analytic desum=",2es16.8)') desum
        call cpu_time(finish)
        gftime = gftime + finish - start
        call cpu_time(start)
!        call add_susc_block(ia2,ja2,suscblock,suscwork,chi0escf)
        call add_susc_block(ia2,ja2,suscblock,suscwork,kssusc0)
        call cpu_time(finish)
        basistime = basistime + finish - start
!     ----------------------------------------------------------------
      end if
!     ----------------------------------------------------------------
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end do            ! atom i
  end do              ! atom j
! Add contributions to KS susc
!  kssusc0 = kssusc0 + chi0escf
  chi0_updn = chi0escf_updn + chi0eb_updn + chi0ef_updn
  chi0_dnup = chi0escf_dnup + chi0eb_dnup + chi0ef_dnup
  write(555,'(6es16.8)') czero, chi0_updn, chi0_dnup
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  where (abs(kssusc0) < susctol) kssusc0 = 0.d0
! ------------------------------------------------------------------
  if (lenhanced) then
    call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc0,ndensum,kernel,ndensum,czero,denominator,ndensum)
    do iq=1,ndensum
      denominator(iq,iq) = denominator(iq,iq) + 1.d0
    end do
    call zgesv(ndensum,ndensum,denominator,ndensum,ipiv,kssusc0,ndensum,info)
    if (info /= 0) stop 'dyn_susc_real2: failure in zgesv'
!    call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc0,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
!    write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
!    write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
!    if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
!    kssusc0 = x
  end if
! ------------------------------------------------------------------
! Symmetrize
!  call symmetrize(nalmsb,kssusc0,susctol)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! multipoles
  suscylm = czero
  do jq=1,ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
    do iq=1,ndensum
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
      suscylm(is,js,ilm,jlm,ia2,ja2) = suscylm(is,js,ilm,jlm,ia2,ja2) + suscnorm(iq)*kssusc0(iq,jq)*suscnorm(jq)
    end do
  end do
! Temporary xc kernel
!  suscylm(1,2,1,1,1,1) = suscylm(1,2,1,1,1,1)/(1.d0 - uparam*suscylm(1,2,1,1,1,1))
!  suscylm(2,1,1,1,1,1) = suscylm(2,1,1,1,1,1)/(1.d0 - uparam*suscylm(2,1,1,1,1,1))
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
            block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*suscylm(is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
          end do
          end do
        end do
        end do
        suscylm(:,:,ilm,jlm,ia2,ja2) = block
      end do
      end do
    end do
    end do
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Final filtering (global)
!  maxelem = maxval(abs(suscylm))
!  where (abs(suscylm) < susctol*maxelem) suscylm = 0.d0
!  maxelem = maxval(abs(suscy00))
!  where (abs(suscy00) < susctol*maxelem) suscy00 = 0.d0
!  write(*,'(/,"KS susc ia, ja, im, il, jm, jl, susc ylm")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm=1,lmmax0
        do ilm=1,lmmax0
          norm = sum(abs(suscylm(:,:,ilm,jlm,ia2,ja2)))
          if (abs(norm) > susctol) then
            do j=1,4
!              write(*,'(6i4,8f16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscylm(i,j,ilm,jlm,ia2,ja2),i=1,4)
            end do
          end if
        end do
      end do
    end do
  end do
  write(iofile,'(1000es16.8)') real(czero), ((((suscylm(i,j,1,1,ia2,ja2),i=1,4),j=1,4),ia2=1,nasusc2),ja2=1,nasusc2)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ********************
! Positive frequencies
  do iw=1,nw
! ********************
    write(*,'(i4,2es16.8)') iw, zw(iw)
!    write(*,'("test element=",2i4)') lmsb2i(1,5,2,1),lmsb2i(1,5,1,1)
    suscylm = czero; suscy00 = czero; maxelem = 0.d0
    chi0_dnup = czero; chi0escf_dnup = czero; chi0eb_dnup = czero; chi0ef_dnup = czero
    chi0_updn = czero; chi0escf_updn = czero; chi0eb_updn = czero; chi0ef_updn = czero
    kssusc0 = czero
!    if (lanalytic) chi0escf = czero
!    if (nonanalytic .and. nweb > 0) chi0eb = czero
!    if (nonanalytic .and. nwef > 0) chi0ef = czero
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    do ja2=1,nasusc2      ! atom j
      ja = iasusc2(ja2)
      jq0 = sum(nalmsbgf(1:ja-1))
      jq1 = sum(nalmsbgf(1:ja))
      do ia2=1,nasusc2    ! atom i
        ia = iasusc2(ia2)
        iq0 = sum(nalmsbgf(1:ia-1))
        iq1 = sum(nalmsbgf(1:ia))
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!       ----------------------------------------------------------------
        if (lanalytic) then
!       ----------------------------------------------------------------
          call cpu_time(start)
          suscblock = czero; desum = czero
!         Integration on usual box contour
          iew0 = nescf + 2*(iw-1)*(nescf+nweb+nwef) 
          do ie=1,nescf     ! energy
!            write(*,*) "ie", ie
            iew = iew0 + nescf + ie
            jew = iew0 + ie
            de = desusc(ie)
            desum = desum + de
!            write(*,'(2i8,8es16.8)') iew, jew, esusc(iew), esusc(jew), esusc(iew) - esusc(jew), de
!           get the projected GFs
!           E + w + i0
            call projected_gf(iew,ia,ja,gfijw,onsite,struct)
            if (lrot) call local_frame(ia,ja,magdir,gfijw,gf)
!           E + i0
            call projected_gf(ie,ja,ia,gfji0,onsite,struct)
            if (lrot) call local_frame(ja,ia,magdir,gfji0,gf)
!           E - i0
            gfij0 = conjg(transpose(gfji0))
!           E - w - i0
            call projected_gf(jew,ia,ja,gfjiw,onsite,struct)
            if (lrot) call local_frame(ia,ja,magdir,gfjiw,gf)
            gfjiw = conjg(transpose(gfjiw))
            do j=1,jq1-jq0
              jq = jq0 + j
              i3 = i2almsbgf(:,jq)
              q3 = i3(1); q4 = i3(2)!; ja = i3(3)
              do i=1,iq1-iq0
                iq = iq0 + i
                i3 = i2almsbgf(:,iq)
                q1 = i3(1); q2 = i3(2)!; ia = i3(3)
                norm = (conjg(de)*gfij0(q2,q3)*gfjiw(q4,q1) - de*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
                suscblock(i,j) = suscblock(i,j) + norm
              end do
            end do
!           E + w + i0
            gfijwup = 1.d0/(esusc(iew) - e0up); gfijwdn = 1.d0/(esusc(iew) - e0dn)
!           E + i0
            gfji0up = 1.d0/(esusc(ie) - e0up);  gfji0dn = 1.d0/(esusc(ie) - e0dn)
!           E - i0
            gfij0up = 1.d0/(esusc(ie) - e0up);  gfij0dn = 1.d0/(esusc(ie) - e0dn)
            gfij0up = conjg(gfij0up); gfij0dn = conjg(gfij0dn)
!           E - w - i0
            gfjiwup = 1.d0/(esusc(jew) - e0up); gfjiwdn = 1.d0/(esusc(jew) - e0dn)
            gfjiwup = conjg(gfjiwup); gfjiwdn = conjg(gfjiwdn)
!            gfij0up = gfij0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfij0dn = gfij0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!            gfji0up = gfji0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfji0dn = gfji0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!            gfijwup = gfijw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfijwdn = gfijw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!            gfjiwup = gfjiw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfjiwdn = gfjiw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
            chi0escf_updn = chi0escf_updn + (conjg(de)*gfij0up*gfjiwdn - de*gfijwup*gfji0dn)/i2pi
            chi0escf_dnup = chi0escf_dnup + (conjg(de)*gfij0dn*gfjiwup - de*gfijwdn*gfji0up)/i2pi
          end do            ! energy
!          write(*,'("Analytic desum=",2es16.8)') desum
          call cpu_time(finish)
          gftime = gftime + finish - start
          call cpu_time(start)
!          call add_susc_block(ia2,ja2,suscblock,suscwork,chi0escf)
          call add_susc_block(ia2,ja2,suscblock,suscwork,kssusc0)
          call cpu_time(finish)
          basistime = basistime + finish - start
!       ----------------------------------------------------------------
        end if
!       ----------------------------------------------------------------
        if (lnonanalytic .and. nweb > 0) then
!       ----------------------------------------------------------------
          call cpu_time(start)
!         Integration on real axis around bottom of box contour
          suscblock = czero; desum = czero
!         ++++++++++++++++
!         loop over panels
          do jw=1,iw
!         ++++++++++++++++
!           panel for Gij(E)
            iew0 = nescf + 2*(jw-1)*(nescf+nweb+nwef) + 2*nescf + nweb
!           panel for Gji(E-w)
            jew0 = nescf + 2*(iw-jw)*(nescf+nweb+nwef) + 2*nescf + nweb
            do ie=0,nweb  ! energy
!              write(*,*) "ie", ie
!             find energies in storage
              iew = iew0 + ie
              jew = jew0 - ie
              de = desusc(iew)
!             handle ebottom (not repeated in energy mesh)
              if (ie == 0) then
                if (jw == 1) then
                  iew = 1
                  de = desusc(iew0+nweb)
                else
                  iew = iew0 - 2*(nescf+nwef) - nweb
                  de = desusc(iew)
                end if
              end if
              if (ie == nweb) then
                jew = jew0 - 2*(nescf+nweb+nwef) 
                if (jw == iw) jew = 1
              end if
              desum = desum + de
!              write(*,'(4i8,8es16.8)') iew0, jew0, iew, jew, esusc(iew), esusc(jew), esusc(iew) - esusc(jew), de
!             get the projected GFs
!             E + i0
              call projected_gf(iew,ia,ja,gfij0,onsite,struct)
              if (lrot) call local_frame(ia,ja,magdir,gfij0,gf)
!             E - w - i0
              call projected_gf(jew,ia,ja,gfjiw,onsite,struct)
              if (lrot) call local_frame(ia,ja,magdir,gfjiw,gf)
              gfjiw = conjg(transpose(gfjiw))
              do j=1,jq1-jq0
                jq = jq0 + j
                i3 = i2almsbgf(:,jq)
                q3 = i3(1); q4 = i3(2)!; ja = i3(3)
                do i=1,iq1-iq0
                  iq = iq0 + i
                  i3 = i2almsbgf(:,iq)
                  q1 = i3(1); q2 = i3(2)!; ia = i3(3)
                  norm = -de*gfij0(q2,q3)*gfjiw(q4,q1)/i2pi
                  suscblock(i,j) = suscblock(i,j) + norm
                end do
              end do
              gfij0up = 1.d0/(esusc(iew) - e0up); gfij0dn = 1.d0/(esusc(iew) - e0dn)
              gfjiwup = 1.d0/(esusc(jew) - e0up); gfjiwdn = 1.d0/(esusc(jew) - e0dn)
              gfjiwup = conjg(gfjiwup); gfjiwdn = conjg(gfjiwdn)
!              gfij0up = gfij0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfij0dn = gfij0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!              gfjiwup = gfjiw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfjiwdn = gfjiw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
              chi0eb_updn = chi0eb_updn - de*gfij0up*gfjiwdn/i2pi
              chi0eb_dnup = chi0eb_dnup - de*gfij0dn*gfjiwup/i2pi
            end do  ! energy
!         ++++++
          end do
!         ++++++
!          write(*,'("Ebottom desum=",2es16.8)') desum
          call cpu_time(finish)
          gftime = gftime + finish - start
          call cpu_time(start)
!          call add_susc_block(ia2,ja2,suscblock,suscwork,chi0eb)
          call add_susc_block(ia2,ja2,suscblock,suscwork,kssusc0)
          call cpu_time(finish)
          basistime = basistime + finish - start
!       ----------------------------------------------------------------
        end if
!       ----------------------------------------------------------------
        if (lnonanalytic .and. nwef > 0) then
!       ----------------------------------------------------------------
          call cpu_time(start)
!         Integration on real axis around Efermi
          suscblock = czero; desum = czero
!         ++++++++++++++++
!         loop over panels
          do jw=1,iw
!         ++++++++++++++++
!           panel for Gij(E)
            iew0 = nescf + 2*(jw-1)*(nescf+nweb+nwef) + 2*(nescf+nweb) + nwef
!           panel for Gji(E-w)
            jew0 = nescf + 2*(iw-jw)*(nescf+nweb+nwef) + 2*(nescf+nweb) + nwef
            do ie=0,nwef  ! energy
!              write(*,*) "ie", ie
!             find energies in storage
              iew = iew0 + ie
              jew = jew0 - ie
              de = desusc(iew)
!             handle efermi (not repeated in energy mesh)
              if (ie == 0) then
                if (jw == 1) then
                  iew = nescf
                  de = desusc(iew0+nwef)
                else
                  iew = iew0 - 2*(nescf+nweb) - nwef
                  de = desusc(iew)
                end if
              end if
              if (ie == nwef) then
                jew = jew0 - 2*(nescf+nweb+nwef)
                if (jw == iw) jew = nescf
              end if
              desum = desum + de
!              write(*,'(4i8,8es16.8)') iew0, jew0, iew, jew, esusc(iew), esusc(jew), esusc(iew) - esusc(jew), de
!             get the projected GFs
!             E + i0
              call projected_gf(iew,ia,ja,gfij0,onsite,struct)
              if (lrot) call local_frame(ia,ja,magdir,gfij0,gf)
!             E - w - i0
              call projected_gf(jew,ia,ja,gfjiw,onsite,struct)
              if (lrot) call local_frame(ia,ja,magdir,gfjiw,gf)
              gfjiw = conjg(transpose(gfjiw))
              do j=1,jq1-jq0
                jq = jq0 + j
                i3 = i2almsbgf(:,jq)
                q3 = i3(1); q4 = i3(2)!; ja = i3(3)
                do i=1,iq1-iq0
                  iq = iq0 + i
                  i3 = i2almsbgf(:,iq)
                  q1 = i3(1); q2 = i3(2)!; ia = i3(3)
                  norm = de*gfij0(q2,q3)*gfjiw(q4,q1)/i2pi
                  suscblock(i,j) = suscblock(i,j) + norm
                end do
              end do
              gfij0up = 1.d0/(esusc(iew) - e0up); gfij0dn = 1.d0/(esusc(iew) - e0dn)
              gfjiwup = 1.d0/(esusc(jew) - e0up); gfjiwdn = 1.d0/(esusc(jew) - e0dn)
              gfjiwup = conjg(gfjiwup); gfjiwdn = conjg(gfjiwdn)
!              gfij0up = gfij0(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfij0dn = gfij0(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
!              gfjiwup = gfjiw(lmsb2i(1,5,2,ia),lmsb2i(1,5,2,ia)); gfjiwdn = gfjiw(lmsb2i(1,5,1,ia),lmsb2i(1,5,1,ia))
              chi0eb_updn = chi0eb_updn + de*gfij0up*gfjiwdn/i2pi
              chi0eb_dnup = chi0eb_dnup + de*gfij0dn*gfjiwup/i2pi
            end do  ! energy
!         ++++++
          end do
!         ++++++
!          write(*,'("EFermi  desum=",2es16.8)') desum
          call cpu_time(finish)
          gftime = gftime + finish - start
          call cpu_time(start)
!          call add_susc_block(ia2,ja2,suscblock,suscwork,chi0ef)
          call add_susc_block(ia2,ja2,suscblock,suscwork,kssusc0)
          call cpu_time(finish)
          basistime = basistime + finish - start
!       ----------------------------------------------------------------
        end if
!       ----------------------------------------------------------------
!       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      end do            ! atom i
    end do              ! atom j
!   Add contributions to KS susc
!    if (lanalytic) kssusc0 = kssusc0 + chi0escf
!    if (nonanalytic .and. nweb > 0) kssusc0 = kssusc0 + chi0eb
!    if (nonanalytic .and. nwef > 0) kssusc0 = kssusc0 + chi0ef
    chi0_updn = chi0escf_updn + chi0eb_updn + chi0ef_updn
    chi0_dnup = chi0escf_dnup + chi0eb_dnup + chi0ef_dnup
    write(555,'(6es16.8)') zw(iw), chi0_updn, chi0_dnup
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    where (abs(kssusc0) < susctol) kssusc0 = 0.d0
!   ------------------------------------------------------------------
    if (lenhanced) then
      call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc0,ndensum,kernel,ndensum,czero,denominator,ndensum)
      do iq=1,ndensum
        denominator(iq,iq) = denominator(iq,iq) + 1.d0
      end do
      call zgesv(ndensum,ndensum,denominator,ndensum,ipiv,kssusc0,ndensum,info)
      if (info /= 0) stop 'dyn_susc_real2: failure in zgesv'
!      call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc0,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
!      write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
!      write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
!      if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
!      kssusc0 = x
    end if
!   ------------------------------------------------------------------
!   Symmetrize
!    call symmetrize(nalmsb,kssusc0,susctol)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   multipoles
    suscylm = czero
    do jq=1,ndensum
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
      do iq=1,ndensum
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
        suscylm(is,js,ilm,jlm,ia2,ja2) = suscylm(is,js,ilm,jlm,ia2,ja2) + suscnorm(iq)*kssusc0(iq,jq)*suscnorm(jq)
      end do
    end do
!   Temporary xc kernel
!    suscylm(1,2,1,1,1,1) = suscylm(1,2,1,1,1,1)/(1.d0 - uparam*suscylm(1,2,1,1,1,1))
!    suscylm(2,1,1,1,1,1) = suscylm(2,1,1,1,1,1)/(1.d0 - uparam*suscylm(2,1,1,1,1,1))
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
              block(i,j) = block(i,j) + ds2c(i,i2is(1,is),i2is(2,is))*suscylm(is,js,ilm,jlm,ia2,ja2)*pc2s(i2is(1,js),i2is(2,js),j)
            end do
            end do
          end do
          end do
          suscylm(:,:,ilm,jlm,ia2,ja2) = block
        end do
        end do
      end do
      end do
    end if
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Final filtering (global)
!    maxelem = maxval(abs(suscylm))
!    where (abs(suscylm) < susctol*maxelem) suscylm = 0.d0
!    maxelem = maxval(abs(suscy00))
!    where (abs(suscy00) < susctol*maxelem) suscy00 = 0.d0
!    write(*,'(/,"KS susc ia, ja, im, il, jm, jl, susc ylm")')
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do ia2=1,nasusc2
        ia = iasusc2(ia2)
        do jlm=1,lmmax0
          do ilm=1,lmmax0
            norm = sum(abs(suscylm(:,:,ilm,jlm,ia2,ja2)))
            if (abs(norm) > susctol) then
              do j=1,4
!                write(*,'(6i4,8f16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscylm(i,j,ilm,jlm,ia2,ja2),i=1,4)
              end do
            end if
          end do
        end do
      end do
    end do
    write(iofile,'(1000es16.8)') real(zw(iw)), ((((suscylm(i,j,1,1,ia2,ja2),i=1,4),j=1,4),ia2=1,nasusc2),ja2=1,nasusc2)
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ******
  end do
! ******
  close(iofile)
  deallocate(suscblock,suscwork,zw,gf)
!  if (lanalytic) deallocate(chi0escf)
!  if (nonanalytic .and. nweb > 0) deallocate(chi0eb)
!  if (nonanalytic .and. nwef > 0) deallocate(chi0ef)
  write(*,'(/," Dynamic KS susc chiGF time=",f10.3," s")') gftime
  write(*,'(" Dynamic KS susc basis time=",f10.3," s")') basistime
! All done!
  end subroutine dynamic_susc2
