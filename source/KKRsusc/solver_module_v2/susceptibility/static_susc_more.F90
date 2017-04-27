  subroutine static_susc_more(onsite,struct)
! Static KS susceptibility, spin and orbital, spherical averages
! Also first and second frequency derivatives
! Using GF basis
  use global

  implicit none

  logical,           intent(in)  :: onsite, struct
! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
! -----------------------------------------------------------------
  complex(kind=c8b), allocatable :: gfijw(:,:), gfjiw(:,:), gfij0(:,:), gfji0(:,:), gf(:,:)
  complex(kind=c8b), allocatable :: gfijdef(:,:), gfjidef(:,:), gfijef(:,:), gfjief(:,:)
  complex(kind=c8b), allocatable :: gfijdeb(:,:), gfjideb(:,:), gfijeb(:,:), gfjieb(:,:)
  complex(kind=c8b), allocatable :: chi0ssw0(:,:,:,:), chi0ssw1(:,:,:,:), chi0ssw2(:,:,:,:)
  complex(kind=c8b), allocatable :: chi0sow0(:,:,:,:), chi0sow1(:,:,:,:), chi0sow2(:,:,:,:)
  complex(kind=c8b), allocatable :: chi0osw0(:,:,:,:), chi0osw1(:,:,:,:), chi0osw2(:,:,:,:)
  complex(kind=c8b), allocatable :: chi0oow0(:,:,:,:), chi0oow1(:,:,:,:), chi0oow2(:,:,:,:)
  complex(kind=c8b) :: de1, de2, e, norm0, norm1, norm2
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4
  integer(kind=i4b) :: i2(2), i3(3), i4(4), ia, ja, ia2, ja2, i, j
  integer(kind=i4b) :: iq, jq, iq0, iq1, jq0, jq1
  real(kind=r8b)    :: maxelem0, maxelem1, maxelem2, start, finish, gftime, re, im
  real(kind=r8b)    :: rea, ima, reb, imb, ram
  integer(kind=i4b) :: ne, ip, ie, je, ie0
  integer(kind=i4b) :: nepan, netot
  integer(kind=i4b) :: ipan(nescf)
  complex(kind=c8b) :: ea(nescf), eb(nescf), z(nescf), w(nescf), w0(nescf), w1(nescf)
  complex(kind=c8b) :: w00(nescf,nescf), w10(nescf,nescf), w01(nescf,nescf), w11(nescf,nescf) 

  ram = 16*(13*nlmsb**2 + 12*16*nasusc2**2)
  ram = ram/1024.d0**3
  write(*,'("RAM for static_susc_more: ",f8.3," GB")') ram
  allocate(gf(nlmsb,nlmsb),gfijw(nlmsb,nlmsb),gfjiw(nlmsb,nlmsb),gfij0(nlmsb,nlmsb),gfji0(nlmsb,nlmsb))
  allocate(gfijdef(nlmsb,nlmsb),gfjidef(nlmsb,nlmsb),gfijef(nlmsb,nlmsb),gfjief(nlmsb,nlmsb))
  allocate(gfijdeb(nlmsb,nlmsb),gfjideb(nlmsb,nlmsb),gfijeb(nlmsb,nlmsb),gfjieb(nlmsb,nlmsb))
  allocate(chi0ssw0(4,4,nasusc2,nasusc2),chi0ssw1(4,4,nasusc2,nasusc2),chi0ssw2(4,4,nasusc2,nasusc2))
  allocate(chi0sow0(4,4,nasusc2,nasusc2),chi0sow1(4,4,nasusc2,nasusc2),chi0sow2(4,4,nasusc2,nasusc2))
  allocate(chi0osw0(4,4,nasusc2,nasusc2),chi0osw1(4,4,nasusc2,nasusc2),chi0osw2(4,4,nasusc2,nasusc2))
  allocate(chi0oow0(4,4,nasusc2,nasusc2),chi0oow1(4,4,nasusc2,nasusc2),chi0oow2(4,4,nasusc2,nasusc2))
  chi0ssw0 = 0.d0; chi0ssw1 = 0.d0; chi0ssw2 = 0.d0
  chi0sow0 = 0.d0; chi0sow1 = 0.d0; chi0sow2 = 0.d0
  chi0osw0 = 0.d0; chi0osw1 = 0.d0; chi0osw2 = 0.d0
  chi0oow0 = 0.d0; chi0oow1 = 0.d0; chi0oow2 = 0.d0
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  open(file='meshpanels.dat',unit=iofile,status='old')
  read(iofile,*) nepan, netot
  if (netot /= nescf) stop 'static_susc_more: netot /= nescf'
  do ie=1,nepan
    read(iofile,*) rea, ima, reb, imb, ipan(ie)
    write(*,'(4f10.6,i8)') rea, ima, reb, imb, ipan(ie)
    ea(ie) = cmplx(rea,ima); eb(ie) = cmplx(reb,imb)
  end do
  if (netot /= sum(ipan(1:nepan))) stop 'static_susc2: netot /= sum(ipan)'
  close(iofile)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  call cpu_time(start)
  do ja2=1,nasusc2      ! atom j
    ja = iasusc2(ja2)
    do ia2=1,nasusc2    ! atom i
      ia = iasusc2(ia2)
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
          if (lrot) call local_frame(ja,ia,magdir,gfji0,gf)
!         Gji(E - i0) = Gij(E + i0)^\dagger
          call projected_gf(ie0+je,ia,ja,gfjiw,onsite,struct)
          if (lrot) call local_frame(ia,ja,magdir,gfjiw,gf)
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
            if (lrot) call local_frame(ia,ja,magdir,gfijw,gf)
!           Gij(E - i0) = Gji(E + i0)^\dagger
            call projected_gf(ie0+ie,ja,ia,gfij0,onsite,struct)
            if (lrot) call local_frame(ja,ia,magdir,gfij0,gf)
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
!           Static susc
            do q4=1,nlmsba(ja)
              i3 = i2lmsb(:,q4,ja)
              p4 = i3(1); lm4 = i3(2); s4 = i3(3)
              i2 = i2lm(:,lm4)
              m4 = i2(1); l4 = i2(2)
            do q3=1,nlmsba(ja)
              i3 = i2lmsb(:,q3,ja)
              p3 = i3(1); lm3 = i3(2); s3 = i3(3)
              i2 = i2lm(:,lm3)
              m3 = i2(1); l3 = i2(2)
            do q2=1,nlmsba(ia)
              i3 = i2lmsb(:,q2,ia)
              p2 = i3(1); lm2 = i3(2); s2 = i3(3)
              i2 = i2lm(:,lm2)
              m2 = i2(1); l2 = i2(2)
            do q1=1,nlmsba(ia)
              i3 = i2lmsb(:,q1,ia)
              p1 = i3(1); lm1 = i3(2); s1 = i3(3)
              i2 = i2lm(:,lm1)
              m1 = i2(1); l1 = i2(2)
!             -----------------------------------------------------------
!             matrix elements for the susceptibility
!             +1/i2pi \int dz^* Tr A Gji(z^*)^\dagger B Gij(z^*)^\dagger - 1/i2pi \int dz Tr A Gij(z) B Gji(z)
              norm0 = (conjg(w00(ie,je))*gfij0(q2,q3)*gfjiw(q4,q1) - w00(ie,je)*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
!             -1/i2pi \int dz^* Tr A Gji(z^*)^\dagger B d/dz^* Gij(z^*)^\dagger - 1/i2pi \int dz Tr A d/dz Gij(z) A Gji(z)
              norm1 = -(conjg(w01(ie,je))*gfij0(q2,q3)*gfjiw(q4,q1) + w10(ie,je)*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
!             -1/i2pi \int dz^* Tr A d/dz^* Gji(z^*)^\dagger B d/dz^* Gij(z^*)^\dagger + 1/i2pi \int dz A d/dz Gij(z) B d/dz Gji(z)
              norm2 = -(conjg(w11(ie,je))*gfij0(q2,q3)*gfjiw(q4,q1) - w11(ie,je)*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
!             -----------------------------------------------------------
!             operators are diagonal in the radial basis functions
              if (p1 == p2 .and. p3 == p4) then
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             spin-spin susceptibility
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if (lm1 == lm2 .and. lm3 == lm4) then
              do j=1,4
                do i=1,4
                  chi0ssw0(i,j,ia2,ja2) = chi0ssw0(i,j,ia2,ja2) + norm0*pauli(s1,s2,i)*pauli(s3,s4,j)
                  chi0ssw1(i,j,ia2,ja2) = chi0ssw1(i,j,ia2,ja2) + norm1*pauli(s1,s2,i)*pauli(s3,s4,j)
                  chi0ssw2(i,j,ia2,ja2) = chi0ssw2(i,j,ia2,ja2) + norm2*pauli(s1,s2,i)*pauli(s3,s4,j)
                end do
              end do
              end if
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             spin-orbital susceptibility
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if (lm1 == lm2 .and. s3 == s4) then
              do j=1,3
                do i=1,4
                  chi0sow0(i,j,ia2,ja2) = chi0sow0(i,j,ia2,ja2) + norm0*pauli(s1,s2,i)*lorb(lm3,lm4,j)
                  chi0sow1(i,j,ia2,ja2) = chi0sow1(i,j,ia2,ja2) + norm1*pauli(s1,s2,i)*lorb(lm3,lm4,j)
                  chi0sow2(i,j,ia2,ja2) = chi0sow2(i,j,ia2,ja2) + norm2*pauli(s1,s2,i)*lorb(lm3,lm4,j)
                end do
              end do
              end if
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             orbital-spin susceptibility
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if (s1 == s2 .and. lm3 == lm4) then
              do j=1,4
                do i=1,3
                  chi0osw0(i,j,ia2,ja2) = chi0osw0(i,j,ia2,ja2) + norm0*lorb(lm1,lm2,i)*pauli(s3,s4,j)
                  chi0osw1(i,j,ia2,ja2) = chi0osw1(i,j,ia2,ja2) + norm1*lorb(lm1,lm2,i)*pauli(s3,s4,j)
                  chi0osw2(i,j,ia2,ja2) = chi0osw2(i,j,ia2,ja2) + norm2*lorb(lm1,lm2,i)*pauli(s3,s4,j)
                end do
              end do
              end if
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!             orbital-orbital susceptibility
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if (s1 == s2 .and. s3 == s4) then
              do j=1,3
                do i=1,3
                  chi0oow0(i,j,ia2,ja2) = chi0oow0(i,j,ia2,ja2) + norm0*lorb(lm1,lm2,i)*lorb(lm3,lm4,j)
                  chi0oow1(i,j,ia2,ja2) = chi0oow1(i,j,ia2,ja2) + norm1*lorb(lm1,lm2,i)*lorb(lm3,lm4,j)
                  chi0oow2(i,j,ia2,ja2) = chi0oow2(i,j,ia2,ja2) + norm2*lorb(lm1,lm2,i)*lorb(lm3,lm4,j)
                end do
              end do
              end if
!             ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              end if
!             -----------------------------------------------------------
            end do
            end do
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
      do q4=1,nlmsba(ja)
        i3 = i2lmsb(:,q4,ja)
        p4 = i3(1); lm4 = i3(2); s4 = i3(3)
        i2 = i2lm(:,lm4)
        m4 = i2(1); l4 = i2(2)
      do q3=1,nlmsba(ja)
        i3 = i2lmsb(:,q3,ja)
        p3 = i3(1); lm3 = i3(2); s3 = i3(3)
        i2 = i2lm(:,lm3)
        m3 = i2(1); l3 = i2(2)
      do q2=1,nlmsba(ia)
        i3 = i2lmsb(:,q2,ia)
        p2 = i3(1); lm2 = i3(2); s2 = i3(3)
        i2 = i2lm(:,lm2)
        m2 = i2(1); l2 = i2(2)
      do q1=1,nlmsba(ia)
        i3 = i2lmsb(:,q1,ia)
        p1 = i3(1); lm1 = i3(2); s1 = i3(3)
        i2 = i2lm(:,lm1)
        m1 = i2(1); l1 = i2(2)
!       Nothing for static susc
        norm0 = czero
!       1/i2pi Gij(EF+i0) Gij(EF+i0)^\dagger - 1/i2pi Gij(EB+i0) Gij(EB+i0)^\dagger
        norm1 = czero
        norm1 = norm1 + (gfijef(q2,q3)*conjg(gfijef(q1,q4)) - gfijeb(q2,q3)*conjg(gfijeb(q1,q4)))/i2pi
!       write something
        norm2 = czero
        norm2 = norm2 - gfijdef(q2,q3)*(gfjief(q4,q1) - conjg(gfijef(q1,q4)))/i2pi
        norm2 = norm2 + gfijdeb(q2,q3)*(gfjieb(q4,q1) - conjg(gfijeb(q1,q4)))/i2pi
        norm2 = norm2 - (gfijef(q2,q3) - conjg(gfjief(q3,q2)))*conjg(gfjidef(q1,q4))/i2pi
        norm2 = norm2 + (gfijeb(q2,q3) - conjg(gfjieb(q3,q2)))*conjg(gfjideb(q1,q4))/i2pi
!       -----------------------------------------------------------
!       operators are diagonal in the radial basis functions
        if (p1 == p2 .and. p3 == p4) then
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       spin-spin susceptibility
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (lm1 == lm2 .and. lm3 == lm4) then
        do j=1,4
          do i=1,4
            chi0ssw0(i,j,ia2,ja2) = chi0ssw0(i,j,ia2,ja2) + norm0*pauli(s1,s2,i)*pauli(s3,s4,j)
            chi0ssw1(i,j,ia2,ja2) = chi0ssw1(i,j,ia2,ja2) + norm1*pauli(s1,s2,i)*pauli(s3,s4,j)
            chi0ssw2(i,j,ia2,ja2) = chi0ssw2(i,j,ia2,ja2) + norm2*pauli(s1,s2,i)*pauli(s3,s4,j)
          end do
        end do
        end if
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       spin-orbital susceptibility
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (lm1 == lm2 .and. s3 == s4) then
        do j=1,3
          do i=1,4
            chi0sow0(i,j,ia2,ja2) = chi0sow0(i,j,ia2,ja2) + norm0*pauli(s1,s2,i)*lorb(lm3,lm4,j)
            chi0sow1(i,j,ia2,ja2) = chi0sow1(i,j,ia2,ja2) + norm1*pauli(s1,s2,i)*lorb(lm3,lm4,j)
            chi0sow2(i,j,ia2,ja2) = chi0sow2(i,j,ia2,ja2) + norm2*pauli(s1,s2,i)*lorb(lm3,lm4,j)
          end do
        end do
        end if
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       orbital-spin susceptibility
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (s1 == s2 .and. lm3 == lm4) then
        do j=1,4
          do i=1,3
            chi0osw0(i,j,ia2,ja2) = chi0osw0(i,j,ia2,ja2) + norm0*lorb(lm1,lm2,i)*pauli(s3,s4,j)
            chi0osw1(i,j,ia2,ja2) = chi0osw1(i,j,ia2,ja2) + norm1*lorb(lm1,lm2,i)*pauli(s3,s4,j)
            chi0osw2(i,j,ia2,ja2) = chi0osw2(i,j,ia2,ja2) + norm2*lorb(lm1,lm2,i)*pauli(s3,s4,j)
          end do
        end do
        end if
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       orbital-orbital susceptibility
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if (s1 == s2 .and. s3 == s4) then
        do j=1,3
          do i=1,3
            chi0oow0(i,j,ia2,ja2) = chi0oow0(i,j,ia2,ja2) + norm0*lorb(lm1,lm2,i)*lorb(lm3,lm4,j)
            chi0oow1(i,j,ia2,ja2) = chi0oow1(i,j,ia2,ja2) + norm1*lorb(lm1,lm2,i)*lorb(lm3,lm4,j)
            chi0oow2(i,j,ia2,ja2) = chi0oow2(i,j,ia2,ja2) + norm2*lorb(lm1,lm2,i)*lorb(lm3,lm4,j)
          end do
        end do
        end if
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        end if
!       -----------------------------------------------------------
      end do
      end do
      end do
      end do
!     ------------------------------------------------------------------
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end do            ! atom i
  end do              ! atom j
  call cpu_time(finish)
  gftime = finish - start
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*,'(/,"chi0ssw0 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0ssw0(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0sow0 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0sow0(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0osw0 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0osw0(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0oow0 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0oow0(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*,'(/,"chi0ssw1 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0ssw1(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0sow1 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0sow1(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0osw1 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0osw1(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0oow1 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0oow1(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(*,'(/,"chi0ssw2 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0ssw2(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0sow2 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0sow2(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0osw2 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0osw2(i,j,ia2,ja2),j=1,4)
      end do
    end do
  end do
  write(*,'(/,"chi0oow2 ia, ja, susc")')
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja2=1,nasusc2
      ja = iasusc2(ja2)
      do i=1,4
        write(*,'(2i4,8es16.8)') ia, ja, (chi0oow2(i,j,ia2,ja2),j=1,3)
      end do
    end do
  end do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  deallocate(chi0ssw0,chi0ssw1,chi0ssw2,chi0sow0,chi0sow1,chi0sow2,chi0osw0,chi0osw1,chi0osw2,chi0oow0,chi0oow1,chi0oow2)
  deallocate(gf,gfijw,gfjiw,gfij0,gfji0,gfijdef,gfjidef,gfijef,gfjief,gfijdeb,gfjideb,gfijeb,gfjieb)
  write(*,'(/," Static  KS susc chiGF time=",f10.3," s")') gftime
! All done!
  end subroutine static_susc_more
