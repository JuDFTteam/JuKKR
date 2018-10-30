  subroutine static_susc(onsite,struct)
! Static KS susceptibility
! Using density basis
  use global

  implicit none

  logical,           intent(in)  :: onsite, struct
! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
! -----------------------------------------------------------------
  complex(kind=c8b) :: suscylm(4,4,lmmax0,lmmax0,nasusc2,nasusc2)
  complex(kind=c8b) :: suscy00(4,4,lmmax,lmmax,nasusc2,nasusc2)
  complex(kind=c8b) :: suscden(4,lmmax,nasusc2), rhoden(ndensum), rhodenylm(4,lmmax0,nasusc2)
  complex(kind=c8b) :: gfijw(nlmsb,nlmsb), gfjiw(nlmsb,nlmsb)
  complex(kind=c8b) :: gfij0(nlmsb,nlmsb), gfji0(nlmsb,nlmsb)
  complex(kind=c8b) :: gf(nlmsb,nlmsb)
  complex(kind=c8b), allocatable :: suscblock(:,:), suscwork(:,:)
  complex(kind=c8b) :: norm, de, e, norm2, dosylm(2,lmmax0,nasusc2)
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4
  integer(kind=i4b) :: ib, jb, is, js
  integer(kind=i4b) :: i2(2), i3(3), i4(4), ie, ia, ja, ia2, ja2, jlm, ilm, j, i, klm
  integer(kind=i4b) :: iq, jq, iq0, iq1, jq0, jq1
  integer(kind=i4b) :: ipiv(ngfsum), info
  real(kind=r8b)    :: gaunti, gauntj, doublegaunt, maxelem, start, finish, gftime, basistime, re, im
  real(kind=r8b)    :: dr(nrmax)
  complex(kind=c8b) :: work(nrmax), block(4,4)
  integer(kind=i4b) :: ne, nr
  complex(kind=c8b), external :: radint
!  real(kind=r8b)    :: r(nalmsb), c(nalmsb), rwork(2*nalmsb), rcond, ferr(nalmsb), berr(nalmsb)
!  complex(kind=c8b) :: work(2*nalmsb), temp(nalmsb,nalmsb)
!  complex(kind=c8b) :: x(nalmsb,nalmsb)
!  character*1       :: equed

  suscylm = 0.d0; suscy00 = 0.d0; maxelem = 0.d0; suscden = 0.d0
  kssusc0 = 0.d0; dosylm = 0.d0
  gftime = 0.d0; basistime = 0.d0
  allocate(suscblock(ngfmax,ngfmax),suscwork(ngfmax,ngfmax))
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ja2=1,nasusc2      ! atom j
    ja = iasusc2(ja2)
    jq0 = sum(nalmsbgf(1:ja-1))
    jq1 = sum(nalmsbgf(1:ja))
    do ia2=1,nasusc2    ! atom i
      ia = iasusc2(ia2)
      iq0 = sum(nalmsbgf(1:ia-1))
      iq1 = sum(nalmsbgf(1:ia))
      suscblock = czero
!     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call cpu_time(start)
!     Integration on usual box contour
      do ie=1,nescf     ! energy
!        write(*,*) "ie", ie
        e  = escf(ie)
        de = descf(ie)
!       get the projected GFs
!       EF + i0
        call projected_gf(ie,ia,ja,gfijw,onsite,struct)
        if (lrot) call local_frame(ia,ja,magdir,gfijw,gf)
!       EF - i0
        gfjiw = conjg(transpose(gfijw))
!       EF + i0
        call projected_gf(ie,ja,ia,gfji0,onsite,struct)
        if (lrot) call local_frame(ja,ia,magdir,gfji0,gf)
!       EF - i0
        gfij0 = conjg(transpose(gfji0))
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
      end do            ! energy
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
                                                               + rgaunt(ilm,jlm,klm)*overlap(iq,jq,ia)*(gfijw(iq,jq) - gfjiw(iq,jq))/i2pi
                if (is == 2 .and. js == 2) dosylm(1,klm,ia2) = dosylm(1,klm,ia2) &
                                                               + rgaunt(ilm,jlm,klm)*overlap(iq,jq,ia)*(gfijw(iq,jq) - gfjiw(iq,jq))/i2pi
              end if
            end do
          end do
        end do
      end if
!     ------------------------------------------------------------------
      call cpu_time(finish)
      gftime = gftime + finish - start
      call cpu_time(start)
      call add_susc_block(ia2,ja2,suscblock,suscwork,kssusc0)
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
  suscylm = 0.d0
  do jq=1,ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
    do iq=1,ndensum
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
      suscylm(is,js,ilm,jlm,ia2,ja2) = suscylm(is,js,ilm,jlm,ia2,ja2) + suscnorm(iq)*kssusc0(iq,jq)*suscnorm(jq)
    end do
  end do
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
  write(*,'(/,"KS susc ia, ja, im, il, jm, jl, susc ylm")')
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm=1,1!lmmax0
        do ilm=1,1!lmmax0
          norm = sum(abs(suscylm(:,:,ilm,jlm,ia2,ja2)))
          if (abs(norm) > susctol) then
            do j=1,4
              write(*,'(6i4,8f16.8)') ia, ja, i2lm(:,ilm), i2lm(:,jlm), (suscylm(i,j,ilm,jlm,ia2,ja2),i=1,4)
            end do
          end if
        end do
      end do
    end do
  end do
! DOS at EF test
  where (abs(dosylm) < susctol) dosylm = 0.d0
  write(*,'(/,"DOS at EF sum rule: ia, dos, susc (up then dn)")')
  do ia2=1,nasusc2
    do ilm=1,lmmax0
      write(*,'(2i4,8es16.8)') ia2, ilm, real(dosylm(1,ilm,ia2)), real(sum(suscylm(3,3:4,ilm,1,ia2,:))),  &
          real(dosylm(2,ilm,ia2)), real(sum(suscylm(4,3:4,ilm,1,ia2,:)))
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
      if (ilm == 1 .and. is == 1) then
        work(1:nr) = suscbasis(1:nr,ib,ilm,is,ia2)*br(1:nr,ia)
        norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
        vlmsbden(iq) = norm*2.d0
        norm2 = norm2 + vlmsbden(iq)*suscnorm(iq)
      end if
    end do
    write(*,'("KS susc ia=",i4," potnorm=",2f16.8)') ia, norm2
  end do
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
  deallocate(suscblock,suscwork)
  write(*,'(/," Static  KS susc chiGF time=",f10.3," s")') gftime
  write(*,'(" Static  KS susc basis time=",f10.3," s")') basistime
! All done!
  end subroutine static_susc
