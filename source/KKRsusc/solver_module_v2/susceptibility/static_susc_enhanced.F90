  subroutine static_susc_enhanced()
! Static KS susceptibility including kernel
! Using density basis
  use global

  implicit none

! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
! -----------------------------------------------------------------
  complex(kind=c8b) :: suscylm(4,4,lmmax0,lmmax0,nasusc2,nasusc2)
  complex(kind=c8b) :: suscy00(4,4,lmmax,lmmax,nasusc2,nasusc2)
  complex(kind=c8b) :: norm, de, e, norm2
  integer(kind=i4b) :: ib, jb, is, js
  integer(kind=i4b) :: i2(2), i3(3), i4(4), ie, ia, ja, ia2, ja2, jlm, ilm, j, i
  integer(kind=i4b) :: iq, jq, iq0, iq1, jq0, jq1
  integer(kind=i4b) :: ipiv(ngfsum), info
  real(kind=r8b)    :: gaunti, gauntj, doublegaunt, maxelem, start, finish, re, im
  real(kind=r8b)    :: dr(nrmax)
  complex(kind=c8b) :: work(nrmax), block(4,4)
  integer(kind=i4b) :: ne, nr
  complex(kind=c8b), external :: radint
!  real(kind=r8b)    :: r(nalmsb), c(nalmsb), rwork(2*nalmsb), rcond, ferr(nalmsb), berr(nalmsb)
!  complex(kind=c8b) :: work(2*nalmsb), temp(nalmsb,nalmsb)
!  complex(kind=c8b) :: x(nalmsb,nalmsb)
!  character*1       :: equed

  call cpu_time(start)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  call zcopy(ndensum*ndensum,kssusc0,1,kssusc,1)
  call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc0,ndensum,kernel,ndensum,czero,denominator,ndensum)
  do iq=1,ndensum
    denominator(iq,iq) = denominator(iq,iq) + cone
  end do
  call zgesv(ndensum,ndensum,denominator,ndensum,ipiv,kssusc,ndensum,info)
  if (info /= 0) stop 'dyn_susc_real2: failure in zgesv'
! call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc0,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
! write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
! write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
! if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
! kssusc0 = x
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! multipoles
  suscylm = czero
  do jq=1,ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja2 = i4(4)
    do iq=1,ndensum
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3); ia2 = i4(4)
      suscylm(is,js,ilm,jlm,ia2,ja2) = suscylm(is,js,ilm,jlm,ia2,ja2) + suscnorm(iq)*kssusc(iq,jq)*suscnorm(jq)
    end do
  end do
! ----------------------------------------------------------------------
  write(*,'(/,"long magnetization+charge")')
  do ia2=1,nasusc2
    do ilm=1,lmmax0
      write(*,'(2i4,4es16.8)') ia2, ilm, sum(suscylm(3,3,ilm,1,ia2,:) + suscylm(4,4,ilm,1,ia2,:) - suscylm(3,4,ilm,1,ia2,:) - suscylm(4,3,ilm,1,ia2,:)), &
          sum(suscylm(3,3,ilm,1,ia2,:) - suscylm(4,4,ilm,1,ia2,:) - suscylm(3,4,ilm,1,ia2,:) + suscylm(4,3,ilm,1,ia2,:))
    end do
  end do
! ----------------------------------------------------------------------
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
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  call cpu_time(finish)
  write(*,'(" Static  KS susc time=",f10.3," s")') finish - start
! All done!
  end subroutine static_susc_enhanced
