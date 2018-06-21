  subroutine build_kxcalda2()
! assembles the xc ALDA kernel in the density basis
! input:  kxc multipoles up to lmmax2 (as GS density)
! output: kxc multipoles up to lmmax0 (as susceptibility)
  use global

  implicit none

! -----------------------------------------------------------------
  complex(kind=c8b) :: kxcylm(4,4,lmmax0,lmmax0,nasusc2)
  complex(kind=c8b) :: kxcy00(4,4,lmmax,lmmax,nasusc2)
  real(kind=r8b),    parameter :: fourpi = 16.d0*atan(1.d0)
  integer(kind=i4b) :: s1, s2, s3, s4, is, js, nden0, nden1
  integer(kind=i4b) :: i2(2), i4(4), ie, ia, ia2, jlm, ilm, klm, ib, jb, j, i, nr, iq, jq
  real(kind=r8b)    :: dr(nrmax), maxnorm, dummy(5)
  real(kind=r8b)    :: maxelem
  complex(kind=c8b) :: work(nrmax), norm, norm2, block(4,4), kxc(nrmax)
  integer(kind=i4b) :: iqmax, jqmax
  complex(kind=c8b), external :: radint

  kxcylm = 0.d0; kxcy00 = 0.d0!; kernel = 0.d0
! use the ASA
! use the fact that the basis was constructed per l channel
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ia2=1,nasusc2   ! atom i
    ia = iasusc2(ia2)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    nr = nrpts(ia)
    dr(1:nr) = drmesh(1:nr,ia)!/rmesh(1:nr,ia)**2
    maxnorm = 0.d0
    nden0 = sum(nalmsbden(1:ia2-1))
    nden1 = sum(nalmsbden(1:ia2))
    do jq=1+nden0,nden1
      i4 = i2almsbden(:,jq)
      jb = i4(1); jlm = i4(2); js = i4(3)!; ia = i4(4)
      i2 = i2is(:,js)
      s3 = i2(1); s4 = i2(2)
      do iq=1+nden0,nden1
        i4 = i2almsbden(:,iq)
        ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
        i2 = i2is(:,is)
        s1 = i2(1); s2 = i2(2)
!       ***********************************
!        if (ib == jb .and. ilm == jlm) then
!       ***********************************
!        write(iodb,'("build_kxcalda2: indices=",7i8)') ia, is, js, ilm, jlm, ib, jb
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        norm = 0.d0
        do klm=1,lmmax0
        if (abs(rgaunt(ilm,jlm,klm)) > ylmtol) then
          kxc = 0.d0
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         convert from cartesian to spin
          do j=1,4
            do i=1,4
              kxc(1:nr) = kxc(1:nr) + pc2s(s1,s2,i)*kxclm(1:nr,klm,i,j,ia)*ds2c(j,s3,s4)*rgaunt(ilm,jlm,klm)
            end do
          end do
!          if (ib == 1 .and. jb == 1 .and. ilm == 1 .and. jlm == 1) then
!            write(*,'("kxclm=",2i4,16f8.3," | ",2f8.3)') is, js, kxclm(nr,1,:,:,ia), kxc(nr)
!          end if
!         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         matrix elements
          work(1:nr) = suscbasis(1:nr,ib,ilm,is,ia2)*suscbasis(1:nr,jb,jlm,js,ia2)*kxc(1:nr)/fourpi
          norm = norm + radint(nr,work,dr,npanat(ia),ircutat(:,ia))
        end if
        end do
!       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!        if (abs(norm) < susctol) norm = 0.d0
!       adds only the longitudinal kernel
        if (is > 2 .and. js > 2) then
          kxcylm(is,js,ilm,jlm,ia2) = kxcylm(is,js,ilm,jlm,ia2) + norm*suscnorm(iq)*suscnorm(jq)
          kernel(iq,jq) = kernel(iq,jq) + norm
        end if
        if (abs(norm) > maxnorm) then
          maxnorm = abs(norm)
          iqmax = iq; jqmax = jq
        end if
!        if (abs(norm) > susctol) write(iodb,'("ia,is,js,ilm,jlm,ib,jb,norm=",7i8,2es16.8)') ia, is, js, ilm, jlm, ib, jb, norm
!       ******
!        end if
!       ******
      end do
    end do
    write(iodb,'("build_kxcalda2: norm=",7i8,2es16.8)') ia, i2almsbden(1:3,iqmax), i2almsbden(1:3,jqmax), maxnorm
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do            ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(loutsusc) then 
    write(*,'(/,"xc kernel ia, im, il, jm, jl, kxc ylm")')
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm=1,1!lmmax0
      do ilm=1,1!lmmax0
        norm = sum(abs(kxcylm(:,:,ilm,jlm,ia2)))
        if (real(norm) > atol) then
          if (lcartesian) then
            block = 0.d0
            do j=1,4
            do i=1,4
              do js=1,4
              do is=1,4
                block(i,j) = block(i,j) + ps2c(i,i2is(1,is),i2is(2,is))*kxcylm(is,js,ilm,jlm,ia2)*dc2s(i2is(1,js),i2is(2,js),j)
              end do
              end do
            end do
            end do
            kxcylm(:,:,ilm,jlm,ia2) = block
          end if
          do i=1,4
            write(*,'(5i4,8f16.8)') ia, i2lm(:,ilm), i2lm(:,jlm), (kxcylm(i,j,ilm,jlm,ia2),j=1,4)
          end do
        end if
      end do
      end do
    end do
  end if 
!  kxcsusc = real(kxcsusc)
! All done
  end subroutine build_kxcalda2
