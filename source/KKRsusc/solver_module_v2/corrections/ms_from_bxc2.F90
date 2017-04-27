  subroutine ms_from_bxc2(ie,vxcdiff,vsocz,vsocpm,msgf,msxc,mssocz,mssocpm)
! check the magnetization sum rule
  use global

  implicit none

! Current energy
  integer(kind=i4b), intent(in)    :: ie
! xc splitting
  complex(kind=c8b), intent(in)    :: vxcdiff(nbmax,nbmax,lmmax,lmmax,nasusc2)
! SOC z splitting
  complex(kind=c8b), intent(in)    :: vsocz(nbmax,nbmax,lmmax,lmmax,nasusc2)
! non-local SOC spin splitting
  complex(kind=c8b), intent(in)    :: vsocpm(nbmax,nbmax,lmmax,lmmax,nasusc2,nasusc2)
! magnetization constructed from site-diagonal block of GF
  complex(kind=c8b), intent(inout) :: msgf(nbmax,nbmax,lmmax,lmmax,nasusc2)
! contribution from the xc potential
  complex(kind=c8b), intent(inout) :: msxc(nbmax,nbmax,lmmax,lmmax,nasusc2,nasusc2)
! contribution from vsocz
  complex(kind=c8b), intent(inout) :: mssocz(nbmax,nbmax,lmmax,lmmax,nasusc2,nasusc2)
! contribution from vsocpm
  complex(kind=c8b), intent(inout) :: mssocpm(nbmax,nbmax,lmmax,lmmax,nasusc2,nasusc2)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0), tol = 1.d-7
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi), iu = (0.d0,1.d0)
  integer(kind=i4b) :: i3(3), ia, ia2, i1, ib1, ilm1, is1, i2, ib2, ilm2, is2, ja, ja2, j1, jb1, jlm1, js1, j2, jb2, jlm2, js2
  complex(kind=c8b) :: gfij(nlmsb,nlmsb), gfji(nlmsb,nlmsb), tmp, msgflm, msxclm(nasusc2), mssoczlm(nasusc2), mssocpmlm(nasusc2)
  real(kind=r8b)    :: re, im


  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      call projected_gf(ie,ia,ja,gfij,.true.,.true.)
      call projected_gf(ie,ja,ia,gfji,.true.,.true.)
!     ------------------------------------------------------------------
!     magnetization from site-diagonal GF
      if (ia == ja) then
        do i2=1,nlmsba(ia)
          i3 = i2lmsb(:,i2,ia)
          ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
          do i1=1,nlmsba(ia)
            i3 = i2lmsb(:,i1,ia)
            ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
            tmp = (conjg(descf(ie)*pauli(is1,is2,3)*gfji(i2,i1)) - descf(ie)*pauli(is2,is1,3)*gfij(i1,i2))/i2pi 
            msgf(ib1,ib2,ilm1,ilm2,ia2) = msgf(ib1,ib2,ilm1,ilm2,ia2) + tmp
          end do
        end do
      end if
!     ------------------------------------------------------------------
      do j2=1,nlmsba(ja)
      do j1=1,nlmsba(ja)
        i3 = i2lmsb(:,j1,ja)
        jb1 = i3(1); jlm1 = i3(2); js1 = i3(3)
        i3 = i2lmsb(:,j2,ja)
        jb2 = i3(1); jlm2 = i3(2); js2 = i3(3)
        do i2=1,nlmsba(ia)
        do i1=1,nlmsba(ia)
          i3 = i2lmsb(:,i1,ia)
          ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
          i3 = i2lmsb(:,i2,ia)
          ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
!         --------------------------------------------------------------
!         spin of ja (dn,up)
          if (js1 == 2 .and. js2 == 1) then
!           spin of ia (up,dn)
            if (is1 == 2 .and. is2 == 1) then
!              write(*,'("+- ",7i4)') ie, ia, ja, i1, i2, j1, j2
!           magnetization from xc splitting
              tmp = (conjg(descf(ie)*gfji(j1,i1)*vxcdiff(jb2,jb1,jlm2,jlm1,ja2)*gfij(i2,j2)) &
                      - descf(ie)*gfij(i1,j1)*vxcdiff(jb1,jb2,jlm1,jlm2,ja2)*gfji(j2,i2))/i2pi
              msxc(ib1,ib2,ilm1,ilm2,ia2,ja2) = msxc(ib1,ib2,ilm1,ilm2,ia2,ja2) + 0.5d0*tmp
!           magnetization from SOC z splitting
              tmp = (conjg(descf(ie)*gfji(j1,i1)*vsocz(jb2,jb1,jlm2,jlm1,ja2)*gfij(i2,j2)) &
                      - descf(ie)*gfij(i1,j1)*vsocz(jb1,jb2,jlm1,jlm2,ja2)*gfji(j2,i2))/i2pi
              mssocz(ib1,ib2,ilm1,ilm2,ia2,ja2) = mssocz(ib1,ib2,ilm1,ilm2,ia2,ja2) + 0.5d0*tmp
!           magnetization from non-local SOC splitting
              tmp = (conjg(descf(ie)*gfji(j1,i1)*vsocpm(jb2,jb1,jlm2,jlm1,ja2,ja2)*gfij(i2,j2)) &
                      - descf(ie)*gfij(i1,j1)*vsocpm(jb1,jb2,jlm1,jlm2,ja2,ja2)*gfji(j2,i2))/i2pi
              mssocpm(ib1,ib2,ilm1,ilm2,ia2,ja2) = mssocpm(ib1,ib2,ilm1,ilm2,ia2,ja2) + 0.5d0*tmp
            end if
          end if
!         --------------------------------------------------------------
!         spin of ja (dn,up)
          if (js1 == 1 .and. js2 == 2) then
!           spin of ia (up,dn)
            if (is1 == 1 .and. is2 == 2) then
!              write(*,'("+- ",7i4)') ie, ia, ja, i1, i2, j1, j2
!           magnetization from xc splitting
              tmp = (conjg(descf(ie)*gfji(j1,i1)*vxcdiff(jb2,jb1,jlm2,jlm1,ja2)*gfij(i2,j2)) &
                      - descf(ie)*gfij(i1,j1)*vxcdiff(jb1,jb2,jlm1,jlm2,ja2)*gfji(j2,i2))/i2pi
              msxc(ib1,ib2,ilm1,ilm2,ia2,ja2) = msxc(ib1,ib2,ilm1,ilm2,ia2,ja2) + 0.5d0*tmp
!           magnetization from SOC z splitting
              tmp = (conjg(descf(ie)*gfji(j1,i1)*vsocz(jb2,jb1,jlm2,jlm1,ja2)*gfij(i2,j2)) &
                      - descf(ie)*gfij(i1,j1)*vsocz(jb1,jb2,jlm1,jlm2,ja2)*gfji(j2,i2))/i2pi
              mssocz(ib1,ib2,ilm1,ilm2,ia2,ja2) = mssocz(ib1,ib2,ilm1,ilm2,ia2,ja2) + 0.5d0*tmp
!           magnetization from non-local SOC splitting
              tmp = (conjg(descf(ie)*gfji(j1,i1)*vsocpm(jb2,jb1,jlm2,jlm1,ja2,ja2)*gfij(i2,j2)) &
                      - descf(ie)*gfij(i1,j1)*vsocpm(jb1,jb2,jlm1,jlm2,ja2,ja2)*gfji(j2,i2))/i2pi
              mssocpm(ib1,ib2,ilm1,ilm2,ia2,ja2) = mssocpm(ib1,ib2,ilm1,ilm2,ia2,ja2) + 0.5d0*tmp
            end if
          end if
!         --------------------------------------------------------------
        end do
        end do
      end do
      end do
    end do
  end do
! Print output
  if (ie == nescf) then
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
!      do ilm2=1,lmmax
!      do ilm1=1,lmmax
!      ilm1 = ilm2
!        do ib1=1,nbmax
!        do ib2=1,nbmax
!          write(*,'(5i4,10es12.4)') ia, ib1, ib2, ilm1, ilm2, msgf(ib1,ib2,ilm1,ilm2,ia2), &
!             sum(msxc(ib1,ib2,ilm1,ilm2,ia2,:) + mssocz(ib1,ib2,ilm1,ilm2,ia2,:) + mssocpm(ib1,ib2,ilm1,ilm2,ia2,:)),  &
!             sum(msxc(ib1,ib2,ilm1,ilm2,ia2,:)), sum(mssocz(ib1,ib2,ilm1,ilm2,ia2,:)), sum(mssocpm(ib1,ib2,ilm1,ilm2,ia2,:))
!        end do
!        end do
!      end do
!      end do
      do jlm1=1,lmmax2
        msgflm = 0.d0; msxclm = 0.d0; mssoczlm = 0.d0; mssocpmlm = 0.d0
        do ilm2=1,lmmax
        do ilm1=1,lmmax
          do ib2=1,iwsusc(i2lm(2,ilm2),1,ia)
          do ib1=1,iwsusc(i2lm(2,ilm1),1,ia)
            i2 = lmsb2i(ib2,ilm2,1,ia)
            i1 = lmsb2i(ib1,ilm1,1,ia)
            msgflm    = msgflm    + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*msgf(ib1,ib2,ilm1,ilm2,ia2)
            msxclm    = msxclm    + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*msxc(ib1,ib2,ilm1,ilm2,ia2,:)
            mssoczlm  = mssoczlm  + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*mssocz(ib1,ib2,ilm1,ilm2,ia2,:)
            mssocpmlm = mssocpmlm + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*mssocpm(ib1,ib2,ilm1,ilm2,ia2,:)
          end do
          end do
        end do
        end do
        re = real(msgflm); im = aimag(msgflm)
        if (abs(re) < tol) re = 0.d0
        if (abs(im) < tol) im = 0.d0
        msgflm = cmplx(re,im)
        do ja2=1,nasusc2
          re = real(msxclm(ja2)); im = aimag(msxclm(ja2))
          if (abs(re) < tol) re = 0.d0
          if (abs(im) < tol) im = 0.d0
          msxclm(ja2) = cmplx(re,im)
          re = real(mssoczlm(ja2)); im = aimag(mssoczlm(ja2))
          if (abs(re) < tol) re = 0.d0
          if (abs(im) < tol) im = 0.d0
          mssoczlm(ja2) = cmplx(re,im)
          re = real(mssocpmlm(ja2)); im = aimag(mssocpmlm(ja2))
          if (abs(re) < tol) re = 0.d0
          if (abs(im) < tol) im = 0.d0
          mssocpmlm(ja2) = cmplx(re,im)
        end do
        write(*,'(3i4,10es12.4)') ia, i2lm(:,jlm1), msgflm, sum(msxclm+mssoczlm+mssocpmlm), sum(msxclm), sum(mssoczlm), sum(mssocpmlm)
      end do
    end do
  end if
! All done!
  end subroutine ms_from_bxc2
