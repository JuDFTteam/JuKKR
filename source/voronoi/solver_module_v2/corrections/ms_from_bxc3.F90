  subroutine ms_from_bxc3(ie,vxcdiff,vsocz,vsocpm,magdir1,msgf,msxc,mssocz,mssocpm)
! check the magnetization sum rule
! the ordering of the labels on the product of GFs was changed to agree with static_susc
! Bext is included in vsocz and vsocpm (check pot_correction)
  use global

  implicit none

! Current energy
  integer(kind=i4b), intent(in)    :: ie
! xc splitting
  complex(kind=c8b), intent(in)    :: vxcdiff(nbmax,nbmax,lmmax,lmmax,nasusc)
! SOC z splitting
  complex(kind=c8b), intent(in)    :: vsocz(nbmax,nbmax,lmmax,lmmax,nasusc)
! non-local SOC spin splitting
  complex(kind=c8b), intent(in)    :: vsocpm(nbmax,nbmax,lmmax,lmmax,nasusc,nasusc)
! spin quantization axes
  real(kind=r8b),    intent(in)    :: magdir1(3,nasusc)
! magnetization constructed from site-diagonal block of GF
  complex(kind=c8b), intent(inout) :: msgf(nbmax,nbmax,lmmax,lmmax,nasusc2)
! contribution from the xc potential
  complex(kind=c8b), intent(inout) :: msxc(nbmax,nbmax,lmmax,lmmax,nasusc2)
! contribution from vsocz
  complex(kind=c8b), intent(inout) :: mssocz(nbmax,nbmax,lmmax,lmmax,nasusc2)
! contribution from vsocpm
  complex(kind=c8b), intent(inout) :: mssocpm(nbmax,nbmax,lmmax,lmmax,nasusc2)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0), tol = 1.d-7
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi), iu = (0.d0,1.d0)
  integer(kind=i4b) :: i3(3), ia, ia2, i1, ib1, ilm1, is1, i2, ib2, ilm2, is2, ja, ja2, j1, jb1, jlm1, js1, j2, jb2, jlm2, js2, ka, ka2
  complex(kind=c8b) :: gfijr(nlmsb,nlmsb), gfkir(nlmsb,nlmsb)
  complex(kind=c8b) :: gfija(nlmsb,nlmsb), gfkia(nlmsb,nlmsb), gf(nlmsb,nlmsb)
  complex(kind=c8b) :: tmp, msgflm, msxclm, mssoczlm, mssocpmlm
  real(kind=r8b)    :: re, im


  do ia2=1,nasusc2
    ia = iasusc2(ia2)
    do ja=1,nasusc
!    if (isoc(ja) /= 0 .or. ibfield(ja) /= 0 .or. ikxc(ja) /= 0) then
!     Gij(E+i0)
      call projected_gf(ie,ia,ja,gfijr,.true.,.true.)
      if (lrot) call local_frame(ia,ja,magdir1,gfijr,gf)
!     Gij(E-i0) = Gji(E+i0)^\dagger
      call projected_gf(ie,ja,ia,gfija,.true.,.true.)
      if (lrot) call local_frame(ja,ia,magdir1,gfija,gf)
      gfija = conjg(transpose(gfija))
!     magnetization from site-diagonal GF
      if (ia == ja) then
        do i2=1,nlmsba(ia)
          i3 = i2lmsb(:,i2,ia)
          ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
          do i1=1,nlmsba(ia)
            i3 = i2lmsb(:,i1,ia)
            ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
            tmp = pauli(is1,is2,3)*(conjg(descf(ie))*gfija(i2,i1) - descf(ie)*gfijr(i2,i1))/i2pi 
            msgf(ib1,ib2,ilm1,ilm2,ia2) = msgf(ib1,ib2,ilm1,ilm2,ia2) + tmp
          end do
        end do
      end if
      do ka=1,nasusc
!      if (isoc(ka) == 0 .and. ibfield(ka) == 0 .and. ja /= ka) cycle
!       Gki(E+i0)
        call projected_gf(ie,ka,ia,gfkir,.true.,.true.)
        if (lrot) call local_frame(ka,ia,magdir1,gfkir,gf)
!       Gki(E-i0) = Gik(E+i0)^\dagger
        call projected_gf(ie,ia,ka,gfkia,.true.,.true.)
        if (lrot) call local_frame(ia,ka,magdir1,gfkia,gf)
        gfkia = conjg(transpose(gfkia))
!       ------------------------------------------------------------------
        do j2=1,nlmsba(ka)
        do j1=1,nlmsba(ja)
          i3 = i2lmsb(:,j1,ja)
          jb1 = i3(1); jlm1 = i3(2); js1 = i3(3)
          i3 = i2lmsb(:,j2,ka)
          jb2 = i3(1); jlm2 = i3(2); js2 = i3(3)
          do i2=1,nlmsba(ia)
          do i1=1,nlmsba(ia)
            i3 = i2lmsb(:,i1,ia)
            ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
            i3 = i2lmsb(:,i2,ia)
            ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
!           --------------------------------------------------------------
!           spin of ja (dn,up)
            if (js1 == 2 .and. js2 == 1) then
!             spin of ia (up,dn)
              if (is1 == 1 .and. is2 == 2) then
!                write(*,'("+- ",7i4)') ie, ia, ja, i1, i2, j1, j2
!               ----------------------------------------------------------
                if (ja == ka) then
!                 magnetization from xc splitting
                  tmp = (conjg(descf(ie))*gfija(i2,j1)*conjg(vxcdiff(jb2,jb1,jlm2,jlm1,ja))*gfkia(j2,i1) &
                              - descf(ie)*gfijr(i2,j1)*vxcdiff(jb1,jb2,jlm1,jlm2,ja)*gfkir(j2,i1))/i2pi
!                 if Bxc is an active kernel
                  if (ikxc(ja) == 1 .or. ikxc(ja) == 3) then
                    msxc(ib1,ib2,ilm1,ilm2,ia2) = msxc(ib1,ib2,ilm1,ilm2,ia2) + 0.5d0*tmp
                  else
                    mssocz(ib1,ib2,ilm1,ilm2,ia2) = mssocz(ib1,ib2,ilm1,ilm2,ia2) + 0.5d0*tmp
                  end if
!                 magnetization from SOC z splitting
                  tmp = (conjg(descf(ie))*gfija(i2,j1)*conjg(vsocz(jb2,jb1,jlm2,jlm1,ja))*gfkia(j2,i1) &
                              - descf(ie)*gfijr(i2,j1)*vsocz(jb1,jb2,jlm1,jlm2,ja)*gfkir(j2,i1))/i2pi
                  mssocz(ib1,ib2,ilm1,ilm2,ia2) = mssocz(ib1,ib2,ilm1,ilm2,ia2) + 0.5d0*tmp
                end if
!               ----------------------------------------------------------
!               magnetization from non-local SOC splitting
                tmp = (conjg(descf(ie))*gfija(i2,j1)*conjg(vsocpm(jb2,jb1,jlm2,jlm1,ka,ja))*gfkia(j2,i1) &
                            - descf(ie)*gfijr(i2,j1)*vsocpm(jb1,jb2,jlm1,jlm2,ja,ka)*gfkir(j2,i1))/i2pi
                mssocpm(ib1,ib2,ilm1,ilm2,ia2) = mssocpm(ib1,ib2,ilm1,ilm2,ia2) + 0.5d0*tmp
!               ----------------------------------------------------------
              end if
            end if
!           ------------------------------------------------------------
!           spin of ja (up,dn)
            if (js1 == 1 .and. js2 == 2) then
!             spin of ia (up,dn)
              if (is1 == 2 .and. is2 == 1) then
!                write(*,'("+- ",7i4)') ie, ia, ja, i1, i2, j1, j2
!               ----------------------------------------------------------
                if (ja == ka) then
!               magnetization from xc splitting
                  tmp = (conjg(descf(ie))*gfija(i2,j1)*conjg(vxcdiff(jb2,jb1,jlm2,jlm1,ja))*gfkia(j2,i1) &
                              - descf(ie)*gfijr(i2,j1)*vxcdiff(jb1,jb2,jlm1,jlm2,ja)*gfkir(j2,i1))/i2pi
!                 if Bxc is an active kernel
                  if (ikxc(ja) == 1 .or. ikxc(ja) == 3) then
                    msxc(ib1,ib2,ilm1,ilm2,ia2) = msxc(ib1,ib2,ilm1,ilm2,ia2) + 0.5d0*tmp
                  else
                    mssocz(ib1,ib2,ilm1,ilm2,ia2) = mssocz(ib1,ib2,ilm1,ilm2,ia2) + 0.5d0*tmp
                  end if
!               magnetization from SOC z splitting
                  tmp = (conjg(descf(ie))*gfija(i2,j1)*conjg(vsocz(jb2,jb1,jlm2,jlm1,ja))*gfkia(j2,i1) &
                              - descf(ie)*gfijr(i2,j1)*vsocz(jb1,jb2,jlm1,jlm2,ja)*gfkir(j2,i1))/i2pi
                  mssocz(ib1,ib2,ilm1,ilm2,ia2) = mssocz(ib1,ib2,ilm1,ilm2,ia2) + 0.5d0*tmp
                end if
!               ----------------------------------------------------------
!               magnetization from non-local SOC splitting
                tmp = (conjg(descf(ie))*gfija(i2,j1)*conjg(vsocpm(jb2,jb1,jlm2,jlm1,ka,ja))*gfkia(j2,i1) &
                            - descf(ie)*gfijr(i2,j1)*vsocpm(jb1,jb2,jlm1,jlm2,ja,ka)*gfkir(j2,i1))/i2pi
                mssocpm(ib1,ib2,ilm1,ilm2,ia2) = mssocpm(ib1,ib2,ilm1,ilm2,ia2) + 0.5d0*tmp
!               ----------------------------------------------------------
              end if
            end if
!           ------------------------------------------------------------
          end do
          end do
        end do
        end do
      end do
!    end if
    end do
  end do
! Print output
  if (ie == nescf) then
!   pass the magnetization to the kernel routine
!    mtotsusc = msgf
!    mxcsusc  = msgf - mssocz - mssocpm
    mtotsusc = msxc + mssocz + mssocpm
    mxcsusc  = msxc
    write(*,'(" Collinear sum rule: msgflm, msxclm+mssoczlm+mssocpmlm, msxclm, mssoczlm, mssocpmlm")')
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do jlm1=1,lmmax0!lmmax2
        msgflm = 0.d0; msxclm = 0.d0; mssoczlm = 0.d0; mssocpmlm = 0.d0
        do ilm2=1,lmmax
        do ilm1=1,lmmax
          do ib2=1,iwsusc(i2lm(2,ilm2),1,ia)
          do ib1=1,iwsusc(i2lm(2,ilm1),1,ia)
            i2 = lmsb2i(ib2,ilm2,1,ia)
            i1 = lmsb2i(ib1,ilm1,1,ia)
            msgflm    = msgflm    + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*msgf(ib1,ib2,ilm1,ilm2,ia2)
            msxclm    = msxclm    + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*msxc(ib1,ib2,ilm1,ilm2,ia2)
            mssoczlm  = mssoczlm  + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*mssocz(ib1,ib2,ilm1,ilm2,ia2)
            mssocpmlm = mssocpmlm + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*mssocpm(ib1,ib2,ilm1,ilm2,ia2)
          end do
          end do
        end do
        end do
        re = real(msgflm); im = aimag(msgflm)
        if (abs(re) < tol) re = 0.d0
        if (abs(im) < tol) im = 0.d0
        msgflm = cmplx(re,im)
        re = real(msxclm); im = aimag(msxclm)
        if (abs(re) < tol) re = 0.d0
        if (abs(im) < tol) im = 0.d0
        msxclm = cmplx(re,im)
        re = real(mssoczlm); im = aimag(mssoczlm)
        if (abs(re) < tol) re = 0.d0
        if (abs(im) < tol) im = 0.d0
        mssoczlm = cmplx(re,im)
        re = real(mssocpmlm); im = aimag(mssocpmlm)
        if (abs(re) < tol) re = 0.d0
        if (abs(im) < tol) im = 0.d0
        mssocpmlm = cmplx(re,im)
        write(*,'("ia=",i4," lm=",2i4,2es12.4," | ",3es12.4)') ia, i2lm(:,jlm1), real(msgflm), real(msxclm+mssoczlm+mssocpmlm), real(msxclm), real(mssoczlm), real(mssocpmlm)
      end do
    end do
  end if
! All done!
  end subroutine ms_from_bxc3
