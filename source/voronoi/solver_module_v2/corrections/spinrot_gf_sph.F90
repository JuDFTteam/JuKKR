  subroutine spinrot_gf_sph(magdir0,magdir1)
! Rotates the entire projected GF according to local spin rotations
! Assumes that the wfns are spherical
  use global

  implicit none

! Initial spin axes
  real(kind=r8b), intent(in) :: magdir0(3,nasusc)
! Final spin axis
  real(kind=r8b), intent(in) :: magdir1(3,nasusc)
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0)
  complex(kind=c8b) :: spinrot(nsmax,nsmax,nasusc)
  complex(kind=c8b), allocatable :: pzrsave(:,:), pzlsave(:,:), gfsave(:,:)
  integer(kind=i4b) :: ie, i3(3), i2(2), ilms, jlms, ia, ja
  integer(kind=i4b) :: i, j, is, js, ib, jb, ilm, jlm
  integer(kind=i4b) :: p, q, ps, qs, pb, qb, plm, qlm
  real(kind=r8b)    :: start, finish


  allocate(pzrsave(nlmsb,nlms),pzlsave(nlmsb,nlms),gfsave(max(nlms,nlmsb),max(nlms,nlmsb)))
! Spin rotation matrices
  do ia=1,nasusc
    call spin_rotation(magdir0(:,ia),magdir1(:,ia),pauli,spinrot(:,:,ia))
!    write(*,'(i4,8f8.4)') ia, spinrot(:,:,ia)
  end do
  call cpu_time(start)
! **************
  do ie=1,nesusc
  do ia=1,nasusc
! **************
!   Rotation of the spherical scattering solutions
!    call zcopy(nlmsb*nlms,pzr(1,1,ia,ie),1,pzrsave,1)
!    call zcopy(nlmsb*nlms,pzl(1,1,ia,ie),1,pzlsave,1)
!    call zscal(nlmsb*nlms,czero,pzr(1,1,ia,ie),1)
!    call zscal(nlmsb*nlms,czero,pzl(1,1,ia,ie),1)
    pzrsave(1:nlmsba(ia),:) = pzr(1:nlmsba(ia),:,ia,ie)
    pzlsave(1:nlmsba(ia),:) = pzl(1:nlmsba(ia),:,ia,ie)
    pzr(:,:,ia,ie) = czero
    pzl(:,:,ia,ie) = czero
    do ilms=1,nlms
      i2 = i2lms(:,ilms)
      ilm = i2(1); is = i2(2)
      do p=1,nlmsba(ia)
        i3 = i2lmsb(:,p,ia)
        pb = i3(1); plm = i3(2); ps = i3(3)
        do jlms=1,nlms
          i2 = i2lms(:,jlms)
          jlm = i2(1); js = i2(2)
          do q=1,nlmsba(ia)
            i3 = i2lmsb(:,q,ia)
            qb = i3(1); qlm = i3(2); qs = i3(3)
!           selection rules
            if (ilm == jlm .and. plm == qlm .and. plm == ilm .and. qlm == jlm .and. pb == qb) then
              pzr(p,ilms,ia,ie) = pzr(p,ilms,ia,ie) + spinrot(ps,qs,ia)*pzrsave(q,jlms)*conjg(spinrot(is,js,ia))
              pzl(p,ilms,ia,ie) = pzl(p,ilms,ia,ie) + conjg(spinrot(ps,qs,ia))*pzlsave(q,jlms)*spinrot(is,js,ia)
            end if
          end do
        end do
      end do
    end do
!   Rotation of the spherical single-site GF
    gfsave(1:nlmsba(ia),1:nlmsba(ia)) = gfpq(1:nlmsba(ia),1:nlmsba(ia),ia,ie)
    gfpq(:,:,ia,ie) = 0.d0
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        do q=1,nlmsba(ia)
          i3 = i2lmsb(:,q,ia)
          qb = i3(1); qlm = i3(2); qs = i3(3)
          do p=1,nlmsba(ia)
            i3 = i2lmsb(:,p,ia)
            pb = i3(1); plm = i3(2); ps = i3(3)
!           selection rules
            if (ilm == jlm .and. plm == qlm .and. plm == ilm .and. qlm == jlm .and. ib == pb .and. jb == qb) then
              gfpq(i,j,ia,ie) = gfpq(i,j,ia,ie) + spinrot(is,ps,ia)*gfsave(p,q)*conjg(spinrot(js,qs,ia))
            end if
          end do
        end do
      end do
    end do
! ******
  end do
  end do
! ******
  call cpu_time(finish)
  write(*,'("spinrot_gf: site diagonal time",f8.3," s")') finish - start
! ----------------------------------------------------------------------
  call cpu_time(start)
! **************
  do ie=1,nesusc
  do ja=1,nasusc
  do ia=1,nasusc
! **************
!   Rotation of the structural GF
!   save current block
    do jlms=1,nlms
      j = alms2i(jlms,ja)
      do ilms=1,nlms
        i = alms2i(ilms,ia)
        gfsave(ilms,jlms) = gstruct(i,j,ie)
        gstruct(i,j,ie) = czero
      end do
    end do
!   now rotate
    do jlms=1,nlms
      j = alms2i(jlms,ja)
      i2 = i2lms(:,jlms)
      jlm = i2(1); js = i2(2)
      do ilms=1,nlms
        i = alms2i(ilms,ia)
        i2 = i2lms(:,ilms)
        ilm = i2(1); is = i2(2)
        do q=1,nlms
          i2 = i2lms(:,q)
          qlm = i2(1); qs = i2(2)
          do p=1,nlms
            i2 = i2lms(:,p)
            plm = i2(1); ps = i2(2)
!           selection rules
            if (ilm == plm .and. jlm == qlm) then
              gstruct(i,j,ie) = gstruct(i,j,ie) + spinrot(is,ps,ia)*gfsave(p,q)*conjg(spinrot(js,qs,ja))
            end if
          end do
        end do
      end do
    end do
! ******
  end do
  end do
  end do
! ******
  call cpu_time(finish)
  write(*,'("spinrot_gf: structural gf time",f8.3," s")') finish - start
  deallocate(pzrsave,pzlsave,gfsave)
! All done!
  end subroutine spinrot_gf_sph
