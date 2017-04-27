  subroutine build_bsoctwist(vsocb,bsoctwist)
! assembles the SOC potential in the projection basis
  use global !, only: i4b, r8b, c8b, nrmax, nrpts, drmesh, phiref, nlmsba, nlmsb, nlms, i2lmsb, i2lms, i2lm, ldots, atol, pauli, lorb

  implicit none

! SOC potential in the basis
  complex(kind=c8b), intent(in)  :: vsocb(nlmsb,nlmsb,nasusc)
! SOC potential for noncollinear sum rule
  complex(kind=c8b), intent(out) :: bsoctwist(nlmsb,nlmsb,3,nasusc)
! -----------------------------------------------------------------
  real(kind=r8b), parameter :: uz(3) = (/0.d0,0.d0,1.d0/)
  integer(kind=i4b) :: ia, is, js, ib, jb, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j, k, l
  complex(kind=c8b) :: proj(nbmax,nbmax,lmmax,lmmax,3)
  real(kind=r8b)    :: rotmat(3,3,nasusc)

  bsoctwist = 0.d0
! ----------------------------------------------------------------------
  do ia=1,nasusc
!   rotation matrices to local frame
!    call rotvec(uz,magdir(:,ia),rotmat(:,:,ia))
    call rotvec(uz,uz,rotmat(:,:,ia))
!   spin-dependent part of the potential
    proj = 0.d0
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      i2 = i2lm(:,jlm)
      jm = i2(1); jl = i2(2)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
!        proj(ib,jb,ilm,jlm,:) = proj(ib,jb,ilm,jlm,:) + 0.5d0*pauli(js,is,1:3)*vsocb(i,j,ia)
        proj(ib,jb,ilm,jlm,:) = proj(ib,jb,ilm,jlm,:) + 0.5d0*matmul(transpose(rotmat(:,:,ia)),pauli(js,is,1:3))*vsocb(i,j,ia)
      end do
    end do
!   tensor product in cartesian components
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      i2 = i2lm(:,jlm)
      jm = i2(1); jl = i2(2)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        i2 = i2lm(:,ilm)
        im = i2(1); il = i2(2)
!       components of \vec{B} \times \vec{\sigma}
        bsoctwist(i,j,1,ia) = bsoctwist(i,j,1,ia) + proj(ib,jb,ilm,jlm,2)*pauli(is,js,3) - proj(ib,jb,ilm,jlm,3)*pauli(is,js,2)
        bsoctwist(i,j,2,ia) = bsoctwist(i,j,2,ia) + proj(ib,jb,ilm,jlm,3)*pauli(is,js,1) - proj(ib,jb,ilm,jlm,1)*pauli(is,js,3)
        bsoctwist(i,j,3,ia) = bsoctwist(i,j,3,ia) + proj(ib,jb,ilm,jlm,1)*pauli(is,js,2) - proj(ib,jb,ilm,jlm,2)*pauli(is,js,1)
      end do
    end do
  end do
! All done
  end subroutine build_bsoctwist
