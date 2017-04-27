  subroutine ms_from_bxc4(ie,bxc,bsoc,magdir1,msgf,msxc,mssoc)
! noncollinear magnetization sum rule
! the ordering of the labels on the product of GFs was changed to agree with static_susc
! Bext is included in Bsoc (check pot_correction)
  use global

  implicit none

! Current energy
  integer(kind=i4b), intent(in)    :: ie
! twist Bxc
  complex(kind=c8b), intent(in)    :: bxc(nlmsb,nlmsb,3,nasusc)
! twist SOC + magnetic field + ...
  complex(kind=c8b), intent(in)    :: bsoc(nlmsb,nlmsb,3,nasusc)
! spin quantization axes
  real(kind=r8b),    intent(in)    :: magdir1(3,nasusc)
! magnetization constructed from site-diagonal block of GF
  complex(kind=c8b), intent(inout) :: msgf(3,nbmax,nbmax,lmmax,lmmax,nasusc2)
! contribution from the xc potential
  complex(kind=c8b), intent(inout) :: msxc(3,nbmax,nbmax,lmmax,lmmax,nasusc,nasusc2)
! contribution from remaining spin-dependent potentials
  complex(kind=c8b), intent(inout) :: mssoc(3,nbmax,nbmax,lmmax,lmmax,nasusc,nasusc2)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0), tol = 1.d-7, uz(3) = (/0.d0,0.d0,1.d0/)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi), iu = (0.d0,1.d0), czero = (0.d0,0.d0), cone = (1.d0,0.d0)
  integer(kind=i4b) :: i3(3), ia, ia2, i1, ib1, ilm1, is1, i2, ib2, ilm2, is2, ja, j1, jb1, jlm1, js1, j2, jb2, jlm2, js2, i, j, k
  complex(kind=c8b), allocatable :: gfij(:,:), gfji(:,:), work(:,:), work1(:,:), work2(:,:)
  real(kind=r8b),    allocatable :: msgflm(:,:), msxclm(:,:), mssoclm(:,:), rotmat(:,:,:)
  complex(kind=c8b) :: tmp1, tmp2, mvec(3)
  real(kind=r8b)    :: rotgen(3,3,3), rvec(3)
  complex(kind=c8b) :: twistpauli(2,2,3,3)
  real(kind=r8b)    :: re, im
  integer(kind=i4b) :: jlms, jb, jlm, js, ilms, ib, ilm, is, iatom

!  write(*,'("ms_from_bxc4: ie=",i4)') ie
! ----------------------------------------------------------------------
! Generators of spatial rotations
! Any similarity with the Levi-Civita symbol is not coincidence
  rotgen(:,:,:) = 0.d0
! x
  rotgen(2,3,1) =  1.d0; rotgen(3,2,1) = -1.d0
! y
  rotgen(1,3,2) = -1.d0; rotgen(3,1,2) =  1.d0
! z
  rotgen(1,2,3) =  1.d0; rotgen(2,1,3) = -1.d0
! ----------------------------------------------------------------------
! Twisted Pauli matrices: for each generator rotate Pauli matrices
  twistpauli(:,:,:,:) = 0.d0
  do j=1,3
    do i=1,3
      do k=1,3
        twistpauli(:,:,i,j) = twistpauli(:,:,i,j) + rotgen(i,k,j)*pauli(:,:,k)
      end do
    end do
  end do
! ----------------------------------------------------------------------
! Loop over all atoms and collect all contributions
  allocate(gfij(nlmsb,nlmsb),gfji(nlmsb,nlmsb),work(nlmsb,nlmsb),work1(nlmsb,nlmsb),work2(nlmsb,nlmsb))
! ***********************************
! Loop over atoms for susceptibility
  suscatoms: do ia2=1,nasusc2
    ia = iasusc2(ia2)
!   Gij(E+i0)
    call projected_gf(ie,ia,ia,gfij,.true.,.true.)
!    if (lrot) call local_frame(ia,ia,magdir,gfij,work)
!   ------------------------------------------------------------------
!   magnetization from site-diagonal GF
    do i2=1,nlmsba(ia)
      i3 = i2lmsb(:,i2,ia)
      ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
      do i1=1,nlmsba(ia)
        i3 = i2lmsb(:,i1,ia)
        ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
!       vector components of spin density
        do i=1,3
          tmp1 = -descf(ie)*pauli(is2,is1,i)*gfij(i1,i2)
          tmp1 = tmp1 + conjg(descf(ie)*pauli(is1,is2,i)*gfij(i2,i1)) 
          msgf(i,ib1,ib2,ilm1,ilm2,ia2) = msgf(i,ib1,ib2,ilm1,ilm2,ia2) + tmp1/i2pi
        end do
      end do
    end do
!   ------------------------------------------------------------------
!   Loop over all atoms in the system
    allatoms: do ja=1,nasusc
! ***********************************
!     Gij(E+i0)
      call projected_gf(ie,ia,ja,gfij,.true.,.true.)
!      if (lrot) call local_frame(ia,ja,magdir,gfij,work)
!     Gji(E+i0)
      call projected_gf(ie,ja,ia,gfji,.true.,.true.)
!      if (lrot) call local_frame(ja,ia,magdir,gfji,work)
!     magnetization from all spin-dependent potentials in the system
!     --------------------------------------------------------------
      gen: do k=1,3
!       ----------------------------------------------------------------
!       Gij(E+i0) (bxc x s) --> work 
        call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,gfij,nlmsb,bxc(:,:,k,ja),nlmsb,czero,work,nlmsb)
!       -dE * ( Gij(E+i0) (bxc x s) ) Gji(E+i0) --> work1
        call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ja),-descf(ie),work,nlmsb,gfji,nlmsb,czero,work1,nlmsb)
!       Gij(E-i0) (bxc x s) --> work 
        call zgemm('C','C',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,gfji,nlmsb,bxc(:,:,k,ja),nlmsb,czero,work,nlmsb)
!       dE^* * ( Gij(E-i0) (bxc x s) ) Gji(E-i0) --> work1
        call zgemm('N','C',nlmsba(ia),nlmsba(ia),nlmsba(ja),conjg(descf(ie)),work,nlmsb,gfij,nlmsb,cone,work1,nlmsb)
!       ----------------------------------------------------------------
!       Gij(E+i0) (bsoc x s) --> work 
        call zgemm('N','N',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,gfij,nlmsb,bsoc(:,:,k,ja),nlmsb,czero,work,nlmsb)
!       -dE * ( Gij(E+i0) (bsoc x s) ) Gji(E+i0) --> work2
        call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ja),-descf(ie),work,nlmsb,gfji,nlmsb,czero,work2,nlmsb)
!       Gij(E-i0) (bsoc x s) --> work 
        call zgemm('C','C',nlmsba(ia),nlmsba(ja),nlmsba(ja),cone,gfji,nlmsb,bsoc(:,:,k,ja),nlmsb,czero,work,nlmsb)
!       dE^* * ( Gij(E-i0) (bsoc x s) ) Gji(E-i0) --> work2
        call zgemm('N','C',nlmsba(ia),nlmsba(ia),nlmsba(ja),conjg(descf(ie)),work,nlmsb,gfij,nlmsb,cone,work2,nlmsb)
!       ----------------------------------------------------------------
!       cartesian components of the result
        ixyz: do i=1,3
!         basis labels for susc atom
          do i2=1,nlmsba(ia)
          do i1=1,nlmsba(ia)
            i3 = i2lmsb(:,i1,ia)
            ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
            i3 = i2lmsb(:,i2,ia)
            ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
            tmp1 = twistpauli(is1,is2,i,k)*work1(i2,i1)
            msxc(i,ib1,ib2,ilm1,ilm2,ja,ia2)  = msxc(i,ib1,ib2,ilm1,ilm2,ja,ia2)  + 0.5d0*tmp1/i2pi
            tmp2 = twistpauli(is1,is2,i,k)*work2(i2,i1)
            mssoc(i,ib1,ib2,ilm1,ilm2,ja,ia2) = mssoc(i,ib1,ib2,ilm1,ilm2,ja,ia2) + 0.5d0*tmp2/i2pi
          end do
          end do
        end do ixyz
      end do gen
!     ----------------------------------------------------------------
! *****************
    end do allatoms
  end do suscatoms
! *****************
! ----------------------------------------------------------------------
! Print output
  if (ie == nescf) then
    write(*,'(" Noncollinear sum rule: msgflm, msxclm+mssoclm, msxclm, mssoclm")')
    allocate(rotmat(3,3,nasusc2),msgflm(3,lmmax0),msxclm(3,lmmax0),mssoclm(3,lmmax0))
    mtotsusc = 0.d0; mxcsusc = 0.d0; msocsusc = 0.d0
    suscatom: do ia2=1,nasusc2
      ia = iasusc2(ia2)
!     rotation matrices to local frame
      call rotvec(uz,magdir(:,ia),rotmat(:,:,ia2))
!      call rotvec(uz,uz,rotmat(:,:,ia2))
      msgflm = 0.d0; msxclm = 0.d0; mssoclm = 0.d0
      lm2: do ilm2=1,lmmax
      lm1: do ilm1=1,lmmax
!       spin-independent basis assumed
        b2: do ib2=1,iwsusc(i2lm(2,ilm2),1,ia)
        b1: do ib1=1,iwsusc(i2lm(2,ilm1),1,ia)
          i2 = lmsb2i(ib2,ilm2,1,ia)
          i1 = lmsb2i(ib1,ilm1,1,ia)
!         pass the z-magnetization in the local frame to the kernel routine:  m_local = R^T m_global
!         ------------------------------------------------------------
!         spin density from onsite GF in local frame
          mvec(:) = msgf(:,ib1,ib2,ilm1,ilm2,ia2)
          msgf(:,ib1,ib2,ilm1,ilm2,ia2) = matmul(transpose(rotmat(:,:,ia2)),mvec(:))
!         ------------------------------------------------------------
!         spin density from all SOC+Bext atoms
          mvec(:) = sum(mssoc(:,ib1,ib2,ilm1,ilm2,:,ia2),dim=2)
!         spin density from SOC+Bext atoms in susc
!          mvec(:) = sum(mssoc(:,ib1,ib2,ilm1,ilm2,iasusc2(1:nasusc2),ia2),dim=2)
          msocsusc(ib1,ib2,ilm1,ilm2,ia2) = sum(rotmat(:,3,ia2)*mvec(:))
!         ------------------------------------------------------------
!         spin density from all Bxc atoms
!          mvec(:) = sum(msxc(:,ib1,ib2,ilm1,ilm2,:,ia2),dim=2)
!         spin density from Bxc atoms in susc
          mvec(:) = sum(msxc(:,ib1,ib2,ilm1,ilm2,iasusc2(1:nasusc2),ia2),dim=2)
          mxcsusc(ib1,ib2,ilm1,ilm2,ia2) = sum(rotmat(:,3,ia2)*mvec(:))
!         ------------------------------------------------------------
!         total spin density from onsite GF
!          mtotsusc(ib1,ib2,ilm1,ilm2,ia2) = msgf(ib1,ib2,ilm1,ilm2,ia2)
!         total spin density from Bxc + SOC + Bext
          mtotsusc(ib1,ib2,ilm1,ilm2,ia2) = mxcsusc(ib1,ib2,ilm1,ilm2,ia2) + msocsusc(ib1,ib2,ilm1,ilm2,ia2)
!         ------------------------------------------------------------
!         expand the spin densities in multipoles
          lmden: do jlm1=1,lmmax0!lmmax2
!           spin density from onsite GF in local frame
            msgflm(:,jlm1)  = msgflm(:,jlm1)  + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*msgf(:,ib1,ib2,ilm1,ilm2,ia2)
!           spin density created by Bxc of all atoms
!            msxclm(:,jlm1)  = msxclm(:,jlm1)  + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*sum(msxc(:,ib1,ib2,ilm1,ilm2,:,ia2),dim=2)
!           spin density created by Bxc of susc atoms
            msxclm(:,jlm1)  = msxclm(:,jlm1)  + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*sum(msxc(:,ib1,ib2,ilm1,ilm2,iasusc2(1:nasusc2),ia2),dim=2)
!           spin density created by SOC and Bext of all atoms
            mssoclm(:,jlm1) = mssoclm(:,jlm1) + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*sum(mssoc(:,ib1,ib2,ilm1,ilm2,:,ia2),dim=2)
!           spin density created by SOC and Bext of susc atoms
!            mssoclm(:,jlm1) = mssoclm(:,jlm1) + rgaunt(ilm1,ilm2,jlm1)*overlap(i1,i2,ia)*sum(mssoc(:,ib1,ib2,ilm1,ilm2,iasusc2(1:nasusc2),ia2),dim=2)
          end do lmden
!         ------------------------------------------------------------
        end do b1
        end do b2
      end do lm1
      end do lm2
!     kill small numbers
      where (abs(msgflm) < tol)  msgflm  = 0.d0
      where (abs(msxclm) < tol)  msxclm  = 0.d0
      where (abs(mssoclm) < tol) mssoclm = 0.d0
!     print contributions to current atom
      do jlm1=1,lmmax0
        msxclm(:,jlm1)  = matmul(transpose(rotmat(:,:,ia2)),msxclm(:,jlm1))
        mssoclm(:,jlm1) = matmul(transpose(rotmat(:,:,ia2)),mssoclm(:,jlm1))
        write(*,'("ia=",i4," lm=",2i4,2(3es12.4,2x)" | ",3(3es12.4,2x))') ia, i2lm(:,jlm1), msgflm(:,jlm1), &
            msxclm(:,jlm1)+mssoclm(:,jlm1), msxclm(:,jlm1), mssoclm(:,jlm1)
      end do
    end do suscatom
    deallocate(rotmat,msgflm,msxclm,mssoclm)
  end if
  deallocate(gfij,gfji,work,work1,work2)
! All done!
  end subroutine ms_from_bxc4
