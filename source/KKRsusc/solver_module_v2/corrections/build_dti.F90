  subroutine build_dti(ie,ia,lldautoo,dti,jii)
! Matrix elements of Bxc with regular scattering solutions
! Eq. 7 of PRB 79, 045209 (2009)
! Augmented with LDA+U contribution

  use global

  implicit none

! Which energy
  integer(kind=i4b), intent(in) :: ie
! Which atom
  integer(kind=i4b), intent(in) :: ia
! Whether to add LDA+U correction to spin splitting
  logical,  intent(in) :: lldautoo
! t-matrix perturbations, Pauli matrices in basis
!  complex(kind=c8b), intent(out) :: dti(nlms,nlms,3)
  complex(kind=c8b), intent(out) :: dti(nlmsb,nlmsb,6)
! magnetic exchange energy
  real(kind=r8b),    intent(out) :: jii
! ----------------------------------------------------------------------
  integer(kind=i4b) :: nr, i2(2), i3(3), i, j, ilms, jlms, ilmsb, jlmsb, is, js, il, im, jl, jm, ilm, jlm, ib, jb
  real(kind=r8b)    :: uminusj, dr(nrmax), bxcl(nbmax,nbmax,0:nlmax), di1(3,3)
  complex(kind=c8b) :: work(nrmax), norm, vldauspin(nbmax,nbmax,lmmax,lmmax), sproj(nsmax,nsmax)
  complex(kind=c8b), external :: radint

! projection on local spin quantization axis
  sproj = 0.d0
  do i=1,3
!    sproj = sproj + 0.5d0*pauli(:,:,i)*magdir(i,ia)
    sproj = sproj + pauli(:,:,i)*magdir(i,ia)
  end do
! derivates of direction of magnetic moment
  do j=1,3
    do i=1,3
      if (i == j) cycle
      di1(i,j) = -magdir(i,ia)*magdir(j,ia)
    end do
    di1(j,j) = (1.d0 - magdir(j,ia))*(1.d0 + magdir(j,ia))
  end do
! Assumption: ASA
! Assumption: basis independent of spin -> making the following code silly
  nr = nrpts(ia)
  dr(1:nr) = drmesh(1:nr,ia)
  bxcl = 0.d0
  do is=1,nsmax
    do il=0,nlmax
      do jb=1,iwsusc(il,is,ia)
        do ib=1,iwsusc(il,is,ia)
          work(1:nr) = phiref(1:nr,ib,il,is,ia)*br(1:nr,ia)*phiref(1:nr,jb,il,is,ia)
!         spin-independent basis: bxcl is doubled, so this is divided here
          bxcl(ib,jb,il) = bxcl(ib,jb,il) + real(radint(nr,work,dr,npanat(ia),ircutat(:,ia)))/nsmax
        end do
      end do
!      if (ie == 1) write(iodb,'(2i4,100es16.8)') ia, il, bxcl(:,:,il)
    end do
  end do
! Magnetic exchange energy
  work(1:nr) = br(1:nr,ia)*new_rho2ns(1:nr,1,2,ia)
  jii = real(radint(nr,work,dr,npanat(ia),ircutat(:,ia)))
! ----------------------------------------------------------------------
! Assumption: Dudarev LDA+U
  vldauspin = 0.d0
  if (lldautoo .and. lldau .and. ildau(ia) /= 0) then
!   contribution from LDA+U to spin splitting
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
!     ----------------------------------------------------------------
!     selection rules
        if (jl == il) then
          uminusj = ueff(il,ia) - jeff(il,ia)
!       revise the basis
          norm = -rhomat(lms2i(ilm,is),lms2i(jlm,js),ia)
          if (ilm == jlm .and. is == js) norm = norm + 0.5d0 + vshift(il,is,ia)
          vldauspin(ib,jb,ilm,jlm) = vldauspin(ib,jb,ilm,jlm) + uminusj*overlap(i,j,ia)*norm*sproj(js,is)
        end if
!     ----------------------------------------------------------------
      end do
    end do
  end if
! ----------------------------------------------------------------------
  dti = 0.d0
  do i=1,3
!    if (i == 1) sproj =  pauli(:,:,1)*magdir(3,ia) - pauli(:,:,3)*magdir(1,ia)
!    if (i == 2) sproj =  pauli(:,:,2)*magdir(3,ia) - pauli(:,:,3)*magdir(2,ia)
!    if (i == 3) sproj =  pauli(:,:,2)*magdir(1,ia) - pauli(:,:,1)*magdir(2,ia)
    do jlmsb=1,nlmsba(ia)
      i3 = i2lmsb(:,jlmsb,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      do ilmsb=1,nlmsba(ia)
        i3 = i2lmsb(:,ilmsb,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
!       assumption: ASA Bxc
        if (ilm == jlm) then
          il = i2lm(2,ilm)
          dti(ilmsb,jlmsb,i) = dti(ilmsb,jlmsb,i) + pauli(is,js,i)*bxcl(ib,jb,il)
!          dti(ilmsb,jlmsb,i) = dti(ilmsb,jlmsb,i) + sum(pauli(is,js,1:3)*di1(1:3,i))*bxcl(ib,jb,il)
!         Pauli matrices in the projection basis
          if (ib == jb) dti(ilmsb,jlmsb,i+3) = pauli(is,js,i)
!          dti(ilmsb,jlmsb,i) = dti(ilmsb,jlmsb,i) + sproj(is,js)*bxcl(ib,jb,il)
        end if
!       LDA+U
        dti(ilmsb,jlmsb,i) = dti(ilmsb,jlmsb,i) + pauli(is,js,i)*vldauspin(ib,jb,ilm,jlm)
!        dti(ilmsb,jlmsb,i) = dti(ilmsb,jlmsb,i) + sum(pauli(is,js,1:3)*di1(1:3,i))*vldauspin(ib,jb,ilm,jlm)
!        dti(ilmsb,jlmsb,i) = dti(ilmsb,jlmsb,i) + sproj(is,js)*vldauspin(ib,jb,ilm,jlm)
      end do
    end do
  end do
!  do i=1,3
!    do jlms=1,nlms
!    do ilms=1,nlms
!      do jlmsb=1,nlmsba(ia)
!        i3 = i2lmsb(:,jlmsb,ia)
!        jb = i3(1); jlm = i3(2); js = i3(3)
!        do ilmsb=1,nlmsba(ia)
!          i3 = i2lmsb(:,ilmsb,ia)
!          ib = i3(1); ilm = i3(2); is = i3(3)
!         assumption: ASA Bxc
!          if (ilm == jlm) then
!            il = i2lm(2,ilm)
!            dti(ilms,jlms,i) = dti(ilms,jlms,i) + pzl(ilmsb,ilms,ia,ie)*pauli(is,js,i)*bxcl(ib,jb,il)*pzr(jlmsb,jlms,ia,ie)
!          end if
!         LDA+U
!          dti(ilms,jlms,i) = dti(ilms,jlms,i) + pzl(ilmsb,ilms,ia,ie)*pauli(is,js,i)*vldauspin(ib,jb,ilm,jlm)*pzr(jlmsb,jlms,ia,ie)
!        end do
!      end do
!    end do
!    end do
!  end do
! All done!
  end subroutine build_dti
