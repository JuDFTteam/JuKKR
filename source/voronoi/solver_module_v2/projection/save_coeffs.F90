  subroutine save_coeffs(ie,is,symmetrize)
! Put projection coefficients in memory
! The onsite term is added to the big GF array
! For the collinear SRA case
  use global

  implicit none

! Current energy point
  integer(kind=i4b), intent(in) :: ie
! Current spin channel
  integer(kind=i4b), intent(in) :: is
! Whether to symmetrize
  logical,           intent(in) :: symmetrize
! -----------------------------------------------------------------
  integer(kind=i4b) :: i, i1, j, il, lm, im, ia, ib, jb, nb, ng, istart, iend
  complex(kind=c8b) :: work

! Symmetrization
  if (symmetrize) then
    iend = 0
    do ng=1,ngroup
      istart = iend + 1
      iend   = iend + igroup(ng)
      do il=0,nlmax
        nb = iwsusc(il,is,istart)
        if (any(iwsusc(il,is,istart:iend) /= nb)) then
          write(*,'("save_coeffs: symmetrization failed for istart,iend=",100i4)') istart, iend, il, is, iwsusc(il,is,istart:iend)
          stop
        end if
!       ****************
        if (nb > 0) then
!       ****************
          do ib=1,nb
            work = sum(pzc(ib,il,is,istart:iend))/igroup(ng)
            pzc(ib,il,is,istart:iend) = work
            if (isra == 1) then
              work = sum(fzc(ib,il,is,istart:iend))/igroup(ng)
              fzc(ib,il,is,istart:iend) = work
            end if
            do jb=1,nb
              work = sum(pqc(ib,jb,il,is,istart:iend))/igroup(ng)
              pqc(ib,jb,il,is,istart:iend) = work
              if (isra == 1) then
                work = sum(fqc(ib,jb,il,is,istart:iend))/igroup(ng)
                fqc(ib,jb,il,is,istart:iend) = work
                work = sum(psc(ib,jb,il,is,istart:iend))/igroup(ng)
                psc(ib,jb,il,is,istart:iend) = work
                work = sum(fsc(ib,jb,il,is,istart:iend))/igroup(ng)
                fsc(ib,jb,il,is,istart:iend) = work
              end if
            end do
          end do
!       ******
        end if
!       ******
      end do
    end do
  end if
! -----------------------------------------------------------------
! Put coefficients in memory
!  pzl(:,:,:,ie) = 0.d0
!  pzr(:,:,:,ie) = 0.d0
!  if (isra == 1) then
!    pzl(:,:,:,ie) = 0.d0
!    pzr(:,:,:,ie) = 0.d0
!  end if
!  gfpq(:,:,:,ie) = 0.d0
!  if (isra == 1) then
!    gffq(:,:,:,ie) = 0.d0
!    gfps(:,:,:,ie) = 0.d0
!    gffs(:,:,:,ie) = 0.d0
!  end if
  do ia=1,nasusc
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
!   Check if the angular momentum block is to keep
      if (nb > 0) then
        do im=-il,il
          lm = lm2i(im,il)
          do ib=1,nb
!       Saving wfn coefficients
            i = lmsb2i(ib,lm,is,ia)
            i1 = lms2i(lm,is)
!            write(*,'(10i4)') ia, is, lm, ib, i
            pzl(i,i1,ia,ie) = pzc(ib,il,is,ia)
            pzr(i,i1,ia,ie) = pzc(ib,il,is,ia)
            if (isra == 1) then
              fzl(i,i1,ia,ie) = fzc(ib,il,is,ia)
              fzr(i,i1,ia,ie) = fzc(ib,il,is,ia)
            end if
!           -------------------------------------------------------
!         Diagonal blocks coming from the RH term in the GF
            do jb=1,nb
              j = lmsb2i(jb,lm,is,ia)
              gfpq(i,j,ia,ie) = pqc(ib,jb,il,is,ia)
!          TEST TEST TEST TEST TEST TEST TEST TEST TEST
!              if (ib /= jb) gfpq(i,j,ia,ie) = 0.d0
!          TEST TEST TEST TEST TEST TEST TEST TEST TEST
              if (isra == 1) then
                gffq(i,j,ia,ie) = fqc(ib,jb,il,is,ia)
                gfps(i,j,ia,ie) = psc(ib,jb,il,is,ia)
                gffs(i,j,ia,ie) = fsc(ib,jb,il,is,ia)
              end if
            end do
!           -------------------------------------------------------
          end do
        end do
      end if
    end do
  end do
! All done
  end subroutine save_coeffs
