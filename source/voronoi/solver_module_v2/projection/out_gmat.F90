  subroutine out_gmat(nhost,lmmaxp,nsec,gmat,ie,is)
! Output structural GF to outsusc.dat
! New: keep it in memory
  use global

  implicit none

! --> number of atoms for host
  integer(kind=i4b), intent(in) :: nhost
! --> size of the angular momentum matrix
  integer(kind=i4b), intent(in) :: lmmaxp
! --> size of the GF matrix
  integer(kind=i4b), intent(in) :: nsec
! --> the structural GF matrix
  complex(kind=c8b), intent(in) :: gmat(nsec,nsec)
! --> current energy and spin channel
  integer(kind=i4b), intent(in) :: ie, is
! -----------------------------------------------------------------
  complex(kind=c8b) :: gtrim(lmmax,lmmax)
  real(kind=r8b)    :: gmax, re, im
  integer(kind=i4b) :: ia, ja, i1, i2, ilm, jlm

  if (lhdio) write(iomain,'(" GF: ia,ja, iasusc(ia),iasusc(ja), i1,i2")')
  do ja=1,nasusc
    i2 = (iasusc(ja)-nhost-1)*lmmaxp
    do ia=1,nasusc
      i1 = (iasusc(ia)-nhost-1)*lmmaxp
      if (lhdio) write(iomain,'(6i8)') ia,ja, iasusc(ia),iasusc(ja), i1,i2
      gmax = maxval(abs(gmat(i1+1:i1+lmmax,i2+1:i2+lmmax)))
      do jlm=1,lmmax
        do ilm=1,lmmax
          re = real(gmat(i1+ilm,i2+jlm))
          im = aimag(gmat(i1+ilm,i2+jlm))
          if (abs(re) < gfilter*gmax) re = 0.d0
          if (abs(im) < gfilter*gmax) im = 0.d0
          gtrim(ilm,jlm) = cmplx(re,im)
        end do
      end do
      if (lhdio) then
!       to disk
        do jlm=1,lmmax
          write(iomain,'(1000es16.8)') gtrim(1:lmmax,jlm)
        end do
      else
!       to ram
        call save_gsij(gtrim,is,ie,ia,ja,lmmaxp)
      end if
    end do
  end do
! All done
  end subroutine out_gmat
