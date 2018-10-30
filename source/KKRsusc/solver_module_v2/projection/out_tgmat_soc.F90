  subroutine out_tgmat_soc(nhost,natyp,lmmaxp,nsec,tmat,gmat,ie,is,write_mpi)
! Output t-mat and structural GF to outsusc.dat
! BIG FUCK-UP: the correct t-matrix from kkrimp is tmatll, not dtmtrx!
! New: keep it in memory
  use global, only: i4b, r8b, c8b, lhdio, iomain, iodb, gfilter, nasusc, iasusc, lmmax, nesusc, magdir

  implicit none

! --> number of atoms for host and total
  integer(kind=i4b), intent(in) :: nhost, natyp
! --> size of the angular momentum matrix
  integer(kind=i4b), intent(in) :: lmmaxp
! --> size of the GF matrix
  integer(kind=i4b), intent(in) :: nsec
! --> the t-matrices
  complex(kind=c8b), intent(in) :: tmat(lmmaxp,lmmaxp,natyp)
! --> the structural GF matrix
  complex(kind=c8b), intent(in) :: gmat(nsec,nsec)
! --> current energy and spin channel
  integer(kind=i4b), intent(in) :: ie, is
! --> write ? 
  logical,           intent(in) :: write_mpi
! -----------------------------------------------------------------
  complex(kind=c8b) :: gtrim(2*lmmax,2*lmmax)
  real(kind=r8b)    :: gmax, re, im
  integer(kind=i4b) :: ia, ja, i1, i2, ilm, jlm
  real(kind=r8b)    :: start, finish, magdir0(3,nasusc), magdir1(3,nasusc)

  do ia=1,nasusc
    i1 = (iasusc(ia)-nhost-1)*lmmaxp
!    write(iodb,'("out_tmat: ie, is, dims=",8i8)') ie, is, nhost, lmmaxp, nsec, ia, i1, lmmax
!    gmax = maxval(abs(tmat(i1+1:i1+lmmax,i1+1:i1+lmmax)))
!    do jlm=1,lmmax
!      do ilm=1,lmmax
!        re = real(tmat(i1+ilm,i1+jlm))
!        im = aimag(tmat(i1+ilm,i1+jlm))
!        if (abs(re) < gfilter*gmax) re = 0.d0
!        if (abs(im) < gfilter*gmax) im = 0.d0
!        gtrim(ilm,jlm) = cmplx(re,im)
!      end do
!    end do
!    gtrim(1:lmmax,1:lmmax) = tmat(i1+1:i1+lmmax,i1+1:i1+lmmax)
    gtrim(1:2*lmmax,1:2*lmmax) = tmat(1:2*lmmax,1:2*lmmax,iasusc(ia))
    if (lhdio .and. write_mpi) then
!     To disk ram for MPI      
      do jlm=1,2*lmmax
        write(iomain,'(1000es16.8)') gtrim(1:2*lmmax,jlm)
      end do
    else
!     to ram
      call save_tmati(gtrim,is,ie,ia,2*lmmax)
    end if
  end do
! ----------------------------------------------------------------------
  do ja=1,nasusc
    i2 = (iasusc(ja)-nhost-1)*lmmaxp
    do ia=1,nasusc
      i1 = (iasusc(ia)-nhost-1)*lmmaxp
!      write(iodb,'("out_gstruct: ie, is, dims=",10i8)') ie, is, nhost, lmmaxp, nsec, ia, ja, i1, i2, lmmax
!      gmax = maxval(abs(gmat(i1+1:i1+lmmax,i2+1:i2+lmmax)))
!      do jlm=1,lmmax
!        do ilm=1,lmmax
!          re = real(gmat(i1+ilm,i2+jlm))
!          im = aimag(gmat(i1+ilm,i2+jlm))
!          if (abs(re) < gfilter*gmax) re = 0.d0
!          if (abs(im) < gfilter*gmax) im = 0.d0
!          gtrim(ilm,jlm) = cmplx(re,im)
!        end do
!      end do
      gtrim(1:2*lmmax,1:2*lmmax) = gmat(i1+1:i1+2*lmmax,i2+1:i2+2*lmmax)
!      write(*,'(1000es16.8)') gtrim(1:2*lmmax,jlm)
      if (lhdio .and. write_mpi) then
!       To disk and ram for MPI
        do jlm=1,2*lmmax
          write(iomain,'(1000es16.8)') gtrim(1:2*lmmax,jlm)
        end do
      else
!       to ram
        call save_gsij(gtrim,is,ie,ia,ja,2*lmmax)
      end if
    end do
  end do
! All done
  end subroutine out_tgmat_soc
