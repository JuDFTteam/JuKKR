  subroutine build_vsigmab(ia,ie,vsigmab)
! reads in the self-energy data and prepares the potential

  implicit none

  integer(kind=i4b), intent(in)  :: ia, ie
  complex(kind=c8b), intent(out) :: vsigmab(nlmsb,nlmsb)
! -----------------------------------------------------------------
  integer(kind=i4b) :: is, js, ib, jb, nbi, nbj, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j, nr, je, k
  real(kind=r8b)    :: dr(nrmax), maxnorm, upot, ere, x(100)
  complex(kind=c8b) :: work(nrmax), norm, sigma(nlms,nlms)
  integer(kind=i4b) :: imax, jmax

  vsigmab = 0.d0; maxnorm = 0.d0; sigma = 0.d0
! Only for first atom
  if (ia > 1) return
! Read data from file
  open(file='Vbias_sigma.dat',unit=iofile,status='old')
! Skip header
  read(iofile,*)
  read(iofile,*)
! Loop over energy
  do je=1,nesusc
!   Read the d-block: 5*5 m-values * 2 spins * 2 reals
    read(iofile,*) ere, x
!   Save desired energy point
    if (je == ie) then
!      if (abs(ere - real(esusc(ie))) > 1.d-6) stop 'check energy'
      k = 0
!      do is=1,nsmax
      do is=nsmax,1,-1
!      is = 2
        do ilm=5,9
          do jlm=5,9
            k = k + 1
            i = lms2i(ilm,is); j = lms2i(jlm,is)
            sigma(i,j) = cmplx(x(2*k-1),x(2*k))
!            sigma(i,j) = cmplx(0.d0,x(2*k))
          end do
        end do
      end do
    end if
  end do
  close(iofile)
! ----------------------------------------------------------------------
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
!   ----------------------------------------------------------------
!   selection rules
      if (jl == 2 .and. il == 2) then
!     revise the basis
        norm = sigma(lms2i(ilm,is),lms2i(jlm,js))
        if (abs(norm) > maxnorm) then
          maxnorm = abs(norm)
          imax = i; jmax = j
        end if
        vsigmab(i,j) = vsigmab(i,j) + overlap(i,j,ia)*norm
      end if
!   ----------------------------------------------------------------
    end do
  end do
!  write(iodb,'(" LDA+U kernel, ia=",i4)') ia
!  write(iodb,'("vldaub norm:",6i4,2es16.8)') i2lmsb(:,imax,ia), i2lmsb(:,jmax,ia), maxnorm
!  write(iodb,'("  ib ilm  is, jb jlm  js, matrix element")')
  do j=1,nlmsba(ia)
    do i=1,nlmsba(ia)
      if (abs(vsigmab(i,j)) > 1.d-8) then
!        write(iodb,'(6i4,2es10.2)') i2lmsb(:,i,ia), i2lmsb(:,j,ia), vldaub(i,j)
      end if
    end do
  end do
! All done
  end subroutine build_vsigmab

