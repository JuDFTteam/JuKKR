  subroutine ratfit_gf2(ia,ja,onsite,struct)
! Rational function fit
! TEST VERSION
  use global

  implicit none

! Block of GF
  integer(kind=i4b), intent(in)  :: ia, ja
! To add onsite or structural parts of GF
  logical,           intent(in)  :: onsite, struct
! -----------------------------------------
  integer(kind=i4b), parameter :: itermaxfit = 100
  real(kind=r8b),    parameter :: tol = 1.d-8
! -----------------------------------------
  integer(kind=i4b) :: ie, i, j, m, n, p, iter, bestiter, bestnum, bestden, maxie, maxi, maxj
  complex(kind=c8b) :: gf(nlmsb,nlmsb), gfdata(nesusc,nlmsb,nlmsb)
  complex(kind=c8b) :: y, dy(numd+dend), efac
  complex(kind=c8b) :: a(nesusc,numd+dend), b(nesusc), numroots(numd), denroots(dend)
  complex(kind=c8b) :: dev(nesusc), fit(numd+dend), bestfit(numd+dend)
  real(kind=r8b)    :: lambda, weights(nesusc), maxdev, maxij, avgij, oldrms, newrms
  integer(kind=i4b) :: ipiv(numd+dend), info

  gffit(:,:,:,ia,ja) = 0.d0
! Fill in array of data
  do ie=1,nesusc
    call projected_gf(ie,ia,ja,gf,onsite,struct)
!    call symmetrize(nlmsb,gf,gfilter)
    gfdata(ie,:,:) = gf
  end do
! Get coefficients of fit
  avgij = 0.d0; maxij = -1.d0
  do j=1,nlmsba(ja)
  do i=1,nlmsba(ia)
    maxdev = maxval(abs(gfdata(1:nesusc,i,j)))
!   **************************
    if (maxdev > gfilter) then
!   **************************
!    write(*,'("fit ia=",i4," lmsb=",6i4)') ia, i2lmsb(:,i,ia), i2lmsb(:,j,ia)
!   --------------------------------------------------------------------
!   initialize fitting parameters
!    lambda = 1.d0
!    bestfit = 1.d0
!    fit = bestfit
!   compute initial value for rms
!    oldrms = 0.d0
!    do ie=1,nesusc
!      weights(ie) = 1.d0
!      call ratval2(numd,dend,bestfit,eshift,esusc(ie),y,dy)
!      oldrms = oldrms + weights(ie)*abs(gfdata(ie,i,j) - y)**2
!    end do
!    oldrms = sqrt(oldrms)/nesusc
!    write(iodb,'("ratfit_gf2: i,j=",2i4," rms=",2es16.8)') i, j, oldrms
!   --------------------------------------------------------------------
!   Now iterate to improve fit parameters
    do iter=1,itermaxfit
      a = 0.d0; b = 0.d0
      do ie=1,nesusc
        weights(ie) = 1.d0  ! unbiased fit
!     rhs: weighted
        b(ie) = weights(ie)*gfdata(ie,i,j)
!     numerator loop
        efac = 1.d0
        do p=1,numd
          a(ie,p) = weights(ie)*efac
          efac = efac*(esusc(ie) - eshift)
        end do
!     denominator loop
        efac = 1.d0
        do p=1,dend-1
          efac = efac*(esusc(ie) - eshift)
          a(ie,numd+p) = -efac*b(ie)
        end do
      end do
!     Solve for new parameters
      call least_squares(nesusc,numd+dend-1,a,b,maxdev)
!      write(iodb,'("SVD err=",2es16.8)') maxdev
!     Was there an improvement?
      fit(1:numd) = b(1:numd)
      fit(numd+1) = 1.d0
      fit(2+numd:numd+dend) = b(1+numd:numd+dend-1)
      newrms = 0.d0
      do ie=1,nesusc
        weights(ie) = 1.d0
        call ratval2(numd,dend,fit,eshift,esusc(ie),y,dy)
        newrms = newrms + weights(ie)*abs(gfdata(ie,i,j) - y)**2
      end do
      newrms = sqrt(newrms)/nesusc
      bestfit = fit
    end do
    write(iodb,'("ratfit_gf2: i,j=",2i4," rms=",2es16.8)') i, j, newrms
!   --------------------------------------------------------------------
    maxdev = 0.d0
    do ie=1,nesusc
      call ratval2(numd,dend,bestfit,eshift,esusc(ie),y,dy)
      lambda = abs(gfdata(ie,i,j) - y)
      if (lambda > maxdev) maxdev = lambda
    end do
    gffit(:,i,j,ia,ja) = bestfit
    write(iodb,'("ratfit_gf2: i,j=",2i4," maxdev",2es16.8)') i, j, maxdev
!   ******
    end if
!   ******
  end do
  end do
!  write(iodb,'("ratfit_gf2: ia,ja, i,j=",8i4," maxij=",2es16.8)') ia, ja, i2lmsb(:,maxi,ia), i2lmsb(:,maxj,ja), maxij
! All done!
  end subroutine ratfit_gf2

