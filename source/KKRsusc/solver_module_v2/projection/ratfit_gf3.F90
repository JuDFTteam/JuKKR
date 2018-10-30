  subroutine ratfit_gf3(ia,ja,onsite,struct)
! Rational function fit
! TEST VERSION

  implicit none

! Block of GF
  integer(kind=i4b), intent(in)  :: ia, ja
! To add onsite or structural parts of GF
  logical,           intent(in)  :: onsite, struct
! -----------------------------------------
  real(kind=r8b),    parameter :: tol = 1.d-8
  integer(kind=i4b), parameter :: itermax = 10
! -----------------------------------------
  integer(kind=i4b) :: ie, i, j, p, iter, bestiter, maxie, maxi, maxj
  complex(kind=c8b) :: gf(nlmsb,nlmsb), gfdata(nesusc,nlmsb,nlmsb)
  complex(kind=c8b) :: y, efac
  complex(kind=c8b) :: a(nesusc,numd+dend), b(nesusc), numroots(numd), denroots(dend)
  complex(kind=c8b) :: dev(nesusc), fit(numd+dend), bestfit(numd+dend), weights(nesusc)
  real(kind=r8b)    :: avgdev, maxdev, mindev, maxij, avgij, rms

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
    rms = 1.d99
!   --------------------------------------------------------------------
!   first pass: all points weighted in the same way
    weights = 1.d0
    do iter=1,itermax
!   fill in design matrix; each equation is weighted
      do ie=1,nesusc
!        weights(ie) = weights(ie)/aimag(esusc(ie))
!        weights(ie) = weights(ie)/(minval(abs(esusc(ie) - escf)) + 1.d-4)
!     rhs: weighted
        b(ie) = weights(ie)*gfdata(ie,i,j)
!     numerator loop (degree numd-1)
        efac = 1.d0
        do p=1,numd
          a(ie,p) = weights(ie)*efac
          efac = efac*(esusc(ie) - eshift)
        end do
!     denominator loop (degree dend-1)
        efac = 1.d0
        do p=1,dend-1
!          efac = efac*(esusc(ie) - eshift)
          a(ie,numd+p) = -efac*b(ie)
          efac = efac*(esusc(ie) - eshift)
        end do
        b(ie) = efac*b(ie)
      end do
!   fit
      call least_squares(nesusc,numd+dend-1,a,b,maxdev)
!      write(iodb,'("SVD err=",2es16.8)') maxdev
      fit = 0.d0
      do p=1,numd
        fit(p) = b(p)
      end do
!      fit(numd+1) = 1.d0
      do p=1,dend-1
!      do p=2,iden
        fit(numd+p) = b(numd+p)
      end do
      fit(numd+dend) = 1.d0
!   compute deviations on fitted points
      do ie=1,nesusc
        call ratval(numd,dend,fit,eshift,esusc(ie),y)
        dev(ie) = y - gfdata(ie,i,j)
      end do
      avgdev = sum(abs(dev))/nesusc
      maxdev = maxval(abs(dev))
!      if (any(aimag(denroots(1:iden)) > 0.d0)) avgdev = 10.d0*avgdev
!   best fit so far?
      if (avgdev < rms) then
        rms = avgdev; bestiter = iter
        bestfit = fit
      end if
!      write(iodb,'("ratfit_gf:",4i4," maxdev, avgdev=",2es16.8)') i, j, inum, iden, maxdev, avgdev
!      if (i == j) then
!        write(iodb,'(2es16.8)') fit
!      end if
!     weights are defined as the reciprocal of the denominator
      do ie=1,nesusc
!     denominator loop (degree dend-1)
        efac = 1.d0; y = 0.d0
        do p=1,dend-1
!          efac = efac*(esusc(ie) - eshift)
          y = y + efac*bestfit(numd+p)
          efac = efac*(esusc(ie) - eshift)
        end do
        weights(ie) = 1.d0/(y + efac)
      end do
    end do
!   best fit
!      if (mindev > tol) write(iodb,'("ratfit_gf: ia,ja, i,j=",8i4," maxdev=",3i4,es16.8)') ia, ja, i2lmsb(:,i), i2lmsb(:,j), bestiter, bestnum, bestden, mindev
!   --------------------------------------------------------------------
    gffit(:,i,j,ia,ja) = bestfit
! check roots
    dev(1:numd) = bestfit(1:numd)
    dev(numd+1:dend) = bestfit(numd+1:dend)
    call rat_roots(eshift,numd,dend,dev,numroots,denroots)
    write(iodb,'("ratfit_gf: ia,ja, i,j=",8i4," num, den, rms=",3i4,2es16.8)') ia, ja, i2lmsb(:,i,ia), i2lmsb(:,j,ja), bestiter, numd, dend, rms
    if (rms > maxij) then
      maxi = i; maxj = j; maxij = rms
    end if
!    avgij = avgij + mindev
!   ******
    end if
!   ******
  end do
  end do
  write(iodb,'("ratfit_gf: ia,ja, i,j=",8i4," maxij=",2es16.8)') ia, ja, i2lmsb(:,maxi,ia), i2lmsb(:,maxj,ja), maxij
! All done!
  end subroutine ratfit_gf3

