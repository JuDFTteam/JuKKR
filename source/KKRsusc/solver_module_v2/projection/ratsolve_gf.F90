  subroutine ratsolve_gf(ia,ja,onsite,struct)
! Rational function approximation
! TEST VERSION

  implicit none

! Block of GF
  integer(kind=i4b), intent(in)  :: ia, ja
! To add onsite or structural parts of GF
  logical,           intent(in)  :: onsite, struct
! -----------------------------------------
  real(kind=r8b),    parameter :: tol = 1.d-8
  integer(kind=i4b), parameter :: itermax = 1
! -----------------------------------------
  integer(kind=i4b) :: ie, i, j, p, iter, bestiter, bestnum, bestden, maxie, maxi, maxj, inum, iden
  complex(kind=c8b) :: gf(nlmsb,nlmsb), gfdata(nesusc,nlmsb,nlmsb)
  complex(kind=c8b) :: phase, y, efac
  complex(kind=c8b) :: dev(nesusc), numroots(numd), denroots(dend)
  real(kind=r8b)    :: weights(nesusc), avgdev, maxdev, mindev, maxij, avgij, rms
  complex(kind=c8b) :: a(nesusc,nesusc), af(nesusc,nesusc), b(nesusc), x(nesusc), work(2*nesusc)
  real(kind=r8b)    :: r(nesusc), c(nesusc), rwork(2*nesusc), rcond, ferr, berr
  integer(kind=i4b) :: info, ipiv(nesusc)
  character*1       :: equed = 'B'

  gffit(:,:,:,ia,ja) = 0.d0
  if (1 + numd + dend /= nesusc) stop '1 + numd + dend /= nesusc'
! Fill in array of data
  do ie=1,nesusc
    call projected_gf(ie,ia,ja,gf,onsite,struct)
!    call symmetrize(nlmsb,gf,gfilter)
    gfdata(ie,:,:) = gf
  end do
! Get coefficients of fit
  maxij = -1.d0
  do j=1,nlmsba(ja)
  do i=1,nlmsba(ia)
    maxdev = maxval(abs(gfdata(1:nesusc,i,j)))
!   **************************
    if (maxdev > gfilter) then
!   **************************
!    write(*,'("fit ia=",i4," lmsb=",6i4)') ia, i2lmsb(:,i,ia), i2lmsb(:,j,ia)
!   --------------------------------------------------------------------
      do ie=1,nesusc
!     rhs:
        b(ie) = gfdata(ie,i,j)
!     numerator loop
        efac = 1.d0
        a(ie,1) = efac
        do p=1,numd
          efac = efac*(esusc(ie) - eshift)
          a(ie,1+p) = efac
        end do
!     denominator loop
        efac = 1.d0
        do p=1,dend
          efac = efac*(esusc(ie) - eshift)
          a(ie,1+numd+p) = -efac*gfdata(ie,i,j)
        end do
      end do
!   solve sys of eqs with iterative improvement
!      call zgesv(nesusc,1,a,nesusc,ipiv,b,nesusc,info)
!      x = b
      call zgesvx('N','N',nesusc,1,a,nesusc,af,nesusc,ipiv,equed,r,c,b,nesusc,x,nesusc,rcond,ferr,berr,work,rwork,info)
!      write(iodb,'("info=",i4,3es10.1)') info, rcond, ferr, berr
      if (info /= 0 .and. info /= nesusc + 1) stop 'ratsolve_gf: failure in zgesvx'
!   compute deviations on fitted points
      do ie=1,nesusc
        call ratval(numd,dend,x,eshift,esusc(ie),y)
        dev(ie) = y - gfdata(ie,i,j)
        weights(ie) = abs(dev(ie))
      end do
      rms = maxval(weights)
!      call rat_roots(eshift,inum,iden,fit,numroots,denroots)
!      if (any(aimag(denroots(1:iden)) > 0.d0)) avgdev = 10.d0*avgdev
!      write(iodb,'("ratfit_gf: inum, iden=",2i4," maxdev, avgdev=",2es16.8)') inum, iden, maxdev, avgdev
!   best fit
!      if (mindev > tol) write(iodb,'("ratfit_gf: ia,ja, i,j=",8i4," maxdev=",3i4,es16.8)') ia, ja, i2lmsb(:,i), i2lmsb(:,j), bestiter, bestnum, bestden, mindev
!   --------------------------------------------------------------------
      gffit(:,i,j,ia,ja) = x
!      write(iodb,'("ratsolve_gf: ia,ja, i,j=",8i4," rms=",2es16.8)') ia, ja, i2lmsb(:,i,ia), i2lmsb(:,j,ja), rms
      if (rms > maxij) then
        maxi = i; maxj = j; maxij = rms
      end if
!   ******
    end if
!   ******
  end do
  end do
  write(iodb,'("ratsolve_gf: ia,ja, i,j=",8i4," maxij=",2es16.8)') ia, ja, i2lmsb(:,maxi,ia), i2lmsb(:,maxj,ja), maxij
! All done!
  end subroutine ratsolve_gf

