  subroutine overlaps_susc(ia,gfsum)
 Analysis of the product basis for the susceptibility

  implicit none

  integer(kind=i4b), intent(in) :: ia
  complex(kind=c8b), intent(in) :: gfsum(nlmsb,nlmsb)
! ---------------------------------------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-3
  integer(kind=i4b) :: i3(3), nb, nr, i, j, lm, ib, ievals
  integer(kind=i4b) :: q1, b1, lm1, s1
  integer(kind=i4b) :: q2, b2, lm2, s2
  integer(kind=i4b) :: q3, b3, lm3, s3
  integer(kind=i4b) :: q4, b4, lm4, s4
  integer(kind=i4b) :: info, lrwork, liwork
  integer(kind=i4b), allocatable :: iwork(:), lmsb2i2(:,:)
  real(kind=r8b),    allocatable :: big_overlaps(:,:), rwork(:), evals(:)
  complex(kind=c8b) :: work(nrmax), norm, rho2, rho2proj, rho, rhoc
  real(kind=r8b)    :: dr(nrmax), doublegaunt, doublegaunt0


  nb = nlmsba(ia)**2
  allocate(big_overlaps(nb,nb),lmsb2i2(2,nb))
  big_overlaps = 0.d0; rho2 = 0.d0; rho2proj = 0.d0
  nr = nrpts(ia); dr = drmesh(:,ia)
  j = 0
  do q4=1,nlmsba(ia)
    i3 = i2lmsb(:,q4,ia)
    b4 = i3(1); lm4 = i3(2); s4 = i3(3)
    do q3=1,nlmsba(ia)
      i3 = i2lmsb(:,q3,ia)
      b3 = i3(1); lm3 = i3(2); s3 = i3(3)
      j = j + 1
      lmsb2i2(:,j) = (/q3,q4/)
      i = 0
      do q2=1,nlmsba(ia)
        i3 = i2lmsb(:,q2,ia)
        b2 = i3(1); lm2 = i3(2); s2 = i3(3)
        do q1=1,nlmsba(ia)
          i3 = i2lmsb(:,q1,ia)
          b1 = i3(1); lm1 = i3(2); s1 = i3(3)
          i = i + 1
!         -----------------------------------------------------------------------------
          if (s1 == s3 .and. s2 == s4) then
            doublegaunt0 = 0.d0
            do lm=1,lmmax0  ! truncated Ylm expansion
              doublegaunt0 = doublegaunt0 + gaunt(lm1,lm2,lm)*gaunt(lm3,lm4,lm)
            end do
            doublegaunt = doublegaunt0
            do lm=lmmax0+1,lmmax2  ! full Ylm expansion
              doublegaunt = doublegaunt + gaunt(lm1,lm2,lm)*gaunt(lm3,lm4,lm)
            end do
            if (abs(doublegaunt) > ylmtol) then
              write(iodb,'("big_overlaps: ia=",i4," nonzero q=",4i8,es16.8)') ia, q1, q2, q3, q4, doublegaunt
              work = phiref(:,b1,i2lm(2,lm1),s1,ia)*phiref(:,b2,i2lm(2,lm2),s2,ia)
              work = work*phiref(:,b3,i2lm(2,lm3),s3,ia)*phiref(:,b4,i2lm(2,lm4),s4,ia)
              work(1:nr) = work(1:nr)/rmesh(1:nr,ia)**2
              norm = radint(nr,work(1:nr),dr(1:nr))
              rho2 = rho2 + real(norm)*doublegaunt*gfsum(q1,q2)*conjg(gfsum(q4,q3))
              if (abs(doublegaunt0) > ylmtol) big_overlaps(i,j) = real(norm)*doublegaunt0
            end if
          end if
!         -----------------------------------------------------------------------------
        end do
      end do
    end do
  end do
! diagonalize; eigenvalues in ascending order
  lrwork = 2*nb*(1 + 3*nb) + 1
  liwork = 5*nb + 3
  allocate(evals(nb),rwork(lrwork),iwork(liwork))
  call dsyevd('V','U',nb,big_overlaps,nb,evals,rwork,lrwork,iwork,liwork,info)
  ievals = 0
  do ib=1,nb
    write(iodb,'("big_overlaps: ia,ib=",i4,i8," evals(ib)=",es16.8)') ia, ib, evals(ib)
    if (abs(evals(ib)) > tol*evals(nb)) then
      write(iodb,'("evec(ib)  ib ilm  is  jb jlm  js")')
      if (sum(big_overlaps(:,ib)) < 0.d0) big_overlaps(:,ib) = -big_overlaps(:,ib)
      do i=1,nb
        if (abs(big_overlaps(i,ib)) > 1.d-4) then
          write(iodb,'(f8.4,1000i4)') big_overlaps(i,ib), i2lmsb(:,lmsb2i2(1,i),ia), i2lmsb(:,lmsb2i2(2,i),ia)
        end if
      end do
      ievals = ievals + 1 
    end if
  end do
  write(iodb,'("big_overlaps: ia=",i4," evals > tol=",es10.1,i8)') ia, tol, ievals
! check how much of the density matrix is lost
  rho2proj = 0.d0
  do ib=nb-ievals+1,nb
    rho = 0.d0; rhoc = 0.d0
    i = 0
    do q2=1,nlmsba(ia)
      do q1=1,nlmsba(ia)
        i = i + 1
        rho  = rho  + big_overlaps(i,ib)*gfsum(q1,q2)
        rhoc = rhoc + big_overlaps(i,ib)*conjg(gfsum(q2,q1))
      end do
    end do
    rho2proj = rho2proj + evals(ib)*rho*rhoc
  end do
  write(iodb,'("big_overlaps: ia=",i4," rho2=",2es16.8,"  rho2proj=",2es16.8,"  diff=",es16.8)') ia, rho2, rho2proj, abs(rho2 - rho2proj)
  deallocate(big_overlaps,evals,rwork,iwork)
! All done!
  end subroutine overlaps_susc

