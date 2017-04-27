  subroutine tmat_proj(ia,ie,is,tmat)
! compute t-matrix using basis
! t_l = \int dr r j_l V R_l
  use global
  use bessel_new

  implicit none

! which atom
  integer(kind=i4b), intent(in)    :: ia
! which energy                   
  integer(kind=i4b), intent(in)    :: ie
! which spin component           
  integer(kind=i4b), intent(in)    :: is
! t-matrix in projection basis
  complex(kind=c8b), intent(inout) :: tmat(lmmax,lmmax)
! ----------------------------------------------------------------------
  real(kind=r8b)    :: dr(nrmax), r(nrmax)
  complex(kind=c8b) :: rjl(nrmax,0:nlmax), work(nrmax), tmatl(0:nlmax), norm, norm2
  integer(kind=i4b) :: nr, i, i3(3), ilm, il, ib, ir, js, jlm
  complex(kind=c8b), external :: radint

! radial mesh info
  nr = nrpts(ia)
  dr(1:nr) = drmesh(1:nr,ia)
  r(1:nr)  = rmesh(1:nr,ia)
! r * J_l(sqrt(E)*r)
  do il=0,nlmax
    do ir=1,nr
      norm = eksusc(ie)*r(ir)
      rjl(ir,il) = r(ir)*bessj(il,norm)
    end do
  end do
! compute t-matrix
  tmatl(:) = 0.d0
  do i=1,nlmsba(ia)
    i3 = i2lmsb(:,i,ia)
    ib = i3(1); ilm = i3(2); js = i3(3)
    if (is == js) then
!     \int dr r J_l(r;E) V(r) R_l(r;E)
      il = i2lm(2,ilm)
      work(1:nr) = rjl(1:nr,il)*(vr(1:nr,ia) + (2*is-3)*br(1:nr,ia) - 2.d0*zat(ia)/r(1:nr))*phiref(1:nr,ib,il,is,ia)
      norm = radint(nr,work(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
      if (i2lm(1,ilm) == 0) tmatl(il) = tmatl(il) + norm*pzr(i,lms2i(ilm,is),ia,ie)
    end if
  end do
! update t-matrix
  do jlm=1,lmmax
    do ilm=1,lmmax
      if (ilm == jlm) then
        norm = tmat(ilm,ilm); norm2 = tmatl(i2lm(2,ilm))
        if (i2lm(1,ilm) == 0) write(*,'("ie,ia,is,ilm=",4i4,"  tmat,tmatl,diff=",4es16.8,f6.1)') ie, ia, is, ilm, norm, norm2, 1.d2*abs((norm2-norm))/abs(norm)
        tmat(ilm,jlm) = tmatl(i2lm(2,ilm))
      else
        tmat(ilm,jlm) = 0.d0
      end if
    end do
  end do
! test output
!  if (ia == 1 .and. is == 1 .and. ie == 1) then
!    open(file='bessel.dat',unit=iofile2,status='replace')
!    do i=1,1001
!      norm = 1.d-2*(i-1.d0)
!      write(iofile2,'(100es16.8)') norm, (bessn(il,norm),il=0,nlmax)
!    end do
!    close(iofile2)
!  end if
! All done  
  end subroutine tmat_proj
