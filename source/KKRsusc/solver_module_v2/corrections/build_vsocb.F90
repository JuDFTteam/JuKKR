  subroutine build_vsocb(ie,ia,vsoc,vsocb,u,ipick)
! assembles the SOC potential in the projection basis
  use global !, only: i4b, r8b, c8b, nrmax, nrpts, drmesh, phiref, nlmsba, nlmsb, nlms, i2lmsb, i2lms, i2lm, ldots, atol, pauli, lorb

  implicit none

! Which energy
  integer(kind=i4b), intent(in)    :: ie
! Which atom
  integer(kind=i4b), intent(in)    :: ia
! Spherical spin-averaged radial SOC potential
  complex(kind=c8b), intent(in)    :: vsoc(nrmax)
! SOC potential in the basis (added to)
  complex(kind=c8b), intent(inout) :: vsocb(nlmsb,nlmsb)
! Direction of spin moment
  real(kind=r8b),    intent(in)    :: u(3)
! Transverse, longitudinal or full L.S operator
  integer(kind=i4b), intent(in)    :: ipick
! -----------------------------------------------------------------
  integer(kind=i4b) :: is, js, ib, jb, nbi, nbj, i3(3), i2(2)
  integer(kind=i4b) :: jm, jl, im, il, jlm, ilm, i, j
  integer(kind=i4b) :: nr
  real(kind=r8b)    :: dr(nrmax), socparam(nbmax,nbmax,nsmax,nsmax,0:nlmax)
  complex(kind=c8b) :: work(nrmax), norm
  integer(kind=i4b) :: imax, jmax
  complex(kind=c8b), external :: radint

! use the structure of L.S
! also use the fact that the basis was constructed per l channel
!  maxnorm = 0.d0
! ----------------------------------------------------------------------
  nr = nrpts(ia)
  dr(1:nr) = drmesh(1:nr,ia)
! SOC parameters
!  if (ie == nescf) then
!  do il=0,nlmax
!    write(*,'("SOC params for ia=",i4,"  il=",i4)') ia, il
!    do js=1,nsmax
!      do is=1,nsmax
!        write(*,'(" is, js=",2i4)') is, js
!        do jb=1,iwsusc(il,js,ia)
!          do ib=1,iwsusc(il,is,ia)
!            work(1:nr) = phiref(1:nr,ib,il,is,ia)*vsoc(1:nr)*phiref(1:nr,jb,il,js,ia)
!            socparam(ib,jb,is,js,il) = real(radint(nr,work,dr,npanat(ia),ircutat(:,ia)))
!            write(*,'(2i4,2f16.8)') ib, jb, socparam(ib,jb,is,js,il), socparam(ib,jb,is,js,il)*13610d0
!          end do
!        end do
!      end do
!    end do
!  end do
!  end if
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
!      if (abs(ldots(i,j,ia)) > atol) then
      if (il == jl) then
!       revise the basis
        work(1:nr) = phiref(1:nr,ib,il,is,ia)*vsoc(1:nr)*phiref(1:nr,jb,il,js,ia)
        norm = radint(nr,work,dr,npanat(ia),ircutat(:,ia))
!        if (abs(norm) > maxnorm) then
!          maxnorm = abs(norm)
!          imax = i; jmax = j
!        end if
!       full L.S
        if (ipick == 1) vsocb(i,j) = vsocb(i,j) + norm*sum(lorb(ilm,jlm,:)*pauli(is,js,1:3))
!       longitudinal part 
        if (ipick == 2) vsocb(i,j) = vsocb(i,j) + norm*sum(lorb(ilm,jlm,:)*u)*sum(pauli(is,js,1:3)*u)
!       transverse part
        if (ipick == 3) vsocb(i,j) = vsocb(i,j) + norm*(sum(lorb(ilm,jlm,:)*pauli(is,js,1:3)) - sum(lorb(ilm,jlm,:)*u)*sum(pauli(is,js,1:3)*u))
      end if
!   ----------------------------------------------------------------
    end do
  end do
!  write(iodb,'(" SOC kernel, ia=",i4)') ia
!  write(iodb,'("vsocb norm:",6i4,2es16.8)') i2lmsb(:,imax,ia), i2lmsb(:,jmax,ia), maxnorm
!  write(iodb,'("  ib ilm  is, jb jlm  js, matrix element")')
!  do j=1,nlmsba(ia)
!    do i=1,nlmsba(ia)
!      if (abs(vsocb(i,j)) > 1.d-8) then
!        write(iodb,'(6i4,2es10.2)') i2lmsb(:,i,ia), i2lmsb(:,j,ia), vsocb(i,j)
!      end if
!    end do
!  end do
! All done

end subroutine build_vsocb
