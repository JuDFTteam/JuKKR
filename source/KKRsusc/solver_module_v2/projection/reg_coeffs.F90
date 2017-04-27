  subroutine reg_coeffs(natomd,lmaxd,nrmaxd,pz,fz,is,ie,ek,tmat)
! Energy dependent projection coefficients of regular scattering wfns
! Expects solutions of the ASA SRA problem, no m-dependence of pz or fz
! BEWARE OF THE SPECIAL TREATMENT OF THE RADIAL MESH, IR=1 NOT USED
  use global
  use bessel_new

  implicit none

! --> dimensions of the wavefunction arrays
  integer(kind=i4b), intent(in)    :: natomd, lmaxd, nrmaxd
! --> regular scattering solutions
  complex(kind=c8b), intent(in)    :: pz(nrmaxd,0:lmaxd,natomd), fz(nrmaxd,0:lmaxd,natomd)
! --> spin channel
  integer(kind=i4b), intent(in)    :: is
! --> current energy point
  integer(kind=i4b), intent(in)    :: ie
! --> defined value of the square-root of the energy
  complex(kind=c8b), intent(in)    :: ek
! --> tmatrix-in the projection basis
  complex(kind=c8b), intent(inout) :: tmat(0:nlmax,0:nlmax,natomd) 
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nr, nr0, nr1, il, nb, ib, ir
  real(kind=r8b)    :: r(nrmax), dr(nrmax)
  complex(kind=c8b) :: rjl(nrmax,0:nlmax), tmatl(0:nlmax), work1(nrmax), work2(nrmax), norm1, norm2
  complex(kind=c8b), external :: radint

  do ia=1,nasusc
    ih  = iasusc(ia)
    nr0 = nrpts0(ia); nr1 = nrpts1(ia); nr = nr1 - nr0 + 1
    dr  = drproj(:,ia)
    r   = rmesh(:,ia)
!   r * J_l(sqrt(E)*r)
    do il=0,nlmax
      do ir=1,nr
        norm1 = ek*r(ir)
        rjl(ir,il) = r(ir)*bessj(il,norm1)
      end do
    end do
!   projection coefficients
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
      tmatl(il) = 0.d0
      if (nobasis(il,is,ia) .and. nb > 0) then
        write(*,'("reg_coeffs: no basis",3i4)') il, is, ia
        stop
      end if
!     projection on each basis function
!     basis functions assumed normalized to 1 with weight dr
      do ib=1,nb
        work1(1:nr) = pz(nr0:nr1,il,ih)*phiref(1:nr,ib,il,is,ia)
        norm1 = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
        pzc(ib,il,is,ia) = norm1
!     -------------------------------------------------------------
        if (isra == 1) then
          work2(1:nr) = fz(nr0:nr1,il,ih)*phiref(1:nr,ib,il,is,ia)
          norm2 = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
          fzc(ib,il,is,ia) = norm2
        end if
!     -------------------------------------------------------------
!       t-matrix in projection basis
        work1(1:nr) = rjl(1:nr,il)*(vr(1:nr,ia) + (2*is-3)*br(1:nr,ia) - 2.d0*zat(ia)/r(1:nr))*phiref(1:nr,ib,il,is,ia)
        tmatl(il) = tmatl(il) + pzc(ib,il,is,ia)*radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
      end do
!     check how much of the function is left
      if (nb > 0) then
        work1(1:nr) = pz(nr0:nr1,il,ih)
        do ib=1,nb
          norm1 = pzc(ib,il,is,ia)
          work1(1:nr) = work1(1:nr) - norm1*phiref(1:nr,ib,il,is,ia)
        end do
        work1 = work1*conjg(work1)
        norm1 = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
        work1(1:nr) = conjg(pz(nr0:nr1,il,ih))*pz(nr0:nr1,il,ih)
        norm1 = norm1/radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!     -------------------------------------------------------------
        norm2 = 0.d0
      if (isra == 1) then
        work2(1:nr) = fz(nr0:nr1,il,ih)
        do ib=1,nb
          norm2 = fzc(ib,il,is,ia)
          work2(1:nr) = work2(1:nr) - norm2*phiref(1:nr,ib,il,is,ia)
        end do
        work2 = work2*conjg(work2)
        norm2 = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
        work2(1:nr) = conjg(fz(nr0:nr1,il,ih))*fz(nr0:nr1,il,ih)
        norm2 = norm2/radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
      end if
!     -------------------------------------------------------------
        if (lhdio) write(iodb,'("reg_coeffs",4i4,2es16.3)') ie, ia, is, il, sqrt(abs(norm1)), sqrt(abs(norm2))
      end if
!     -------------------------------------------------------------
!     t-matrix test
!      norm1 = tmat(il,ih); norm2 = tmatl(il)
!      write(*,'("ie,ia,is,il=",4i4,"  tmat,tmatl,diff=",4es16.8,f6.1)') ie, ia, is, il, norm1, norm2, 1.d2*abs(norm1-norm2)/abs(norm2)
!     Save the t-matrix in projection basis
      if (ltmatproj)  tmat(il,il,ih) = tmatl(il)
!     -------------------------------------------------------------
    end do
  end do
! coefficients computed
  noregcoeffs(is,ie) = .false.
! All done
  end subroutine reg_coeffs
