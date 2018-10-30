  subroutine irr_coeffs(natomd,lmaxd,nrmaxd,pz,qz,fz,sz,is,ie,ek)
! Energy dependent projection coefficients of the on-site part of the GF
! Expects solutions of the ASA SRA problem, no m-dependence of pz or fz
! BEWARE OF THE SPECIAL TREATMENT OF THE RADIAL MESH, IR=1 NOT USED
  use global

  implicit none

! --> dimensions of the wavefunction arrays
  integer(kind=i4b), intent(in) :: natomd, lmaxd, nrmaxd
! --> regular scattering solutions
  complex(kind=c8b), intent(in) :: pz(nrmaxd,0:lmaxd,natomd), fz(nrmaxd,0:lmaxd,natomd)
! --> irregular scattering solutions
  complex(kind=c8b), intent(in) :: qz(nrmaxd,0:lmaxd,natomd), sz(nrmaxd,0:lmaxd,natomd)
! --> spin channel
  integer(kind=i4b), intent(in) :: is
! --> current energy point
  integer(kind=i4b), intent(in) :: ie
! --> defined value of the square-root of the energy
  complex(kind=c8b), intent(in) :: ek
! -----------------------------------------------------------------
  logical, parameter :: separable = .false.
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nr, nr0, nr1, il, nb
  integer(kind=i4b) :: ib, jb, ir
  real(kind=r8b)    :: dr(nrmax)
  complex(kind=c8b) :: work1(nrmax), work2(nrmax), norm, norm2
  complex(kind=c8b), external :: radint

  do ia=1,nasusc
    ih = iasusc(ia)
    nr0 = nrpts0(ia); nr1 = nrpts1(ia); nr = nr1 - nr0 + 1
    dr = drproj(:,ia)
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
      if (nobasis(il,is,ia) .and. nb > 0) then
        write(*,'("irr_coeffs: no basis",3i4)') il, is, ia
        stop
      end if
!   ====================================================================
      if (separable) then
!   expansion of the irregular solutions
      do jb=1,nb
        work1(1:nr) = qz(nr0:nr1,il,ih)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
        if (isra == 1) then
          work2(1:nr) = sz(nr0:nr1,il,ih)*phiref(1:nr,jb,il,is,ia)
          norm2 = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
        end if
        do ib=1,nb
          pqc(ib,jb,il,is,ia) = norm*ek*pzc(ib,il,is,ia)
          if (isra == 1) then
            psc(ib,jb,il,is,ia) = norm2*ek*pzc(ib,il,is,ia)
            fqc(ib,jb,il,is,ia) = norm*ek*fzc(ib,il,is,ia)
            fsc(ib,jb,il,is,ia) = norm2*ek*fzc(ib,il,is,ia)
          end if
        end do
      end do
!   symmetrize
!      pqc(1:nb,1:nb,il,is,ia) = 0.5d0*(pqc(1:nb,1:nb,il,is,ia) - transpose(pqc(1:nb,1:nb,il,is,ia)))
!   ====================================================================
      else
!   ====================================================================
!   double expansion of the on-site product
!   basis functions assumed normalized to 1 with weight dr
      do jb=1,nb
      do ib=1,nb
!     integral over r' with switching between P and Q
!     ir <-> r; work2 stores data for r'
        do ir=1,nr
!       r' <= r --> Q(r)P(r')
          work2(1:ir) = qz(nr0+ir-1,il,ih)*pz(nr0:nr0+ir-1,il,ih)
!       r' > r --> P(r)Q(r')
          work2(ir+1:nr) = pz(nr0+ir-1,il,ih)*qz(nr0+ir:nr1,il,ih)
          work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
!       integrate over r'
          norm = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
          work1(ir) = norm
        end do
!     multiply by basis function and integrate over r
        work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!     include the sqrt of the energy
        pqc(ib,jb,il,is,ia) = norm*ek
!     -------------------------------------------------------------
      if (isra == 1) then
!     integral over r' with switching between P and S
        do ir=1,nr
!       r' <= r --> S(r)P(r')
          work2(1:ir) = sz(nr0+ir-1,il,ih)*pz(nr0:nr0+ir-1,il,ih)
!       r' > r --> P(r)S(r')
          work2(ir+1:nr) = pz(nr0+ir-1,il,ih)*sz(nr0+ir:nr,il,ih)
          work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
          norm = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
          work1(ir) = norm
        end do
!     multiply by basis function and integrate over r'
        work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!     include the sqrt of the energy
        psc(ib,jb,il,is,ia) = norm*ek
!     -------------------------------------------------------------
!     integral over r' with switching between F and Q
        do ir=1,nr
!       r' <= r --> Q(r)F(r')
          work2(1:ir) = qz(nr0+ir-1,il,ih)*fz(nr0:nr0+ir-1,il,ih)
!       r' > r --> F(r)Q(r')
          work2(ir+1:nr) = fz(nr0+ir-1,il,ih)*qz(nr0+ir:nr,il,ih)
          work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
          norm = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
          work1(ir) = norm
        end do
!     multiply by basis function and integrate over r'
        work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!     include the sqrt of the energy
        fqc(ib,jb,il,is,ia) = norm*ek
!     -------------------------------------------------------------
!     integral over r' with switching between F and S
        do ir=1,nr
!       r' <= r --> S(r)F(r')
          work2(1:ir) = sz(nr0+ir-1,il,ih)*fz(nr0:nr0+ir-1,il,ih)
!       r' > r --> F(r)S(r')
          work2(ir+1:nr) = fz(nr0+ir-1,il,ih)*sz(nr0+ir:nr,il,ih)
          work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
          norm = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
          work1(ir) = norm
        end do
!     multiply by basis function and integrate over r'
        work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
        norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!     include the sqrt of the energy
        fsc(ib,jb,il,is,ia) = norm*ek
      end if
!     -------------------------------------------------------------
      end do
      end do
!     ==================================================================
      end if
!     ==================================================================
    end do
  end do
! coefficients computed
  noirrcoeffs = .false.
! All done
  end subroutine irr_coeffs
