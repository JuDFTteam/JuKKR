  subroutine irr_coeffs_soc(natomd,lmsize,nrmaxd,pz,qz,fz,sz,pzleft,qzleft,fzleft,szleft,ie,ek)
! Energy dependent projection coefficients of the on-site part of the GF
! Expects solutions of the ASA SRA problem, no m-dependence of pz or fz
! BEWARE OF THE SPECIAL TREATMENT OF THE RADIAL MESH, IR=1 NOT USED
! CAREFUL --> lmsize = 2*(lmaxd+1)**2
! Separable option is deleted here
  use global

  implicit none

! --> dimensions of the wavefunction arrays
  integer(kind=i4b), intent(in) :: natomd, lmsize, nrmaxd
! --> regular scattering solutions right
  complex(kind=c8b), intent(in) :: pz(nrmaxd,lmsize,lmsize,natomd), fz(nrmaxd,lmsize,lmsize,natomd)
! --> regular scattering solution left
  complex(kind=c8b), intent(in) :: pzleft(nrmaxd,lmsize,lmsize,natomd), fzleft(nrmaxd,lmsize,lmsize,natomd)
! --> irregular scattering solutions
  complex(kind=c8b), intent(in) :: qz(nrmaxd,lmsize,lmsize,natomd), sz(nrmaxd,lmsize,lmsize,natomd)
! --> irregular scattering solution left
  complex(kind=c8b), intent(in) :: qzleft(nrmaxd,lmsize,lmsize,natomd), szleft(nrmaxd,lmsize,lmsize,natomd)
! --> current energy point
  integer(kind=i4b), intent(in) :: ie
! --> defined value of the square-root of the energy
  complex(kind=c8b), intent(in) :: ek
! ----------------------------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nr, nr0, nr1, nb1, nb2
  integer(kind=i4b) :: ib, jb, ir, i, j, k, m
  integer(kind=i4b) :: ilmsn, ilmsn2, ilmsn3, ilm, ilm2, ilm3, il, il2, il3, im, im2
  integer(kind=i4b) :: is, is2, is3 ! Spin channels
  real(kind=r8b)    :: dr(nrmax), mom_loc(3,nasusc)
  complex(kind=c8b) :: work1(nrmax), work2(nrmax), norm, norm2
  logical           :: rotate_onsite, rotate_regular 
  complex(kind=c8b), external :: radint

  ! Atoms loop
  do ia=1,nasusc
    ih = iasusc(ia)
    nr0 = nrpts0(ia); nr1 = nrpts1(ia); nr = nr1 - nr0 + 1
    dr = drproj(:,ia)
    do ilmsn2 = 1, lmsize 
      is2  = i2lms_new(2,ilmsn2)
      ilm2 = i2lms_new(1,ilmsn2)  
      il2  = i2lm(2,ilm2)
      im2  = i2lm(1,ilm2)
      do ilmsn = 1, lmsize
        is   = i2lms_new(2,ilmsn)
        ilm  = i2lms_new(1,ilmsn)
        il   = i2lm(2,ilm)
        im   = i2lm(1,ilm)
        nb1 = iwsusc(il,is,ia)
        nb2 = iwsusc(il2,is2,ia)
        if (nobasis(il,is,ia) .and. nb1 > 0) then
          write(*,'("reg_coeffs: no basis",6i4)') il, im, is, im2, is2, ia
          stop
        end if
!       SOC ---> Radial solution Rlm,m'(r) and Phil(r) basis
!       SOC does not couple l --> l' 
        if (il == il2) then

!         Construct the onsite green function                         
          ! ====================================================================
          ! double expansion of the on-site product
          ! basis functions assumed normalized to 1 with weight dr
          do jb=1,nb1
            do ib=1,nb1
              ! integral over r' with switching between P and Q
              ! ir <-> r; work2 stores data for r'
              do ir=1,nr
                ! initialize array for the loop
                work2 = 0.d0
                do ilmsn3 = 1, lmsize  
                  ilm3 = i2lms_new(1,ilmsn3)
                  il3  = i2lm(2,ilm3)
                  if (il3 == il .and. il3 == il2) then
                    ! r' <= r --> Q(r)P(r') and sum over s and l (SOC) 
                    work2(1:ir) = work2(1:ir) + qz(nr0+ir-1,ilmsn,ilmsn3,ih)*pzleft(nr0:nr0+ir-1,ilmsn2,ilmsn3,ih)
                    ! r' > r --> P(r)Q(r') and sum over s and l (SOC)
                    work2(ir+1:nr) = work2(ir+1:nr) + pz(nr0+ir-1,ilmsn,ilmsn3,ih)*qzleft(nr0+ir:nr,ilmsn2,ilmsn3,ih)
                  endif
                enddo ! ilmsn3 
                work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
                ! integrate over r'
                norm = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
                work1(ir) = norm
              enddo ! ir
              ! multiply by basis function and integrate over r
              work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is2,ia)
              norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
              i = lmsb2i(ib,ilm,is,ia)
              j = lmsb2i(jb,ilm2,is2,ia)           
              ! Store the onsite green function and include the sqrt of the energy
              gfpq(i,j,ia,ie) = norm*ek        
            enddo ! ib 
          enddo ! jb
!                 ! -------------------------------------------------------------
!                 if (isra == 1) then
!                   ! integral over r' with switching between P and S
!                   do ir=1,nr
!                     ! r' <= r --> S(r)P(r')
!                     work2(1:ir) = sz(nr0+ir-1,il,ih)*pz(nr0:nr0+ir-1,il,ih)
!                     ! r' > r --> P(r)S(r')
!                     work2(ir+1:nr) = pz(nr0+ir-1,il,ih)*sz(nr0+ir:nr,il,ih)
!                     work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
!                     norm = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!                     work1(ir) = norm
!                   end do
!                   ! multiply by basis function and integrate over r'
!                   work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
!                   norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!                   ! include the sqrt of the energy
!                   psc(ib,jb,il,is,ia) = norm*ek
!                   ! -------------------------------------------------------------
!                   ! integral over r' with switching between F and Q
!                   do ir=1,nr
!                     ! r' <= r --> Q(r)F(r')
!                     work2(1:ir) = qz(nr0+ir-1,il,ih)*fz(nr0:nr0+ir-1,il,ih)
!                     ! r' > r --> F(r)Q(r')
!                     work2(ir+1:nr) = fz(nr0+ir-1,il,ih)*qz(nr0+ir:nr,il,ih)
!                     work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
!                     norm = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!                     work1(ir) = norm
!                   end do
!                   ! multiply by basis function and integrate over r'
!                   work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
!                   norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!                   ! include the sqrt of the energy
!                   fqc(ib,jb,il,is,ia) = norm*ek
!                   ! -------------------------------------------------------------
!                   ! integral over r' with switching between F and S
!                   do ir=1,nr
!                     ! r' <= r --> S(r)F(r')
!                     work2(1:ir) = sz(nr0+ir-1,il,ih)*fz(nr0:nr0+ir-1,il,ih)
!                     ! r' > r --> F(r)S(r')
!                     work2(ir+1:nr) = fz(nr0+ir-1,il,ih)*sz(nr0+ir:nr,il,ih)
!                     work2(1:nr) = work2(1:nr)*phiref(1:nr,ib,il,is,ia)
!                     norm = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!                     work1(ir) = norm
!                   end do
!                   ! multiply by basis function and integrate over r'
!                   work1(1:nr) = work1(1:nr)*phiref(1:nr,jb,il,is,ia)
!                   norm = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!                   ! include the sqrt of the energy
!                   fsc(ib,jb,il,is,ia) = norm*ek
!                 end if
!                 ! -------------------------------------------------------------
!               end do
!             end do

        else ! ll' coupling 

            do jb=1,nb2              
              do ib=1,nb1
                i   = lmsb2i(ib,ilm,is,ia)
                j   = lmsb2i(jb,ilm2,is2,ia)         
                gfpq(i,j,ia,ie) = 0.d0        
              enddo ! ib
            enddo ! jb
 
        endif ! ll' coupling

      enddo ! ilmsn1
    enddo ! ilmsn2
  enddo ! iatom

  ! coefficients computed
  noirrcoeffs_soc(:,:,ie) = .false.

! Rotate the projected onsite Green function to the global frame
!  rotate_onsite  = .true.
!  rotate_regular = .false.
! Moment along z in the local frame
!  do ia = 1, nasusc
!    mom_loc(:,ia) = (/0.d0,0.d0,1.d0/)
!  end do 
!  call spinrot_scatt_sol(ie,mom_loc,magdir,rotate_onsite,rotate_regular)

! All done
  end subroutine irr_coeffs_soc
