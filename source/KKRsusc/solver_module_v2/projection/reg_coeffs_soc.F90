  subroutine reg_coeffs_soc(natomd,lmsize,nrmaxd,pz,fz,ie,ek,rightleft,tmat)
! Energy dependent projection coefficients of regular scattering wfns
! Expects solutions of the ASA SRA problem, no m-dependence of pz or fz
! BEWARE OF THE SPECIAL TREATMENT OF THE RADIAL MESH, IR=1 NOT USED
! CAREFUL --> lmsize = 2*(lmaxd+1)**2
  use global
  use bessel_new 

  implicit none

! --> dimensions of the wavefunction arrays
  integer(kind=i4b), intent(in) :: natomd, lmsize, nrmaxd
! --> regular scattering solutions
  complex(kind=c8b), intent(in) :: pz(nrmaxd,lmsize,lmsize,natomd), fz(nrmaxd,lmsize,lmsize,natomd)
! --> current energy point
  integer(kind=i4b), intent(in) :: ie
! --> defined value of the square-root of the energy
  complex(kind=c8b), intent(in) :: ek
! --> tmatrix-in the projection basis
  complex(kind=c8b), intent(inout) :: tmat(lmsize,lmsize,natomd) 
! ------------------------------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nr, nr0, nr1, il, nb, ib, ir
  integer(kind=i4b) :: ilmsn, ilmsn2, ilm, ilm2, il2, im, im2, i, j
  integer(kind=i4b) :: is, is2  ! Spin channels 
  real(kind=r8b)    :: r(nrmax), dr(nrmax), mom_loc(3,nasusc)
  complex(kind=c8b) :: work1(nrmax), work2(nrmax), norm1, norm2
  logical           :: rightleft, rotate_onsite, rotate_regular
  complex(kind=c8b) :: rjl(nrmax,0:nlmax), tmatl(0:nlmax), vpot(2,2,natomd)
  complex(kind=c8b), external :: radint
! ------------------------------------------------------------------------------------

  ! Atoms loop  
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
    enddo ! il 
    do ilmsn2 = 1, lmsize 
      is2  = i2lms_new(2,ilmsn2)
      ilm2 = i2lms_new(1,ilmsn2)  
      il2  = i2lm(2,ilm2)
      im2  = i2lm(1,ilm2)
      do ilmsn = 1, lmsize
        is  = i2lms_new(2,ilmsn)
        ilm = i2lms_new(1,ilmsn)
        il  = i2lm(2,ilm)
        im  = i2lm(1,ilm)
        nb  = iwsusc(il,is,ia)
        tmatl(il) = 0.d0
        if (nobasis(il,is,ia) .and. nb > 0) then
          write(*,'("reg_coeffs: no basis",6i4)') il, im, is, im2, is2, ia
          stop
        end if
!       SOC ---> Radial solution Rlm,m'(r) and Phil(r) basis
!       SOC does not couple l --> l' 
        if (il == il2) then
!         projection on each basis function
!         basis functions assumed normalized to 1 with weight dr
          do ib=1,nb          
            work1(1:nr) = pz(nr0:nr1,ilmsn,ilmsn2,ih)*phiref(1:nr,ib,il,is,ia)
            norm1 = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
            pzc(ib,il,is,ia) = norm1
!           ----------------------------------------------------------------
            if (isra == 1) then
              work2(1:nr) = fz(nr0:nr1,ilmsn,ilmsn2,ih)*phiref(1:nr,ib,il,is,ia)
              norm2 = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
              fzc(ib,il,is,ia) = norm2
            end if
!           ----------------------------------------------------------------
!           t-matrix in projection basis
!            work1(1:nr)  = rjl(1:nr,il)*(vpot(is,is2,ih) - 2.d0*zat(ia)/r(1:nr))*phiref(1:nr,ib,il,is,ia)
!            tmatl(il) = tmatl(il) + pzc(ib,il,is,ia)*radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
          end do 
!         check how much of the function is left
          if (nb > 0) then
            work1(1:nr) = pz(nr0:nr1,ilmsn,ilmsn2,ih)
            do ib=1,nb
              norm1 = pzc(ib,il,is,ia)
              work1(1:nr) = work1(1:nr) - norm1*phiref(1:nr,ib,il,is,ia)
            end do
!           To correct
            work1(1:nr) = work1(1:nr)*conjg(work1(1:nr))
            norm1 = radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
            work1(1:nr) = conjg(pz(nr0:nr1,ilmsn,ilmsn2,ih))*pz(nr0:nr1,ilmsn,ilmsn2,ih)
            norm1 = norm1/radint(nr,work1(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
!           ----------------------------------------------------------------
            norm2 = 0.d0
            if (isra == 1) then
              work2(1:nr) = fz(nr0:nr1,ilmsn,ilmsn2,ih)
              do ib=1,nb
                norm2 = fzc(ib,il,is,ia)
                work2(1:nr) = work2(1:nr) - norm2*phiref(1:nr,ib,il,is,ia)
              end do
              work2 = work2*conjg(work2)
              norm2 = radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
              work2(1:nr) = conjg(fz(nr0:nr1,ilmsn,ilmsn2,ih))*fz(nr0:nr1,ilmsn,ilmsn2,ih)
              norm2 = norm2/radint(nr,work2(1:nr),dr(1:nr),npanat(ia),ircutat(:,ia))
            end if           
!           ----------------------------------------------------------------
!            if (lhdio) write(iodb,'("reg_coeffs",7i4,2es16.3)') ie, ia, il, im, is, im2, is2, sqrt(abs(norm1)), sqrt(abs(norm2))
          end if ! nb
!         Store the coefficients  
          j = lms2i(ilm2,is2)
          if (rightleft) then
            do ib = 1, nb
              i = lmsb2i(ib,ilm,is,ia)
              pzr(i,j,ia,ie) = pzc(ib,il,is,ia) 
              if (isra == 1) fzr(i,j,ia,ie) = fzc(ib,il,is,ia)
            end do
          else ! Left solutions 
            do ib = 1, nb
              i = lmsb2i(ib,ilm,is,ia)
              pzl(i,j,ia,ie) = pzc(ib,il,is,ia)
              if (isra == 1) fzl(i,j,ia,ie) = fzc(ib,il,is,ia)
            end do   
          endif ! right left

        else ! ll' coupling 

          j = lms2i(ilm2,is2)
          if (rightleft) then
            do ib = 1, nb       
              i = lmsb2i(ib,ilm,is,ia)
              pzr(i,j,ia,ie) = 0.d0 
              if (isra == 1) fzr(i,j,ia,ie) = 0.d0
            end do
          else ! Left solutions 
            do ib = 1, nb       
              i = lmsb2i(ib,ilm,is,ia)
              pzl(i,j,ia,ie) = 0.d0 
              if (isra == 1) fzl(i,j,ia,ie) = 0.d0
            end do
          endif ! right left

        endif ! ll' coupling    

!       -------------------------------------------------------------
!       t-matrix in the projection basis ?
        if (ltmatproj) tmat(ilmsn,ilmsn2,ia) = tmatl(il)
!       ------------------------------------------------------------- 

      end do ! ilmsn1
    end do ! ilmsn2
  end do ! iatom
  ! coefficients computed   
  noregcoeffs_soc (:,:,ie) = .false.

! Rotate the projected regular solutions to the global frame
!  rotate_onsite  = .false. 
!  rotate_regular = .true. 
! Moment along z in the local frame
!  do ia = 1, nasusc 
!    mom_loc(:,ia) = (/0.d0,0.d0,1.d0/)
!  end do
!  call spinrot_scatt_sol(ie,mom_loc,magdir,rotate_onsite,rotate_regular)

! All done
  end subroutine reg_coeffs_soc
