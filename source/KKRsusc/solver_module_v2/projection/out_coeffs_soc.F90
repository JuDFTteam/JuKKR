  subroutine out_coeffs_soc(lmsize,ie,e,ek,de)
! Output projection coefficients to outsusc.dat file
  use global

  implicit none

! --> energy point label
  integer(kind=i4b), intent(in) :: ie
! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(in) :: e, ek, de
! -----------------------------------------------------------------------
  real(kind=r8b), parameter :: pi = 4.d0*atan(1.d0)
  integer(kind=i4b) :: ia, ih, il, il2, im, im2, ilm, ilm2, ilmsn, ilmsn2
  integer(kind=i4b) :: nb, ib, jb, i, j  ! basis
  integer(kind=i4b) :: is, is2 ! spin channels 
  integer(kind=i4b) :: lmsize
! -----------------------------------------------------------------------

  if (.not.lhdio) then
    esusc(ie) = e; eksusc(ie) = ek
    if (ikkr == 1 ) then  ! old impurity code
      desusc(ie) = -pi*de
    else                  ! KKRFLEX
      desusc(ie) = -pi*de/nsmax
    end if
    if (.not.lscfmesh) then
      escf(ie) = esusc(ie); descf(ie) = desusc(ie)
      if (ie == nesusc) then
        write(*,'(/," de sums to=",2f12.6,/)') sum(descf(1:nescf))
        efscf = real(escf(ie))  ! Fermi energy
      end if
    end if
!    Not used anymore --> Already done (see reg_coeff_soc.F90 and irr_coef_soc.F90)
!    call save_coeffs(ie,is,.false.)
    return
  end if

  write(iomain,'(" ie, e, ek, de=",i4)') ie
  write(iomain,'(6es16.8)') e, ek, de
  write(iomain,'(" Projection coefficients:")')

  do ia=1,nasusc
    ih = iasusc(ia)
    write(iomain,'(" ia=",i4,"  i1=",i4)') ia, ih
    write(iomain,*)'SOC --> Coefficients diagonal in l'
!   Loop over LL' and ss'
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
!       SOC does not couple l --> l' 
        if (il == il2) then    
          nb = iwsusc(il,is,ia)
          if (nb > 0) then
            write(iomain,'(6i4)') il, im, im2, is, is2, nb
            ! Right coeff
            write(iomain,*) " pz right coeff "                
            j = lms2i(ilm2,is2)
            do ib=1,nb
              i = lmsb2i(ib,ilm,is,ia) 
              write(iomain,'(100es16.8)')pzr(i,j,ia,ie)
            end do 
            ! Left coeff
            write(iomain,*) " pz left coeff "                
            do ib=1,nb
              i = lmsb2i(ib,ilm,is,ia) 
              write(iomain,'(100es16.8)')pzl(i,j,ia,ie)
            end do
!           ---------------------------------------------------------------
            if (isra == 1) then
              ! Right coeff
              write(iomain,*) " fz right coeff "                
              j = lms2i(ilm2,is2)
              do ib=1,nb
                i = lmsb2i(ib,ilm,is,ia) 
                write(iomain,'(100es16.8)')fzr(i,j,ia,ie)
              end do 
              ! Left coeff
              write(iomain,*) " fz left coeff "                
              do ib=1,nb
                i = lmsb2i(ib,ilm,is,ia) 
                write(iomain,'(100es16.8)')fzl(i,j,ia,ie)
              end do
            end if ! isra
!           ---------------------------------------------------------------
          endif ! nb
        end if ! ll' coupling
      end do ! ilmsn
    end do ! ilmsn2
  end do ! ia

  write(iomain,'(" Onsite coefficients:")')
  do ia=1,nasusc
    ih = iasusc(ia)
    write(iomain,'(" ia=",i4,"  i1=",i4)') ia, ih
    write(iomain,*)'SOC --> Coefficients diagonal in l'
!   Loop over LL' and ss'
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
        nb = iwsusc(il,is,ia)
!       SOC does not couple l --> l' 
        if (il == il2) then    
          if (nb > 0) then
            write(iomain,'(6i4)') il, im, im2, is, is2, nb
            write(iomain,*) " pq coeff"
            do ib= 1, nb
              do jb = 1, nb 
                i = lmsb2i(ib,ilm,is,ia)
                j = lmsb2i(jb,ilm2,is2,ia)           
                write(iomain,'(100es16.8)') gfpq(i,j,ia,ie)
              end do ! jb
            end do ! ib
!           ---------------------------------------------------------------
            if (isra == 1) then
              write(iomain,*) " ps coeff"
              do ib= 1, nb
                do jb = 1, nb 
                  i = lmsb2i(ib,ilm,is,ia)
                  j = lmsb2i(jb,ilm2,is2,ia)           
                  write(iomain,'(100es16.8)') gfps(i,j,ia,ie)
                end do ! jb
              end do ! ib
              write(iomain,*) " fq coeff"
              do ib= 1, nb
                do jb = 1, nb 
                  i = lmsb2i(ib,ilm,is,ia)
                  j = lmsb2i(jb,ilm2,is2,ia)           
                  write(iomain,'(100es16.8)') gffq(i,j,ia,ie)
                end do ! jb
              end do ! ib
              write(iomain,*) " fs coeff"
              do ib= 1, nb
                do jb = 1, nb 
                  i = lmsb2i(ib,ilm,is,ia)
                  j = lmsb2i(jb,ilm2,is2,ia)           
                  write(iomain,'(100es16.8)') gffs(i,j,ia,ie)
                end do ! jb
              end do ! ib
!             ---------------------------------------------------------------
            endif ! isra         
          end if ! nb
        end if ! ll' coupling 
      end do ! ilmsn
    end do ! ilmsn2
  end do ! ia

! All done
  end subroutine out_coeffs_soc
