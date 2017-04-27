  subroutine in_coeffs_soc(lmsize,ie,e,ek,de)
! Output projection coefficients to outsusc.dat file
  use global

  implicit none

! --> energy point label
  integer(kind=i4b), intent(in) :: ie, lmsize
! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(out):: e, ek, de
! -----------------------------------------------------------------------
  real(kind=r8b),   parameter :: pi  = 4.d0*atan(1.d0)
  real(kind=r8b),   parameter :: tol = 1.d-6
  complex(kind=c8b)             :: er, ekr, der
  integer(kind=i4b) :: ia, ih, il, il2, im, im2, ilm, ilm2, ilmsn, ilmsn2
  integer(kind=i4b) :: nb, ib, jb, i, j  ! basis
  integer(kind=i4b) :: is, is2 ! spin channels   
  real(kind=r8b)    :: x(1000)
  integer(kind=i4b) :: n(1000)
  character*60      :: header
! -----------------------------------------------------------------------

  read(iomain,'(a)') header ! ie, e, ek, de
  read(iomain,'(6es16.8)') x(1:6)
  ! Read energies and energy weights
  e  = cmplx(x(1),x(2))
  ek = cmplx(x(3),x(4))
  de = cmplx(x(5),x(6))     

  ! Test if ordering of energies is correct
!  if (abs(real(er) - real(e)) > tol) then 
!   write(*,*) abs(real(er) - real(e)) 
!   stop 'in_coeffs_soc --> Problem in the energy parallelization'
!  end if 
!  if (abs(aimag(er) - aimag(e)) > tol) then 
!    write(*,*) abs(aimag(er) - aimag(e))
!    stop 'in_coeffs_soc --> Problem in the energy parallelization'    
!  end if 

  read(iomain,'(a)') header ! proj coeffs
  do ia=1,nasusc
    ih = iasusc(ia)
    read(iomain,'(a)') header ! ia, i1   
    read(iomain,'(a)') header ! 'SOC --> Coefficients diagonal in l'
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
            read(iomain,*) n(1:6)!il, im, im2, is, is2, nb
            ! Right coeff
            read(iomain,'(a)') header ! pz right coeff
            j = lms2i(ilm2,is2)
            do ib=1,nb
              i = lmsb2i(ib,ilm,is,ia)
              read(iomain,*) x(1:2)  
              pzr(i,j,ia,ie) = cmplx(x(1),x(2))
            end do 
            ! Left coeff
            read(iomain,'(a)') header ! pz left coeff 
            do ib=1,nb
              i = lmsb2i(ib,ilm,is,ia) 
              read(iomain,*) x(1:2) 
              pzl(i,j,ia,ie) = cmplx(x(1),x(2))
            end do
!           ---------------------------------------------------------------
            if (isra == 1) then
              ! Right coeff
              read(iomain,'(a)') header !" fz right coeff "
              j = lms2i(ilm2,is2)
              do ib=1,nb
                i = lmsb2i(ib,ilm,is,ia)
                read(iomain,*) x(1:2) 
                fzr(i,j,ia,ie) = cmplx(x(1),x(2)) 
              end do 
              ! Left coeff
              read(iomain,'(a)') header !" fz left coeff "                
              do ib=1,nb
                i = lmsb2i(ib,ilm,is,ia) 
                read(iomain,*) x(1:2)
                fzl(i,j,ia,ie) = cmplx(x(1),x(2)) 
              end do
            end if ! isra
!           ---------------------------------------------------------------
          endif ! nb
        end if ! ll' coupling
      end do ! ilmsn
    end do ! ilmsn2
  end do ! ia

  read(iomain,'(a)') header ! Onsite coefficients 
  do ia=1,nasusc
    ih = iasusc(ia)
    read(iomain,'(a)') header ! ia, ih
    read(iomain,'(a)') header ! SOC --> Coefficients diagonal in l
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
            read(iomain,*) n(1:6) ! il, im, im2, is, is2, nb
            read(iomain,'(a)') header !" pq coeff"
            do ib= 1, nb
              do jb = 1, nb 
                i = lmsb2i(ib,ilm,is,ia)
                j = lmsb2i(jb,ilm2,is2,ia)        
                read(iomain,*) x(1:2) 
                gfpq(i,j,ia,ie) = cmplx(x(1),x(2))  
              end do ! jb
            end do ! ib
!           ---------------------------------------------------------------
            if (isra == 1) then
              read(iomain,'(a)') header !" ps coeff"
              do ib= 1, nb
                do jb = 1, nb 
                  i = lmsb2i(ib,ilm,is,ia)
                  j = lmsb2i(jb,ilm2,is2,ia)           
                  read(iomain,*) x(1:2)
                  gfps(i,j,ia,ie) = cmplx(x(1),x(2))
                end do ! jb
              end do ! ib
              read(iomain,'(a)') header !" fq coeff"
              do ib= 1, nb
                do jb = 1, nb 
                  i = lmsb2i(ib,ilm,is,ia)
                  j = lmsb2i(jb,ilm2,is2,ia)     
                  read(iomain,*) x(1:2)         
                  gffq(i,j,ia,ie) = cmplx(x(1),x(2))
                end do ! jb
              end do ! ib
              read(iomain,'(a)') header ! "fs coeff"
              do ib= 1, nb
                do jb = 1, nb 
                  i = lmsb2i(ib,ilm,is,ia)
                  j = lmsb2i(jb,ilm2,is2,ia)      
                  read(iomain,*) x(1:2)     
                  gffs(i,j,ia,ie) = cmplx(x(1),x(2))
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
  end subroutine in_coeffs_soc
