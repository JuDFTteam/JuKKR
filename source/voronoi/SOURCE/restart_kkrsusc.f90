  subroutine restart_kkrsusc(itc,lmaxd,config,ielast,ez,wz,my_rank)

! --> modules
  use global
  use type_config
  use mod_config
  use mod_physic_params, only: cvlight
! ---------------------------------------------------------------------------------

  implicit none

! --> energy point label, mpi_stuff
  integer(kind=i4b), intent(in) :: itc, lmaxd, ielast
! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(in) :: ez(ielast), wz(ielast)
! ---------------------------------------------------------------------------------
  type(config_type)            :: config
  integer(kind=i4b), parameter :: nloutsusc =  40
  integer(kind=i4b) :: ia, ih, il, il2, im, im2, ilm, ilm2, ilmsn, ilmsn2, ie
  integer(kind=i4b) :: nb, ib, jb, i, j, k, l  ! basis
  integer(kind=i4b) :: lmsize, my_rank
  complex(kind=c8b) :: ek_susc
  logical           :: exists
! ---------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------
! Check if outsusc.dat files are written 
! ---------------------------------------------------------------------------------
! Projection coeff and onsite GF
  inquire(file='outsusc.dat',exist=exists)
  if (exists) then
    write(*,'(" projection coeff, gfpq, gmat and tmat read from outsusc.dat")')
  else 
    stop 'to run in restart mode provide outsusc.dat first'
  end if
! Structural GF
  inquire(file='outsusc_green.dat',exist=exists)
  if (exists) then
    write(*,'(" structural green function, gmat is read from outsusc_green.dat")')
  else 
    stop 'to run in restart mode provide outsusc_green.dat first'
  end if
! T-matrix
   inquire(file='outsusc_tmat.dat',exist=exists)
  if (exists) then
    write(*,'(" t-matrix, gmat is read from outsusc_tmat.dat")')
  else 
    stop 'to run in restart mode provide outsusc_tmat.dat first'
  end if
! ---------------------------------------------------------------------------------
! Read kkrsusc options
!  open(file='outsusc.dat',unit=iomain,status='old')
  call outsusc_in_soc(itc,my_rank)
! ---------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------
!  lmsize = 2*lmmax
! Energy loop 
  
! First read projection coeff and then onsite gf
!  do ie = 1, ielast
!    if (config%nsra == 1) ek_susc = sqrt(ez(ie))                                          
!    if (config%nsra == 2) ek_susc = sqrt(ez(ie) + (ez(ie)/cvlight)*(ez(ie)/cvlight))      
!    call in_coeffs_soc(lmsize,ie,ez(ie),ek_susc,wz(ie)) 
!  end do   
! ---------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------
! Second read gmat and tmat
!  do ie = 1, ielast
!    ! Read in gmat
!    call in_gmat_soc(ie)
!    ! Read in tmat  
!    call in_tmat_soc(ie)
!  end do
! ---------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------
! The basis function are recalculated (Fast)
! ---------------------------------------------------------------------------------
 
! All done
  end subroutine restart_kkrsusc
