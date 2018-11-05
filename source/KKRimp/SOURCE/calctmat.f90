module mod_calctmat
!-------------------------------------------------------------------------------
!> Summary: Calculate the t-matrices for the actual system 
!> Author: Phivos Mavropoulos, Hubert Ebert, Voicu Popescu
!> Category: KKRimp, single-site 
!>           
!-------------------------------------------------------------------------------

contains
!-------------------------------------------------------------------------------
!> Summary: Calculate the t-matrices for the actual system 
!> Author: Phivos Mavropoulos, Hubert Ebert, Voicu Popescu
!> Category: KKRimp, single-site 
!>           
!-------------------------------------------------------------------------------

SUBROUTINE CALCTMAT(ERYD,VPOT,CELL,ZATOM,LMAXATOM,TMATLL,config,ispin,nspin)!     &
! C!           C                IDOLDAU,LOPT,WLDAU)

! 
! 
! C
! C *********************************************************************
! C * For KREL = 1 (relativistic mode)                                  *
! C *                                                                   *
! C *  NPOTD = 2 * NATYPD                                               *
! C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! C *  NSPIND = 1                                                       *
! C *                                                                   *
! C *  LDA+U implementation     Mar. 2002-Dec.2004                      *
! C *                           ph.mavropoulos, h. ebert, v. popescu    *
! C * Notes:                                                            *
! C *  average WLDAU for spherical wavefunctions:                       *
! C *  The spherical part of the d or f wavefunction is found by adding *
! C *  the average interaction potential WLDAUAV to the spherical       *
! C *  potential. Then the non-spherical parts are found by using only  *
! C *  the deviation of WLDAU from the average. This speeds up the      *
! C *  convergence of the Born series. See also subroutines             *
! C *  regsol, pnstmat and pnsqns                                       *
! C *                                                                   *
! C *********************************************************************
  use nrtype
  use mod_physic_params, only: cvlight
!   use configparams, only: icst,ins,nsra
  use mod_gauntharmonics, only: gauntcoeff
  use mod_pnstmat
  use mod_wfmesh
  use mod_cradwf
  use type_cell
  use type_config
  use mod_config, only: config_testflag

  implicit none
!interface variables
  complex(kind=dpc),intent(in)              ::  eryd
  type(cell_type),intent(in)                ::  cell
  real(kind=dp),intent(in)                  ::  vpot(1:cell%nrmaxd,(2*lmaxatom+1)**2)
  real(kind=dp),intent(in)                  ::  zatom
  integer,intent(in)                        ::  lmaxatom
  complex(kind=dpc)                         ::  tmatll((lmaxatom+1)**2,(lmaxatom+1)**2)
  type(config_type),intent(in)              ::  config
  integer,intent(in)                        ::  ispin
  integer,intent(in)                        ::  nspin

!local variables
  real(kind=dp),allocatable                 ::  vpottemp(:,:)
  integer                                   ::  lm1,lm2
  complex(kind=dpc)                         ::  ek
  integer                                   ::  lmmaxatom,lmpotatom
  logical, parameter                        ::  lirrsol=.true.
  complex(kind=dpc),allocatable             ::  alpha(:), &
                                                fz(:,:), & 
                                                pns(:,:,:,:), &
                                                pz(:,:), &
                                                qz(:,:), & 
                                                sz(:,:), &
                                                tmat(:)
  real(kind=dp),allocatable                 ::  rs(:,:), &
                                                s(:)
!   complex(kind=dpc)                         ::  alpha(0:lmaxatom), &
!                                                 fz(cell%nrmax,0:lmaxatom), & 
!                                                 pns((lmaxatom+1)**2,(lmaxatom+1)**2,cell%nrmin_ns:cell%nrmax,2), &
!                                                 pz(cell%nrmax,0:lmaxatom), &
!                                                 qz(cell%nrmax,0:lmaxatom), & 
!                                                 sz(cell%nrmax,0:lmaxatom), &
!                                                 tmat(0:lmaxatom)
!   real(kind=dp)                             ::  rs(cell%nrmax,0:lmaxatom), &
!                                                 s(0:lmaxatom)



  integer                                   ::  idoldau,lmlo,lmhi,mmax,m1,lopt
  real(kind=dp)                             ::  wldauav
  real(kind=dp),allocatable                 ::  cutoff(:),wldau(:,:,:)

lmmaxatom = (lmaxatom+1)**2
lmpotatom = (2*lmaxatom+1)**2

  allocate (                                    alpha(0:lmaxatom), &
                                                fz(cell%nrmax,0:lmaxatom), & 
                                                pns((lmaxatom+1)**2,(lmaxatom+1)**2,cell%nrmin_ns:cell%nrmax,2), &
                                                pz(cell%nrmax,0:lmaxatom), &
                                                qz(cell%nrmax,0:lmaxatom), & 
                                                sz(cell%nrmax,0:lmaxatom), &
                                                tmat(0:lmaxatom) )
  allocate (                                    rs(cell%nrmax,0:lmaxatom), &
                                                s(0:lmaxatom) )



! write(*,*) 'sdf'


! write(*,*) 'tmatll alloc?',allocated(tmatll)
! stop
! allocate( tmatll((lmaxatom+1)**2, (lmaxatom+1)**2) )
allocate(cutoff(cell%nrmax), wldau(2*lmaxatom+1,2*lmaxatom+1,nspin))

! allocate(vpottemp(cell%nrmax-cell%nrmin_ns,(2*lmaxatom+1)**2))
allocate(vpottemp(cell%nrmin_ns:cell%nrmax,(2*lmaxatom+1)**2))

tmatll=(0.0_dp,0.0_dp)

!-------------------------------------------------------
!-- jet to be implemented LDA+U
!-------------------------------------------------------
idoldau=0
if ( idoldau.eq.1 ) then
  wldauav = 0.d0                                        
  lmlo = lopt*lopt + 1
  lmhi = (lopt+1)*(lopt+1)
  mmax = lmhi - lmlo + 1
  do m1 = 1,mmax                                        
    wldauav = wldauav + wldau(m1,m1,ispin)         
  enddo                                                 
  wldauav = wldauav/dble(mmax)                        
  do m1 = 1,cell%nrmax
    cutoff(m1) = 1.d0
  end do
end if                                                    


call wfmesh(eryd,ek,cvlight,config%nsra,zatom,cell%rmesh,s,rs, &
            cell%nrmax,lmaxatom)

call cradwf(eryd,ek,config%nsra,alpha,cell%npan,cell%nrcut,cvlight,rs,s,              &
            pz,fz,qz,sz,tmat,vpot(:,1),cell%drmeshdi,cell%rmesh,zatom,lirrsol, &
            idoldau,lopt,wldauav,cutoff,lmaxatom,lmaxatom+1,cell%nrmax) 


if ( config_testflag('tmatdebug') ) then
  write(1999,'(50000F)') cell%rmesh
!   write(5001,'(50000F)') cell%rmesh
  do lm1=0,lmaxatom
      write(2000,'(50000E)') pz(:,lm1)
      write(2001,'(50000E)') fz(:,lm1)
      write(2002,'(50000E)') qz(:,lm1)
      write(2003,'(50000E)') sz(:,lm1)
  end do
!   stop
end if




if ( config%ins.eq.0 ) then
  do lm1 = 1,lmmaxatom
    tmatll(lm1,lm1) = tmat(gauntcoeff(lmaxatom)%loflm(lm1))
  end do
else
  vpottemp=vpot(cell%nrmin_ns:cell%nrmax,:)
  call pnstmat(cell%drmeshdi,ek,config%icst,pz,qz,fz,sz,pns,tmatll,vpottemp,& !vpot(cell%nrmin_ns:cell%nrmax,:), &
               cell%npan, cell%nrcut,config%nsra, &
               gauntcoeff(lmaxatom)%cleb,gauntcoeff(lmaxatom)%icleb,gauntcoeff(lmaxatom)%iend,gauntcoeff(lmaxatom)%loflm, &
               gauntcoeff(lmaxatom)%ncleb, &
               tmat,lmaxatom, &
               idoldau,lopt,lmlo,lmhi,wldau,wldauav,cutoff, &
               cell%nrmax, cell%nrmin_ns,lmaxatom,2*lmaxatom+1,lmmaxatom,lmpotatom) 
end if

if ( config_testflag('tmatdebug') ) then
!   write(1999,'(50000F)') cell%rmesh
  write(5003,'(5000F)')  cell%rmesh(cell%nrmin_ns:cell%nrmax)
!   write(5001,'(50000F)') cell%rmesh
  do lm1=1,(lmaxatom+1)**2
    do lm2=1,(lmaxatom+1)**2
      write(2100,'(50000E)') pns(lm2,lm1,:,1)
      write(2101,'(50000E)') pns(lm2,lm1,:,2)
    end do
  end do
!   stop
end if



! do lm1=1,(lmaxatom+1)**2
!   write(3001,'(5000F)') tmatll(:,lm1)
! end do
! stop

end subroutine calctmat


end module mod_calctmat
