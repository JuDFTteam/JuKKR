!------------------------------------------------------------------------------------
!> Summary: Module handling the relativistic exchange interactions
!> Author:
!> For details see Ebert and Mankovsky, PRB 79, 045209 (2009) 
!------------------------------------------------------------------------------------
module mod_calccouplingconstants

contains 

!-------------------------------------------------------------------------------
!> Summary: Relativistic exchange interactions
!> Author:
!> Category: physical-observables, input-output, KKRimp
!> Deprecated: False 
!> For details see Ebert and Mankovsky, PRB 79, 045209 (2009) 
!-------------------------------------------------------------------------------
subroutine calcJijmatrix(gmat,tmat,natom,wez,Jijmatrix,Aimatrix)
use type_tmat, only: tmat_type
use type_gmat, only: gmat_type
use mod_mathtools, only: matmat
use mod_config, only: config_testflag
use mod_types, only: t_inc

implicit none
type(gmat_type)         :: gmat
type(tmat_type)         :: tmat(natom,1)
integer                 :: natom
double complex          :: wez
double precision        :: Jijmatrix(3,3,natom,natom)
double precision        :: Aimatrix(3,natom)

!local
integer :: iatom,jatom,kalpha,lalpha
integer :: istart,istop,jstart,jstop
integer :: ilmsize,jlmsize
integer :: ilm1 
double precision :: Jscal

double complex,allocatable :: Gij(:,:),Gji(:,:)
double complex,allocatable :: Atemp(:,:),Btemp(:,:),Ctemp(:,:)
double complex,allocatable :: deltaTik(:,:),deltaTjl(:,:)
double complex             :: traceC
double precision           :: Jijmatrixtemp(3,3)
double complex             :: Jijmatrixtemp_complex(3,3)
double precision           :: Aimatrixtemp(3)
double complex             :: Aimatrixtemp_complex(3)

! double precision           :: pi
if (t_inc%i_write>0) write(1337,*) '##########################################'
if (t_inc%i_write>0) write(1337,*) '  starting calculation of the Jij-matrix'
if (t_inc%i_write>0) write(1337,*) '##########################################'

if (.not. allocated (gmat%gmat) ) then
  stop '[calcJijmatrix] gmat not allocated'
end if
! pi=3.1415926535897931D0
! print *,pi


Jscal=1.0D0/2.0D0
! The factor -2/pi is already contained in wez, where the 2 is due to the double
! occupancy with with spin up/down in a non-magnetic calculation. Since we do 
! magnetic calculations, we need to correct for that.
! factor 1/2.0D0 -> wez correction
! factor of 4.0D0 -> check Ebert!



do iatom=1,natom

  do jatom=1,natom
 
    istart=gmat%iatom2nlmindex(1,iatom)
    istop =gmat%iatom2nlmindex(3,iatom)
    jstart=gmat%iatom2nlmindex(1,jatom)
    jstop =gmat%iatom2nlmindex(3,jatom)
  
    ilmsize=istop-istart+1
    jlmsize=jstop-jstart+1
  
  !         write(*,*) 'istart',istart
  !         write(*,*) 'istop',istop
  !         write(*,*) 'istart',istart
  !         write(*,*) 'jstop',jstop
  !         write(*,*) 'ilmsize',ilmsize
  !         write(*,*) 'jlmsize',jlmsize
  
    allocate( Gij(ilmsize,jlmsize) )
    allocate( Gji(jlmsize,ilmsize) )
    allocate( deltaTik(ilmsize,ilmsize) )
    allocate( deltaTjl(jlmsize,jlmsize) )
  
    allocate( Atemp(ilmsize,jlmsize))
    allocate( Btemp(jlmsize,ilmsize))
    allocate( Ctemp(ilmsize,ilmsize))
  
    Gij = gmat%Gmat(istart:istop, jstart:jstop )
    Gji = gmat%Gmat(jstart:jstop, istart:istop )
  
    do kalpha=1,3
      deltaTik = tmat(iatom,1)%deltaT_Jij(:,:,kalpha)

      if (iatom==jatom) then
        Atemp = matmat(deltaTik,Gij)
        traceC=(0.0D0,0.0D0)
        do ilm1=1,ilmsize
          traceC=traceC+Atemp(ilm1,ilm1)
        end do
        Aimatrixtemp_complex(kalpha)=      traceC*Jscal
        Aimatrixtemp        (kalpha)=dimag(traceC*Jscal)
      end if
    
      do lalpha=1,3
      

!         write(92921,*) Gij
!         write(92922,*) Gji
!               write(*,*) 'deltaT_Jij(:,:,kalpha)',ubound(tmat(iatom,1)%deltaT_Jij)
! stop
        deltaTjl = tmat(jatom,1)%deltaT_Jij(:,:,lalpha)

!         write(92923,*) deltaTik
!         write(92924,*) deltaTjl

      
        Atemp = matmat(deltaTik,Gij)
!         write(92925,*) Atemp

        Btemp = matmat(deltaTjl,Gji)
!         write(92926,*) Btemp
        Ctemp = matmat(Atemp,Btemp)
!         write(92927,*) Ctemp
      
        traceC=(0.0D0,0.0D0)
        do ilm1=1,ilmsize
          traceC=traceC+Ctemp(ilm1,ilm1)
        end do
      
        Jijmatrixtemp_complex(kalpha,lalpha)=      traceC*Jscal
        Jijmatrixtemp        (kalpha,lalpha)=dimag(traceC*Jscal)

      
      end do !lalpha
    
    end do !kalpha

    deallocate( deltaTik, deltaTjl )
    deallocate( Gij,Gji )
    deallocate( Atemp,Btemp,Ctemp)

    if (iatom==jatom) then
      Aimatrix(:,iatom)=Aimatrix(:,iatom)+dimag(wez*Aimatrixtemp_complex)
    end if
    Jijmatrix(:,:,iatom,jatom)=Jijmatrix(:,:,iatom,jatom)+dimag(wez*Jijmatrixtemp_complex)

    if ( config_testflag('Jij(E)') ) then
      write(234932875,'(2I5,500E25.14)') iatom, jatom, wez, dimag(wez*Jijmatrixtemp_complex)
      write(234932876,'(2I5,500E25.14)') iatom, jatom,Jijmatrixtemp
    end if 
  end do !jatom

end do !iatom

end subroutine calcJijmatrix


!-------------------------------------------------------------------------------
!> Summary: Matrix elements of Bxc for the exchange interactions
!> Author:
!> Category: physical-observables, KKRimp
!> Deprecated: False 
!> The matrix elements of Bxc with Pauli matrices between two regular scattering wavefunctions are computed
!> This is then transformed from the local frame to the global frame
!> For details see Ebert and Mankovsky, PRB 79, 045209 (2009) 
!-------------------------------------------------------------------------------
subroutine calccouplingdeltat(wavefunction,deltaTmat,cellnew,gauntcoeff,theta,phi,lmmax,lmsize,lmax,lmpot,nrmax)
use nrtype, only: dp
use type_wavefunction, only: wavefunction_type 
use type_cellnew, only: cell_typenew
use type_gauntcoeff, only: gauntcoeff_type
use mod_mathtools, only: matmat, matmat1T,matmatT1
use mod_vllmat, only: vllmat
use mod_intcheb_cell, only: intcheb_cell
use mod_rotatespinframe, only: rotatematrix

implicit none
type(wavefunction_type)                   :: wavefunction
type(cell_typenew)                        :: cellnew
type(gauntcoeff_type)                     :: gauntcoeff
integer                                   :: lmmax
integer                                   :: lmsize
integer                                   :: lmax
integer                                   :: lmpot
integer                                   :: nrmax
!local
double precision                          :: theta,phi
integer                                   :: leftsol
double precision                          :: bpot(nrmax,lmpot,1)
double complex                          :: Bpotll(lmmax,lmmax,nrmax)
double complex                            :: lambda(2,2,3)
integer                                   :: ir,ilm1,ilm2,kspin,mspin,lspin,mspinshift,lspinshift,lm1,lm2
double complex                            :: bpot2(lmsize,lmsize),rtemp(lmsize,lmsize),atemp(lmsize,lmsize)
double complex                            :: RBpotR(lmsize,lmsize,nrmax),RBpotR_integrated(lmsize,lmsize,3),fntemp(nrmax)
double complex                            :: deltaTmat(lmsize,lmsize,3)

if (allocated(wavefunction%rllleft) ) then
  leftsol=1
else
  leftsol=0
end if

bpot(:,:,1) = (cellnew%vpotnew(:,:,1)-cellnew%vpotnew(:,:,2)) *0.5D0 ! not quite sure if this factor should appear here
                                                              !        check the reference of Ebert

call vllmat(1, cellnew%nrmaxnew, cellnew%nrmaxnew, lmmax, lmmax, bpotll, bpot, lmpot, gauntcoeff%cleb, gauntcoeff%icleb, gauntcoeff%iend, 2, 0.0_dp, cellnew%rmeshnew, 0, gauntcoeff%ncleb)

call calclambda(lambda,theta,phi)


! ##########################################################################################3
! calculate  deltaT_alpha = int ( R^T(r) (U sigma_alpha U^T -sigma_z) * B(r) R(r) )
! where all functions are matricies in L
! ##########################################################################################3
do kspin=1,3

  do ir=1,nrmax
  
    do lspin=1,2
      lspinshift=lmmax*(lspin-1)
      do mspin=1,2
        mspinshift=lmmax*(mspin-1)
        do lm1=1,lmmax
          do lm2=1,lmmax
            bpot2(mspinshift+lm2,lspinshift+lm1)=bpotll(lm2,lm1,ir)*lambda(mspin,lspin,kspin)
          end do 
        end do 
      end do !mspin
    end do !lspin

    rtemp = wavefunction%rll(:lmsize,:lmsize,ir,1)
    atemp= matmat(bpot2,rtemp)
!     leftsol=1
    if (leftsol==1) then
      rtemp = wavefunction%rllleft(:lmsize,:lmsize,ir,1)
    else
      rtemp = wavefunction%rll(:lmsize,:lmsize,ir,1)
    end if
  
    rbpotr(:,:,ir) = matmatT1(rtemp,atemp)
  
  end do !nrmax
  
  
  do ilm1=1,lmsize
    do ilm2=1,lmsize
!       write(*,*) ilm1,ilm2,lmsize
      fntemp = RBpotR(ilm1,ilm2,:)
      call intcheb_cell(fntemp, RBpotR_integrated(ilm1,ilm2,kspin), cellnew%rpan_intervall, cellnew%ipan_intervall, cellnew%npan_tot, cellnew%ncheb, cellnew%nrmaxnew)
    end do
  end do
  call rotatematrix(RBpotR_integrated(:,:,kspin), theta, phi, lmmax,0 ) ! 'loc->glob')
! stop
end do !kspin

deltaTmat=RBpotR_integrated

! stop
end subroutine calccouplingdeltat


!-------------------------------------------------------------------------------
!> Summary: Pauli matrices transformed from the global to the local frame
!> Author:
!> Category: special-functions, physical-observables, KKRimp
!> Deprecated: False 
!> 
!-------------------------------------------------------------------------------
subroutine calclambda(lambda,theta,phi)
use mod_rotatespinframe, only: rotatematrix
implicit none
double complex :: lambda(2,2,3)
double precision :: theta, phi
double complex :: sigmatemp(2,2,3),sigma(2,2,3)
integer        :: ispin
call calc_sigma(sigma)
sigmatemp=sigma
do ispin=1,3
  call rotatematrix(sigmatemp(:,:,ispin), theta, phi, 1, 1) !'glob->loc')
!   lambda(:,:,ispin) = sigmatemp(:,:,ispin)- sigma(:,:,3)
  lambda(:,:,ispin) = sigmatemp(:,:,ispin) !test
end do !ispin
end subroutine calclambda


!-------------------------------------------------------------------------------
!> Summary: Initialization of Pauli matrices
!> Author:
!> Category: initialization, special-functions, physical-observables, KKRimp
!> Deprecated: False 
!> The assignment of the spin indices is provided in two ways
!-------------------------------------------------------------------------------
subroutine calc_sigma(sigma)
implicit none
double complex :: sigma(2,2,3)
integer        :: verbose
character(len=*),parameter :: conventionmode='kkr'
verbose=0

if (conventionmode=='normal') then
  sigma(1,1,1)=( 0.0D0, 0.0D0)
  sigma(1,2,1)=( 1.0D0, 0.0D0)
  sigma(2,1,1)=( 1.0D0, 0.0D0)
  sigma(2,2,1)=( 0.0D0, 0.0D0)
  
  sigma(1,1,2)=( 0.0D0, 0.0D0)
  sigma(1,2,2)=( 0.0D0,-1.0D0)
  sigma(2,1,2)=( 0.0D0, 1.0D0)
  sigma(2,2,2)=( 0.0D0, 0.0D0)
  
  sigma(1,1,3)=( 1.0D0, 0.0D0)
  sigma(1,2,3)=( 0.0D0, 0.0D0)
  sigma(2,1,3)=( 0.0D0, 0.0D0)
  sigma(2,2,3)=(-1.0D0, 0.0D0)
elseif (conventionmode=='kkr') then
  sigma(1,1,1)=( 0.0D0, 0.0D0)
  sigma(1,2,1)=( 1.0D0, 0.0D0)
  sigma(2,1,1)=( 1.0D0, 0.0D0)
  sigma(2,2,1)=( 0.0D0, 0.0D0)
  
  sigma(1,1,2)=( 0.0D0, 0.0D0)
  sigma(1,2,2)=( 0.0D0, 1.0D0)
  sigma(2,1,2)=( 0.0D0,-1.0D0)
  sigma(2,2,2)=( 0.0D0, 0.0D0)
  
  sigma(1,1,3)=(-1.0D0, 0.0D0)
  sigma(1,2,3)=( 0.0D0, 0.0D0)
  sigma(2,1,3)=( 0.0D0, 0.0D0)
  sigma(2,2,3)=( 1.0D0, 0.0D0)
else
  stop '[calc_sigma] wrong mode'
end if


if (verbose==1) then
  write(*,*) '#################################'
  write(*,*) 'calculation OF Pauli matricies'
  write(*,*) '#################################'
  write(*,*) 'sigma_x'
  write(*,'(4f6.2)') sigma(1,1,1), sigma(1,2,1)
  write(*,'(4f6.2)') sigma(2,1,1), sigma(2,2,1)
  
  write(*,*) 'sigma_y'
  write(*,'(4f6.2)') sigma(1,1,2), sigma(1,2,2)
  write(*,'(4f6.2)') sigma(2,1,2), sigma(2,2,2)
  
  write(*,*) 'sigma_z'
  write(*,'(4f6.2)') sigma(1,1,3), sigma(1,2,3)
  write(*,'(4f6.2)') sigma(2,1,3), sigma(2,2,3)
  write(*,*) '#################################'
end if

end subroutine calc_sigma


!-------------------------------------------------------------------------------
!> Summary: Output of exchange interactions to files
!> Author:
!> Category: input-output, physical-observables, KKRimp
!> Deprecated: False 
!> 
!-------------------------------------------------------------------------------
subroutine calccouplingconstants_writeoutJij(natom,Jijmatrix,Aimatrix,density,ITSCF)
use type_density, only: density_type
use  rotaterealspace, only: rotaterealspace_matrix2
implicit none
integer :: natom,ITSCF
double precision :: Jijmatrix(3,3,natom,natom)
double precision :: Aimatrix(3,natom)
type(density_type),allocatable  ::  density(:)
!local
integer :: iatom,jatom
double precision :: Jijmatrix_temp(3,3)


  write(34536254,*) '# '
  write(34536254,*) '################################################################'
  write(34536254,*) '# Iteration number ',ITSCF
  write(34536254,*) '################################################################'
  write(34536256,'(A)') '#iatom, iatom, Jij, Dij_x, Dij_y, Dij_z'

  write(34536258,*) '# '
  write(34536258,*) '################################################################'
  write(34536258,*) '# Iteration number ',ITSCF
  write(34536258,*) '################################################################'
  write(34536259,'(A)') '#iatom, iatom, Jij, Dij_x, Dij_y, Dij_z'

  write(34536268,*) '# '
  write(34536268,*) '################################################################'
  write(34536268,*) '# Iteration number ',ITSCF
  write(34536268,*) '################################################################'
  write(34536268,'(A)') '#iatom, Ai(x) Ai(y) Ai(z)'


  do iatom=1,natom

  write(34536268,'(I5,20E25.14)') iatom,Aimatrix(:,iatom)



    do jatom=1,natom
      write(34536254,'(A,2I3)') '# coupling between',iatom,jatom
      write(34536254,'(A,I7,A,2F8.2)') '# direction of atom',iatom,' theta/phi:', density(iatom)%theta,density(iatom)%phi
      write(34536254,'(A,I7,A,2F8.2)') '# direction of atom',jatom,' theta/phi:', density(jatom)%theta,density(jatom)%phi
      write(34536254,'(A,E12.3)') '# Jij is ',(Jijmatrix(1,1,iatom,jatom)+Jijmatrix(2,2,iatom,jatom)+Jijmatrix(3,3,iatom,jatom))/3.0D0
      write(34536254,'(A,3E12.3)') '# Dij is ',(Jijmatrix(2,3,iatom,jatom)-Jijmatrix(3,2,iatom,jatom))/2.0D0, &
                                          (Jijmatrix(3,1,iatom,jatom)-Jijmatrix(1,3,iatom,jatom))/2.0D0, &
                                          (Jijmatrix(1,2,iatom,jatom)-Jijmatrix(2,1,iatom,jatom))/2.0D0
      write(34536254,'(A,3E12.3)') '# Sij is ',(Jijmatrix(2,3,iatom,jatom)+Jijmatrix(3,2,iatom,jatom))/2.0D0, &
                                          (Jijmatrix(3,1,iatom,jatom)+Jijmatrix(1,3,iatom,jatom))/2.0D0, &
                                          (Jijmatrix(1,2,iatom,jatom)+Jijmatrix(2,1,iatom,jatom))/2.0D0
      write(34536254,'(A)') '# gen. Jij matrix is'
      write(34536254,*) Jijmatrix(1,1,iatom,jatom),Jijmatrix(1,2,iatom,jatom),Jijmatrix(1,3,iatom,jatom)
      write(34536254,*) Jijmatrix(2,1,iatom,jatom),Jijmatrix(2,2,iatom,jatom),Jijmatrix(2,3,iatom,jatom)
      write(34536254,*) Jijmatrix(3,1,iatom,jatom),Jijmatrix(3,2,iatom,jatom),Jijmatrix(3,3,iatom,jatom)
      write(34536254,*) '################################################################'

      write(34536256,'(2I6,12E12.3)') iatom,jatom, &
                                    (Jijmatrix(1,1,iatom,jatom)+Jijmatrix(2,2,iatom,jatom)+Jijmatrix(3,3,iatom,jatom))/3.0D0, &
                                    (Jijmatrix(2,3,iatom,jatom)-Jijmatrix(3,2,iatom,jatom))/2.0D0, &
                                    (Jijmatrix(3,1,iatom,jatom)-Jijmatrix(1,3,iatom,jatom))/2.0D0, &
                                    (Jijmatrix(1,2,iatom,jatom)-Jijmatrix(2,1,iatom,jatom))/2.0D0

       Jijmatrix_temp=Jijmatrix(:,:,iatom,jatom)
       call rotaterealspace_matrix2(Jijmatrix_temp,density(iatom)%theta,density(iatom)%phi,& 
                                                                density(jatom)%theta,density(jatom)%phi,'glob->loc')


      write(34536258,'(A,2I3)') '# coupling between',iatom,jatom
      write(34536258,'(A,I7,A,2F8.2)') '# direction of atom',iatom,' theta/phi:', density(iatom)%theta,density(iatom)%phi
      write(34536258,'(A,I7,A,2F8.2)') '# direction of atom',jatom,' theta/phi:', density(jatom)%theta,density(jatom)%phi
      write(34536258,'(A,E12.3)') '# Jij is ',(Jijmatrix_temp(1,1)+Jijmatrix_temp(2,2)+Jijmatrix_temp(3,3))/3.0D0
      write(34536258,'(A,3E12.3)') '# Dij is ',(Jijmatrix_temp(2,3)-Jijmatrix_temp(3,2))/2.0D0, &
                                          (Jijmatrix_temp(3,1)-Jijmatrix_temp(1,3))/2.0D0, &
                                          (Jijmatrix_temp(1,2)-Jijmatrix_temp(2,1))/2.0D0
      write(34536258,'(A,3E12.3)') '# Sij is ',(Jijmatrix_temp(2,3)+Jijmatrix_temp(3,2))/2.0D0, &
                                          (Jijmatrix_temp(3,1)+Jijmatrix_temp(1,3))/2.0D0, &
                                          (Jijmatrix_temp(1,2)+Jijmatrix_temp(2,1))/2.0D0
      write(34536258,'(A)') '# gen. Jij matrix is'
      write(34536258,*) Jijmatrix_temp(1,1),Jijmatrix_temp(1,2),Jijmatrix_temp(1,3)
      write(34536258,*) Jijmatrix_temp(2,1),Jijmatrix_temp(2,2),Jijmatrix_temp(2,3)
      write(34536258,*) Jijmatrix_temp(3,1),Jijmatrix_temp(3,2),Jijmatrix_temp(3,3)
      write(34536258,*) '################################################################'

      write(34536259,'(2I6,12E12.3)') iatom,jatom, &
                                    (Jijmatrix_temp(1,1)+Jijmatrix_temp(2,2)+Jijmatrix_temp(3,3))/3.0D0, &
                                    (Jijmatrix_temp(2,3)-Jijmatrix_temp(3,2))/2.0D0, &
                                    (Jijmatrix_temp(3,1)-Jijmatrix_temp(1,3))/2.0D0, &
                                    (Jijmatrix_temp(1,2)-Jijmatrix_temp(2,1))/2.0D0


    end do
  end do

end subroutine calccouplingconstants_writeoutJij




end module mod_calccouplingconstants
