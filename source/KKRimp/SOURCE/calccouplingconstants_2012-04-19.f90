module mod_calccouplingconstants
contains 
subroutine calcJijmatrix(gmat,tmat,natom,wez,Jijmatrix)
use type_tmat
use type_gmat
use mod_mathtools

implicit none
type(gmat_type)         :: gmat
type(tmat_type)         :: tmat(natom,1)
integer                 :: natom
double complex          :: wez
double precision        :: Jijmatrix(3,3,natom,natom)

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

! double precision           :: pi
write(1337,*) '##########################################'
write(1337,*) '  starting calculation of the Jij-matrix'
write(1337,*) '##########################################'

if (.not. allocated (gmat%gmat) ) then
  stop '[calcJijmatrix] gmat not allocated'
end if
! pi=3.1415926535897931D0
! print *,pi
Jscal=-1.0D0/pi

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
    
      do lalpha=1,3
      

!         write(92921,*) Gij
!         write(92922,*) Gji
!               write(*,*) 'deltaT_Jij(:,:,kalpha)',ubound(tmat(iatom,1)%deltaT_Jij)
! stop
        deltaTik = tmat(iatom,1)%deltaT_Jij(:,:,kalpha)
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

    Jijmatrix(:,:,iatom,jatom)=Jijmatrix(:,:,iatom,jatom)+dimag(wez*Jijmatrixtemp_complex)

    write(2349875,*) iatom, jatom,Jijmatrixtemp
  end do !jatom

end do !iatom

end subroutine !calcJijmatrix



subroutine calccouplingdeltat(wavefunction,deltaTmat,cellnew,gauntcoeff,theta,phi,lmmax,lmsize,lmax,lmpot,nrmax)
use type_wavefunction
use type_cellnew
use type_gauntcoeff
use mod_mathtools, only: matmat, matmat1T,matmatT1
use mod_vllmat
use mod_chebyshev, only: intcheb_cell
use mod_rotatespinframe  
implicit none
type(wavefunction_type)                   :: wavefunction
type(cell_typenew)                        :: cellnew
type(gauntcoeff_type)                     :: gauntcoeff
double precision                          :: zatom
integer                                   :: lmmax
integer                                   :: lmsize
integer                                   :: lmax
integer                                   :: lmpot
integer                                   :: nrmax
!local
double precision                          :: vecn
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

! lmsize= ubound(wavefunction%Rll,1)
! print *,'lmsize, lmmax',lmsize,lmmax

bpot(:,:,1) = (cellnew%vpotnew(:,:,1)-cellnew%vpotnew(:,:,2))*0.5D0

! do ir=1,nrmax
! write(7776,'(50000E)') bpot(ir,:,1)
! end do

call vllmat(bpotll,bpot,lmax,lmmax,lmpot,1, &
            cellnew%nrmaxnew,gauntcoeff,0.0D0,cellnew%rmeshnew,lmmax,0,2,1,'NS')

!   do ir=1,nrmax
!  write(7777,'(50000E)') bpotll(:,:,ir)
! end do
call calclambda(lambda,theta,phi)
! print *,theta,phi
! print *,lambda(1,:,1)
! print *,lambda(2,:,1)
! 
! print *,lambda(1,:,2)
! print *,lambda(2,:,2)
! 
! print *,lambda(1,:,3)
! print *,lambda(2,:,3)


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
      call intcheb_cell(fntemp,cellnew,RBpotR_integrated(ilm1,ilm2,kspin))
    end do
  end do
  call rotatematrix(RBpotR_integrated(:,:,kspin),theta, phi,lmmax,'loc->glob')
! stop
end do !kspin

!   do kspin=1,3
! !     write(*,*) 'spin ',kspin
!     do ilm1=1,lmsize
!       write(7700-1+kspin,'(50000E)') RBpotR_integrated(ilm1,:,kspin)
!     end do !ilm1
!   end do !kspin
deltaTmat=RBpotR_integrated

! stop
end subroutine !calcBpotspin



subroutine calclambda(lambda,theta,phi)
use mod_rotatespinframe
implicit none
double complex :: lambda(2,2,3)
double precision :: theta, phi
double complex :: sigmatemp(2,2,3),sigma(2,2,3)
integer        :: ispin
call calc_sigma(sigma)
sigmatemp=sigma
do ispin=1,3
  call rotatematrix(sigmatemp(:,:,ispin),theta, phi,1,'glob->loc')
!   lambda(:,:,ispin) = sigmatemp(:,:,ispin)- sigma(:,:,3)
  lambda(:,:,ispin) = sigmatemp(:,:,ispin) !test
end do !ispin
end subroutine calclambda

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
  stop'[calc_sigma] wrong mode'
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




end module mod_calccouplingconstants
