module SourceTerms

contains

subroutine SourceTermSuperVector(lcut,kinenergy,r,J,H,J2,H2)
! calculates the 2-entry vector source term functions J,H,B used in the relativistic Lippmann-Schwinger equations
! input: lcut: maximal l value,  W: relativistic energy, z: function argument
! output: arrays with entries from 1 to Lambdamax*2 for the given argument z

use SpinSphericals ! requires SpinSphericals.f90
use Constants ! requires Constants.f90
use Lebedev ! requires Lebedev.f90
use DiracConfig ! requires DiracConfig.f90
use MOD_BESHANK ! requires beshank.f90

implicit none

! input variables
double precision                :: r
double complex                  :: kinenergy
integer                         :: lcut

! output variables
double complex,dimension(:)     :: J,H,J2,H2


! parameters
integer                         :: Lambdacut, Lambda, kappa, l,m

double complex                  :: energy
double complex                  :: k,prefactor
double complex,allocatable,dimension(:):: HL,JL,NL

integer,allocatable,dimension(:):: KappaArray, LambdabarArray
real,allocatable,dimension(:)   :: MuArray

Lambdacut = getLambdacut(lcut)
call KappaMuArray(lcut,KappaArray,MuArray)
call makeLambdabarArray(lcut,LambdabarArray)

allocate(JL(0:lcut+1))
allocate(HL(0:lcut+1))
allocate(NL(0:lcut+1))

! allocate(J(Lambdacut*2))
! allocate(H(Lambdacut*2))
! allocate(N(Lambdacut*2))
! allocate(J2(Lambdacut*2))
! allocate(H2(Lambdacut*2))
! allocate(N2(Lambdacut*2))

! write(*,*) "lcut=",lcut
! write(*,*) "Lambdacut=",Lambdacut

energy = kinenergy + emass * lightspeed**2
k = sqrt(DCMPLX(energy**2 - emass**2 * lightspeed**4)) / lightspeed / hbar
! write(*,*) "k=",k

if(verbose == 2) then
    write(*,*) "source term r=",r
end if
call BESHANK(HL,JL,DCMPLX(k*r),lcut+1)

m=1
! 1st component
do Lambda = 1,Lambdacut
  kappa = KappaArray(Lambda)
  l = getL(kappa)
!  write(*,*) "Lambda=",Lambda,"l=",l,"  Bessel function J=",JL(l)
  J(m) = JL(l)
  H(m) = HL(l)
!   N(m) = NL(l)
  J2(m) = J(m)
  H2(m) = H(m)
!   N2(m) = N(m)
  m = m+1
end do
! 2nd component
prefactor = i * lightspeed * hbar * k / (energy + emass * lightspeed**2)
! theoretically the expression is equal to prefactor = i * sqrt((energy - emass**2 * lightspeed) / (energy + emass**2 * lightspeed))
! however, this is numerically unstable, so the method used for the calculation is preferable
! write(*,*) 
! write(*,*) "second component with prefactor=",prefactor
! write(*,*) 
do Lambda = 1,Lambdacut
  kappa = KappaArray(Lambda)
  l = getLbar(kappa)
  ! write(*,*) "Lambda=",Lambda,"l=",l,"  Bessel function J=",JL(l)
  if (l < 0) then
    J(m) = 0
    H(m) = 0
!     N(m) = 0
  else
    J(m) = prefactor * SIGN(1,kappa) * JL(l)
    H(m) = prefactor * SIGN(1,kappa) * HL(l)
!     N(m) = prefactor * SIGN(1,kappa) * NL(l)
    J2(m) = -J(m)
    H2(m) = -H(m)
!     N2(m) = -N(m)
  end if
  m = m+1
end do

H  = -i*H
H2 = -i*H2
end subroutine SourceTermSuperVector

end module SourceTerms