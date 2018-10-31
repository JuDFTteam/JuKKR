module Potential
! Version mit Speedup durch Vermeiden doppelter Integrationen
! zusätzlicher Speedup durch Berücksichtigung von verschwindenden rel. Gaunt-Koeffizienten

use mod_datatypes, only: dp

contains

subroutine PotentialMatrixArray(lcut,lcut_input,zatom,meshpoints,nrmax,kinenergy,VLLin,PotMatrixArray)
use SpinSphericals
use Lebedev
use DiracConfig
implicit none

! input 
integer                                     :: lcut, lcut_input, nrmax
complex (kind=dp)                              :: kinenergy
double precision                            :: meshpoints(nrmax)
double precision                            :: zatom,VLLin(nrmax,(2*lcut+1)**2,2)

! output
complex (kind=dp)                              :: PotMatrixArray(:,:,:)

! other variables
integer                                     :: k
complex (kind=dp), allocatable, dimension(:,:)   :: PotMatrix
complex (kind=dp), allocatable                  :: chi1array(:,:), chi2array(:,:)
double precision                            :: theta_array(integrationpoints),phi_array(integrationpoints), weight_array(integrationpoints)
integer, allocatable                         :: DcoeffIndexListA(:,:), DcoeffIndexListB(:,:)
complex (kind=dp), allocatable                  :: DcoeffListA(:), DcoeffListB(:)

if(verbose > 0) then
    write(*,*) ""
    write(*,*) ""
    write(*,*) "  #############################################################"
    write(*,*) "  #                                                           #"
    write(*,*) "  #          FULL POTENTIAL SINGLE-SITE DIRAC SOLVER          #"
    write(*,*) "  #          Pascal Kordt, July 2010 - September 2011         #"
    write(*,*) "  #                                                           #"
    write(*,*) "  #############################################################"
    write(*,*) ""
    write(*,*) "kinenergy=",kinenergy
    write(*,*) "lcut=",lcut
    write(*,*) "zatom=",zatom
    write(*,*) "nrmax=",nrmax
    write(*,*) ""
end if

call makeSpinSphericalArray(lcut,chi1array,chi2array)
call lebedevarray(theta_array,phi_array,weight_array)


! wPot-Version
if(spherical_only /= 1) then
    if(verbose > 0) then
        write(*,*) "******* full potential mode *******"
    end if
    call readDcoeff(lcut,DcoeffIndexListA,DcoeffIndexListB,DcoeffListA,DcoeffListB)
else
    if(verbose > 0) then
        write(*,*) "*** quick mode for spherical potentials without magnetic field B ***"
    end if
end if

write(*,*) ""
write(*,*) "Starting to calculate the w-coefficients for each r-value of the mesh"
do k=1,nrmax
    if(verbose >1 ) then
        write(*,*)
        write(*,*) "MESH POINT ",k," OF ",nrmax
     end if
!    r = meshpoints(k)
    call PotentialSuperMatrix(lcut,lcut_input,zatom,meshpoints,k,kinenergy,chi1array,chi2array,theta_array,phi_array,weight_array,DcoeffIndexListA,DcoeffIndexListB,DcoeffListA,DcoeffListB,VLLin,PotMatrix)
    PotMatrixArray(:,:,k) = PotMatrix(:,:)
end do
write(*,*) "done. (w-coefficients)"

end subroutine PotentialMatrixArray


subroutine PotentialSuperMatrix(lcut,lcut_input,zatom,meshpoints,nr,kinenergy,chi1array,chi2array,theta_array,phi_array,weight_array,DcoeffIndexListA,DcoeffIndexListB,DcoeffListA,DcoeffListB,VLLin,PotMatrix)
use Constants
use SpinSphericals
use Lebedev
use DiracConfig
implicit none

! l cut-off for the expansion in spin spherical harmonics
integer                         :: lcut
double precision                :: zatom,meshpoints(:)
integer                         :: nr

! l cut-off for the input potential, expanded in spherical harmonics
integer                         :: lcut_input
double precision                :: r
complex (kind=dp)                  :: kinenergy
complex (kind=dp), allocatable      :: chi1array(:,:), chi2array(:,:)
double precision                :: VLLin(:,:,:)


complex (kind=dp)                  :: energy
integer                         :: Lambdacut,Lambda1,Lambda2
complex (kind=dp), allocatable, dimension(:,:,:) :: wcoeff
complex (kind=dp), allocatable, dimension(:,:) :: PotMatrix

double precision                 :: theta_array(integrationpoints), phi_array(integrationpoints), weight_array(integrationpoints)
integer, allocatable              :: DcoeffIndexListA(:,:), DcoeffIndexListB(:,:)
complex (kind=dp), allocatable       :: DcoeffListA(:), DcoeffListB(:)

r = meshpoints(nr)

Lambdacut = getLambdacut(lcut)
if(.NOT. allocated(PotMatrix)) then
    allocate(PotMatrix(2*Lambdacut,2*Lambdacut))
end if

energy = kinenergy + emass * lightspeed**2

if(spherical_only == 1) then
    ! accelerated calculation for spherical potentials with B=0
    call wPotentialExpansionNoB(lcut,lcut_input,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,DcoeffIndexListA,DcoeffIndexListB,DcoeffListA,DcoeffListB,VLLin,wcoeff)
else
    ! full-potential calculation
    call wPotentialExpansion(lcut,lcut_input,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,DcoeffIndexListA,DcoeffIndexListB,DcoeffListA,DcoeffListB,VLLin,wcoeff)
end if


! Block a
do Lambda1=1,Lambdacut
do Lambda2=1,Lambdacut
    PotMatrix(Lambda1,Lambda2) = wcoeff(Lambda1,Lambda2,1)
end do
end do

! Block b
do Lambda1=1,Lambdacut
do Lambda2=1,Lambdacut
    PotMatrix(Lambda1,Lambda2+Lambdacut) = wcoeff(Lambda1,Lambda2,2)
end do
end do

! Block c
do Lambda1=1,Lambdacut
do Lambda2=1,Lambdacut
   PotMatrix(Lambda1+Lambdacut,Lambda2) = wcoeff(Lambda1,Lambda2,3)
end do
end do

! Block d
do Lambda1=1,Lambdacut
do Lambda2=1,Lambdacut
    PotMatrix(Lambda1+Lambdacut,Lambda2+Lambdacut) = wcoeff(Lambda1,Lambda2,4)
end do
end do

! write(*,*) "(W+mc²)/c²/hbar²=",(energy + emass * lightspeed**2) / lightspeed**2 / hbar**2
PotMatrix = PotMatrix * (energy + emass * lightspeed**2) / lightspeed**2 / hbar**2

end subroutine PotentialSuperMatrix


subroutine readDcoeff(lcut,DcoeffIndexListA,DcoeffIndexListB,DcoeffListA,DcoeffListB)
use SpinSphericals ! requires SpinSphericals.f90
use Constants ! requires Constants.f90
use Lebedev ! requires Lebedev.f90
use DiracConfig ! requires DiracConfig.f90
use RelativisticGauntCoefficients ! requires RelativisticGauntCoefficients.f90

implicit none

! input
! l cut-off for the expansion in spin spherical harmonics
integer                         :: lcut
integer                         :: Lambdacut
integer, allocatable             :: DcoeffIndexListA(:,:), DcoeffIndexListB(:,:)
complex (kind=dp), allocatable      :: DcoeffListA(:), DcoeffListB(:)
double precision                :: temp_realpart, temp_imagpart
character(len=100)              :: filename
character(len=2)                :: lcutstring
integer                         :: currentline, filelengthA, filelengthB

integer, allocatable, dimension(:):: KappaArray, LambdabarArray
real, allocatable, dimension(:)   :: MuArray

call KappaMuArray(lcut,KappaArray,MuArray)
call makeLambdabarArray(lcut,LambdabarArray)
Lambdacut = LambdabarArray(getLambdacut(lcut))

! read the relativistic Gaunt-coefficients D from the files
write(lcutstring,'(I02.2)') lcut
filename = 'RelativisticGauntCoeffA_lcut'//trim(lcutstring)//'.txt'
if(CheckExistence(filename) /= 1) then
    write(*,*) ""
    write(*,*) "The files containing the relativistic equivalence of the Gaunt coefficients"
    write(*,*) "could not be found in your working directory and will be created now."
    write(*,*) "This can take a while..."
    write(*,*) "Don't worry, this procedure has to be done only once. In future calculations"
    write(*,*) "the file created now will be used."
    write(*,*) ""
    call writeDcoeff(lcut)
end if

if(verbose > 0) then
    write(*,*) "Reading the 1st set of D-coeffients from ", filename
end if
filelengthA = GetNumberOfLines(filename)
open(unit=1,file=filename)

filename = 'RelativisticGauntCoeffB_lcut'//trim(lcutstring)//'.txt'
if(CheckExistence(filename) /= 1) then
    write(*,*) ""
    write(*,*) "The files containing the relativistic equivalence of the Gaunt coefficients"
    write(*,*) "could not be found in your working directory and will be created now."
    write(*,*) "This can take a while..."
    write(*,*) "Don't worry, this procedure has to be done only once. In future calculations"
    write(*,*) "the file created now will be used."
    write(*,*) ""
    call writeDcoeff(lcut)
end if
if(verbose > 0) then
    write(*,*) "Reading the 2nd set of D-coeffients from ", filename
end if
filelengthB = GetNumberOfLines(filename)
open(unit=2,file=filename)

! allocate arrays
allocate(DcoeffIndexListA(4,filelengthA))
allocate(DcoeffIndexListB(4,filelengthB))
allocate(DcoeffListA(filelengthA))
allocate(DcoeffListB(filelengthB))

if(verbose > 0) then
    write(*,*) "Preparing the D-coefficients array..."
end if

do currentline=1,filelengthA
    read(1,'(I3.3, " ",  I3.3, " ", I3.3, " ", I3.3, " ", E23.16, " ", E23.16)') DcoeffIndexListA(1,currentline), DcoeffIndexListA(2,currentline), DcoeffIndexListA(3,currentline), DcoeffIndexListA(4,currentline), temp_realpart, temp_imagpart
    DcoeffListA(currentline) = temp_realpart + i*temp_imagpart
end do

do currentline=1,filelengthB
    read(2,'(I3.3, " ",  I3.3, " ", I3.3, " ", I3.3, " ", E23.16, " ", E23.16)') DcoeffIndexListB(1,currentline), DcoeffIndexListB(2,currentline), DcoeffIndexListB(3,currentline), DcoeffIndexListB(4,currentline), temp_realpart, temp_imagpart
    DcoeffListB(currentline) = temp_realpart + i*temp_imagpart
end do

close(1)
close(2)

if(verbose > 0) then
    write(*,*) "done."
end if
end subroutine readDcoeff


subroutine wPotentialExpansionNoB(lcut,lcut_input,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,DcoeffIndexListA,DcoeffIndexListB,DcoeffListA,DcoeffListB,VLLin,wcoeff)
use Lebedev
use SpinSphericals
use Constants
use DiracConfig

implicit none

! input
! l cut-off for the expansion in spin spherical harmonics
integer                         :: lcut
integer, allocatable             :: DcoeffIndexListA(:,:), DcoeffIndexListB(:,:)
complex (kind=dp), allocatable      :: DcoeffListA(:), DcoeffListB(:)
double precision                :: zatom,VLLin(:,:,:)
double precision                :: meshpoints(:)
integer                         :: nr


! l cut-off for the input potential, expanded in spherical harmonics
integer                         :: lcut_input
double precision                :: r
complex (kind=dp), allocatable      :: chi1array(:,:), chi2array(:,:)
double precision                :: theta_array(integrationpoints), phi_array(integrationpoints), weight_array(integrationpoints)

double precision                :: potPhi
integer                         :: Lambda1, Lambdacut
complex (kind=dp), allocatable      :: wcoeff(:,:,:)
integer, allocatable, dimension(:):: KappaArray, LambdabarArray
real, allocatable, dimension(:)   :: MuArray

r = meshpoints(nr)

Lambdacut = getLambdacut(lcut)
call KappaMuArray(lcut,KappaArray,MuArray)
call makeLambdabarArray(lcut,LambdabarArray)

call getPotPhi(zatom,meshpoints,nr,0.0d0,0.0d0,VLLin,potPhi)

allocate(wcoeff(Lambdacut,Lambdacut,4))
wcoeff(:,:,:) = 0.0d0

if(verbose==2) then
    write(*,*) "computing the w-coefficients"
end if

do Lambda1=1,Lambdacut
    ! only diagonal entries are /= 0
    wcoeff(Lambda1,Lambda1,1) = echarge * potPhi 
    wcoeff(Lambda1,Lambda1,4) = echarge * potPhi 
end do
end subroutine wPotentialExpansionNoB


subroutine wPotentialExpansion(lcut,lcut_input,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,DcoeffIndexListA,DcoeffIndexListB,DcoeffListA,DcoeffListB,VLLin,wcoeff)
use SpinSphericals ! requires SpinSphericals.f90
use Constants ! requires Constants.f90
use Lebedev ! requires Lebedev.f90
use DiracConfig ! requires DiracConfig.f90

implicit none

! input
! l cut-off for the expansion in spin spherical harmonics
integer                         :: lcut
integer, allocatable             :: DcoeffIndexListA(:,:), DcoeffIndexListB(:,:)
complex (kind=dp), allocatable      :: DcoeffListA(:), DcoeffListB(:)

double precision                :: VLLin(:,:,:)
double precision                :: zatom,meshpoints(:)
integer                         :: nr


! l cut-off for the input potential, expanded in spherical harmonics
integer                         :: lcut_input
double precision                :: r
complex (kind=dp), allocatable      :: chi1array(:,:), chi2array(:,:)
double precision                :: theta_array(integrationpoints), phi_array(integrationpoints), weight_array(integrationpoints)

integer                         :: Lambdacut, Lambda1, Lambda2, Lambda3, Lambda4, test_zero, test_nonzero
integer                         :: linecount, filelengthA, filelengthB
complex (kind=dp), allocatable      :: vcoeff(:,:,:), wcoeff(:,:,:)


integer, allocatable, dimension(:):: KappaArray, LambdabarArray
real, allocatable, dimension(:)   :: MuArray

r = meshpoints(nr)

Lambdacut = getLambdacut(lcut)
call KappaMuArray(lcut,KappaArray,MuArray)
call makeLambdabarArray(lcut,LambdabarArray)

allocate(wcoeff(Lambdacut,Lambdacut,4))
wcoeff(:,:,:) = 0.0d0

call PotentialExpansion(lcut,lcut_input,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,VLLin,vcoeff)
if(verbose==2) then
    write(*,*) "computing sums for the the w-coefficients..."
end if

filelengthA = size(DcoeffListA(:))
filelengthB = size(DcoeffListB(:))

do linecount = 1, filelengthA
    Lambda1 = DcoeffIndexListA(1,linecount)
    Lambda2 = DcoeffIndexListA(2,linecount)
    Lambda3 = DcoeffIndexListA(3,linecount)
    Lambda4 = DcoeffIndexListA(4,linecount)
    wcoeff(Lambda1,Lambda4,1) = wcoeff(Lambda1,Lambda4,1) + DcoeffListA(linecount) * vcoeff(Lambda2,Lambda3,1)
end do

do linecount = 1, filelengthB
    Lambda1 = DcoeffIndexListB(1,linecount)
    Lambda2 = DcoeffIndexListB(2,linecount)
    Lambda3 = DcoeffIndexListB(3,linecount)
    Lambda4 = DcoeffIndexListB(4,linecount)
    wcoeff(Lambda1,Lambda4,4) = wcoeff(Lambda1,Lambda4,4) + DcoeffListB(linecount) * vcoeff(Lambda2,Lambda3,4)
end do

! clean up
test_zero=0
test_nonzero=0
do Lambda1=1,Lambdacut
do Lambda4=1,Lambdacut
    test_nonzero=test_nonzero+2
    if(w_eps > 0 .and. abs(wcoeff(Lambda1,Lambda4,1)) < w_eps) then
        wcoeff(Lambda1,Lambda4,1)=0.0d0
        test_zero=test_zero+1
        test_nonzero=test_nonzero-1
    end if
    if(w_eps > 0 .and. abs(wcoeff(Lambda1,Lambda4,4)) < w_eps) then
        wcoeff(Lambda1,Lambda4,4)=0.0d0
        test_zero=test_zero+1
        test_nonzero=test_nonzero-1
    end if
end do
end do

if(verbose==2) then
    write(*,*) "done."
end if

end subroutine wPotentialExpansion


subroutine PotentialExpansion(lcut,lcut_input,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,VLLin,vcoeff)
use SpinSphericals ! requires SpinSphericals.f90
use Constants ! requires Constants.f90
use Lebedev ! requires Lebedev.f90
use DiracConfig ! requires DiracConfig.f90

implicit none

! input
! l cut-off for the expansion in spin spherical harmonics
integer                         :: lcut
double precision                :: VLLin(:,:,:)
double precision                :: zatom,meshpoints(:)
integer                         :: nr

! l cut-off for the input potential, expanded in spherical harmonics
integer                         :: lcut_input
double precision                :: r
complex (kind=dp), allocatable      :: chi1array(:,:), chi2array(:,:)
double precision                :: theta_array(integrationpoints), phi_array(integrationpoints), weight_array(integrationpoints)

! variables
double precision                :: theta,phi


! variables
integer                         :: j,k
double precision                :: deviation
double precision                :: potPhi,potBx,potBy,potBz
double precision                :: eval_phi,eval_theta
integer                         :: out_potprecision
complex (kind=dp), dimension(4,4)  :: potMatrixIn,potMatrix,differenceMatrix ! 4x4 potential matrix

integer                         :: Lambdacut, Lambda1, Lambda2
complex (kind=dp), dimension(4)    :: v
complex (kind=dp), allocatable, dimension(:,:,:) :: vcoeff

integer, allocatable, dimension(:):: KappaArray, LambdabarArray
real, allocatable, dimension(:)   :: MuArray
complex (kind=dp)                  :: chi11,chi12,chi21,chi22
complex (kind=dp), allocatable        :: nuL(:,:,:), nuR(:,:,:)


! test options
! 1: write the precision of the potential expansion for a sample point, 0: don't write it
out_potprecision = 1
! sample point to evaluate
eval_phi   = 0.3d0*Pi
eval_theta = 0.1d0*Pi

r = meshpoints(nr)

Lambdacut = getLambdacut(lcut)
call KappaMuArray(lcut,KappaArray,MuArray)
call makeLambdabarArray(lcut,LambdabarArray)
allocate(vcoeff(Lambdacut,Lambdacut,4))

if(verbose > 1) then
    write(*,*) "computing potential expansion..."
end if

call getPotPhi(zatom,meshpoints,nr,eval_theta,eval_phi,VLLin,potPhi)
call getPotB(zatom,meshpoints,nr,eval_theta,phi,VLLin,potBx,potBy,potBz)
! potPhi now contains the value for the scalar potential Phi and potBx,potBy,potBz the B-Field values


! calculate the potential matrix (this has to be changed in case of a full vector poential A)
potMatrixIn(1,1) = echarge * potPhi - muB * potBz
potMatrixIn(1,2) = - muB * potBx  + i * muB * potBy
potMatrixIn(2,1) = - muB * potBx - i * muB * potBy
potMatrixIn(2,2) = echarge * potPhi + muB * potBz

potMatrixIn(1,3) = 0.0d0
potMatrixIn(1,4) = 0.0d0
potMatrixIn(2,3) = 0.0d0
potMatrixIn(2,4) = 0.0d0

potMatrixIn(3,1) = 0.0d0
potMatrixIn(3,2) = 0.0d0
potMatrixIn(4,1) = 0.0d0
potMatrixIn(4,2) = 0.0d0

potMatrixIn(3,3) = echarge * potPhi + muB * potBz
potMatrixIn(3,4) = muB * potBx - i * muB * potBy
potMatrixIn(4,3) = muB * potBx + i * muB * potBy
potMatrixIn(4,4) = echarge * potPhi - muB * potBz

call nuCoefficients(lcut,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,VLLin,nuL,nuR)

do Lambda1=1,Lambdacut
do Lambda2=1,Lambdacut
  call ExpansionCoefficients(Lambda1,Lambda2,lcut,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,nuL,nuR,v)

  do j=1,4
    ! clean up
    if(v_eps > 0 .and. abs(v(j)) < v_eps) then
        v(j)=0.0d0
    end if
  vcoeff(Lambda1,Lambda2,j) = v(j)
  end do
end do
end do


if(verbose>1) then
! calculate the potential matrix and check accuracy
  potMatrix(:,:) = 0
  do Lambda1=1,Lambdacut
  do Lambda2=1,Lambdacut

    call SpinSphericalHarmonic(KappaArray(Lambda1),MuArray(Lambda1),theta,phi,chi11,chi12)
    call SpinSphericalHarmonic(KappaArray(Lambda2),MuArray(Lambda2),theta,phi,chi21,chi22)
    potMatrix(1,1) = potMatrix(1,1) + vcoeff(Lambda1,Lambda2,1) * chi11 * conjg(chi21)
    potMatrix(1,2) = potMatrix(1,2) + vcoeff(Lambda1,Lambda2,1) * chi11 * conjg(chi22)
    potMatrix(2,1) = potMatrix(2,1) + vcoeff(Lambda1,Lambda2,1) * chi12 * conjg(chi21)
    potMatrix(2,2) = potMatrix(2,2) + vcoeff(Lambda1,Lambda2,1) * chi12 * conjg(chi22)

    call SpinSphericalHarmonic(KappaArray(Lambda1),MuArray(Lambda1),theta,phi,chi11,chi12)
    call SpinSphericalHarmonic(KappaArray(LambdabarArray(Lambda2)),MuArray(LambdabarArray(Lambda2)),theta,phi,chi21,chi22)
    potMatrix(1,3) = potMatrix(1,3) + vcoeff(Lambda1,Lambda2,2) * chi11 * conjg(chi21)
    potMatrix(1,4) = potMatrix(1,4) + vcoeff(Lambda1,Lambda2,2) * chi11 * conjg(chi22)
    potMatrix(2,3) = potMatrix(2,3) + vcoeff(Lambda1,Lambda2,2) * chi12 * conjg(chi21)
    potMatrix(2,4) = potMatrix(2,4) + vcoeff(Lambda1,Lambda2,2) * chi12 * conjg(chi22)

    call SpinSphericalHarmonic(KappaArray(LambdabarArray(Lambda1)),MuArray(LambdabarArray(Lambda1)),theta,phi,chi11,chi12)
    call SpinSphericalHarmonic(KappaArray(Lambda2),MuArray(Lambda2),theta,phi,chi21,chi22)
    potMatrix(3,1) = potMatrix(3,1) + vcoeff(Lambda1,Lambda2,3) * chi11 * conjg(chi21)
    potMatrix(3,2) = potMatrix(3,2) + vcoeff(Lambda1,Lambda2,3) * chi11 * conjg(chi22)
    potMatrix(4,1) = potMatrix(4,1) + vcoeff(Lambda1,Lambda2,3) * chi12 * conjg(chi21)
    potMatrix(4,2) = potMatrix(4,2) + vcoeff(Lambda1,Lambda2,3) * chi12 * conjg(chi22)

    call SpinSphericalHarmonic(KappaArray(LambdabarArray(Lambda1)),MuArray(LambdabarArray(Lambda1)),theta,phi,chi11,chi12)
    call SpinSphericalHarmonic(KappaArray(LambdabarArray(Lambda2)),MuArray(LambdabarArray(Lambda2)),theta,phi,chi21,chi22)
    potMatrix(3,3) = potMatrix(3,3) + vcoeff(Lambda1,Lambda2,4) * chi11 * conjg(chi21)
    potMatrix(3,4) = potMatrix(3,4) + vcoeff(Lambda1,Lambda2,4) * chi11 * conjg(chi22)
    potMatrix(4,3) = potMatrix(4,3) + vcoeff(Lambda1,Lambda2,4) * chi12 * conjg(chi21)
    potMatrix(4,4) = potMatrix(4,4) + vcoeff(Lambda1,Lambda2,4) * chi12 * conjg(chi22)
  ! write(*,*) "Lambda1=",Lambda1,"Lambda2=",Lambda2,"v4-coefficient",vcoeff(4,Lambda1,Lambda2)

  end do
  end do

  ! calculate the sum of all absolute values of the difference matrix (matrix 1-norm)
  deviation = 0
  differenceMatrix = potMatrixIn-potMatrix
  do j=1,4
  do k=1,4
    deviation = deviation + abs(differenceMatrix(j,k))
  end do
  end do
  write(*,*) "done. (expansion accuracy: ", deviation, ")"
end if

end subroutine PotentialExpansion


subroutine PotentialExpansionNoB(lcut,lcut_input,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,VLLin,vcoeff)
! obsolete
use SpinSphericals ! requires SpinSphericals.f90
use Constants ! requires Constants.f90
use Lebedev ! requires Lebedev.f90

implicit none

! input
! l cut-off for the expansion in spin spherical harmonics
integer                         :: lcut
double precision                :: zatom,meshpoints(:)
integer                         :: nr

! l cut-off for the input potential, expanded in spherical harmonics
integer                         :: lcut_input
double precision                :: r
complex (kind=dp), allocatable      :: chi1array(:,:), chi2array(:,:)
double precision                :: theta_array(integrationpoints), phi_array(integrationpoints), weight_array(integrationpoints)
double precision                :: VLLin(:,:,:)

! variables
double precision                :: theta,phi

! variables
integer                         :: j,k
double precision                :: deviation
double precision                :: potPhi
double precision                :: eval_phi,eval_theta
complex (kind=dp), dimension(4,4)  :: potMatrixIn,potMatrix,differenceMatrix ! 4x4 potential matrix

integer                         :: Lambdacut, Lambda1, Lambda2
complex (kind=dp), allocatable, dimension(:,:,:) :: vcoeff

integer, allocatable, dimension(:):: KappaArray, LambdabarArray
real, allocatable, dimension(:)   :: MuArray
complex (kind=dp)                  :: chi11,chi12,chi21,chi22

r = meshpoints(nr)

Lambdacut = getLambdacut(lcut)
call KappaMuArray(lcut,KappaArray,MuArray)
call makeLambdabarArray(lcut,LambdabarArray)
allocate(vcoeff(Lambdacut,Lambdacut,4))

vcoeff(:,:,:) = 0.0d0

write(*,*) "Starting to calculate the potential expansion for r=",r
write(*,*) "!!! Using the PotentialExpansionNoB subroutine. This subroutine is obsolete and should only be used for testing! !!!"

! sample point to evaluate
eval_phi   = 0.3d0*Pi
eval_theta = 0.1d0*Pi

call getPotPhi(zatom,meshpoints,nr,eval_theta,eval_phi,VLLin,potPhi)

! calculate the potential matrix (this has to be changed in case of a full vector poential A)
potMatrixIn(1,1) = echarge * PotPhi
potMatrixIn(1,2) = 0.0d0
potMatrixIn(2,1) = 0.0d0
potMatrixIn(2,2) = echarge * PotPhi

potMatrixIn(1,3) = 0.0d0
potMatrixIn(1,4) = 0.0d0
potMatrixIn(2,3) = 0.0d0
potMatrixIn(2,4) = 0.0d0

potMatrixIn(3,1) = 0.0d0
potMatrixIn(3,2) = 0.0d0
potMatrixIn(4,1) = 0.0d0
potMatrixIn(4,2) = 0.0d0

potMatrixIn(3,3) = echarge * PotPhi
potMatrixIn(3,4) = 0.0d0
potMatrixIn(4,3) = 0.0d0
potMatrixIn(4,4) = echarge * PotPhi

! The only entries /= 0 are for Lambda1=Lambda2=1,2 and j=1,4(a block and d block)
! Lambda1=Lambda2=1
vcoeff(1,1,1) = 4*Pi*echarge*PotPhi 
vcoeff(1,1,4) = 4*Pi*echarge*PotPhi 
! Lambda1=Lambda2=2
vcoeff(2,2,1) = 4*Pi*echarge*PotPhi 
vcoeff(2,2,4) = 4*Pi*echarge*PotPhi 

! calculate the potential matrix
potMatrix(:,:) = 0
do Lambda1=1,Lambdacut
do Lambda2=1,Lambdacut

  call SpinSphericalHarmonic(KappaArray(Lambda1),MuArray(Lambda1),theta,phi,chi11,chi12)
  call SpinSphericalHarmonic(KappaArray(Lambda2),MuArray(Lambda2),theta,phi,chi21,chi22)
  potMatrix(1,1) = potMatrix(1,1) + vcoeff(Lambda1,Lambda2,1) * chi11 * conjg(chi21)
  potMatrix(1,2) = potMatrix(1,2) + vcoeff(Lambda1,Lambda2,1) * chi11 * conjg(chi22)
  potMatrix(2,1) = potMatrix(2,1) + vcoeff(Lambda1,Lambda2,1) * chi12 * conjg(chi21)
  potMatrix(2,2) = potMatrix(2,2) + vcoeff(Lambda1,Lambda2,1) * chi12 * conjg(chi22)

  call SpinSphericalHarmonic(KappaArray(Lambda1),MuArray(Lambda1),theta,phi,chi11,chi12)
  call SpinSphericalHarmonic(KappaArray(LambdabarArray(Lambda2)),MuArray(LambdabarArray(Lambda2)),theta,phi,chi21,chi22)
  potMatrix(1,3) = potMatrix(1,3) + vcoeff(Lambda1,Lambda2,2) * chi11 * conjg(chi21)
  potMatrix(1,4) = potMatrix(1,4) + vcoeff(Lambda1,Lambda2,2) * chi11 * conjg(chi22)
  potMatrix(2,3) = potMatrix(2,3) + vcoeff(Lambda1,Lambda2,2) * chi12 * conjg(chi21)
  potMatrix(2,4) = potMatrix(2,4) + vcoeff(Lambda1,Lambda2,2) * chi12 * conjg(chi22)

  call SpinSphericalHarmonic(KappaArray(LambdabarArray(Lambda1)),MuArray(LambdabarArray(Lambda1)),theta,phi,chi11,chi12)
  call SpinSphericalHarmonic(KappaArray(Lambda2),MuArray(Lambda2),theta,phi,chi21,chi22)
  potMatrix(3,1) = potMatrix(3,1) + vcoeff(Lambda1,Lambda2,3) * chi11 * conjg(chi21)
  potMatrix(3,2) = potMatrix(3,2) + vcoeff(Lambda1,Lambda2,3) * chi11 * conjg(chi22)
  potMatrix(4,1) = potMatrix(4,1) + vcoeff(Lambda1,Lambda2,3) * chi12 * conjg(chi21)
  potMatrix(4,2) = potMatrix(4,2) + vcoeff(Lambda1,Lambda2,3) * chi12 * conjg(chi22)

  call SpinSphericalHarmonic(KappaArray(LambdabarArray(Lambda1)),MuArray(LambdabarArray(Lambda1)),theta,phi,chi11,chi12)
  call SpinSphericalHarmonic(KappaArray(LambdabarArray(Lambda2)),MuArray(LambdabarArray(Lambda2)),theta,phi,chi21,chi22)
  potMatrix(3,3) = potMatrix(3,3) + vcoeff(Lambda1,Lambda2,4) * chi11 * conjg(chi21)
  potMatrix(3,4) = potMatrix(3,4) + vcoeff(Lambda1,Lambda2,4) * chi11 * conjg(chi22)
  potMatrix(4,3) = potMatrix(4,3) + vcoeff(Lambda1,Lambda2,4) * chi12 * conjg(chi21)
  potMatrix(4,4) = potMatrix(4,4) + vcoeff(Lambda1,Lambda2,4) * chi12 * conjg(chi22)

end do
end do

! calculate the sum of all absolute values of the difference matrix (matrix 1-norm)
deviation = 0
differenceMatrix = potMatrixIn-potMatrix
do j=1,4
do k=1,4
  deviation = deviation + abs(differenceMatrix(j,k))
end do
end do
write(*,*) "accuracy of the potential expansion (sum of the absolute deviations):", deviation
write(*,*) "...finished calculating the potential expansion for r=",r

end subroutine PotentialExpansionNoB


subroutine makeSpinSphericalArray(lcut,chi1array,chi2array)
  use SpinSphericals
  use Lebedev
  implicit none
  ! input: lcut
  integer                               :: lcut
  ! output: chi1array, chi2array
  complex (kind=dp), allocatable            :: chi1array(:,:), chi2array(:,:)
  
  complex (kind=dp)                        :: chi1,chi2
  double precision                      :: x,y,z,r,theta,phi,weight
  integer                               :: kappa
  real                                  :: mu
  
  integer                               :: Lambdacut, j, Lambda
  integer, allocatable, dimension(:)      :: KappaArray
  real, allocatable, dimension(:)         :: MuArray
  
  Lambdacut = getLambdacut(lcut+1)
  call KappaMuArray(lcut,KappaArray,MuArray)
  
  allocate(chi1array(integrationpoints,Lambdacut))
  allocate(chi2array(integrationpoints,Lambdacut))

  do j=1,integrationpoints
      call lebedevpoint(j,x,y,z,weight)
      call Cartesian2Spherical(x,y,z,r,theta,phi)

      do Lambda=1,Lambdacut
          kappa = KappaArray(Lambda)
          mu    = MuArray(Lambda)
          call SpinSphericalHarmonic(kappa,mu,theta,phi,chi1,chi2)
          chi1array(j,Lambda) = chi1
          chi2array(j,Lambda) = chi2
      end do
  end do
end subroutine makeSpinSphericalArray


subroutine nuCoefficients(lcut,zatom,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,VLLin,nuL,nuR)
! input: lcut,position r,SphiSpherical-Arrays chi1array,chi2array
! output nuL(Lambda,matrixindex,eigenvalueindex)
!        nuR(Lambda,matrixindex,eigenvalueindex)
  use SpinSphericals
  use Lebedev
  use DiracConfig
  implicit none
  integer                :: Lambda
  double precision                      :: zatom,meshpoints(:)
  integer                               :: nr

  double precision            :: r_fix
  double precision                      :: theta_array(integrationpoints), phi_array(integrationpoints), weight_array(integrationpoints)
  double precision                      :: VLLin(:,:,:)


  complex (kind=dp), allocatable        :: nuL(:,:,:), nuR(:,:,:)
  integer                :: j
  ! ^ index 1..4 corresponds to a,b,c,d i.e. the different sub-matrices, index i=1,2 corresponds to the different eigenvalues of the respective sub-matrix 
  complex (kind=dp), dimension(4,2)    :: eigenvalue
  complex (kind=dp), dimension(4,2,2)    :: eigenvector
  complex (kind=dp)            :: chi1,chi2
  double precision            :: potPhi,potBx,potBy,potBz
  integer                :: matrixindex,eigenvalueindex
 
  integer                :: lcut, Lambdacut
  integer, allocatable, dimension(:)    :: KappaArray, LambdabarArray
  real, allocatable, dimension(:)        :: MuArray
  complex (kind=dp), allocatable            :: chi1array(:,:), chi2array(:,:)

  r_fix = meshpoints(nr)

  Lambdacut = getLambdacut(lcut)
  call makeLambdabarArray(lcut,LambdabarArray)
  call KappaMuArray(lcut,KappaArray,MuArray)

  allocate(nuL(Lambdacut,4,2))
  allocate(nuR(Lambdacut,4,2))

  nuL(:,:,:)  = 0
  nuR(:,:,:)  = 0

  do j=1,integrationpoints

     call getPotPhi(zatom,meshpoints,nr,theta_array(j),phi_array(j),VLLin,potPhi)
     call getPotB(zatom,meshpoints,nr,theta_array(j),phi_array(j),VLLin,potBx,potBy,potBz) 

    ! get eigenvalue entries of the 2x2 sub-matrices of the 4x4 potential matrix
    call SubMatrixEigenvalue(potPhi,potBx,potBy,potBz,eigenvalue)
    ! get eigenvector entries of the 2x2 sub-matrices of the 4x4 potential matrix
    call SubMatrixEigenvector(potPhi,potBx,potBy,potBz,eigenvector)

    do Lambda=1,Lambdacut
        ! write entries of the Spin Spherical Harmonics into chi1 (1st component) and chi2 (2nd component)
        do matrixindex=1,4
        if(matrixindex==1 .OR. matrixindex==2) then
            chi1 = chi1array(j,Lambda)
            chi2 = chi2array(j,Lambda)
        else if (matrixindex==3 .OR. matrixindex==4) then
            chi1 = chi1array(j,LambdabarArray(Lambda))
            chi2 = chi2array(j,LambdabarArray(Lambda))
        end if
        if(spherical_only/=1 .or. (matrixindex == 1 .or. matrixindex==4)) then
            do eigenvalueindex=1,2
            nuL(Lambda,matrixindex,eigenvalueindex) = nuL(Lambda,matrixindex,eigenvalueindex) + &
                eigenvalue(matrixindex,eigenvalueindex) * (conjg(chi1) * eigenvector(matrixindex,eigenvalueindex,1) + conjg(chi2) * eigenvector(matrixindex,eigenvalueindex,2)) &
                * weight_array(j)
            end do
        end if
        end do ! matrixindex loop

        do matrixindex=1,4
        if(matrixindex==1 .OR. matrixindex==3) then
            chi1 = chi1array(j,Lambda)
            chi2 = chi2array(j,Lambda)
        else if (matrixindex==2 .OR. matrixindex==4) then
            chi1 = chi1array(j,LambdabarArray(Lambda))
            chi2 = chi2array(j,LambdabarArray(Lambda))
        end if
        if(spherical_only/=1 .or. (matrixindex == 1 .or. matrixindex==4)) then
            do eigenvalueindex=1,2
            nuR(Lambda,matrixindex,eigenvalueindex) = nuR(Lambda,matrixindex,eigenvalueindex) + &
                (conjg(eigenvector(matrixindex,eigenvalueindex,1)) * chi1 + conjg(eigenvector(matrixindex,eigenvalueindex,2)) * chi2 ) &
                * weight_array(j)
            end do
        end if
        end do ! matrixindex loop
    end do ! Lambda


  end do ! integrationpoints loop 

end subroutine nuCoefficients



subroutine ExpansionCoefficients(Lambda1,Lambda2,lcut,meshpoints,nr,chi1array,chi2array,theta_array,phi_array,weight_array,nuL,nuR,vcoeff)
! input: Lambda1,Lambda2,lcut,position r,SphiSpherical-Arrays chi1array,chi2array
! output v(1),v(2),v(3),v(4) expansion coefficients
  use SpinSphericals
  use Lebedev
  use DiracConfig
  implicit none
  complex (kind=dp), dimension(4)        :: vcoeff ! index 1..4 corresponds to a,b,c,d
  integer                :: Lambda1, Lambda2
  double precision                      :: meshpoints(:)
  integer                               :: nr

  double precision            :: r_fix
  double precision                      :: theta_array(integrationpoints), phi_array(integrationpoints), weight_array(integrationpoints)

  complex (kind=dp), allocatable        :: nuL(:,:,:), nuR(:,:,:)

  ! ^ index 1..4 corresponds to a,b,c,d i.e. the different sub-matrices, index i=1,2 corresponds to the different eigenvalues of the respective sub-matrix 
  integer                :: matrixindex
 
  integer                :: lcut, Lambdacut, printtimer
  integer, allocatable, dimension(:)    :: KappaArray, LambdabarArray
  real, allocatable, dimension(:)        :: MuArray
  complex (kind=dp), allocatable            :: chi1array(:,:), chi2array(:,:)
  printtimer = 1 ! 1=print computing step times 0=don't
  Lambdacut = getLambdacut(lcut)
  call makeLambdabarArray(lcut,LambdabarArray)
  call KappaMuArray(lcut,KappaArray,MuArray)

  r_fix = meshpoints(nr)

  vcoeff(:) = 0

  do matrixindex=1,4
      vcoeff(matrixindex) = nuL(Lambda1,matrixindex,1) * nuR(Lambda2,matrixindex,1) + nuL(Lambda1,matrixindex,2) * nuR(Lambda2,matrixindex,2)
  end do
end subroutine ExpansionCoefficients



subroutine SubMatrixEigenvalue(potPhi,potBx,potBy,potBz,eigenvalue)
  ! Input:
  ! eigenvalueindex = 1,2 corresponds to the two different eigenvalues of the corresponding sub-matrix
  ! matrixindex = 1,2,3,4 corresponds to a,b,c,d numbering the four different sub-matrices
  ! potPhi: scalar potential
  ! potBx,potBy,potBz: entries of the B-Field
  ! Output
  ! eigenvalue(matrixindex,eigenvalueindex)
  ! matrixindex = 1,2,3,4 corresponds to a,b,c,d numbering the four different sub-matrices
  ! eigenvalueindex = 1,2 corresponds to the two different eigenvalues of the corresponding sub-matrix
  use Constants
  implicit none
  double precision            :: potPhi,potBx,potBy,potBz
  complex (kind=dp), dimension(4,2)    :: eigenvalue
  eigenvalue(1,1) = echarge * potPhi + muB * sqrt(potBx**2 + potBy**2 + potBz**2)
  eigenvalue(1,2) = echarge * potPhi - muB * sqrt(potBx**2 + potBy**2 + potBz**2)
  eigenvalue(2,1) = 0
  eigenvalue(2,2) = 0
  eigenvalue(3,1) = 0
  eigenvalue(3,2) = 0
  eigenvalue(4,1) = echarge * potPhi + muB * sqrt(potBx**2 + potBy**2 + potBz**2)
  eigenvalue(4,2) = echarge * potPhi - muB * sqrt(potBx**2 + potBy**2 + potBz**2)
end subroutine SubMatrixEigenvalue


subroutine SubMatrixEigenvector(potPhi,potBx,potBy,potBz,eigenvector)
  ! Input:
  ! potPhi: scalar potential
  ! potBx,potBy,potBz: entries of the B-Field
  ! Output
  ! eigenvector(matrixindex,eigenvalueindex,entryindex)
  ! eigenvalueindex = 1,2 corresponds to the two different eigenvalues of the corresponding sub-matrix
  ! matrixindex = 1,2,3,4 corresponds to a,b,c,d numbering the four different sub-matrices
  ! ventryindex = 1,2 corresponds to the 1st and 2nd entry of the eigenvector
  use Constants
  implicit none
  double precision            :: potPhi,potBx,potBy,potBz
  complex (kind=dp), dimension(4,2,2)    :: eigenvector
  double precision                      :: norm
  integer                               :: matrixindex, eigenvalueindex
  if(potBx==0 .AND. potBy==0) then

    eigenvector(1,1,1) = 0
    eigenvector(1,1,2) = 1

    eigenvector(1,2,1) = 1
    eigenvector(1,2,2) = 0

    eigenvector(2,1,1) = 0
    eigenvector(2,1,2) = 1

    eigenvector(2,2,1) = 1
    eigenvector(2,2,2) = 0

    eigenvector(3,1,1) = 0
    eigenvector(3,1,2) = 1

    eigenvector(3,2,1) = 1
    eigenvector(3,2,2) = 0

    eigenvector(4,1,1) = 1
    eigenvector(4,1,2) = 0

    eigenvector(4,2,1) = 0
    eigenvector(4,2,2) = 1

else

    if(sqrt(potBx**2 + potBy**2 + potBz**2) + potBz == 0) then
      eigenvector(1,1,1) = 0
    else 
      eigenvector(1,1,1) = (-potBx + i * potBy) / ( sqrt(potBx**2 + potBy**2 + potBz**2) + potBz)
    end if
    eigenvector(1,1,2) = 1

    if(- sqrt(potBx**2 + potBy**2 + potBz**2) + potBz == 0) then
      eigenvector(1,2,1) = 0
    else
      eigenvector(1,2,1) = (-potBx + i * potBy) / (-sqrt(potBx**2 + potBy**2 + potBz**2) + potBz)
    end if
    eigenvector(1,2,2) = 1

    eigenvector(2,1,1) = 0
    eigenvector(2,1,2) = 1

    eigenvector(2,2,1) = 1
    eigenvector(2,2,2) = 0

    eigenvector(3,1,1) = 0
    eigenvector(3,1,2) = 1

    eigenvector(3,2,1) = 1
    eigenvector(3,2,2) = 0

    if(sqrt(potBx**2 + potBy**2 + potBz**2) - potBz == 0) then
      eigenvector(4,1,1) = 0
    else
      eigenvector(4,1,1) = ( potBx - i * potBy) / ( sqrt(potBx**2 + potBy**2 + potBz**2) - potBz)
    end if
    eigenvector(4,1,2) = 1

    if(-sqrt(potBx**2 + potBy**2 + potBz**2) - potBz == 0) then
      eigenvector(4,2,1) = 0
    else
      eigenvector(4,2,1) = ( potBx - i * potBy) / (-sqrt(potBx**2 + potBy**2 + potBz**2) - potBz)
    end if
    eigenvector(4,2,2) = 1
  end if

  ! These eigenvectors are orthogonal, but not yet normalised
  ! normalisation:
  do matrixindex=1,4
  do eigenvalueindex = 1,2
    norm = 1.0d0
    if(eigenvector(matrixindex,eigenvalueindex,1) /= 0 .OR. eigenvector(matrixindex,eigenvalueindex,2) /= 0) then
      norm = sqrt(abs(eigenvector(matrixindex,eigenvalueindex,1))**2.0d0 + abs(eigenvector(matrixindex,eigenvalueindex,2))**2.0d0)
      eigenvector(matrixindex,eigenvalueindex,1) = eigenvector(matrixindex,eigenvalueindex,1) / norm
      eigenvector(matrixindex,eigenvalueindex,2) = eigenvector(matrixindex,eigenvalueindex,2) / norm
    end if
  end do
  end do
end subroutine SubMatrixEigenvector



subroutine getPotPhi(zatom,meshpoints,nr,theta,phi,VLLin,potPhi) 
! Sample potential
  use SpinSphericals
  use Constants
  implicit none
  double precision :: r,theta,phi,potPhi,V_spinup,V_spindown
  double precision :: zatom,meshpoints(:)
  integer          :: nr,lcut,l,m
  double precision :: VLLin(:,:,:) !VLLin(nrmax,(2*lcut+1)**2,2)
  integer          :: lmmax,lm,nrmax

  r = meshpoints(nr)
  lmmax=size(VLLin(1,:,1),1)
  nrmax=size(VLLin(:,1,1),1)
  lcut = int(1.0d0/2.0d0 * (sqrt(real(lmmax)) -1)) ! lmmax=(2*lcut+1)**2
  V_spindown = 0.0d0
  lm=1
  do l=0,2*lcut,1
  do m=-l,l,1
      if (lm==1) then
          ! Due to a strange convention the V00 component has a prefactor that the values for other l,m indices don't have
          V_spindown = V_spindown + sqrt(4.0d0*Pi) * VLLin(nr,lm,1) * RealSphericalHarmonic(l,m,theta,phi)
      else
          V_spindown = V_spindown +                  VLLin(nr,lm,1) * RealSphericalHarmonic(l,m,theta,phi)
      end if
      lm = lm+1
  end do
  end do
  V_spinup = 0.0d0
  lm=1
  do l=0,lcut,1
  do m=-l,l,1
      if (lm==1) then
          ! Due to a strange convention the V00 component has a prefactor that the values for other l,m indices don't have
          V_spinup = V_spinup + sqrt(4.0d0*Pi) * VLLin(nr,lm,2) * RealSphericalHarmonic(l,m,theta,phi)
      else
          V_spinup = V_spinup +                  VLLin(nr,lm,2) * RealSphericalHarmonic(l,m,theta,phi)
      end if
      lm = lm+1
  end do
  end do
  ! account for different convention in the Dirac solver:
  V_spindown = V_spindown / echarge 
  V_spinup   = V_spinup   / echarge
  ! calculate total potential:
  potPhi = - echarge / 4.0d0 / Pi / eps0 / r * zatom + 0.5d0 * (V_spindown + V_spinup)
end subroutine getPotPhi 

subroutine getPotB(zatom,meshpoints,nr,theta,phi,VLLin,potBx,potBy,potBz)
! Sample potential
  use SpinSphericals
  use Constants
  implicit none
  double precision                :: r,theta,phi,potBx,potBy,potBz,V_spinup,V_spindown
  double precision                :: zatom,meshpoints(:)
  integer                         :: nr,lcut,l,m
  double precision                :: VLLin(:,:,:)
  integer                         :: lmmax,lm,nrmax

  r = meshpoints(nr)
  lmmax=size(VLLin(1,:,1),1)
  nrmax=size(VLLin(:,1,1),1)
  lcut = int(1/2 * (sqrt(real(lmmax)) -1)) ! lmmax=(2*lcut+1)**2

  V_spindown = 0.0d0
  lm=1
  do l=0,lcut,1
  do m=-l,l,1
      if (lm==1) then
          ! Due to a strange convention the V00 component has a prefactor that the values for other l,m indices don't have
          V_spindown = V_spindown + sqrt(4.0d0*Pi) * VLLin(nr,lm,1) * RealSphericalHarmonic(l,m,theta,phi)
      else
          V_spindown = V_spindown +                  VLLin(nr,lm,1) * RealSphericalHarmonic(l,m,theta,phi)
      end if
      lm = lm+1
  end do
  end do
  V_spinup = 0.0d0
  lm=1
  do l=0,2*lcut,1
  do m=-l,l,1
      if (lm==1) then
          ! Due to a strange convention the V00 component has a prefactor that the values for other l,m indices don't have
          V_spinup = V_spinup + sqrt(4.0d0*Pi) * VLLin(nr,lm,2) * RealSphericalHarmonic(l,m,theta,phi)
      else
          V_spinup = V_spinup +                  VLLin(nr,lm,2) * RealSphericalHarmonic(l,m,theta,phi)
      end if
      lm = lm+1
  end do
  end do
  ! account for different convention in the Dirac solver:
  V_spindown = V_spindown / echarge 
  V_spinup   = V_spinup   / echarge

  potBx = 0
  potBy = 0
  potBz = 0.5d0 * (- V_spindown + V_spinup)
end subroutine getPotB 

integer function CheckExistence(filename)
  character(len=100), intent(in)  :: filename
  integer                         :: status
  open(unit=78, file=filename, status='old', action='read', iostat=status)
  if(status == 0) then
      ! file exists
      CheckExistence = 1
  else
      ! file does not exist
      CheckExistence = 0
  end if
  close(unit=78)
end function CheckExistence

complex (kind=dp) function GetNumberOfLines(filename)
! opens a file and finds out the number of lines
  use DiracConfig
  implicit none
  character(len=100), intent(in)  :: filename
  integer                         :: nvals
  integer                         :: status
  double precision                :: value

  nvals=0
  open(unit=77, file=filename, status='old', action='read', iostat=status)
  if(status == 0) then
      if(verbose > 1) then
          write(*,*) "successfully opened ", filename
      end if
      do
          read(77,*,iostat=status) value
          if(status /= 0) exit
          nvals = nvals + 1
      end do
      if(status > 0) then
          write(*,*) "an error occured while reading from ", filename
      else
          if(verbose > 1) then
              write(*,*) "the file ", filename, " has ", nvals, " lines "
          end if
      end if
  else
      write(*,*) "ERROR in GetNumberOfLines():"
      write(*,*) "could not open file ", filename
  end if
  GetNumberOfLines = nvals

  close(unit=77)
end function GetNumberOfLines

end module Potential
