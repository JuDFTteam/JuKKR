module RelativisticGauntCoefficients

contains

subroutine writeDcoeff(lcut)
use SpinSphericals
use Lebedev
implicit none

integer                         :: lcut

character(len=40)               :: filename
character(len=2)                :: lcutstring

double complex,allocatable      :: Dcoeff(:,:,:,:)
double complex,allocatable      :: chi(:,:)
integer                         :: Lambda1,Lambda2,Lambda3,Lambda4,Lambdacut,Lambdacut_total,kappa
real                            :: mu

integer,allocatable             :: KappaArray(:), LambdabarArray(:)
real,allocatable                :: MuArray(:)


double complex                  :: func_value, integral_contribution
integer                         :: j, totalentries, nonzeroentriesA, nonzeroentriesB
double precision                :: weight,x,y,z,r,theta,phi
double precision                :: eps

eps=1d-16
write(*,*) "values smaller than ",eps," are set to 0"

call KappaMuArray(lcut+1,KappaArray,MuArray)
call makeLambdabarArray(lcut+1,LambdabarArray)
Lambdacut=getLambdacut(lcut)
Lambdacut_total=LambdabarArray(getLambdacut(lcut))
write(*,*) "Lambdacut=",Lambdacut

allocate(chi(Lambdacut_total,2))
allocate(Dcoeff(Lambdacut_total,Lambdacut_total,Lambdacut_total,Lambdacut_total))
 chi(:,:) = 0.0d0
 Dcoeff(:,:,:,:) = 0.0d0

write(lcutstring,'(I02.2)') lcut
filename = 'RelativisticGauntCoeffA_lcut'//trim(lcutstring)//'.txt'
open(unit=1,file=filename)
write(*,*) "Storing coefficients into file ", filename
filename = 'RelativisticGauntCoeffB_lcut'//trim(lcutstring)//'.txt'
open(unit=2,file=filename)
write(*,*) "Storing coefficients for Lambda-bar into file ", filename

write(*,*) "progress:"
do j=1,integrationpoints
    write(*,*) j,"/",integrationpoints
    call lebedevpoint(j,x,y,z,weight)
    call Cartesian2Spherical(x,y,z,r,theta,phi)

    do Lambda1=1,Lambdacut_total
        ! calculate values of the spin spherical harmonics
        kappa=KappaArray(Lambda1)
        mu=MuArray(Lambda1)
        call SpinSphericalHarmonic(kappa,mu,theta,phi,chi(Lambda1,1),chi(Lambda1,2))
    end do

    do Lambda1=1,Lambdacut_total
    do Lambda2=1,Lambdacut_total
    do Lambda3=1,Lambdacut_total
    do Lambda4=1,Lambdacut_total
        func_value =  (CONJG(chi(Lambda1,1))*chi(Lambda2,1) + CONJG(chi(Lambda1,2))*chi(Lambda2,2)) &
                    * (CONJG(chi(Lambda3,1))*chi(Lambda4,1) + CONJG(chi(Lambda3,2))*chi(Lambda4,2))
        integral_contribution = func_value*weight
        Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4) = Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4) + integral_contribution
    end do
    end do
    end do
    end do
end do


totalentries = Lambdacut**4
nonzeroentriesA = 0
nonzeroentriesB = 0
do Lambda1=1,Lambdacut
do Lambda2=1,Lambdacut
do Lambda3=1,Lambdacut
do Lambda4=1,Lambdacut
    ! file A for Lambda-coefficients
    if(DBLE(Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4)) > eps .OR. DIMAG(Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4)) > eps) then
        nonzeroentriesA = nonzeroentriesA + 1
        if(DBLE(Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4)) < eps) then
            Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4) = DIMAG(Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4)) 
        elseif(DIMAG(Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4)) < eps) then
            Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4) = DBLE(Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4)) 
        end if
        write(1,'(I3.3, " ",  I3.3, " ", I3.3, " ", I3.3, " ", E23.16, " ", E23.16)') Lambda1, Lambda2, Lambda3, Lambda4, DBLE(Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4)), DIMAG(Dcoeff(Lambda1,Lambda2,Lambda3,Lambda4))
    end if
    ! file B for (overlined) Lambda-bar-coefficients
    if(DBLE(Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4))) > eps .OR. DIMAG(Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4))) > eps) then
        nonzeroentriesB = nonzeroentriesB + 1
        if(DBLE(Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4))) < eps) then
            Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4)) = DIMAG(Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4))) 
        elseif(DIMAG(Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4))) < eps) then
            Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4)) = DBLE(Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4))) 
        end if
        write(2,'(I3.3, " ",  I3.3, " ", I3.3, " ", I3.3, " ", E23.16, " ", E23.16)') Lambda1, Lambda2, Lambda3, Lambda4, DBLE(Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4))), DIMAG(Dcoeff(LambdabarArray(Lambda1),LambdabarArray(Lambda2),LambdabarArray(Lambda3),LambdabarArray(Lambda4)))
    end if

end do
end do
end do
end do

close(1)
close(2)

write(*,*) "Out of ", totalentries, " entries ", nonzeroentriesA, " are non-zero for file A"
write(*,*) "(i.e. a rate of ", DBLE(nonzeroentriesA)/DBLE(totalentries), ")"
write(*,*) "and ", nonzeroentriesB, " are non-zero for file B"
write(*,*) "(i.e. a rate of ", DBLE(nonzeroentriesB)/DBLE(totalentries), ")"
write(*,*) ""
write(*,*) "done."

end subroutine writeDcoeff

end module RelativisticGauntCoefficients