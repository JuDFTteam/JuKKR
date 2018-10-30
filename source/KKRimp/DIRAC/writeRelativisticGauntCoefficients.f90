program writeRelativisticGauntCoefficients
use SpinSphericals
use Lebedev
use RelativisticGauntCoefficients
implicit none
integer                         :: lcut
write(*,*) "This programme computes the Relativistic Gaunt Coefficients and writes them into a file (new version)."
write(*,*) "Enter l-cutoff value:"
read(*,*) lcut
call writeDcoeff(lcut)
end program writeRelativisticGauntCoefficients