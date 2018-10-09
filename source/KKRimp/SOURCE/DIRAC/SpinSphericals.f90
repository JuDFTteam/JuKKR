module SpinSphericals

  use Constants

contains

  subroutine SpinSphericalHarmonic(kappa,mu,theta,phi,chi1,chi2)
    ! 1st component of the spin spherical harmonic chi for given values of kappa and mu
    ! at a given coordinate specified by theta and phi
    ! input: kappa,mu,theta,phi
    ! output: chi1 (first component), chi2 (second component)
    ! requires AssociatedLegendrePolynomials.f90
    implicit none
    integer :: kappa,l
    real :: mu,ms
    double precision :: theta,phi
    double complex :: chi1,chi2

    l=getL(kappa)
    ! first component
    ms = -.5e0
    if(ClebshGordan(kappa,mu,ms)/=0) then
      chi1 = dcmplx(ClebshGordan(kappa,mu,ms)) * ComplexSphericalHarmonic(l,int(mu-ms),theta,phi)
     else
       chi1 = 0
     end if
    ! (mu-ms should be an integer, the int function converts it to the correct data type)
    ! second component
    ms = .5e0
    if(ClebshGordan(kappa,mu,ms)/=0) then
      chi2 = dcmplx(ClebshGordan(kappa,mu,ms)) * ComplexSphericalHarmonic(l,int(mu-ms),theta,phi)
    else
      chi2 = 0
     end if
  end subroutine SpinSphericalHarmonic


  double complex function ComplexSphericalHarmonic(l,m,theta,phi)
    ! calculates the value of the Complex Spherical Harmonics
    ! input l=0..10, m=0..l
    implicit none
    integer, intent(in) :: l,m
    double precision, intent(in) :: theta, phi
  !  ComplexSphericalHarmonic = AssociatedLegendrePolynomial(l,abs(m),cos(theta))
    ComplexSphericalHarmonic = sqrt(((2.0d0*l+1.0d0)*factorial(l-abs(m)))/(4*Pi*factorial(l+abs(m)))) * exp(i*m*phi) * AssociatedLegendrePolynomial(l,m,cos(theta))
  end function ComplexSphericalHarmonic


  double complex function RealSphericalHarmonic(l,m,theta,phi)
    ! calculates the value of the Real Spherical Harmonics
    ! input l=0..10, m=0..l
    implicit none
    integer, intent(in) :: l,m
    double precision, intent(in) :: theta, phi
    if(m == 0) then
        RealSphericalHarmonic = ComplexSphericalHarmonic(l,m,theta,phi)
    else if (m>0) then
        RealSphericalHarmonic = 1.0d0/sqrt(2.0d0) *(  ComplexSphericalHarmonic(l, m,theta,phi)+ComplexSphericalHarmonic(l,-m,theta,phi) )
    else if (m<0) then
        RealSphericalHarmonic = i/sqrt(2.0d0)     *(  ComplexSphericalHarmonic(l,-m,theta,phi)-ComplexSphericalHarmonic(l, m,theta,phi) )
    end if
  end function RealSphericalHarmonic


  double precision function AssociatedLegendrePolynomial(l,mm,x)
    ! calculates the value of the Associated Legendre Polynomials
    ! input l=0..8, m=-l..l
    ! the expressions for the polynomials are stored in an external file
    ! that is called AssociatedLegendrePolynomialsExplicitely.f90
    implicit none
    integer, intent(in) :: l,mm
    integer                       :: m
    double precision, intent(in) :: x
    double precision :: func_value, prefactor
    if(mm>=0) then
      prefactor = 1.0d0
    else
      prefactor = (-1.0d0)**mm 
    end if
    m = abs(mm)
    include 'AssociatedLegendrePolynomialsExplicitely.f90'
    AssociatedLegendrePolynomial = prefactor * func_value
  end function AssociatedLegendrePolynomial


  integer function factorial(n)
    implicit none
    integer, intent(in) :: n
    integer :: temp, j
    temp = 1
    if(n==0) then
      factorial = 1
    else
      do j=1,n
         temp = temp*j
      end do
    factorial = temp
    end if
  end function factorial


  double precision function ClebshGordan(kappa,mu,ms)
    ! computes the Clebsh-Gordan coefficients
    implicit none
    integer, intent(in) :: kappa
    real, intent(in) :: mu, ms
    integer :: l
    real :: j
    j = getJ(kappa)
    l = getL(kappa)
    if(j==l+.5e0 .AND. ms==.5e0) then
      ClebshGordan = sqrt((l+mu+.5d0)/(2.0d0*l+1.0d0))
    else if(j==l+.5e0 .AND. ms==-.5e0) then
      ClebshGordan = sqrt((l-mu+.5d0)/(2.0d0*l+1.0d0))
    else if(j==l-.5e0 .AND. ms==.5e0) then
      ClebshGordan = -sqrt((l-mu+.5d0)/(2.0d0*l+1.0d0))
    else if(j==l-.5e0 .AND. ms==-.5e0) then
      ClebshGordan = sqrt((l+mu+.5d0)/(2.0d0*l+1.0d0))
    else
      write(*,*) 'Error in module SpinSpherical, function ClebshGordan: incorrect j value.'
      stop
    end if  
  end function ClebshGordan


  integer function getL(kappa)
    ! computes l from a given value of kappa
    ! see e.g. E.M. Rose, Rel. Elektronentheorie I, S. 43
    implicit none
      integer, intent(in) :: kappa
      if (kappa >= 0) then
         getL = kappa
      else
         getL = -kappa-1
      end if
  end function getL


  integer function getLbar(kappa)
    ! computes l_bar from a given value of kappa
    ! see e.g. E.M. Rose, Rel. Elektronentheorie I, S. 44
    implicit none
      integer, intent(in) :: kappa
      if (kappa > 0) then
         getLbar = kappa-1
      else
         getLbar = -kappa
      end if
  end function getLbar


  real function getJ(kappa)
    ! computes j from a given value of kappa
    ! see e.g. E.M. Rose, Rel. Elektronentheorie I, S. 43
    implicit none
      integer, intent(in) :: kappa
      getJ = Abs(kappa)-.5d0
  end function getJ


  subroutine makeLambdabarArray(lcut,LambdabarArray)
    implicit none
    ! input variable
    integer   :: lcut
    ! outupt
    integer,allocatable,dimension(:) :: LambdabarArray
    ! other variables
    integer                               :: Lambda,Lambdacut
    integer,allocatable,dimension(:)      :: KappaArray
    real,allocatable,dimension(:)         :: MuArray
    
    Lambdacut = getLambdacut(lcut)
    allocate(LambdabarArray(Lambdacut))
    call KappaMuArray(lcut,KappaArray,MuArray)

    do Lambda=1,Lambdacut
      LambdabarArray(Lambda) = Lambda-2*KappaArray(Lambda)
    end do
  end subroutine makeLambdabarArray
    

  subroutine KappaMuArray(lcut,KappaArray,MuArray)
    implicit none
    ! input variable
    integer :: lcut
    ! output variables
    integer,allocatable,dimension(:) :: KappaArray
    real,allocatable,dimension(:) :: MuArray
    ! other variables
    integer :: k,kk
    integer :: absolutkappa,kappa
    real :: j, mu
    allocate(KappaArray(getLambdacut(lcut+2)))
    allocate(MuArray(getLambdacut(lcut+2)))
    k=1
    do absolutkappa=1,lcut+2
      j = absolutkappa - 0.5e0
      ! case j=l+1/2
      kappa = - absolutkappa
      do kk=1,absolutkappa*2
        mu = -j + (real(kk)-1.0e0)
        KappaArray(k) = kappa
        MuArray(k) = mu
        k = k+1
      end do

      ! case j=l-1/2
      if(absolutkappa /= lcut+1) then
        kappa = absolutkappa
        do kk=1,absolutkappa*2
          mu = -j + (real(kk)-1.0e0)
          KappaArray(k) = kappa
          MuArray(k) = mu
          k = k+1
        end do
      end if

    end do
  end subroutine KappaMuArray


  integer function getLambdacut(lcut)
    implicit none
    integer, intent(in) :: lcut
    getLambdacut=2*(lcut+1)**2
  end function getLambdacut


  subroutine Cartesian2Spherical(x,y,z,r,theta,phi)
    implicit none
    double precision :: x,y,z
    double precision :: r,theta,phi
    r = sqrt(x**2+y**2+z**2)
    if(x>0) then
      phi = atan(y/x)
    else if (x==0) then
      phi = sign(1.0d0,y)*Pi/2
    else if (x<0 .AND. y>=0) then
      phi = atan(y/x)+Pi
    else ! x<0 .AND. y<0
      phi = atan(y/x)-Pi
    end if
    if(r>0) then
      theta = acos(z/r)
    else
      theta = 0
    end if
  end subroutine Cartesian2Spherical


  subroutine Spherical2Cartesian(x,y,z,r,theta,phi)
    implicit none
    double precision :: x,y,z
    double precision :: r,theta,phi
    x=r*sin(theta)*cos(phi)
    y=r*sin(theta)*sin(phi)
    z=r*cos(theta)
  end subroutine Spherical2Cartesian


end module SpinSphericals

