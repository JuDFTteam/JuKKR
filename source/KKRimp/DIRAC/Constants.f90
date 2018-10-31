module Constants

double precision, parameter     :: Pi = 3.14159265358979323846d0
complex, parameter              :: i  = (0.0d0,1.0d0)

! ! SI UNITS
! ! electron charge
! double precision, parameter     :: echarge = 1.60217648740d-19 ! C
! ! reduced Planck constant
! double precision, parameter     :: hbar = 1.05457162853d-34 ! Js
! ! electron rest mass
! double precision, parameter     :: emass = 9.1093821545d-31 ! kg
! ! speed of light
! double precision, parameter     :: lightspeed = 299792458d0 ! m/s
! ! Bohr magneton
! double precision, parameter     :: muB = echarge * hbar / 2.0d0 / emass / lightspeed ! = 9.2740091523*10**(-24) J/T #
! ! vacuum permittivity
! double precision, parameter     :: eps0 = 8.854187817d-12 ! F·m−1

! ! ATOMIC HARTREE UNITS
! ! electron charge
! double precision, parameter     :: echarge    = 1.0d0
! ! reduced Planck constant
! double precision, parameter     :: hbar       = 1.0d0
! ! electron rest mass
! double precision, parameter     :: emass      = 1.0d0
! ! speed of light
! double precision, parameter     :: lightspeed = 137.035999074d0 ! 2010 CODATA value
! ! Bohr magneton
! double precision, parameter     :: muB        = echarge * hbar / 2.0d0 / emass / lightspeed
! ! vacuum permittivity
! double precision, parameter     :: eps0       = 1.0d0 / 4.0d0 / Pi

! ATOMIC RYDBERG UNITS
! electron charge
double precision, parameter     :: echarge    = sqrt(2.0d0)
! reduced Planck constant
double precision, parameter     :: hbar       = 1.0d0
! electron rest mass
double precision, parameter     :: emass      = 0.5d0
! speed of light
double precision, parameter     :: lightspeed = 274.0720442d0 !echarge**2 / hbar * 137.035999074d0 ! =2*alpha=274...
! Bohr magneton
double precision, parameter     :: muB        = echarge * hbar / 2.0d0 / emass / lightspeed
! vacuum permittivity
double precision, parameter     :: eps0       = 1.0d0 / 4.0d0 / Pi

end module Constants