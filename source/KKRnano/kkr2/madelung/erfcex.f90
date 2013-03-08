double precision function erfcex(z)
implicit none
!-----------------------------------------------------------------------
!
!     calculates complementary errorfunction times sqrt(pi)
!      times exp(z*z)  by continued fractions
!
!-----------------------------------------------------------------------
!     .. scalar arguments ..
double precision z
!     ..
!     .. local scalars ..
double precision bound,erf1,exzz,f,fa,q,ratio,sqrtpi,term,u,ua, &
                 v,x,xa,y,z2,zz
!     ..
!     .. intrinsic functions ..
intrinsic abs,atan,exp,sqrt
!     ..
bound=3.d-11
sqrtpi = sqrt(4.0d0*atan(1.0d0))
zz = z*z
!
!---> choose algorithm
!
if (z.lt.1.5d0) then

!     this exponential was outside the if statement
!     but for large arguments the exponent blow up
!     changes made 21/10/99

   exzz = exp(zz)

   z2 = 2.0d0*zz
   erf1 = z
   ratio = 1.0d0
   term = z
10    continue
   ratio = ratio + 2.0d0
   term = term*z2/ratio
   erf1 = erf1 + term
   if (term.gt.bound) go to 10
   erfcex = sqrtpi*exzz - 2.0d0*erf1

else
!
!---> continued fraction expansion : abramowitz p. 298, eq. (7.1.14)
!
   u = 1.0d0
   v = 0.0d0
   x = z
   y = 1.0d0
   q = 0.5d0
   f = (u+v*q)/ (x+y*q)
20    continue
   ua = u
   u = u*z + v*q
   v = ua
   xa = x
   x = x*z + y*q
   y = xa
   q = q + 0.5d0
   fa = f
   f = (u+v*q)/ (x+y*q)
   if (abs(fa-f).gt.bound*f) go to 20
   erfcex = f
end if

end
