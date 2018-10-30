! Please note that for l>6 calculated values become very inaccurate
! and there are large deviations from the orthonormality relation!
if (l==0 .AND. m==0) then 
 func_value = 1.0d0
else if (l==1 .AND. m==0) then 
 func_value = x
else if (l==1 .AND. m==1) then 
 func_value = -(1.0d0-x**2.0d0)**(1.0d0/2.0d0)
else if (l==2 .AND. m==0) then 
 func_value = 3d0/2d0*x**2d0-1d0/2d0
else if (l==2 .AND. m==1) then 
 func_value = -3d0*(1d0-x**2d0)**(1d0/2d0)*x
else if (l==2 .AND. m==2) then 
 func_value = 3d0-3d0*x**2d0
else if (l==3 .AND. m==0) then 
 func_value = x**3d0+3d0/2d0*(x**2d0-1d0)*x
else if (l==3 .AND. m==1) then 
 func_value = -1d0/48d0*(1d0-x**2d0)**(1d0/2d0)*(360d0*x**2d0-72d0)
else if (l==3 .AND. m==2) then 
 func_value = 15d0*(1d0-x**2d0)*x
else if (l==3 .AND. m==3) then 
 func_value = -15d0*(1d0-x**2d0)**(3d0/2d0)
else if (l==4 .AND. m==0) then 
 func_value = x**4d0+3d0*(x**2d0-1d0)*x**2d0+3d0/8d0*(x**2d0-1d0)**2d0
else if (l==4 .AND. m==1) then 
 func_value = -1d0/384d0*(1d0-x**2d0)**(1d0/2d0)*(3840d0*x**3d0+2880d0*(x**2d0-1d0)*x)
else if (l==4 .AND. m==2) then 
 func_value = 1d0/384d0*(1d0-x**2d0)*(20160d0*x**2d0-2880d0)
else if (l==4 .AND. m==3) then 
 func_value = -105d0*(1d0-x**2d0)**(3d0/2d0)*x
else if (l==4 .AND. m==4) then 
 func_value = 105d0*(1d0-x**2d0)**2d0
else if (l==5 .AND. m==0) then 
 func_value = x**5d0+5d0*(x**2d0-1d0)*x**3d0+15d0/8d0*(x**2d0-1d0)**2d0*x
else if (l==5 .AND. m==1) then 
 func_value = -1d0/3840d0*(1d0-x**2d0)**(1d0/2d0)*(57600d0*x**4d0+86400d0*(x**2d0-1d0)*x**2d0+7200d0*(x**2d0-1d0)**2d0)
else if (l==5 .AND. m==2) then 
 func_value = 1d0/3840d0*(1d0-x**2d0)*(403200d0*x**3d0+201600d0*(x**2d0-1d0)*x)
else if (l==5 .AND. m==3) then 
 func_value = -1d0/3840d0*(1d0-x**2d0)**(3d0/2d0)*(1814400d0*x**2d0-201600d0)
else if (l==5 .AND. m==4) then 
 func_value = 945d0*(1d0-x**2d0)**2d0*x
else if (l==5 .AND. m==5) then 
 func_value = -945d0*(1d0-x**2d0)**(5d0/2d0)
else if (l==6 .AND. m==0) then 
 func_value = x**6d0+15d0/2d0*(x**2d0-1d0)*x**4d0+45d0/8d0*(x**2d0-1d0)**2d0*x**2d0+5d0/16d0*(x**2d0-1d0)**3d0
else if (l==6 .AND. m==1) then 
 func_value = -1d0/46080d0*(1d0-x**2d0)**(1d0/2d0)*(967680d0*x**5d0+2419200d0*(x**2d0-1d0)*x**3d0+604800d0*(x**2d0-1d0)**2d0*x)
else if (l==6 .AND. m==2) then 
 func_value = 1d0/46080d0*(1d0-x**2d0)*(9676800d0*x**4d0+9676800d0*(x**2d0-1d0)*x**2d0+604800d0*(x**2d0-1d0)**2d0)
else if (l==6 .AND. m==3) then 
 func_value = -1d0/46080d0*(1d0-x**2d0)**(3d0/2d0)*(58060800d0*x**3d0+21772800d0*(x**2d0-1d0)*x)
else if (l==6 .AND. m==4) then 
 func_value = 1d0/46080d0*(1d0-x**2d0)**2d0*(239500800d0*x**2d0-21772800d0)
else if (l==6 .AND. m==5) then 
 func_value = -10395d0*(1d0-x**2d0)**(5d0/2d0)*x
else if (l==6 .AND. m==6) then 
 func_value = 10395d0*(1d0-x**2d0)**3d0
else if (l==7 .AND. m==0) then 
 func_value = x**7d0+21d0/2d0*(x**2d0-1d0)*x**5d0+105d0/8d0*(x**2d0-1d0)**2d0*x**3d0+35d0/16d0*(x**2d0-1d0)**3d0*x
else if (l==7 .AND. m==1) then 
 func_value = -1d0/645120d0*(1d0-x**2d0)**(1d0/2d0)*(18063360d0*x**6d0+67737600d0*(x**2d0-1d0)*x**4d0+33868800d0*(x**2d0-1d0)**2d0*x**2d0+1411200d0*(x**2d0-1d0)**3d0)
else if (l==7 .AND. m==2) then 
 func_value = 1d0/645120d0*(1d0-x**2d0)*(243855360d0*x**5d0+406425600d0*(x**2d0-1d0)*x**3d0+76204800d0*(x**2d0-1d0)**2d0*x)
else if (l==7 .AND. m==3) then 
 func_value = -1d0/645120d0*(1d0-x**2d0)**(3d0/2d0)*(2032128000d0*x**4d0+1524096000d0*(x**2d0-1d0)*x**2d0+76204800d0*(x**2d0-1d0)**2d0)
else if (l==7 .AND. m==4) then 
 func_value = 1d0/645120d0*(1d0-x**2d0)**2d0*(11176704000d0*x**3d0+3353011200d0*(x**2d0-1d0)*x)
else if (l==7 .AND. m==5) then 
 func_value = -1d0/645120d0*(1d0-x**2d0)**(5d0/2d0)*(43589145600d0*x**2d0-3353011200d0)
else if (l==7 .AND. m==6) then 
 func_value = 135135d0*(1d0-x**2d0)**3d0*x
else if (l==7 .AND. m==7) then 
 func_value = -135135d0*(1d0-x**2d0)**(7d0/2d0)
else if (l==8 .AND. m==0) then 
 func_value = x**8d0+14d0*(x**2d0-1d0)*x**6d0+105d0/4d0*(x**2d0-1d0)**2d0*x**4d0+35d0/4d0*(x**2d0-1d0)**3d0*x**2d0+35d0/128d0*(x**2d0-1d0)**4d0
else if (l==8 .AND. m==1) then 
 func_value = -1d0/10321920d0*(1d0-x**2d0)**(1d0/2d0)*(371589120d0*x**7d0+1950842880d0*(x**2d0-1d0)*x**5d0+1625702400d0*(x**2d0-1d0)**2d0*x**3d0+203212800d0*(x**2d0-1d0)**3d0*x)
else if (l==8 .AND. m==2) then 
 func_value = 1d0/10321920d0*(1d0-x**2d0)*(6502809600d0*x**6d0+16257024000d0*(x**2d0-1d0)*x**4d0+6096384000d0*(x**2d0-1d0)**2d0*x**2d0+203212800d0*(x**2d0-1d0)**3d0)
else if (l==8 .AND. m==3) then 
 func_value = -1d0/10321920d0*(1d0-x**2d0)**(3d0/2d0)*(71530905600d0*x**5d0+89413632000d0*(x**2d0-1d0)*x**3d0+13412044800d0*(x**2d0-1d0)**2d0*x)
else if (l==8 .AND. m==4) then 
 func_value = 1d0/10321920d0*(1d0-x**2d0)**2d0*(536481792000d0*x**4d0+321889075200d0*(x**2d0-1d0)*x**2d0+13412044800d0*(x**2d0-1d0)**2d0)
else if (l==8 .AND. m==5) then 
 func_value = -1d0/10321920d0*(1d0-x**2d0)**(5d0/2d0)*(2789705318400d0*x**3d0+697426329600d0*(x**2d0-1d0)*x)
else if (l==8 .AND. m==6) then 
 func_value = 1d0/10321920d0*(1d0-x**2d0)**3d0*(10461394944000d0*x**2d0-697426329600d0)
else if (l==8 .AND. m==7) then 
 func_value = -2027025d0*(1d0-x**2d0)**(7d0/2d0)*x
else if (l==8 .AND. m==8) then 
 func_value = 2027025d0*(1d0-x**2d0)**4d0
else if (m>l) then 
 func_value = 0
else
write(*,*) 'Error in module SpinSpherical, function AssociatedLegendrePolynomials: incorrect l,m value. Choose not more than l=8.'
stop
end if
