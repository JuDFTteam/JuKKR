c from the specfun package as included in scipy

        subroutine sphj(n,x,nm,sj,dj)
c       modified to allow n=0 case (also in csphjy, sphy)
c
c       =======================================================
c       purpose: compute spherical bessel functions jn(x) and
c                their derivatives
c       input :  x --- argument of jn(x)
c                n --- order of jn(x)  ( n = 0,1,… )
c       output:  sj(n) --- jn(x)
c                dj(n) --- jn'(x)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the starting
c                point for backward recurrence
c       =======================================================
c
        implicit none ! double precision (a-h,o-z)
        integer, intent(in) :: n
        integer, intent(out) :: nm
        double precision, intent(in) :: x
        double precision, intent(out) :: sj(0:n), dj(0:n)
        
        integer, external :: msta1, msta2
        double precision :: f, f0, f1, cs, sa, sb
        integer :: k, m
        nm=n
        if (dabs(x).lt.1.0d-100) then
           do 10 k=0,n
              sj(k)=0.0d0
10            dj(k)=0.0d0
           sj(0)=1.0d0
           if (n.gt.0) dj(1)=.3333333333333333d0
           return
        endif
        sj(0)=dsin(x)/x
        dj(0)=(dcos(x)-dsin(x)/x)/x
        if (n.lt.1) return
        sj(1)=(sj(0)-dcos(x))/x
        if (n.ge.2) then
           sa=sj(0)
           sb=sj(1)
           m=msta1(x,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(x,n,15)
           endif
           f=0.0d0
           f0=0.0d0
           f1=1.0d0-100
           do 15 k=m,0,-1
              f=(2.0d0*k+3.0d0)*f1/x-f0
              if (k.le.nm) sj(k)=f
              f0=f1
15            f1=f
           cs=0.0d0
           if (dabs(sa).gt.dabs(sb)) cs=sa/f
           if (dabs(sa).le.dabs(sb)) cs=sb/f0
           do 20 k=0,nm
20            sj(k)=cs*sj(k)
        endif
        do 25 k=1,nm
25         dj(k)=sj(k-1)-(k+1.0d0)*sj(k)/x
        endsubroutine

        integer function msta1(x,mp)
c       ===================================================
c       purpose: determine the starting point for backward
c                recurrence such that the magnitude of
c                jn(x) at that point is about 10^(-mp)
c       input :  x     --- argument of jn(x)
c                mp    --- value of magnitude
c       output:  msta1 --- starting point
c       ===================================================
        implicit none ! double precision (a-h,o-z)
        double precision, intent(in) :: x
        integer, intent(in) :: mp
        double precision, external :: envj
        double precision :: a0, f0, f1, f
        integer :: n0, n1, nn, it
        
        a0=dabs(x)
        n0=int(1.1d0*a0)+1
        f0=envj(n0,a0)-mp
        n1=n0+5
        f1=envj(n1,a0)-mp
        do 10 it=1,20
           nn=n1-(n1-n0)/(1.0d0-f0/f1)
           f=envj(nn,a0)-mp
           if(abs(nn-n1).lt.1) go to 20
           n0=n1
           f0=f1
           n1=nn
 10        f1=f
 20     msta1=nn
        endfunction


        integer function msta2(x,n,mp)
c       ===================================================
c       purpose: determine the starting point for backward
c                recurrence such that all jn(x) has mp
c                significant digits
c       input :  x  --- argument of jn(x)
c                n  --- order of jn(x)
c                mp --- significant digit
c       output:  msta2 --- starting point
c       ===================================================
        implicit none ! double precision (a-h,o-z)
        double precision, intent(in) :: x
        integer, intent(in) :: n, mp
        double precision, external :: envj
        double precision :: a0, hmp, ejn, obj, f0, f1, f
        integer :: n0, n1, nn, it
        a0=dabs(x)
        hmp=0.5d0*mp
        ejn=envj(n,a0)
        if (ejn.le.hmp) then
           obj=mp
           n0=int(1.1*a0)+1
        else
           obj=hmp+ejn
           n0=n
        endif
        f0=envj(n0,a0)-obj
        n1=n0+5
        f1=envj(n1,a0)-obj
        do 10 it=1,20
           nn=n1-(n1-n0)/(1.0d0-f0/f1)
           f=envj(nn,a0)-obj
           if (abs(nn-n1).lt.1) go to 20
           n0=n1
           f0=f1
           n1=nn
10         f1=f
20      msta2=nn+10
        endfunction

        real*8 function envj(n,x)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: x
        envj=0.5d0*dlog10(6.28d0*n)-n*dlog10(1.36d0*x/n)
        endfunction

c       **********************************
