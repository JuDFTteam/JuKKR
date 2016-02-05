      SUBROUTINE SPHERE(LMAX,YR,WTYR,RIJ,IJEND,IJD1)
      implicit none
c-----------------------------------------------------------------------
C     IJD=434
c     generate an angular mesh and spherical harmonics at those
c     mesh points. For an angular integration the weights are ge-
c     rated .
c
c     R. Zeller      Feb. 1996
c     Small change for GGA implementation
c     Nikos          Dec. 1996 
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      include 'inc.p'
      INTEGER IJEND,LMAX,IJD1
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI,R,R1,R2,R3
      INTEGER IJ,LM1,NP,N,NN
C     ..
C     .. External Subroutines ..
      EXTERNAL YMY,CYLM02
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*)
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION W(1000),Y(1000)
      DOUBLE PRECISION COSX(IJD),COSFI(IJD),SINFI(IJD)
      DOUBLE PRECISION FAI(IJD),f0,nd(3,3),dx1,dx2,dx3
      DOUBLE PRECISION ND2(3,3),f1
C
C     ..
C     .. Data statements ..
      DATA PI/3.14159265359D0/
C     ..
      PI=4.d0*DATAN(1.d0)
      write (6,*) 'SPHERE FOR ***GGA***  : read lebedev.mesh' 
c      READ (89,*) IJEND
      IF (IJD.GT.1000) STOP 'SPHERE'
c
c
      DO 30 IJ = 1,IJD
        CALL LEBEDEV(IJ,R1,R2,R3,W(IJ))
c        READ (89,*) R1,R2,R3,W(IJ)
c
c      make a small rotation
c 
          f0 = 0.08d0
c         f1 = 0.05d0
          ND(1,1) = COS(f0)
          ND(1,2) = 0
          ND(1,3) = SIN(f0)
          ND(2,1) = 0
          ND(2,2) = 1
          ND(2,3) = 0.0
          ND(3,1) = -SIN(f0)
          ND(3,2) = 0
          ND(3,3) = COS(f0)
c
c         ND2(1,1) = dcos(f1)
c         ND2(1,2) = -dsin(f1)
c         ND2(1,3) = 0
c         ND2(2,1) = dsin(f1)
c         ND2(2,2) = dcos(f1)
c         ND2(2,3) = 0.0
c         ND2(3,1) = 0
c         ND2(3,2) = 0
c         ND2(3,3) = 1.
c
          DX1=ND(1,1)*R1 + ND(2,1)*R2 +
     +          ND(3,1)*R3
          DX2=ND(1,2)*R1 + ND(2,2)*R2 +
     +          ND(3,2)*R3
          DX3=ND(1,3)*R1 + ND(2,3)*R2 +
     +          ND(3,3)*R3
c         R1 = DX1
c         R2 = DX2
c         R3 = DX3
c         DX1=ND2(1,1)*R1 + ND2(2,1)*R2 +
c    +          ND2(3,1)*R3
c         DX2=ND2(1,2)*R1 + ND2(2,2)*R2 +
c    +          ND2(3,2)*R3
c         DX3=ND2(1,3)*R1 + ND2(2,3)*R2 +
c    +          ND2(3,3)*R3
        R1 = DX1
        R2 = DX2
        R3 = DX3
           
        RIJ(IJ,1) = R1
        RIJ(IJ,2) = R2
        RIJ(IJ,3) = R3
        CALL YMY(R1,R2,R3,R,Y,LMAX)
        DO 10 LM1 = 1, (LMAX+1)**2
          YR(IJ,LM1) = Y(LM1)
   10   CONTINUE
c
c---> multiply the spherical harmonics with the weights
c
        DO 20 LM1 = 1, (LMAX+1)**2
          WTYR(IJ,LM1) = YR(IJ,LM1)*W(IJ)*PI*4.D0
   20   CONTINUE
c produce what is neaded for GGA
        COSX(IJ) =  R3
        IF (ABS(R3).NE.1.d0) THEN 
        COSFI(IJ) =R1/SQRT(1.d0-R3*R3)
        SINFI(IJ) =R2/SQRT(1.d0-R3*R3)
        IF (ABS(COSFI(IJ)).GT.1.d0) 
     &               COSFI(IJ)=COSFI(IJ)/ABS(COSFI(IJ))
        IF (ABS(SINFI(IJ)).GT.1.d0) 
     &               SINFI(IJ)=SINFI(IJ)/ABS(SINFI(IJ))
        FAI(IJ) = ACOS(COSFI(IJ))
          ELSE IF (SINFI(IJ).EQ.0.d0) THEN
             FAI(IJ) = pi/2.d0
          ELSE
        COSFI(IJ) =R1
        SINFI(IJ) =R2
        IF (ABS(COSFI(IJ)).GT.1.d0) 
     &                COSFI(IJ)=COSFI(IJ)/ABS(COSFI(IJ))
          END IF
        FAI(IJ) = ACOS(COSFI(IJ))
!       WRITE(6,*)IJ,COSX(IJ),COSFI(IJ),SINFI(IJ),FAI(IJ) 
   30 CONTINUE 

       NP = IJEND
      write(6,9000) IJEND,N,NN,LMAX
 9000 format(1x,' IJEND N NN LMAX',4I6)
 9001 format(4D20.12)
      CALL CYLM02(LMAX,COSX,FAI,SINFI,COSFI)

      END
cdeck
      SUBROUTINE CYLM02(lmax,COSX,FAI,SINFI,COSFI)
c.....------------------------------------------------------------------
c     preparation of cylm0(=ylm(ip,i)), cylmt1(=dylm/dtheta),
c     cylmt2(=d2ylm/dt2), 
c     cylmf1, cylmf2 are for fai.
c     cylmtf=d2ylm/dfdt
c     i=1,2,....,(lmax+1)**2                  
c.....------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      include 'inc.p'
      INTEGER N,LMMAXD,LMD0
      PARAMETER (N=2* (LPOTD+1),LMMAXD=(LPOTD+1)**2)
      PARAMETER (LMD0=2*LPOTD+1)

c.....------------------------------------------------------------------
      DOUBLE COMPLEX cylm0(LMMAXD)
      DOUBLE COMPLEX cylmt1(LMMAXD),cylmt2(LMMAXD)
      DOUBLE COMPLEX cylmf1(LMMAXD),cylmf2(LMMAXD),cylmtf(LMMAXD)
      DOUBLE COMPLEX ep1f,ep2f,em1f,em2f,ci
      DOUBLE PRECISION bb1(lmmaxd)
      INTEGER llmax,np,ip  
      DOUBLE PRECISION COSX(IJD),FAI(IJD),COSFI(IJD),SINFI(IJD)
c.....------------------------------------------------------------------
      common/cylmdnp/thet(IJD),ylm(IJD,LMMAXD),dylmt1(IJD,LMMAXD),
     &  dylmt2(IJD,LMMAXD),dylmf1(IJD,LMMAXD),dylmf2(IJD,LMMAXD),
     &  dylmtf(IJD,LMMAXD)
c.....------------------------------------------------------------------
      DOUBLE PRECISION yl(21)                                                  
c.....------------------------------------------------------------------
c                                                                       
      ci=CMPLX(0.d0,1.d0)
      one=1.d0
      pi=4.d0*DATAN(one)                                                
      llmax=(lmax+1)**2

        nph=ijd/2
        np = ijd
       write(6,*) ' in cylm02 np=',np
        
      do 11 ip=1,np                                                     

c
c
c     if(ip.gt.nph) cosx(ip)=-cosx(ip)
c     this is changed (new lebedev points)
c    
       thet(ip)=ACOS(cosx(ip))
       fi=fai(ip)
       di=2*fai(ip)
c      SIN2F = 2.d0*SINFI(IP)*COSFI(IP)
c      COS2F = COSFI(IP)*COSFI(IP) - SINFI(IP)*SINFI(IP)
c      ep1f=CMPLX(COSFI(IP),SINFI(IP))
       ep1f=cmplx(COS(fi),SIN(fi))
       em1f=conjg(ep1f)
c      ep2f=CMPLX(COS2F,SIN2F)
       ep2f=cmplx(COS(di),SIN(di))
       em2f=conjg(ep2f)

       do 21 l=0,lmax                                                  
c
          call spher(yl,l,cosx(ip))
        do 20 m=-l,l                                                    
          mm=l+m+1                                                      
          i=(l+1)**2-l+m   
          aaa=m*fai(ip)
          ccc=COS(aaa)
          sss=SIN(aaa)
          cylm0(i) = yl(mm) * CMPLX(ccc,sss)
c de moivre formula
c        cylm0(i) = yl(mm) * (CMPLX(COSFI(IP),SINFI(IP)))**M
c9010  format(1x,' i cylm0',i5,2f10.5)
   20   continue 

         do 22 m=-l,l                                                   
          i=(l+1)**2-l+m                                                
          cylmt1(i)=0.
          cylmt2(i)=0.
          cylmtf(i)=0.
   22    continue 

        do 23 m=-l,l                                                    
          i=(l+1)**2-l+m                                                

          lmm1m=l-m-1
          lmm=l-m
          lmm1=l-m+1
          lmm2=l-m+2
          lm1m=l+m-1
          lm=l+m
          lm1=l+m+1
          lm2=l+m+2

          cylmt2(i)=cylmt2(i)-(lmm*lm1+lmm1*lm)/4.*cylm0(i)

          if(m+2.le.l) cylmt2(i)=cylmt2(i)+
     &      SQRT(FLOAT(lmm1m*lmm*lm1*lm2))/4*cylm0(i+2)*em2f

          if(m+1.le.l) cylmt1(i)=cylmt1(i)+
     &        SQRT(FLOAT(lmm*lm1))/2*cylm0(i+1)*em1f

          if(m-1.ge.-l) cylmt1(i)=cylmt1(i)-
     &        SQRT(FLOAT(lm*lmm1))/2*cylm0(i-1)*ep1f

          if(m-2.ge.-l) cylmt2(i)=cylmt2(i)+
     &      SQRT(FLOAT(lmm1*lmm2*lm1m*lm))/4*cylm0(i-2)*ep2f

   23   continue 

         do 24 m=-l,l                                                   
          i=(l+1)**2-l+m                                                
          cylmf1(i)=ci*m*cylm0(i)
          cylmf2(i)=-m*m*cylm0(i)
          cylmtf(i)=ci*m*cylmt1(i)
   24    continue 

   21    continue 
c  
c        calculate real spherical harmonics differenciated
c
c
c        write(6,9005) (cylm0(i),i=1,5)
C9005 format(1x,' cylm0',4f10.5)
         call trarea(cylm0,bb1,lmax) 

         DO 31 M=1,llmax
   31    ylm(ip,m)=bb1(m)
c
c        write(6,9006) (ylm(ip,i),i=1,5)
c9006 format(1x,' ylm',10f10.5)
c
c
         call trarea(cylmt1,bb1,lmax)
         DO 32 M=1,llmax
   32    dylmt1(ip,m)=bb1(m)
C
         call trarea(cylmt2,bb1,lmax)
         DO 33 M=1,llmax
   33    dylmt2(ip,m)=bb1(m)
c
         call trarea(cylmf1,bb1,lmax)
         DO 34 M=1,llmax
   34    dylmf1(ip,m)=bb1(m)
C
         call trarea(cylmf2,bb1,lmax)
         DO 35 m=1,llmax
   35    dylmf2(ip,m)=bb1(m)
C
         call trarea(cylmtf,bb1,lmax)
         DO 36 m=1,llmax
   36    dylmtf(ip,m)=bb1(m)
C
C
   !      if(ip.eq.196.or.ip.eq.np) then
   !      write(6,9007) thet(ip),fai(ip),
   !  &         (ylm(ip,m),dylmt1(ip,m),dylmt2(ip,m),
   !  &         dylmf1(ip,m),dylmf2(ip,m),dylmtf(ip,m),m=1,5)
 ! 9007    format(1x,' thet fai ylm  ',7f9.5)
c        stop99
c        
   !     end if

   11 continue                                                          
         return                                                         
         end   
cdeck
c subroutine trarea
c
c     from complex to real  (differenciated spherical harmonics)
c
      SUBROUTINE TRAREA(a,b,lmax)
      implicit DOUBLE PRECISION (a-h,o-z)
      include 'inc.p' 
      INTEGER LMMAXD,lmax
      parameter(LMMAXD=(LPOTD+1)**2)
      INTEGER NN
      parameter(NN=2*(LPOTD+1))
      DOUBLE COMPLEX a(LMMAXD),ci
      DOUBLE PRECISION b(LMMAXD)
      data rtwo/1.414213562373d0/
      data pi,ci/3.14159265359d0,(0.d0,1.d0)/


C
C    calculate real the spherical harmonics derivetived
C
         i=0
         do 60 l=0,lmax
         i=i+l +1
!         b(i)=REAL(a(i))    ! Linux specific
         b(i)=DBLE(a(i)) 
c        write(6,9000) a(i),b(i)
 9000    format(1x,' a=',4f10.5,' b=',f10.5)
         sgm=-1.e0
         do 70 m=1,l
!         b(i-m)= REAL( ci*( a(i-m)-conjg(a(i-m)) ) )/rtwo
!         b(i+m)= sgm*REAL((a(i+m)+conjg(a(i+m))))/rtwo
c Linux specific
         b(i-m)= DBLE( ci*( a(i-m)-conjg(a(i-m)) ) )/rtwo
         b(i+m)= sgm*DBLE((a(i+m)+conjg(a(i+m))))/rtwo
c        write(6,9000) a(i-m),conjg(a(i-m)),b(i-m)
         sgm=-sgm
   70 continue
      i=i+l
   60 continue
      return
      end

cdeck
c subroutine spher  ====*====3====*====4====*====5====*====6====*====7  
c                                                                       
c      spherical haminics exept the facter exp(i*m*phi)                 
c                                                                       
c      m=-l to l , for given l.                                         
c      x=cos(theta)                                                     
c                                                                       
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7  
c                                                                       
      SUBROUTINE SPHER(ylm,l,x)                                         
c                                                                       
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      include 'inc.p' 
      INTEGER LMMAXD 
      PARAMETER (LMD0=2*LPOTD+1)                          
c                                                                       
      DOUBLE PRECISION ylm(2*l+1)                                              
c     REAL*8 ylm(LMD0)                                              
c                                                                       
      pi = 4.0*DATAN(1.0d0) 
      IF (L.GT.LPOTD) THEN  
      WRITE(6,*) 'ERROR IN SPHER : ',L,LPOTD 
      STOP 
      END IF                                               
c                                                                       
c                                                                       
      ovr1=ABS(x)-1.d0                                                  
      if(ovr1.gt.0.1d-12) then                                          
        write(6,200) x                                                  
  200   format(//3x,'==invalid argument for spher; x=',e24.16,' ==')    
        stop                                                            
      else if(abs(ovr1).lt.1.d-10) then                                 
        if(x.gt.0.0) then                                               
          fac=1.0                                                       
        else                                                            
          fac=(-1)**l                                                   
        end if                                                          
        l2=2*l+1                                                        
        do 10 i = 1,l2                                                  
   10   ylm(i) = 0.0                                                    
        ylm(l+1) = SQRT(DFLOAT(l2)/(4.0*pi))*fac                        
        return                                                          
      end if                                                            
c                                                                       
c l<0                                                                   
      if(l.lt.0) then                                                   
        write(6,*) ' === l=',l,' < 0  : in sub.spher. ==='              
        stop '=== stop in sub.spher. (l<0) ==='                         
c l=0                                                                   
      else if(l.eq.0) then                                              
        ylm(1) = SQRT(1.0/(4.0*pi))                                     
c l=1                                                                   
      else if(l.eq.1) then                                              
        fac = SQRT(3.0/(4.0*pi))                                        
        ylm(1) = fac*SQRT((1.0-x*x)/2.0)                                
        ylm(2) = fac*x                                                  
        ylm(3) = -ylm(1)                                                
c l>1                                                                   
      else                                                              
        ylm(1) = 1.0                                                    
        ylm(2) = x                                                      
        do 20 i = 2,l                                                   
   20   ylm(i+1) = ((2*i-1)*x*ylm(i)-(i-1)*ylm(i-1))/i                  
        fac = 1.0d0/SQRT(1.0d0- x*x)                                    
        do 30 m = 1,l                                                   
          lm = l + m                                                    
          ylm(lm+1)=fac*(-(l-m+1)*x*ylm(lm)+(lm-1)*ylm(l))              
          if (m.lt.l) then                                              
            nn = m +1                                                   
            do 32 i = nn,l                                              
            ii = l-i + nn                                               
   32       ylm(ii)=fac*(-(ii-m)*x*ylm(ii)+(ii+m-2)*ylm(ii-1))          
          end if                                                        
   30   continue                                                        
        fac = SQRT((2*l+1)/(4.0d0*pi))                                  
        ylm(l+1) = fac * ylm(l+1)                                       
        do 40 m = 1,l                                                   
          fac = -fac/SQRT(DFLOAT((l+m)*(l-m+1)))                        
          lm = l + 1 + m                                                
          ln = l + 1 - m                                                
          qq = ylm(lm)                                                  
          ylm(lm) = fac * qq                                            
          ylm(ln) = ABS(fac) * qq                                       
   40   continue                                                        
      end if                                                            
c                                                                       
      return                                                            
      end    












