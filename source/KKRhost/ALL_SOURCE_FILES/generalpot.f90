! **********************************************************************
SUBROUTINE generalpot(ifile,natps,natyp,nspin,z,alat,rmt,rmtnew,  &
    rws,r,drdi,vm2z,irws,a,b,ins,irns,  &
    lpot,vins,qbound,irc,kshape,efermi,vbc,ecore, lcore,ncore,lmpotd,irmd,irmind)
! **************************************************
! * The subroutine writes out the potential cards
! * in a standard r-mesh that can be read in and
! * interpolated to a different r-mesh from subroutine
! * start No shape function information is needed
! * and all nessecery data are stored in the potential
! * card.
! *                                      ver. 18.5.2000
! ***************************************************
!     ..
IMPLICIT NONE
!..
!.. Scalar Arguments ..
INTEGER LMPOTD,IRMD,IRMIND
DOUBLE PRECISION ALAT,QBOUND
INTEGER IFILE,INS,KSHAPE,LPOT,NATPS,NATYP,NSPIN
!..
!.. Array Arguments ..
DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),ECORE(20,*),EFERMI, &
                 R(IRMD,*),RMT(*),RMTNEW(*),RWS(*),VBC(2), &
                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
INTEGER IRC(*),IRNS(*),IRWS(*),LCORE(20,*),NCORE(*)
!..
!.. Local Scalars ..
DOUBLE PRECISION A1,B1,RMAX,RMT1,RMTNW1,RV,SUM,Z1,PARSUM, &
                 PARSUMDERIV,R0,RINTER,DR,MAXA
INTEGER I,ICORE,IH,IP,IR,IRMIN,IRNS1,IS,ISAVE,J,LM,LMNR, &
        LMPOT,NCORE1,NR,NZ1,NR_U,IRMIN_U,IRNS_U, &
        IMT1,LM1,IRNSTOT
!..
!.. Local Arrays ..
DOUBLE PRECISION DRADI(IRMD),ECORE1(20),RA(IRMD),VM2ZA(IRMD), &
                 RR_U(IRMD),DRDI_U(IRMD)
DOUBLE PRECISION VM2ZB(IRMD),VM2Z_U(IRMD), &
            VINS_U(IRMIND:IRMD,LMPOTD),VINSA(IRMIND:IRMD,LMPOTD), &
            VINSB(IRMIND:IRMD,LMPOTD)
INTEGER LCORE1(20)
 CHARACTER*4 ELEMNAME(0:113)
!..
!.. Intrinsic Functions ..
INTRINSIC SQRT
!        1      2      3      4      5      6      7      8      9    
DATA ELEMNAME/'VAC', &
  'H   ','He  ','Li  ','Be  ','B   ','C   ','N   ','O   ','F   ', &
  'Ne  ', &
  'Na  ','Mg  ','Al  ','Si  ','P   ','S   ','Cl  ','Ar  ','K   ', &
  'Ca  ', &
  'Sc  ','Ti  ','V   ','Cr  ','Mn  ','Fe  ','Co  ','Ni  ','Cu  ', &
  'Zn  ', &
  'Ga  ','Ge  ','As  ','Se  ','Br  ','Kr  ','Rb  ','Sr  ','Y   ', &
  'Zr  ', &
  'Nb  ','Mo  ','Tc  ','Ru  ','Rh  ','Pd  ','Ag  ','Cd  ','In  ', &
  'Sn  ', &
  'Sb  ','Te  ','I   ','Xe  ','Cs  ','Ba  ','La  ','Ce  ','Pr  ', &
  'Nd  ', &
  'Pm  ','Sm  ','Eu  ','Gd  ','Tb  ','Dy  ','Ho  ','Er  ','Tm  ', &
  'Yb  ', &
  'Lu  ','Hf  ','Ta  ','W   ','Re  ','Os  ','Ir  ','Pt  ','Au  ', &
  'Hg  ', &
  'Tl  ','Pb  ','Bi  ','Po  ','At  ','Rn  ','Fr  ','Ra  ','Ac  ', &
  'Th  ', &
  'Pa  ','U   ','Np  ','Pu  ','Am  ','Cm  ','Bk  ','Cf  ','Es  ', &
  'Fm  ', &
  'Md  ','No  ','Lr  ','Rf  ','Db  ','Sg  ','Bh  ','Hs  ','Mt  ', &
  'Uun ','Uuu ','Uub ','NoE '/

isave = 1
lmpot = (lpot+1)* (lpot+1)

DO  ih = 1,natyp
  DO  is = 1,nspin
    DO lm =1,lmpotd
      DO ir =irmind,irmd
        vinsa(ir,lm) = 0.d0
        vinsb(ir,lm) = 0.d0
      END DO
    END DO
    ip = nspin* (ih-1) + is
    rmt1 = rmt(ih)
    rmtnw1 = rmtnew(ih)
    z1 = z(ih)
    rmax = rws(ih)
    IF (kshape == 0) THEN
      nr = irws(ih)
      irns1 = 0
    ELSE
      nr = irc(ih)
      irns1 = irns(ih)
    END IF
    
    irmin = nr - irns1
    a1 = a(ih)
    b1 = b(ih)
    ncore1 = ncore(ip)
    
    DO  j = 1,nr
      ra(j) = r(j,ih)
      dradi(j) = drdi(j,ih)
      vm2za(j) = vm2z(j,ip)
    END DO
    DO lm1=1,lmpot
      DO j=irmind,irmd
        vinsa(j,lm1) = vins(j,lm1,ip)
      END DO
    END DO
    
    IF (ncore1 >= 1) THEN
      
      DO  j = 1,ncore1
        lcore1(j) = lcore(j,ip)
        ecore1(j) = ecore(j,ip)
      END DO
    END IF
    
!  Generate uniform mesh RUNI
    
    nr_u = nr
    irns_u = irns1
    irmin_u = nr_u
    IF (ins > 0) irmin_u = nr_u - irns_u
    
    IF (ins == 0) THEN
      DO i=1,nr_u
        rr_u(i) = ra(i)
        drdi_u(i) = dradi(i)
      END DO
      imt1 = 0
    ELSE
      imt1 = ANINT(LOG(rmtnw1/b1+1.0D0)/a1) + 1
      DO i=1,imt1
        rr_u(i) = ra(i)
        drdi_u(i) = dradi(i)
      END DO
      rinter =  rmax - rmtnw1
      dr = rinter/FLOAT(nr-imt1)
      DO i=1,nr-imt1
        drdi_u(imt1+i) = dr
        rr_u(imt1+i) = rr_u(imt1) + dr*FLOAT(i)
      END DO
      CALL doubleraus1(nr,irmin,lmpot,ra,dradi,vm2za,vinsa,  &
          irmd,irmind,lmpotd)
      
!     After this sub the arrays are rearanged and nr is not
!     the same anymore in the case of FP. If ins.eq.0 there is
!     no nead for doubleraus1. IRMIN should remain the same
      
    END IF
! ----------------------------------------------------------------
! Now the new mesh is generated
    
! test
!  write(6,*) nr_u,imt1,irns_u
!  do i=1,nr_u
!    write(6,*) i,ra(i),rr_u(i)
!  end do
! test
    
    maxa = 1.d35
    CALL spline(irmd,ra,vm2za,nr,maxa,maxa,vm2zb)
    IF (ins > 0) THEN
      DO lm1=1,lmpot
        irnstot = nr - irmin ! nr has changed irmin is the same
!write(6,*) ' Testing ',nr,irmin,irnstot,irmind
        CALL spline(irmd-irmind,ra(irmind),  &
            vinsa(irmind,lm1),irnstot,maxa,maxa, vinsb(irmind,lm1))
      END DO              ! LM1
    END IF
    
! OK with spline
    
    DO ir = 1,nr_u
      r0 = rr_u(ir)
      CALL splint(ra,vm2za,vm2zb,nr,r0,parsum,parsumderiv)
      vm2z_u(ir) = parsum
    END DO
    IF (ins > 0) THEN
!IRNSTOT = NR_U - IRMIN_U
      DO lm1=1,lmpot
        DO ir = irmin_u,nr_u
          r0 = rr_u(ir)
          CALL splint(ra(irmind),vinsa(irmind,lm1),  &
              vinsb(irmind,lm1),irnstot,r0, parsum,parsumderiv)
          vins_u(ir,lm1) = parsum
        END DO
      END DO
    END IF
!write(6,*) ' All interpolation ok now write'
!     --------------------------------------------------------------
    WRITE (ifile,FMT=8000)
    nz1 = z1
    IF (nspin == 1) THEN
      WRITE (ifile,FMT=8010) elemname(nz1),z1
    ELSE IF (is == 1) THEN
      WRITE (ifile,FMT=8012) elemname(nz1),z1
    ELSE IF (is == 2) THEN
      WRITE (ifile,FMT=8011) elemname(nz1),z1
    END IF
    WRITE (ifile,FMT=8020)
!          write (ifile,*) ALAT,RMAX,RMTNW1,RMT1
    WRITE (ifile,FMT=8030) alat,rmax,rmtnw1,rmt1
    WRITE (ifile,FMT=8040) nr_u,imt1,irns1
    WRITE (ifile,FMT=8050) a1,b1
    WRITE (ifile,FMT=8060) efermi,vbc(is)
    WRITE (ifile,FMT=8070) ncore1,lmpot
    IF (ncore1 >= 1) WRITE (ifile,FMT=9040) (lcore1(icore),  &
        ecore1(icore),icore=1,ncore1)
    
    IF (ins == 0 .OR. (ih < natps.AND.ins <= 2)) THEN
      
!---  >       store only the spherically averaged potential
!     (in mt or as - case)
!     this is done always for the host
      
      WRITE (ifile,FMT=9051) (vm2z_u(ir),ir=1,nr_u)
    ELSE
      
!---  >     store the full potential , but the non spherical contribution
!     only from irns1 up to irws1 ;
!     remember that the lm = 1 contribution is multiplied
!     by a factor 1/sqrt(4 pi)
      
      WRITE (ifile,FMT=9060) nr_u,irns1,lmpot,isave
      WRITE (ifile,FMT=9070) (vm2z_u(ir),ir=1,nr_u)
      IF (lpot > 0) THEN
        lmnr = 1
        DO  lm = 2,lmpot
          sum = 0.0D0
          DO  ir = irmin,nr_u
            rv = vins_u(ir,lm)*rr_u(ir)
            sum = sum + rv*rv*dradi(ir)
          END DO
          
          IF (SQRT(sum) > qbound) THEN
            lmnr = lmnr + 1
            WRITE (ifile,FMT=9060) lm
            WRITE (ifile,FMT=9070) (vins_u(ir,lm),ir=irmin,nr_u)
          END IF
          
        END DO
        
!---  >         write a one to mark the end
        
        IF (lmnr < lmpot) WRITE (ifile,FMT=9060) isave
      END IF
      
    END IF
    
  END DO
END DO


8000 FORMAT (' GENERAL POTENTIAL MESH             exc:')
8010 FORMAT ('#  ',a4,'POTENTIAL             Z = ',f8.3)
8011 FORMAT ('#  ',a4,'POTENTIAL SPIN UP     Z=  ',f8.3)
8012 FORMAT ('#  ',a4,'POTENTIAL SPIN DOWN   Z=  ',f8.3)
8020 FORMAT ('#')
8030 FORMAT (4F12.8, '   # alat, rmax, rmaxlog, rmt')
8040 FORMAT (1P,3I6,31X,'  # IRWS, IRMT, IRNS ')
8050 FORMAT (2D15.8,19X,'  # A , B ')
8060 FORMAT (3F12.8,13X,'  # Ef, vbc ')
8070 FORMAT (1P,2I5,39X,'  # NCORE, LMPOT' )
9000 FORMAT (7A4,6X,'  exc:',a24,3X,a10)
9010 FORMAT (3F12.8)
9020 FORMAT (f10.5,/,f10.5,2F15.10)
9030 FORMAT (i3,/,2D15.8,/,2I2)
9040 FORMAT (i5,1P,d20.11)
9050 FORMAT (1P,2D15.6,1P,d15.8)
9051 FORMAT (1P,4D20.12)
9060 FORMAT (10I5)
9070 FORMAT (1P,4D20.13)
END SUBROUTINE generalpot
! **********************************************************************

!***********************************************************************

SUBROUTINE spline(nmax,x,y,n,yp1,ypn,y2)
! Given arrays x(1:n) and  y(1:n) containing a tabulated function,
! i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn
! for the rst derivative of the interpolating function at points
! 1 and n, respectively, this routine returns an array y2(1:n) of
! length n which contains the second derivatives of the interpolating
! function at the tabulated points xi.
! If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
! signaled to set the corresponding boundary condition for a natural
! spline, with zero second derivative on that boundary.
! Parameter: NMAX is the largest anticipated value of n.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      implicit none
      integer n,nmax 
      double precision yp1,ypn,x(nmax),y(nmax),y2(nmax) 
      INTEGER I,K 
      DOUBLE PRECISION P,QN,SIG,UN,U(NMAX) 

IF (n > nmax) STOP 'SPLINE: N > NMAX.'
IF (yp1 > 0.99D30) THEN
! The lower boundary condition is set either to be "natural"
  y2(1) = 0.d0
  u(1) = 0.d0
ELSE
! or else to have a specified first derivative.
  y2(1) = -0.5D0
  u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
END IF

DO i = 2,n-1
! This is the decomposition loop of the tridiagonal algorithm. y2 and u
! are used for temporary storage of the decomposed factors.
  sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
  p = sig * y2(i-1) + 2.d0
  y2(i) = (sig-1.d0)/p
  u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))  &
      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1)) / p
END DO

IF (ypn > .99D30) THEN
! The upper boundary condition is set either to be "natural"
  qn = 0.d0
  un = 0.d0
ELSE
! or else to have a specified 1rst derivative.
  qn = 0.5D0
  un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
END IF
y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.d0)
DO k = n-1,1,-1
! This is the backsubstitution loop of the tridiagonal algorithm.
  y2(k)=y2(k)*y2(k+1)+u(k)
END DO

RETURN
END SUBROUTINE spline
! **********************************************************************

! **********************************************************************
SUBROUTINE splint(xa,ya,y2a,n,x,y,yderiv)
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the xai's in order), and given the array y2a(1:n), which
! is the output from spline above, and given a value of x, this routine
! returns a cubic-spline interpolated value y and the derivative yderiv.
! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X,Y,YDERIV,XA(*),YA(*),Y2A(*)
      INTEGER K,KHI,KLO
      DOUBLE PRECISION A,B,H
! We will find the right place in the table by means of bisection.
! This is optimal if sequential calls to this routine are at random
! values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of
! klo and khi and test if they remain appropriate on the
! next call.
klo=1
khi=n
1    IF (khi-klo > 1) THEN
  k=(khi+klo)/2
  IF(xa(k) > x)THEN
    khi=k
  ELSE
    klo=k
  END IF
  GO TO 1
END IF
! klo and khi now bracket the input value of x.
h=xa(khi)-xa(klo)
! The xa's must be distinct.
IF (h == 0.d0) STOP 'Bad XA input in SPLINT'
! Cubic spline polynomial is now evaluated.
a = (xa(khi)-x)/h
b = (x-xa(klo))/h
y = a*ya(klo) + b*ya(khi) +  &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * (h**2)/6.d0
yderiv = (ya(khi)-ya(klo))/h -  &
    ((3.d0*a*a-1.d0)*y2a(klo) - (3.d0*b*b-1.d0)*y2a(khi))*h/6.d0

RETURN
END SUBROUTINE splint
! **********************************************************************

! **********************************************************************
!************************************************************************

SUBROUTINE doubleraus1(irmax,irmin,lmpot,rr,drdi,vpot,vins,  &
    irmd,irmind,lmpotd)
! Gets rid of the double-points in the radial mesh, i.e. the points
! where RR(I) = RR(I-1). Returns the "new" mesh in the same array,
! and rearranges accordingly the WAVEF defined at the same mesh.
! IRMAX is also altered to the new value.
      IMPLICIT NONE
      INTEGER IRMD,LMPOTD,IRMIND
      INTEGER NCOUNTMAX
      PARAMETER(NCOUNTMAX=500)
! Input and output:
      INTEGER IRMAX
      DOUBLE PRECISION RR(IRMD),DRDI(IRMD),VPOT(IRMD), &
     &                 VINS(IRMIND:IRMD,LMPOTD)
! Inside:
      INTEGER IR,ICOUNT,NC,NCOUNT,LMPOT,IRMIN
      INTEGER LM1
      INTEGER IDOUBLE(NCOUNTMAX)

! Find double-points:
ncount = 0
DO ir = 2,irmax
  IF ((rr(ir)-rr(ir-1)) < 1.d-20) THEN
    ncount = ncount + 1
    idouble(ncount) = ir
  END IF
END DO
IF (ncount+1 > ncountmax) STOP 'DOUBLERAUS2: Too many double-points.'
idouble(ncount+1) = irmax+1   ! To be used below.

! Rearrange the arrays.
DO icount = 1,ncount
  DO ir = idouble(icount)-icount+1,idouble(icount+1)-icount
    IF((ir+icount) <= irmax) THEN
      rr(ir) = rr(ir+icount)
      drdi(ir) = drdi(ir+icount)
      vpot(ir) = vpot(ir+icount)
    END IF
  END DO
END DO
irmax = irmax - ncount
ncount = 0
DO nc = 1,ncountmax
  idouble(nc) = 0
END DO

DO ir = irmin,irmax
  IF ((rr(ir)-rr(ir-1)) < 1.d-20) THEN
    ncount = ncount + 1
    idouble(ncount) = ir
  END IF
END DO
! Rearrange the arrays.
DO icount = 1,ncount
  DO ir = idouble(icount)-icount+1,idouble(icount+1)-icount
    DO lm1=1,lmpot
      vins(ir,lm1) = vins(ir+icount,lm1)
    END DO
  END DO
END DO
RETURN
END SUBROUTINE doubleraus1
! **********************************************************************
