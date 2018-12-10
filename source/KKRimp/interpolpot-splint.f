!------------------------------------------------------------------------------------
!> Summary: Interpolation of the potential using cubic spline from an old mesh to a new mesh 
!> Author: People who wrote it
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!>
!> One can write Latex comments like this \(i\hbar\frac{\partial \psi}{\partial t}=-\mathcal{H}\psi\)
!> or add labeled equations using the standard latex way
!> \begin{equation}
!> \mathbf{A} \mathbf{v}= \eta\mathbf{v}
!> \end{equation}
!> **FORd** also accepts markdown style so you can _write with style_ 
!> 
!> **IMPORTANT**
!> The JM-KKR follows the coding conventions noted in this example, one that is
!> not obvious is that **each level of indentation consists of two spaces**. Please keep this:
!> _do it for the children_.
!> So please keep the conventions.
! These are special boxes for ford, notice that this comment does not appear in the html file.
! These boxes contain important information and should be added when necessary. ALWAYS remember to close the box
! BEFORE oppening a new one or they will be nested.
!------------------------------------------------------------------------------------

      module mod_interpolspline

      contains

!-------------------------------------------------------------------------------
!> Summary: Cubic-spline interpolation of the potential using Spline and Splint subroutines
!> Author: Who wrote this subroutine
!> Category: Interpolation, Numerical-tools, potential, old-mesh, new-mesh
!> Deprecated: TRUE ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------


      subroutine interpolspline(rmesh,rmeshnew,vpot,vpotnew,
     +                          nrmax,nrmaxnew)
      use type_cell
      use type_cellnew
      implicit none
!interface
      double precision              :: rmesh(nrmax)
      double precision              :: rmeshnew(nrmaxnew)
      double precision              :: vpot(nrmax)
      double precision              :: vpotnew(nrmaxnew)
      integer                       :: nrmax
      integer                       :: nrmaxnew
!local
      double precision              :: maxa
      double precision              :: splinetemp(nrmax)
      double precision              :: PARSUM, PARSUMDERIV,r0
      integer                       :: ir
! Vpot,cell,lmaxatom,cellnew,ispin,nspin)
! implicit none
! double precision               :: Vpot(cell%nrmaxd,(lmaxatom+1)**2)
! type(cell_type)                :: cell
! integer                        :: lmaxatom
! type(cell_typenew)             :: cellnew
! integer                        :: ispin
! integer                        :: nspin


!     call interpolpot(cell%rmesh(irmin:irmax),cellnew%rmeshnew(irminnew:irmaxnew), &
!                      vpot(irmin:irmax,lm1),vpotnew(irminnew:irmaxnew,lm1),&
!                     irmax-irmin+1,irmaxnew-irminnew+1,8)

!            call interpolatecell(VPOT(:,:,ISPIN,IATOM),CELL(IATOM),lmaxatom(iatom),cellnew(iatom),ispin,nspin)

      maxa = 1.d35
      CALL SPLINE(nrmax,Rmesh,vpot,nrmax,maxa,maxa,splinetemp)  
!           CALL SPLINE(IRMDJJ,R,VM2Z,NR,maxa,maxa,VM2ZB)  


      DO IR = 1,nrmaxnew
         R0 = Rmeshnew(IR)
         CALL SPLINT(Rmesh,vpot,splinetemp,nrmax,R0,PARSUM,PARSUMDERIV)
         vpotnew(IR) = PARSUM
      END DO


!       end subroutine 

!             interpolatecell(VPOT(:,:,ISPIN,IATOM),CELL(IATOM),lmaxatom(iatom),cellnew(iatom),ispin,nspin)






      end subroutine  interpolspline






c***********************************************************************

!-------------------------------------------------------------------------------
!> Summary:  Spline subroutine is take from "Numerical Recipes in Fortran 77" and is used for interpolation
!> Given arrays x(1:n) and  y(1:n) containing a tabulated function, 
!> i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn 
!> for the 1rst derivative of the interpolating function at points 
!> 1 and n, respectively, this routine returns an array y2(1:n) of 
!> length n which contains the second derivatives of the interpolating 
!> function at the tabulated points xi. 
!> If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
!> signaled to set the corresponding boundary condition for a natural
!> spline, with zero second derivative on that boundary. 
!> Parameter: NMAX is the largest anticipated value of n. 
!> Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
!>
!> Author: Who wrote this subroutine
!> Category: Numerical-tools, Interpolation
!> Deprecated: TRUE ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------



      SUBROUTINE spline(NMAX,x,y,n,yp1,ypn,y2) 
      implicit none
      INTEGER n,NMAX 
      REAL*8          yp1,ypn,x(NMAX),y(NMAX),y2(NMAX) 
c Given arrays x(1:n) and  y(1:n) containing a tabulated function, 
c i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn 
c for the 1rst derivative of the interpolating function at points 
c 1 and n, respectively, this routine returns an array y2(1:n) of 
c length n which contains the second derivatives of the interpolating 
c function at the tabulated points xi. 
c If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
c signaled to set the corresponding boundary condition for a natural
c spline, with zero second derivative on that boundary. 
c Parameter: NMAX is the largest anticipated value of n. 
c Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      INTEGER i,k 
      REAL*8          p,qn,sig,un,u(NMAX) 

      if (n.gt.nmax) stop 'SPLINE: n > NMAX.'
      if (yp1.gt.0.99d30) then
c The lower boundary condition is set either to be "natural" 
         y2(1) = 0.d0
         u(1) = 0.d0
      else
c or else to have a specified first derivative. 
         y2(1) = -0.5d0
         u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      endif 

      do i = 2,n-1  
c This is the decomposition loop of the tridiagonal algorithm. y2 and u 
c are used for temporary storage of the decomposed factors. 
         sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
         p = sig * y2(i-1) + 2.d0 
         y2(i) = (sig-1.d0)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) 
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1)) / p 
      enddo

      if (ypn.gt.0.99d30) then   
c The upper boundary condition is set either to be "natural"
         qn = 0.d0
         un = 0.d0
      else
c or else to have a specified 1rst derivative. 
         qn = 0.5d0
         un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1)) / (qn*y2(n-1)+1.d0) 
      do k = n-1,1,-1 
c This is the backsubstitution loop of the tridiagonal algorithm. 
         y2(k)=y2(k)*y2(k+1)+u(k) 
      enddo

      return 
      END subroutine



c***********************************************************************

!-------------------------------------------------------------------------------
!> Summary:  Splint is taken from the "Numerical Recipes in Fortran 77" and is used for interpolation
!> Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
!> function (with the xai's in order), and given the array y2a(1:n), which
!> is the output from spline above, and given a value of x, this routine
!> returns a cubic-spline interpolated value y and the derivative yderiv.
!> Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
!> Author: Who wrote this subroutine
!> Category: Numerical-tools, Interpolation
!> Deprecated: TRUE ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------


      SUBROUTINE splint(xa,ya,y2a,n,x,y,yderiv)
      implicit none
      INTEGER n
      REAL*8         x,y,yderiv,xa(*),ya(*),y2a(*)
c Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
c function (with the xai's in order), and given the array y2a(1:n), which
c is the output from spline above, and given a value of x, this routine
c returns a cubic-spline interpolated value y and the derivative yderiv.
c Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
      INTEGER k,khi,klo
      REAL*8         a,b,h
c We will  nd the right place in the table by means of bisection.
c This is optimal if sequential calls to this routine are at random
c values of x. If sequential calls are in order, and closely
c spaced, one would do better to store previous values of
c klo and khi and test if they remain appropriate on the
c next call.
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
         else
            klo=k
         endif
         goto 1
      endif
c klo and khi now bracket the input value of x.
      h=xa(khi)-xa(klo)
c The xa's must be distinct.
      if (h.eq.0.d0) pause 'bad xa input in splint'
c Cubic spline polynomial is now evaluated.
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo) + b*ya(khi) +
     &     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi)) * (h**2)/6.d0
      yderiv = (ya(khi)-ya(klo))/h - 
     &     ((3.d0*a*a-1.d0)*y2a(klo) - (3.d0*b*b-1.d0)*y2a(khi))*h/6.d0

      return
      END SUBROUTINE









      end module mod_interpolspline
