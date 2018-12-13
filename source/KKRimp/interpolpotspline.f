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
      use mod_splint, only: splint_real
      use mod_spline, only: spline_real
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
      CALL spline_real(nrmax,Rmesh,vpot,nrmax,maxa,maxa,splinetemp)  
!           CALL SPLINE(IRMDJJ,R,VM2Z,NR,maxa,maxa,VM2ZB)  


      DO IR = 1,nrmaxnew
         R0 = Rmeshnew(IR)
         CALL splint_real(Rmesh,vpot,splinetemp,nrmax,R0,PARSUM,
     &       PARSUMDERIV)
         vpotnew(IR) = PARSUM
      END DO


!       end subroutine 

!             interpolatecell(VPOT(:,:,ISPIN,IATOM),CELL(IATOM),lmaxatom(iatom),cellnew(iatom),ispin,nspin)



      end subroutine  interpolspline




      end module mod_interpolspline
