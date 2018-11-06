!------------------------------------------------------------------------------------
!> Summary: Module handling common array dimensions
!> Author:
!> Array dimensions are defined according to lmax; used in main program
!------------------------------------------------------------------------------------
!> @note Notes on the code
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!------------------------------------------------------------------------------------
module arrayparams
use nrtype
! initial global variables
!  INTEGER                          ::  NRMAXD                      ! maximum number of rad. points for allocationg
!  INTEGER                          ::  NPAND  ,NRSHAPED,  NLMSHAPED
!  INTEGER                          ::  NCORESTATED                ! maximum number of core states (for allocating)
   INTEGER                          ::  LMAXD                      ! maximum lmax from atominfo
   INTEGER                          ::  LPOTD
!       PARAMETER ( NSPIND = KREL + (1-KREL)*(KSP+1) )
! ( LPOTD = 2*LMAXD )
!       PARAMETER ( NCLEB = (LMAXD*2+1)**2 * (LMAXD+1)**2 )


!derived global variables
   INTEGER        ::     LMMAXD
!      PARAMETER (LMMAXD= (KREL+1) * (LMAXD+1)**2)
   INTEGER        ::     LMPOTD
!      PARAMETER (LMPOTD= (LPOTD+1)**2)
!    INTEGER        ::     IRMIND
!       PARAMETER (IRMIND=IRMD-IRNSD)
   INTEGER        ::     MMAXD
!       PARAMETER ( MMAXD = 2*LMAXD+1 )
!    INTEGER        ::     LM2D
!       PARAMETER (LMMAXD= (2*LMAXD+1)**2)
!    REAL(KIND=DP)  ::     CVLIGHT
!       PARAMETER (CVLIGHT=274.0720442D0)
   INTEGER LM2D
!       PARAMETER (LM2D= (2*LMAXD+1)**2)

integer :: IRMIND
integer :: IRMD,IRMAXD
integer :: IEMXD

integer :: INS
integer :: IRMTD
integer :: NCORED
contains

!-------------------------------------------------------------------------------
!> Summary: Define common array dimensions based on lmax
!> Author:
!> Category: initialization, KKRimp
!> Deprecated: False 
!> Array dimensions are defined according to lmax
!-------------------------------------------------------------------------------
!> @note Notes on the code
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!-------------------------------------------------------------------------------
subroutine arrayparams_set(LMAXD1)
integer,save          :: first = 1

integer :: lmaxd1

IF (first/=1) stop '[array_params] Trying to change the array parameters more then one time not permitted'
first=0

LMAXD=LMAXD1
LPOTD=2*LMAXD
LMMAXD = (LMAXD+1)**2
LMPOTD = (LPOTD+1)**2
MMAXD  =  2*LMAXD+1
LM2D= (2*LMAXD+1)**2

end subroutine arrayparams_set
end module arrayparams
