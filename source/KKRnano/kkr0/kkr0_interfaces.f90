module kkr0_interfaces
  interface

subroutine TESTDIM( &
  NSPIN, NAEZ, LMAX, IRM, NREF, IRNS, NCLS &
  )
  integer, intent(in) :: NSPIN
  integer, intent(in) :: NAEZ
  integer, intent(in) :: LMAX
  integer, intent(in) :: IRM
  integer, intent(in) :: NREF
  integer, intent(in) :: NCLS
  integer, dimension(*), intent(in) :: IRNS

end subroutine TESTDIM

! reads the initial potential
! reads shapefunctions
! creates radial mesh
!subroutine STARTB1( &
!IFILE, IPF, IPFE, IPE, KVREL, KHFELD, LMAX, NBEG, NEND, RMTNEW, RMT,  &
!ITITLE, HFIELD, IMT, IRC, VCONST, IRNS, LPOT, NSPIN, IRMIN, NTCELL,  &
!IRCUT, IPAN, THETAS, IFUNM, NFU, LLMSP, LMSP, EFERMI, LRECPOT, VBC,  &
!RWS, LCORE, NCORE, DRDI, R, ZAT, A, B, IRWS, INIPOL, IINFO &
!)
!  integer, intent(in) :: IFILE
!  integer, intent(in)  :: IPF
!  integer, intent(in)  :: IPFE
!  integer, intent(in)  :: IPE
!  integer, intent(in)  :: KVREL
!  integer, intent(in)  :: KHFELD
!  integer, intent(in)  :: LMAX
!  integer, intent(in) :: NBEG
!  integer, intent(in) :: NEND
!  real(kind=8), intent(in)  :: HFIELD
!  real(kind=8), intent(in)  :: VCONST
!  integer, intent(in)  :: LPOT
!  integer, intent(in)  :: NSPIN
!  real(kind=8), intent(inout)  :: EFERMI
!  integer, intent(in)  :: LRECPOT
!  integer, intent(in) :: IINFO
!
!  real(kind=8), dimension(*), intent(inout) :: RMTNEW
!  real(kind=8), dimension(*), intent(inout) :: RMT
!  integer, dimension(20,*), intent(inout) :: ITITLE
!  integer, dimension(*), intent(inout) :: IMT
!  integer, dimension(*), intent(inout) :: IRC
!
!  integer, dimension(*), intent(in) :: IRNS
!  integer, dimension(*), intent(inout) :: IRMIN
!
!  integer, dimension(*), intent(inout) :: NTCELL
!  !integer, dimension(?,?), intent(inout) :: IRCUT
!  integer, dimension(*), intent(inout) :: IPAN
!
!  !real(kind=8), dimension(?,?,?), intent(inout) :: THETAS
!  !integer, dimension(?,?), intent(inout) :: IFUNM
!
!  integer, dimension(*), intent(inout) :: NFU
!  !integer, dimension(?,?), intent(inout) :: LLMSP
!  !integer, dimension(?,?), intent(inout) :: LMSP
!  real(kind=8), dimension(*), intent(inout) :: VBC
!  real(kind=8), dimension(*), intent(inout) :: RWS
!  integer, dimension(20,*), intent(inout) :: LCORE
!  integer, dimension(*), intent(inout) :: NCORE
!  !real(kind=8), dimension(?,?), intent(inout) :: DRDI
!  !real(kind=8), dimension(?,?), intent(inout) :: R
!  real(kind=8), dimension(*), intent(inout) :: ZAT
!  real(kind=8), dimension(*), intent(inout) :: A
!  real(kind=8), dimension(*), intent(inout) :: B
!  integer, dimension(*), intent(inout) :: IRWS
!
!  integer, dimension(*), intent(in) :: INIPOL
!
!end subroutine STARTB1

function TEST(STRING)

  logical :: TEST
  character(len=8), intent(in) :: STRING

end function TEST


! This routine only writes some output and checks if IMT is odd
subroutine CALRMT( &
IPF, IPFE, IPE, IMT, Z, RMT, RWS, RMTNEW, ALAT, DRDI, A, B, IRWS, R,  &
IFILE &
)
  integer, intent(in) :: IPF
  integer, intent(in) :: IPFE
  integer, intent(in) :: IPE
  integer, intent(in) :: IMT
  real(kind=8), intent(in) :: Z
  real(kind=8), intent(in) :: RMT
  real(kind=8), intent(in) :: RWS
  real(kind=8), intent(in) :: RMTNEW
  real(kind=8), intent(in) :: ALAT
  real(kind=8), intent(in) :: A
  real(kind=8), intent(in) :: B
  integer, intent(in) :: IRWS
  integer, intent(in) :: IFILE
  real(kind=8), dimension(*), intent(in) :: DRDI
  real(kind=8), dimension(*), intent(in) :: R

end subroutine CALRMT

subroutine RCSTOP(C)
  character(len=8), intent(in) :: C
end subroutine RCSTOP

  end interface
end module kkr0_interfaces
