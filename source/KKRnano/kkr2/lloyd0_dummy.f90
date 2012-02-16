!> Dummy routine for lloyd0 - for debugging purposes.
subroutine LLOYD0_DUMMY(EZ,WEZ,CLEB,DRDI,R,IRMIN, &
                  VINS,VISP,THETAS,ZAT,ICLEB, &
                  IFUNM1,IPAN,IRCUT,LMSP1,JEND,LOFLM,NTCELL,ICST, &
                  IELAST,IEND,NAEZ,NSPIN,NSRA, &
                  WEZRN,RNORM, &                                     ! <
                  GMATN, &                                           ! >
                  LLY_GRDT, &                                        ! >
                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &                 ! >
                  DMATLDAU, &                                        ! <
                  LMPIC,MYLRANK, &                                   ! >
                  LCOMM,LSIZE, &                                     ! >
                  !   new input parameters after inc.p removal
                  prod_lmpid_smpid_empid, lmax, irmd, irnsd, iemxd, &
                  irid, nfund, ncelld, ipand, ncleb)

  implicit none
  include 'mpif.h'

  integer :: lmax
  integer :: irmd
  integer :: irnsd
  integer :: iemxd
  integer :: prod_lmpid_smpid_empid
  integer :: irid
  integer :: nfund
  integer :: ncelld
  integer :: ipand
  integer :: ncleb

  double complex :: EZ(IEMXD)  ! in
  double complex :: WEZ(IEMXD) ! in
  double precision::CLEB(NCLEB,2)   !in
  double precision::DRDI(IRMD,NAEZ) !in
  double precision::R(IRMD,NAEZ)    !in
  integer::IRMIN(NAEZ) !in
  double precision::VINS((IRMD-IRNSD):IRMD,(2*LMAX+1)**2, 2)  !in?
  double precision::VISP(IRMD,2) ! in?
  double precision::THETAS(IRID,NFUND,NCELLD) !in
  double precision::ZAT(NAEZ) !in
  integer::ICLEB(NCLEB,3)
  !integer::IFUNM1(LMXSPD,NAEZ) (2*LPOTD+1)**2
  integer::IFUNM1((4*LMAX+1)**2,NAEZ) !in?
  integer::IPAN(NAEZ) !in
  integer::IRCUT(0:IPAND,NAEZ) !in
  !integer::LMSP1(LMXSPD,NAEZD)
  integer::LMSP1((4*LMAX+1)**2, NAEZ) !in
  integer::JEND((LMAX+1)**2, 0:LMAX, 0:LMAX) !in
  integer::LOFLM((2*LMAX+1)**2) !in
  integer::NTCELL(NAEZ) !in
  integer::ICST
  integer::IELAST
  integer::IEND
  integer::NAEZ  !in
  integer::NSPIN !in
  integer::NSRA  !in
  double complex :: WEZRN(IEMXD,2)  ! out
  double precision::RNORM(IEMXD,2)  !out
  double complex :: GMATN((LMAX+1)**2, (LMAX+1)**2, IEMXD, NSPIN) ! inout?
  double complex :: LLY_GRDT(IEMXD,NSPIN) ! in

  ! -------------- LDA+U --------------------------------------------
  !double complex :: DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  double complex :: DMATLDAU(2*LMAX+1,2*LMAX+1,NSPIN,LMAX+1)  ! inout?
  double complex :: PHILDAU(IRMD,LMAX+1) !in?
  double precision::WMLDAU(2*LMAX+1,2*LMAX+1,NSPIN,LMAX+1) !in?
  integer::NLDAU !in
  logical::LDAU  !in
  integer::LLDAU(LMAX+1) !in?
  !------------------------------------------------------------------

  integer::MYLRANK(prod_lmpid_smpid_empid) !in
  integer::LCOMM(prod_lmpid_smpid_empid)   !in
  integer::LSIZE(prod_lmpid_smpid_empid)   !in
  integer::LMPIC                           !in

!---------------- Local variables -----------------------------------

  integer :: ii, jj

  do jj = 1, 2
    do ii = 1, iemxd
      WEZRN(ii, jj) = WEZ(ii)
      RNORM(ii, jj) = 1.0
    end do
  end do

end subroutine LLOYD0_DUMMY
