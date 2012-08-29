! JUST FOR TESTING purposes
! replace by proper implementation
module TEST_lcutoff_mod
  implicit none

  SAVE

  integer lm_low
  double precision cutoff_radius
  integer, dimension(:), allocatable :: lmarray
  integer :: cutoffmode

  CONTAINS

  subroutine initLcutoff(rbasis, bravais, lmmaxd, atomindex)
    use lcutoff_mod
    implicit none
    double precision, dimension(:,:), intent(in) :: rbasis
    double precision, dimension(3,3), intent(in) :: bravais
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: atomindex

    allocate(lmarray(size(rbasis,2)))

    open(91, file='lcutoff', form='formatted')
      read(91,*) cutoff_radius
      read(91,*) lm_low
      read(91,*) cutoffmode
    close(91)

    call getLMarray(lmarray, rbasis, rbasis(:,atomindex), bravais, cutoff_radius, lmmaxd, lm_low)

  end subroutine

  subroutine cropGLLH(GLLH, lmmaxd, naclsd, naezd, lmarray, numn0, indn0)
    implicit none
    double complex GLLH(LMMAXD,NACLSD*LMMAXD,NAEZD)
    integer lmmaxd
    integer naclsd
    integer naezd
    integer, dimension(naezd) :: lmarray
    integer, dimension(naezd) :: numn0
    integer, dimension(naezd, naclsd) :: indn0
    !-----
    integer ii, jj
    integer lmmax1, lmmax2, lm1, lm2
    integer clustersitelm
    double complex, parameter :: CZERO = (0.0d0, 0.0d0)

    do ii = 1, naezd
      do jj = 1, numn0(ii)

      lmmax1 = lmarray(ii)
      lmmax2 = lmarray(indn0(ii, jj))

        do lm2 = 1, lmmaxd
          do lm1 = 1, lmmaxd

          clustersitelm = (jj-1)*lmmaxd + lm2

          if (lm1 > lmmax1 .or. lm2 > lmmax2) then
            GLLH(lm1, clustersitelm, ii) = CZERO
          end if

          end do
        end do


      end do
    end do
  end subroutine


!------------------------------------------------------------------------------
!> Generates matrix (\Delta T G_ref - 1) BUT PERFORMS L-CUTOFF BY TRUNCATING
!> T-MATRIX
!> on input: GLLH contains G_ref, on output: GLLH contains coefficient matrix
subroutine generateCoeffMatrixCROPPED(GLLH, NUMN0, INDN0, TMATLL, NAEZ, lmmaxd, naclsd, lmarray)
  implicit none

  double complex, parameter :: CONE  = ( 1.0D0,0.0D0)
  double complex, parameter :: CZERO = ( 0.0D0,0.0D0)

  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd
  integer :: NAEZ
  double complex :: GLLH(LMMAXD,NACLSD*LMMAXD,NAEZ)
  integer :: INDN0(NAEZ,NACLSD)
  integer :: NUMN0(NAEZ)
  doublecomplex :: TMATLL(lmmaxd,lmmaxd,NAEZ)
  integer, dimension(:) :: lmarray

  !---------- local --------------
  double complex :: TGH(lmmaxd)
  integer :: IL1B
  integer :: IL2B
  integer :: LM1
  integer :: LM2
  integer :: LM3
  integer :: site_index
  integer :: site_lm_index
  integer :: cluster_site_index
  integer :: cluster_site_lm_index

  integer :: lmmax1, lmmax2, lmmax3

  ! -------------- Calculation of (Delta_t * G_ref - 1) ---------------
  !
  !
  ! NUMN0(site_index) is the number of atoms in the reference cluster
  ! of atom/site 'site_index' (inequivalent atoms only!)
  ! INDN0 stores the index of the atom in the basis corresponding to
  ! the reference cluster atom (inequivalent atoms only!)
  ! -------------------------------------------------------------------

  !$omp parallel do private(site_index, site_lm_index, cluster_site_index, &
  !$                        cluster_site_lm_index, IL1B, IL2B, &
  !$                        LM1, LM2, LM3, TGH)
  do site_index=1,NAEZ
    IL1B=LMMAXD*(site_index-1)
    do cluster_site_index=1,NUMN0(site_index)

    lmmax2 = lmarray(INDN0(site_index,cluster_site_index))
    lmmax1 = lmarray(site_index)
    lmmax3 = lmmax1

      do LM2=1,LMMAXD
        cluster_site_lm_index=LMMAXD*(cluster_site_index-1)+LM2
        IL2B=LMMAXD*(INDN0(site_index,cluster_site_index)-1)+LM2

        TGH = CZERO
        do LM1=1,lmmax1
          do LM3=1,lmmax3
            TGH(LM1)=TGH(LM1)+TMATLL(LM1,LM3,site_index)*GLLH(LM3,cluster_site_lm_index,site_index)
          enddo !lm3
        enddo !lm1

        do LM1=1,LMMAXD
          site_lm_index=IL1B+LM1
          GLLH(LM1,cluster_site_lm_index,site_index) = TGH(LM1)

          if (site_lm_index == IL2B) then
            ! substract 1 only at the 'diagonal'
            GLLH(LM1,cluster_site_lm_index,site_index) = GLLH(LM1,cluster_site_lm_index,site_index) - CONE
          endif

          if (lm1 > lmmax1 .or. lm2 > lmmax2) then
            GLLH(LM1,cluster_site_lm_index,site_index) = CZERO
          end if

        enddo !lm1
      enddo !lm2

    enddo
  enddo
  !$omp end parallel do
end subroutine

end module TEST_lcutoff_mod
