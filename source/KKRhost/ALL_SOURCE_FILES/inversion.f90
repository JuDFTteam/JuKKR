module mod_inversion

contains

! ************************************************************************
subroutine inversion(gllke, invmod, icheck)
  ! ************************************************************************
  ! This subroutine calculates the inversion of a matrix
  ! in 4 different ways depending on the form of the matrix

  ! INVMOD = 0  ----> total inversion scheme
  ! INVMOD = 1  ----> band matrix inversion scheme
  ! INVMOD = 2  ----> corner band matrix inversion scheme
  ! INVMOD = 3  ----> godfrin module

  ! ------------------------------------------------------------------------
  use global_variables
  use mod_datatypes, only: dp
  use godfrin
   use mod_invslab
   use mod_invsupercell
  implicit none

  complex (kind=dp) :: ci, czero, cone
  parameter (ci=(0.e0_dp,1.e0_dp), czero=(0.e0_dp,0.e0_dp), &
    cone=(1.e0_dp,0.e0_dp))

  complex (kind=dp) :: gllke(alm, alm), gdi(ndim_slabinv, ndim_slabinv, nlayerd), &
    gup(ndim_slabinv, ndim_slabinv, nlayerd), gdow(ndim_slabinv, ndim_slabinv, nlayerd)
  complex (kind=dp), allocatable :: gtemp(:, :)
  integer :: i, i1, ip1, ii1, il1, ldi1, ip2, ii2, il2, ldi2, j, invmod
  integer :: lm1, lm2, info, ipvt(alm), nlayer
  integer :: icheck(naezd/nprincd, naezd/nprincd)
  ! total matrix inversion

  allocate (gtemp(alm,alm))

  ! Naive number of layers in each principal layer
  nlayer = naezd/nprincd

  ! ---------------------------------------------------------------------
  !                       full matrix inversion
  ! ---------------------------------------------------------------------
  if (invmod==0) then

    ! initialize unit matrix
    do i = 1, alm
      do j = 1, alm
        gtemp(i, j) = czero
        if (i==j) then
          gtemp(i, j) = cone
        end if
      end do
    end do

    call zgetrf(alm, alm, gllke, alm, ipvt, info)
    call zgetrs('N', alm, alm, gllke, alm, ipvt, gtemp, alm, info)
    call zcopy(alm*alm, gtemp, 1, gllke, 1)

  ! ---------------------------------------------------------------------
  !           block tridiagonal inversion (slab or supercell)
  ! ---------------------------------------------------------------------
  else if ((invmod==1) .or. (invmod==2)) then
    ! copy blocks of assumed tridiagonal matrix
    do i1 = 1, nlayerd
      do ip1 = 1, nprincd
        do ip2 = 1, nprincd
          ii1 = (i1-1)*nprincd + ip1
          ii2 = (i1-1)*nprincd + ip2
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              ldi1 = lmmaxd*(ip1-1) + lm1
              il1 = lmmaxd*(ii1-1) + lm1
              ldi2 = lmmaxd*(ip2-1) + lm2
              il2 = lmmaxd*(ii2-1) + lm2
              gdi(ldi1, ldi2, i1) = gllke(il1, il2)
            end do
          end do
        end do
      end do
    end do
    ! this part now is correct also for    ! changes 20/10/99
    ! supercell geometry : 20/10/99
    !--> upper linear part
    DO I1=1,NLAYERD
      DO IP1=1,NPRINCD
      DO IP2=1,NPRINCD
        DO LM1=1,LMMAXD
        DO LM2=1,LMMAXD
          LDI1 = LMMAXD*(IP1-1)+LM1
          LDI2 = LMMAXD*(IP2-1)+LM2
          IF (I1.LE.(NLAYERD-1)) THEN
            II1 = (I1-1)*NPRINCD+IP1
            II2 = I1*NPRINCD+IP2
            IL1 = LMMAXD*(II1-1)+LM1
            IL2 = LMMAXD*(II2-1)+LM2
            GUP(LDI1,LDI2,I1) =  GLLKE(IL1,IL2)
          ELSE
            II1 = IP1
            II2 = (NLAYERD-1)*NPRINCD+IP2
            IL1 = LMMAXD*(II1-1)+LM1
            IL2 = LMMAXD*(II2-1)+LM2
            GDOW(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
          END IF
        END DO
        END DO
      END DO
      END DO
    END DO
    !--> lower linear part
    DO I1=1,NLAYERD 
      DO IP1=1,NPRINCD
      DO IP2=1,NPRINCD    
        DO LM1=1,LMMAXD
        DO LM2=1,LMMAXD
          LDI1 = LMMAXD*(IP1-1)+LM1
          LDI2 = LMMAXD*(IP2-1)+LM2
          IF (I1.LE.(NLAYERD-1)) THEN
            II1 = I1*NPRINCD+IP1
            II2 = (I1-1)*NPRINCD+IP2
            IL1 = LMMAXD*(II1-1)+LM1
            IL2 = LMMAXD*(II2-1)+LM2
            GDOW(LDI1,LDI2,I1) =  GLLKE(IL1,IL2)
          ELSE
            II1 = (NLAYERD-1)*NPRINCD+IP1
            II2 = IP2
            IL1 = LMMAXD*(II1-1)+LM1
            IL2 = LMMAXD*(II2-1)+LM2
            GUP(LDI1,LDI2,I1) = GLLKE(IL1,IL2)
          END IF
        END DO
        END DO
      END DO
      END DO
    END DO
    ! end of the corrected part  20/10/99

!   slab: matrix is block tridiagonal
    if (invmod==1) then
      ! write (6,*) '-------slab calculation--------'
      call invslab(gdi, gup, gdow, gllke, icheck)
    else if (invmod==2) then
      ! write (6,*) '-------supercell calculation--------'
      call invsupercell(gdi, gup, gdow, gllke, icheck)
    end if

  ! -----------------------------------------------------------------
  !                       godfrin module
  ! -----------------------------------------------------------------
  else if (invmod==3) then
    call sparse_inverse(gllke,t_godfrin%na,t_godfrin%nb, &
      t_godfrin%bdims,t_godfrin%ldiag,t_godfrin%lper,    &
      t_godfrin%lpardiso)  ! GODFRIN Flaviano
  !------------------------------------------------------------------
  else
    ! if it gets here, did you have a coffee before running the code?
    stop 'UNKNOWN INVERSION MODE !!!'
  endif
  !------------------------------------------------------------------

  deallocate (gtemp)

  return

end subroutine inversion

end module mod_inversion
