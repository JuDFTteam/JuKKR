! It is planned to move some blocks of Lloyd's formula
! calculations here.

! Dependencies: kkr_helpers_mod
module lloyds_formula_mod
  implicit none

  contains

  !----------------------------------------------------------------------------
  !> Calculates TRACE(X) where
  !> X = (1 - g0 * \Delta t_ref)^(-1) * d/dE (1 - g0 * \Delta t_ref).
  !> This is needed for Lloyd's formula
  !> The matrix X can also be retrieved but it is not needed for
  !> Lloyd's formula.
  !> previously in GREFSY
  !> @param[in,out] LUfact   the LU-factorisation of (1 - g0 t_ref) as obtained in GREFSY
  !!                         dim (NGD,NGD)
  !> @param[in,out] DGTDE    INPUT: the following energy derivative:
  !!                         d/dE (1 - g0 * \Delta t_ref)
  !!                         Only for central atom: first LMGF0D columns !
  !!                         CHANGED ON OUTPUT: content: matrix X
  !> @param[in,out] IPVT     integer work array of dimension NGD
  !> @param[out]    LLY_G0TR Contains result: TRACE(X)
  !> @param[in]     NGD      leading dimension of matrix LUfact, NGD=lmmaxd*naclsd
  !> @param[in]     NDIM     logical dimension of matrix NDIM=lmmaxd * (#atoms in cluster)
  !> @param[in]     LMGF0D   same as lmmaxd ???
  subroutine calcLloydTraceX(LUfact, DGTDE, IPVT, LLY_G0TR, &
                             NDIM, LMGF0D, NGD)
    implicit none

    double complex, intent(inout), dimension(NGD,NGD) :: LUfact
    double complex, intent(inout), dimension(NGD,LMGF0D) :: DGTDE
    integer, intent(in) :: NDIM  ! actual dimension
    integer, intent(in) :: LMGF0D
    integer, intent(in) :: NGD   ! lmmaxd*naclsd  leading dimension

    double complex, parameter :: CZERO = (0.0d0, 0.0d0)
    integer :: I
    integer :: info

    CALL ZGETRS('N',NDIM,LMGF0D,LUfact,NGD,IPVT,DGTDE,NGD,info)

    LLY_G0TR = CZERO

    DO I = 1,LMGF0D
       LLY_G0TR = LLY_G0TR - DGTDE(I,I)  ! why minus?
    ENDDO

  end subroutine

  !----------------------------------------------------------------------------
  !> Calculates the energy derivative of the reference Green's function.
  !> previously in gll95
  !> This is needed for Lloyd's formula
  !> @param[in,out] LUfact  the LU-factorisation of (1 - g0 t_ref) as obtained in GREFSY
  !!                        dim (NGD,NGD)
  !> @param[in,out] GREF0   the reference Green's function
  !!                        dim (NGD,LMGF0D)
  !> @param[in,out] DGDE    on input:  the derivative of the free space Green's function
  !!                        on output: derivative of the reference Green's function
  !!                        dim (NGD,LMGF0D) -> RESULT
  !> @param[in,out] DGTDE0  the following energy derivative: d/dE (1 - g0 * \Delta t_ref)
  !!                        unchanged on output
  !> @param[in,out] IPVT    integer work array of dimension NGD
  !> @param[in]     NGD     leading dimension of matrix LUfact, NGD=lmmaxd*naclsd
  !> @param[in]     NDIM    logical dimension of matrix NDIM=lmmaxd * (#atoms in cluster)
  !> @param[in]     LMGF0D  same as lmmaxd ???

  subroutine calcDerivativeGref(LUfact, GREF0, DGTDE0, DGDE, IPVT, &
                                NDIM, LMGF0D, NGD)
    implicit none

    integer, intent(in) :: NDIM  ! actual dimension
    integer, intent(in) :: LMGF0D
    integer, intent(in) :: NGD   ! lmmaxd*naclsd  leading dimension
    integer, intent(inout), dimension(NGD) :: IPVT
    double complex, intent(inout), dimension(NGD,NGD) :: DGTDE0
    double complex, intent(inout), dimension(NGD,LMGF0D) :: GREF0
    double complex, intent(inout), dimension(NGD,NGD) :: DGDE ! This is the output!
    double complex, intent(inout), dimension(NGD,NGD) :: LUfact

    double complex, parameter :: CONE= (1.D0,0.D0)
    integer :: INFO

    !debug
    if (NDIM > NGD) then
      write(*,*) "calcDerivativeGref: NDIM > NGD"
      stop
    end if

    call ZGEMM('N','N',NDIM,LMGF0D,NDIM,-CONE,DGTDE0,NGD, &
                GREF0,NGD,CONE,DGDE,NGD)

    ! use LU-factorisation stored in LUfact
    call ZGETRS('N',NDIM,LMGF0D,LUfact,NGD,IPVT,DGDE,NGD,INFO)

     !do N2 = 1,LMGF0D
     !  do N1 = 1,NGD
     !    DGDEOUT(N1,N2) = DGDE(N1,N2)
     !  enddo
     !enddo

  end subroutine

  !----------------------------------------------------------------------------
  ! previously in gll95
  ! TODO: lmmaxd and LMGF0D are the same????
  subroutine calcDerivativeFreeGreens(energy, RATOM, NATOM, ALAT, &
                                      DGDE, NGD, &
                                      lmmaxd, LMGF0D, &
                                      CLEB, ICLEB, LOFLM, IEND, ncleb) ! Gaunt coeffs
    use kkr_helpers_mod
    implicit none

    double precision, intent(in) :: RATOM(3,NATOM)
    double precision, intent(in) :: ALAT
    double complex, intent(in) :: energy
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: NATOM
    integer, intent(in) :: LMGF0D
    integer, intent(in) :: NGD

    double complex, dimension(NGD,NGD), intent(inout) :: DGDE

    double precision, intent(inout) :: CLEB(NCLEB)
    integer, intent(inout) :: ICLEB(NCLEB,3)
    integer, intent(inout) :: LOFLM(*)  ! WHAT IS THE DIMENSION ???!!!
    integer, intent(inout) :: IEND
    integer, intent(in) :: ncleb

    ! local variables
    double complex, dimension(LMGF0D,LMGF0D) :: DGLLDE

    double complex, parameter :: CZERO = (0.0d0, 0.0d0)
    double precision, dimension(3) :: RDIFF
    integer :: N1
    integer :: N2
    integer :: LM1
    integer :: LM2
    integer :: ind
    integer :: site_lm_index1
    integer :: site_lm_index2
    integer :: lmaxd
     !
     ! ---> construct derivative of free Green's function
     !
     if (lmmaxd /= lmgf0d) then
       write(*,*) "calcDerivativeFreeGreens: lmmaxd /= lmgf0d"
       stop
     end if

     !lmaxd = sqrt(lmmaxd - 1) ! TODO: replace by something better
     lmaxd = lmmaxTolmax(lmmaxd)

     do N1 = 1,NATOM
       do N2 = 1,NATOM
         do ind = 1,3
           RDIFF(ind) = - (RATOM(ind,N1)-RATOM(ind,N2))*ALAT
         end do

         if (N1/=N2) then

           call DGFREE(RDIFF,energy,DGLLDE,CLEB,ICLEB,LOFLM,IEND,lmaxd, ncleb)

           do LM2 = 1,LMGF0D
             site_lm_index2 = (N2-1)*LMGF0D + LM2
             do LM1 = 1,LMGF0D
               site_lm_index1 = (N1-1)*LMGF0D + LM1
               DGDE(site_lm_index1,site_lm_index2) = DGLLDE(LM1,LM2)
             end do
           end do
         else
           do LM2 = 1,LMGF0D
             site_lm_index2 = (N2-1)*LMGF0D + LM2
             do LM1 = 1,LMGF0D
               site_lm_index1 = (N1-1)*LMGF0D + LM1
               DGDE(site_lm_index1,site_lm_index2) = CZERO
             end do
           end do

         end if  ! (N1 /= N2)

       end do
     end do

  end subroutine

end module lloyds_formula_mod
