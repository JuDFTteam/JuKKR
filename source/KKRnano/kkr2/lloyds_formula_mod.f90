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
    integer, intent(inout), dimension(NGD) :: IPVT
    double complex, intent(out) :: LLY_G0TR

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

  !----------------------------------------------------------------------------
  ! This routine should be useable unchanged for Gref, as well as for its
  ! Derivative (used for Lloyd's formula)
  ! NEEDS TO KNOW ALL REFERENCE CLUSTERS AND ALL REFERENCE GREEN'S FUNCS !!!!
  ! HOW TO PARALLELISE ???

!  subroutine fourierTransformReferenceGreen(ALAT, NACLS, RR, EZOA, k_vector )
!
!    ! GLLH: temporary array - changed on output
!    ! DGDE: output -> RESULT
!    ! DGINP: the input
!
!    use kkr_helpers_mod
!    implicit none
!
!    double precision, dimension(3), intent(in) :: k_vector
!
!    ! local vars
!    integer :: site_index
!    integer :: ref_cluster_index
!    integer :: cluster_site_index
!    integer :: cluster_site_lm_index
!    integer :: LM1
!    integer :: LM2
!
!    ! local
!
!    GLLH = (0.0d0, 0.0d0)
!
!    do site_index = 1,NAEZ
!
!      ref_cluster_index = CLS(site_index)
!
!      call DLKE1(ALAT,NACLS,RR,EZOA(1,site_index), &
!                 k_vector, ref_cluster_index, EIKRM, EIKRP, &
!                 nrd, naclsd)
!
!      call DLKE0(site_index,GLLH, EIKRP, EIKRM, &
!                 ref_cluster_index,NACLS,ATOM(1,site_index),NUMN0,INDN0,DGINP(1,1,1,ref_cluster_index), &
!                 naez, lmax, naclsd)
!    end do
!
!    do site_index=1,NAEZ
!
!      do cluster_site_index=1,NUMN0(site_index)
!        do LM2=1,LMMAXD
!          cluster_site_lm_index=LMMAXD*(cluster_site_index-1)+LM2
!
!          if (INDN0(site_index,cluster_site_index) == IAT) then
!            do LM1=1,LMMAXD
!              site_lm_index=LMMAXD*(site_index-1)+LM1
!              DGDE(site_lm_index,LM2)= GLLH(LM1,cluster_site_lm_index,site_index)
!            enddo
!          endif
!
!        enddo
!      enddo
!
!    enddo
!
!end subroutine


  !> Calculates the energy derivative of the inverse of the
  !> scattering path operator.
  !>
  !> @param[in]     site_lm_size   number of sites * lmmaxd (ALM)
  !> @param[in]     lmmaxd
  !> @param[in]     alat           lattice parameter
  !> @param[in,out] DPDE_LOCAL
  !> @param[in,out] GLLKE_X        reference Greens-Function at (k,E)
  !> @param[in,out] DGDE           derivative of reference Greens-Function at (k,E)
  !> @param[in,out] DTmatDE_LOCAL  derivative of Delta T-Matrix
  !> @param[in,out] Tmat_local                   Delta T-Matrix
  !!                unchanged on output???
  subroutine calcDerivativeP(site_lm_size, lmmaxd, alat, &
                             DPDE_LOCAL, GLLKE_X, DGDE, DTmatDE_LOCAL, Tmat_local)
    ! calculate the following expression:

    ! dP(E,k)   dGref(E,k)                             d \Delta T(E)
    ! ------- = ---------- * \Delta T(E) + Gref(E,k) * -------------
    !   dE        dE                                    dE

     integer, intent(in) :: site_lm_size
     integer, intent(in) :: lmmaxd
     double precision, intent(in) :: alat

     double complex, dimension(site_lm_size, lmmaxd), intent(inout) :: DPDE_LOCAL
     double complex, dimension(site_lm_size, lmmaxd), intent(inout) :: GLLKE_X
     double complex, dimension(site_lm_size, lmmaxd), intent(inout) :: DGDE
     double complex, dimension(lmmaxd,lmmaxd), intent(inout) :: DTmatDE_LOCAL
     double complex, dimension(lmmaxd,lmmaxd), intent(inout) :: Tmat_local

     double complex, parameter :: CONE = ( 1.0D0,0.0D0)
     double complex, parameter :: CZERO= ( 0.0D0,0.0D0)

     double precision :: TWO_PI
     double complex   :: CFCTORINV

     TWO_PI = 8.D0*ATAN(1.D0)
     CFCTORINV = (CONE*TWO_PI)/ALAT

     DPDE_LOCAL = CZERO

     call ZGEMM('N','N',site_lm_size,LMMAXD,LMMAXD,CONE, &
                DGDE,site_lm_size, &
                Tmat_local,LMMAXD,CZERO, &
                DPDE_LOCAL,site_lm_size)

     ! WHY CFCTORINV ??? - must be a remainder from Fourier-transform
     call ZGEMM('N','N',site_lm_size,LMMAXD,LMMAXD,CFCTORINV, &
                GLLKE_X,site_lm_size, &
                DTmatDE_LOCAL,LMMAXD,CONE,DPDE_LOCAL,site_lm_size)

   end subroutine calcDerivativeP


   !---------------------------------------------------------------------------
   !> Calculates TRACE(X) for the real system. Only the local contribution for
   !> one k-point is calculated.
   !
   !                /  -1    dM  \
   ! calculate  Tr  | M   * ---- |
   !                \        dE  /

   ! NOTE: Later, don't forget to integrate over k! BZTR2 = BZTR2 + TRACEK*VOLCUB(k_point_index)

   !> @param[in]     site_lm_size   number of sites * lmmaxd (ALM)
   !> @param[in]     lmmaxd
   !> @param[in]     DPDE_LOCAL
   !> @param[in]     GLLKE1             scattering path operator
   !> @param[in]     inv_Tmat           inverse of local Delta T-Matrix (MSSQ) WHY???
   !> @param[out]    TRACEK             resulting trace

   subroutine calcLloydTraceXRealSystem(DPDE_LOCAL, GLLKE1, inv_Tmat, TRACEK, site_lm_size, lmmaxd)

     integer, intent(in) :: site_lm_size
     integer, intent(in) :: lmmaxd
     double complex, dimension(site_lm_size, lmmaxd), intent(in) :: DPDE_LOCAL
     double complex, dimension(site_lm_size, lmmaxd), intent(in) :: GLLKE1
     double complex, dimension(lmmaxd, lmmaxd), intent(in) :: inv_Tmat
     double complex, intent(out) :: TRACEK


     double complex, parameter :: CZERO= ( 0.0D0,0.0D0)
     double complex :: GTDPDE
     integer :: LM1
     integer :: LM2
     integer :: site_lm_index

     TRACEK=CZERO

     do LM1=1,LMMAXD
       do LM2=1,LMMAXD
         GTDPDE = CZERO
         do site_lm_index = 1, site_lm_size
           GTDPDE = GTDPDE + GLLKE1(site_lm_index,LM2)*DPDE_LOCAL(site_lm_index,LM1)
         enddo
         TRACEK = TRACEK + inv_Tmat(LM1,LM2)*GTDPDE   ! ?????? why
       enddo
     enddo
   ! NOTE: Later, don't forget to integrate over k! BZTR2 = BZTR2 + TRACEK*VOLCUB(k_point_index)
   ! TODO:
   ! after integration in k-space do the following:
   ! *) multiply trace with norm. factor (from k-space integration) and add TR_ALPH
   ! *) the contributions to the trace of all the processes have to be summed up (communication!)
   ! This is the original code that does this:
   !if(LLY == 1)  then
   !  BZTR2 = BZTR2*NSYMAT/VOLBZ + TR_ALPH
   !  TRACE=CZERO
   !  call MPI_ALLREDUCE(BZTR2,TRACE,1, &
   !  MPI_DOUBLE_COMPLEX,MPI_SUM, &
   !  LCOMM(LMPIC),IERR)
   !  LLY_GRDT = TRACE
   !endif

   end subroutine

!------------------------------------------------------------------------------
! Renormalise calculated density of states.
! previously in main2
! @param[in] LMAXD1 lmax + 1
! @param[in] IELAST number of energy points
! @param[in] IEMXD stupid dimension parameter for energy points IELAST <= IEMXD
! @param[in] RNORM normalisation factors: WARNING: 2nd dimension always 2 :-(
! WARNING: included spin loop!!! was originally in spin loop - move out!!!

  subroutine renormalizeDOS(DEN,RNORM,LMAXD1,IELAST,NSPIN,IEMXD)
    implicit none

    integer, intent(in) :: LMAXD1
    integer, intent(in) :: IELAST
    integer, intent(in) :: NSPIN
    integer, intent(in) :: IEMXD

    double complex, dimension(0:LMAXD1,IEMXD,NSPIN), intent(inout) :: DEN
    double precision, dimension(IEMXD,2), intent(in) :: RNORM

    integer :: IE
    integer :: ISPIN
    integer :: L_ind

    !DEBUG
    if (IELAST > IEMXD) then
      write(*,*) "renormalizeDOS: IELAST > IEMXD"
      write(*,*) IELAST, IEMXD
      stop
    end if

    if (NSPIN < 1 .or. NSPIN > 2) then
      write(*,*) "renormalizeDOS: invalid NSPIN value", NSPIN
      stop
    end if

    do ISPIN = 1,NSPIN
      do IE=1,IELAST
        do L_ind=0,LMAXD1
          DEN(L_ind,IE,ISPIN)=DEN(L_ind,IE,ISPIN)*RNORM(IE,ISPIN)
        end do
      end do
    end do

  end subroutine

end module lloyds_formula_mod
