  subroutine build_dti_v2(ie,ia,lldautoo,dti)
! Commutators of t-matrix with generators of spin and orbital rotations

  use global

  implicit none

! Which energy
  integer(kind=i4b), intent(in) :: ie
! Which atom
  integer(kind=i4b), intent(in) :: ia
! Whether to add LDA+U correction to spin splitting
  logical,  intent(in) :: lldautoo
! t-matrix perturbations
  complex(kind=c8b), intent(inout) :: dti(nlms,nlms,3)
! ----------------------------------------------------------------------
  real(kind=r8b), parameter    :: tol = 1.d-4
  complex(kind=c8b), parameter :: iu = (0.d0,1.d0), cone = (1.d0,0.d0), czero = (0.d0,0.d0)
  integer(kind=i4b) :: i2(2), i, j, ilms, jlms, is, js, ilm, jlm
  real(kind=r8b)    :: axis(3), axislen, udir(3,3)
  complex(kind=c8b) :: genrot(nlms,nlms,3), fac


! Assumption: ASA
! Construct generators for rotations
  udir(:,1) = (/1.d0,0.d0,0.d0/)
  udir(:,2) = (/0.d0,1.d0,0.d0/)
  udir(:,3) = (/0.d0,0.d0,1.d0/)
  genrot = czero
  do i=1,3
    do jlms=1,nlms
      i2 = i2lms(:,jlms)
      jlm = i2(1); js = i2(2)
      do ilms=1,nlms
        i2 = i2lms(:,ilms)
        ilm = i2(1); is = i2(2)
        if (is  == js ) genrot(ilms,jlms,i) = genrot(ilms,jlms,i) + lorb(ilm,jlm,i)
        if (ilm == jlm) genrot(ilms,jlms,i) = genrot(ilms,jlms,i) + pauli(is,js,i)/2.d0
      end do
    end do
  end do
! ----------------------------------------------------------------------
  do j=1,3
    dti(:,:,j) = czero
!   magdir x udir
    axis(1) = magdir(2,ia)*udir(3,j) - magdir(3,ia)*udir(2,j)
    axis(2) = magdir(3,ia)*udir(1,j) - magdir(1,ia)*udir(3,j)
    axis(3) = magdir(1,ia)*udir(2,j) - magdir(2,ia)*udir(1,j)
    axislen = sqrt(dot_product(axis,axis))
    if (axislen > tol) then
      axis = axis/axislen
      if (ie == nescf) write(*,'("axis=",3f8.4)') axis
!     i [ tmat, (L + S).axis ]
      do i=1,3 
        fac = +iu*axis(i)
        call zgemm('N','N',nlms,nlms,nlms,fac,tmatrix(:,:,ia,ie),nlms,genrot(:,:,i),nlms,cone,dti(:,:,j),nlms)
        fac = -iu*axis(i)
        call zgemm('N','N',nlms,nlms,nlms,fac,genrot(:,:,i),nlms,tmatrix(:,:,ia,ie),nlms,cone,dti(:,:,j),nlms)
      end do
    end if
  end do
! All done!
  end subroutine build_dti_v2
