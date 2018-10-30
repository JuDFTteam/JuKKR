  subroutine orb_rotation(u0,u1,orbrot)
! Constructs the orb rotation from u0 to u1
! See Rodrigues rotation formula e.g. on wikipedia
! CHECK: Y_L(R^-1 u) = \sum_L' R_LL' Y_L'(u) i.e. should the argument be inversely rotated?
! REVISE: mirroring or inversion?
  use global, only: iodb, i4b, r8b, c8b, nll, ull, wll, ylm, nlmax, lmmax, i2lm

  implicit none

! Initial direction
  real(kind=r8b), intent(in)  :: u0(3)
! Final direction
  real(kind=r8b), intent(in)  :: u1(3)
! Orb rotation matrix
  real(kind=r8b), intent(out) :: orbrot(lmmax,lmmax)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: tol = 1.d-6, rtol = 1.d-10
  real(kind=r8b)    :: axis(3), u0len, u1len, cosa, sina, urot(3)
  integer(kind=i4b) :: ilm, jlm, ill
  real(kind=c8b)    :: rotylm(lmmax) !, test1(lmmax,lmmax)
  logical           :: mirror

! Lengths of the vectors (if not unit vectors)
  u0len = sqrt(dot_product(u0,u0))
  u1len = sqrt(dot_product(u1,u1))
  if (u0len < tol .or. u1len < tol) stop 'orb_rotation: check u0 and u1'
  cosa = dot_product(u0,u1)/(u0len*u1len)
  sina = sqrt(1.d0 + cosa)*sqrt(1.d0 - cosa)
! Rotation axis: axis = u0 x u1 / | u0 x u1 |
  axis(1) = u0(2)*u1(3) - u0(3)*u1(2)
  axis(2) = u0(3)*u1(1) - u0(1)*u1(3)
  axis(3) = u0(1)*u1(2) - u0(2)*u1(1)
  u0len = sqrt(dot_product(axis,axis))
  write(iodb,'("axis=",3es16.8)') axis
! Collinear?
  if (u0len < tol) then
    if (cosa < 0.d0) then
!     mirroring
      mirror = .true.
      axis = u1/u1len
    else
!     identity
      orbrot = 0.d0
      do ilm=1,lmmax
        orbrot(ilm,ilm) = 1.d0
      end do
      return
    end if
! Rotation
  else
    mirror = .false.
    axis = axis/u0len
  end if
  write(iodb,'("mirror, rot axis=",l4,3es16.8)') mirror, axis
!  open(file='ymy.dat',status='replace',unit=1234)
!  do ill=0,1000
!    alpha = 4.d0*atan(1.d0)*ill/1001
!    urot(1) = cos(1.d0/3.d0)*sin(alpha)
!    urot(2) = sin(1.d0/3.d0)*sin(alpha)
!    urot(3) = cos(alpha)
!    call rymy(urot,nlmax,lmmax,rotylm)
!    write(1234,'(1000es16.8)') alpha, rotylm
!  end do
!  close(1234)
! Rotation matrix
  do jlm=1,lmmax
    do ilm=1,lmmax
      orbrot(ilm,jlm) = 0.d0
!     ------------------------------------------------------------------
!     rotation matrix is block diagonal
      if (i2lm(2,ilm) == i2lm(2,jlm)) then
!     ------------------------------------------------------------------
!      write(iodb,'("rot ilm, jlm=",2i4)') ilm, jlm
      do ill=1,nll
!       projection of the original vector on the axis of rotation
        u0len = dot_product(ull(:,ill),axis)
        if (mirror) then
!       mirrored unit vector
          urot(1) = ull(1,ill) - 2.d0*u0len*axis(1)
          urot(2) = ull(2,ill) - 2.d0*u0len*axis(2)
          urot(3) = ull(3,ill) - 2.d0*u0len*axis(3)
        else
!       rotated spherical unit vector: coordinates transformed by inverse rotation -> change sign of rotation angle
          urot(1) = cosa*ull(1,ill) - sina*(axis(2)*ull(3,ill) - axis(3)*ull(2,ill)) + (1.d0 - cosa)*u0len*axis(1)
          urot(2) = cosa*ull(2,ill) - sina*(axis(3)*ull(1,ill) - axis(1)*ull(3,ill)) + (1.d0 - cosa)*u0len*axis(2)
          urot(3) = cosa*ull(3,ill) - sina*(axis(1)*ull(2,ill) - axis(2)*ull(1,ill)) + (1.d0 - cosa)*u0len*axis(3)
        end if
        if (ilm == 1 .and. jlm == 1 .and. ill < 0) write(iodb,'("ull, urot=",6es16.8)') ull(:,ill), urot
!       Ylm for rotated angles
        call rymy(urot,nlmax,lmmax,rotylm)
!       matrix element
        orbrot(ilm,jlm) = orbrot(ilm,jlm) + wll(ill)*rotylm(ilm)*ylm(ill,jlm)
      end do
!     ------------------------------------------------------------------
      end if
!     ------------------------------------------------------------------
    end do
  end do
  where (abs(orbrot) < rtol) orbrot = 0.d0
!  write(iodb,'("Orbital rotation matrix times its transpose")')
!  test1 = matmul(orbrot,transpose(orbrot))
!  do ilm=1,lmmax
!    write(iodb,'(100f8.4)') test1(ilm,1:lmmax)
!  end do
! All done!
  end subroutine orb_rotation
