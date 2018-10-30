  subroutine spinrot_rhomat(magdir0,magdir1)
! Rotates the site-diagonal density matrix according to local spin rotations
  use global

  implicit none

! Initial spin axes
  real(kind=r8b), intent(in) :: magdir0(3,nasusc)
! Final spin axis
  real(kind=r8b), intent(in) :: magdir1(3,nasusc)
! ----------------------------------------------------------------------
  logical, parameter :: lhund = .true.
  complex(kind=c8b) :: spinrot(nsmax,nsmax)
  real(kind=r8b)    :: orbrot(lmmax,lmmax)
  complex(kind=c8b), allocatable :: rhomatsave(:,:)
  integer(kind=i4b) :: i2(2), ia
  integer(kind=i4b) :: i, j, is, js, ilm, jlm
  integer(kind=i4b) :: p, q, ps, qs, plm, qlm


  allocate(rhomatsave(nlms,nlms))
! **************
  do ia=1,nasusc
! **************
! Spin rotation matrices
    call spin_rotation(magdir0(:,ia),magdir1(:,ia),pauli,spinrot)
!    write(*,'(i4,8f8.4)') ia, spinrot(:,:,ia)
! Orbital rotation matrices
    call orb_rotation(magdir0(:,ia),magdir1(:,ia),orbrot)
!   Rotation of the density matrix
    rhomatsave = rhomat(:,:,ia)
    rhomat(:,:,ia) = 0.d0
    do j=1,nlms
      i2 = i2lms(:,j)
      jlm = i2(1); js = i2(2)
      do q=1,nlms
        i2 = i2lms(:,q)
        qlm = i2(1); qs = i2(2)
        do i=1,nlms
          i2 = i2lms(:,i)
          ilm = i2(1); is = i2(2)
          do p=1,nlms
            i2 = i2lms(:,p)
            plm = i2(1); ps = i2(2)
!           rotation of the Hund density matrix to a new axis
            if (lhund) then
              rhomat(i,j,ia) = rhomat(i,j,ia) + orbrot(ilm,plm)*spinrot(is,ps)*rhomatsave(p,q)*conjg(spinrot(js,qs))*orbrot(jlm,qlm)
!           spin rotation only
            else
!             selection rules
              if (ilm == plm .and. jlm == qlm) then
                rhomat(i,j,ia) = rhomat(i,j,ia) + spinrot(is,ps)*rhomatsave(p,q)*conjg(spinrot(js,qs))
              end if
            end if
          end do
        end do
      end do
    end do
! ******
  end do
! ******
  deallocate(rhomatsave)
! All done!
  end subroutine spinrot_rhomat
