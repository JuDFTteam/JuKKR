  subroutine get_rhoden(ia,lmmax,gfsum,rhoden)
! components of the magnetization in the basis

  implicit none

! which atom, angular momentum cutoff
  integer(kind=i4b), intent(in)  :: ia, lmmax
! density matrix
  complex(kind=c8b), intent(in)  :: gfsum(nlmsb,nlmsb)
! density matrix the density basis
  complex(kind=c8b), intent(out) :: rhoden(ndenmax)
! ---------------------------------------------------
  real(kind=r8b), parameter :: fourpi = 16.d0*atan(1.d0)
  integer(kind=i4b) :: nr, ir, i2(2), i3(3), i4(4), k, n, ill, ngf0, ngf1, nden0, nden1
  integer(kind=i4b) :: i, ilm, il, im, ib, is, iq
  integer(kind=i4b) :: j, jlm, jl, jm, jb, js, jq
  complex(kind=c8b) :: rotation(ndenmax,ndenmax), rhodenp(ndenmax), multipoles(nsmax2,lmmax0)


  if (lmmax < lmmax0) stop 'get_rhoden: lmmax < lmmax0'

! ----------------------------------------------------------------------
! Convert density matrix to the density basis
  rhodenp = 0.d0; multipoles = 0.d0
  nden0 = sum(nalmsbden(1:ia-1))
  nden1 = sum(nalmsbden(1:ia))
  do jq=1,nden1-nden0
    i4 = i2almsbden(:,jq+nden0)
    ib = i4(1); ilm = i4(2); is = i4(3)
    ngf0 = sum(nalmsbgf(1:ia-1))
    ngf1 = sum(nalmsbgf(1:ia))
    do iq=1,ngf1-ngf0
      if (abs(dengaunt(iq,jq,ia)) > ylmtol) then
        i3 = i2almsbgf(:,iq+ngf0)
        i = i3(1); j = i3(2)!; ia = i3(3)
        rhodenp(jq) = rhodenp(jq) + dengaunt(iq,jq,ia)*gfsum(i,j)
      end if
    end do
!   Set transverse components to zero
    if (is < 3) rhodenp(jq) = 0.d0
  end do
! Rotation matrix
  rotation = 0.d0
  do jq=1,nden1-nden0
    i4 = i2almsbden(:,jq+nden0)
    jb = i4(1); jlm = i4(2); js = i4(3)
    do iq=1,nden1-nden0
      i4 = i2almsbden(:,iq+nden0)
      ib = i4(1); ilm = i4(2); is = i4(3)
      if (ilm == jlm .and. ib == jb) then
      if (is == 1 .and. js == 1) rotation(iq,jq) =  0.5d0
      if (is == 2 .and. js == 1) rotation(iq,jq) = -0.5d0
      if (is == 3 .and. js == 1) rotation(iq,jq) = -0.5d0
      if (is == 4 .and. js == 1) rotation(iq,jq) =  0.5d0
      if (is == 1 .and. js == 2) rotation(iq,jq) = -0.5d0
      if (is == 2 .and. js == 2) rotation(iq,jq) =  0.5d0
      if (is == 3 .and. js == 2) rotation(iq,jq) = -0.5d0
      if (is == 4 .and. js == 2) rotation(iq,jq) =  0.5d0
      if (is == 1 .and. js == 3) rotation(iq,jq) =  0.5d0
      if (is == 2 .and. js == 3) rotation(iq,jq) =  0.5d0
      if (is == 3 .and. js == 3) rotation(iq,jq) =  0.5d0
      if (is == 4 .and. js == 3) rotation(iq,jq) =  0.5d0
      if (is == 1 .and. js == 4) rotation(iq,jq) = -0.5d0
      if (is == 2 .and. js == 4) rotation(iq,jq) = -0.5d0
      if (is == 3 .and. js == 4) rotation(iq,jq) =  0.5d0
      if (is == 4 .and. js == 4) rotation(iq,jq) =  0.5d0
      end if
    end do
  end do
! Apply rotation to bring the magnetization from z to x
  rhoden = matmul(rotation,rhodenp)
! Test: multipoles
  do jq=1,nden1-nden0
    i4 = i2almsbden(:,jq+nden0)
    ib = i4(1); ilm = i4(2); is = i4(3)
    if (is > 2) rhoden(jq) = 0.d0
    multipoles(is,ilm) = multipoles(is,ilm) + rhoden(jq)*suscnorm(jq+nden0)
  end do
  do ilm=1,lmmax0
    write(*,'("get_rhoden: multipoles=",2i4,8f12.8)') i2lm(:,ilm), multipoles(:,ilm)
  end do
! All done!
  end subroutine get_rhoden

