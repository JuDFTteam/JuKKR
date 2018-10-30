  subroutine mix_rhomat(mix)
! Mix density matrix (for LDA+U)
  use global

  implicit none

  real(kind=r8b), intent(in) :: mix
! ---------------------------------------------------------------------
  real(kind=r8b), parameter :: tol = 1.d-8
  complex(kind=c8b) :: newrhomat(nlms,nlms,nasusc)
  integer(kind=i4b) :: ia, il, ilms
  real(kind=r8b)    :: uminusj, qin, qout

! Save current rhomat
  newrhomat = rhomat
! Read rhomat from previous iteration
  call read_rhomat
! Compute vshift for frozen density matrix
  if (abs(ldaumix) < tol) then
    do ia=1,nasusc
      if (ildau(ia) == 1) then
        do il=0,nlmax
!         --------------------------------------------------------------
          uminusj = ueff(il,ia) - jeff(il,ia)
          if (abs(uminusj) > tol) then
            qin = 0.d0; qout = 0.d0
            do ilms=1,nlms
              if (i2lm(2,i2lms(1,ilms)) == il) then
                qin  = qin  + rhomat(ilms,ilms,ia)
                qout = qout + newrhomat(ilms,ilms,ia)
              end if
            end do
!           Penalty for wrong occupancy
!            vshift(il,1:nsmax,ia) = 0.d0
!            if (abs(uminusj*(qout-qin)) < 0.1d0) vshift(il,1:nsmax,ia) = uminusj*(qout-qin)
            write(iodb,'("mix_rhomat: ia, il, qin, qout, vshift=",2i4,3f12.6)') ia, il, qin, qout, vshift(il,1,ia)
          end if
!         --------------------------------------------------------------
        end do
      end if
    end do
  end if
! Update rhomat:
  do ia=1,nasusc
    if (abs(sum(rhomat(:,:,ia))) < tol) then
      rhomat(:,:,ia) = newrhomat(:,:,ia)
    else
      rhomat(:,:,ia) = (1.d0 - mix)*rhomat(:,:,ia) + mix*newrhomat(:,:,ia)
    end if
  end do
! All done!
  end subroutine mix_rhomat
