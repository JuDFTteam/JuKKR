module mod_symlat

contains

subroutine symlat(nsymop, platcp, symopm)
  ! - Supplies the point symmetry operations of the lattice
  ! ----------------------------------------------------------------------
  ! i Inputs:
  ! i   platcp:lattice vectors of most compact primitive unit cell
  ! o Outputs:
  ! o   nsymop:number of allowed symmetry operations
  ! o   symopm:symmetry operation matrix
  ! r Remarks:
  ! r   symlat analyzes the primitive translations of the bravais
  ! r   lattice in order to supply the symmetry operations of the lattice.
  ! r   It gives the number nsymop of allowed operations as well as
  ! r   these operations themselves.
  ! ----------------------------------------------------------------------

  use :: mod_datatypes, only: dp
  use mod_dinv33
  use mod_dmpy
  use mod_latvec
  use mod_rotmat
  implicit none
  ! Passed parameters:
  integer :: nsymop
  real (kind=dp) :: platcp(3, 3), symopm(9, *)
  ! Local parameters:
  integer :: i, iprint, ltmax, ll1, m, m1, m2, m3, mm, nrot(4)
  parameter (ltmax=3, ll1=ltmax*2+1, iprint=20)
  real (kind=dp) :: platt(9), qlatcp(3, 3), mat(9), vecg(3), vol
  logical :: lirr

  data nrot/2, 3, 4, 6/

  mm(i, m) = ltmax - (mod(i,ll1**m)-mod(i,ll1**(m-1)))/ll1**(m-1)

  call dinv33(platcp, 1, qlatcp, vol)
  call rotmat(-1, .false., 1, symopm(1,1), [0.0_dp, 0.0_dp, 0.0_dp])
  call rotmat(-1, .true., 1, symopm(1,2), [0.0_dp, 0.0_dp, 0.0_dp])
  nsymop = 2
  ! --- find all possible rotation axis
  do i = 0, (ll1**3-1)/2 - 1
    m1 = mm(i, 1)
    m2 = mm(i, 2)
    m3 = mm(i, 3)
    lirr = .true.
    do m = 2, ll1
      lirr = lirr .and. (mod(m1,m)/=0 .or. mod(m2,m)/=0 .or. mod(m3,m)/=0)
    end do
    if (lirr) then
      do m = 1, 3
        vecg(m) = m1*platcp(m, 1) + m2*platcp(m, 2) + m3*platcp(m, 3)
      end do
      do m = 1, 4
        ! --------- create the matrix of the symmetry operation
        call rotmat(-1, .false., nrot(m), mat, vecg)
        call dmpy(mat, 3, 1, platcp, 3, 1, platt, 3, 1, 3, 3, 3)
        ! --------- check the primitive translations for the symmetry
        ! operations
        if (latvec(3,qlatcp,platt)) then
          call rotmat(-1, .false., nrot(m), symopm(1,nsymop+1), vecg)
          call rotmat(-1, .true., nrot(m), symopm(1,nsymop+2), vecg)
          nsymop = nsymop + 2
          if (m/=1) then
            call rotmat(-1, .false., -nrot(m), symopm(1,nsymop+1), vecg)
            call rotmat(-1, .true., -nrot(m), symopm(1,nsymop+2), vecg)
            nsymop = nsymop + 2
          end if
        end if
      end do
    end if
  end do
  if (iprint>=30) write (1337, 100) nsymop
100 format (/, ' SYMLAT: lattice invariant under ', i2, &
    ' symmetry operations.')
end subroutine symlat

end module mod_symlat
