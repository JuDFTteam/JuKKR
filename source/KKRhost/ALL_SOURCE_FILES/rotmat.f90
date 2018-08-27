module mod_rotmat

contains

subroutine rotmat(iopt, li, nrot, symopm, vecg)
  ! - Converts rotation/rotoinversion matrix <-> (nrot,vecg,li)
  ! ----------------------------------------------------------------------
  ! i Inputs:
  ! i   iopt  := -1 to convert (nrot,vecg,li) in symopm
  ! i          =  1 to convert symopm in (nrot,vecg,li)
  ! i Inputs/Outputs:
  ! io  li    :if T: inversion or rotoinversion
  ! io  nrot  :rotation angle = 2*pi/nrot
  ! io  symopm:symmetry operation matrix
  ! io  vecg  :rotation axis
  ! ----------------------------------------------------------------------
  use :: mod_datatypes, only: dp
  use mod_ddet33
  use mod_errmsg
  use mod_nrmliz
  use mod_rinit
  implicit none
  ! Passed parameters:
  integer :: iopt, nrot
  real (kind=dp) :: vecg(3), symopm(3, 3)
  logical :: li
  ! Local parameters:
  integer :: i, idamax, in, j
  real (kind=dp) :: costbn, detop, dnrm2, omcos, sintbn, sinpb3, tiny, &
    twopi, vfac
  character (len=144) :: messg
  parameter (twopi=6.28318530717958648e0_dp)
  parameter (tiny=1.0e-3_dp)

  if (iopt==-1) then
    call rinit(9, symopm)
    in = abs(nrot)
    if (in==1) then
      call dcopy(3, 1.e0_dp, 0, symopm, 4)
    else if (in==2 .or. in==3 .or. in==4 .or. in==6) then
      sintbn = sin(twopi/nrot)
      costbn = cos(twopi/nrot)
      omcos = 1.e0_dp - costbn
      if (dnrm2(3,vecg,1)<tiny) call errmsg(' ROTMAT: zero rotation vector.$', &
        4)
      call nrmliz(1, vecg, vecg)
      symopm(1, 1) = omcos*vecg(1)*vecg(1) + costbn
      symopm(1, 2) = omcos*vecg(1)*vecg(2) - sintbn*vecg(3)
      symopm(1, 3) = omcos*vecg(1)*vecg(3) + sintbn*vecg(2)
      symopm(2, 1) = omcos*vecg(2)*vecg(1) + sintbn*vecg(3)
      symopm(2, 2) = omcos*vecg(2)*vecg(2) + costbn
      symopm(2, 3) = omcos*vecg(2)*vecg(3) - sintbn*vecg(1)
      symopm(3, 1) = omcos*vecg(3)*vecg(1) - sintbn*vecg(2)
      symopm(3, 2) = omcos*vecg(3)*vecg(2) + sintbn*vecg(1)
      symopm(3, 3) = omcos*vecg(3)*vecg(3) + costbn
    else
      call errmsg(' ROTMAT: bad nrot.$', 3)
    end if
    if (li) call dscal(9, -1.e0_dp, symopm(1,1), 1)

  else if (iopt==1) then
    ! ----- First calculate determinant.
    detop = ddet33(symopm)
    if (abs(abs(detop)-1.0e0_dp)>tiny) call errmsg( &
      ' ROTMAT: determinant is not +/- 1$', 4)
    detop = sign(1.e0_dp, detop)
    li = detop < 0.e0_dp
    ! ----- multiply operation symopm with detop
    call dscal(9, detop, symopm(1,1), 1)
    ! ----- For the rotation angle we have due to the normalization of v:
    ! ----- sum_i symopm(i,i) = sum_i (1-cos) v_i*v_i+3*cos = 1 + 2 * cos,
    costbn = -0.5e0_dp
    call daxpy(3, 0.5e0_dp, symopm(1,1), 4, costbn, 0)
    if (abs(costbn-1.e0_dp)<tiny) then
      nrot = 1
      call rinit(3, vecg)
    else
      nrot = nint(twopi/acos(max(-1.e0_dp,costbn)))
      ! ------- for nrot > 2 the matrix is non-symmetric and the rotation
      ! ------- axis can be calculated from the antisymmetric part.
      ! ------- for nrot = 2 this not possible. However, the squared vector
      ! ------- components are given by:  mat(i,i) = 2 v_i * v_i - 1.
      ! ------- This is used for the largest component. The others are taken
      ! ------- from: mat(i,j) = 2 v_i * v_j for i ne j. This way we also
      ! ------- get the right phases between the components.
      if (nrot==2) then
        do i = 1, 3
          vecg(i) = 0.5e0_dp*(symopm(i,i)+1.0e0_dp)
        end do
        j = idamax(3, vecg, 1)
        if (vecg(j)<0.0e0_dp) then
          write (messg, 100) j, symopm(j, j)
          call errmsg(messg, 4)
        end if
        vecg(j) = sqrt(vecg(j))
        vfac = 0.5e0_dp/vecg(j)
        do i = 1, 3
          if (i/=j) vecg(i) = vfac*symopm(i, j)
        end do
      else
        vecg(1) = symopm(3, 2) - symopm(2, 3)
        vecg(2) = symopm(1, 3) - symopm(3, 1)
        vecg(3) = symopm(2, 1) - symopm(1, 2)
      end if
      ! ------- next renormalize at least one component to 1 in order to
      ! ------- allow for abbreviations as 'D', 'X', 'Y' or 'Z'
      sinpb3 = sqrt(.75e0_dp)
      if (abs((sinpb3-abs(vecg(1)))*(sinpb3-abs(vecg(2)))*(sinpb3- &
        abs(vecg(3))))>tiny) then
        do j = 3, 1, -1
          vfac = abs(vecg(j))
          if (vfac>tiny) call dscal(3, 1.e0_dp/vfac, vecg, 1)
        end do
      end if
    end if
    call dscal(9, detop, symopm(1,1), 1)
  end if

100 format (' ROTMAT: Bad component ', i1, ' of operation ', &
    '. Diagonal element =', f9.5, '$')
end subroutine rotmat

end module mod_rotmat
