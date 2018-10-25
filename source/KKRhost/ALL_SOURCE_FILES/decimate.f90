module mod_decimate

contains

  !-------------------------------------------------------------------------------
  !> Summary: Decimation method
  !> Author: 
  !> Category: KKRhost, structural-greensfunction, k-points
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Performs decimation step using t-matrices of left and right hosts to create structural GF
  !>
  !> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to
  ! Fortran90
  !-------------------------------------------------------------------------------
  subroutine decimate(gllke, naez, tinvbup, tinvbdown, vacflag, factl, nlbasis, nrbasis)

    use :: global_variables, only: lmmaxd, alm, ndim_slabinv, nprincd 
    use :: mod_datatypes, only: dp
    use :: mod_bofm, only: bofm
    use :: mod_surfgf, only: surfgf

    implicit none

    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: nlbasis !! Number of basis layers of left host (repeated units)
    integer, intent (in) :: nrbasis !! Number of basis layers of right host (repeated units)
    logical, dimension (2), intent (in) :: vacflag !! flag indicating if outside is vacuum or not
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: factl
    complex (kind=dp), dimension (lmmaxd, lmmaxd, *), intent (in) :: tinvbup !! left host t-matrix
    complex (kind=dp), dimension (lmmaxd, lmmaxd, *), intent (in) :: tinvbdown !! right host t-matrix
    complex (kind=dp), dimension (alm, alm), intent (inout) :: gllke
    ! .. Local Scalars
    integer :: ldi1, ldi1t, ldi2, ldi2t, lm1, lm2, nlayer, icall
    integer :: ichck, ihost, ii1, ii2, il1, il2, ip1, ip1t, ip2, ip2t, itermax
    real (kind=dp) :: errmax
    ! .. Local Arrays
    complex (kind=dp), dimension (ndim_slabinv, ndim_slabinv) :: a1, an, b1, bn, c1, cn, x1, xn
    ! ..
    data icall/0/
    ! .. Save statement
    save :: icall, nlayer, itermax, errmax, ichck
    ! .. External Functions
    logical :: opt
    external :: opt

    icall = icall + 1
    ! ----------------------------------------------------------------------------
    if (icall==1) then
      nlayer = naez/nprincd
      ! Parameters for the "decimation" technique.
      itermax = 300
      errmax = 1.0e-180_dp
      ichck = 1
    end if
    ! ----------------------------------------------------------------------------
    if (.not. vacflag(1)) then
      ! -------------------------------------------------------------------------
      ! Get the matrix B1
      ! -------------------------------------------------------------------------
      call bofm(1, 1, b1, ndim_slabinv, gllke, alm)

      ! Now Subtract t-mat of left host
      do ip1 = 1, nprincd
        ihost = mod(ip1-1, nlbasis) + 1
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            il1 = lmmaxd*(ip1-1) + lm1
            il2 = lmmaxd*(ip1-1) + lm2
            b1(il1, il2) = (b1(il1,il2)-tinvbup(lm1,lm2,ihost))
          end do
        end do
      end do

      call bofm(1, 2, c1, ndim_slabinv, gllke, alm)
      call bofm(2, 1, a1, ndim_slabinv, gllke, alm)

      ! It performs the 'space decimation' iterative procedure.
      call surfgf(ndim_slabinv, a1, b1, c1, x1, itermax, errmax, ichck)
      ! Adds to the matrix GLLKE the elements that couples the
      ! interface to the two half-spaces.
      do ip1 = 1, nprincd
        do ip2 = 1, nprincd
          ii1 = ip1
          ii2 = ip2
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              ldi1 = lmmaxd*(ip1-1) + lm1
              il1 = lmmaxd*(ii1-1) + lm1
              ldi2 = lmmaxd*(ip2-1) + lm2
              il2 = lmmaxd*(ii2-1) + lm2
              gllke(il1, il2) = gllke(il1, il2) - x1(ldi1, ldi2)
            end do
          end do
        end do
      end do
    end if

    if (.not. vacflag(2)) then

      ! If 'ONEBULK' is activated then it calculates the xn decimated element
      ! from the x1 element: this is just in the case of equal bulks on the

      if (.not. opt('ONEBULK ')) then

        ! ----------------------------------------------------------------------
        ! Get the matrix BN
        ! ----------------------------------------------------------------------
        call bofm(nlayer, nlayer, bn, ndim_slabinv, gllke, alm)

        ! Now Substract t-mat right host
        ! Notes : the indexing is easier like that
        do ip1 = 1, nprincd
          ihost = nrbasis - mod(ip1, nrbasis)
          ihost = mod(ip1-1, nrbasis) + 1
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              il1 = lmmaxd*(ip1-1) + lm1
              il2 = lmmaxd*(ip1-1) + lm2
              bn(il1, il2) = (bn(il1,il2)-tinvbdown(lm1,lm2,ihost))
            end do
          end do
        end do

        call bofm(nlayer, nlayer-1, an, ndim_slabinv, gllke, alm)
        call bofm(nlayer-1, nlayer, cn, ndim_slabinv, gllke, alm)

        ! It performs the 'space decimation' iterative procedure.
        call surfgf(ndim_slabinv, cn, bn, an, xn, itermax, errmax, ichck)

      else

        do ip1 = 1, nprincd
          do ip2 = 1, nprincd
            ip1t = (nprincd+1) - ip2
            ip2t = (nprincd+1) - ip1
            do lm1 = 1, lmmaxd
              do lm2 = 1, lmmaxd
                ldi1 = lmmaxd*(ip1-1) + lm1
                ldi2 = lmmaxd*(ip2-1) + lm2
                ldi1t = lmmaxd*(ip1t-1) + lm2
                ldi2t = lmmaxd*(ip2t-1) + lm1
                xn(ldi1t, ldi2t) = factl(lm1, lm2)*x1(ldi1, ldi2)
              end do
            end do
          end do
        end do
      end if
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Added on 1.02.2000
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Adds to the matrix GLLKE the elements that couples the
      ! interface to the two half-spaces.
      do ip1 = 1, nprincd
        do ip2 = 1, nprincd
          ii1 = (nlayer-1)*nprincd + ip1
          ii2 = (nlayer-1)*nprincd + ip2
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              ldi1 = lmmaxd*(ip1-1) + lm1
              il1 = lmmaxd*(ii1-1) + lm1
              ldi2 = lmmaxd*(ip2-1) + lm2
              il2 = lmmaxd*(ii2-1) + lm2
              gllke(il1, il2) = gllke(il1, il2) - xn(ldi1, ldi2)
            end do
          end do
        end do
      end do
    end if

    return

  end subroutine decimate

end module mod_decimate
