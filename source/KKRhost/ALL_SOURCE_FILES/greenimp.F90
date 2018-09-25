module mod_greenimp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate impurity GF by solving dyson equation
  !> Author: N. H. Long
  !> Date: 05.2013
  !> Category: KKRhost, input-output
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine greenimp(natomimp, dtmtrx, e)
    use :: mod_version_info
    use :: global_variables
    use :: mod_datatypes, only: dp

    implicit none

    complex (kind=dp) :: e, e1
    complex (kind=dp) :: cone, czero
    parameter (cone=(1d0,0d0), czero=(0d0,0d0))
    integer :: natomimp, ndim, info
    integer :: i, j, lm1, lm2, ilm, jlm, ilm1, jlm1
    integer :: ipvt1(natomimp*lmmaxso)
    complex (kind=dp) :: dtmtrx(lmmaxso*natomimp, lmmaxso*natomimp)
    complex (kind=dp), allocatable :: gimp(:, :), gi(:, :)

    ! read in GF of the host
    allocate (gimp(natomimp*lmmaxso,natomimp*lmmaxso))

    gimp = czero
    write (6, *) 'read in Green for host'
    read (60, '(2e17.9)') e1
    ! READ(60,*) E1
    do j = 1, natomimp
      do lm2 = 1, lmmaxso
        jlm = (j-1)*lmmaxso + lm2
        do i = 1, natomimp
          do lm1 = 1, lmmaxso
            ilm = (i-1)*lmmaxso + lm1
            ! READ(60,*) JLM1,ILM1,GIMP(ILM,JLM)
            read (60, '((2I5),(2e17.9))') jlm1, ilm1, gimp(ilm, jlm)
          end do
        end do
      end do
    end do

    ! calculate impurity GF
    allocate (gi(natomimp*lmmaxso,natomimp*lmmaxso))
    gi = czero
    ndim = natomimp*lmmaxso
    ! -G_host * delta t
    call zgemm('N', 'N', ndim, ndim, ndim, -cone, gimp, ndim, dtmtrx, ndim, czero, gi, ndim)
    do i = 1, ndim
      gi(i, i) = cone + gi(i, i)
    end do

    ! solve linear equation
    call zgetrf(ndim, ndim, gi, ndim, ipvt1, info)
    call zgetrs('N', ndim, ndim, gi, ndim, ipvt1, gimp, ndim, info)

    ! write down to the file GMATLL_GES
    write (59, '(2(e17.9,1X))') e
    do lm1 = 1, ndim
      do lm2 = 1, ndim
        write (59, '((2I5),(2e17.9))') lm2, lm1, gimp(lm2, lm1)
      end do
    end do
    deallocate (gimp)
    deallocate (gi)
  end subroutine greenimp

end module mod_greenimp
