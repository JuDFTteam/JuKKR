subroutine drvbastrans(rc, crel, rrel, srrel, nrrel, irrel, nlmax, nkmmax, &
  nmuemax, nkmpmax, nkmax, linmax)
  use :: mod_datatypes, only: dp
  ! ********************************************************************
  ! *                                                                  *
  ! *                                                                  *
  ! ********************************************************************
  implicit none

  ! Dummy arguments
  integer :: linmax, nkmax, nkmmax, nkmpmax, nlmax, nmuemax
  complex (kind=dp) :: crel(nkmmax, nkmmax), rc(nkmmax, nkmmax), &
    rrel(nkmmax, nkmmax), srrel(2, 2, nkmmax)
  integer :: irrel(2, 2, nkmmax), nrrel(2, nkmmax)

  ! Local variables
  real (kind=dp) :: cgc(nkmpmax, 2)
  integer :: i, ikm1lin(linmax), ikm2lin(linmax), il, imue, iprint, &
    kaptab(nmuemax), ltab(nmuemax), mmax, nmuetab(nmuemax), &
    nsollm(nlmax, nmuemax)

  if (nkmmax/=2*nlmax**2) stop ' Check NLMAX,NKMMAX in < DRVBASTRANS > '
  if (nmuemax/=2*nlmax) stop ' Check NLMAX,NMUEMAX in < DRVBASTRANS > '
  if (nkmpmax/=(nkmmax+2*nlmax)) stop &
    ' Check NLMAX,NKMMAX,NKMPMAX in < DRVBASTRANS > '
  if (nkmax/=2*nlmax-1) stop ' Check NLMAX,NKMAX in < DRVBASTRANS > '
  if (linmax/=(2*nlmax*(2*nlmax-1))) stop &
    ' Check NLMAX,LINMAX in < DRVBASTRANS > '

  iprint = 0

  do i = 1, nmuemax
    ltab(i) = i/2
    if (2*ltab(i)==i) then
      kaptab(i) = ltab(i)
    else
      kaptab(i) = -ltab(i) - 1
    end if
    nmuetab(i) = 2*abs(kaptab(i))
  end do

  do il = 1, nlmax
    mmax = 2*il
    do imue = 1, mmax
      if ((imue==1) .or. (imue==mmax)) then
        nsollm(il, imue) = 1
      else
        nsollm(il, imue) = 2
      end if
    end do
  end do

  call ikmlin(iprint, nsollm, ikm1lin, ikm2lin, nlmax, nmuemax, linmax, nlmax)

  call calccgc(ltab, kaptab, nmuetab, cgc, nkmax, nmuemax, nkmpmax)

  ! ---------------------------- now calculate the transformation matrices

  call strsmat(nlmax-1, cgc, srrel, nrrel, irrel, nkmmax, nkmpmax)

  call bastrmat(nlmax-1, cgc, rc, crel, rrel, nkmmax, nkmpmax)

  return
end subroutine drvbastrans
