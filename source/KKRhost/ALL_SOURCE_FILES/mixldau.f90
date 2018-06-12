subroutine mixldau(mmaxd, nspind, natypd, natyp, nspin, lopt, wldauold, wldau)
  use :: mod_datatypes, only: dp
  implicit none
  ! Input:
  integer :: natypd, nspind, mmaxd
  integer :: lopt(natypd)
  integer :: natyp, nspin
  real (kind=dp) :: wldauold(mmaxd, mmaxd, nspind, natypd)
  ! Input/Output:
  real (kind=dp) :: wldau(mmaxd, mmaxd, nspind, natypd)
  ! Inside:
  integer :: iat, is, m1, m2, mmax
  integer :: ier
  real (kind=dp) :: xmix, xmix2, rmserr
  character (len=256) :: uio                               ! NCOLIO=256



  external :: ioinput


  ! First calculate rms error in interaction matrix

  do iat = 1, natyp
    rmserr = 0.e0_dp
    if (lopt(iat)>=0) then
      mmax = 2*lopt(iat) + 1
      do is = 1, nspin
        do m2 = 1, mmax
          do m1 = 1, mmax
            rmserr = rmserr + (wldau(m1,m2,is,iat)-wldauold(m1,m2,is,iat))**2
          end do
        end do
      end do
      rmserr = sqrt(rmserr)
      write (1337, 100) iat, rmserr
100   format ('LDA+U interaction matrix rms error for atom', i6, ' = ', e10.2)
    end if
  end do

  ! Now mix old/new interaction matrices
  ier = 0
  call ioinput('MIXLDAU         ', uio, 1, 7, ier)
  if (ier/=0) then
    write (*, *) 'MIXLDAU not found, setting to 1.'
    return
  else
    read (unit=uio, fmt=*) xmix
    write (1337, *) 'Using MIXLDAU = ', xmix
  end if

  xmix2 = 1.e0_dp - xmix

  do iat = 1, natyp
    if (lopt(iat)>=0) then
      do is = 1, nspin
        do m2 = 1, mmaxd
          do m1 = 1, mmaxd
            wldau(m1, m2, is, iat) = xmix*wldau(m1, m2, is, iat) + &
              xmix2*wldauold(m1, m2, is, iat)
          end do
        end do
      end do
    end if
  end do

end subroutine mixldau
