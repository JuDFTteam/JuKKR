!< constructs potential including big/small components and with relativistic mass terms etc included
subroutine vllmatsra(vll0, vll, rmesh, lmsize, nrmax, nrmaxd, eryd, lmax, &
  lval_in, cmode)

  use :: constants
  use :: mod_datatypes, only: dp
  ! ************************************************************************************
  ! The perturbation matrix for the SRA-equations are set up
  ! ************************************************************************************
  implicit none

  ! inputs
  integer, intent (in) :: lmax     ! < Maximum l component in wave function expansion
  integer, intent (in) :: nrmax    ! < NTOTD*(NCHEBD+1)
  integer, intent (in) :: nrmaxd   ! < dimension for rmesh (maximum of nrmax values of all atoms)
  integer, intent (in) :: lmsize   ! < (lmax+2)^2
  integer, intent (in) :: lval_in  ! < l-value used in spherical calculation of calcpsh (lmsize=1)
  complex (kind=dp), intent (in) :: eryd  ! < energy (used in rel-mass factor)
  character (len=*), intent (in) :: cmode ! < either 'Ref=0' or 'Ref=Vsph' which determines the used reference system (for calcsph trick or direct evaluation)
  real (kind=dp), dimension (nrmaxd), intent (in) :: rmesh ! < radial mesh
  complex (kind=dp), dimension (lmsize, lmsize, nrmax), intent (in) :: vll0 ! < input potential in (l,m,s) basis

  ! outputs
  complex (kind=dp), dimension (2*lmsize, 2*lmsize, nrmax), intent (out) :: vll ! < output potential in (l,m,s) basis with big/small components

  ! locals
  integer :: ilm, lval, mval, ival, ir
  integer, dimension (lmsize) :: loflm
  complex (kind=dp) :: mass, mass0

  ! ************************************************************************************
  ! determine the bounds of the matricies to get the lm-expansion and the max.
  ! number of radial points
  ! ************************************************************************************

  ! ************************************************************************************
  ! calculate the index array to determine the L value of an LM index
  ! in case of spin-orbit coupling 2*(LMAX+1)**2 are used instead of
  ! (LMAX+1)**2
  ! the second half refers to the second spin and has the the same L value
  ! ************************************************************************************
  ilm = 0

  if (lmsize==1) then
    loflm(1) = lval_in
  else if ((lmax+1)**2==lmsize) then
    do lval = 0, lmax
      do mval = -lval, lval
        ilm = ilm + 1
        loflm(ilm) = lval
      end do
    end do
  else if (2*(lmax+1)**2==lmsize) then
    do ival = 1, 2
      do lval = 0, lmax
        do mval = -lval, lval
          ilm = ilm + 1
          loflm(ilm) = lval
        end do
      end do
    end do
  else
    stop '[vllmatsra] error'
  end if

  vll(:,:,:) = czero

  if (cmode=='Ref=0') then
    ! without SRA-trick the full matrix is constructed and the proper mass terms (as described in the PhD thesis of Bauer) have to be included
    vll(1:lmsize, 1:lmsize, :) = vll0(1:lmsize,1:lmsize,:) ! /cvlight

    do ir = 1, nrmax
      do ival = 1, lmsize
        lval = loflm(ival)
        mass = cone + (eryd-vll0(ival,ival,ir))/cvlight**2
        mass0 = cone + eryd/cvlight**2

        ! ************************************************************************************
        ! Conventional potential matrix
        ! ************************************************************************************

        vll(lmsize+ival, lmsize+ival, ir) = -vll0(ival, ival, ir)/cvlight**2 ! TEST 9/22/2011
        vll(ival, ival, ir) = vll(ival, ival, ir) + (1.0e0_dp/mass-1.0e0_dp/mass0)*lval*(lval+1)/rmesh(ir)**2

        ! ************************************************************************************
        ! The pertubation matrix is changed in the following way

        ! from  / V11  V12 \   to    / V21  V22 \
        !       \ V21  V22 /         \-V11 -V12 /
        ! because of the convention used for the left solution
        ! ************************************************************************************
      end do                       ! ival

    end do                         ! ir
  else if (cmode=='Ref=Vsph') then
    ! for SRA-trick the small component of vll vanishes
    vll(lmsize+1:2*lmsize, 1:lmsize, :) = vll0(1:lmsize,1:lmsize, :)
  end if

end subroutine vllmatsra
