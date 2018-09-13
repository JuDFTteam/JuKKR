module mod_rites

contains

! -------------------------------------------------------------------------------
! SUBROUTINE: RITES
!> @brief this subroutine stores in 'ifile' the necessary results
!> (potentials etc.) to start self-consistency iterations

!> @ details Modified for the full potential case - if ins .gt. 0 there
!> is written a different potential card
!> if the sum of absolute values of an lm component of vins (non
!> spher. potential) is less than the given rms error qbound this
!> component will not be stored .

!> see to subroutine start , where most of the arrays are described)

!> @note modified by B. Drittler  aug. 1988
!> @note Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to
! Fortran90
! -------------------------------------------------------------------------------
subroutine rites(ifile, natps, natyp, nspin, z, alat, rmt, rmtnew, rws, &
  ititle, r, drdi, vm2z, irws, a, b, txc, kxc, ins, irns, lpot, vins, qbound, &
  irc, kshape, efermi, vbc, ecore, lcore, ncore, ecorerel, nkcore, kapcore, &
  lmpot)

  use :: global_variables
  use :: mod_datatypes, only: dp

  ! .. Scalar Arguments
  integer, intent (in) :: ins      !! 0 (MT), 1(ASA), 2(Full Potential)
  integer, intent (in) :: kxc      !! Type of xc-potential 0=vBH 1=MJW 2=VWN
                                   ! 3=PW91
  integer, intent (in) :: lpot     !! Maximum l component in potential
                                   ! expansion
  integer, intent (in) :: lmpot    !! (LPOT+1)**2
  integer, intent (in) :: ifile    !! Unit specifier for potential card
  integer, intent (in) :: natps
  integer, intent (in) :: natyp    !! Number of kinds of atoms in unit cell
  integer, intent (in) :: nspin    !! Counter for spin directions
  integer, intent (in) :: kshape   !! Exact treatment of WS cell
  real (kind=dp), intent (in) :: alat !! Lattice constant in a.u.
  real (kind=dp), intent (in) :: qbound !! Convergence parameter for the
                                        ! potential
  real (kind=dp), intent (in) :: efermi !! Fermi energy
  ! .. Array Arguments
  real (kind=dp), dimension (*), intent (in) :: a !! Constants for
                                                  ! exponential R mesh
  real (kind=dp), dimension (*), intent (in) :: b !! Constants for
                                                  ! exponential R mesh
  real (kind=dp), dimension (*), intent (in) :: z
  real (kind=dp), dimension (*), intent (in) :: rws !! Wigner Seitz radius
  real (kind=dp), dimension (2), intent (in) :: vbc !! Potential constants
  real (kind=dp), dimension (*), intent (in) :: rmt !! Muffin-tin radius of
                                                    ! true system
  real (kind=dp), dimension (*), intent (in) :: rmtnew !! Adapted muffin-tin
                                                       ! radius
  real (kind=dp), dimension (irmd, *), intent (in) :: r !! Radial mesh ( in
                                                        ! units a Bohr)
  real (kind=dp), dimension (irmd, *), intent (in) :: vm2z
  real (kind=dp), dimension (irmd, *), intent (in) :: drdi !! Derivative
                                                           ! dr/di
  real (kind=dp), dimension (20, *), intent (in) :: ecore !! Core energies
                                                          ! !(2), 22.5,2000
  real (kind=dp), dimension (irmind:irmd, lmpot, *), intent (in) :: vins !!
                                                                         ! Non-spherical
                                                                         ! part
                                                                         ! of
                                                                         ! the
                                                                         ! potential
  ! ----------------------------------------------------------------------------
  integer, dimension (20, natyp), intent (in) :: nkcore
  integer, dimension (20, 2*natyp), intent (in) :: kapcore
  real (kind=dp), dimension (krel*20+(1-krel), 2*natyp), &
    intent (in) :: ecorerel        ! relativistic core energies
  ! ----------------------------------------------------------------------------
  integer, dimension (*), intent (in) :: irc !! R point for potential cutting
  integer, dimension (*), intent (in) :: irns !! Position of atoms in the
                                              ! unit cell in units of bravais
                                              ! vectors
  integer, dimension (*), intent (in) :: irws !! R point at WS radius
  integer, dimension (*), intent (in) :: ncore !! Number of core states
  integer, dimension (20, *), intent (in) :: ititle
  integer, dimension (20, *), intent (in) :: lcore !! Angular momentum of
                                                   ! core states
  character (len=124), dimension (*), intent (in) :: txc
  ! .. Local Scalars
  integer :: i, icore, ih, inew, ip, ir, irmin, irns1, is, isave, j, lm, lmnr, &
    ncore1, nr
  real (kind=dp) :: a1, b1, rmax, rmt1, rmtnw1, rv, sign, sum, z1
  ! .. Local Arrays
  integer, dimension (20) :: lcore1
  real (kind=dp), dimension (20) :: ecore1
  real (kind=dp), dimension (irmd) :: dradi
  real (kind=dp), dimension (irmd) :: ra
  real (kind=dp), dimension (irmd) :: vm2za
  real (kind=dp), dimension (20, 2) :: ecore2
  character (len=3), dimension (4) :: txtk
  character (len=1), dimension (0:3) :: txtl
  ! ..
  data txtl/'s', 'p', 'd', 'f'/
  data txtk/'1/2', '3/2', '5/2', '7/2'/
  ! ..
  ! -------------------------------------------------------------------
  isave = 1
  inew = 1


  do ih = 1, natyp
    do is = 1, nspin
      if (is==nspin) then
        sign = 1.0e0_dp
      else
        sign = -1.0e0_dp
      end if
      ip = nspin*(ih-1) + is

      rmt1 = rmt(ih)
      rmtnw1 = rmtnew(ih)
      z1 = z(ih)
      rmax = rws(ih)
      if (kshape==0) then
        nr = irws(ih)
      else
        nr = irc(ih)
      end if

      irns1 = irns(ih)
      irmin = nr - irns1
      a1 = a(ih)
      b1 = b(ih)
      ncore1 = ncore(ip)

      do j = 1, nr
        ra(j) = r(j, ih)
        dradi(j) = drdi(j, ih)
        ! -------------------------------------------------------------------
        ! Store only lm=1 component of the potential
        ! -------------------------------------------------------------------
        vm2za(j) = vm2z(j, ip)
      end do                       ! J

      open (ifile, file='out_potential', form='formatted')
      write (ifile, fmt=100)(ititle(i,ip), i=1, 7), txc(kxc+1)
      write (ifile, fmt=110) rmt1, alat, rmtnw1
      write (ifile, fmt=120) z1, rmax, efermi, vbc(is)
      write (ifile, fmt=130) nr, a1, b1, ncore1, inew

      if (ncore1>=1) then

        if (krel==0) then
          do j = 1, ncore1
            lcore1(j) = lcore(j, ip)
            ecore1(j) = ecore(j, ip)
          end do
          write (ifile, fmt=140)(lcore1(icore), ecore1(icore), icore=1, &
            ncore1)
        else
          do j = 1, ncore1
            lcore1(j) = lcore(j, ip)
            ecore2(j, 1) = ecorerel(j, 2*ih-1)
            ecore2(j, 2) = ecorerel(j, 2*ih)
          end do
          ! ----------------------------------------------------------------
          ! independent of spin, the \mu-averaged relativistic core energies
          ! are written out for \kappa = -l-1,l
          ! format compatible with the non-(scalar) relativistic mode
          ! however, the next read in has no meaning for the REL core-solver
          ! a detailed output of the core energies is supplied by < CORE >
          ! ----------------------------------------------------------------
          do icore = 1, ncore1
            write (ifile, fmt=150) lcore1(icore), (ecore2(icore,i+1), txtl( &
              lcore1(icore)), txtk(abs(kapcore(icore,2*ih-1+i))), i=0, &
              nkcore(icore,ih)-1)
          end do
        end if

      end if

      if (ins==0 .or. (ih<natps .and. ins<=2)) then
        ! -------------------------------------------------------------------
        ! store only the spherically averaged potential
        ! (in mt or as - case)
        ! this is done always for the host
        ! -------------------------------------------------------------------
        if (inew==0) then
          write (ifile, fmt=160)(ra(ir), dradi(ir), vm2za(ir), ir=1, nr)
        else
          write (ifile, fmt=170)(vm2za(ir), ir=1, nr)
        end if
      else
        ! -------------------------------------------------------------------
        ! store the full potential , but the non spherical contribution
        ! only from irns1 up to irws1 ;
        ! remember that the lm = 1 contribution is multiplied
        ! by a factor 1/sqrt(4 pi)
        ! -------------------------------------------------------------------
        write (ifile, fmt=180) nr, irns1, lmpot, isave
        write (ifile, fmt=190)(vm2za(ir), ir=1, nr)
        if (lpot>0) then
          lmnr = 1
          do lm = 2, lmpot
            sum = 0.0e0_dp
            do ir = irmin, nr
              rv = vins(ir, lm, ip)*ra(ir)
              sum = sum + rv*rv*dradi(ir)
            end do                 ! IR
            if (sqrt(sum)>qbound) then
              lmnr = lmnr + 1
              write (ifile, fmt=180) lm
              write (ifile, fmt=190)(vins(ir,lm,ip), ir=irmin, nr)
            end if
          end do                   ! LM
          ! ----------------------------------------------------------------
          ! Write a one to mark the end
          ! ----------------------------------------------------------------
          if (lmnr<lmpot) write (ifile, fmt=180) isave
        end if
      end if
    end do                         ! IS
  end do                           ! IH

  close (ifile)

100 format (7a4, 6x, '  exc:', a124, 3x, a10)
110 format (3f20.15)
  ! 9010 FORMAT (3F) !f12.8) maybe change to higher accuracy in writeout?
120 format (f10.5, /, f10.5, 2f20.15)
130 format (i3, /, 2d15.8, /, 2i2)
140 format (i5, 1p, d20.11)
150 format (i5, 2(1p,d20.11,2x,a1,a3))
160 format (1p, 2d15.6, 1p, d15.8)
170 format (1p, 4d20.12)
180 format (10i5)
190 format (1p, 4d20.13)
end subroutine rites

end module mod_rites
