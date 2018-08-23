module mod_madelung2d

contains

subroutine madelung2d(lpot, yrg, wg, naez, alat, vol, bravais, recbv, rbasis, &
  rmax, gmax, nlbasis, nleft, zperleft, tleft, nrbasis, nright, zperight, &
  tright, lmxspd, lassld, lpotd, lmpotd, nmaxd, ishld, nembd1, wlength)
  ! **********************************************************************
  ! *                                                                    *
  ! * This subroutine calculates the Madelung potential coefficients     *
  ! * in the 2D case and stores them in the DA-file abmad.unformatted    *
  ! * For each layer in the slab, the summation is split into three      *
  ! * parts (see also VINTERFACE):                                       *
  ! * within the slab, over the NLEFT*NLBASIS left host sites and over   *
  ! * the NRIGHT*NRBASIS right host sites, the last two steps only in    *
  ! * case of decimation run                                             *
  ! *                                                                    *
  ! * all positions must be scaled with ALAT to get them correct         *
  ! * (done in EWALD2D)                                                  *
  ! *                                                                    *
  ! * The record index is:                                               *
  ! *   (IQ1-1)*NAEZ + IQ2                 for (IQ1,IQ2) within the slab *
  ! *   NAEZ*NAEZ + (IQ1-1)*NLEFT*NLBASIS  for (IQ1,(IL,IBL)), IQ1 in    *
  ! *                  + (IL-1)*NLEFT+IBL  slab, (IL,IBL) in the left    *
  ! *   NAEZ*NAEZ + NAEZ*NLEFT*NLBASIS                                   *
  ! *             + (IQ1-1)*NRIGHT*NRBASIS for (IQ1,(IR,IBR)), IQ1 in    *
  ! *             + (IR-1)*NRIGHT+IBR      slab, (IR,IBR) in the right   *
  ! *                                                                    *
  ! **********************************************************************
  use :: mod_datatypes, only: dp
   use mod_madelgaunt
   use mod_madelcoef
   use mod_madel2out
   use mod_ewald2d
   use mod_lattice2d
  implicit none
  ! ..
  ! .. Scalar Arguments ..
  integer :: lpot, naez, wlength
  integer :: nlbasis, nleft, nrbasis, nright
  integer :: lassld, lpotd, lmpotd, lmxspd, nmaxd, ishld, nembd1
  real (kind=dp) :: alat, vol, rmax, gmax
  ! ..
  ! .. Array Arguments ..
  real (kind=dp) :: bravais(3, 3), recbv(3, 3)
  real (kind=dp) :: yrg(lassld, 0:lassld, 0:lassld), wg(lassld)
  real (kind=dp) :: rbasis(3, *)
  real (kind=dp) :: zperight(3), zperleft(3)
  real (kind=dp) :: tleft(3, nembd1), tright(3, nembd1)
  ! ..
  ! .. Local Scalars ..
  integer :: iq1, iq2, iend, nclebd, iprint
  integer :: i, ib, ih, ileft, iright
  integer :: lrecamad, irec, nleftoff, nrightoff, nleftall, nrightall
  integer :: ngmax, nrmax, nshlg, nshlr
  logical :: opt
  ! ..
  ! .. Local Arrays ..
  ! .. Attention: LMXSPD*LMPOTD appears as NCLEB1 in other routines
  real (kind=dp) :: cleb(lmxspd*lmpotd)
  real (kind=dp) :: bm(lmpotd), vec2(3), sum(lmxspd)
  real (kind=dp) :: gn2(2, nmaxd), rm2(2, nmaxd)
  real (kind=dp) :: avmad(lmpotd, lmpotd)
  integer :: nsg(ishld), nsr(ishld)
  integer :: icleb(lmxspd*lmpotd, 3)
  ! ..
  ! .. External Functions/Subroutines
  ! ......................................................................
  iprint = 0
  nclebd = lmxspd*lmpotd

  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  write (1337, '(79("="))')
  write (1337, '(18X,A)') 'MADELUNG2D: setting 2D Madelung coefficients'
  write (1337, '(79("="))')
  write (1337, *)
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  ! ======================================================================
  call lattice2d(alat, bravais, recbv, ngmax, nrmax, nshlg, nshlr, nsg, nsr, &
    gn2, rm2, rmax, gmax, iprint, nmaxd, ishld)
  ! ======================================================================

  lrecamad = wlength*2*lmpotd*lmpotd
  open (69, access='direct', recl=lrecamad, file='avmad.unformatted', &
    form='unformatted')

  ! --> calculate the gaunt coefs

  call madelgaunt(lpot, yrg, wg, cleb, icleb, iend, lassld, nclebd)

  ! --> calculate the madelung coefficients to be used for VMAD

  ! **********************************************************************
  ! ********************************************** loop over atoms in slab
  do iq1 = 1, naez

    ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    ! 1.  Summation in all layers in the slab
    ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    ! ++++++++++++++++++++++++++++++++ loop over all other sites in the slab
    do iq2 = 1, naez

      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      if (iq1==1 .and. iq2==1) then
        write (1337, '(5X,2A,/)') '< EWALD2D > : calculating 2D-lattice sums ' &
          , 'inside the slab'
        if (iprint>=2) write (1337, 100)
      end if
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT


      ! make ewald sumation in plane and inverse space
      ! sum if rz<>0 (out of plane)

      ! WRITE(99,*) 'Layer pair:',IQ1,IQ2
      call ewald2d(lpot, alat, rbasis(1,iq1), rbasis(1,iq2), iq1, iq2, rm2, &
        nrmax, nshlr, nsr, gn2, ngmax, nshlg, nsg, sum, vol, lassld, lmxspd)
      ! WRITE(99,*) 'SUM: ',SUM

      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      if (iprint>=2) then
        write (1337, 110) iq1, iq2, sum(1)
        if (iq2==naez .and. iq1/=naez) write (1337, '(20X,20("-"))')
      end if
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      call madelcoef(.true., lpot, avmad, bm, sum, cleb, icleb, iend, lpotd, &
        lmpotd, lmxspd, nclebd)

      irec = iq2 + naez*(iq1-1)
      write (69, rec=irec) avmad
    end do
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  if (iprint>=2) write (1337, '(18X,22("-"),/)')
  ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
  ! ********************************************** loop over atoms in slab

  ! ######################################################################
  if (opt('DECIMATE')) then

    nleftoff = naez*naez           ! record offsets
    nrightoff = nleftoff + naez*nleft*nlbasis ! left and right
    nleftall = nleft*nlbasis
    nrightall = nright*nrbasis

    ! ********************************************** loop over atoms in slab
    do iq1 = 1, naez
      ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      ! 2.  Summation in the LEFT bulk side
      ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      ileft = 0
      ! ++++++++++++++++++++++++++++++++ loop over all sites in the left host
      do ih = 1, nleft
        do ib = 1, nlbasis
          do i = 1, 3
            vec2(i) = (tleft(i,ib)+(ih-1)*zperleft(i))
          end do
          ileft = ileft + 1

          ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
          ! OUTPUT
          if (iq1==1 .and. ileft==1) then
            write (1337, '(5X,2A,/)') &
              '< EWALD2D > : calculating 2D-lattice sums ', 'slab - left host'
            if (iprint>=2) write (1337, 100)
          end if
          ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
          ! OUTPUT


          ! -->  make ewald sumation for m= 0 l<5 rz=0 (in plane) and
          ! Inverse space sum if rz<>0 (out of plane)

          call ewald2d(lpot, alat, rbasis(1,iq1), vec2, iq1, ih, rm2, nrmax, &
            nshlr, nsr, gn2, ngmax, nshlg, nsg, sum, vol, lassld, lmxspd)

          ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
          ! OUTPUT
          if (iprint>=2) then
            write (1337, 110) iq1, ileft, sum(1)
            if (ileft==nleftall .and. iq1/=naez) write (1337, '(20X,20("-"))')
          end if
          ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
          ! OUTPUT

          call madelcoef(.true., lpot, avmad, bm, sum, cleb, icleb, iend, &
            lpotd, lmpotd, lmxspd, nclebd)

          irec = ileft + nleftall*(iq1-1) + nleftoff
          write (69, rec=irec) avmad
        end do                     ! ib loop in left host basis
      end do                       ! ih loop in layers to get convergence
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (ileft/=nleftall) then
        write (6, *) ' < MADELUNG2D > : index error ', &
          'ILEFT <> NLEFT*NLBASIS'
        stop
      end if
    end do                         ! ILAY1 loop
    ! **********************************************************************

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    if (iprint>=2) write (1337, '(18X,22("-"),/)')
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

    ! ********************************************** loop over atoms in slab
    do iq1 = 1, naez
      ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      ! 3.  Summation in the RIGHT bulk side
      ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      ! ++++++++++++++++++++++++++++++++ loop over all sites in the right host
      iright = 0
      do ih = 1, nright
        do ib = 1, nrbasis
          do i = 1, 3
            vec2(i) = (tright(i,ib)+(ih-1)*zperight(i))
          end do
          iright = iright + 1

          ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
          ! OUTPUT
          if (iq1==1 .and. iright==1) then
            write (1337, '(5X,2A,/)') &
              '< EWALD2D > : calculating 2D-lattice sums ', &
              'slab - right host'
            if (iprint>=2) write (1337, 100)
          end if
          ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
          ! OUTPUT

          ! -->  make ewald sumation (in plane) and
          ! Inverse space sum if rz<>0 (out of plane)

          call ewald2d(lpot, alat, rbasis(1,iq1), vec2, iq1, ih, rm2, nrmax, &
            nshlr, nsr, gn2, ngmax, nshlg, nsg, sum, vol, lassld, lmxspd)

          ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
          ! OUTPUT
          if (iprint>=2) then
            write (1337, 110) iq1, iright, sum(1)
            if (iright==nrightall .and. iq1/=naez) write (1337, &
              '(20X,20("-"))')
          end if
          ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
          ! OUTPUT

          call madelcoef(.true., lpot, avmad, bm, sum, cleb, icleb, iend, &
            lpotd, lmpotd, lmxspd, nclebd)

          irec = iright + nrightall*(iq1-1) + nrightoff
          write (69, rec=irec) avmad
        end do                     ! ib loop in right host basis
      end do                       ! ih loop in layers to get convergence
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (iright/=nrightall) then
        write (6, *) ' < MADELUNG2D > : index error ', &
          'IRIGHT <> NRIGHT*NRBASIS'
        stop
      end if
    end do                         ! ILAY1 loop
    ! **********************************************************************

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    if (iprint>=2) write (1337, '(18X,22("-"),/)')
    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

  end if
  ! ######################################################################
  close (69)

  if (iprint<1) return
  ! ======================================================================

  call madel2out(iprint, naez, lrecamad, lmpotd, nleftoff, nrightoff, &
    nleftall, nrightall)

100 format (8x, '2D Lattice sum (LMXSP = 1)', /, 18x, '  IQ1  IQ2  SUM', /, &
    18x, 23('-'))
110 format (18x, 2i5, d12.4)
end subroutine madelung2d

end module mod_madelung2d
