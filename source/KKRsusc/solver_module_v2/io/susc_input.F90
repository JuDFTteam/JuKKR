  subroutine susc_input(lmax,natyp,nspin,nrad,itc)
! Handles the input file and sets main parameters
  use global

  implicit none

! --> values given in the main program
  integer(kind=i4b), intent(in) :: lmax, natyp, nspin, nrad, itc
! -----------------------------------------------------------------
  integer(kind=i4b) :: istart, iend, ng, nbtmp
  integer(kind=i4b) :: i, ie, is, ib, il, iw(0:lmax), ntot, itemp1, itemp2, itemp3
  real(kind=r8b)    :: ere1, eim1, ere2, eim2, utemp(0:lmax), jtemp(0:lmax), socs, bf, dir(3)
  real(kind=r8b)    :: er, ei, der, dei

  open(file='inpsusc.dat',unit=iomain,status='old')
  read(iomain,*) ! ---------------------
  read(iomain,*) ! Data for ...
  read(iomain,*) ! ---------------------
  read(iomain,*) ! ikkr, ldos, lsoc, lldau, lbfield, lsusc, lrot, lfit=
  read(iomain,*) ikkr, ldos, lsoc, lldau, lbfield, lsusc, lrot, lfit
  read(iomain,*) ! ibasis, method, basistol, freeze=
  read(iomain,*) ibasis, method, basistol, freeze
  read(iomain,*) ! idos, nedos, e0dos, e1dos, eimdos=
  read(iomain,*) idos, nedos, e0dos, e1dos, eimdos
  read(iomain,*) ! lrhomat, ldaumix, ldauitc=
  read(iomain,*) lrhomat, ldaumix, ldauitc
  read(iomain,*) ! ldynsusc, nomega, omegamin, omegamax, domega, nlmax0=
  read(iomain,*) ldynsusc, nomega, omegamin, omegamax, domega, nlmax0
  read(iomain,*) ! lkha, lkxc, cartesian, analytic, nonanalytic, enhanced, itermax, lambdamix=
  read(iomain,*) lkha, lkxc, cartesian, analytic, nonanalytic, enhanced, itermax, lambdamix
  read(iomain,*) ! ispinrot, urot(1:3), dirmix=
  read(iomain,*) ispinrot, urot(1:3), dirmix
  read(iomain,*) ! ifit, numd, dend, eshift, lregf, fudge=
  read(iomain,*) ifit, numd, dend, ere1, eim1, lregf, fudge
  eshift = cmplx(ere1,eim1)
  read(iomain,*) ! ---------------------
  read(iomain,*) ! nesusc, nasusc, nbmax, nlmax, sra, nrmax, nsmax=
  read(iomain,*) nesusc, nasusc, nbmax, nlmax, isra, nrmax, nsmax
!  if (ne /= nepts) stop 'nesusc /= nepts'   ! messy in fpimpu
  if (lmax /= nlmax)  stop 'susc_input: lmax /= nlmax'
  if (natyp < nasusc) stop 'susc_input: natyp < nasusc'
  if (nrad > nrmax)   stop 'susc_input: nrad > nrmax'
  if (nsmax /= nspin) stop 'susc_input: nspin /= nsmax'
! -----------------------------------------------------------------
! inconsistencies?
  call check_global_options(itc)
! -----------------------------------------------------------------
! initialize remaining parameters and allocate arrays needed to read input
  call init_param
! -----------------------------------------------------------------
! read atom labels in groups
! if a set of atoms is to be treated in the same way, only one read is needed
  read(iomain,*) ! ngroup, igroup(1:ngroup)=
  read(iomain,*) ngroup
  read(iomain,*) igroup(1:ngroup)
  ntot = sum(igroup(1:ngroup))
  if (ntot /= nasusc) stop 'susc_input: atom groups and nasusc'
! loop over atom groups
  iend = 0
! **************
  do ng=1,ngroup
! **************
    istart = iend + 1
    iend   = iend + igroup(ng)
    read(iomain,*) ! atomic labels, atom numbers in main program
    read(iomain,*) iasusc(istart:iend)
    read(iomain,*) ! isoc, socscaling, ildau, ueff(0:lmax), jeff(0:lmax)=
    read(iomain,*) itemp1, socs, itemp2, utemp(0:nlmax), jtemp(0:nlmax)
!   SOC options
    isoc(istart:iend) = itemp1
    if (lsoc .and. itemp1 /= 0) then
      socscaling(istart:iend) = socs
    else
      socscaling(istart:iend) = 0.d0
    end if
!   LDA+U options
    ildau(istart:iend) = itemp2
    if (lldau .and. itemp2 == 1) then
      do il=0,nlmax
        ueff(il,istart:iend) = utemp(il)
        jeff(il,istart:iend) = jtemp(il)
      end do
    else
      ueff(0:nlmax,istart:iend) = 0.d0
      jeff(0:nlmax,istart:iend) = 0.d0
    end if
    read(iomain,*) ! ibfield, blen, bdir=
    read(iomain,*) itemp1, bf, dir
!   external field options
    ibfield(istart:iend) = itemp1
    if (lbfield .and. itemp1 /= 0) then
      blen(istart:iend) = bf
      do i=1,3
        bdir(i,istart:iend) = dir(i)
      end do
    else
      blen(istart:iend) = 0.d0
      bdir(:,istart:iend) = 0.d0
    end if
    read(iomain,*) ! isusc, ikha, ikxc=
    read(iomain,*) itemp1, itemp2, itemp3
!   Susceptibility and kernel options
    if (lsusc) isusc(istart:iend) = itemp1
    if (lkha)  ikha(istart:iend)  = itemp2
    if (lkxc)  ikxc(istart:iend)  = itemp3
!   Basis options
    read(iomain,*) ! nspin, iwsusc(0:lmax), ewsusc(1:nbmax,0:lmax,1:nspin)
    read(iomain,*) is, iw(0:nlmax)
!   which atoms are magnetic (is == 2) or not (is == 1)
    issusc(istart:iend) = is
    if (is > nsmax) stop 'susc_input: is > nsmax'
    do il=0,nlmax
!     info on how many basis functions for each l channel
      iwsusc(il,1:nsmax,istart:iend) = iw(il)
!     read the projection energies for each l channel
!     different energies for each spin are possible
      do ib=1,iw(il)
        read(iomain,*) i, ere1, eim1, ere2, eim2
        ewsusc(ib,il,1:nsmax,istart:iend) = cmplx(ere1,eim1)
        if (is == 2) ewsusc(ib,il,2,istart:iend) = cmplx(ere2,eim2)
      end do
    end do
! ******
  end do
! ******
  close(iomain)
  if (any(iasusc(1:nasusc) > natyp)) stop 'susc_input: atom labels'
! -----------------------------------------------------------------
! inconsistencies?
  call check_atomic_options(itc)
! -----------------------------------------------------------------
! allocate global arrays according to input
  call init_arrays
  call init_basis
! ----------------------------------------------------------------------
  if (.not.lhdio) then
!   Set up real spherical harmonics and Gaunt coefficients
    call ymy_gaunts
!   Read SCF mesh from file
    write(*,*) "Reading E-mesh for integration"
    open(file='emesh.scf',unit=iofile,status='old')
    read(iofile,*) nescf
    write(*,*) "nescf=", nescf
    do ie=1,nescf
      read(iofile,*) er, ei, der, dei
      escf(ie)  = cmplx(er,ei)
      ekscf(ie) = sqrt(escf(ie))
      descf(ie) = cmplx(der,dei)
!      descf(ie) = -pi*cmplx(der,dei)
!      write(*,'("e,ek,de=",6es16.8)') escf(ie), ekscf(ie), descf(ie)
    end do
    efscf = real(escf(nescf))
    write(*,'("de sums to=",2f12.6,/)') sum(descf(1:nescf))
    close(iofile)
  end if
! All done!
  end subroutine susc_input
