  subroutine  outsusc_in_soc(itc,my_rank)
! read all info from outsusc.dat
  use global

  implicit none

  integer(kind=i4b), intent(in) :: itc
! -------------------------------------
  real(kind=r8b), parameter :: pi = 4.d0*atan(1.d0)
  integer(kind=i4b) :: istart, iend, ng
  integer(kind=i4b) :: i, is, ib, il, iw(0:4), ntot, itemp1, itemp2, itemp3
  real(kind=r8b)    :: ere1, eim1, ere2, eim2, utemp(0:4), jtemp(0:4), socs, bf, dir(3)
  integer(kind=i4b) :: ib1, ie, ia, nb, nsum
  integer(kind=i4b) :: icount, im
  real(kind=r8b)    :: ram, er, ei, der, dei
  complex(kind=c8b) :: e, ek, de
  complex(kind=c8b), allocatable :: gmat(:,:,:,:), tmat(:,:,:)
  character*60      :: header
  integer(kind=i4b) :: bounds(2)
! -------------------------------------
! No mpi here (my_rank = o from kkrflex) 
  integer(kind=i4b) ::  my_rank 
! -------------------------------------
  write(*,'(/" Opening outsusc.dat"/)')
  open(file='outsusc.dat',unit=iomain,status='old')
  read(iomain,*) ! ---------------------
  read(iomain,*) ! Data for ...
  read(iomain,*) ! ---------------------
  read(iomain,*) ! ikkr, ldos, lsoc, lldau, lbfield, lsusc, lrot, lfit=
  read(iomain,*) ikkr, ldos, lsoc, lldau, lbfield, lsusc, lrot, lfit
  read(iomain,*) ! ibasis, ibasismethod, basistol, lbasisfreeze=
  read(iomain,*) ibasis, ibasismethod, basistol, lbasisfreeze
  read(iomain,*) ! idos, nedos, e0dos, e1dos=
  read(iomain,*) idos, nedos, ere1, eim1, ere2, eim2
  dose0 = cmplx(ere1,eim1); dose1 = cmplx(ere1,eim1)
  read(iomain,*) ! lrhomat, ldaumix, ldauitc=
  read(iomain,*) lrhomat, ldaumix, ldauitc
  read(iomain,*) ! ldynsusc, nomega, omegamin, omegamax, domega, nlmax0=
  read(iomain,*) ldynsusc, nomega, ere1, eim1, ere2, eim2, domega, nlmax0
  omegamin = cmplx(ere1,eim1); omegamax = cmplx(ere2,eim2)
  read(iomain,*) ! lkha, lkxc, lcartesian, lanalytic, lnonanalytic, lenhanced, itermax, lambdamix=
  read(iomain,*) lkha, lkxc, lcartesian, lanalytic, lnonanalytic, lenhanced, itermax, lambdamix
  read(iomain,*) ! ispinrot, urot(1:3), dirmix=
  read(iomain,*) ispinrot, urot(1:3), dirmix3
  read(iomain,*) ! ifit, numd, dend, eshift, lregf, fudge=
  read(iomain,*) ifit, numd, dend, ere1, eim1, lregf, fudge
  eshift = cmplx(ere1,eim1)
  read(iomain,*) ! ---------------------
  read(iomain,*) ! nesusc, nasusc, nbmax, nlmax, isra, nrmax, nsmax=
  read(iomain,*) nesusc, nasusc, nbmax, nlmax, isra, nrmax, nsmax
! -----------------------------------------------------------------
! inconsistencies?
!  call check_global_options(itc)
! -----------------------------------------------------------------
! initialize remaining parameters and allocate arrays needed to read input
!  call init_param
! -----------------------------------------------------------------
! read atom labels in groups
! if a set of atoms is to be treated in the same way, only one read is needed
  read(iomain,*) ! ngroup, igroup(1:ngroup)=
  read(iomain,*) ngroup
  read(iomain,*) igroup(1:ngroup)
  ntot = sum(igroup(1:ngroup))
  if (ntot /= nasusc) stop 'outsusc_in: atom groups and nasusc'
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
    if (lldau .and. itemp2 == 1) then
      ildau(istart:iend) = itemp2
      do il=0,nlmax
        ueff(il,istart:iend) = utemp(il)
        jeff(il,istart:iend) = jtemp(il)
      end do
    else
      ildau(istart:iend) = 0
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
    if (lsusc) then
      isusc(istart:iend) = itemp1
      if (lkha) ikha(istart:iend) = itemp2
      if (lkxc) ikxc(istart:iend) = itemp3
    else
      isusc(istart:iend) = 0
      ikha(istart:iend) = 0
      ikxc(istart:iend) = 0
    end if
!   Basis options
    read(iomain,*) ! nspin, iwsusc(0:lmax), ewsusc(1:nbmax,0:lmax,1:nspin)
    read(iomain,*) is, iw(0:nlmax)
!   which atoms are magnetic (is == 2) or not (is == 1)
    issusc(istart:iend) = is
    if (is > nsmax) stop 'outsusc_in: is > nsmax'
    do il=0,nlmax
!     info on how many basis functions for each l channel
      iwsusc(il,1:nsmax,istart:iend) = iw(il)
    end do
! ******
  end do
! ******
! -----------------------------------------------------------------
! inconsistencies?
!  call check_atomic_options(itc)
! -----------------------------------------------------------------
! *******************************************
! allocate coeff, onsite Gf and structural Gf
! ******************************************* 
  bounds(1) = 1
  bounds(2) = nesusc
  call init_arrays_tgmat(my_rank,1,bounds)  
  call init_gfpartsfit_coeff_gfpq(my_rank,1,bounds)
! allocate global arrays according to input
!  call init_arrays(my_rank)
! -----------------------------------------------------------------
! Read radial mesh, potentials, core densities and basis functions
! reallocate these arrays, in case nbmax changed
!  call init_basis(my_rank)
  nlmsb = 0
  do ia=1,nasusc
    call read_rmesh(ia)
    nsum = 0
    do is=1,nsmax
      do il=0,nlmax
        nb = iwsusc(il,is,ia) 
        nsum = nsum + nb*(2*il+1)
        if (nb > 0) call read_wfns(ia,il,is)
      end do
    end do 
    nlmsb = max(nlmsb,nsum)
  end do
! allocate parts of the GF according to input
!  call init_gfpartsfit(my_rank)
! -----------------------------------------------------------------
! Set up overlaps for the product basis
  call overlaps_gf
! Set up real spherical harmonics and Gaunt coefficients
  call ymy_gaunts
! Set up the angular momentum matrices
  call orbmoment
! Read groundstate density
! ********************
  call read_rho2ns
! ********************
! Set up overlaps for the density basis
  call overlaps_susc2
! Put the SCF potential in the basis --> magdir is needed, moved to pot_correction
!  call build_vscfb
! ----------------------------------------------------------------------
! Now read projection coefficients, t-matrices and structural GFs
! Big loop to read everything from outsusc.dat
! auxiliary storage for reads
! ******************************************************************
! allocate(tmat(lmmax,lmmax,nasusc),gmat(lmmax,lmmax,nasusc,nasusc))
! ******************************************************************
  if (.not.lfit) nescf = nesusc
  write(iodb,'(/"Reading projection coeff and gfpq"/)')
  do ie=1,nesusc
    write(iodb,*) "ie =", ie
!   ****************************************************
!   do is=1,nsmax
!   read(iomain,'(a)') header ! spin separator
!   write(iodb,*) header
!   write(iodb,'(/,"ispin=",i4)') is
!   ****************************************************
!   Read coefficients
    call in_coeffs_soc(2*lmmax,ie,e,ek,de) 
!   ****************************************************
!   call in_coeffs(is,ie,e,ek,de)
!   ****************************************************
!   Save energy mesh; de should be dummy
    if (ikkr == 1 ) then  ! old impurity code
      desusc(ie) = -pi*de
    else                  ! KKRFLEX
      desusc(ie) = -pi*de/nsmax
    end if
    esusc(ie) = e; eksusc(ie) = ek
    if (.not.lfit) then
      escf(ie) = esusc(ie); descf(ie) = desusc(ie)
    end if
    write(iodb,'("e,ek,de=",6es16.8)') e, ek, de
!   ***************************************************
!   Save projection coefficients
!   call save_coeffs(ie,is,.false.)
!   ***************************************************
  end do
  write(iodb,'(/"Reading gmat and tmat"/)') 
  do ie=1,nesusc
    write(iodb,*) "ie =", ie
!   Read t-matrices
    call in_gmat_soc(ie)  
!   Read structural GF
    call in_tmat_soc(ie)
  end do ! ie
! **************************************
!     Put it somewhere
!     This is the collinear case
!      write(iodb,'("Saving t-mat and GF in RAM",2i4)') ie, is
!      call save_tmcoll(tmat,is,ie)
!      call save_gscoll(gmat,is,ie)
!    end do
!  end do
! deallocate(tmat,gmat)
! *************************************
  write(*,'(/" Closing outsusc.dat"/)')
  close(iomain)
! ----------------------------------------------------------------------
! Whether to read SCF energy mesh from file
  if (lscfmesh) then
    write(*,'(" Reading E-mesh for integration")')
    open(file='emesh.scf',unit=iofile,status='old')
     read(iofile,*) nescf
    write(*,'(" nescf=",i8)') nescf
    do ie=1,nescf
      read(iofile,*) er, ei, der, dei
      escf(ie)  = cmplx(er,ei)
      ekscf(ie) = sqrt(escf(ie))
      descf(ie) = cmplx(der,dei)
!      descf(ie) = -pi*cmplx(der,dei)
      write(iodb,'(" e,ek,de=",6es16.8)') escf(ie), ekscf(ie), descf(ie)
    end do
    efscf = real(escf(nescf))
    write(iodb,'(" de sums to=",2f12.6,/)') sum(descf(1:nescf))
    close(iofile)
  end if
! All done
  end subroutine outsusc_in_soc
