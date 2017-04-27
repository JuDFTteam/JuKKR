  subroutine new_input(lmax,natyp,nspin,nrad,itc,my_rank)
! Handles the input file and sets main parameters
  use global

  implicit none

! --> values given in the main program
  integer(kind=i4b), intent(in) :: lmax, natyp, nspin, nrad, itc, my_rank
! -------------------------------------------------------------------------------------------
! line number, column number
  integer(kind=i4b) :: iline, ipos
! was key found?
  logical           :: found
! storage for a single line
  character(len=nchars) :: line
! -------------------------------------------------------------------------------------------
  integer(kind=i4b) :: istart, iend, ig, ng, nbtmp
  integer(kind=i4b) :: i, ie, ia, is, ib, il, iw(0:lmax), ntot, itemp1, itemp2, itemp3
  real(kind=r8b)    :: ere1, eim1, ere2, eim2, utemp(0:lmax), jtemp(0:lmax), socs, bf, dir(3)
  real(kind=r8b)    :: er, ei, der, dei, alat


! Load input file to memory
  call load_input('newinpsusc.dat',nchars,nlines,my_rank)
! Scan input file for global options
  call set_global_options(lmax,natyp,nspin,nrad,itc,my_rank)
! -----------------------------------------------------------------
! Initialize remaining parameters and allocate arrays to read atom info
  call init_param
! -----------------------------------------------------------------
! Combine atoms into groups
  call set_groups(itc,my_rank)
  if (any(iasusc(1:nasusc) > natyp)) stop 'new_input: atom labels'
! Scan input file for basis options
  call set_basis_options(itc,my_rank)
! Scan input file for groundstate options
  call set_groundstate_options(itc,my_rank)
! Scan input file for DOS options
  call set_dos_options(itc,my_rank)
! Scan input file for Jij options
  call set_jij_options(itc,my_rank)
! Scan input file for spin rotation options
  call set_spinrot_options(itc,my_rank)
! Scan input file for magnetic field options
  call set_bfield_options(itc,my_rank)
! Scan input file for self-energy options
  call set_selfe_options(itc,my_rank)
! Scan input file for SOC options
  call set_soc_options(itc,my_rank)
! Scan input file for LDA+U options
  call set_ldau_options(itc,my_rank)
! Scan input file for GF fit options
  call set_fit_options(itc,my_rank)
! Scan input file for susceptibility options
  call set_susc_options(itc,my_rank)
! -----------------------------------------------------------------
! allocate global arrays according to input
  call init_arrays(my_rank)
  call init_basis(my_rank)
! ----------------------------------------------------------------------
  if (.not.lhdio) then
!   Set up real spherical harmonics and Gaunt coefficients
    call ymy_gaunts
!   Whether to read the energy mesh from a separate file
    if (lscfmesh) then
      if (my_rank == 0) then
        write(*,'(" Reading E-mesh for integration")')
      end if ! my_rank
      open(file='emesh.scf',unit=iofile,status='old')
      read(iofile,*) nescf
      if (my_rank ==0) then
        write(*,'(" nescf=",i8)') nescf
      end if ! my_rank
      do ie=1,nescf
        read(iofile,*) er, ei, der, dei
        escf(ie)  = cmplx(er,ei)
        ekscf(ie) = sqrt(escf(ie))
        descf(ie) = cmplx(der,dei)
!        descf(ie) = -pi*cmplx(der,dei)
!        write(*,'("e,ek,de=",6es16.8)') escf(ie), ekscf(ie), descf(ie)
      end do
      efscf = real(escf(nescf))
      if (my_rank ==0) then  
        write(*,'(" de sums to=",2f12.6,/)') sum(descf(1:nescf))
      end if ! my_rank
      close(iofile)
    else
      if (my_rank ==0) then
        write(*,'(" E-mesh for integration the same as in outsusc")')
      end if ! my_rank
    end if
!   Whether to read atomic positions from a separate file
    if (lpositions) then
      if (my_rank ==0) then
        write(*,'(" Reading atomic positions")')
        write(iodb,'(" Atomic positions:")')
      end if ! my_rank 
!     Read atomic positions
      open(file='positions.dat',unit=iofile,status='old')
      read(iofile,*) alat ! nasusc
      do ia=1,nasusc
        read(iofile,*) ri(:,ia)
        ri(:,ia) = ri(:,ia)*alat
        if (my_rank ==0) then
          write(iodb,'(i4,3f12.8)') ia, ri(:,ia)
        end if ! my_rank
      end do
      close(iofile)
    else
      if (my_rank ==0) then
        write(*,'(" WARNING: Write interface to get positions from host KKR code")')
      end if ! my_rank
    end if
  end if
! All done!
  end subroutine new_input
