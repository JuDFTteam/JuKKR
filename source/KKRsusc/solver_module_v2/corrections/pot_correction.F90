  subroutine pot_correction(itc)
! Driver routine for including potential corrections
  use global

  implicit none

! Current SCF iteration
  integer(kind=i4b), intent(in) :: itc
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: iu = (0.d0,1.d0), cone = (1.d0,0.d0)
! Speed of light
  real(kind=r8b), parameter :: c = 274.0720442d0
! temporary parameters
  logical, parameter :: onsite = .true., struct = .true.
  logical, parameter :: intra = .true., inter = .true.
! Spherical spin-averaged radial SOC potential
  complex(kind=c8b) :: vsoc(nrmax,nasusc,nesusc)
! Potential correction
  complex(kind=c8b), allocatable :: deltapot(:,:,:)
! Spin moment from old susc sum rule
  complex(kind=c8b), allocatable :: msgfz(:,:,:,:,:), msxcz(:,:,:,:,:), mssocz(:,:,:,:,:), mssocpm(:,:,:,:,:)
! Different pieces of the spin splitting for old sum rule
  complex(kind=c8b), allocatable :: vxcdiff(:,:,:,:,:), vsocz(:,:,:,:,:), vsocpm(:,:,:,:,:,:)
! Spin moment from susc new sum rule
  complex(kind=c8b), allocatable :: msgftwist(:,:,:,:,:,:), msxctwist(:,:,:,:,:,:,:), mssoctwist(:,:,:,:,:,:,:)
! Spin-dependent potentials for new sum rule
  complex(kind=c8b), allocatable :: bxctwist(:,:,:,:), bsoctwist(:,:,:,:)
! Change in t-matrix
  complex(kind=c8b), allocatable :: dtmat(:,:,:)
! Split for the energy loop
  integer(kind=i4b) :: nesplit
! Misc
  integer(kind=i4b) :: ie, ia, i, j, ia2, is, ilms
  complex(kind=c8b) :: spinrot(nsmax,nsmax)
  real(kind=r8b)    :: start, finish, magdir0(3,nasusc), magdir1(3,nasusc)
! Variables for handling model self-energy
  integer(kind=i4b) :: naselfe, norbmax, npmax, norb
  integer(kind=i4b), allocatable :: iaselfe(:), iselfe(:), iorb(:), senumd(:), sedend(:)
  real(kind=r8b),    allocatable :: efshift(:), se_infty(:,:,:), se_fit(:,:,:,:)
  integer(kind=i4b) :: i3(3), i1, ib1, ilm1, is1, i2, ib2, ilm2, is2, ja, j1, jb1, jlm1, js1, j2, jb2, jlm2, js2, k, iatom, jlms, js, ib, jb, ilm, jlm
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Check if there is anything to do here:
  write(*,'(/,"Entering pot_correction",/)')
 
  if (.not.lsoc .and. .not.lldau .and. .not.lbfield .and. .not.lrot .and. .not.lsusc .and. .not.lsemodel) then
! NEW: Lichtenstein Jij
    if (ljij) then
      write(*,'(/,"Computing Jij only",/)')
!      call symmetrize_tgf
!     the actual calculation of Jij's
      if (ljijtensor) then
        call anisotropic_jij(.false.)
      else
        call lichtenstein_jij
      end if
    else
      write(*,'(/,"Nothing to do in pot_correction",/)')
    end if
!   that's it
    return
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Input is collinear ASA -> enforce this symmetry
!  call symmetrize_tgf
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (lsemodel) then
    write(*,'(/,"**** Handling model self-energy ****",/)')
    open(file='se_model.dat',unit=iofile,status='old')
!   number of atoms to apply model self-energy to, maximum number of orbitals
    read(iofile,*) ! skip line
    read(iofile,*) naselfe, norbmax, npmax
    write(*,'("naselfe=",i4," norbmax=",i4," npmax=",i4)') naselfe, norbmax, npmax
    if (naselfe > nasusc) stop 'SE model: naselfe > nasusc'
    if (norbmax > lmmax)  stop 'SE model: norbmax > lmmax'
    allocate(iaselfe(nasusc),iselfe(naselfe))
    allocate(iorb(norbmax),senumd(nasusc),sedend(nasusc))
    allocate(efshift(nasusc),se_infty(2,nlms,nasusc),se_fit(2,nlms,npmax,nasusc))
!   list of atoms numbered as in nasusc
    read(iofile,*) ! skip line
    read(iofile,*) iselfe(1:naselfe)
    write(*,'("iselfe=",100i4)') iselfe(1:naselfe)
!   initialize key controlling whether SE is applied to an atom
    iaselfe(1:nasusc) = 0
    do ia2=1,naselfe
      iaselfe(iselfe(ia2)) = 1
    end do
    write(*,'("iaselfe=",100i4)') iaselfe(1:nasusc)
!   now read model data for each atom
    do ia2=1,naselfe
      ia = iselfe(ia2)
      read(iofile,*) ! skip line
      read(iofile,*) i, norb, senumd(ia), sedend(ia), efshift(ia)
      write(*,'("ia=",i4," norb=",i4," numd, dend=",2i4," efshift=",f16.8)') i, norb, senumd(ia), sedend(ia), efshift(ia)
      if (norb > norbmax) stop 'SE model: norb(ia) > norbmax'
      if (1+senumd(ia)+sedend(ia) > npmax) stop 'SE model: 1+senumd(ia)+sedend(ia) > npmax'
      read(iofile,*) ! skip line
      read(iofile,*) iorb(1:norb)
      write(*,'("iorb=",100i4)') iorb(1:norb)
!     set to which orbitals the self-energy applies
!     read parameters for each orbital
      read(iofile,*) ! skip line
      se_infty(:,:,ia) = 0.d0; se_fit(:,:,:,ia) = 0.d0
      do i=1,norb
        do is=1,nsmax
          ilms = lms2i(iorb(i),is)
          read(iofile,*) se_infty(:,ilms,ia), se_fit(:,ilms,1:1+senumd(ia)+sedend(ia),ia)
          write(*,'("SE params=",1000es8.1)') se_infty(:,ilms,ia), se_fit(:,ilms,1:1+senumd(ia)+sedend(ia),ia)
        end do
      end do
    end do
    close(iofile)
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  allocate(deltapot(nlmsb,nlmsb,nasusc),dtmat(nlms,nlms,nasusc))
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! A global spin rotation is done by simple spin matrices
  call cpu_time(start)
! Set spin axes
  call spin_directions(magdir0,magdir1)
!  if (lrot) call spinrot_gf(magdir0,magdir1)
! with the input GF being ASA this can be made faster
  if (lrot) then
    if (lgrefsph) then
      call spinrot_tgf_sph(magdir0,magdir1)
!   this should handle full-potential like input
    else
      call spinrot_tgf(magdir0,magdir1)
    end if
  end if
!  write(*,*) "after spinrot_gf"
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Prepare for the magnetization sum rule
  if (lsusc .and. lsumrule) then
!   noncollinear sum rule
    if (lnewsumrule) then
      nesplit = 0
      allocate(bxctwist(nlmsb,nlmsb,3,nasusc),msxctwist(3,nbmax,nbmax,lmmax,lmmax,nasusc,nasusc2))
      allocate(bsoctwist(nlmsb,nlmsb,3,nasusc),mssoctwist(3,nbmax,nbmax,lmmax,lmmax,nasusc,nasusc2))
      allocate(msgftwist(3,nbmax,nbmax,lmmax,lmmax,nasusc2))
      bxctwist = 0.d0; msxctwist = 0.d0; bsoctwist = 0.d0; mssoctwist = 0.d0; msgftwist = 0.d0
!     Potential in the projection basis
      call build_vscfb2(bxctwist)
!   old sum rule
    else
      nesplit = nescf
      allocate(msgfz(nbmax,nbmax,lmmax,lmmax,nasusc2),msxcz(nbmax,nbmax,lmmax,lmmax,nasusc2))
      allocate(vxcdiff(nbmax,nbmax,lmmax,lmmax,nasusc))
      msgfz = 0.d0; msxcz = 0.d0
      allocate(mssocz(nbmax,nbmax,lmmax,lmmax,nasusc2),mssocpm(nbmax,nbmax,lmmax,lmmax,nasusc2))
      allocate(vsocz(nbmax,nbmax,lmmax,lmmax,nasusc),vsocpm(nbmax,nbmax,lmmax,lmmax,nasusc,nasusc))
      mssocz = 0.d0; mssocpm = 0.d0
!     Potential in the projection basis
      call build_vscfb
!     Splitting coming from xc potential
      call build_vxcdiff2(vlmsbgf,magdir1,vxcdiff)
    end if
  else
!   Potential in the projection basis
    call build_vscfb
    nesplit = 0
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! First iteration: save density matrix to file
  if (lnewdensity .and. itc == 1 .and. .not.lrhomat) call save_rhomat
! Read density matrix from file
  if (lldau) call read_rhomat
! Rotate density matrix
  if (lrot) call spinrot_rhomat(magdir0,magdir1)
!  write(*,*) "after spinrot_rhomat"
! Radial part of SOC potential
  if (lsoc) then
    do ie=1,nesusc
      do ia=1,nasusc
        if (isoc(ia) > 0) call build_vsoc(c,esusc(ie),zat(ia),nrpts(ia),rmesh(:,ia),vr(:,ia),socscaling(ia),vsoc(:,ia,ie),npanat(ia),ircutat(:,ia)) 
      end do
    end do
  end if
! calculate gradient of scalar relativistic mass (test calc for susc)
! added by Sascha
! allocate grad_mass
  if(lscalarcorr) then
    if(allocated(grad_mass)) deallocate(grad_mass)
    allocate(grad_mass(1:3,1:nrmax,1:lmmax,1:nasusc))
    grad_mass(:,:,:,:)=0.d0
    do ie=1,nesusc
      do ia=1,nasusc
        call build_grad_mass(c,esusc(ie),zat(ia),nrpts(ia),rmesh(:,ia),vr(:,ia),socscaling(ia),ia)
      end do
    end do
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop over energies for susc sum rule
! ***************
  do ie=1,nesplit
! ***************
    deltapot = 0.d0; dtmat = 0.d0
!   --------------------------------------------------------------------
    if (lsemodel) then
!     Construct the model self-energy
      do ia=1,nasusc
!       SE correction
        if (iaselfe(ia) == 1) call build_selfe(ie,ia,senumd(ia),sedend(ia),efshift(ia),se_infty(:,:,ia),se_fit(:,:,:,ia),deltapot(:,:,ia))
      end do
!      write(*,*) "after build_pot", ie
!     --------------------------------------------------------------------
!     Now solve the Lippmann-Schwinger equation in the basis
!     This updates the regular solutions and the single-site GF
      call lippmann_schwinger(ie,nasusc,nlms,nlmsb,nlmsba,deltapot,pzl(:,:,:,ie),pzr(:,:,:,ie),dtmat,gfpq(:,:,:,ie))
!      write(*,*) "after lippmann_schwinger", ie
!     Next the structural GF has to be updated
      call structural_gf(ie,nasusc,nlms,nalms,alms2i,dtmat,gstruct(:,:,ie))
!      write(*,*) "after structural_gf", ie
!     Update total t-matrix
      call zaxpy(nlms*nlms*nasusc,cone,dtmat,1,tmatrix(:,:,:,ie),1)
    end if
!   --------------------------------------------------------------------
    if (lldau) then
!     First compute the LDA+U correction
      do ia=1,nasusc
!       LDA+U correction
        if (ildau(ia) == 1) call build_vldaub(ie,ia,deltapot(:,:,ia))
      end do
!      write(*,*) "after build_pot", ie
!     --------------------------------------------------------------------
!     Now solve the Lippmann-Schwinger equation in the basis
!     This updates the regular solutions and the single-site GF
      call lippmann_schwinger(ie,nasusc,nlms,nlmsb,nlmsba,deltapot,pzl(:,:,:,ie),pzr(:,:,:,ie),dtmat,gfpq(:,:,:,ie))
!      write(*,*) "after lippmann_schwinger", ie
!     Next the structural GF has to be updated
      call structural_gf(ie,nasusc,nlms,nalms,alms2i,dtmat,gstruct(:,:,ie))
!      write(*,*) "after structural_gf", ie
!     Update total t-matrix
      call zaxpy(nlms*nlms*nasusc,cone,dtmat,1,tmatrix(:,:,:,ie),1)
    end if
!   --------------------------------------------------------------------
!   First compute the potential corrections that don't flip spin
    vsocz = 0.d0; deltapot = 0.d0; dtmat = 0.d0
    if (lsoc .or. lbfield ) then
!     SOC Lz.Sz if not applied in the host
      if (.not.lsoc_new) then
        if (lsoc) then
          do ia=1,nasusc
            if (isoc(ia) == 1 .or. isoc(ia) == 2) then
              call build_vsocb(ie,ia,vsoc(:,ia,ie),deltapot(:,:,ia),magdir1(:,ia),2)
!             call torque_soc2(ie,ia,vsoc(:,ia,ie),magdir0(:,ia),magdir1(:,ia))
            end if
          end do
        end if
      end if ! SOC host
!     Zeeman Bz.(Lz+Sz)
      if (lbfield) then
        do ia=1,nasusc
          if (ibfield(ia) > 0) call build_zeeman(ia,magdir1(:,ia),deltapot(:,:,ia),2)
        end do
      end if
!     Splitting coming from SOC z and Bz
      call build_vsocz2(deltapot,magdir1,vsocz)
!      write(*,*) "after build_pot", ie
!     --------------------------------------------------------------------
!     Now solve the Lippmann-Schwinger equation in the basis
!     This updates the regular solutions and the single-site GF
      call lippmann_schwinger(ie,nasusc,nlms,nlmsb,nlmsba,deltapot,pzl(:,:,:,ie),pzr(:,:,:,ie),dtmat,gfpq(:,:,:,ie))
!      write(*,*) "after lippmann_schwinger", ie
!     Next the structural GF has to be updated
      call structural_gf(ie,nasusc,nlms,nalms,alms2i,dtmat,gstruct(:,:,ie))
!      write(*,*) "after structural_gf", ie
!     Update total t-matrix
      call zaxpy(nlms*nlms*nasusc,cone,dtmat,1,tmatrix(:,:,:,ie),1)
    end if
!   --------------------------------------------------------------------
!   Now compute the potential corrections that flip spin
    vsocpm = 0.d0; deltapot = 0.d0; dtmat = 0.d0
    if (lsoc .or. lbfield) then
!     SOC +-
      if (.not.lsoc_new) then
        if (lsoc) then
          do ia=1,nasusc
            if (isoc(ia) == 1 .or. isoc(ia) == 3) then
              call build_vsocb(ie,ia,vsoc(:,ia,ie),deltapot(:,:,ia),magdir1(:,ia),3)
            end if
          end do
        end if
      end if ! SOC host
!     Zeeman +-
      if (lbfield) then
        do ia=1,nasusc
          if (ibfield(ia) > 0) call build_zeeman(ia,magdir1(:,ia),deltapot(:,:,ia),3)
        end do
      end if
!     Splitting coming from SOC +- and external B+-
      call build_vsocpm2(ie,deltapot,magdir1,vsocpm)
!     --------------------------------------------------------------------
!     Now solve the Lippmann-Schwinger equation in the basis
!     This updates the regular solutions and the single-site GF
      call lippmann_schwinger(ie,nasusc,nlms,nlmsb,nlmsba,deltapot,pzl(:,:,:,ie),pzr(:,:,:,ie),dtmat,gfpq(:,:,:,ie))
!      write(*,*) "after lippmann_schwinger", ie
!     Next the structural GF has to be updated
      call structural_gf(ie,nasusc,nlms,nalms,alms2i,dtmat,gstruct(:,:,ie))
!      write(*,*) "after structural_gf", ie
!     Update total t-matrix
      call zaxpy(nlms*nlms*nasusc,cone,dtmat,1,tmatrix(:,:,:,ie),1)
    end if
!   --------------------------------------------------------------------
!   Full SOC added to the GF: check sum rule
    call ms_from_bxc3(ie,vxcdiff,vsocz,vsocpm,magdir1,msgfz,msxcz,mssocz,mssocpm)
!   --------------------------------------------------------------------
! ******
  end do
! ******
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop over remaining energies
  if (lsoc .or. lldau .or. lbfield .or. lsemodel .or. lsusc) then
! ********************
  do ie=nesplit+1,nesusc
! ********************
    deltapot = 0.d0; dtmat = 0.d0
!   First compute the potential corrections
    do ia=1,nasusc
!   --------------------------------------------------------------------
!     SE correction
      if (lsemodel) then 
        if (iaselfe(ia) == 1) call build_selfe(ie,ia,senumd(ia),sedend(ia),efshift(ia),se_infty(:,:,ia),se_fit(:,:,:,ia),deltapot(:,:,ia))
      end if
!     LDA+U correction
      if (lldau) then
        if (ildau(ia) == 1) call build_vldaub(ie,ia,deltapot(:,:,ia))
      end if
!     SOC if not applied in the host
      if (.not.lsoc_new) then
        if (lsoc) then
          if (isoc(ia) > 0) call build_vsocb(ie,ia,vsoc(:,ia,ie),deltapot(:,:,ia),magdir1(:,ia),isoc(ia))
!         call torque_soc2(ie,ia,vsoc(:,ia,ie),magdir0(:,ia),magdir1(:,ia))
        end if
      end if ! SOC host
!     Zeeman
      if (lbfield) then
        if (ibfield(ia) > 0) call build_zeeman(ia,magdir1(:,ia),deltapot(:,:,ia),1)
      end if
    end do
!    write(*,*) "after build_pot", ie
!   Now solve the Lippmann-Schwinger equation in the basis
!   This updates the regular solutions and the single-site GF
    call lippmann_schwinger(ie,nasusc,nlms,nlmsb,nlmsba,deltapot,pzl(:,:,:,ie),pzr(:,:,:,ie),dtmat,gfpq(:,:,:,ie))
!    write(*,*) "after lippmann_schwinger", ie
!   Next the structural GF has to be updated
    call structural_gf(ie,nasusc,nlms,nalms,alms2i,dtmat,gstruct(:,:,ie))
!    write(*,*) "after structural_gf", ie
!   Update total t-matrix
    call zaxpy(nlms*nlms*nasusc,cone,dtmat,1,tmatrix(:,:,:,ie),1)

!   Include the SOC host in sumrule    
    if (lsusc .and. lsumrule .and. lnewsumrule) then
      if (lsoc_new .and. lsoc) then 
        do ia = 1, nasusc
          ! compute vsoc in the basis (careful might be inconsistent)
          call build_vsocb(ie,ia,vsoc(:,ia,ie),deltapot(:,:,ia),magdir1(:,ia),isoc(ia))
          ! recompute the bfield contribution since it was overwritten by build_vsoc
          if (lbfield) then
            if (ibfield(ia) > 0) call build_zeeman(ia,magdir1(:,ia),deltapot(:,:,ia),1)
          end if                                                                 
        end do ! nasusc
      end if ! SOC host
      call build_bsoctwist(deltapot,bsoctwist)

!     noncollinear sum rule
      call ms_from_bxc4(ie,bxctwist,bsoctwist,magdir1,msgftwist,msxctwist,mssoctwist)
    end if
! ******
  end do
! ******
  end if
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  call cpu_time(finish)
  write(*,'("***********************************")')
  write(*,'(" Pot correction time:",f12.4," s")') finish - start
  write(*,'("***********************************"/)')
! If the fit is being used, it also has to be recomputed
! Recalculate the groundstate properties
!  call ms_from_bxc
  call groundstate_new
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Compute the static susceptibility and the kernels
!  if (lsusc) call static_susc_more(onsite,struct)
! Added by Sascha:
  if (lsusc .and. lcurrcorr) call gradient_susc_basis
  if (lsusc) call static_susc2(onsite,struct)
  if (lsusc .and. lkxc .and. lsumrule) call kxc_sumrule
  if (lsusc .and. lkxc) call build_kxcalda2
  if (lsusc .and. lkha) call build_khartree(intra,inter)
  if (lsusc .and. (lkxc .or. lkha) .and. lenhanced) then
    call susc_denominator
    call static_susc_enhanced
  end if
  if (lsusc .and. ldynsusc) call dyn_susc_expansion
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! If there was a spin rotation the spin direction has to be updated
! NEW: Lichtenstein Jij
  if (ljij) then
    if (ljijtensor) then
      call anisotropic_jij(lldau)
    else
      call lichtenstein_jij
    end if
  end if
! This might have to move before spinrot_rhomat
  if (lrot) call new_directions(magdir0,magdir1)
! rotate back the density matrix
  if (lrot) then
    do ia=1,nasusc
!      magdir1(:,ia) = (/0.d0,0.d0,1.d0/)
      magdir1(:,ia) = magdir0(:,ia)
      magdir0(:,ia) = magdir(:,ia)
    end do
    call spinrot_rhomat(magdir0,magdir1)
  end if
! Mix and save density matrix
  if (lldau .and. itc > ldauitc) then
!   Mix density matrix; compute vshift
    call mix_rhomat(ldaumix)
!   Save density matrix
    if (lnewdensity) call save_rhomat
  end if
! Free memory
  deallocate(deltapot,dtmat)
  if (lsusc .and. lsumrule) then
    if (lnewsumrule) then
      deallocate(msgftwist,msxctwist,bxctwist,mssoctwist,bsoctwist)
    else
      deallocate(msgfz,msxcz,vxcdiff,mssocz,mssocpm,vsocz,vsocpm)
    end if
  end if
  if (lsemodel) deallocate(iaselfe,iselfe,iorb,senumd,sedend,efshift,se_infty,se_fit)
  write(*,'(/,"Leaving pot_correction",/)')
! All done!
  end subroutine pot_correction
