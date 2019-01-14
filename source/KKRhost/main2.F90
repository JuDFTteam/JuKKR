!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!-----------------------------------------------------------------------------------
!> Summary: Wrapper module for the calculation of the DFT quantities for the JM-KKR package
!> Author: Philipp Ruessmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,       
!> and many others ...
!> The code uses the information obtained in the main0 module, this is
!> mostly done via the `get_params_2()` call, that obtains parameters of the type
!> `t_params` and passes them to local variables
!-----------------------------------------------------------------------------------
module mod_main2


  use :: mod_datatypes, only: dp
  private
  public :: main2

contains

  !-------------------------------------------------------------------------------
  !> Summary: Main wrapper routine dealing with the calculation of the DFT quantities
  !> Author: Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
  !> and many others ...
  !> Category: potential, xc-potential, total-energy, KKRhost
  !> Deprecated: False 
  !> Calculates the potential from density, exc-potential, calculate total energy, ...
  !-------------------------------------------------------------------------------
  !> @note JC: there seems to be an array called sum, this is dangerous as it 
  !> can get confuder by the FORTRAN intrinsic function sum(). The name should
  !> be changed.
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine main2()

    use :: mod_constants, only: pi
    use :: mod_runoptions, only: disable_charge_neutrality, no_madelung, print_program_flow, relax_SpinAngle_Dirac, &
      search_Efermi, simulate_asa, slow_mixing_Efermi, symmetrize_potential_cubic, symmetrize_potential_madelung, &
      use_decimation, use_rigid_Efermi, use_semicore, use_spherical_potential_only, write_deci_tmat, write_kkrimp_input, &
      write_madelung_file, write_potential_tests, write_rho2ns
    use :: global_variables, only: krel, ipand, npotd, natomimpd, lmxspd, iemxd, nspotd, irid, ngshd, linterface, &
      nfund, ncelld, irmd, nembd1, nembd, irmind, lmmaxd, wlength, natypd, naezd, lmpotd, lpotd, lmaxd, nspind, nspotd, &
      ipand, ngshd, irid, nfund, ncelld
    use :: mod_main0, only: lcore, ncore, ircut, ipan, ntcell, lpot, nlbasis, nrbasis, nright, nleft, natomimp, atomimp, &
      natyp, naez, lly, lmpot, nsra, ins, nspin, lmax, imix, qbound, fcm, itdbry, irns, kpre, kshape, kte, kvmad, kxc, &
      icc, ishift, ixipol, kforce, ifunm, lmsp, imt, irc, irmin, irws, llmsp, ititle, nfu, hostimp, ilm_map, imaxsh, &
      ielast, npol, npnt1, npnt2, npnt3, itscf, scfsteps, iesemicore, kaoez, iqat, noq, npolsemi, n1semi, n2semi, n3semi, &
      zrel, jwsrel, irshift, mixing, lambda_xc, a, b, thetas, drdi, rmesh, zat, rmt, rmtnew, rws, emin, emax, tk, alat, &
      cmomhost, conc, gsh, ebotsemi, emusemi, tksemi, vins, visp, rmrel, drdirel, vbc, r2drdirel, ecore, ez, wez, txc, lly, &
      lrhosym, idoldau, lopt, nshell, nemb, fsemicore, qmgam, fact, qmphi, qmtet, ipf, idosemicore, thesme, dez, vtrel, btrel
    use :: mod_types, only: t_inc
    use :: mod_wunfiles, only: t_params, get_params_2, read_density, save_emesh, save_scfinfo
    use :: mod_profiling, only: memocc
    use :: mod_brydbm, only: brydbm
    use :: mod_ecoub, only: ecoub
    use :: mod_epathtb, only: epathtb
    use :: mod_epotinb, only: epotinb
    use :: mod_espcb, only: espcb
    use :: mod_etotb1, only: etotb1
    use :: mod_force, only: force
    use :: mod_forceh, only: forceh
    use :: mod_forcxc, only: forcxc
    use :: mod_mtzero, only: mtzero
    use :: mod_mdirnewang, only: mdirnewang
    use :: mod_mixstr, only: mixstr
    use :: mod_rhosymm, only: rhosymm
    use :: mod_relpotcvt, only: relpotcvt
    use :: mod_rhototb, only: rhototb
    use :: mod_vmadelblk, only: vmadelblk
    use :: mod_vintras, only: vintras
    use :: mod_vinterface, only: vinterface
    use :: mod_scfiterang, only: scfiterang
    use :: mod_rites, only: rites
    use :: mod_writekkrflex, only: writekkrflex
    use :: mod_vxcdrv, only: vxcdrv
    use :: mod_rinit, only: rinit 
    use :: mod_convol, only: convol
    use :: mod_types, only: t_madel

    implicit none

    real (kind=dp), parameter :: eps = 1.0d-12

    integer :: iobroy
    parameter (iobroy=20)
    integer :: nmvecmax
    parameter (nmvecmax=4)
    ! ..
    ! .. Local scalars
    integer :: nk                  !! ITERMDIR variables
    integer :: irc1
    integer :: ipot
    integer :: nmvec               !! ITERMDIR variables
    integer :: icont
    integer :: ispin
    integer :: irmin1
    integer :: lmaxd1
    integer :: lsmear
    integer :: lrecabmad
    integer :: i_stat, i_all
    integer :: i, j, ie, i1, i2, ih, it, io, lm, ir, irec
    integer :: special_straight_mixing !!id to specify modified straight mixing scheme: 0=normal, 1=alternating mixing factor (i.e. reduced mixing factor in every odd iteration), 2=charge-neurality based mixing factor (former: 'alt mix' and 'spec mix')
    real (kind=dp) :: df
    real (kind=dp) :: rv
    real (kind=dp) :: mix
    real (kind=dp) :: sum
    real (kind=dp) :: fpi        !! 4\(\pi\)
    real (kind=dp) :: rfpi       !! \(\sqrt{4\pi}\)
    real (kind=dp) :: efold
    real (kind=dp) :: efnew
    real (kind=dp) :: denef
    real (kind=dp) :: fsold
    real (kind=dp) :: vshift       ! fxf
    real (kind=dp) :: rmsavm
    real (kind=dp) :: rmsavq
    real (kind=dp) :: rmsav0
    real (kind=dp) :: chrgnt
    real (kind=dp) :: chrgold
    real (kind=dp) :: excdiff      !! Scale magn. part of xc-potential
    real (kind=dp) :: e2shift
    real (kind=dp) :: erravang     !! ITERMDIR variables
    real (kind=dp) :: chrgsemicore
    ! .. Local Arrays
    integer, dimension (natypd) :: lcoremax
    integer, dimension (natypd, naezd) :: itoq
    integer, dimension (20, natypd) :: nkcore
    integer, dimension (20, npotd) :: kapcore
    real (kind=dp), dimension (natypd) :: eu !! LDA+U
    real (kind=dp), dimension (natypd) :: edc !! LDA+U
    real (kind=dp), dimension (lmpotd) :: c00
    real (kind=dp), dimension (lmpotd) :: bvmad
    real (kind=dp), dimension (natypd) :: denefat
    real (kind=dp), dimension (2) :: vmt_init
    real (kind=dp), dimension (irmd, npotd) :: rhoc !! core charge density
    real (kind=dp), dimension (lmpotd, lmpotd) :: avmad
    real (kind=dp), dimension (0:lpotd, natypd) :: excnm !! Scale magn. part of xc-potential
    real (kind=dp), dimension (lmpotd, naezd) :: vinters
    real (kind=dp), dimension (irmd, npotd) :: vspsmdum
    logical, dimension (natypd, lmpot) :: lpotsymm
    ! -------------------------------------------------------------------------
    ! ITERMDIR variables
    ! -------------------------------------------------------------------------
    real (kind=dp), dimension (natypd, nmvecmax) :: mvgam
    real (kind=dp), dimension (natypd, nmvecmax) :: mvphi
    real (kind=dp), dimension (natypd, nmvecmax) :: mvtet
    complex (kind=dp), dimension (natypd, 3, nmvecmax) :: mvevi
    complex (kind=dp), dimension (natypd, 3, nmvecmax) :: mvevief
    ! -------------------------------------------------------------------------
    ! ECOU(0:LPOT,NATYP)      ! Coulomb energy
    ! EPOTIN(NATYP),          ! energy of input potential (EPOTINB
    ! ESPC(0:3,NPOTD),        ! energy single particle core
    ! ESPV(0:LMAXD1,NPOTD)    ! energy single particle valence
    ! ! both changed for the relativistic case
    ! EXC(0:LPOT,NATYP),      ! E_xc
    ! -------------------------------------------------------------------------
    real (kind=dp), dimension (natypd) :: epotin !! energy of input potential (EPOTINB
    real (kind=dp), dimension (0:3, npotd) :: espc !! energy single particle core
    real (kind=dp), dimension (0:lpotd, natypd) :: exc !! exchange correlation energy
    real (kind=dp), dimension (0:lpotd, natypd) :: ecou !! Coulomb energy
    real (kind=dp), dimension (0:lmaxd+1, npotd) :: espv !! energy single particle valence both changed for the relativistic case
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd) :: rhoorb
    real (kind=dp), dimension (krel*20+(1-krel), npotd) :: ecorerel
    ! -------------------------------------------------------------------------
    ! CMINST(LMPOT,NATYP)            ! charge moment of interstitial
    ! CMOM(LMPOT,NATYP)              ! LM moment of total charge
    ! CHRGATOM(NATYP,
    ! 2*KREL+(1-KREL)*NSPIND) ! total charge per atom
    ! -------------------------------------------------------------------------
    real (kind=dp), dimension (lmpotd, natypd) :: cmom !! LM moment of total charge
    real (kind=dp), dimension (lmpotd, natypd) :: cminst !! charge moment of interstitial
    real (kind=dp), dimension (natypd, 2*krel+(1-krel)*nspind) :: chrgatom !! total charge per atom
    ! -------------------------------------------------------------------------
    ! FORCES
    ! -------------------------------------------------------------------------
    real (kind=dp), dimension (-1:1, natypd) :: flm !! Forces
    real (kind=dp), dimension (-1:1, natypd) :: flmc !! Forces
    ! -------------------------------------------------------------------------
    ! For SIMULASA
    ! -------------------------------------------------------------------------
    integer :: ipos, ilm_mapp, ias

    ! .. Allocatable arrays
    real (kind=dp), dimension (:, :, :), allocatable :: vons !! output potential (nonspherical VONS)

    ! -------------------------------------------------------------------------
    ! R2NEF (IRMD,LMPOT,NATYP,2)  ! rho at FERMI energy
    ! RHO2NS(IRMD,LMPOT,NATYP,2)  ! radial density
    ! nspin=1            : (*,*,*,1) radial charge density
    ! nspin=2 or krel=1  : (*,*,*,1) rho(2) + rho(1) -> charge
    ! (*,*,*,2) rho(2) - rho(1) -> mag. moment
    ! RHOC(IRMD,NPOTD)              ! core charge density
    ! -------------------------------------------------------------------------
    real (kind=dp), dimension (:, :, :, :), allocatable :: r2nef !! rho at FERMI energy
    real (kind=dp), dimension (:, :, :, :), allocatable :: rho2ns !! radial density
    ! -------------------------------------------------------------------------
    ! Scale magn. part of xc-potential:
    real (kind=dp), dimension (:, :, :), allocatable :: vxcm
    real (kind=dp), dimension (:, :, :), allocatable :: vxcnm
    real (kind=dp), dimension (:, :, :, :), allocatable :: rho2nsnm

    ! .. External Subroutines


    lmaxd1 = lmax + 1

    ! Allocations
    ! allocate(THETAS(IRID,NFUND,NCELLD),stat=i_stat)
    ! call memocc(i_stat,product(shape(THETAS))*kind(THETAS),'THETAS','main2')
    ! allocate(THESME(IRID,NFUND,NCELLD),stat=i_stat)
    ! call memocc(i_stat,product(shape(THESME))*kind(THESME),'THESME','main2')
    allocate (vons(irmd,lmpot,npotd), stat=i_stat)
    call memocc(i_stat, product(shape(vons))*kind(vons), 'VONS', 'main2')
    ! allocate(VINS(IRMIND:IRMD,LMPOT,NSPOTD),stat=i_stat)
    ! call memocc(i_stat,product(shape(VINS))*kind(VINS),'VINS','main2')
    allocate (vxcm(irmd,lmpot,npotd), stat=i_stat)
    call memocc(i_stat, product(shape(vxcm))*kind(vxcm), 'VXCM', 'main2')
    allocate (vxcnm(irmd,lmpot,npotd), stat=i_stat)
    call memocc(i_stat, product(shape(vxcnm))*kind(vxcnm), 'VXCNM', 'main2')
    allocate (r2nef(irmd,lmpot,natyp,2), stat=i_stat)
    call memocc(i_stat, product(shape(r2nef))*kind(r2nef), 'R2NEF', 'main2')
    allocate (rho2ns(irmd,lmpot,natyp,2), stat=i_stat)
    call memocc(i_stat, product(shape(rho2ns))*kind(rho2ns), 'RHO2NS', 'main2')
    allocate (rho2nsnm(irmd,lmpot,natyp,2), stat=i_stat)
    call memocc(i_stat, product(shape(rho2nsnm))*kind(rho2nsnm), 'RHO2NSNM', 'main2')

    ! Consistency check
    if ((krel<0) .or. (krel>1)) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
    if ((krel==1) .and. (nspind==2)) stop ' set NSPIN = 1 for KREL = 1 in the inputcard'
    ! -------------------------------------------------------------------------
    ! This routine previously used to read from unformatted files created by
    ! the main0 module, now  instead of unformatted files take parameters from
    ! types defined in wunfiles.F90
    ! -------------------------------------------------------------------------
    call get_params_2(t_params, krel, natyp, ipand, npotd, natomimpd, lmxspd, nfund, &
      lmpot, ncelld, irmd, nembd1, nembd, irmind, nsra, ins, nspin, ipan, ircut, lcore, &
      ncore, lmax, ntcell, lpot, nlbasis, nrbasis, nright, nleft, natomimp, atomimp, &
      imix, qbound, fcm, itdbry, irns, kpre, kshape, kte, kvmad, kxc, icc, ishift, &
      ixipol, kforce, ifunm, lmsp, imt, irc, irmin, irws, llmsp, ititle, nfu, hostimp, &
      ilm_map, imaxsh, ielast, npol, npnt1, npnt2, npnt3, itscf, scfsteps, iesemicore, &
      kaoez, iqat, noq, lly, npolsemi, n1semi, n2semi, n3semi, zrel, jwsrel, irshift, &
      mixing, lambda_xc, a, b, thetas, drdi, rmesh, zat, rmt, rmtnew, rws, emin, emax, &
      tk, alat, efold, chrgold, cmomhost, conc, gsh, ebotsemi, emusemi, tksemi, vins, &
      visp, rmrel, drdirel, vbc, fsold, r2drdirel, ecore, ez, wez, txc, linterface, &
      lrhosym, ngshd, naez, irid, nspotd, iemxd, special_straight_mixing)


    ! -------------------------------------------------------------------------
    ! Reading the density parameters stored in t_params
    ! -------------------------------------------------------------------------
    call read_density(t_params, rho2ns, r2nef, rhoc, denef, denefat, espv, ecore, &
      idoldau, lopt, eu, edc, chrgsemicore, rhoorb, ecorerel, nkcore, kapcore, krel, &
      natyp, npotd, irmd, lmpot, lmaxd1)

    ! -------------------------------------------------------------------------
    ! End read in variables
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    ! Setting up constants
    ! -------------------------------------------------------------------------
    fpi = 4.0d0*pi
    rfpi = sqrt(fpi)
    rmsav0 = 1.0d10
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Setting dummy argument LSMEAR to allow compatibility with IMPURITY
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lsmear = 0

    icont = 1
    ipf = 1337
    nspin = 2*krel + (1-krel)*nspin
    idosemicore = 0
    if (use_semicore) idosemicore = 1
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!  ITERATION BEGIN  !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    itscf = itscf + 1              ! initialised to 0 in main0
    t_inc%i_iteration = itscf

    write (1337, '(/,79("*"))')
    write (1337, '(19X,A,I3,A,I3,A)') '****** ITERATION : ', itscf, ' OUT OF ', scfsteps, ' ******'
    write (1337, '(79("*"),/)')

    if (imix>=3) then
      open (iobroy, file='broy_io.unformatted', form='unformatted', status='unknown')
      open (iobroy+2, file='broy_io2.unformatted', form='unformatted', status='unknown')
    end if
    ! -------------------------------------------------------------------------
    ! the next four lines may not always work
    ! -------------------------------------------------------------------------
    nshell(0) = natyp
    do i1 = 1, natyp
      nshell(i1) = 1
    end do
    ! -------------------------------------------------------------------------
    ! Determine total charge density expanded in spherical harmonics
    ! -------------------------------------------------------------------------
    if (print_program_flow) write (1337, *) '>>> RHOTOTB'
    call rhototb(ipf, natyp, naez, nspin, rho2ns, rhoc, rhoorb, zat, drdi, irws, ircut, &
      nfu, llmsp, thetas, ntcell, kshape, ipan, chrgnt, itscf, nshell, noq, conc, kaoez,&
      chrgatom, irmd, nemb, lmpot)

    if (print_program_flow) write (1337, *) '<<< RHOTOTB'

    if (write_rho2ns) then     ! Bauer
      open (unit=324234, file='out_rhotot')
      do i1 = 1, natyp
        write (324234, *) '#IATOM', i1
        write (324234, '(50000F14.7)') rho2ns(:, :, i1, 1)
        if (nspin==2) write (324234, '(50000F14.7)') rho2ns(:, :, i1, 2)
      end do
      close (unit=324234)
    end if

    ! -------------------------------------------------------------------------
    ! Determine new Fermi level due to valence charge up to old Fermi level
    ! EMAX and density of states DENEF
    ! -------------------------------------------------------------------------
    if (itscf>1 .and. chrgnt*chrgold<0.0_dp .and. abs(chrgnt)>5.d-2) then
      e2shift = chrgnt/(chrgnt-chrgold)*(emax-efold)
    else
      e2shift = chrgnt/denef
    end if

    e2shift = min(abs(e2shift), 0.05_dp)*sign(1.0_dp, e2shift)
    efold = emax
    chrgold = chrgnt
    if (disable_charge_neutrality) then
      write (1337, *) 'test-opt no-neutr: Setting FERMI level shift to zero'
      e2shift = 0.0_dp
    end if
    if (slow_mixing_Efermi) then
      write (1337, *) 'test-opt slow-neu: FERMI level shift * STRMIX'
      e2shift = e2shift*mixing
    end if
    if (ishift==0) emax = emax - e2shift
    ! -------------------------------------------------------------------------
    fsemicore = 0d0
    if (idosemicore==1) then
      ! ----------------------------------------------------------------------
      ! Semicore treatment, recalculate the normalisation factor
      ! ----------------------------------------------------------------------
      if (chrgsemicore<1d-10) chrgsemicore = 1d-10
      ! Number of semicore bands
      i1 = nint(chrgsemicore)
      fsemicore = real(i1, kind=dp)/chrgsemicore*fsold
      write (1337, '(6X,"< SEMICORE > : ",/,21X,"charge found in semicore :",F10.6,/,21X,"new normalisation factor :",F20.16,/)') chrgsemicore, fsemicore
    end if
    ! -------------------------------------------------------------------------
    ! write (6,FMT=9020) EFOLD,E2SHIFT
    write (1337, fmt=110) efold, e2shift
    ! -------------------------------------------------------------------------
    ! Divided by NAEZ because the weight of each atom has been already
    ! taken into account in 1c
    ! -------------------------------------------------------------------------
    write (1337, fmt=120) emax, denef/real(naez, kind=dp)
    write (6, fmt=120) emax, denef/real(naez, kind=dp)
    write (1337, '(79("+"),/)')
    ! -------------------------------------------------------------------------
    df = 2.0d0/pi*e2shift/real(nspin, kind=dp)
    ! -------------------------------------------------------------------------
    ! ISPIN LOOP
    ! -------------------------------------------------------------------------
    do ispin = 1, nspin
      ! ----------------------------------------------------------------------
      if (kte==1) then
        do i1 = 1, natyp
          ipot = (i1-1)*nspin + ispin
          espv(0, ipot) = espv(0, ipot) - efold*chrgnt/real(nspin*naez, kind=dp)
        end do
      end if                       ! (kte.eq.1)
      ! ----------------------------------------------------------------------
      ! Get correct density
      ! ----------------------------------------------------------------------
      if (.not. (use_decimation)) then
        do i1 = 1, natyp
          do lm = 1, lmpot
            call daxpy(irc(i1), df, r2nef(1,lm,i1,ispin), 1, rho2ns(1,lm,i1,ispin), 1)
          end do
        end do
      end if
      ! ----------------------------------------------------------------------
    end do
    ! -------------------------------------------------------------------------
    ! End of ISPIN loop
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------
    ! ITERMDIR
    ! -------------------------------------------------------------------------
    if ((krel==1) .and. (relax_SpinAngle_Dirac)) then
      mvevi = t_params%mvevi
      mvevief = t_params%mvevief

      call rinit(naez, qmgam)
      do i1 = 1, naez
        itoq(1, i1) = kaoez(1, i1)
      end do
      nk = 2*lmax + 1
      nmvec = 2

      fact(0) = 1.0d0
      do i = 1, 100
        fact(i) = fact(i-1)*real(i, kind=dp)
      end do
      ! ----------------------------------------------------------------------
      if (.not. (use_decimation)) then
        do i1 = 1, natyp
          do lm = 1, nmvec
            do it = 1, 3
              mvevi(i1, it, lm) = mvevi(i1, it, lm) + e2shift*mvevief(i1, it, lm)
            end do
          end do
        end do
      end if
      ! ----------------------------------------------------------------------
      do i1 = 1, natyp
        call mdirnewang(i1, nmvec, mvevi, mvphi, mvtet, mvgam, natyp, lmax, nmvecmax)
      end do

      open (67, file='itermdir.unformatted', form='unformatted')
      call scfiterang(itscf,itoq,fact,mvphi,mvtet,mvgam,qmphi,qmtet,qmgam,naez,nk,  &
        erravang,naez,natyp,nmvecmax,lmmaxd)
      t_params%mvevi = mvevi
      t_params%mvevief = mvevief
    end if
    ! -------------------------------------------------------------------------
    ! End of ITERMDIR
    ! -------------------------------------------------------------------------

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! POTENTIAL PART
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (lrhosym) then
      call rhosymm(lmpot,nspin,1,natyp,rho2ns,ixipol,irws,ircut,ipan,kshape,natyp,  &
        irmd)
    end if

    cminst(:, :) = 0.0_dp
    cmom(:, :) = 0.0_dp
    vons(:, :, :) = 0.0_dp
    call vintras(cmom, cminst, lpot, nspin, 1, natyp, rho2ns, vons, rmesh, drdi, irws, &
      ircut, ipan, kshape, ntcell, ilm_map, ifunm, imaxsh, gsh, thetas, lmsp, lmpot, natyp)

    if (write_potential_tests) then     ! Bauer
      open (unit=786785, file='test_vintraspot')
      do i1 = 1, nspin*natyp
        write (786785, *) '# atom/spin index ', i1
        write (786785, '(50000E25.16)') vons(:, :, i1)
      end do                       ! iatom
      close (786785)
    end if

    ! -------------------------------------------------------------------------
    ! fivos     IF ( .NOT.no_madelung.AND. ( SCFSTEPS.GT.1 )
    ! fivos     &     .OR. (ICC .GT. 0 ) )THEN
    ! -------------------------------------------------------------------------
    if (linterface) then
      call vinterface(cmom, cminst, lpot, nspin, naez, natyp, vons, zat, rmesh, irws, &
        ircut, ipan, kshape, noq, kaoez, iqat, conc, chrgatom(1,1), icc, hostimp, &
        nlbasis, nleft, nrbasis, nright, cmomhost, chrgnt, vinters, naez, lmpot)
      ! ----------------------------------------------------------------------
    else
      ! ----------------------------------------------------------------------
      call vmadelblk(cmom, cminst, lpot, nspin, naez, vons, zat, rmesh, irws, ircut, &
        ipan, kshape, noq, kaoez, conc, chrgatom(1,1), icc, hostimp, vinters, nemb, &
        lmpot, natyp)
    end if

    if (write_kkrimp_input) then
      call writekkrflex(natomimp, nspin, ielast, (lpot+1)**2, alat, natyp, kshape, vbc, &
        atomimp, hostimp, noq, zat, kaoez, conc, cmom, cminst, vinters, nemb, naez)
    end if

    ! -------------------------------------------------------------------------
    ! fivos      END IF
    ! -------------------------------------------------------------------------
    if (use_spherical_potential_only) vons(1:irmd, 2:lmpot, 1:npotd) = 0.0_dp
    if (write_potential_tests) then     ! bauer
      open (unit=54633163, file='test_vpotout_inter')
      do i1 = 1, natyp*nspin
        write (54633163, *) '# atom ', i1
        write (54633163, '(50000E25.16)') vons(:, :, i1)
      end do                       ! iatom
      close (54633163)
    end if                         ! config_testflag('write_gmatonsite')

    ! -------------------------------------------------------------------------
    ! Write the CMOMS to a file
    ! -------------------------------------------------------------------------
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! In case of DECIMATION output, we store ONLY information connected
    ! with the effective CPA-medium (for the case of NO-CPA NAEZ=NATYP)
    ! hence the CMOMS are calculated site-dependent. In the same format
    ! are read in by <MAIN0> -- < CMOMSREAD >     v.popescu 01/02/2002
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (write_deci_tmat .and. (itscf==1)) then
      open (37, file='decifile', form='formatted', position='append')
      write (37, 150) naez, lmpot
      do ih = 1, naez
        write (37, *) ih
        do lm = 1, lmpot
          c00(lm) = 0.0d0
          ! ----------------------------------------------------------------
          ! Store the charge on SITE IH
          ! ----------------------------------------------------------------
          do io = 1, noq(ih)
            it = kaoez(io, ih)
            c00(lm) = c00(lm) + cmom(lm, it)*conc(it)
            if (ins/=0) c00(lm) = c00(lm) + cminst(lm, it)*conc(it)
            if (lm==1) c00(1) = c00(1) - zat(it)/rfpi*conc(it)
          end do
          ! ----------------------------------------------------------------
        end do
        write (37, 160)(c00(lm), lm=1, lmpot)
      end do
      close (37)
    end if
    ! -------------------------------------------------------------------------
    ! FORCES
    ! -------------------------------------------------------------------------
    if ((kforce==1) .and. (krel/=1)) then
      if (ins==0) then
        call forceh(cmom,flm,lpot,nspin,1,natyp,rho2ns,vons,rmesh,drdi,irws,zat)
        call force(flm,flmc,lpot,nspin,1,natyp,rhoc,vons,rmesh,drdi,irws)
      else
        call forceh(cmom,flm,lpot,nspin,1,natyp,rho2ns,vons,rmesh,drdi,imt,zat)
        call force(flm,flmc,lpot,nspin,1,natyp,rhoc,vons,rmesh,drdi,imt)
      end if
    end if
    ! -------------------------------------------------------------------------
    ! Force Calculation stops here look after VXCDRV
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    ! ENERGIES
    ! -------------------------------------------------------------------------
    if (kte==1) then
      ! Single-particle core energy
      call espcb(espc, nspin, natyp, ecore, lcore, lcoremax, ncore)
      ! Energy of the input potential: Int V(r) rho(r) d^3r
      call epotinb(epotin, nspin, natyp, rho2ns, visp, rmesh, drdi, ins, irmin, irws, &
        lpot, vins, ircut, ipan, zat)
      ! Coulomb hartree energy
      call ecoub(cmom, ecou, lpot, nspin, natyp, rho2ns, vons, zat, rmesh, drdi, irws, &
        kvmad, kshape, ircut, ipan, imaxsh, ifunm, ilm_map, ntcell, gsh, thetas, lmsp, &
        lpot)
    end if
    ! -------------------------------------------------------------------------
    ! End of calculation of the energy
    ! -------------------------------------------------------------------------
    vxcm(:, :, :) = 0.0_dp
    exc(:, :) = 0.0_dp
    call vxcdrv(exc, kte, kxc, lpot, nspin, 1, natyp, rho2ns, vxcm, rmesh, drdi, a, &
      irws, ircut, ipan, ntcell, kshape, gsh, ilm_map, imaxsh, ifunm, thetas, lmsp)

    if (use_spherical_potential_only) vons(1:irmd, 2:lmpot, 1:npotd) = 0.0_dp

    ! Recalculate XC-potential with zero spin density for magn. moment scaling
    vxcnm(:, :, :) = 0.0_dp          ! Initialize
    excnm(:, :) = 0.0_dp
    if (abs(lambda_xc-1.0_dp)>eps .and. nspin==2) then
      rho2nsnm(:, :, :, 1) = rho2ns(:, :, :, 1) ! Copy charge density
      rho2nsnm(:, :, :, 2) = 0.0_dp  ! Set spin density to zero
      call vxcdrv(excnm, kte, kxc, lpot, nspin, 1, natyp, rho2nsnm, vxcnm, rmesh, drdi, &
        a, irws, ircut, ipan, ntcell, kshape, gsh, ilm_map, imaxsh, ifunm, thetas, lmsp)
      ! Compute the EXC-difference
      excdiff = 0.0_dp
      do i1 = 1, natyp
        do lm = 0, lpot
          excdiff = excdiff + exc(lm, i1) - excnm(lm, i1)
        end do
      end do
      write (1337, *) 'LAMBDA_XC=', lambda_xc, 'EXCDIF=', excdiff
    end if

    ! Add xc-potential with magn. part weighted by lambda_xc
    vons(:, :, :) = vons(:, :, :) + lambda_xc*vxcm(:, :, :) + (1.0_dp-lambda_xc)*vxcnm(:, :, :)
    exc(:, :) = lambda_xc*exc(:, :) + (1.0_dp-lambda_xc)*excnm(:, :)


    if (write_potential_tests) then     ! bauer
      open (unit=57633263, file='test_vpotout_xc')
      do i1 = 1, natyp*nspin
        write (57633263, *) '# atom ', i1
        write (57633263, '(50000E25.16)') vons(:, :, i1)
      end do                       ! iatom
      close (57633263)
    end if                         ! config_testflag('write_gmatonsite')

    ! -------------------------------------------------------------------------
    ! FORCES
    ! -------------------------------------------------------------------------
    ! Force calculation continues here
    ! -------------------------------------------------------------------------
    if ((kforce==1) .and. (krel/=1)) then
      if (kshape==0) then
        call forcxc(flm,flmc,lpot,nspin,1,natyp,rhoc,vons,rmesh,alat,drdi,irws,0)
      else
        call forcxc(flm,flmc,lpot,nspin,1,natyp,rhoc,vons,rmesh,alat,drdi,imt,0)
      end if
    end if
    ! -------------------------------------------------------------------------
    ! Force calculation ends
    ! -------------------------------------------------------------------------
    if (ishift==2) then            ! fxf
      ! OPEN (67,FILE='vmtzero_init',FORM='formatted')                 ! fxf
      ! READ (67,*) VMT_INIT(1)                                        ! fxf
      ! CLOSE(67)                                                      ! fxf
      vmt_init(1) = 0.0_dp
      vmt_init(2) = vmt_init(1)    ! read initial muffin-tin zero                ! fxf
      ! vmt_init is needed as a common reference for potential mixing if more
      ! iterations are to be remembered, e.g., in Anderson mixing.
      ! ----------------------------------------------------------------------
      ! Shift new potential to initial muffin-tin zero                        ! fxf
      call mtzero(lmpot, natyp, conc, nspin, vons, vmt_init, zat, rmesh, drdi, imt, &
        ircut, ipan, ntcell, lmsp, ifunm, thetas, irws, e2shift, ishift, nshell, &
        linterface)      ! fxf
      open (67, file='vmtzero', form='formatted') ! fxf
      write (67, *) vmt_init(1)    ! fxf
      close (67)                   ! fxf
      ! Shift old potential to initial muffin-tin zero for correct mixing     ! fxf
      vshift = -vbc(1)             ! fxf
      call potenshift(visp, vins, natyp, nspin, ircut, irc, irmin, ntcell, imaxsh, &
        ilm_map, ifunm, lmsp, lmpot, gsh, thetas, thesme, rfpi, rmesh, kshape, vshift, &
        irmd, npotd, irmind, lmxspd) ! fxf
    else if (ishift==1) then
      ! Shift new potential to old MT-zero for correct mixing
      ! (convolution with shapes is done later)
      do ispin = 1, nspin
        do ih = 1, natyp
          ipot = nspin*(ih-1) + ispin
          irc1 = irc(ih)
          vshift = rfpi*vbc(ispin)
          vons(1:irc1, 1, ipot) = vons(1:irc1, 1, ipot) + vshift
        end do
      end do

    else                           ! fxf
      ! Before fxf, only the following call was present.
      call mtzero(lmpot, natyp, conc, nspin, vons, vbc, zat, rmesh, drdi, imt, ircut, &
        ipan, ntcell, lmsp, ifunm, thetas, irws, e2shift, ishift, nshell, linterface)
    end if                         ! fxf
    ! -------------------------------------------------------------------------
    write (1337, '(79("="),/)')
    ! -------------------------------------------------------------------------
    ! Convolute potential with shape function for next iteration
    ! -------------------------------------------------------------------------

    if (write_potential_tests) then     ! bauer
      open (unit=12633269, file='test_vpotout_shift')
      do i1 = 1, natyp*nspin
        write (12633269, *) '# atom ', i1
        write (12633269, '(50000E25.16)') vons(:, :, i1)
      end do                       ! iatom
      close (12633269)
    end if                         ! config_testflag('write_gmatonsite')

    if (kshape/=0) then
      do ispin = 1, nspin
        do i1 = 1, natyp
          ipot = nspin*(i1-1) + ispin

          if (write_potential_tests) then ! bauer
            open (unit=12642269, file='test_convol')
            write (12642269, *) '# atom ', i1

            write (12642269, *) ircut(1, i1), irc(i1), imaxsh(lmpot), ilm_map, ifunm, &
              lmpot, gsh, thetas, zat(i1), rfpi, rmesh(:, i1), vons(:, :, ipot), lmsp
            close (12642269)
          end if                   ! config_testflag('write_gmatonsite')

          call convol(ircut(1,i1), irc(i1), ntcell(i1), imaxsh(lmpot), ilm_map, ifunm, &
            lmpot, gsh, thetas, thesme, zat(i1), rfpi, rmesh(:,i1), vons(:,:,ipot), &
            vspsmdum(1,1), lmsp)
        end do
      end do
    end if

    if (write_potential_tests) then     ! bauer
      open (unit=57633269, file='test_vpotout_conv')
      do i1 = 1, natyp*nspin
        write (57633269, *) '# atom ', i1
        write (57633269, '(50000E25.16)') vons(:, :, i1)
      end do                       ! iatom
      close (57633269)
    end if                         ! config_testflag('write_gmatonsite')

    ! -------------------------------------------------------------------------
    ! Symmetrisation of the potentials
    ! -------------------------------------------------------------------------
    ! Keep only symmetric part of the potential
    if (symmetrize_potential_cubic) then
      write (1337, *) 'Keeping only symmetric part of potential:'
      write (1337, *) 'Components L = 1, 11, 21, 25, 43, 47.'
      do ipot = 1, npotd
        do lm = 1, lmpot
          if (lm/=1 .and. lm/=11 .and. lm/=21 .and. lm/=25 .and. lm/=43 .and. lm/=47) then
            do i = 1, irmd
              vons(i, lm, ipot) = 0.0_dp
            end do
          end if
        end do
      end do
    end if

    if (symmetrize_potential_madelung) then
      ! declarations needed:
      ! real (kind=dp) AVMAD(LMPOT,LMPOT),BVMAD(LMPOT)
      ! INTEGER LRECABMAD,I2,IREC
      ! LOGICAL LPOTSYMM(NATYP,LMPOT)
      lrecabmad = wlength*2*lmpot*lmpot + wlength*2*lmpot
      if (write_madelung_file) open (69, access='direct', recl=lrecabmad, file='abvmad.unformatted', form='unformatted')
      do i1 = 1, natyp
        do lm = 1, lmpot
          lpotsymm(i1, lm) = .false.
        end do
        do i2 = 1, natyp
          irec = i2 + naez*(i1-1)
          if (write_madelung_file) then
            read (69, rec=irec) avmad, bvmad
          else
            bvmad = t_madel%bvmad(irec,:)
          end if
          do lm = 1, lmpot
            if (abs(bvmad(lm))>1d-10) lpotsymm(i1, lm) = .true.
          end do
        end do
        do lm = 1, lmpot
          if (lpotsymm(i1,lm)) then
            write (1337, *) 'atom ', i1, 'lm = ', lm, ' contribution used'
          else
            do ispin = 1, nspin
              ipot = nspin*(i1-1) + ispin
              do ir = 1, irmd
                vons(ir, lm, ipot) = 0.0d0
              end do
            end do
          end if
        end do
      end do
      if (write_madelung_file) close (69)
    end if

    if (write_potential_tests) then     ! bauer
      open (unit=54633563, file='test_vpotout')
      do i1 = 1, natyp*nspin
        write (54633563, *) '# atom ', i1
        write (54633563, '(50000E25.16)') vons(:, :, i1)
      end do                       ! iatom
      close (54633563)
    end if                         ! config_testflag('write_gmatonsite')

    ! for simulasa:
    if (simulate_asa) then
      do ias = 1, npotd
        do ilm_mapp = 1, lmpot
          do ipos = 1, irmd
            if (ilm_mapp/=1) then
              vons(ipos, ilm_mapp, ias) = 0.0_dp
            end if
          end do
        end do
      end do
    end if

    ! -------------------------------------------------------------------------
    ! Final construction of the potentials (straight mixing)
    ! -------------------------------------------------------------------------
    mix = mixing
    if (special_straight_mixing==1) mix = mixing/real(1+mod(itscf,2), kind=dp)
    if (special_straight_mixing==2) then
      mix = mixing/(1.0d0+1.0d+3*abs(chrgnt)/real(naez*nspin,kind=dp))
    end if
    write (1337, *) 'MIXSTR', mix
    call mixstr(rmsavq, rmsavm, ins, lpot, lmpot, 0, nshell, 1, natyp, conc, nspin, &
      itscf, rfpi, fpi, ipf, mix, fcm, irc, irmin, rmesh, drdi, vons, visp, vins, &
      vspsmdum, vspsmdum, lsmear)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of  POTENTIAL PART
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! -------------------------------------------------------------------------
    if (itscf/=1) rmsav0 = 1.0d2*max(rmsavq, rmsavm)

    write (1337, fmt=140) mix
    write (1337, '(79("="),/)')
    ! -------------------------------------------------------------------------
    if (max(rmsavq,rmsavm)<qbound) then
      t_inc%i_iteration = t_inc%n_iteration
    else
      ! ----------------------------------------------------------------------
      ! Potential mixing procedures: Broyden or Andersen updating schemes
      ! ----------------------------------------------------------------------
      if (imix>=3) then
        call brydbm(visp, vons, vins, vspsmdum, vspsmdum, ins, lmpot, rmesh, drdi, mix, &
          conc, irc, irmin, nspin, 1, natyp, itdbry, imix, iobroy, ipf, lsmear)
      end if
      ! ----------------------------------------------------------------------
      ! Reset to start new iteration
      ! ----------------------------------------------------------------------
      do i = 1, nspin*natyp
        it = i
        if (nspin==2) it = (i+1)/2

        irc1 = irc(it)
        call dcopy(irc1, vons(1,1,i), 1, visp(1,i), 1)

        if ((ins/=0) .and. (lpot>0)) then
          irmin1 = irmin(it)
          do lm = 2, lmpot
            do j = irmin1, irc1
              vins(j, lm, i) = vons(j, lm, i)
            end do
            sum = 0.0d0
            do ir = irmin1, irc1
              rv = vins(ir, lm, i)*rmesh(ir, it)
              sum = sum + rv*rv*drdi(ir, it)
            end do
            if (sqrt(sum)<qbound) then
              do j = irmin1, irc1
                vins(j, lm, i) = 0.0d0
              end do
            end if
          end do
        end if
      end do
      ! ----------------------------------------------------------------------
    end if
    ! -------------------------------------------------------------------------
    rewind 11

    efnew = emax
    if (use_rigid_Efermi .or. use_decimation) efnew = efold

    if (ishift==2) then            ! Shift mixed potential to new muffin-tin zero       ! fxf
      vbc(1) = vbc(1) + e2shift    ! fxf
      vbc(2) = vbc(1)              ! fxf
      vshift = vbc(1)              ! fxf
      call potenshift(visp, vins, natyp, nspin, ircut, irc, irmin, ntcell, imaxsh, &
        ilm_map, ifunm, lmsp, lmpot, gsh, thetas, thesme, rfpi, rmesh, kshape, vshift, &
        irmd, npotd, irmind, lmxspd) ! fxf
      write (1337, *) 'New VMT ZERO:', vbc(1) ! fxf
    end if                         ! fxf

    call rites(11, 1, natyp, nspin, zat, alat, rmt, rmtnew, rws, ititle, rmesh, drdi, &
      visp, irws, a, b, txc, kxc, ins, irns, lpot, vins, qbound, irc, kshape, efnew, &
      vbc, ecore, lcore, ncore, ecorerel, nkcore, kapcore, lmpot)
    close (11)
    ! -------------------------------------------------------------------------
    ! ENERGIES calculation
    ! -------------------------------------------------------------------------
    if (kte==1) then
      call etotb1(ecou, epotin, espc, espv, exc, kpre, lmax, lpot, lcoremax, nspin, &
        natyp, nshell(1), conc, idoldau, lopt, eu, edc)
    end if
    ! -------------------------------------------------------------------------
    ! End of ENERGIES calculation
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    ! CONVERGENCY TESTS
    ! -------------------------------------------------------------------------
    if (search_Efermi .and. (abs(e2shift)<1d-8)) then
      t_inc%i_iteration = t_inc%n_iteration
      icont = 0
      go to 100
    end if
    ! -------------------------------------------------------------------------
    if (max(rmsavq,rmsavm)<qbound) then
      write (6, '(17X,A)') '++++++ SCF ITERATION CONVERGED ++++++'
      write (6, '(79("*"))')
      icont = 0
      go to 100
      ! ----------------------------------------------------------------------
    else
      ! ---------------------------------------------------------------------
      if (max(rmsavq,rmsavm)>rmsav0) then
        write (6, *) 'ITERATION DIVERGED ---'
        icont = 0
        go to 100
      end if
      ! ----------------------------------------------------------------------
      if (itscf>=scfsteps) then
        t_inc%i_iteration = t_inc%n_iteration
        write (6, '(12X,A)') '++++++ NUMBER OF SCF STEPS EXHAUSTED ++++++'
        write (6, '(79("*"))')
        icont = 0
        go to 100
      end if
    end if
    ! -------------------------------------------------------------------------
100 continue                       ! jump mark
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!    ITERATION END    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! -------------------------------------------------------------------------
    ! Update energy contour
    ! -------------------------------------------------------------------------
    if (icont==1) then
      call epathtb(ez, dez, emax, ielast, iesemicore, idosemicore, emin, emax, tk, &
        npol, npnt1, npnt2, npnt3, ebotsemi, emusemi, tksemi, npolsemi, n1semi, n2semi, &
        n3semi, iemxd)
      do ie = 1, ielast
        wez(ie) = -2.0_dp/pi*dez(ie)
        if (ie<=iesemicore) wez(ie) = wez(ie)*fsemicore
      end do
      write (1337, '(79("="))')
    end if
    ! -------------------------------------------------------------------------
    ! Convert VISP potential to the relativistic form VTREL,BTREL.
    ! -------------------------------------------------------------------------
    if (krel==1) then
      call relpotcvt(2, visp, zat, rmesh, drdi, ircut, vtrel, btrel, zrel, rmrel, &
        jwsrel, drdirel, r2drdirel, irshift, ipand, irmd, npotd, natyp)
    end if

    ! -------------------------------------------------------------------------
    ! Write out information for the next iteration
    ! -------------------------------------------------------------------------
    ! -------------------------------------------------------------------------
    ! New_energy_mesh
    ! -------------------------------------------------------------------------
    call save_emesh(ielast, ez, wez, emin, emax, iesemicore, fsemicore, npol, tk, &
      npnt1, npnt2, npnt3, ebotsemi, emusemi, tksemi, npolsemi, n1semi, n2semi, &
      n3semi, iemxd, t_params)
    ! -------------------------------------------------------------------------
    ! Output_potential
    ! -------------------------------------------------------------------------
    call save_scfinfo(t_params, vins, visp, ecore, vbc, rmrel, drdirel, r2drdirel, zrel, & 
      jwsrel, irshift, vtrel, btrel, itscf, scfsteps, efold, chrgold, cmomhost, krel, &
      irmind, irmd, lmpot, nspotd, natyp, npotd, nembd1)
    ! -------------------------------------------------------------------------
110 format ('                old', ' E Fermi ', f14.10, ' Delta E_F = ', e16.8)
120 format ('                new', ' E FERMI ', f14.10, '  DOS(E_F) = ', f12.6)
130 format (19x, ' TIME IN ITERATION : ', f9.2, /, 79('*'))
140 format (20x, 'mixing factor used : ', 1p, d12.2)
150 format ('CMOMC', 2i6)
160 format (4d22.14)

    ! Deallocate arrays
    i_all = -product(shape(vxcm))*kind(vxcm)
    deallocate (vxcm, stat=i_stat)
    call memocc(i_stat, i_all, 'vxcm', 'main2')
    i_all = -product(shape(vxcnm))*kind(vxcnm)
    deallocate (vxcnm, stat=i_stat)
    call memocc(i_stat, i_all, 'vxcnm', 'main2')
    i_all = -product(shape(thetas))*kind(thetas)
    deallocate (r2nef, stat=i_stat)
    call memocc(i_stat, i_all, 'R2NEF', 'main2')
    i_all = -product(shape(rho2ns))*kind(rho2ns)
    deallocate (rho2ns, stat=i_stat)
    call memocc(i_stat, i_all, 'RHO2NS', 'main2')
    i_all = -product(shape(rho2nsnm))*kind(rho2nsnm)
    deallocate (rho2nsnm, stat=i_stat)
    call memocc(i_stat, i_all, 'RHO2NSNM', 'main2')
    deallocate (vons, stat=i_stat)
    i_all = product(shape(vons))*kind(vons)
    call memocc(i_stat, i_all, 'VONS', 'main2')

  end subroutine main2


  ! ----------------------------------------------------------------------------
  !> Summary: Adds a constant (=VSHIFT) to the potentials of atoms                   
  !> Author:                                                                         
  !> Category: potential, KKRhost                                                    
  !> Deprecated: False                                                               
  !> Adds a constant (=VSHIFT) to the potentials of atoms                            
  ! ----------------------------------------------------------------------------
  subroutine potenshift(visp, vins, natyp, nspin, ircut, irc, irmin, ntcell, imaxsh, &
    ilm_map, ifunm, lmsp, lmpot, gsh, thetas, thesme, rfpi, rmesh, kshape, vshift, &
    irmd, npotd, irmind, lmxspd)

    use :: global_variables, only: nspotd, ipand, ngshd, irid, nfund, ncelld
    use :: mod_rinit, only: rinit 
    use :: mod_convol, only: convol
    implicit none

    ! .. Input variables
    integer, intent (in) :: irmd   !! Maximum number of radial points
    integer, intent (in) :: lmpot  !! (LPOT+1)**2
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: npotd  !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: lmxspd !! (2*LPOT+1)**2
    integer, intent (in) :: irmind !! IRMD-IRNSD
    integer, intent (in) :: kshape !! Exact treatment of WS cell
    real (kind=dp), intent (in) :: rfpi
    real (kind=dp), intent (in) :: vshift
    integer, dimension (natyp), intent (in) :: irc !! R point for potential cutting
    integer, dimension (natyp), intent (in) :: irmin !! Max R for spherical treatment
    integer, dimension (natyp), intent (in) :: ntcell !! Index for WS cell
    integer, dimension (0:lmpot), intent (in) :: imaxsh
    integer, dimension (ngshd, 3), intent (in) :: ilm_map
    integer, dimension (natyp, lmxspd), intent (in) :: lmsp
    integer, dimension (natyp, lmxspd), intent (in) :: ifunm
    integer, dimension (0:ipand, natyp), intent (in) :: ircut !! R points of panel borders
    real (kind=dp), dimension (ngshd), intent (in) :: gsh
    real (kind=dp), dimension (irmd, natyp), intent (in) :: rmesh !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (irid, nfund, ncelld), intent (in) :: thesme
    real (kind=dp), dimension (irid, nfund, ncelld), intent (in) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    ! .. Input/Output:
    real (kind=dp), dimension (irmd, npotd), intent (inout) :: visp !! Spherical part of the potential
    real (kind=dp), dimension (irmind:irmd, lmpot, nspotd), intent (inout) :: vins !! Non-spherical part of the potential
    ! .. Local variables
    integer :: ispin, ih, ipot, ir, lm, imt1, irc1, irmin1
    real (kind=dp), dimension (irmd) :: pshiftr
    real (kind=dp), dimension (irmd, lmpot) :: pshiftlmr

    do ih = 1, natyp
      imt1 = ircut(1, ih)
      irc1 = irc(ih)
      irmin1 = irmin(ih)
      do ispin = 1, nspin
        write (1337, *) 'SHIFTING OF THE POTENTIALS OF ATOM', ih, ' BY', vshift, 'RY.'
        ipot = nspin*(ih-1) + ispin

        call rinit(irmd*lmpot, pshiftlmr)
        call rinit(irmd, pshiftr)
        do ir = 1, irc1
          pshiftlmr(ir, 1) = vshift
        end do

        if (kshape==0) then        ! ASA
          do ir = 1, irc1
            visp(ir, ipot) = visp(ir, ipot) + pshiftlmr(ir, 1)
          end do
        else                       ! Full-potential
          call convol(imt1,irc1,ntcell(ih),imaxsh(lmpot),ilm_map,ifunm,lmpot,gsh,   &
            thetas,thesme,0.0_dp,rfpi,rmesh(1,ih),pshiftlmr,pshiftr,lmsp)

          do ir = 1, irc1
            visp(ir, ipot) = visp(ir, ipot) + pshiftlmr(ir, 1)
          end do

          do lm = 2, lmpot
            do ir = irmin1, irc1
              vins(ir, lm, ipot) = vins(ir, lm, ipot) + pshiftlmr(ir, lm)*rfpi
            end do
          end do
        end if                     ! (kshape.eq.0)
      end do
    end do

  end subroutine potenshift

end module mod_main2
