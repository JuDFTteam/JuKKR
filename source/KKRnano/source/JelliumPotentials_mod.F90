!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module JelliumPotentials_mod
!-------------------------------------------------------------------------------
!> Summary: Start potentials are interpolated from files
!> Author: Phivos Mavropoulos, Bernhard H Drittler, Paul F Baumeister, Marcel Bornemann
!> Category: KKRnano, initialization, potential, input-output, single-site, core-electrons
!-------------------------------------------------------------------------------
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: jellstart12
  
  contains
  

  !---------- Routines for creation of Jellium potentials -----------

  subroutine jellstart12(nspin, ins, natoms, z, idshape, rwscl, rmtcl, &
                         meshn, xrn, drn, irws, irns, alatnew, qbound, &
                         lpot, irmd, irnsd, nspind, &
                         atom_index, elementdatabasepath)
  
  
  ! Date: 2016-01-21  Time: 16:35:08
  
  ! ******************************************************
  ! * This subroutine reads a jellium potential from the database
  ! * file. and interpolates to the new mesh
  ! ******************************************************

!   use DimParams_mod, only: DimParams

  integer, intent(in)             :: nspin
  integer, intent(in)             :: ins
  integer, intent(in)             :: natoms
  double precision, intent(in)    :: z(:)
  integer, intent(in)             :: idshape(*) ! translates atom indices into shape indices
  double precision, intent(in)    :: rwscl(*)
  double precision, intent(in)    :: rmtcl(*)
  integer, intent(in)             :: meshn(:)
  double precision, intent(in)    :: xrn(:,:)
  double precision, intent(in)    :: drn(:,:)
  integer, intent(in)             :: irws(:)
  integer, intent(in)             :: irns(:)
  double precision, intent(in)    :: alatnew
  double precision, intent(inout) :: qbound
  integer, intent(in)             :: lpot, irmd, irnsd, nspind ! former members of 
  integer, intent(in)             :: atom_index ! ==atom_id?
  character(len=*), intent(in)    :: elementdatabasepath

  ! parameters that are no longer taken from 'inc.geometry' but depend on dims
  integer :: npotd, lmpotd, irmind, inslpd, lmxspd, irmdjj

  double precision :: efermi, rws0, br, za, zvali, einf, ar, amsh
  integer :: kshape, lcore(20), ncore
  double precision :: ea, s1, vbc(2), ecore1(20), bout, rmtout, r0, rmaxout, rmtnew
  integer :: i, ir, iri, ispin, ipot, id, lm1, irnsout, irmtout, irwsout, nr, iat, lcore1(20), irc, nz, irs1, nsec, nzvali, nc, ios
  logical, allocatable :: potlm(:)
  double precision, allocatable :: u(:), drdi(:), ecore(:), rmesh(:), vins(:,:), vm2z(:), vinsout(:,:), vm2zout(:), vm2zb(:), rout(:), vinsb(:,:), drdiout(:), work(:,:), ra(:)
  character(len=40) :: baner
  character(len=4) :: aaaa, tran
  character(len=2) :: txtc(20), suffix
  character(len=256) :: atompot
  character(len=32) :: filename
  double precision, parameter :: maxa = 1.d35, aout = 0.025d0

  character(len=4), parameter :: elem_file(0:113) = [ &
      'Vac0', 'H_01', 'He02', 'Li03', 'Be04', 'B_05', 'C_06', 'N_07', 'O_08', 'F_09', &
      'Ne10', 'Na11', 'Mg12', 'Al13', 'Si14', 'P_15', 'S_16', 'Cl17', 'Ar18', 'K_19', &
      'Ca20', 'Sc21', 'Ti22', 'V_23', 'Cr24', 'Mn25', 'Fe26', 'Co27', 'Ni28', 'Cu29', &
      'Zn30', 'Ga31', 'Ge32', 'As33', 'Se34', 'Br35', 'Kr36', 'Rb37', 'Sr38', 'Y_39', &
      'Zr40', 'Nb41', 'Mo42', 'Tc43', 'Ru44', 'Rh45', 'Pd46', 'Ag47', 'Cd48', 'In49', &
      'Sn50', 'Sb51', 'Te52', 'I_53', 'Xe54', 'Cs55', 'Ba56', 'La57', 'Ce58', 'Pr59', &
      'Nd60', 'Pm61', 'Sm62', 'Eu63', 'Gd64', 'Tb65', 'Dy66', 'Ho67', 'Er68', 'Tm69', &
      'Yb70', 'Lu71', 'Hf72', 'Ta73', 'W_74', 'Re75', 'Os76', 'Ir77', 'Pt78', 'Au79', &
      'Hg80', 'Tl81', 'Pb82', 'Bi83', 'Po84', 'At85', 'Rn68', 'Fr87', 'Ra88', 'Ac89', &
      'Th90', 'Pa91', 'U_92', 'Np93', 'Pu94', 'Am95', 'Cm96', 'Bk97', 'Cf98', 'Es99', &
      'Fm__', 'Md__', 'No__', 'Lr__', 'Rf__', 'Db__', 'Sg__', 'Bh__', 'Hs__', 'Mt__', &
      'Uun_', 'Uuu_', 'Uub_', 'NoE_' ]
  !
  ! Periodic System of elements according to Aco Z. Muradjan
  !                                                              H He  1s
  !                                                              LiBe  2s
  !                                                  B C N O F Ne      2p
  !                                                              NaMg  3s
  !                                                  AlSiP S ClAr      3p
  !                                                              K Ca  4s
  !                              ScTiV CrMnFeCoNiCuZn                  3d
  !                                                  GaGeAsSeBrKr      4p
  !                                                              RbSr  5s
  !                              Y ZrNbMoTcRuRhPdAgCd                  4d
  !                                                  InSnSbTeI Xe      5p
  !                                                              CsBa  6s
  !  LaCePrNdPmSmEuGdTbDyHoErTmYb                                      4f
  !                              LuHfTaW ReOsIrPtAuHg                  5d
  !                                                  TlPbBiPoAtRn      6p
  !                                                              FrRa  7s
  !  AcThPaU NpPuAmCmBkCfEsFmMdNo                                      5f
  !                              LrRfDbSgBhHsMtDsRgCn                  6d
  !                                                  utuqupuhusuo      7p
  !                                                              unud  8s
  !   
  !     --------------------------------------------------------------
  character(len=*), parameter :: F142="(7x,f8.4,7x,f8.4,7x,i5,7x,f8.4)"

  npotd = nspind*natoms
  lmpotd = (lpot+1)**2
  irmind = irmd - irnsd
  inslpd = (irnsd+1)*lmpotd
  lmxspd = (2*lpot+1)**2
  irmdjj = 1501
  kshape = 2 ! always full-pot calculations
  write(*, *) 'irmd', irmd

  ! allocate arrays
  allocate(u(irmd))
  allocate(drdi(irmdjj))
  allocate(ecore(20))
  allocate(rmesh(irmdjj))
  allocate(vins(irmind:irmd, lmpotd))
  allocate(vm2z(irmdjj))
  allocate(vinsout(irmind:irmd, lmpotd))
  allocate(vm2zout(irmd))
  allocate(vm2zb(irmdjj))
  allocate(rout(irmd))
  allocate(vinsb(irmd, lmpotd))
  allocate(drdiout(irmd))
  allocate(work(irmd, lmpotd))
  allocate(ra(irmd))
  allocate(potlm(lmpotd))

  write(6, *) ' ****  READING  POTENTIAL  **** '

  write(filename, fmt="(a,i7.7)") "potential.", atom_index
  open(19, file=filename, form="formatted", action='write', iostat=ios) ; assert(ios == 0)

  vins(:,:) = 0.d0
  vm2z(:) = 0.d0

  do iat = 1, natoms
    do ispin = 1, nspin
    
      potlm(:) = .false.
      ipot = nspin*(iat - 1) + ispin
      
  ! find out what atom is needed
      
      nz = int(z(iat))
      suffix = '  ' ; if (((nz >= 24 .and. nz <= 28) .or. (nz >= 57 .and. nz <= 70)) .and. ispin == 2) suffix = 's2'
      atompot = elementdatabasepath-'/'-elem_file(nz)-'.pot'-suffix
      write(6, "(9a)") ' Using database ....: ',atompot
      open(21, file=atompot, status='old', iostat=ios)
      if (ios /= 0) then
        write(6, *) ' Error in JELLSTART '
        write(6, *) ' Potential.............', elem_file(nz)
        write(6, *) ' Does not exist in the database'
        die_here("Unable to find"+elem_file(nz)+"in the database!")
      endif
      
  !           irws1 = nr
      
  ! --------------------------------------------------------------------
      
      efermi = .409241d+00 ! why?
      vbc(1:2) = .5d0
  !           read potential from jellium
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
      read(21, fmt="(3x,a40,3x,a4,3x,a4)") baner, aaaa, tran
      read(21, fmt=F142) rws0, s1, irs1, br
      read(21, fmt=F142) za, zvali, nsec, einf ! nsec is not used here
  !     calculate number of core states
      nz = int(za)
      nzvali = int(zvali) ! integer number of valence states
      nc = nz - nzvali
      selectcase (nc)
      case ( 0); ncore = 0 ! vacuum
      case ( 2); ncore = 1 ! 1s
      case ( 4); ncore = 2 ! 1s2s
      case (10); ncore = 3 ! 1s2s2p
      case (12); ncore = 4 ! 1s2s2p3s
      case (18); ncore = 5 ! 1s2s2p3s3p
      case (28); ncore = 6 ! 1s2s2p3s3p3d
      case (30); ncore = 7 ! 1s2s2p3s3p3d4s
      case (36); ncore = 8 ! 1s2s2p3s3p3d4s4p
      case (46); ncore = 9 ! 1s2s2p3s3p3d4s4p4d
      case (48); ncore = 10 ! 1s2s2p3s3p3d4s4p4d4s
      case (54); ncore = 11 ! 1s2s2p3s3p3d4s4p4d4s4p
      case (68); ncore = 12 ! 1s2s2p3s3p3d4s4p4d4s4p4f
      case (78); ncore = 13 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d
      case (80); ncore = 14 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d6s
      case (86); ncore = 15 ! 1s2s2p3s3p3d4s4p4d4s4p4f5d6s4p
      case default ! ToDo: error handling: what if nc is not in the list?
        ncore = 0
        warn(6, "no matching number of core electrons found! nc ="+nc)
      endselect ! nc
      write(6, *) '*************************************'
      write(6, *) '   Potential Interpolation Program   '
      write(6, *) '   Using the Jellium Database v1.0   '
      write(6, *) '*************************************'
      write(6, fmt="('Jellium Fermi Energy :',f10.5,' Ry')") efermi
      write(6, fmt="('Potential Atomic Number :',f7.2)") za
      write(6, fmt="('Number of Core   States :',i4)") ncore
      read(21, fmt="(i3,a2,d15.8)") (lcore(nc), txtc(nc), ecore(nc), nc=1, ncore)
      write(6, *) ' ** Position of the Core States ** '
      do i = 1, ncore
        write(6, fmt="(i3,a2,f15.6, ' Ry')") lcore(i), txtc(i), ecore(i)
        selectcase (txtc(i))
        case ('s '); lcore(i) = 0
        case ('p '); lcore(i) = 1
        case ('d '); lcore(i) = 2
        case ('f '); lcore(i) = 3
        case default ! ToDo: error handling: what if txtc(i) is not in the list?
          warn(6, "unable to convert '"-txtc(i)-"' into an ell-quantum number!")
        endselect ! txtc(i)
      enddo ! i
      ncore = ncore
      lcore1(1:ncore) = lcore(1:ncore)
      ecore1(1:ncore) = ecore(1:ncore)
      write(6, fmt="('All other states are above :',f8.4,' Ry in Energy')") einf
      write(6, *) '**********************************************'
      
      read(21, fmt="(4d15.8)") vm2z(1:irs1)
      read(21, fmt="(a4)") tran
      close (21)
  !     make mesh r0
      ar = log(s1/br + 1.d0)/(irs1 - 1.d0)
      ea = exp(ar)
      amsh = 1.d0
      rmesh(1) = 0.d0
      drdi(1) = br*ar
      do i = 2, irs1
        amsh = amsh*ea
        rmesh(i) = br*(amsh - 1.d0)
        drdi(i) = drdi(1)*amsh
      enddo ! i
      write(6, *) 'Jellium Potential Read In ',irs1,' points'
      nr = irs1
      
  ! --------------------------------------------------------------------
      
  !     the input mesh has been constructed. now construct the output mesh.
      
      id = idshape(iat) ! translate atom index into shape index
      rout(1) = 0.d0
      
      rmaxout = rwscl(id)
      rmtout  = rmtcl(id)
      irwsout = irws(iat)
      irmtout = irws(iat) - meshn(id)
      irnsout = irns(iat)  ! 22.1.12 changed from irns(id) to irns(iat)
      
      
      if (ins == 0) then
        bout = rmaxout/(exp(aout*(irwsout - 1.d0)) - 1.d0)
        do ir = 2, irwsout
          ea = exp(aout*(ir - 1.d0))
          rout(ir) = bout*(ea - 1.d0)
          drdiout(ir) = aout*bout*ea
        enddo ! ir
        do i = 1, irwsout
          if (rout(i) < rmtout) irmtout = i
        enddo ! i
        if (mod(irmtout, 2) == 0) irmtout = irmtout + 1
        rmtnew = rout(irmtout)
        rmaxout = rout(irwsout)
      else
        bout = rmtout/(exp(aout*(irmtout - 1.d0)) - 1.d0)
        do ir = 2, irmtout
          ea = exp(aout*(ir - 1.d0))
          rout(ir) = bout*(ea - 1.d0)
          drdiout(ir) = aout*bout*ea
        enddo ! ir
        do iri = 1, meshn(id)
          ir = iri + irmtout
          rout(ir) = alatnew*xrn(iri, id) ! scaling is always 1.d0
          drdiout(ir) = alatnew*drn(iri, id)
        enddo ! iri
        rmtnew = rout(irmtout)
        rmaxout = rout(irwsout)
      endif
      
  !  ok now interpolate
      
      call spline(irmdjj, rmesh, vm2z, nr, maxa, maxa, vm2zb)
      
  ! ok with spline
      
      vm2zout(1) = vm2z(1)
      do ir = 2, irwsout
        r0 = rout(ir)
        vm2zout(ir) = splint(rmesh, vm2z, vm2zb, nr, r0)
      enddo ! ir
      
      if (ins > 0) then
        irc = irwsout - irnsout
        do lm1 = 1, lmpotd
          do ir = irc, irwsout
            vinsout(ir, lm1) = 0.d0
          enddo ! ir
        enddo ! lm
      endif
      
      call ritesone12(19, ispin, za, alatnew, rmtout, rmtnew, rmaxout,  &
          rout, drdiout, vm2zout, irwsout, aout, bout, ins, irnsout,  &
          vinsout, qbound, irwsout, kshape, efermi, vbc,  &
          ecore1, lcore1, ncore, elem_file(nz), nspin, lpot, irmd, irnsd)

    enddo ! ispin
  enddo ! iat


  endsubroutine ! jellstart12



  subroutine ritesone12(ifile, is, z, alat, rmt, rmtnew, rws, r, drdi, vm2z, irws, a, b, ins, irns, vins, &
               qbound, irc, kshape, efermi, vbc, ecore, lcore, ncore, elem_name, nspin, lpot, irmd, irnsd)
  ! compare routine rites in PotentialConverter_mod
  
  
  ! Date: 2016-01-21  Time: 16:59:42
  
  ! ************************************************************************
  !      this subroutine stores in 'ifile' the necessary results
  !      (potentials e.t.c.) to start self-consistency iterations
  !      modified for the full potential case - if ins .gt. 0 there
  !       is written a different potential card
  !       if the sum of absolute values of an lm component of vins (non
  !       spher. potential) is less than the given rms error qbound this
  !       component will not be stored .
  !        (see to subroutine start , where most of the arrays are
  !         described)
  !                            modified by b. drittler  aug. 1988
  !-----------------------------------------------------------------------
  integer, intent(in)             :: ifile
  integer, intent(in)             :: is
  double precision, intent(in)    :: z
  double precision, intent(in)    :: alat
  double precision, intent(in)    :: rmt
  double precision, intent(in)    :: rmtnew
  double precision, intent(in)    :: rws
  double precision, intent(in)    :: r(:)
  double precision, intent(in)    :: drdi(:)
  double precision, intent(in)    :: vm2z(:)
  integer, intent(in)             :: irws
  double precision, intent(in)    :: a
  double precision, intent(in)    :: b
  integer, intent(in)             :: ins
  integer, intent(in)             :: irns
  double precision, intent(in)    :: qbound
  integer, intent(in)             :: irc
  integer, intent(in)             :: kshape
  double precision, intent(in)    :: efermi
  double precision, intent(in)    :: vbc(2)
  double precision, intent(in)    :: ecore(20)
  integer, intent(in)             :: lcore(20)
  integer, intent(in)             :: ncore
  character(len=4), intent(in)    :: elem_name
  integer, intent(in)             :: nspin, lpot, irmd, irnsd
  double precision, intent(in)    :: vins(irmd-irnsd:irmd,(lpot+1)**2)

  !     .. locals ..
  double precision :: rv, sm
  integer :: ic, ir, irmin, lm, lmnr, lmpot, nr
  integer, parameter :: inew=1, isave=1
  character(len=*), parameter :: F9060="(10i5)", F9070="(1p,4d20.13)"
  character(len=4) :: spin, updn  

  lmpot = (lpot+1)**2

  nr = irc ; if (kshape == 0) nr = irws

  irmin = nr - irns

  spin = ''; updn = ''
  if (nspin > 1) then
    spin = 'SPIN'
    updn = 'UP'; if (is == 1) updn = 'DOWN'
  endif ! spin
  write(ifile, fmt="(a4,' POTENTIAL ',a4,' ',a4,10x,'  exc:',a24)") elem_name, spin, updn, ''
  
  !     write(ifile, fmt="(7a4,6x,'  exc:',a24,3x,a10)") ititle(1:7), txc(kxc+1)
  write(ifile, fmt="(3f12.8)") rmt, alat, rmtnew
  write(ifile, fmt="(f10.5)") z
  write(ifile, fmt="(f10.5,2f15.10)") rws, efermi, vbc(is)
  write(ifile, fmt="(i0)") nr
  write(ifile, fmt="(2d15.8)") a, b
  write(ifile, fmt="(2i2)") ncore, inew
  do ic = 1, ncore
    write(ifile, fmt="(i5,1p,d20.11)") lcore(ic), ecore(ic)
  enddo ! ic

  if (ins == 0) then
  !--->       store only the spherically averaged potential
  !           (in mt or as - case)
  !           this is done always for the host
    if (inew == 0) then
      write(ifile, fmt="(1p,2d15.6,1p,d15.8)") (r(ir), drdi(ir), vm2z(ir), ir=1,nr)
    else  ! inew
      write(ifile, fmt="(1p,4d20.12)") vm2z(1:nr)
    endif ! inew

  else  ! ins
    
  !--->     store the full potential , but the non spherical contribution
  !         only from irns1 up to irws1 ;
  !         remember that the lm = 1 contribution is multiplied
  !         by a factor 1/sqrt(4 pi)
    
    write(ifile, fmt=F9060) nr, irns, lmpot, isave
    write(ifile, fmt=F9070) vm2z(1:nr)
    if (lpot > 0) then
      lmnr = 1
      do lm = 2, lmpot
        sm = 0.d0
        do ir = irmin, nr
          rv = vins(ir,lm)*r(ir)
          sm = sm + rv*rv*drdi(ir)
        enddo ! ir
        if (sm > qbound**2) then
          lmnr = lmnr + 1
          write(ifile, fmt=F9060) lm
          write(ifile, fmt=F9070) vins(irmin:nr,lm)
        endif
      enddo ! lm
      if (lmnr < lmpot) write(ifile, fmt=F9060) isave ! write a one to mark the end
    endif
    
  endif ! ins
  
  endsubroutine ! ritesone12


  !***********************************************************************
  
  
  ! Date: 2016-01-12  Time: 14:48:44

  subroutine spline(nmax, x, y, n, yp1, ypn, y2)
    integer, intent(in)             :: nmax
    double precision, intent(in)    :: x(nmax)
    double precision, intent(in)    :: y(nmax)
    integer, intent(in)             :: n
    double precision, intent(in)    :: yp1
    double precision, intent(in)    :: ypn
    double precision, intent(out)   :: y2(nmax)

    ! Given arrays x(1:n) and  y(1:n) containing a tabulated function, 
    ! i.e., y i = f(xi), with x1<x2<...<xN , and given values yp1 and ypn
    ! for the 1rst derivative of the interpolating function at points
    ! 1 and n, respectively, this routine returns an array y2(1:n) of
    ! length n which contains the second derivatives of the interpolating
    ! function at the tabulated points xi.
    ! If yp1 and/or ypn are equal to 1.e30 or larger, the routine is
    ! signaled to set the corresponding boundary condition for a natural
    ! spline, with zero second derivative on that boundary.
    ! Parameter: NMAX is the largest anticipated value of n.
    ! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
    integer :: i, k
    double precision :: p, qn, sig, un, u(nmax)
    double precision, parameter :: boundary = 0.99d30

    if (n > nmax) die_here("spline: impossible "+n+"= n > nmax ="+nmax)
    if (yp1 > boundary) then
    ! the lower boundary condition is set either to be "natural"
      y2(1) = 0.d0
      u(1) = 0.d0
    else
    ! or else to have a specified first derivative.
      y2(1) = -0.5d0
      u(1)=(3.d0/(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
    endif

    do i = 2, n-1
    ! this is the decomposition loop of the tridiagonal algorithm. y2 and u
    ! are used for temporary storage of the decomposed factors.
      sig = (x(i) - x(i-1))/(x(i+1) - x(i-1))
      p = sig*y2(i-1) + 2.d0
      y2(i) = (sig - 1.d0)/p
      u(i) = (6.d0*((y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1))/(x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig*u(i-1))/p
    enddo ! i

    if (ypn > boundary) then
    ! the upper boundary condition is set either to be "natural"
      qn = 0.d0
      un = 0.d0
    else
    ! or else to have a specified first derivative.
      qn = 0.5d0
      un = (3.d0/(x(n) - x(n-1)))*(ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))
    endif
    y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.d0)
    do k = n-1, 1, -1
    ! this is the backsubstitution loop of the tridiagonal algorithm.
      y2(k) = y2(k)*y2(k+1) + u(k)
    enddo ! k

  endsubroutine ! spline


  !***********************************************************************
  
  
  ! Date: 2016-01-12  Time: 14:48:43

  double precision function splint(xa, ya, y2a, n, x, yderiv) result(y)
    double precision, intent(in)    :: xa(*)
    double precision, intent(in)    :: ya(*)
    double precision, intent(in)    :: y2a(*)
    integer, intent(in)             :: n
    double precision, intent(in)    :: x
    double precision, intent(out), optional :: yderiv

    ! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
    ! function (with the xai's in order), and given the array y2a(1:n), which
    ! is the output from spline above, and given a value of x, this routine
    ! returns a cubic-spline interpolated value y and the derivative yderiv.
    ! Taken from "Numerical Recipes in Fortran 77", W.H.Press et al.
    integer :: k, khi, klo
    double precision :: a, b, h
    ! We will  nd the right place in the table by means of bisection.
    ! This is optimal if sequential calls to this routine are at random
    ! values of x. If sequential calls are in order, and closely
    ! spaced, one would do better to store previous values of
    ! klo and khi and test if they remain appropriate on the
    ! next call.
    klo = 1
    khi = n
    do while (khi - klo > 1)
      k = (khi + klo)/2
      if (xa(k) > x) then
        khi = k
      else
        klo = k
      endif
    enddo ! while
    ! klo and khi now bracket the input value of x.
    h = xa(khi) - xa(klo)
    ! the xa's must be distinct.
    if (h == 0.d0) die_here("bad xa input in splint")
    ! cubic spline polynomial is now evaluated.
    a = (xa(khi) - x)/h
    b = (x - xa(klo))/h
    y = a*ya(klo) + b*ya(khi) + ((a*(a*a - 1.d0))*y2a(klo) + (b*(b*b - 1.d0))*y2a(khi))*(h*h)/6.d0
    if (present(yderiv)) &
    yderiv = (ya(khi) - ya(klo))/h - ((3.d0*a*a - 1.d0)*y2a(klo) - (3.d0*b*b - 1.d0)*y2a(khi))*h/6.d0

  endfunction ! splint

endmodule ! JelliumPotentials_mod
