!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Subroutine to handle allocation/deallocation of arrays
!> Author: Jonathan Chico 
!> Module to handle the allocation of arrays to be later distributed it aims to bring 
!> modularity to the memory management. Of this way it _should_ be easier to ensure
!> that only arrays that are needed for a given task are allocated/deallocated
!------------------------------------------------------------------------------------
!> @todo The number of arrays in the misc section should be reduced, and they should
!> be located in the appropriate routines
!> @endtodo 
!------------------------------------------------------------------------------------
module memoryhandling

  use :: mod_profiling
  use :: mod_datatypes, only: dp

  implicit none
  private :: dp

contains


  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the unit cell.
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, geometry, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe the
  !> unit cell. 
  !-------------------------------------------------------------------------------  
  subroutine allocate_cell(flag,naez,nemb,natyp,cls,imt,irws,irns,ntcell,refpot,kfg,&
    kaoez,rmt,zat,rws,mtfac,rmtref,rmtrefat,rmtnew,rbasis,lmxc,fpradius)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: naez   !! number of atoms in unit cell
    integer, intent (in) :: nemb   !! number of 'embedding' positions
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    integer, dimension (:), allocatable, intent (inout) :: cls !! Cluster around atomic sites
    integer, dimension (:), allocatable, intent (inout) :: imt !! R point at MT radius
    integer, dimension (:), allocatable, intent (inout) :: irws !! R point at WS radius
    integer, dimension (:), allocatable, intent (inout) :: irns !! Position of atoms in the unit cell in units of bravais vectors
    integer, dimension (:), allocatable, intent (inout) :: lmxc
    integer, dimension (:), allocatable, intent (inout) :: ntcell !! Index for WS cell
    integer, dimension (:), allocatable, intent (inout) :: refpot !! Ref. pot. card at position
    integer, dimension (:, :), allocatable, intent (inout) :: kfg
    integer, dimension (:, :), allocatable, intent (inout) :: kaoez !! atom types located at a given site
    real (kind=dp), dimension (:), allocatable, intent (inout) :: rmt !! Muffin-tin radius of true system
    real (kind=dp), dimension (:), allocatable, intent (inout) :: zat !! Nuclear charge
    real (kind=dp), dimension (:), allocatable, intent (inout) :: rws !! Wigner Seitz radius
    real (kind=dp), dimension (:), allocatable, intent (inout) :: mtfac !! Scaling factor for radius  MT
    real (kind=dp), dimension (:), allocatable, intent (inout) :: rmtref !! Muffin-tin radius of reference system
    real (kind=dp), dimension (:), allocatable, intent (inout) :: rmtrefat
    real (kind=dp), dimension (:), allocatable, intent (inout) :: rmtnew !! Adapted muffin-tin radius
    real (kind=dp), dimension (:), allocatable, intent (inout) :: fpradius !! R point at which full-potential treatment starts
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: rbasis !! Position of atoms in the unit cell in units of bravais vectors


    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then
      allocate (refpot(naez+nemb), stat=i_stat)
      call memocc(i_stat, product(shape(refpot))*kind(refpot), 'REFPOT', 'allocate_cell')
      refpot = 0
      allocate (rmtrefat(naez+nemb), stat=i_stat)
      call memocc(i_stat, product(shape(rmtrefat))*kind(rmtrefat), 'RMTREFAT', 'allocate_cell')
      rmtrefat = -1.e0_dp          ! Signals the need for later calculation
      allocate (kaoez(natyp,naez+nemb), stat=i_stat)
      call memocc(i_stat, product(shape(kaoez))*kind(kaoez), 'KAOEZ', 'allocate_cell')
      kaoez = 0
      allocate (cls(naez+nemb), stat=i_stat)
      call memocc(i_stat, product(shape(cls))*kind(cls), 'CLS', 'allocate_cell')
      cls = 1
      allocate (rbasis(3,naez+nemb), stat=i_stat)
      call memocc(i_stat, product(shape(rbasis))*kind(rbasis), 'RBASIS', 'allocate_cell')
      rbasis = 0.e0_dp
      allocate (mtfac(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(mtfac))*kind(mtfac), 'MTFAC', 'allocate_cell')
      mtfac = 0.e0_dp
      allocate (zat(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(zat))*kind(zat), 'ZAT', 'allocate_cell')
      zat = -1.e0_dp               ! Negative value signals read-in from pot-file
      allocate (ntcell(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(ntcell))*kind(ntcell), 'NTCELL', 'allocate_cell')
      ntcell = 0
      allocate (rmtref(naez), stat=i_stat)
      call memocc(i_stat, product(shape(rmtref))*kind(rmtref), 'RMTREF', 'allocate_cell')
      rmtref = -1.e0_dp
      allocate (irns(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(irns))*kind(irns), 'IRNS', 'allocate_cell')
      irns = -1                    ! Negative value signals to use FPRADIUS
      allocate (rws(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rws))*kind(rws), 'RWS', 'allocate_cell')
      rws = 0.e0_dp
      allocate (irws(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(irws))*kind(irws), 'IRWS', 'allocate_cell')
      irws = 0
      allocate (rmt(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rmt))*kind(rmt), 'RMT', 'allocate_cell')
      rmt = 0.e0_dp
      allocate (rmtnew(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rmtnew))*kind(rmtnew), 'RMTNEW', 'allocate_cell')
      rmtnew = 0.e0_dp
      allocate (imt(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(imt))*kind(imt), 'IMT', 'allocate_cell')
      imt = 0
      allocate (kfg(4,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(kfg))*kind(kfg), 'KFG', 'allocate_cell')
      kfg = 0
      allocate (lmxc(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(lmxc))*kind(lmxc), 'LMXC', 'allocate_cell')
      lmxc = 0
      allocate (fpradius(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(fpradius))*kind(fpradius), 'FPRADIUS', 'allocate_cell')
      fpradius = -1.e0_dp          ! Negative value signals to use IRNS from pot-file (sub. startb1)
    else
      if (allocated(zat)) then
        i_all = -product(shape(zat))*kind(zat)
        deallocate (zat, stat=i_stat)
        call memocc(i_stat, i_all, 'ZAT', 'allocate_cell')
      end if
      if (allocated(ntcell)) then
        i_all = -product(shape(ntcell))*kind(ntcell)
        deallocate (ntcell, stat=i_stat)
        call memocc(i_stat, i_all, 'NTCELL', 'allocate_cell')
      end if
      if (allocated(ntcell)) then
        i_all = -product(shape(ntcell))*kind(ntcell)
        deallocate (ntcell, stat=i_stat)
        call memocc(i_stat, i_all, 'NTCELL', 'allocate_cell')
      end if
      if (allocated(rmtref)) then
        i_all = -product(shape(rmtref))*kind(rmtref)
        deallocate (rmtref, stat=i_stat)
        call memocc(i_stat, i_all, 'RMTREF', 'allocate_cell')
      end if
      if (allocated(irns)) then
        i_all = -product(shape(irns))*kind(irns)
        deallocate (irns, stat=i_stat)
        call memocc(i_stat, i_all, 'IRNS', 'allocate_cell')
      end if
      if (allocated(irws)) then
        i_all = -product(shape(irws))*kind(irws)
        deallocate (irws, stat=i_stat)
        call memocc(i_stat, i_all, 'IRWS', 'allocate_cell')
      end if
      if (allocated(rws)) then
        i_all = -product(shape(rws))*kind(rws)
        deallocate (rws, stat=i_stat)
        call memocc(i_stat, i_all, 'RWS', 'allocate_cell')
      end if
      if (allocated(rmt)) then
        i_all = -product(shape(rmt))*kind(rmt)
        deallocate (rmt, stat=i_stat)
        call memocc(i_stat, i_all, 'RMT', 'allocate_cell')
      end if
      if (allocated(rmtnew)) then
        i_all = -product(shape(rmtnew))*kind(rmtnew)
        deallocate (rmtnew, stat=i_stat)
        call memocc(i_stat, i_all, 'RMTNEW', 'allocate_cell')
      end if
      if (allocated(imt)) then
        i_all = -product(shape(imt))*kind(imt)
        deallocate (imt, stat=i_stat)
        call memocc(i_stat, i_all, 'IMT', 'allocate_cell')
      end if
      if (allocated(kfg)) then
        i_all = -product(shape(kfg))*kind(kfg)
        deallocate (kfg, stat=i_stat)
        call memocc(i_stat, i_all, 'KFG', 'allocate_cell')
      end if
      if (allocated(kaoez)) then
        i_all = -product(shape(kaoez))*kind(kaoez)
        deallocate (kaoez, stat=i_stat)
        call memocc(i_stat, i_all, 'KAOEZ', 'allocate_cell')
      end if
      if (allocated(refpot)) then
        i_all = -product(shape(refpot))*kind(refpot)
        deallocate (refpot, stat=i_stat)
        call memocc(i_stat, i_all, 'REFPOT', 'allocate_cell')
      end if
      if (allocated(rmtrefat)) then
        i_all = -product(shape(rmtrefat))*kind(rmtrefat)
        deallocate (rmtrefat, stat=i_stat)
        call memocc(i_stat, i_all, 'RMTREFAT', 'allocate_cell')
      end if
      if (allocated(lmxc)) then
        i_all = -product(shape(lmxc))*kind(lmxc)
        deallocate (lmxc, stat=i_stat)
        call memocc(i_stat, i_all, 'LMXC', 'allocate_cell')
      end if
      if (allocated(fpradius)) then
        i_all = -product(shape(fpradius))*kind(fpradius)
        deallocate (fpradius, stat=i_stat)
        call memocc(i_stat, i_all, 'FPRADIUS', 'allocate_cell')
      end if
    end if

  end subroutine allocate_cell

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the left and right host for the calculation of slabs    
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, geometry, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe the 
  !> left and right host for the calculation of slabs  
  !-------------------------------------------------------------------------------  
  subroutine allocate_semi_inf_host(flag, nemb, tleft, tright)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: nemb   !! number of 'embedding' positions
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: tleft  !! Vectors of the basis for the left host
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: tright !! vectors of the basis for the right host

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then

      allocate (tleft(3,nemb+1), stat=i_stat)
      call memocc(i_stat, product(shape(tleft))*kind(tleft), 'TLEFT', 'allocate_semi_inf_host')
      tleft = 0.e0_dp
      allocate (tright(3,nemb+1), stat=i_stat)
      call memocc(i_stat, product(shape(tright))*kind(tright), 'TRIGHT', 'allocate_semi_inf_host')
      tright = 0.e0_dp
    else
      if (allocated(tleft)) then
        i_all = -product(shape(tleft))*kind(tleft)
        deallocate (tleft, stat=i_stat)
        call memocc(i_stat, i_all, 'TLEFT', 'allocate_semi_inf_host')
      end if
      if (allocated(tright)) then
        i_all = -product(shape(tright))*kind(tright)
        deallocate (tright, stat=i_stat)
        call memocc(i_stat, i_all, 'TRIGHT', 'allocate_semi_inf_host')
      end if
    end if


  end subroutine allocate_semi_inf_host

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the potential
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, potential, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe the 
  !> potential
  !-------------------------------------------------------------------------------  
  subroutine allocate_potential(flag,irm,natyp,npotd,ipand,nfund,lmxspd,lmpot,      &
    irmind,nspotd,nfu,irc,ncore,irmin,lmsp,lmsp1,ircut,lcore,llmsp,ititle,      &
    visp,ecore,vins)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: irm
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    integer, intent (in) :: npotd  !! 2*NATYP
    integer, intent (in) :: ipand
    integer, intent (in) :: nfund
    integer, intent (in) :: lmxspd
    integer, intent (in) :: lmpot
    integer, intent (in) :: irmind
    integer, intent (in) :: nspotd
    integer, dimension (:), allocatable, intent (inout) :: nfu
    integer, dimension (:), allocatable, intent (inout) :: irc !! R point for potential cutting
    integer, dimension (:), allocatable, intent (inout) :: ncore !! Number of core states
    integer, dimension (:), allocatable, intent (inout) :: irmin !! Max R for spherical treatment
    integer, dimension (:, :), allocatable, intent (inout) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (:, :), allocatable, intent (inout) :: lmsp1
    integer, dimension (:, :), allocatable, intent (inout) :: ircut !! R points of panel borders
    integer, dimension (:, :), allocatable, intent (inout) :: lcore !! Angular momentum of core states
    integer, dimension (:, :), allocatable, intent (inout) :: llmsp !! lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
    integer, dimension (:, :), allocatable, intent (inout) :: ititle
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: visp !! Spherical part of the potential
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: ecore  !! Core energies
    real (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: vins  !! Non-spherical part of the potential

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then

      allocate (irc(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(irc))*kind(irc), 'IRC', 'allocate_potential')
      irc = 0
      allocate (ircut(0:ipand,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(ircut))*kind(ircut), 'IRCUT', 'allocate_potential')
      ircut = 0
      allocate (irmin(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(irmin))*kind(irmin), 'IRMIN', 'allocate_potential')
      irmin = 0
      allocate (lcore(20,npotd), stat=i_stat)
      call memocc(i_stat, product(shape(lcore))*kind(lcore), 'LCORE', 'allocate_potential')
      lcore = 0
      allocate (ecore(20,npotd), stat=i_stat)
      call memocc(i_stat, product(shape(ecore))*kind(ecore), 'ECORE', 'allocate_potential')
      ecore = 0.e0_dp
      allocate (ncore(npotd), stat=i_stat)
      call memocc(i_stat, product(shape(ncore))*kind(ncore), 'NCORE', 'allocate_potential')
      ncore = 0
      allocate (lmsp1(lmxspd,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(lmsp1))*kind(lmsp1), 'LMSP1', 'allocate_potential')
      lmsp1 = 0
      allocate (llmsp(natyp,nfund), stat=i_stat)
      call memocc(i_stat, product(shape(llmsp))*kind(llmsp), 'LLMSP', 'allocate_potential')
      llmsp = 0
      allocate (lmsp(natyp,lmxspd), stat=i_stat)
      call memocc(i_stat, product(shape(lmsp))*kind(lmsp), 'LMSP', 'allocate_potential')
      lmsp = 0
      allocate (ititle(20,npotd), stat=i_stat)
      call memocc(i_stat, product(shape(ititle))*kind(ititle), 'ITITLE', 'allocate_potential')
      ititle = 0
      allocate (nfu(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(nfu))*kind(nfu), 'NFU', 'allocate_potential')
      nfu = 0
      allocate (vins(irmind:irm,lmpot,nspotd), stat=i_stat)
      call memocc(i_stat, product(shape(vins))*kind(vins), 'VINS', 'allocate_potential')
      vins = 0.e0_dp
      allocate (visp(irm,npotd), stat=i_stat)
      call memocc(i_stat, product(shape(visp))*kind(visp), 'VISP', 'allocate_potential')
      visp = 0.e0_dp

    else
      if (allocated(irc)) then
        i_all = -product(shape(irc))*kind(irc)
        deallocate (irc, stat=i_stat)
        call memocc(i_stat, i_all, 'IRC', 'allocate_potential')
      end if
      if (allocated(ircut)) then
        i_all = -product(shape(ircut))*kind(ircut)
        deallocate (ircut, stat=i_stat)
        call memocc(i_stat, i_all, 'IRCUT', 'allocate_potential')
      end if
      if (allocated(irmin)) then
        i_all = -product(shape(irmin))*kind(irmin)
        deallocate (irmin, stat=i_stat)
        call memocc(i_stat, i_all, 'IRMIN', 'allocate_potential')
      end if
      if (allocated(lcore)) then
        i_all = -product(shape(lcore))*kind(lcore)
        deallocate (lcore, stat=i_stat)
        call memocc(i_stat, i_all, 'LCORE', 'allocate_potential')
      end if
      if (allocated(lmsp1)) then
        i_all = -product(shape(lmsp1))*kind(lmsp1)
        deallocate (lmsp1, stat=i_stat)
        call memocc(i_stat, i_all, 'LMSP1', 'allocate_potential')
      end if
      if (allocated(llmsp)) then
        i_all = -product(shape(llmsp))*kind(llmsp)
        deallocate (llmsp, stat=i_stat)
        call memocc(i_stat, i_all, 'LLMSP', 'allocate_potential')
      end if
      if (allocated(lmsp)) then
        i_all = -product(shape(lmsp))*kind(lmsp)
        deallocate (lmsp, stat=i_stat)
        call memocc(i_stat, i_all, 'LMSP', 'allocate_potential')
      end if
      if (allocated(ititle)) then
        i_all = -product(shape(ititle))*kind(ititle)
        deallocate (ititle, stat=i_stat)
        call memocc(i_stat, i_all, 'ITITLE', 'allocate_potential')
      end if
      if (allocated(nfu)) then
        i_all = -product(shape(nfu))*kind(nfu)
        deallocate (nfu, stat=i_stat)
        call memocc(i_stat, i_all, 'NFU', 'allocate_potential')
      end if
      if (allocated(vins)) then
        i_all = -product(shape(vins))*kind(vins)
        deallocate (vins, stat=i_stat)
        call memocc(i_stat, i_all, 'VINS', 'allocate_potential')
      end if
      if (allocated(visp)) then
        i_all = -product(shape(visp))*kind(visp)
        deallocate (visp, stat=i_stat)
        call memocc(i_stat, i_all, 'VISP', 'allocate_potential')
      end if

    end if

  end subroutine allocate_potential

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the CPA treatment 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, coherent-potential-approx, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe the 
  !> CPA treatment
  !-------------------------------------------------------------------------------  
  subroutine allocate_cpa(flag, naez, natyp, noq, icpa, iqat, hostimp, conc)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: naez   !! number of atoms in unit cell
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell

    integer, dimension (:), allocatable, intent (inout) :: noq !! Number of diff. atom types located
    integer, dimension (:), allocatable, intent (inout) :: icpa !! ICPA = 0/1 site-dependent  CPA flag
    integer, dimension (:), allocatable, intent (inout) :: iqat !! the site on which an atom is located on a given site
    integer, dimension (:), allocatable, intent (inout) :: hostimp
    real (kind=dp), dimension (:), allocatable, intent (inout) :: conc !! concentration of a given atom

    ! .. Local variables
    integer :: i_stat, i_all


    if (flag>0) then

      allocate (noq(naez), stat=i_stat)
      call memocc(i_stat, product(shape(noq))*kind(noq), 'NOQ', 'allocate_cpa')
      noq = 1
      allocate (icpa(naez), stat=i_stat)
      call memocc(i_stat, product(shape(icpa))*kind(icpa), 'ICPA', 'allocate_cpa')
      icpa = 0
      allocate (iqat(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(iqat))*kind(iqat), 'IQAT', 'allocate_cpa')
      iqat = 0
      allocate (conc(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(conc))*kind(conc), 'CONC', 'allocate_cpa')
      conc = 1.e0_dp
      allocate (hostimp(0:natyp), stat=i_stat)
      call memocc(i_stat, product(shape(hostimp))*kind(hostimp), 'HOSTIMP', 'allocate_cpa')
      hostimp = 0

    else

      if (allocated(noq)) then
        i_all = -product(shape(noq))*kind(noq)
        deallocate (noq, stat=i_stat)
        call memocc(i_stat, i_all, 'NOQ', 'allocate_cpa')
      end if
      if (allocated(icpa)) then
        i_all = -product(shape(icpa))*kind(icpa)
        deallocate (icpa, stat=i_stat)
        call memocc(i_stat, i_all, 'ICPA', 'allocate_cpa')
      end if
      if (allocated(iqat)) then
        i_all = -product(shape(iqat))*kind(iqat)
        deallocate (iqat, stat=i_stat)
        call memocc(i_stat, i_all, 'IQAT', 'allocate_cpa')
      end if
      if (allocated(conc)) then
        i_all = -product(shape(conc))*kind(conc)
        deallocate (conc, stat=i_stat)
        call memocc(i_stat, i_all, 'CONC', 'allocate_cpa')
      end if
      if (allocated(hostimp)) then
        i_all = -product(shape(hostimp))*kind(hostimp)
        deallocate (hostimp, stat=i_stat)
        call memocc(i_stat, i_all, 'HOSTIMP', 'allocate_cpa')
      end if

    end if

  end subroutine allocate_cpa

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the LDA+U approach 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe the 
  !> LDA+U approach
  !-------------------------------------------------------------------------------  
  subroutine allocate_ldau(flag, natyp, lopt, ueff, jeff, erefldau)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell

    integer, dimension (:), allocatable, intent (inout) :: lopt !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    real (kind=dp), dimension (:), allocatable, intent (inout) :: ueff !! input U parameter for each atom
    real (kind=dp), dimension (:), allocatable, intent (inout) :: jeff !! input J parameter for each atom
    real (kind=dp), dimension (:), allocatable, intent (inout) :: erefldau !! the energies of the projector's wave functions (REAL)

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then
      allocate (lopt(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(lopt))*kind(lopt), 'LOPT', 'allocate_ldau')
      lopt = -1                    ! not perform lda+u (default)
      allocate (ueff(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(ueff))*kind(ueff), 'UEFF', 'allocate_ldau')
      ueff = 0.e0_dp
      allocate (jeff(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(jeff))*kind(jeff), 'JEFF', 'allocate_ldau')
      jeff = 0.e0_dp
      allocate (erefldau(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(erefldau))*kind(erefldau), 'EREFLDAU', 'allocate_ldau')
      erefldau = 0.5e0_dp
    else
      if (allocated(lopt)) then
        i_all = -product(shape(lopt))*kind(lopt)
        deallocate (lopt, stat=i_stat)
        call memocc(i_stat, i_all, 'LOPT', 'allocate_ldau')
      end if
      if (allocated(ueff)) then
        i_all = -product(shape(ueff))*kind(ueff)
        deallocate (ueff, stat=i_stat)
        call memocc(i_stat, i_all, 'UEFF', 'allocate_ldau')
      end if
      if (allocated(jeff)) then
        i_all = -product(shape(jeff))*kind(jeff)
        deallocate (jeff, stat=i_stat)
        call memocc(i_stat, i_all, 'JEFF', 'allocate_ldau')
      end if
      if (allocated(erefldau)) then
        i_all = -product(shape(erefldau))*kind(erefldau)
        deallocate (erefldau, stat=i_stat)
        call memocc(i_stat, i_all, 'EREFLDAU', 'allocate_ldau')
      end if

    end if

  end subroutine allocate_ldau

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the potentials for the LDA+U approach 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, potential, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe the 
  !> potentials for the LDA+U approach
  !-------------------------------------------------------------------------------  
  subroutine allocate_ldau_potential(flag, irm, natyp, mmaxd, nspind, itldau, wldau,&
     uldau, phildau)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: irm
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    integer, intent (in) :: mmaxd
    integer, intent (in) :: nspind !! Counter for spin directions (KREL+(1-KREL)*(KSP+1))
    integer, dimension (:), allocatable, intent (inout) :: itldau !! integer pointer connecting the NTLDAU atoms to their corresponding index in the unit cell
    real (kind=dp), dimension (:, :, :, :), allocatable, intent (inout) :: wldau !! potential matrix
    real (kind=dp), dimension (:, :, :, :, :), allocatable, intent (inout) :: uldau !! calculated Coulomb matrix elements (EREFLDAU)
    complex (kind=dp), dimension (:, :), allocatable, intent (inout) :: phildau

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then

      allocate (itldau(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(itldau))*kind(itldau), 'ITLDAU', 'allocate_ldau_potential')
      itldau = 0
      allocate (uldau(mmaxd,mmaxd,mmaxd,mmaxd,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(uldau))*kind(uldau), 'ULDAU', 'allocate_ldau_potential')
      uldau = 0.e0_dp
      allocate (wldau(mmaxd,mmaxd,nspind,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(wldau))*kind(wldau), 'WLDAU', 'allocate_ldau_potential')
      wldau = 0.e0_dp
      allocate (phildau(irm,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(phildau))*kind(phildau), 'PHILDAU', 'allocate_ldau_potential')
      phildau = (0.e0_dp, 0.e0_dp)

    else
      if (allocated(itldau)) then
        i_all = -product(shape(itldau))*kind(itldau)
        deallocate (itldau, stat=i_stat)
        call memocc(i_stat, i_all, 'ITLDAU', 'allocate_ldau_potential')
      end if
      if (allocated(uldau)) then
        i_all = -product(shape(uldau))*kind(uldau)
        deallocate (uldau, stat=i_stat)
        call memocc(i_stat, i_all, 'ULDAU', 'allocate_ldau_potential')
      end if
      if (allocated(wldau)) then
        i_all = -product(shape(wldau))*kind(wldau)
        deallocate (wldau, stat=i_stat)
        call memocc(i_stat, i_all, 'WLDAU', 'allocate_ldau_potential')
      end if
      if (allocated(phildau)) then
        i_all = -product(shape(phildau))*kind(phildau)
        deallocate (phildau, stat=i_stat)
        call memocc(i_stat, i_all, 'PHILDAU', 'allocate_ldau_potential')
      end if

    end if

  end subroutine allocate_ldau_potential

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the magnetization 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe the 
  !> magnetization
  !-------------------------------------------------------------------------------  
  subroutine allocate_magnetization(flag, naez, natyp, lmmaxd, inipol, ixipol, qmtet, qmphi, drotq)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: naez   !! number of atoms in unit cell
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    integer, intent (in) :: lmmaxd
    integer, dimension (:), allocatable, intent (inout) :: inipol !! Initial spin polarisation
    integer, dimension (:), allocatable, intent (inout) :: ixipol !! Constraint of spin pol.
    real (kind=dp), dimension (:), allocatable, intent (inout) :: qmtet
    real (kind=dp), dimension (:), allocatable, intent (inout) :: qmphi
    complex (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: drotq !! Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then
      allocate (qmtet(naez), stat=i_stat)
      call memocc(i_stat, product(shape(qmtet))*kind(qmtet), 'QMTET', 'allocate_magnetization')
      qmtet = 0.e0_dp
      allocate (qmphi(naez), stat=i_stat)
      call memocc(i_stat, product(shape(qmphi))*kind(qmphi), 'QMPHI', 'allocate_magnetization')
      qmphi = 0.e0_dp
      allocate (inipol(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(inipol))*kind(inipol), 'INIPOL', 'allocate_magnetization')
      inipol = 0
      allocate (ixipol(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(ixipol))*kind(ixipol), 'IXIPOL', 'allocate_magnetization')
      ixipol = 0
      allocate (drotq(lmmaxd,lmmaxd,naez), stat=i_stat)
      call memocc(i_stat, product(shape(drotq))*kind(drotq), 'DROTQ', 'allocate_magnetization')
      drotq = (0.e0_dp, 0.e0_dp)
    else
      if (allocated(qmtet)) then
        i_all = -product(shape(qmtet))*kind(qmtet)
        deallocate (qmtet, stat=i_stat)
        call memocc(i_stat, i_all, 'QMTET', 'allocate_magnetization')
      end if
      if (allocated(qmphi)) then
        i_all = -product(shape(qmphi))*kind(qmphi)
        deallocate (qmphi, stat=i_stat)
        call memocc(i_stat, i_all, 'QMPHI', 'allocate_magnetization')
      end if
      if (allocated(inipol)) then
        i_all = -product(shape(inipol))*kind(inipol)
        deallocate (inipol, stat=i_stat)
        call memocc(i_stat, i_all, 'INIPOL', 'allocate_magnetization')
      end if
      if (allocated(ixipol)) then
        i_all = -product(shape(ixipol))*kind(ixipol)
        deallocate (ixipol, stat=i_stat)
        call memocc(i_stat, i_all, 'IXIPOL', 'allocate_magnetization')
      end if
      if (allocated(drotq)) then
        i_all = -product(shape(drotq))*kind(drotq)
        deallocate (drotq, stat=i_stat)
        call memocc(i_stat, i_all, 'DROTQ', 'allocate_magnetization')
      end if

    end if

  end subroutine allocate_magnetization

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the spin-orbit coupling (SOC) 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, spin-orbit-coupling, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe the 
  !> spin-orbit coupling (SOC)
  !-------------------------------------------------------------------------------  
  subroutine allocate_soc(flag, krel, natyp, lmax, socscale, cscl, socscl)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: krel
    integer, intent (in) :: lmax   !! Maximum l component in wave function  expansion
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    real (kind=dp), dimension (:), allocatable, intent (inout) :: socscale !! Spin-orbit scaling
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: cscl !! Speed of light scaling
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: socscl 

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then
      allocate (socscl(krel*lmax+1,krel*natyp+(1-krel)), stat=i_stat)
      call memocc(i_stat, product(shape(socscl))*kind(socscl), 'SOCSCL', 'allocate_SOC')
      socscl = 1.e0_dp
      allocate (cscl(krel*lmax+1,krel*natyp+(1-krel)), stat=i_stat)
      call memocc(i_stat, product(shape(cscl))*kind(cscl), 'CSCL', 'allocate_SOC')
      cscl = 0.e0_dp
      allocate (socscale(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(socscale))*kind(socscale), 'SOCSCALE', 'allocate_SOC')
      socscale = 1.e0_dp           ! Spin-orbit scaling
    else
      if (allocated(socscl)) then
        i_all = -product(shape(socscl))*kind(socscl)
        deallocate (socscl, stat=i_stat)
        call memocc(i_stat, i_all, 'SOCSCL', 'allocate_SOC')
      end if
      if (allocated(socscale)) then
        i_all = -product(shape(socscale))*kind(socscale)
        deallocate (socscale, stat=i_stat)
        call memocc(i_stat, i_all, 'SOCSCALE', 'allocate_SOC')
      end if
      if (allocated(cscl)) then
        i_all = -product(shape(cscl))*kind(cscl)
        deallocate (cscl, stat=i_stat)
        call memocc(i_stat, i_all, 'CSCL', 'allocate_SOC')
      end if
    end if

  end subroutine allocate_soc

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe energies 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, total-energy, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe 
  !> energies 
  !-------------------------------------------------------------------------------  
  subroutine allocate_energies(flag, iemxd, ez, dez, wez)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: iemxd

    complex (kind=dp), dimension (:), allocatable, intent (inout) :: ez
    complex (kind=dp), dimension (:), allocatable, intent (inout) :: dez
    complex (kind=dp), dimension (:), allocatable, intent (inout) :: wez

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then
      allocate (ez(iemxd), stat=i_stat)
      call memocc(i_stat, product(shape(ez))*kind(ez), 'EZ', 'allocate_energies')
      ez = (0.e0_dp, 0.e0_dp)
      allocate (dez(iemxd), stat=i_stat)
      call memocc(i_stat, product(shape(dez))*kind(dez), 'DEZ', 'allocate_energies')
      dez = (0.e0_dp, 0.e0_dp)
      allocate (wez(iemxd), stat=i_stat)
      call memocc(i_stat, product(shape(wez))*kind(wez), 'WEZ', 'allocate_energies')
      wez = (0.e0_dp, 0.e0_dp)
    else
      if (allocated(ez)) then
        i_all = -product(shape(ez))*kind(ez)
        deallocate (ez, stat=i_stat)
        call memocc(i_stat, i_all, 'EZ', 'allocate_energies')
      end if
      if (allocated(dez)) then
        i_all = -product(shape(dez))*kind(dez)
        deallocate (dez, stat=i_stat)
        call memocc(i_stat, i_all, 'DEZ', 'allocate_energies')
      end if
      if (allocated(wez)) then
        i_all = -product(shape(wez))*kind(wez)
        deallocate (wez, stat=i_stat)
        call memocc(i_stat, i_all, 'WEZ', 'allocate_energies')
      end if

    end if

  end subroutine allocate_energies

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe relativistic corrections
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe 
  !> relativistic corrections
  !-------------------------------------------------------------------------------  
  subroutine allocate_relativistic(flag,krel,irm,naez,natyp,zrel,jwsrel,irshift,    &
    vtrel,btrel,rmrel,drdirel,r2drdirel,qmgam,qmgamtab,qmphitab,qmtettab)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: krel
    integer, intent (in) :: irm
    integer, intent (in) :: naez   !! number of atoms in unit cell
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    integer, dimension (:), allocatable, intent (inout) :: zrel !! atomic number (cast integer)
    integer, dimension (:), allocatable, intent (inout) :: jwsrel !! index of the WS radius
    integer, dimension (:), allocatable, intent (inout) :: irshift !! shift of the REL radial mesh with respect no NREL
    real (kind=dp), dimension (:), allocatable, intent (inout) :: qmgam
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: vtrel  !! potential (spherical part)
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: btrel  !! magnetic field
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: rmrel  !! radial mesh
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: drdirel !! derivative of radial mesh
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: r2drdirel !! r**2 * drdi
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: qmgamtab
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: qmphitab
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: qmtettab

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then
      allocate (vtrel(irm*krel+(1-krel),natyp), stat=i_stat)
      call memocc(i_stat, product(shape(vtrel))*kind(vtrel), 'VTREL', 'allocate_relativistic')
      vtrel = 0.e0_dp
      allocate (btrel(irm*krel+(1-krel),natyp), stat=i_stat)
      call memocc(i_stat, product(shape(btrel))*kind(btrel), 'BTREL', 'allocate_relativistic')
      btrel = 0.e0_dp
      allocate (drdirel(irm*krel+(1-krel),natyp), stat=i_stat)
      call memocc(i_stat, product(shape(drdirel))*kind(drdirel), 'DRDIREL', 'allocate_relativistic')
      drdirel = 0.e0_dp
      allocate (r2drdirel(irm*krel+(1-krel),natyp), stat=i_stat)
      call memocc(i_stat, product(shape(r2drdirel))*kind(r2drdirel), 'R2DRDIREL', 'allocate_relativistic')
      r2drdirel = 0.e0_dp
      allocate (rmrel(irm*krel+(1-krel),natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rmrel))*kind(rmrel), 'RMREL', 'allocate_relativistic')
      rmrel = 0.e0_dp
      allocate (irshift(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(irshift))*kind(irshift), 'IRSHIFT', 'allocate_relativistic')
      irshift = 0
      allocate (jwsrel(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(jwsrel))*kind(jwsrel), 'JWSREL', 'allocate_relativistic')
      jwsrel = 0
      allocate (zrel(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(zrel))*kind(zrel), 'ZREL', 'allocate_relativistic')
      zrel = 0
      allocate (qmgam(naez), stat=i_stat)
      call memocc(i_stat, product(shape(qmgam))*kind(qmgam), 'QMGAM', 'allocate_relativistic')
      qmgam = 0.e0_dp
      allocate (qmgamtab(naez,3), stat=i_stat)
      call memocc(i_stat, product(shape(qmgamtab))*kind(qmgamtab), 'QMGAMTAB', 'allocate_relativistic')
      qmgamtab = 0.e0_dp
      allocate (qmphitab(naez,3), stat=i_stat)
      call memocc(i_stat, product(shape(qmphitab))*kind(qmphitab), 'QMPHITAB', 'allocate_relativistic')
      qmphitab = 0.e0_dp
      allocate (qmtettab(naez,3), stat=i_stat)
      call memocc(i_stat, product(shape(qmtettab))*kind(qmtettab), 'QMTETTAB', 'allocate_relativistic')
      qmtettab = 0.e0_dp

    else
      if (allocated(vtrel)) then
        i_all = -product(shape(vtrel))*kind(vtrel)
        deallocate (vtrel, stat=i_stat)
        call memocc(i_stat, i_all, 'VTREL', 'allocate_relativistic')
      end if
      if (allocated(btrel)) then
        i_all = -product(shape(btrel))*kind(btrel)
        deallocate (btrel, stat=i_stat)
        call memocc(i_stat, i_all, 'BTREL', 'allocate_relativistic')
      end if
      if (allocated(drdirel)) then
        i_all = -product(shape(drdirel))*kind(drdirel)
        deallocate (drdirel, stat=i_stat)
        call memocc(i_stat, i_all, 'DRDIREL', 'allocate_relativistic')
      end if
      if (allocated(r2drdirel)) then
        i_all = -product(shape(r2drdirel))*kind(r2drdirel)
        deallocate (r2drdirel, stat=i_stat)
        call memocc(i_stat, i_all, 'R2DRDIREL', 'allocate_relativistic')
      end if
      if (allocated(rmrel)) then
        i_all = -product(shape(rmrel))*kind(rmrel)
        deallocate (rmrel, stat=i_stat)
        call memocc(i_stat, i_all, 'RMREL', 'allocate_relativistic')
      end if
      if (allocated(irshift)) then
        i_all = -product(shape(irshift))*kind(irshift)
        deallocate (irshift, stat=i_stat)
        call memocc(i_stat, i_all, 'IRSHIFT', 'allocate_relativistic')
      end if
      if (allocated(jwsrel)) then
        i_all = -product(shape(jwsrel))*kind(jwsrel)
        deallocate (jwsrel, stat=i_stat)
        call memocc(i_stat, i_all, 'JWSREL', 'allocate_relativistic')
      end if
      if (allocated(zrel)) then
        i_all = -product(shape(zrel))*kind(zrel)
        deallocate (zrel, stat=i_stat)
        call memocc(i_stat, i_all, 'ZREL', 'allocate_relativistic')
      end if
      if (allocated(qmgam)) then
        i_all = -product(shape(qmgam))*kind(qmgam)
        deallocate (qmgam, stat=i_stat)
        call memocc(i_stat, i_all, 'QMGAM', 'allocate_relativistic')
      end if
      if (allocated(qmgamtab)) then
        i_all = -product(shape(qmgamtab))*kind(qmgamtab)
        deallocate (qmgamtab, stat=i_stat)
        call memocc(i_stat, i_all, 'QMGAMTAB', 'allocate_relativistic')
      end if
      if (allocated(qmphitab)) then
        i_all = -product(shape(qmphitab))*kind(qmphitab)
        deallocate (qmphitab, stat=i_stat)
        call memocc(i_stat, i_all, 'QMPHITAB', 'allocate_relativistic')
      end if
      if (allocated(qmtettab)) then
        i_all = -product(shape(qmtettab))*kind(qmtettab)
        deallocate (qmtettab, stat=i_stat)
        call memocc(i_stat, i_all, 'QMTETTAB', 'allocate_relativistic')
      end if

    end if

  end subroutine allocate_relativistic

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe relativistic transformations 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, dirac, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe 
  !> relativistic transformations
  !-------------------------------------------------------------------------------  
  subroutine allocate_rel_transformations(flag,lmmaxd,nrrel,irrel,rc,crel,rrel,srrel)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: lmmaxd
    integer, dimension (:, :), allocatable, intent (inout) :: nrrel
    integer, dimension (:, :, :), allocatable, intent (inout) :: irrel
    complex (kind=dp), dimension (:, :), allocatable, intent (inout) :: rc !! NREL REAL spher. harm. >  CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
    complex (kind=dp), dimension (:, :), allocatable, intent (inout) :: crel !! Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
    complex (kind=dp), dimension (:, :), allocatable, intent (inout) :: rrel  !! Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
    complex (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: srrel

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then

      allocate (rrel(lmmaxd,lmmaxd), stat=i_stat)
      call memocc(i_stat, product(shape(rrel))*kind(rrel), 'RREL', 'allocate_rel_transformations')
      rrel = (0.e0_dp, 0.e0_dp)
      allocate (srrel(2,2,lmmaxd), stat=i_stat)
      call memocc(i_stat, product(shape(srrel))*kind(srrel), 'SRREL', 'allocate_rel_transformations')
      srrel = (0.e0_dp, 0.e0_dp)
      allocate (irrel(2,2,lmmaxd), stat=i_stat)
      call memocc(i_stat, product(shape(irrel))*kind(irrel), 'IRREL', 'allocate_rel_transformations')
      irrel = 0
      allocate (nrrel(2,lmmaxd), stat=i_stat)
      call memocc(i_stat, product(shape(nrrel))*kind(nrrel), 'NRREL', 'allocate_rel_transformations')
      nrrel = 0
      allocate (crel(lmmaxd,lmmaxd), stat=i_stat)
      call memocc(i_stat, product(shape(crel))*kind(crel), 'CREL', 'allocate_rel_transformations')
      crel = (0.e0_dp, 0.e0_dp)
      allocate (rc(lmmaxd,lmmaxd), stat=i_stat)
      call memocc(i_stat, product(shape(rc))*kind(rc), 'RC', 'allocate_rel_transformations')
      rc = (0.e0_dp, 0.e0_dp)
    else
      if (allocated(rrel)) then
        i_all = -product(shape(rrel))*kind(rrel)
        deallocate (rrel, stat=i_stat)
        call memocc(i_stat, i_all, 'RREL', 'allocate_rel_transformations')
      end if
      if (allocated(srrel)) then
        i_all = -product(shape(srrel))*kind(srrel)
        deallocate (srrel, stat=i_stat)
        call memocc(i_stat, i_all, 'SRREL', 'allocate_rel_transformations')
      end if
      if (allocated(irrel)) then
        i_all = -product(shape(irrel))*kind(irrel)
        deallocate (irrel, stat=i_stat)
        call memocc(i_stat, i_all, 'IRREL', 'allocate_rel_transformations')
      end if
      if (allocated(nrrel)) then
        i_all = -product(shape(nrrel))*kind(nrrel)
        deallocate (nrrel, stat=i_stat)
        call memocc(i_stat, i_all, 'NRREL', 'allocate_rel_transformations')
      end if
      if (allocated(crel)) then
        i_all = -product(shape(crel))*kind(crel)
        deallocate (crel, stat=i_stat)
        call memocc(i_stat, i_all, 'CREL', 'allocate_rel_transformations')
      end if
      if (allocated(rc)) then
        i_all = -product(shape(rc))*kind(rc)
        deallocate (rc, stat=i_stat)
        call memocc(i_stat, i_all, 'RC', 'allocate_rel_transformations')
      end if

    end if

  end subroutine allocate_rel_transformations

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe clusters 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, geometry, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe 
  !> clusters
  !-------------------------------------------------------------------------------  
  subroutine allocate_clusters(flag,naez,lmax,ncleb,nclsd,nembd1,nsheld,naclsd,     &
    lmpot,natomimpd,nsh1,nsh2,nacls,nshell,atomimp,atom,ezoa,icleb,jend,ratom,      &
    rclsimp,cmomhost,rcls)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: naez   !! number of atoms in unit cell
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: ncleb
    integer, intent (in) :: nclsd
    integer, intent (in) :: nembd1
    integer, intent (in) :: nsheld
    integer, intent (in) :: naclsd
    integer, intent (in) :: lmpot
    integer, intent (in) :: natomimpd
    integer, dimension (:), allocatable, intent (inout) :: nsh1 !! Corresponding index of the sites I/J in (NSH1/2) in the unit cell in a shell
    integer, dimension (:), allocatable, intent (inout) :: nsh2 !! Corresponding index of the sites I/J in (NSH1/2) in the unit cell in a shell
    integer, dimension (:), allocatable, intent (inout) :: nacls !! Number of atoms in cluster
    integer, dimension (:), allocatable, intent (inout) :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
    integer, dimension (:), allocatable, intent (inout) :: atomimp
    integer, dimension (:, :), allocatable, intent (inout) :: atom !! Atom at site in cluster
    integer, dimension (:, :), allocatable, intent (inout) :: ezoa !! EZ of atom at site in cluster
    integer, dimension (:, :), allocatable, intent (inout) :: icleb !! Pointer array
    integer, dimension (:, :, :), allocatable, intent (inout) :: jend !! Pointer array for icleb()
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: ratom
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: rclsimp
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: cmomhost !! Charge moments of each atom of the (left/right) host
    real (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: rcls  !! Real space position of atom in cluster

    integer :: i_stat, i_all

    if (flag>0) then

      allocate (atom(naclsd,naez+(nembd1-1)), stat=i_stat)
      call memocc(i_stat, product(shape(atom))*kind(atom), 'ATOM', 'allocate_clusters')
      atom = 0
      allocate (ratom(3,nsheld), stat=i_stat)
      call memocc(i_stat, product(shape(ratom))*kind(ratom), 'RATOM', 'allocate_clusters')
      ratom = 0.e0_dp
      allocate (rcls(3,naclsd,nclsd), stat=i_stat)
      call memocc(i_stat, product(shape(rcls))*kind(rcls), 'RCLS', 'allocate_clusters')
      rcls = 0.e0_dp
      allocate (rclsimp(3,natomimpd), stat=i_stat)
      call memocc(i_stat, product(shape(rclsimp))*kind(rclsimp), 'RCLSIMP', 'allocate_clusters')
      rclsimp = 0.e0_dp
      allocate (nacls(nclsd), stat=i_stat)
      call memocc(i_stat, product(shape(nacls))*kind(nacls), 'NACLS', 'allocate_clusters')
      nacls = 0
      allocate (ezoa(naclsd,naez+(nembd1-1)), stat=i_stat)
      call memocc(i_stat, product(shape(ezoa))*kind(ezoa), 'EZOA', 'allocate_clusters')
      ezoa = 0
      allocate (atomimp(natomimpd), stat=i_stat)
      call memocc(i_stat, product(shape(atomimp))*kind(atomimp), 'ATOMIMP', 'allocate_clusters')
      atomimp = 0
      allocate (icleb(ncleb,4), stat=i_stat)
      call memocc(i_stat, product(shape(icleb))*kind(icleb), 'ICLEB', 'allocate_clusters')
      icleb = 0
      allocate (nsh1(nsheld), stat=i_stat)
      call memocc(i_stat, product(shape(nsh1))*kind(nsh1), 'NSH1', 'allocate_clusters')
      nsh1 = 0
      allocate (nsh2(nsheld), stat=i_stat)
      call memocc(i_stat, product(shape(nsh2))*kind(nsh2), 'NSH2', 'allocate_clusters')
      nsh2 = 0
      allocate (nshell(0:nsheld), stat=i_stat)
      call memocc(i_stat, product(shape(nshell))*kind(nshell), 'NSHELL', 'allocate_clusters')
      nshell = 0
      allocate (cmomhost(lmpot,nembd1), stat=i_stat)
      call memocc(i_stat, product(shape(cmomhost))*kind(cmomhost), 'CMOMHOST', 'allocate_clusters')
      cmomhost = 0.e0_dp
      allocate (jend(lmpot,0:lmax,0:lmax), stat=i_stat)
      call memocc(i_stat, product(shape(jend))*kind(jend), 'JEND', 'allocate_clusters')
      jend = 0

    else
      if (allocated(atom)) then
        i_all = -product(shape(atom))*kind(atom)
        deallocate (atom, stat=i_stat)
        call memocc(i_stat, i_all, 'ATOM', 'allocate_clusters')
      end if
      if (allocated(ratom)) then
        i_all = -product(shape(ratom))*kind(ratom)
        deallocate (ratom, stat=i_stat)
        call memocc(i_stat, i_all, 'RATOM', 'allocate_clusters')
      end if
      if (allocated(rcls)) then
        i_all = -product(shape(rcls))*kind(rcls)
        deallocate (rcls, stat=i_stat)
        call memocc(i_stat, i_all, 'RCLS', 'allocate_clusters')
      end if
      if (allocated(rclsimp)) then
        i_all = -product(shape(rclsimp))*kind(rclsimp)
        deallocate (rclsimp, stat=i_stat)
        call memocc(i_stat, i_all, 'RCLSIMP', 'allocate_clusters')
      end if
      if (allocated(nacls)) then
        i_all = -product(shape(nacls))*kind(nacls)
        deallocate (nacls, stat=i_stat)
        call memocc(i_stat, i_all, 'NACLS', 'allocate_clusters')
      end if
      if (allocated(ezoa)) then
        i_all = -product(shape(ezoa))*kind(ezoa)
        deallocate (ezoa, stat=i_stat)
        call memocc(i_stat, i_all, 'EZOA', 'allocate_clusters')
      end if
      if (allocated(atomimp)) then
        i_all = -product(shape(atomimp))*kind(atomimp)
        deallocate (atomimp, stat=i_stat)
        call memocc(i_stat, i_all, 'ATOMIMP', 'allocate_clusters')
      end if
      if (allocated(icleb)) then
        i_all = -product(shape(icleb))*kind(icleb)
        deallocate (icleb, stat=i_stat)
        call memocc(i_stat, i_all, 'ICLEB', 'allocate_clusters')
      end if
      if (allocated(nsh1)) then
        i_all = -product(shape(nsh1))*kind(nsh1)
        deallocate (nsh1, stat=i_stat)
        call memocc(i_stat, i_all, 'NSH1', 'allocate_clusters')
      end if
      if (allocated(nsh2)) then
        i_all = -product(shape(nsh2))*kind(nsh2)
        deallocate (nsh2, stat=i_stat)
        call memocc(i_stat, i_all, 'NSH2', 'allocate_clusters')
      end if
      if (allocated(nshell)) then
        i_all = -product(shape(nshell))*kind(nshell)
        deallocate (nshell, stat=i_stat)
        call memocc(i_stat, i_all, 'NSHELL', 'allocate_clusters')
      end if
      if (allocated(cmomhost)) then
        i_all = -product(shape(cmomhost))*kind(cmomhost)
        deallocate (cmomhost, stat=i_stat)
        call memocc(i_stat, i_all, 'CMOMHOST', 'allocate_clusters')
      end if
      if (allocated(jend)) then
        i_all = -product(shape(jend))*kind(jend)
        deallocate (jend, stat=i_stat)
        call memocc(i_stat, i_all, 'JEND', 'allocate_clusters')
      end if

    end if

  end subroutine allocate_clusters

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the functions for the expansion of the Green function 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, special-functions, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe 
  !> the functions for the expansion of the Green function
  !-------------------------------------------------------------------------------  
  subroutine allocate_expansion(flag,lm2d,irid,nfund,ntotd,ncleb,lassld,ncelld,     &
    nchebd,loflm,wg,cleb,yrg,thetas,thetasnew)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: lm2d
    integer, intent (in) :: irid
    integer, intent (in) :: nfund
    integer, intent (in) :: ntotd
    integer, intent (in) :: ncleb
    integer, intent (in) :: lassld
    integer, intent (in) :: ncelld
    integer, intent (in) :: nchebd
    integer, dimension (:), allocatable, intent (inout) :: loflm !! l of lm=(l,m) (GAUNT)
    real (kind=dp), dimension (:), allocatable, intent (inout) :: wg !! Integr. weights for Legendre polynomials
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: yrg  !! Spherical harmonics (GAUNT2)
    real (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: thetas !! shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    real (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: thetasnew

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then

      allocate (wg(lassld), stat=i_stat)
      call memocc(i_stat, product(shape(wg))*kind(wg), 'WG', 'allocate_expansion')
      wg = 0.e0_dp
      allocate (yrg(lassld,0:lassld,0:lassld), stat=i_stat)
      call memocc(i_stat, product(shape(yrg))*kind(yrg), 'YRG', 'allocate_expansion')
      yrg = 0.e0_dp
      allocate (thetas(irid,nfund,ncelld), stat=i_stat)
      call memocc(i_stat, product(shape(thetas))*kind(thetas), 'THETAS', 'allocate_expansion')
      thetas = 0.e0_dp
      allocate (thetasnew(ntotd*(nchebd+1),nfund,ncelld), stat=i_stat)
      call memocc(i_stat, product(shape(thetasnew))*kind(thetasnew), 'THETASNEW', 'allocate_expansion')
      thetasnew = 0.e0_dp
      allocate (cleb(ncleb,2), stat=i_stat)
      call memocc(i_stat, product(shape(cleb))*kind(cleb), 'CLEB', 'allocate_expansion')
      cleb = 0.e0_dp
      allocate (loflm(lm2d), stat=i_stat)
      call memocc(i_stat, product(shape(loflm))*kind(loflm), 'LOFLM', 'allocate_expansion')
      loflm = 0

    else
      if (allocated(wg)) then
        i_all = -product(shape(wg))*kind(wg)
        deallocate (wg, stat=i_stat)
        call memocc(i_stat, i_all, 'WG', 'allocate_expansion')
      end if
      if (allocated(yrg)) then
        i_all = -product(shape(yrg))*kind(yrg)
        deallocate (yrg, stat=i_stat)
        call memocc(i_stat, i_all, 'YRG', 'allocate_expansion')
      end if
      if (allocated(thetas)) then
        i_all = -product(shape(thetas))*kind(thetas)
        deallocate (thetas, stat=i_stat)
        call memocc(i_stat, i_all, 'THETAS', 'allocate_expansion')
      end if
      if (allocated(thetasnew)) then
        i_all = -product(shape(thetasnew))*kind(thetasnew)
        deallocate (thetasnew, stat=i_stat)
        call memocc(i_stat, i_all, 'THETASNEW', 'allocate_expansion')
      end if
      if (allocated(cleb)) then
        i_all = -product(shape(cleb))*kind(cleb)
        deallocate (cleb, stat=i_stat)
        call memocc(i_stat, i_all, 'CLEB', 'allocate_expansion')
      end if
      if (allocated(loflm)) then
        i_all = -product(shape(loflm))*kind(loflm)
        deallocate (loflm, stat=i_stat)
        call memocc(i_stat, i_all, 'LOFLM', 'allocate_expansion')
      end if

    end if

  end subroutine allocate_expansion

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the integration mesh 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, radial-grid, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe 
  !> the integration mesh
  !-------------------------------------------------------------------------------  
  subroutine allocate_mesh(flag, irm, natyp, a, b, rmesh, drdi)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: irm
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    real (kind=dp), dimension (:), allocatable, intent (inout) :: a !! Constants for exponential R mesh
    real (kind=dp), dimension (:), allocatable, intent (inout) :: b
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: rmesh !! Radial mesh ( in units a Bohr)
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: drdi !! Derivative dr/di
    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then

      allocate (drdi(irm,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(drdi))*kind(drdi), 'DRDI', 'allocate_mesh')
      drdi = 0.e0_dp
      allocate (rmesh(irm,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rmesh))*kind(rmesh), 'RMESH', 'allocate_mesh')
      rmesh = 0.e0_dp
      allocate (a(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(a))*kind(a), 'A', 'allocate_mesh')
      a = 0.e0_dp
      allocate (b(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(b))*kind(b), 'B', 'allocate_mesh')
      b = 0.e0_dp

    else
      if (allocated(drdi)) then
        i_all = -product(shape(drdi))*kind(drdi)
        deallocate (drdi, stat=i_stat)
        call memocc(i_stat, i_all, 'DRDI', 'allocate_mesh')
      end if
      if (allocated(rmesh)) then
        i_all = -product(shape(rmesh))*kind(rmesh)
        deallocate (rmesh, stat=i_stat)
        call memocc(i_stat, i_all, 'RMESH', 'allocate_mesh')
      end if
      if (allocated(a)) then
        i_all = -product(shape(a))*kind(a)
        deallocate (a, stat=i_stat)
        call memocc(i_stat, i_all, 'A', 'allocate_mesh')
      end if
      if (allocated(b)) then
        i_all = -product(shape(b))*kind(b)
        deallocate (b, stat=i_stat)
        call memocc(i_stat, i_all, 'B', 'allocate_mesh')
      end if

    end if

  end subroutine allocate_mesh

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays that describe the panels 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, radial-grid, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays that describe 
  !> the panels
  !-------------------------------------------------------------------------------  
  subroutine allocate_pannels(flag,natyp,ntotd,ipan,npan_tot,npan_eq_at,npan_log_at,&
    ipan_intervall, rpan_intervall)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    integer, intent (in) :: ntotd  !! IPAND+30
    integer, dimension (:), allocatable, intent (inout) :: ipan !! Number of panels in non-MT-region
    integer, dimension (:), allocatable, intent (inout) :: npan_tot
    integer, dimension (:), allocatable, intent (inout) :: npan_eq_at
    integer, dimension (:), allocatable, intent (inout) :: npan_log_at
    integer, dimension (:, :), allocatable, intent (inout) :: ipan_intervall
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: rpan_intervall

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then
      allocate (ipan(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(ipan))*kind(ipan), 'IPAN', 'allocate_pannels')
      ipan = 0
      allocate (npan_tot(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(npan_tot))*kind(npan_tot), 'NPAN_TOT', 'allocate_pannels')
      npan_tot = 0
      allocate (npan_eq_at(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(npan_eq_at))*kind(npan_eq_at), 'NPAN_EQ_AT', 'allocate_pannels')
      npan_eq_at = 0
      allocate (npan_log_at(natyp), stat=i_stat)
      call memocc(i_stat, product(shape(npan_log_at))*kind(npan_log_at), 'NPAN_LOG_AT', 'allocate_pannels')
      npan_log_at = 0
      allocate (rpan_intervall(0:ntotd,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rpan_intervall))*kind(rpan_intervall), 'RPAN_INTERVALL', 'allocate_pannels')
      rpan_intervall = 0.e0_dp
      allocate (ipan_intervall(0:ntotd,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(ipan_intervall))*kind(ipan_intervall), 'IPAN_INTERVALL', 'allocate_pannels')
      ipan_intervall = 0

    else
      if (allocated(ipan)) then
        i_all = -product(shape(ipan))*kind(ipan)
        deallocate (ipan, stat=i_stat)
        call memocc(i_stat, i_all, 'IPAN', 'allocate_pannels')
      end if
      if (allocated(npan_tot)) then
        i_all = -product(shape(npan_tot))*kind(npan_tot)
        deallocate (npan_tot, stat=i_stat)
        call memocc(i_stat, i_all, 'NPAN_TOT', 'allocate_pannels')
      end if
      if (allocated(npan_eq_at)) then
        i_all = -product(shape(npan_eq_at))*kind(npan_eq_at)
        deallocate (npan_eq_at, stat=i_stat)
        call memocc(i_stat, i_all, 'NPAN_EQ_AT', 'allocate_pannels')
      end if
      if (allocated(npan_log_at)) then
        i_all = -product(shape(npan_log_at))*kind(npan_log_at)
        deallocate (npan_log_at, stat=i_stat)
        call memocc(i_stat, i_all, 'NPAN_LOG_AT', 'allocate_pannels')
      end if
      if (allocated(rpan_intervall)) then
        i_all = -product(shape(rpan_intervall))*kind(rpan_intervall)
        deallocate (rpan_intervall, stat=i_stat)
        call memocc(i_stat, i_all, 'RPAN_INTERVALL', 'allocate_pannels')
      end if
      if (allocated(ipan_intervall)) then
        i_all = -product(shape(ipan_intervall))*kind(ipan_intervall)
        deallocate (ipan_intervall, stat=i_stat)
        call memocc(i_stat, i_all, 'IPAN_INTERVALL', 'allocate_pannels')
      end if

    end if

  end subroutine allocate_pannels

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of misc arrays 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of misc arrays 
  !-------------------------------------------------------------------------------  
  subroutine allocate_misc(flag,nr,irm,irid,lmax,naez,natyp,nfund,nrefd,iemxd,ntotd,&
    nsheld,lmmaxd,nembd1,nchebd,ncelld,lmxspd,nspindd,nsymaxd,nprincd,ifunm,ifunm1, &
    icheck,vref,s,rr,dror,rnew,rs,rrot,thesme,dsymll,dsymll1,lefttinvll,righttinvll)

    implicit none

    integer, intent (in) :: nr
    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: irm
    integer, intent (in) :: irid
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: naez   !! number of atoms in unit cell
    integer, intent (in) :: natyp  !! number of kinds of atoms in unit cell
    integer, intent (in) :: nfund
    integer, intent (in) :: nrefd
    integer, intent (in) :: iemxd
    integer, intent (in) :: ntotd
    integer, intent (in) :: lmmaxd
    integer, intent (in) :: nsheld
    integer, intent (in) :: nembd1
    integer, intent (in) :: nchebd
    integer, intent (in) :: ncelld
    integer, intent (in) :: lmxspd
    integer, intent (in) :: nspindd
    integer, intent (in) :: nsymaxd
    integer, intent (in) :: nprincd
    integer, dimension (:, :), allocatable, intent (inout) :: ifunm
    integer, dimension (:, :), allocatable, intent (inout) :: ifunm1
    integer, dimension (:, :), allocatable, intent (inout) :: icheck
    real (kind=dp), dimension (:), allocatable, intent (inout) :: vref
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: s
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: rr !! Set of real space vectors (in a.u.)
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: dror
    real (kind=dp), dimension (:, :), allocatable, intent (inout) :: rnew
    real (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: rs
    real (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: rrot
    real (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: thesme
    complex (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: dsymll
    complex (kind=dp), dimension (:, :, :), allocatable, intent (inout) :: dsymll1
    complex (kind=dp), dimension (:, :, :, :, :), allocatable, intent (inout) :: lefttinvll
    complex (kind=dp), dimension (:, :, :, :, :), allocatable, intent (inout) :: righttinvll

    ! .. Local variables
    integer :: i_stat, i_all

    if (flag>0) then
      allocate (rr(3,0:nr), stat=i_stat)
      call memocc(i_stat, product(shape(rr))*kind(rr), 'RR', 'allocate_misc')
      rr = 0.e0_dp
      allocate (vref(nrefd), stat=i_stat)
      call memocc(i_stat, product(shape(vref))*kind(vref), 'VREF', 'allocate_misc')
      vref = 0.e0_dp
      allocate (dror(irm,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(dror))*kind(dror), 'DROR', 'allocate_misc')
      dror = 0.e0_dp
      allocate (s(0:lmax,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(s))*kind(s), 'S', 'allocate_misc')
      s = 0.e0_dp
      allocate (rrot(48,3,nsheld), stat=i_stat)
      call memocc(i_stat, product(shape(rrot))*kind(rrot), 'RROT', 'allocate_misc')
      rrot = 0.e0_dp
      allocate (ifunm(natyp,lmxspd), stat=i_stat)
      call memocc(i_stat, product(shape(ifunm))*kind(ifunm), 'IFUNM', 'allocate_misc')
      ifunm = 0
      allocate (ifunm1(lmxspd,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(ifunm1))*kind(ifunm1), 'IFUNM1', 'allocate_misc')
      ifunm1 = 0
      allocate (rs(irm,0:lmax,natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rs))*kind(rs), 'RS', 'allocate_misc')
      rs = 0.e0_dp
      allocate (thesme(irid,nfund,ncelld), stat=i_stat)
      call memocc(i_stat, product(shape(thesme))*kind(thesme), 'THESME', 'allocate_misc')
      thesme = 0.e0_dp
      allocate (rnew(ntotd*(nchebd+1),natyp), stat=i_stat)
      call memocc(i_stat, product(shape(rnew))*kind(rnew), 'RNEW', 'allocate_misc')
      rnew = 0.e0_dp
      allocate (dsymll(lmmaxd,lmmaxd,nsymaxd), stat=i_stat)
      call memocc(i_stat, product(shape(dsymll))*kind(dsymll), 'DSYMLL', 'allocate_misc')
      dsymll = (0.e0_dp, 0.e0_dp)
      allocate (dsymll1(lmmaxd,lmmaxd,nsymaxd), stat=i_stat)
      call memocc(i_stat, product(shape(dsymll1))*kind(dsymll1), 'DSYMLL1', 'allocate_misc')
      dsymll1 = (0.e0_dp, 0.e0_dp)
      allocate (icheck(naez/nprincd,naez/nprincd), stat=i_stat)
      call memocc(i_stat, product(shape(icheck))*kind(icheck), 'ICHECK', 'allocate_misc')
      icheck = 0
      allocate (lefttinvll(lmmaxd,lmmaxd,nembd1,nspindd,iemxd), stat=i_stat)
      call memocc(i_stat, product(shape(lefttinvll))*kind(lefttinvll), 'LEFTTINVLL', 'allocate_misc')
      lefttinvll = (0.e0_dp, 0.e0_dp)
      allocate (righttinvll(lmmaxd,lmmaxd,nembd1,nspindd,iemxd), stat=i_stat)
      call memocc(i_stat, product(shape(righttinvll))*kind(righttinvll), 'RIGHTTINVLL', 'allocate_misc')
      righttinvll = (0.e0_dp, 0.e0_dp)

    else
      if (allocated(rr)) then
        i_all = -product(shape(rr))*kind(rr)
        deallocate (rr, stat=i_stat)
        call memocc(i_stat, i_all, 'RR', 'allocate_misc')
      end if
      if (allocated(vref)) then
        i_all = -product(shape(vref))*kind(vref)
        deallocate (vref, stat=i_stat)
        call memocc(i_stat, i_all, 'VREF', 'allocate_misc')
      end if
      if (allocated(dror)) then
        i_all = -product(shape(dror))*kind(dror)
        deallocate (dror, stat=i_stat)
        call memocc(i_stat, i_all, 'DROR', 'allocate_misc')
      end if
      if (allocated(s)) then
        i_all = -product(shape(s))*kind(s)
        deallocate (s, stat=i_stat)
        call memocc(i_stat, i_all, 'S', 'allocate_misc')
      end if
      if (allocated(rrot)) then
        i_all = -product(shape(rrot))*kind(rrot)
        deallocate (rrot, stat=i_stat)
        call memocc(i_stat, i_all, 'RROT', 'allocate_misc')
      end if
      if (allocated(ifunm)) then
        i_all = -product(shape(ifunm))*kind(ifunm)
        deallocate (ifunm, stat=i_stat)
        call memocc(i_stat, i_all, 'IFUNM', 'allocate_misc')
      end if
      if (allocated(ifunm1)) then
        i_all = -product(shape(ifunm1))*kind(ifunm1)
        deallocate (ifunm1, stat=i_stat)
        call memocc(i_stat, i_all, 'IFUNM1', 'allocate_misc')
      end if
      if (allocated(rs)) then
        i_all = -product(shape(rs))*kind(rs)
        deallocate (rs, stat=i_stat)
        call memocc(i_stat, i_all, 'RS', 'allocate_misc')
      end if
      if (allocated(thesme)) then
        i_all = -product(shape(thesme))*kind(thesme)
        deallocate (thesme, stat=i_stat)
        call memocc(i_stat, i_all, 'THESME', 'allocate_misc')
      end if
      if (allocated(rnew)) then
        i_all = -product(shape(rnew))*kind(rnew)
        deallocate (rnew, stat=i_stat)
        call memocc(i_stat, i_all, 'RNEW', 'allocate_misc')
      end if
      if (allocated(dsymll)) then
        i_all = -product(shape(dsymll))*kind(dsymll)
        deallocate (dsymll, stat=i_stat)
        call memocc(i_stat, i_all, 'DSYMLL', 'allocate_misc')
      end if
      if (allocated(dsymll1)) then
        i_all = -product(shape(dsymll1))*kind(dsymll1)
        deallocate (dsymll1, stat=i_stat)
        call memocc(i_stat, i_all, 'DSYMLL1', 'allocate_misc')
      end if
      if (allocated(icheck)) then
        i_all = -product(shape(icheck))*kind(icheck)
        deallocate (icheck, stat=i_stat)
        call memocc(i_stat, i_all, 'ICHECK', 'allocate_misc')
      end if
      if (allocated(lefttinvll)) then
        i_all = -product(shape(lefttinvll))*kind(lefttinvll)
        deallocate (lefttinvll, stat=i_stat)
        call memocc(i_stat, i_all, 'LEFTTINVLL', 'allocate_misc')
      end if
      if (allocated(righttinvll)) then
        i_all = -product(shape(righttinvll))*kind(righttinvll)
        deallocate (righttinvll, stat=i_stat)
        call memocc(i_stat, i_all, 'RIGHTTINVLL', 'allocate_misc')
      end if

    end if

  end subroutine allocate_misc

  !-------------------------------------------------------------------------------  
  !> Summary: subroutine handling the allocation/deallocation of arrays handling the Green functions 
  !> Author: Jonathan Chico 
  !> Category: profiling, profiling, KKRhost 
  !> Deprecated: False 
  !> Subroutine handling the allocation/deallocation of arrays handling the
  !> Green functions
  !-------------------------------------------------------------------------------  
  subroutine allocate_green(flag,naez,iemxd,ngshd,nsheld,lmpot,nofgijd,ish,jsh,     &
    kmesh,imaxsh,iqcalc,iofgij,jofgij,ijtabsh,ijtabsym,ijtabcalc,ijtabcalc_i,       &
    ilm_map,gsh)

    implicit none

    integer, intent (in) :: flag   ! Allocate/deallocate (1/-1) arrays
    integer, intent (in) :: naez   !! number of atoms in unit cell
    integer, intent (in) :: iemxd
    integer, intent (in) :: ngshd
    integer, intent (in) :: nsheld
    integer, intent (in) :: lmpot
    integer, intent (in) :: nofgijd
    integer, dimension (:,:), allocatable, intent (inout) :: ish
    integer, dimension (:,:), allocatable, intent (inout) :: jsh
    integer, dimension (:), allocatable, intent (inout) :: kmesh
    integer, dimension (:), allocatable, intent (inout) :: imaxsh
    integer, dimension (:), allocatable, intent (inout) :: iqcalc
    integer, dimension (:), allocatable, intent (inout) :: iofgij !! Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
    integer, dimension (:), allocatable, intent (inout) :: jofgij !! Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
    integer, dimension (:), allocatable, intent (inout) :: ijtabsh !! Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
    integer, dimension (:), allocatable, intent (inout) :: ijtabsym !! Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
    integer, dimension (:), allocatable, intent (inout) :: ijtabcalc !! Linear pointer,specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
    integer, dimension (:), allocatable, intent (inout) :: ijtabcalc_i
    integer, dimension (:,:), allocatable, intent (inout) :: ilm_map
    real (kind=dp), dimension (:), allocatable, intent (inout) :: gsh

    integer :: i_stat, i_all

    if (flag>0) then

      allocate (gsh(ngshd), stat=i_stat)
      call memocc(i_stat, product(shape(gsh))*kind(gsh), 'GSH', 'allocate_green')
      gsh = 0.e0_dp
      allocate (kmesh(iemxd), stat=i_stat)
      call memocc(i_stat, product(shape(kmesh))*kind(kmesh), 'KMESH', 'allocate_green')
      kmesh = 0
      allocate (ilm_map(ngshd,3), stat=i_stat)
      call memocc(i_stat, product(shape(ilm_map))*kind(ilm_map), 'ILM_MAP', 'allocate_green')
      ilm_map = 0
      allocate (iqcalc(naez), stat=i_stat)
      call memocc(i_stat, product(shape(iqcalc))*kind(iqcalc), 'IQCALC', 'allocate_green')
      iqcalc = 0
      allocate (jofgij(nofgijd), stat=i_stat)
      call memocc(i_stat, product(shape(jofgij))*kind(jofgij), 'JOFGIJ', 'allocate_green')
      jofgij = 0
      allocate (iofgij(nofgijd), stat=i_stat)
      call memocc(i_stat, product(shape(iofgij))*kind(iofgij), 'IOFGIJ', 'allocate_green')
      iofgij = 0
      allocate (imaxsh(0:lmpot), stat=i_stat)
      call memocc(i_stat, product(shape(imaxsh))*kind(imaxsh), 'IMAXSH', 'allocate_green')
      imaxsh = 0
      allocate (ijtabsh(nofgijd), stat=i_stat)
      call memocc(i_stat, product(shape(ijtabsh))*kind(ijtabsh), 'IJTABSH', 'allocate_green')
      ijtabsh = 0
      allocate (ijtabsym(nofgijd), stat=i_stat)
      call memocc(i_stat, product(shape(ijtabsym))*kind(ijtabsym), 'IJTABSYM', 'allocate_green')
      ijtabsym = 0
      allocate (ijtabcalc(nofgijd), stat=i_stat)
      call memocc(i_stat, product(shape(ijtabcalc))*kind(ijtabcalc), 'IJTABCALC', 'allocate_green')
      ijtabcalc = 0
      allocate (ish(nsheld,nofgijd), stat=i_stat)
      call memocc(i_stat, product(shape(ish))*kind(ish), 'ISH', 'allocate_green')
      ish = 0
      allocate (jsh(nsheld,nofgijd), stat=i_stat)
      call memocc(i_stat, product(shape(jsh))*kind(jsh), 'JSH', 'allocate_green')
      jsh = 0
      allocate (ijtabcalc_i(nofgijd), stat=i_stat)
      call memocc(i_stat, product(shape(ijtabcalc_i))*kind(ijtabcalc_i), 'IJTABCALC_I', 'allocate_green')
      ijtabcalc_i = 0

    else

      if (allocated(gsh)) then
        i_all = -product(shape(gsh))*kind(gsh)
        deallocate (gsh, stat=i_stat)
        call memocc(i_stat, i_all, 'GSH', 'allocate_misc')
      end if
      if (allocated(kmesh)) then
        i_all = -product(shape(kmesh))*kind(kmesh)
        deallocate (kmesh, stat=i_stat)
        call memocc(i_stat, i_all, 'KMESH', 'allocate_misc')
      end if
      if (allocated(ilm_map)) then
        i_all = -product(shape(ilm_map))*kind(ilm_map)
        deallocate (ilm_map, stat=i_stat)
        call memocc(i_stat, i_all, 'ILM_MAP', 'allocate_misc')
      end if
      if (allocated(iqcalc)) then
        i_all = -product(shape(iqcalc))*kind(iqcalc)
        deallocate (iqcalc, stat=i_stat)
        call memocc(i_stat, i_all, 'IQCALC', 'allocate_misc')
      end if
      if (allocated(jofgij)) then
        i_all = -product(shape(jofgij))*kind(jofgij)
        deallocate (jofgij, stat=i_stat)
        call memocc(i_stat, i_all, 'JOFGIJ', 'allocate_misc')
      end if
      if (allocated(iofgij)) then
        i_all = -product(shape(iofgij))*kind(iofgij)
        deallocate (iofgij, stat=i_stat)
        call memocc(i_stat, i_all, 'IOFGIJ', 'allocate_misc')
      end if
      if (allocated(imaxsh)) then
        i_all = -product(shape(imaxsh))*kind(imaxsh)
        deallocate (imaxsh, stat=i_stat)
        call memocc(i_stat, i_all, 'IMAXSH', 'allocate_misc')
      end if
      if (allocated(ijtabsh)) then
        i_all = -product(shape(ijtabsh))*kind(ijtabsh)
        deallocate (ijtabsh, stat=i_stat)
        call memocc(i_stat, i_all, 'IJTABSH', 'allocate_misc')
      end if
      if (allocated(ijtabsym)) then
        i_all = -product(shape(ijtabsym))*kind(ijtabsym)
        deallocate (ijtabsym, stat=i_stat)
        call memocc(i_stat, i_all, 'IJTABSYM', 'allocate_misc')
      end if
      if (allocated(ijtabcalc)) then
        i_all = -product(shape(ijtabcalc))*kind(ijtabcalc)
        deallocate (ijtabcalc, stat=i_stat)
        call memocc(i_stat, i_all, 'IJTABCALC', 'allocate_misc')
      end if
      if (allocated(ish)) then
        i_all = -product(shape(ish))*kind(ish)
        deallocate (ish, stat=i_stat)
        call memocc(i_stat, i_all, 'ISH', 'allocate_misc')
      end if
      if (allocated(jsh)) then
        i_all = -product(shape(jsh))*kind(jsh)
        deallocate (jsh, stat=i_stat)
        call memocc(i_stat, i_all, 'JSH', 'allocate_misc')
      end if
      if (allocated(ijtabcalc_i)) then
        i_all = -product(shape(ijtabcalc_i))*kind(ijtabcalc_i)
        deallocate (ijtabcalc_i, stat=i_stat)
        call memocc(i_stat, i_all, 'IJTABCALC_I', 'allocate_misc')
      end if
    end if

  end subroutine allocate_green

end module memoryhandling
