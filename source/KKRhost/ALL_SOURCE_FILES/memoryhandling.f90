!-------------------------------------------------------------------------------
! MODULE: memoryhandling
!
! DESCRIPTION:
!> @brief Subroutine to handle allocation/deallocation of arrays
!> @details Module to handle the allocation of arrays to be later distributed
!> it aims to bring modularity to the memory management
!
!> @author
!> Jonathan Chico
!> @date 14.11.2017
!> @todo The number of arrays in the misc section should be reduced, and they should
!> be located in the appropriate routines
!-------------------------------------------------------------------------------
    Module memoryhandling

      Use profiling
      Use mod_datatypes, Only: dp

      Implicit None
      Private :: dp

    Contains

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_cell
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the unit cell.
!
!> @author
!> Jonathan Chico
!> @date 14.11.2017
!----------------------------------------------------------------------------
      Subroutine allocate_cell(flag, naez, nemb, natyp, cls, imt, irws, irns, &
        ntcell, refpot, kfg, kaoez, rmt, zat, rws, mtfac, rmtref, rmtrefat, &
        rmtnew, rbasis, lmxc)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: naez !< number of atoms in unit cell
        Integer, Intent (In) :: nemb !< number of 'embedding' positions
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Integer, Dimension (:), Allocatable, Intent (Inout) :: cls !< Cluster around atomic sites
        Integer, Dimension (:), Allocatable, Intent (Inout) :: imt !< R point at MT radius
        Integer, Dimension (:), Allocatable, Intent (Inout) :: irws !< R point at WS radius
        Integer, Dimension (:), Allocatable, Intent (Inout) :: irns !< Position of atoms in the unit cell in units of bravais vectors
        Integer, Dimension (:), Allocatable, Intent (Inout) :: lmxc
        Integer, Dimension (:), Allocatable, Intent (Inout) :: ntcell !< Index for WS cell
        Integer, Dimension (:), Allocatable, Intent (Inout) :: refpot !< Ref. pot. card  at position
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: kfg
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: kaoez !< atom types located at a given site
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: rmt !< Muffin-tin radius of true system
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: zat !< Nuclear charge
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: rws !< Wigner Seitz radius
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: mtfac !< Scaling factor for radius MT
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: rmtref !< Muffin-tin radius of reference system
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: rmtrefat
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: rmtnew !< Adapted muffin-tin radius
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: rbasis !< Position of atoms in the unit cell in units of bravais vectors

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then
          Allocate (refpot(naez+nemb), Stat=i_stat)
          Call memocc(i_stat, product(shape(refpot))*kind(refpot), 'REFPOT', &
            'allocate_cell')
          refpot = 0
          Allocate (rmtrefat(naez+nemb), Stat=i_stat)
          Call memocc(i_stat, product(shape(rmtrefat))*kind(rmtrefat), &
            'RMTREFAT', 'allocate_cell')
          rmtrefat = -1.E0_dp ! Signals the need for later calculation
          Allocate (kaoez(natyp,naez+nemb), Stat=i_stat)
          Call memocc(i_stat, product(shape(kaoez))*kind(kaoez), 'KAOEZ', &
            'allocate_cell')
          kaoez = 0
          Allocate (cls(naez+nemb), Stat=i_stat)
          Call memocc(i_stat, product(shape(cls))*kind(cls), 'CLS', &
            'allocate_cell')
          cls = 1
          Allocate (rbasis(3,naez+nemb), Stat=i_stat)
          Call memocc(i_stat, product(shape(rbasis))*kind(rbasis), 'RBASIS', &
            'allocate_cell')
          rbasis = 0.E0_dp
          Allocate (mtfac(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(mtfac))*kind(mtfac), 'MTFAC', &
            'allocate_cell')
          mtfac = 0.E0_dp
          Allocate (zat(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(zat))*kind(zat), 'ZAT', &
            'allocate_cell')
          zat = -1.E0_dp ! Negative value signals read-in from pot-file
          Allocate (ntcell(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(ntcell))*kind(ntcell), 'NTCELL', &
            'allocate_cell')
          ntcell = 0
          Allocate (rmtref(naez), Stat=i_stat)
          Call memocc(i_stat, product(shape(rmtref))*kind(rmtref), 'RMTREF', &
            'allocate_cell')
          rmtref = -1.E0_dp
          Allocate (irns(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(irns))*kind(irns), 'IRNS', &
            'allocate_cell')
          irns = -1 ! Negative value signals to use FPRADIUS
          Allocate (rws(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(rws))*kind(rws), 'RWS', &
            'allocate_cell')
          rws = 0.E0_dp
          Allocate (irws(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(irws))*kind(irws), 'IRWS', &
            'allocate_cell')
          irws = 0
          Allocate (rmt(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(rmt))*kind(rmt), 'RMT', &
            'allocate_cell')
          rmt = 0.E0_dp
          Allocate (rmtnew(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(rmtnew))*kind(rmtnew), 'RMTNEW', &
            'allocate_cell')
          rmtnew = 0.E0_dp
          Allocate (imt(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(imt))*kind(imt), 'IMT', &
            'allocate_cell')
          imt = 0
          Allocate (kfg(4,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(kfg))*kind(kfg), 'KFG', &
            'allocate_cell')
          kfg = 0
          Allocate (lmxc(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(lmxc))*kind(lmxc), 'LMXC', &
            'allocate_cell')
          lmxc = 0
        Else
          If (allocated(zat)) Then
            i_all = -product(shape(zat))*kind(zat)
            Deallocate (zat, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ZAT', 'allocate_cell')
          End If
          If (allocated(ntcell)) Then
            i_all = -product(shape(ntcell))*kind(ntcell)
            Deallocate (ntcell, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NTCELL', 'allocate_cell')
          End If
          If (allocated(ntcell)) Then
            i_all = -product(shape(ntcell))*kind(ntcell)
            Deallocate (ntcell, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NTCELL', 'allocate_cell')
          End If
          If (allocated(rmtref)) Then
            i_all = -product(shape(rmtref))*kind(rmtref)
            Deallocate (rmtref, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RMTREF', 'allocate_cell')
          End If
          If (allocated(irns)) Then
            i_all = -product(shape(irns))*kind(irns)
            Deallocate (irns, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IRNS', 'allocate_cell')
          End If
          If (allocated(irws)) Then
            i_all = -product(shape(irws))*kind(irws)
            Deallocate (irws, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IRWS', 'allocate_cell')
          End If
          If (allocated(rws)) Then
            i_all = -product(shape(rws))*kind(rws)
            Deallocate (rws, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RWS', 'allocate_cell')
          End If
          If (allocated(rmt)) Then
            i_all = -product(shape(rmt))*kind(rmt)
            Deallocate (rmt, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RMT', 'allocate_cell')
          End If
          If (allocated(rmtnew)) Then
            i_all = -product(shape(rmtnew))*kind(rmtnew)
            Deallocate (rmtnew, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RMTNEW', 'allocate_cell')
          End If
          If (allocated(imt)) Then
            i_all = -product(shape(imt))*kind(imt)
            Deallocate (imt, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IMT', 'allocate_cell')
          End If
          If (allocated(kfg)) Then
            i_all = -product(shape(kfg))*kind(kfg)
            Deallocate (kfg, Stat=i_stat)
            Call memocc(i_stat, i_all, 'KFG', 'allocate_cell')
          End If
          If (allocated(kaoez)) Then
            i_all = -product(shape(kaoez))*kind(kaoez)
            Deallocate (kaoez, Stat=i_stat)
            Call memocc(i_stat, i_all, 'KAOEZ', 'allocate_cell')
          End If
          If (allocated(refpot)) Then
            i_all = -product(shape(refpot))*kind(refpot)
            Deallocate (refpot, Stat=i_stat)
            Call memocc(i_stat, i_all, 'REFPOT', 'allocate_cell')
          End If
          If (allocated(rmtrefat)) Then
            i_all = -product(shape(rmtrefat))*kind(rmtrefat)
            Deallocate (rmtrefat, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RMTREFAT', 'allocate_cell')
          End If
          If (allocated(lmxc)) Then
            i_all = -product(shape(lmxc))*kind(lmxc)
            Deallocate (lmxc, Stat=i_stat)
            Call memocc(i_stat, i_all, 'LMXC', 'allocate_cell')
          End If
        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_semi_inf_host
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the left and right host for the calculation of slabs
!
!> @author
!> Jonathan Chico
!> @date 20.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_semi_inf_host(flag, nemb, tleft, tright)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: nemb !< number of 'embedding' positions
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: tleft !< Vectors of the basis for the left host
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: tright !< vectors of the basis for the right host

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then

          Allocate (tleft(3,nemb+1), Stat=i_stat)
          Call memocc(i_stat, product(shape(tleft))*kind(tleft), 'TLEFT', &
            'allocate_cell')
          tleft = 0.E0_dp
          Allocate (tright(3,nemb+1), Stat=i_stat)
          Call memocc(i_stat, product(shape(tright))*kind(tright), 'TRIGHT', &
            'allocate_cell')
          tright = 0.E0_dp
        Else
          If (allocated(tleft)) Then
            i_all = -product(shape(tleft))*kind(tleft)
            Deallocate (tleft, Stat=i_stat)
            Call memocc(i_stat, i_all, 'TLEFT', 'allocate_cell')
          End If
          If (allocated(tright)) Then
            i_all = -product(shape(tright))*kind(tright)
            Deallocate (tright, Stat=i_stat)
            Call memocc(i_stat, i_all, 'TRIGHT', 'allocate_cell')
          End If
        End If


      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_potential
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the potential
!
!> @author
!> Jonathan Chico
!> @date 14.11.2017
!----------------------------------------------------------------------------
      Subroutine allocate_potential(flag, naez, nemb, irm, natyp, npotd, &
        ipand, nfund, lmxspd, lmpot, irmind, nspotd, nfu, irc, ncore, irmin, &
        lmsp, lmsp1, ircut, lcore, llmsp, ititle, fpradius, visp, ecore, vins)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: naez !< number of atoms in unit cell
        Integer, Intent (In) :: nemb !< number of 'embedding' positions
        Integer, Intent (In) :: irm
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Integer, Intent (In) :: npotd !< 2*NATYP
        Integer, Intent (In) :: ipand
        Integer, Intent (In) :: nfund
        Integer, Intent (In) :: lmxspd
        Integer, Intent (In) :: lmpot
        Integer, Intent (In) :: irmind
        Integer, Intent (In) :: nspotd
        Integer, Dimension (:), Allocatable, Intent (Inout) :: nfu
        Integer, Dimension (:), Allocatable, Intent (Inout) :: irc !< R point for potential cutting
        Integer, Dimension (:), Allocatable, Intent (Inout) :: ncore !< Number of core states
        Integer, Dimension (:), Allocatable, Intent (Inout) :: irmin !< Max R for spherical treatment
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: lmsp !< 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: lmsp1
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: ircut !< R points of panel borders
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: lcore !< Angular momentum of core states
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: llmsp !< lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: ititle
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: fpradius !< R point at which full-potential treatment starts
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: visp !< Spherical part of the potential
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: ecore !< Core energies
        Real (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: vins !< Non-spherical part of the potential

! .. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then

          Allocate (fpradius(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(fpradius))*kind(fpradius), &
            'FPRADIUS', 'allocate_potential')
          fpradius = -1.E0_dp ! Negative value signals to use IRNS from pot-file (sub. startb1)
          Allocate (irc(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(irc))*kind(irc), 'IRC', &
            'allocate_potential')
          irc = 0
          Allocate (ircut(0:ipand,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(ircut))*kind(ircut), 'IRCUT', &
            'allocate_potential')
          ircut = 0
          Allocate (irmin(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(irmin))*kind(irmin), 'IRMIN', &
            'allocate_potential')
          irmin = 0
          Allocate (lcore(20,npotd), Stat=i_stat)
          Call memocc(i_stat, product(shape(lcore))*kind(lcore), 'LCORE', &
            'allocate_potential')
          lcore = 0
          Allocate (ecore(20,npotd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ecore))*kind(ecore), 'ECORE', &
            'allocate_potential')
          ecore = 0.E0_dp
          Allocate (ncore(npotd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ncore))*kind(ncore), 'NCORE', &
            'allocate_potential')
          ncore = 0
          Allocate (lmsp1(lmxspd,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(lmsp1))*kind(lmsp1), 'LMSP1', &
            'allocate_potential')
          lmsp1 = 0
          Allocate (llmsp(natyp,nfund), Stat=i_stat)
          Call memocc(i_stat, product(shape(llmsp))*kind(llmsp), 'LLMSP', &
            'allocate_potential')
          llmsp = 0
          Allocate (lmsp(natyp,lmxspd), Stat=i_stat)
          Call memocc(i_stat, product(shape(lmsp))*kind(lmsp), 'LMSP', &
            'allocate_potential')
          lmsp = 0
          Allocate (ititle(20,npotd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ititle))*kind(ititle), 'ITITLE', &
            'allocate_potential')
          ititle = 0
          Allocate (nfu(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(nfu))*kind(nfu), 'NFU', &
            'allocate_potential')
          nfu = 0
          Allocate (vins(irmind:irm,lmpot,nspotd), Stat=i_stat)
          Call memocc(i_stat, product(shape(vins))*kind(vins), 'VINS', &
            'allocate_misc')
          vins = 0.E0_dp
          Allocate (visp(irm,npotd), Stat=i_stat)
          Call memocc(i_stat, product(shape(visp))*kind(visp), 'VISP', &
            'allocate_misc')
          visp = 0.E0_dp

        Else
          If (allocated(fpradius)) Then
            i_all = -product(shape(fpradius))*kind(fpradius)
            Deallocate (fpradius, Stat=i_stat)
            Call memocc(i_stat, i_all, 'FPRADIUS', 'allocate_potential')
          End If
          If (allocated(irc)) Then
            i_all = -product(shape(irc))*kind(irc)
            Deallocate (irc, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IRC', 'allocate_potential')
          End If
          If (allocated(ircut)) Then
            i_all = -product(shape(ircut))*kind(ircut)
            Deallocate (ircut, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IRCUT', 'allocate_potential')
          End If
          If (allocated(irmin)) Then
            i_all = -product(shape(irmin))*kind(irmin)
            Deallocate (irmin, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IRMIN', 'allocate_potential')
          End If
          If (allocated(lcore)) Then
            i_all = -product(shape(lcore))*kind(lcore)
            Deallocate (lcore, Stat=i_stat)
            Call memocc(i_stat, i_all, 'LCORE', 'allocate_potential')
          End If
          If (allocated(lmsp1)) Then
            i_all = -product(shape(lmsp1))*kind(lmsp1)
            Deallocate (lmsp1, Stat=i_stat)
            Call memocc(i_stat, i_all, 'LMSP1', 'allocate_potential')
          End If
          If (allocated(llmsp)) Then
            i_all = -product(shape(llmsp))*kind(llmsp)
            Deallocate (llmsp, Stat=i_stat)
            Call memocc(i_stat, i_all, 'LLMSP', 'allocate_potential')
          End If
          If (allocated(lmsp)) Then
            i_all = -product(shape(lmsp))*kind(lmsp)
            Deallocate (lmsp, Stat=i_stat)
            Call memocc(i_stat, i_all, 'LMSP', 'allocate_potential')
          End If
          If (allocated(ititle)) Then
            i_all = -product(shape(ititle))*kind(ititle)
            Deallocate (ititle, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ITITLE', 'allocate_potential')
          End If
          If (allocated(nfu)) Then
            i_all = -product(shape(nfu))*kind(nfu)
            Deallocate (nfu, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NFU', 'allocate_potential')
          End If
          If (allocated(vins)) Then
            i_all = -product(shape(vins))*kind(vins)
            Deallocate (vins, Stat=i_stat)
            Call memocc(i_stat, i_all, 'VINS', 'allocate_misc')
          End If
          If (allocated(visp)) Then
            i_all = -product(shape(visp))*kind(visp)
            Deallocate (visp, Stat=i_stat)
            Call memocc(i_stat, i_all, 'VISP', 'allocate_misc')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_cpa
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the CPA treatment
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_cpa(flag, naez, nemb, natyp, noq, icpa, iqat, &
        hostimp, conc)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: naez !< number of atoms in unit cell
        Integer, Intent (In) :: nemb !< number of 'embedding' positions
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell

        Integer, Dimension (:), Allocatable, Intent (Inout) :: noq !< Number of diff. atom types located
        Integer, Dimension (:), Allocatable, Intent (Inout) :: icpa !< ICPA = 0/1 site-dependent CPA flag
        Integer, Dimension (:), Allocatable, Intent (Inout) :: iqat !< the site on which an atom is located on a given site
        Integer, Dimension (:), Allocatable, Intent (Inout) :: hostimp
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: conc !< concentration of a given atom

! .. Local variables
        Integer :: i_stat, i_all


        If (flag>0) Then

          Allocate (noq(naez), Stat=i_stat)
          Call memocc(i_stat, product(shape(noq))*kind(noq), 'NOQ', &
            'allocate_cpa')
          noq = 1
          Allocate (icpa(naez), Stat=i_stat)
          Call memocc(i_stat, product(shape(icpa))*kind(icpa), 'ICPA', &
            'allocate_cpa')
          icpa = 0
          Allocate (iqat(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(iqat))*kind(iqat), 'IQAT', &
            'allocate_cpa')
          iqat = 0
          Allocate (conc(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(conc))*kind(conc), 'CONC', &
            'allocate_cpa')
          conc = 1.E0_dp
          Allocate (hostimp(0:natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(hostimp))*kind(hostimp), &
            'HOSTIMP', 'allocate_cpa')
          hostimp = 0

        Else

          If (allocated(noq)) Then
            i_all = -product(shape(noq))*kind(noq)
            Deallocate (noq, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NOQ', 'allocate_cpa')
          End If
          If (allocated(icpa)) Then
            i_all = -product(shape(icpa))*kind(icpa)
            Deallocate (icpa, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ICPA', 'allocate_cpa')
          End If
          If (allocated(iqat)) Then
            i_all = -product(shape(iqat))*kind(iqat)
            Deallocate (iqat, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IQAT', 'allocate_cpa')
          End If
          If (allocated(conc)) Then
            i_all = -product(shape(conc))*kind(conc)
            Deallocate (conc, Stat=i_stat)
            Call memocc(i_stat, i_all, 'CONC', 'allocate_cpa')
          End If
          If (allocated(hostimp)) Then
            i_all = -product(shape(hostimp))*kind(hostimp)
            Deallocate (hostimp, Stat=i_stat)
            Call memocc(i_stat, i_all, 'HOSTIMP', 'allocate_cpa')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_ldau
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the LDA+U approach
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_ldau(flag, natyp, lopt, ueff, jeff, erefldau)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell

        Integer, Dimension (:), Allocatable, Intent (Inout) :: lopt !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: ueff !< input U parameter for each atom
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: jeff !< input J parameter for each atom
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: erefldau !< the energies of the projector's wave functions (REAL)

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then
          Allocate (lopt(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(lopt))*kind(lopt), 'LOPT', &
            'allocate_ldau')
          lopt = -1 !  not perform lda+u (default)
          Allocate (ueff(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(ueff))*kind(ueff), 'UEFF', &
            'allocate_ldau')
          ueff = 0.E0_dp
          Allocate (jeff(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(jeff))*kind(jeff), 'JEFF', &
            'allocate_ldau')
          jeff = 0.E0_dp
          Allocate (erefldau(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(erefldau))*kind(erefldau), &
            'EREFLDAU', 'allocate_ldau')
          erefldau = 0.5E0_dp
        Else
          If (allocated(lopt)) Then
            i_all = -product(shape(lopt))*kind(lopt)
            Deallocate (lopt, Stat=i_stat)
            Call memocc(i_stat, i_all, 'LOPT', 'allocate_ldau')
          End If
          If (allocated(ueff)) Then
            i_all = -product(shape(ueff))*kind(ueff)
            Deallocate (ueff, Stat=i_stat)
            Call memocc(i_stat, i_all, 'UEFF', 'allocate_ldau')
          End If
          If (allocated(jeff)) Then
            i_all = -product(shape(jeff))*kind(jeff)
            Deallocate (jeff, Stat=i_stat)
            Call memocc(i_stat, i_all, 'JEFF', 'allocate_ldau')
          End If
          If (allocated(erefldau)) Then
            i_all = -product(shape(erefldau))*kind(erefldau)
            Deallocate (erefldau, Stat=i_stat)
            Call memocc(i_stat, i_all, 'EREFLDAU', 'allocate_ldau')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_ldau_potential
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the potentials for the LDA+U approach
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_ldau_potential(flag, irm, natyp, mmaxd, nspind, &
        itldau, wldau, uldau, phildau)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: irm
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Integer, Intent (In) :: mmaxd
        Integer, Intent (In) :: nspind !< Counter for spin directions (KREL+(1-KREL)*(KSP+1))
        Integer, Dimension (:), Allocatable, Intent (Inout) :: itldau !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
        Real (Kind=dp), Dimension (:, :, :, :), Allocatable, &
          Intent (Inout) :: wldau !< potential matrix
        Real (Kind=dp), Dimension (:, :, :, :, :), Allocatable, &
          Intent (Inout) :: uldau !< calculated Coulomb matrix elements (EREFLDAU)
        Complex (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: phildau

! .. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then

          Allocate (itldau(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(itldau))*kind(itldau), 'ITLDAU', &
            'allocate_ldau_potential')
          itldau = 0
          Allocate (uldau(mmaxd,mmaxd,mmaxd,mmaxd,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(uldau))*kind(uldau), 'ULDAU', &
            'allocate_ldau_potential')
          uldau = 0.E0_dp
          Allocate (wldau(mmaxd,mmaxd,nspind,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(wldau))*kind(wldau), 'WLDAU', &
            'allocate_ldau_potential')
          wldau = 0.E0_dp
          Allocate (phildau(irm,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(phildau))*kind(phildau), &
            'PHILDAU', 'allocate_ldau_potential')
          phildau = (0.E0_dp, 0.E0_dp)

        Else
          If (allocated(itldau)) Then
            i_all = -product(shape(itldau))*kind(itldau)
            Deallocate (itldau, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ITLDAU', 'allocate_ldau_potential')
          End If
          If (allocated(uldau)) Then
            i_all = -product(shape(uldau))*kind(uldau)
            Deallocate (uldau, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ULDAU', 'allocate_ldau_potential')
          End If
          If (allocated(wldau)) Then
            i_all = -product(shape(wldau))*kind(wldau)
            Deallocate (wldau, Stat=i_stat)
            Call memocc(i_stat, i_all, 'WLDAU', 'allocate_ldau_potential')
          End If
          If (allocated(phildau)) Then
            i_all = -product(shape(phildau))*kind(phildau)
            Deallocate (phildau, Stat=i_stat)
            Call memocc(i_stat, i_all, 'PHILDAU', 'allocate_ldau_potential')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_magnetization
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the magnetisation
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_magnetization(flag, naez, natyp, lmmaxd, inipol, &
        ixipol, qmtet, qmphi, drotq)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: naez !< number of atoms in unit cell
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Integer, Intent (In) :: lmmaxd
        Integer, Dimension (:), Allocatable, Intent (Inout) :: inipol !< Initial spin polarisation
        Integer, Dimension (:), Allocatable, Intent (Inout) :: ixipol !< Constraint of spin pol.
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: qmtet
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: qmphi
        Complex (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: drotq !< Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then
          Allocate (qmtet(naez), Stat=i_stat)
          Call memocc(i_stat, product(shape(qmtet))*kind(qmtet), 'QMTET', &
            'allocate_magnetization')
          qmtet = 0.E0_dp
          Allocate (qmphi(naez), Stat=i_stat)
          Call memocc(i_stat, product(shape(qmphi))*kind(qmphi), 'QMPHI', &
            'allocate_magnetization')
          qmphi = 0.E0_dp
          Allocate (inipol(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(inipol))*kind(inipol), 'INIPOL', &
            'allocate_magnetization')
          inipol = 0
          Allocate (ixipol(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(ixipol))*kind(ixipol), 'IXIPOL', &
            'allocate_magnetization')
          ixipol = 0
          Allocate (drotq(lmmaxd,lmmaxd,naez), Stat=i_stat)
          Call memocc(i_stat, product(shape(drotq))*kind(drotq), 'DROTQ', &
            'allocate_magnetization')
          drotq = (0.E0_dp, 0.E0_dp)
        Else
          If (allocated(qmtet)) Then
            i_all = -product(shape(qmtet))*kind(qmtet)
            Deallocate (qmtet, Stat=i_stat)
            Call memocc(i_stat, i_all, 'QMTET', 'allocate_magnetization')
          End If
          If (allocated(qmphi)) Then
            i_all = -product(shape(qmphi))*kind(qmphi)
            Deallocate (qmphi, Stat=i_stat)
            Call memocc(i_stat, i_all, 'QMPHI', 'allocate_magnetization')
          End If
          If (allocated(inipol)) Then
            i_all = -product(shape(inipol))*kind(inipol)
            Deallocate (inipol, Stat=i_stat)
            Call memocc(i_stat, i_all, 'INIPOL', 'allocate_magnetization')
          End If
          If (allocated(ixipol)) Then
            i_all = -product(shape(ixipol))*kind(ixipol)
            Deallocate (ixipol, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IXIPOL', 'allocate_magnetization')
          End If
          If (allocated(drotq)) Then
            i_all = -product(shape(drotq))*kind(drotq)
            Deallocate (drotq, Stat=i_stat)
            Call memocc(i_stat, i_all, 'DROTQ', 'allocate_magnetization')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_SOC
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the spin-orbit coupling (SOC)
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_soc(flag, krel, natyp, lmax, socscale, cscl, socscl)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: krel
        Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: socscale !< Spin-orbit scaling
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: cscl !< Speed of light scaling
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: socscl

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then
          Allocate (socscl(krel*lmax+1,krel*natyp+(1-krel)), Stat=i_stat)
          Call memocc(i_stat, product(shape(socscl))*kind(socscl), 'SOCSCL', &
            'allocate_SOC')
          socscl = 1.E0_dp
          Allocate (cscl(krel*lmax+1,krel*natyp+(1-krel)), Stat=i_stat)
          Call memocc(i_stat, product(shape(cscl))*kind(cscl), 'CSCL', &
            'allocate_SOC')
          cscl = 0.E0_dp
          Allocate (socscale(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(socscale))*kind(socscale), &
            'SOCSCALE', 'allocate_SOC')
          socscale = 1.E0_dp ! Spin-orbit scaling
        Else
          If (allocated(socscl)) Then
            i_all = -product(shape(socscl))*kind(socscl)
            Deallocate (socscl, Stat=i_stat)
            Call memocc(i_stat, i_all, 'SOCSCL', 'allocate_SOC')
          End If
          If (allocated(socscale)) Then
            i_all = -product(shape(socscale))*kind(socscale)
            Deallocate (socscale, Stat=i_stat)
            Call memocc(i_stat, i_all, 'SOCSCALE', 'allocate_SOC')
          End If
          If (allocated(cscl)) Then
            i_all = -product(shape(cscl))*kind(cscl)
            Deallocate (cscl, Stat=i_stat)
            Call memocc(i_stat, i_all, 'CSCL', 'allocate_SOC')
          End If
        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_energies
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe energies
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_energies(flag, iemxd, ez, dez, wez)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: iemxd

        Complex (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: ez
        Complex (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: dez
        Complex (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: wez

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then
          Allocate (ez(iemxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ez))*kind(ez), 'EZ', &
            'allocate_energies')
          ez = (0.E0_dp, 0.E0_dp)
          Allocate (dez(iemxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(dez))*kind(dez), 'DEZ', &
            'allocate_energies')
          dez = (0.E0_dp, 0.E0_dp)
          Allocate (wez(iemxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(wez))*kind(wez), 'WEZ', &
            'allocate_energies')
          wez = (0.E0_dp, 0.E0_dp)
        Else
          If (allocated(ez)) Then
            i_all = -product(shape(ez))*kind(ez)
            Deallocate (ez, Stat=i_stat)
            Call memocc(i_stat, i_all, 'EZ', 'allocate_energies')
          End If
          If (allocated(dez)) Then
            i_all = -product(shape(dez))*kind(dez)
            Deallocate (dez, Stat=i_stat)
            Call memocc(i_stat, i_all, 'DEZ', 'allocate_energies')
          End If
          If (allocated(wez)) Then
            i_all = -product(shape(wez))*kind(wez)
            Deallocate (wez, Stat=i_stat)
            Call memocc(i_stat, i_all, 'WEZ', 'allocate_energies')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_relativistic
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe relativistic corrections
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_relativistic(flag, krel, irm, naez, natyp, zrel, &
        jwsrel, irshift, vtrel, btrel, rmrel, drdirel, r2drdirel, qmgam, &
        qmgamtab, qmphitab, qmtettab)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: krel
        Integer, Intent (In) :: irm
        Integer, Intent (In) :: naez !< number of atoms in unit cell
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Integer, Dimension (:), Allocatable, Intent (Inout) :: zrel !< atomic number (cast integer)
        Integer, Dimension (:), Allocatable, Intent (Inout) :: jwsrel !< index of the WS radius
        Integer, Dimension (:), Allocatable, Intent (Inout) :: irshift !< shift of the REL radial mesh with respect no NREL
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: qmgam
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: vtrel !< potential (spherical part)
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: btrel !< magnetic field
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: rmrel !< radial mesh
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: drdirel !< derivative of radial mesh
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: r2drdirel !< r**2 * drdi
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: qmgamtab
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: qmphitab
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: qmtettab

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then
          Allocate (vtrel(irm*krel+(1-krel),natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(vtrel))*kind(vtrel), 'VTREL', &
            'allocate_relativistic')
          vtrel = 0.E0_dp
          Allocate (btrel(irm*krel+(1-krel),natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(btrel))*kind(btrel), 'BTREL', &
            'allocate_relativistic')
          btrel = 0.E0_dp
          Allocate (drdirel(irm*krel+(1-krel),natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(drdirel))*kind(drdirel), &
            'DRDIREL', 'allocate_relativistic')
          drdirel = 0.E0_dp
          Allocate (r2drdirel(irm*krel+(1-krel),natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(r2drdirel))*kind(r2drdirel), &
            'R2DRDIREL', 'allocate_relativistic')
          r2drdirel = 0.E0_dp
          Allocate (rmrel(irm*krel+(1-krel),natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(rmrel))*kind(rmrel), 'RMREL', &
            'allocate_relativistic')
          rmrel = 0.E0_dp
          Allocate (irshift(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(irshift))*kind(irshift), &
            'IRSHIFT', 'allocate_relativistic')
          irshift = 0
          Allocate (jwsrel(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(jwsrel))*kind(jwsrel), 'JWSREL', &
            'allocate_relativistic')
          jwsrel = 0
          Allocate (zrel(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(zrel))*kind(zrel), 'ZREL', &
            'allocate_relativistic')
          zrel = 0
          Allocate (qmgam(naez), Stat=i_stat)
          Call memocc(i_stat, product(shape(qmgam))*kind(qmgam), 'QMGAM', &
            'allocate_relativistic')
          qmgam = 0.E0_dp
          Allocate (qmgamtab(naez,3), Stat=i_stat)
          Call memocc(i_stat, product(shape(qmgamtab))*kind(qmgamtab), &
            'QMGAMTAB', 'allocate_relativistic')
          qmgamtab = 0.E0_dp
          Allocate (qmphitab(naez,3), Stat=i_stat)
          Call memocc(i_stat, product(shape(qmphitab))*kind(qmphitab), &
            'QMPHITAB', 'allocate_relativistic')
          qmphitab = 0.E0_dp
          Allocate (qmtettab(naez,3), Stat=i_stat)
          Call memocc(i_stat, product(shape(qmtettab))*kind(qmtettab), &
            'QMTETTAB', 'allocate_relativistic')
          qmtettab = 0.E0_dp

        Else
          If (allocated(vtrel)) Then
            i_all = -product(shape(vtrel))*kind(vtrel)
            Deallocate (vtrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'VTREL', 'allocate_relativistic')
          End If
          If (allocated(btrel)) Then
            i_all = -product(shape(btrel))*kind(btrel)
            Deallocate (btrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'BTREL', 'allocate_relativistic')
          End If
          If (allocated(drdirel)) Then
            i_all = -product(shape(drdirel))*kind(drdirel)
            Deallocate (drdirel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'DRDIREL', 'allocate_relativistic')
          End If
          If (allocated(r2drdirel)) Then
            i_all = -product(shape(r2drdirel))*kind(r2drdirel)
            Deallocate (r2drdirel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'R2DRDIREL', 'allocate_relativistic')
          End If
          If (allocated(rmrel)) Then
            i_all = -product(shape(rmrel))*kind(rmrel)
            Deallocate (rmrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RMREL', 'allocate_relativistic')
          End If
          If (allocated(irshift)) Then
            i_all = -product(shape(irshift))*kind(irshift)
            Deallocate (irshift, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IRSHIFT', 'allocate_relativistic')
          End If
          If (allocated(jwsrel)) Then
            i_all = -product(shape(jwsrel))*kind(jwsrel)
            Deallocate (jwsrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'JWSREL', 'allocate_relativistic')
          End If
          If (allocated(zrel)) Then
            i_all = -product(shape(zrel))*kind(zrel)
            Deallocate (zrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ZREL', 'allocate_relativistic')
          End If
          If (allocated(qmgam)) Then
            i_all = -product(shape(qmgam))*kind(qmgam)
            Deallocate (qmgam, Stat=i_stat)
            Call memocc(i_stat, i_all, 'QMGAM', 'allocate_relativistic')
          End If
          If (allocated(qmgamtab)) Then
            i_all = -product(shape(qmgamtab))*kind(qmgamtab)
            Deallocate (qmgamtab, Stat=i_stat)
            Call memocc(i_stat, i_all, 'QMGAMTAB', 'allocate_relativistic')
          End If
          If (allocated(qmphitab)) Then
            i_all = -product(shape(qmphitab))*kind(qmphitab)
            Deallocate (qmphitab, Stat=i_stat)
            Call memocc(i_stat, i_all, 'QMPHITAB', 'allocate_relativistic')
          End If
          If (allocated(qmtettab)) Then
            i_all = -product(shape(qmtettab))*kind(qmtettab)
            Deallocate (qmtettab, Stat=i_stat)
            Call memocc(i_stat, i_all, 'QMTETTAB', 'allocate_relativistic')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_rel_transformations
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe relativistic transformations
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_rel_transformations(flag, lmmaxd, nrrel, irrel, rc, &
        crel, rrel, srrel)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: lmmaxd
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: nrrel
        Integer, Dimension (:, :, :), Allocatable, Intent (Inout) :: irrel
        Complex (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: rc !< NREL REAL spher. harm. >  CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
        Complex (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: crel !< Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
        Complex (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: rrel !< Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
        Complex (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: srrel

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then

          Allocate (rrel(lmmaxd,lmmaxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(rrel))*kind(rrel), 'RREL', &
            'allocate_rel_transformations')
          rrel = (0.E0_dp, 0.E0_dp)
          Allocate (srrel(2,2,lmmaxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(srrel))*kind(srrel), 'SRREL', &
            'allocate_rel_transformations')
          srrel = (0.E0_dp, 0.E0_dp)
          Allocate (irrel(2,2,lmmaxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(irrel))*kind(irrel), 'IRREL', &
            'allocate_rel_transformations')
          irrel = 0
          Allocate (nrrel(2,lmmaxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(nrrel))*kind(nrrel), 'NRREL', &
            'allocate_rel_transformations')
          nrrel = 0
          Allocate (crel(lmmaxd,lmmaxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(crel))*kind(crel), 'CREL', &
            'allocate_rel_transformations')
          crel = (0.E0_dp, 0.E0_dp)
          Allocate (rc(lmmaxd,lmmaxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(rc))*kind(rc), 'RC', &
            'allocate_rel_transformations')
          rc = (0.E0_dp, 0.E0_dp)
        Else
          If (allocated(rrel)) Then
            i_all = -product(shape(rrel))*kind(rrel)
            Deallocate (rrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RREL', 'allocate_rel_transformations')
          End If
          If (allocated(srrel)) Then
            i_all = -product(shape(srrel))*kind(srrel)
            Deallocate (srrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'SRREL', 'allocate_rel_transformations' &
              )
          End If
          If (allocated(irrel)) Then
            i_all = -product(shape(irrel))*kind(irrel)
            Deallocate (irrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IRREL', 'allocate_rel_transformations' &
              )
          End If
          If (allocated(nrrel)) Then
            i_all = -product(shape(nrrel))*kind(nrrel)
            Deallocate (nrrel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NRREL', 'allocate_rel_transformations' &
              )
          End If
          If (allocated(crel)) Then
            i_all = -product(shape(crel))*kind(crel)
            Deallocate (crel, Stat=i_stat)
            Call memocc(i_stat, i_all, 'CREL', 'allocate_rel_transformations')
          End If
          If (allocated(rc)) Then
            i_all = -product(shape(rc))*kind(rc)
            Deallocate (rc, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RC', 'allocate_rel_transformations')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_clusters
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe clusters
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_clusters(flag, naez, lmax, ncleb, nclsd, nembd1, &
        nsheld, naclsd, lmpot, natomimpd, nsh1, nsh2, nacls, nshell, atomimp, &
        atom, ezoa, icleb, jend, ratom, rclsimp, cmomhost, rcls)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: naez !< number of atoms in unit cell
        Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
        Integer, Intent (In) :: ncleb
        Integer, Intent (In) :: nclsd
        Integer, Intent (In) :: nembd1
        Integer, Intent (In) :: nsheld
        Integer, Intent (In) :: naclsd
        Integer, Intent (In) :: lmpot
        Integer, Intent (In) :: natomimpd
        Integer, Dimension (:), Allocatable, Intent (Inout) :: nsh1 !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
        Integer, Dimension (:), Allocatable, Intent (Inout) :: nsh2 !< Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
        Integer, Dimension (:), Allocatable, Intent (Inout) :: nacls !< Number of atoms in cluster
        Integer, Dimension (:), Allocatable, Intent (Inout) :: nshell !< Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
        Integer, Dimension (:), Allocatable, Intent (Inout) :: atomimp
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: atom !< Atom at site in cluster
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: ezoa !< EZ of atom at site in cluster
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: icleb !< Pointer array
        Integer, Dimension (:, :, :), Allocatable, Intent (Inout) :: jend !< Pointer array for icleb()
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: ratom
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: rclsimp
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: cmomhost !< Charge moments of each atom of the (left/right) host
        Real (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: rcls !< Real space position of atom in cluster

        Integer :: i_stat, i_all

        If (flag>0) Then

          Allocate (atom(naclsd,naez+(nembd1-1)), Stat=i_stat)
          Call memocc(i_stat, product(shape(atom))*kind(atom), 'ATOM', &
            'allocate_clusters')
          atom = 0
          Allocate (ratom(3,nsheld), Stat=i_stat)
          Call memocc(i_stat, product(shape(ratom))*kind(ratom), 'RATOM', &
            'allocate_clusters')
          ratom = 0.E0_dp
          Allocate (rcls(3,naclsd,nclsd), Stat=i_stat)
          Call memocc(i_stat, product(shape(rcls))*kind(rcls), 'RCLS', &
            'allocate_clusters')
          rcls = 0.E0_dp
          Allocate (rclsimp(3,natomimpd), Stat=i_stat)
          Call memocc(i_stat, product(shape(rclsimp))*kind(rclsimp), &
            'RCLSIMP', 'allocate_clusters')
          rclsimp = 0.E0_dp
          Allocate (nacls(nclsd), Stat=i_stat)
          Call memocc(i_stat, product(shape(nacls))*kind(nacls), 'NACLS', &
            'allocate_clusters')
          nacls = 0
          Allocate (ezoa(naclsd,naez+(nembd1-1)), Stat=i_stat)
          Call memocc(i_stat, product(shape(ezoa))*kind(ezoa), 'EZOA', &
            'allocate_clusters')
          ezoa = 0
          Allocate (atomimp(natomimpd), Stat=i_stat)
          Call memocc(i_stat, product(shape(atomimp))*kind(atomimp), &
            'ATOMIMP', 'allocate_clusters')
          atomimp = 0
          Allocate (icleb(ncleb,4), Stat=i_stat)
          Call memocc(i_stat, product(shape(icleb))*kind(icleb), 'ICLEB', &
            'allocate_clusters')
          icleb = 0
          Allocate (nsh1(nsheld), Stat=i_stat)
          Call memocc(i_stat, product(shape(nsh1))*kind(nsh1), 'NSH1', &
            'allocate_clusters')
          nsh1 = 0
          Allocate (nsh2(nsheld), Stat=i_stat)
          Call memocc(i_stat, product(shape(nsh2))*kind(nsh2), 'NSH2', &
            'allocate_clusters')
          nsh2 = 0
          Allocate (nshell(0:nsheld), Stat=i_stat)
          Call memocc(i_stat, product(shape(nshell))*kind(nshell), 'NSHELL', &
            'allocate_clusters')
          nshell = 0
          Allocate (cmomhost(lmpot,nembd1), Stat=i_stat)
          Call memocc(i_stat, product(shape(cmomhost))*kind(cmomhost), &
            'CMOMHOST', 'allocate_clusters')
          cmomhost = 0.E0_dp
          Allocate (jend(lmpot,0:lmax,0:lmax), Stat=i_stat)
          Call memocc(i_stat, product(shape(jend))*kind(jend), 'JEND', &
            'allocate_clusters')
          jend = 0

        Else
          If (allocated(atom)) Then
            i_all = -product(shape(atom))*kind(atom)
            Deallocate (atom, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ATOM', 'allocate_clusters')
          End If
          If (allocated(ratom)) Then
            i_all = -product(shape(ratom))*kind(ratom)
            Deallocate (ratom, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RATOM', 'allocate_clusters')
          End If
          If (allocated(rcls)) Then
            i_all = -product(shape(rcls))*kind(rcls)
            Deallocate (rcls, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RCLS', 'allocate_clusters')
          End If
          If (allocated(rclsimp)) Then
            i_all = -product(shape(rclsimp))*kind(rclsimp)
            Deallocate (rclsimp, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RCLSIMP', 'allocate_clusters')
          End If
          If (allocated(nacls)) Then
            i_all = -product(shape(nacls))*kind(nacls)
            Deallocate (nacls, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NACLS', 'allocate_clusters')
          End If
          If (allocated(ezoa)) Then
            i_all = -product(shape(ezoa))*kind(ezoa)
            Deallocate (ezoa, Stat=i_stat)
            Call memocc(i_stat, i_all, 'EZOA', 'allocate_clusters')
          End If
          If (allocated(atomimp)) Then
            i_all = -product(shape(atomimp))*kind(atomimp)
            Deallocate (atomimp, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ATOMIMP', 'allocate_clusters')
          End If
          If (allocated(icleb)) Then
            i_all = -product(shape(icleb))*kind(icleb)
            Deallocate (icleb, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ICLEB', 'allocate_clusters')
          End If
          If (allocated(nsh1)) Then
            i_all = -product(shape(nsh1))*kind(nsh1)
            Deallocate (nsh1, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NSH1', 'allocate_clusters')
          End If
          If (allocated(nsh2)) Then
            i_all = -product(shape(nsh2))*kind(nsh2)
            Deallocate (nsh2, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NSH2', 'allocate_clusters')
          End If
          If (allocated(nshell)) Then
            i_all = -product(shape(nshell))*kind(nshell)
            Deallocate (nshell, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NSHELL', 'allocate_clusters')
          End If
          If (allocated(cmomhost)) Then
            i_all = -product(shape(cmomhost))*kind(cmomhost)
            Deallocate (cmomhost, Stat=i_stat)
            Call memocc(i_stat, i_all, 'CMOMHOST', 'allocate_clusters')
          End If
          If (allocated(jend)) Then
            i_all = -product(shape(jend))*kind(jend)
            Deallocate (jend, Stat=i_stat)
            Call memocc(i_stat, i_all, 'JEND', 'allocate_clusters')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_expansion
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the functions for the expansion of the Green function
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_expansion(flag, lm2d, irid, nfund, ntotd, ncleb, &
        lassld, ncelld, nchebd, loflm, wg, cleb, yrg, thetas, thetasnew)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: lm2d
        Integer, Intent (In) :: irid
        Integer, Intent (In) :: nfund
        Integer, Intent (In) :: ntotd
        Integer, Intent (In) :: ncleb
        Integer, Intent (In) :: lassld
        Integer, Intent (In) :: ncelld
        Integer, Intent (In) :: nchebd
        Integer, Dimension (:), Allocatable, Intent (Inout) :: loflm !< l of lm=(l,m) (GAUNT)
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: wg !< Integr. weights for Legendre polynomials
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: cleb !< GAUNT coefficients (GAUNT)
        Real (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: yrg !< Spherical harmonics (GAUNT2)
        Real (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: thetas !< shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
        Real (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: thetasnew

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then

          Allocate (wg(lassld), Stat=i_stat)
          Call memocc(i_stat, product(shape(wg))*kind(wg), 'WG', &
            'allocate_expansion')
          wg = 0.E0_dp
          Allocate (yrg(lassld,0:lassld,0:lassld), Stat=i_stat)
          Call memocc(i_stat, product(shape(yrg))*kind(yrg), 'YRG', &
            'allocate_expansion')
          yrg = 0.E0_dp
          Allocate (thetas(irid,nfund,ncelld), Stat=i_stat)
          Call memocc(i_stat, product(shape(thetas))*kind(thetas), 'THETAS', &
            'allocate_expansion')
          thetas = 0.E0_dp
          Allocate (thetasnew(ntotd*(nchebd+1),nfund,ncelld), Stat=i_stat)
          Call memocc(i_stat, product(shape(thetasnew))*kind(thetasnew), &
            'THETASNEW', 'allocate_expansion')
          thetasnew = 0.E0_dp
          Allocate (cleb(ncleb,2), Stat=i_stat)
          Call memocc(i_stat, product(shape(cleb))*kind(cleb), 'CLEB', &
            'allocate_expansion')
          cleb = 0.E0_dp
          Allocate (loflm(lm2d), Stat=i_stat)
          Call memocc(i_stat, product(shape(loflm))*kind(loflm), 'LOFLM', &
            'allocate_expansion')
          loflm = 0

        Else
          If (allocated(wg)) Then
            i_all = -product(shape(wg))*kind(wg)
            Deallocate (wg, Stat=i_stat)
            Call memocc(i_stat, i_all, 'WG', 'allocate_expansion')
          End If
          If (allocated(yrg)) Then
            i_all = -product(shape(yrg))*kind(yrg)
            Deallocate (yrg, Stat=i_stat)
            Call memocc(i_stat, i_all, 'YRG', 'allocate_expansion')
          End If
          If (allocated(thetas)) Then
            i_all = -product(shape(thetas))*kind(thetas)
            Deallocate (thetas, Stat=i_stat)
            Call memocc(i_stat, i_all, 'THETAS', 'allocate_expansion')
          End If
          If (allocated(thetasnew)) Then
            i_all = -product(shape(thetasnew))*kind(thetasnew)
            Deallocate (thetasnew, Stat=i_stat)
            Call memocc(i_stat, i_all, 'THETASNEW', 'allocate_expansion')
          End If
          If (allocated(cleb)) Then
            i_all = -product(shape(cleb))*kind(cleb)
            Deallocate (cleb, Stat=i_stat)
            Call memocc(i_stat, i_all, 'CLEB', 'allocate_expansion')
          End If
          If (allocated(loflm)) Then
            i_all = -product(shape(loflm))*kind(loflm)
            Deallocate (loflm, Stat=i_stat)
            Call memocc(i_stat, i_all, 'LOFLM', 'allocate_expansion')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_mesh
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the integration mesh
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_mesh(flag, irm, natyp, a, b, r, drdi)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: irm
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: a !< Constants for exponential R mesh
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: b
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: r !< Radial mesh ( in units a Bohr)
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: drdi !< Derivative dr/di

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then

          Allocate (drdi(irm,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(drdi))*kind(drdi), 'DRDI', &
            'allocate_mesh')
          drdi = 0.E0_dp
          Allocate (r(irm,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(r))*kind(r), 'R', 'allocate_mesh')
          r = 0.E0_dp
          Allocate (a(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(a))*kind(a), 'A', 'allocate_mesh')
          a = 0.E0_dp
          Allocate (b(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(b))*kind(b), 'B', 'allocate_mesh')
          b = 0.E0_dp

        Else
          If (allocated(drdi)) Then
            i_all = -product(shape(drdi))*kind(drdi)
            Deallocate (drdi, Stat=i_stat)
            Call memocc(i_stat, i_all, 'DRDI', 'allocate_mesh')
          End If
          If (allocated(r)) Then
            i_all = -product(shape(r))*kind(r)
            Deallocate (r, Stat=i_stat)
            Call memocc(i_stat, i_all, 'R', 'allocate_mesh')
          End If
          If (allocated(a)) Then
            i_all = -product(shape(a))*kind(a)
            Deallocate (a, Stat=i_stat)
            Call memocc(i_stat, i_all, 'A', 'allocate_mesh')
          End If
          If (allocated(b)) Then
            i_all = -product(shape(b))*kind(b)
            Deallocate (b, Stat=i_stat)
            Call memocc(i_stat, i_all, 'B', 'allocate_mesh')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_pannels
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays that
!> describe the pannels
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_pannels(flag, natyp, ntotd, ipan, npan_tot, &
        npan_eq_at, npan_log_at, ipan_intervall, rpan_intervall)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Integer, Intent (In) :: ntotd
        Integer, Dimension (:), Allocatable, Intent (Inout) :: ipan !< Number of panels in non-MT-region
        Integer, Dimension (:), Allocatable, Intent (Inout) :: npan_tot
        Integer, Dimension (:), Allocatable, Intent (Inout) :: npan_eq_at
        Integer, Dimension (:), Allocatable, Intent (Inout) :: npan_log_at
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: &
          ipan_intervall
        Real (Kind=dp), Dimension (:, :), Allocatable, &
          Intent (Inout) :: rpan_intervall

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then
          Allocate (ipan(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(ipan))*kind(ipan), 'IPAN', &
            'allocate_pannels')
          ipan = 0
          Allocate (npan_tot(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(npan_tot))*kind(npan_tot), &
            'NPAN_TOT', 'allocate_pannels')
          npan_tot = 0
          Allocate (npan_eq_at(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(npan_eq_at))*kind(npan_eq_at), &
            'NPAN_EQ_AT', 'allocate_pannels')
          npan_eq_at = 0
          Allocate (npan_log_at(natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(npan_log_at))*kind(npan_log_at), &
            'NPAN_LOG_AT', 'allocate_pannels')
          npan_log_at = 0
          Allocate (rpan_intervall(0:ntotd,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(rpan_intervall))*kind( &
            rpan_intervall), 'RPAN_INTERVALL', 'allocate_pannels')
          rpan_intervall = 0.E0_dp
          Allocate (ipan_intervall(0:ntotd,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(ipan_intervall))*kind( &
            ipan_intervall), 'IPAN_INTERVALL', 'allocate_pannels')
          ipan_intervall = 0

        Else
          If (allocated(ipan)) Then
            i_all = -product(shape(ipan))*kind(ipan)
            Deallocate (ipan, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IPAN', 'allocate_pannels')
          End If
          If (allocated(npan_tot)) Then
            i_all = -product(shape(npan_tot))*kind(npan_tot)
            Deallocate (npan_tot, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NPAN_TOT', 'allocate_pannels')
          End If
          If (allocated(npan_eq_at)) Then
            i_all = -product(shape(npan_eq_at))*kind(npan_eq_at)
            Deallocate (npan_eq_at, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NPAN_EQ_AT', 'allocate_pannels')
          End If
          If (allocated(npan_log_at)) Then
            i_all = -product(shape(npan_log_at))*kind(npan_log_at)
            Deallocate (npan_log_at, Stat=i_stat)
            Call memocc(i_stat, i_all, 'NPAN_LOG_AT', 'allocate_pannels')
          End If
          If (allocated(rpan_intervall)) Then
            i_all = -product(shape(rpan_intervall))*kind(rpan_intervall)
            Deallocate (rpan_intervall, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RPAN_INTERVALL', 'allocate_pannels')
          End If
          If (allocated(ipan_intervall)) Then
            i_all = -product(shape(ipan_intervall))*kind(ipan_intervall)
            Deallocate (ipan_intervall, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IPAN_INTERVALL', 'allocate_pannels')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_misc
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of misc arrays
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_misc(flag, nr, irm, irid, lmax, naez, natyp, nfund, &
        nrefd, iemxd, ntotd, nsheld, lmmaxd, nembd1, nchebd, ncelld, lmxspd, &
        nspindd, nsymaxd, nprincd, ifunm, ifunm1, icheck, vref, s, rr, dror, &
        rnew, rs, rrot, thesme, dsymll, dsymll1, lefttinvll, righttinvll)

        Implicit None

        Integer, Intent (In) :: nr
        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: irm
        Integer, Intent (In) :: irid
        Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
        Integer, Intent (In) :: naez !< number of atoms in unit cell
        Integer, Intent (In) :: natyp !< number of kinds of atoms in unit cell
        Integer, Intent (In) :: nfund
        Integer, Intent (In) :: nrefd
        Integer, Intent (In) :: iemxd
        Integer, Intent (In) :: ntotd
        Integer, Intent (In) :: lmmaxd
        Integer, Intent (In) :: nsheld
        Integer, Intent (In) :: nembd1
        Integer, Intent (In) :: nchebd
        Integer, Intent (In) :: ncelld
        Integer, Intent (In) :: lmxspd
        Integer, Intent (In) :: nspindd
        Integer, Intent (In) :: nsymaxd
        Integer, Intent (In) :: nprincd
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: ifunm
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: ifunm1
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: icheck
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: vref
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: s
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: rr !< Set of real space vectors (in a.u.)
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: dror
        Real (Kind=dp), Dimension (:, :), Allocatable, Intent (Inout) :: rnew
        Real (Kind=dp), Dimension (:, :, :), Allocatable, Intent (Inout) :: rs
        Real (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: rrot
        Real (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: thesme
        Complex (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: dsymll
        Complex (Kind=dp), Dimension (:, :, :), Allocatable, &
          Intent (Inout) :: dsymll1
        Complex (Kind=dp), Dimension (:, :, :, :, :), Allocatable, &
          Intent (Inout) :: lefttinvll
        Complex (Kind=dp), Dimension (:, :, :, :, :), Allocatable, &
          Intent (Inout) :: righttinvll

!.. Local variables
        Integer :: i_stat, i_all

        If (flag>0) Then
          Allocate (rr(3,0:nr), Stat=i_stat)
          Call memocc(i_stat, product(shape(rr))*kind(rr), 'RR', &
            'allocate_misc')
          rr = 0.E0_dp
          Allocate (vref(nrefd), Stat=i_stat)
          Call memocc(i_stat, product(shape(vref))*kind(vref), 'VREF', &
            'allocate_misc')
          vref = 0.E0_dp
          Allocate (dror(irm,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(dror))*kind(dror), 'DROR', &
            'allocate_misc')
          dror = 0.E0_dp
          Allocate (s(0:lmax,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(s))*kind(s), 'S', 'allocate_misc')
          s = 0.E0_dp
          Allocate (rrot(48,3,nsheld), Stat=i_stat)
          Call memocc(i_stat, product(shape(rrot))*kind(rrot), 'RROT', &
            'allocate_misc')
          rrot = 0.E0_dp
          Allocate (ifunm(natyp,lmxspd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ifunm))*kind(ifunm), 'IFUNM', &
            'allocate_misc')
          ifunm = 0
          Allocate (ifunm1(lmxspd,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(ifunm1))*kind(ifunm1), 'IFUNM1', &
            'allocate_misc')
          ifunm1 = 0
          Allocate (rs(irm,0:lmax,natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(rs))*kind(rs), 'RS', &
            'allocate_misc')
          rs = 0.E0_dp
          Allocate (thesme(irid,nfund,ncelld), Stat=i_stat)
          Call memocc(i_stat, product(shape(thesme))*kind(thesme), 'THESME', &
            'allocate_misc')
          thesme = 0.E0_dp
          Allocate (rnew(ntotd*(nchebd+1),natyp), Stat=i_stat)
          Call memocc(i_stat, product(shape(rnew))*kind(rnew), 'RNEW', &
            'allocate_misc')
          rnew = 0.E0_dp
          Allocate (dsymll(lmmaxd,lmmaxd,nsymaxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(dsymll))*kind(dsymll), 'DSYMLL', &
            'allocate_misc')
          dsymll = (0.E0_dp, 0.E0_dp)
          Allocate (dsymll1(lmmaxd,lmmaxd,nsymaxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(dsymll1))*kind(dsymll1), &
            'DSYMLL1', 'allocate_misc')
          dsymll1 = (0.E0_dp, 0.E0_dp)
          Allocate (icheck(naez/nprincd,naez/nprincd), Stat=i_stat)
          Call memocc(i_stat, product(shape(icheck))*kind(icheck), 'ICHECK', &
            'allocate_misc')
          icheck = 0
          Allocate (lefttinvll(lmmaxd,lmmaxd,nembd1,nspindd,iemxd), &
            Stat=i_stat)
          Call memocc(i_stat, product(shape(lefttinvll))*kind(lefttinvll), &
            'LEFTTINVLL', 'allocate_misc')
          lefttinvll = (0.E0_dp, 0.E0_dp)
          Allocate (righttinvll(lmmaxd,lmmaxd,nembd1,nspindd,iemxd), &
            Stat=i_stat)
          Call memocc(i_stat, product(shape(righttinvll))*kind(righttinvll), &
            'RIGHTTINVLL', 'allocate_misc')
          righttinvll = (0.E0_dp, 0.E0_dp)

        Else
          If (allocated(rr)) Then
            i_all = -product(shape(rr))*kind(rr)
            Deallocate (rr, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RR', 'allocate_misc')
          End If
          If (allocated(vref)) Then
            i_all = -product(shape(vref))*kind(vref)
            Deallocate (vref, Stat=i_stat)
            Call memocc(i_stat, i_all, 'VREF', 'allocate_misc')
          End If
          If (allocated(dror)) Then
            i_all = -product(shape(dror))*kind(dror)
            Deallocate (dror, Stat=i_stat)
            Call memocc(i_stat, i_all, 'DROR', 'allocate_misc')
          End If
          If (allocated(s)) Then
            i_all = -product(shape(s))*kind(s)
            Deallocate (s, Stat=i_stat)
            Call memocc(i_stat, i_all, 'S', 'allocate_misc')
          End If
          If (allocated(rrot)) Then
            i_all = -product(shape(rrot))*kind(rrot)
            Deallocate (rrot, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RROT', 'allocate_misc')
          End If
          If (allocated(ifunm)) Then
            i_all = -product(shape(ifunm))*kind(ifunm)
            Deallocate (ifunm, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IFUNM', 'allocate_misc')
          End If
          If (allocated(ifunm1)) Then
            i_all = -product(shape(ifunm1))*kind(ifunm1)
            Deallocate (ifunm1, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IFUNM1', 'allocate_misc')
          End If
          If (allocated(rs)) Then
            i_all = -product(shape(rs))*kind(rs)
            Deallocate (rs, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RS', 'allocate_misc')
          End If
          If (allocated(thesme)) Then
            i_all = -product(shape(thesme))*kind(thesme)
            Deallocate (thesme, Stat=i_stat)
            Call memocc(i_stat, i_all, 'THESME', 'allocate_misc')
          End If
          If (allocated(rnew)) Then
            i_all = -product(shape(rnew))*kind(rnew)
            Deallocate (rnew, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RNEW', 'allocate_misc')
          End If
          If (allocated(dsymll)) Then
            i_all = -product(shape(dsymll))*kind(dsymll)
            Deallocate (dsymll, Stat=i_stat)
            Call memocc(i_stat, i_all, 'DSYMLL', 'allocate_misc')
          End If
          If (allocated(dsymll1)) Then
            i_all = -product(shape(dsymll1))*kind(dsymll1)
            Deallocate (dsymll1, Stat=i_stat)
            Call memocc(i_stat, i_all, 'DSYMLL1', 'allocate_misc')
          End If
          If (allocated(icheck)) Then
            i_all = -product(shape(icheck))*kind(icheck)
            Deallocate (icheck, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ICHECK', 'allocate_misc')
          End If
          If (allocated(lefttinvll)) Then
            i_all = -product(shape(lefttinvll))*kind(lefttinvll)
            Deallocate (lefttinvll, Stat=i_stat)
            Call memocc(i_stat, i_all, 'LEFTTINVLL', 'allocate_misc')
          End If
          If (allocated(righttinvll)) Then
            i_all = -product(shape(righttinvll))*kind(righttinvll)
            Deallocate (righttinvll, Stat=i_stat)
            Call memocc(i_stat, i_all, 'RIGHTTINVLL', 'allocate_misc')
          End If

        End If

      End Subroutine

!----------------------------------------------------------------------------
! SUBROUTINE: allocate_green
!
! DESCRIPTION:
!> @brief subroutine handling the allocation/deallocation of arrays handling
!> the Green functions
!
!> @author
!> Jonathan Chico
!> @date 19.12.2017
!----------------------------------------------------------------------------
      Subroutine allocate_green(flag, naez, iemxd, ngshd, nsheld, lmpot, &
        nofgijd, ish, jsh, kmesh, imaxsh, iqcalc, iofgij, jofgij, ijtabsh, &
        ijtabsym, ijtabcalc, ijtabcalc_i, ilm_map, gsh)

        Implicit None

        Integer, Intent (In) :: flag ! Allocate/deallocate (1/-1) arrays
        Integer, Intent (In) :: naez !< number of atoms in unit cell
        Integer, Intent (In) :: iemxd
        Integer, Intent (In) :: ngshd
        Integer, Intent (In) :: nsheld
        Integer, Intent (In) :: lmpot
        Integer, Intent (In) :: nofgijd
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: ish
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: jsh
        Integer, Dimension (:), Allocatable, Intent (Inout) :: kmesh
        Integer, Dimension (:), Allocatable, Intent (Inout) :: imaxsh
        Integer, Dimension (:), Allocatable, Intent (Inout) :: iqcalc
        Integer, Dimension (:), Allocatable, Intent (Inout) :: iofgij !< Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
        Integer, Dimension (:), Allocatable, Intent (Inout) :: jofgij !< Linear pointers, similar to NSH1/NSH2 but giving the actual index of sites I,J = 1,NATOMIMP in the cluster
        Integer, Dimension (:), Allocatable, Intent (Inout) :: ijtabsh !< Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
        Integer, Dimension (:), Allocatable, Intent (Inout) :: ijtabsym !< Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
        Integer, Dimension (:), Allocatable, Intent (Inout) :: ijtabcalc !< Linear pointer, specifying whether the block (i,j) has to be calculated needs set up for ICC=-1, not used for ICC=1
        Integer, Dimension (:), Allocatable, Intent (Inout) :: ijtabcalc_i
        Integer, Dimension (:, :), Allocatable, Intent (Inout) :: ilm_map
        Real (Kind=dp), Dimension (:), Allocatable, Intent (Inout) :: gsh

        Integer :: i_stat, i_all

        If (flag>0) Then

          Allocate (gsh(ngshd), Stat=i_stat)
          Call memocc(i_stat, product(shape(gsh))*kind(gsh), 'GSH', &
            'allocate_green')
          gsh = 0.E0_dp
          Allocate (kmesh(iemxd), Stat=i_stat)
          Call memocc(i_stat, product(shape(kmesh))*kind(kmesh), 'KMESH', &
            'allocate_green')
          kmesh = 0
          Allocate (ilm_map(ngshd,3), Stat=i_stat)
          Call memocc(i_stat, product(shape(ilm_map))*kind(ilm_map), &
            'ILM_MAP', 'allocate_green')
          ilm_map = 0
          Allocate (iqcalc(naez), Stat=i_stat)
          Call memocc(i_stat, product(shape(iqcalc))*kind(iqcalc), 'IQCALC', &
            'allocate_green')
          iqcalc = 0
          Allocate (jofgij(nofgijd), Stat=i_stat)
          Call memocc(i_stat, product(shape(jofgij))*kind(jofgij), 'JOFGIJ', &
            'allocate_green')
          jofgij = 0
          Allocate (iofgij(nofgijd), Stat=i_stat)
          Call memocc(i_stat, product(shape(iofgij))*kind(iofgij), 'IOFGIJ', &
            'allocate_green')
          iofgij = 0
          Allocate (imaxsh(0:lmpot), Stat=i_stat)
          Call memocc(i_stat, product(shape(imaxsh))*kind(imaxsh), 'IMAXSH', &
            'allocate_green')
          imaxsh = 0
          Allocate (ijtabsh(nofgijd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ijtabsh))*kind(ijtabsh), &
            'IJTABSH', 'allocate_green')
          ijtabsh = 0
          Allocate (ijtabsym(nofgijd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ijtabsym))*kind(ijtabsym), &
            'IJTABSYM', 'allocate_green')
          ijtabsym = 0
          Allocate (ijtabcalc(nofgijd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ijtabcalc))*kind(ijtabcalc), &
            'IJTABCALC', 'allocate_green')
          ijtabcalc = 0
          Allocate (ish(nsheld,nofgijd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ish))*kind(ish), 'ISH', &
            'allocate_green')
          ish = 0
          Allocate (jsh(nsheld,nofgijd), Stat=i_stat)
          Call memocc(i_stat, product(shape(jsh))*kind(jsh), 'JSH', &
            'allocate_green')
          jsh = 0
          Allocate (ijtabcalc_i(nofgijd), Stat=i_stat)
          Call memocc(i_stat, product(shape(ijtabcalc_i))*kind(ijtabcalc_i), &
            'IJTABCALC_I', 'allocate_green')
          ijtabcalc_i = 0

        Else

          If (allocated(gsh)) Then
            i_all = -product(shape(gsh))*kind(gsh)
            Deallocate (gsh, Stat=i_stat)
            Call memocc(i_stat, i_all, 'GSH', 'allocate_misc')
          End If
          If (allocated(kmesh)) Then
            i_all = -product(shape(kmesh))*kind(kmesh)
            Deallocate (kmesh, Stat=i_stat)
            Call memocc(i_stat, i_all, 'KMESH', 'allocate_misc')
          End If
          If (allocated(ilm_map)) Then
            i_all = -product(shape(ilm_map))*kind(ilm_map)
            Deallocate (ilm_map, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ILM_MAP', 'allocate_misc')
          End If
          If (allocated(iqcalc)) Then
            i_all = -product(shape(iqcalc))*kind(iqcalc)
            Deallocate (iqcalc, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IQCALC', 'allocate_misc')
          End If
          If (allocated(jofgij)) Then
            i_all = -product(shape(jofgij))*kind(jofgij)
            Deallocate (jofgij, Stat=i_stat)
            Call memocc(i_stat, i_all, 'JOFGIJ', 'allocate_misc')
          End If
          If (allocated(iofgij)) Then
            i_all = -product(shape(iofgij))*kind(iofgij)
            Deallocate (iofgij, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IOFGIJ', 'allocate_misc')
          End If
          If (allocated(imaxsh)) Then
            i_all = -product(shape(imaxsh))*kind(imaxsh)
            Deallocate (imaxsh, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IMAXSH', 'allocate_misc')
          End If
          If (allocated(ijtabsh)) Then
            i_all = -product(shape(ijtabsh))*kind(ijtabsh)
            Deallocate (ijtabsh, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IJTABSH', 'allocate_misc')
          End If
          If (allocated(ijtabsym)) Then
            i_all = -product(shape(ijtabsym))*kind(ijtabsym)
            Deallocate (ijtabsym, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IJTABSYM', 'allocate_misc')
          End If
          If (allocated(ijtabcalc)) Then
            i_all = -product(shape(ijtabcalc))*kind(ijtabcalc)
            Deallocate (ijtabcalc, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IJTABCALC', 'allocate_misc')
          End If
          If (allocated(ish)) Then
            i_all = -product(shape(ish))*kind(ish)
            Deallocate (ish, Stat=i_stat)
            Call memocc(i_stat, i_all, 'ISH', 'allocate_misc')
          End If
          If (allocated(jsh)) Then
            i_all = -product(shape(jsh))*kind(jsh)
            Deallocate (jsh, Stat=i_stat)
            Call memocc(i_stat, i_all, 'JSH', 'allocate_misc')
          End If
          If (allocated(ijtabcalc_i)) Then
            i_all = -product(shape(ijtabcalc_i))*kind(ijtabcalc_i)
            Deallocate (ijtabcalc_i, Stat=i_stat)
            Call memocc(i_stat, i_all, 'IJTABCALC_I', 'allocate_misc')
          End If
        End If

      End Subroutine

    End Module
