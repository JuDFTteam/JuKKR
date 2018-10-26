!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_brydbm

  private
  public :: brydbm

contains

  !-------------------------------------------------------------------------------
  !> Summary: Broyden and Anderson mixing
  !> Author: S. Bluegel, B. Drittler
  !> Date: May 1987
  !> Category: KKRhost, mixing
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> imix :
  !>     3      broyden's            f i r s t  m e t h o d
  !>     4      broyden's          s e c o n d  m e t h o d
  !>     5      anderson's     g e n e r a l i z e d   m e t h o d
  !>
  !> Implemented here according to notes of S.B.
  !> broyden's iteration scheme following the papers of :
  !> Srivastava, J. Phys. , 17 (1984) , pp l317
  !> C.G. Broyden in Math.Comput., 19 , pp 577, 1965
  !> C.G. Broyden in IBID, 21 ,pp 368 ,1967
  !>
  !> The method has been generalized to include a metric. The
  !> definition of the necessary inner products are similar to the
  !> discription given in the notes of M.Weinert. The algorithm
  !> discribed in the paper Srivastava has been simplified
  !> (see notes of S.B.)
  !>
  !> The files ui,vi are stored on high speed ssd memory.
  !> Broyden's update treats charge and spin on the same footing
  !> S. Bluegel , KFA , May 1987
  !> 
  !> The anderson method (d.g. Anderson, J. ACM 12, 547 (1964)) has
  !> been generalized and reformulated as an improvement of Broyden's
  !> second method. Successive linesearch is replaced by successive
  !> search on hyperplanes. ( see notes of S.B. )
  !>
  !> Modified for non spherical potential
  !> B. Drittler , Aug. 1988
  !-------------------------------------------------------------------------------
  subroutine brydbm(visp, v, vins, vspsme, vspsmo, ins, lmpot, r, drdi, alpha, atwght, irc, irmin, nspin, natps, natyp, itdept, imix, iobroy, ipf, lsmear)

    use :: mod_types, only: t_inc
    use :: mod_datatypes, only: dp
    use :: global_variables, only: irmd, lmpotd, irmind, ntperd, nspindd, natypd, irnsd
    use :: mod_brysh1, only: brysh1
    use :: mod_brysh2, only: brysh2
    use :: mod_brysh3, only: brysh3
    use :: mod_rcstop, only: rcstop

    implicit none

    ! .. Parameters ..
    integer :: ntird
    integer :: itdthd
    parameter (itdthd=40)
    ! ..
    ! .. External Functions ..
    real (kind=dp) :: atwght(*), drdi(irmd, *), r(irmd, *), v(irmd, lmpotd, *), vins(irmind:irmd, lmpotd, *), visp(irmd, *), vspsmo(irmd, *), vspsme(irmd, *)
    integer :: irc(*), irmin(*)
    ! ..
    ! .. External Subroutines ..
    real (kind=dp) :: cmm, one, rmixiv, smnorm, vmdeno, vmnorm, volinv, zero
    integer :: ia, ij, imap, ir, irc1, irmin1, isp, it, lm, mit, lsmear
    ! ..
    ! .. Intrinsic Functions ..
    real (kind=dp) :: ddot
    external :: ddot
    ! ..
    ! .. Local Arrays ..
    intrinsic :: abs
    ! ..
    ! .. Scalar Arguments ..
    save :: mit, zero, one, wit
    ! ..
    ! .. Data statements ..
    real (kind=dp), allocatable :: am(:), bm(:), fm(:), fm1(:), g(:), sm(:), sm1(:), vi3(:), wit(:), ui2(:), ui3(:), vi2(:)


    real (kind=dp) :: alpha
    integer :: imix, ins, iobroy, ipf, itdept, lmpot, natps, natyp, nspin


    data mit/1/, zero, one/0.0d0, 1.0d0/

    ntird = (irmd*ntperd+(irnsd+1)*(lmpotd-1)*natypd)*nspindd

    allocate (am(2:itdthd-1), bm(2:itdthd-1), fm(ntird), fm1(ntird), g(ntird), sm(ntird), sm1(ntird), vi3(ntird), wit(2:200), ui2(ntird), ui3(ntird), vi2(ntird))

    mit = t_inc%mit_bry

    if (itdept>itdthd .or. itdthd>200) call rcstop('itdbry  ')
    if (imix<=2 .or. imix>5) call rcstop('IMIXD   ')

    if (mit>itdept) mit = 1
    if (imix==3) write (ipf, fmt='('' BROYDEN"S 1ST METHOD USED '')')
    if (imix==4) write (ipf, fmt='('' BROYDEN"S 2ND METHOD USED '')')
    if (imix==5) write (ipf, fmt='('' GENERALIZED ANDERSON METHOD USED '')')
    write (ipf, '(A,i4)') ' Iteration index (read in):', mit

    rmixiv = one/alpha

    ! ---->  the following block is activated only one iteration before
    ! broyden iteration scheme is used
    ! ---->  set up of : sm1 = rho(1) ; fm1=fm[1]=f(rho(1)) - rho(1) ;
    ! metric  g := r*r*drdi
    ! ---->  map data of all muffin-tin spheres into one single vector

    call brysh3(sm1, visp, vins, vspsme, ins, irmin, irc, natps, natyp, nspin, imap, lmpot, lsmear)
    call brysh1(fm1, v, vspsmo, ins, irmin, irc, natps, natyp, nspin, imap, lmpot, lsmear)

    if (imap>ntird) call rcstop('NIRDBRY ')

    do ij = 1, imap
      fm1(ij) = rmixiv*(fm1(ij)-sm1(ij))
    end do

    ij = 0
    do isp = 1, nspin
      do ia = natps, natyp

        irc1 = irc(ia)
        volinv = 3.0d0/(r(irc1,ia)**3)

        do ir = 1, irc1
          ij = ij + 1
          g(ij) = atwght(ia)*volinv*r(ir, ia)*r(ir, ia)*drdi(ir, ia)
        end do

        ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
        ! Next for SMEARED spherical potential
        if (lsmear>0) then
          do ir = 1, irc1
            ij = ij + 1
            g(ij) = atwght(ia)*volinv*r(ir, ia)*r(ir, ia)*drdi(ir, ia)
          end do
        end if
        ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

        if (ins>=1 .and. lmpot>1) then
          irmin1 = irmin(ia)
          do lm = 2, lmpot
            do ir = irmin1, irc1
              ij = ij + 1
              g(ij) = atwght(ia)*volinv*r(ir, ia)*r(ir, ia)*drdi(ir, ia)
            end do
          end do
        end if

      end do
    end do

    if (mit>1) then
      rewind iobroy + 2
      read (iobroy+2)(sm1(ij), ij=1, imap), (fm1(ij), ij=1, imap)

      ! ----> map rho(m) of all mt-spheres into one single vector
      call brysh3(sm, visp, vins, vspsme, ins, irmin, irc, natps, natyp, nspin, imap, lmpot, lsmear)

      ! ----> map f[m] = f(m) - rho(m) = f(rho(m)) - rho(m) of all mt-spheres
      ! into one single vector
      call brysh1(fm, v, vspsmo, ins, irmin, irc, natps, natyp, nspin, imap, lmpot, lsmear)
      do ij = 1, imap
        fm(ij) = rmixiv*(fm(ij)-sm(ij))
      end do

      ! ----> calculate  sm = rho(m) - rho(m-1)
      ! ----> calculate dfm = f[m] - f[m-1]
      do ij = 1, imap
        sm1(ij) = sm(ij) - sm1(ij)
        fm1(ij) = fm(ij) - fm1(ij)
      end do

      ! ----> loop to generate u[m] = u(ij,mit)
      do ij = 1, imap
        ui3(ij) = alpha*fm1(ij) + sm1(ij)
      end do
      rewind iobroy
      do it = 2, mit - 1
        read (iobroy)(ui2(ij), ij=1, imap), (vi2(ij), ij=1, imap), wit

        am(it) = ddot(imap, fm1, 1, vi2, 1)
        call daxpy(imap, -am(it), ui2, 1, ui3, 1)
      end do

      ! ----> print amj = the importance of the history of ui
      write (ipf, fmt='(5x,'' AMJ , ---> J=2,'',I3,/,(9X,1P,7D10.2))') mit - 1, (am(it), it=2, mit-1)

      ! -------->     b r o y d e n ' s   f i r s t   m e t h o d
      if (imix==3) then

        ! ----> calculate dsmnorm
        smnorm = zero
        do ij = 1, imap
          smnorm = smnorm + sm1(ij)*g(ij)*sm1(ij)
        end do

        ! ----> convolute dsm with the metric g
        do ij = 1, imap
          sm1(ij) = g(ij)*sm1(ij)
        end do

        ! ----> loop to generate v[m] = v(ij,mit)
        do ij = 1, imap
          vi3(ij) = alpha*sm1(ij)
        end do
        rewind iobroy
        do it = 2, mit - 1
          read (iobroy)(ui2(ij), ij=1, imap), (vi2(ij), ij=1, imap), wit

          bm(it) = ddot(imap, sm1, 1, ui2, 1)
          call daxpy(imap, -bm(it), vi2, 1, vi3, 1)
        end do

        ! ----> complete the evaluation of v[m]
        vmdeno = ddot(imap, sm1, 1, ui3, 1) - smnorm

        if (abs(vmdeno)<1d-70) call rcstop('BRY0SN  ')

        call dscal(imap, one/vmdeno, vi3, 1)

        ! ----> print bmj = the importance of the history of vi
        write (ipf, fmt='(5x,'' BMJ , ---> J=2,'',I3,/,(9X,1P,7D10.2))') mit - 1, (bm(it), it=2, mit-1)

      ! -------->     b r o y d e n ' s   s e c o n d    m e t h o d
      else if (imix==4) then

        ! ----> calculate v[m] ; convoluted with the metric g
        do ij = 1, imap
          vi3(ij) = g(ij)*fm1(ij)
        end do

        ! ----> calculate #vm# and normalize v[m]
        vmnorm = ddot(imap, vi3, 1, fm1, 1)
        call dscal(imap, one/vmnorm, vi3, 1)

      ! -------->     g e n e r a l i z e d   a n d e r s o n   m e t h o d
      else if (imix==5) then

        ! ----> calculate v[m] ; convoluted with the metric g
        do ij = 1, imap
          vi3(ij) = g(ij)*fm1(ij)
        end do
        rewind iobroy
        do it = 2, mit - 1
          read (iobroy)(ui2(ij), ij=1, imap), (vi2(ij), ij=1, imap), wit

          call daxpy(imap, -am(it)*wit(it), vi2, 1, vi3, 1)
        end do

        ! ----> complete the evaluation of v[m]
        vmdeno = ddot(imap, fm1, 1, vi3, 1)

        if (abs(vmdeno)<1d-70) call rcstop('BRY1SN  ')

        call dscal(imap, one/vmdeno, vi3, 1)

        ! ----> save wit(mit) for next iteration
        wit(mit) = vmdeno

      end if

      ! ----> write u3(ij) and v3(ij) on disk
      write (iobroy)(ui3(ij), ij=1, imap), (vi3(ij), ij=1, imap), wit

      ! ----> update f[m-1] = f[m]  ; rho(m) = rho(m-1)
      ! WRITE (IPF,FMT='(5X,'' CMM = '',1P,D12.4)') CMM
      do ij = 1, imap
        fm1(ij) = fm(ij)
        sm1(ij) = sm(ij)
      end do

      ! ----> calculate cmm
      cmm = ddot(imap, fm, 1, vi3, 1)

      ! ----> update rho(m+1)
      call daxpy(imap, one-cmm, ui3, 1, sm, 1)

      ! ----> map solution back into each mt-sphere
      call brysh2(sm, v, vspsmo, ins, irmin, irc, natps, natyp, nspin, imap, lmpot, lsmear)

    end if
    mit = mit + 1
    t_inc%mit_bry = mit

    rewind iobroy + 2
    write (iobroy+2)(sm1(ij), ij=1, imap), (fm1(ij), ij=1, imap)
    ! *********************************************************************
    deallocate (am, bm, fm, fm1, g, sm, sm1, vi3, wit, ui2, ui3, vi2)
    ! *********************************************************************
    return

  end subroutine brydbm

end module mod_brydbm
