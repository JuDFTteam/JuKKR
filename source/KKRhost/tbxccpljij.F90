!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculates the site-off diagonal XC-coupling parameters \(J_{ij}\)
!> Author: 
!> Calculates the site-off diagonal XC-coupling parameters \(J_{ij}\). According to
!> Lichtenstein et al. JMMM 67, 65 (1987)
!------------------------------------------------------------------------------------
module mod_tbxccpljij

  use :: mod_datatypes

  implicit none

  private
  public :: tbxccpljij

contains


  !-------------------------------------------------------------------------------
  !> Summary: Calculates the site-off diagonal XC-coupling parameters \(J_{ij}\)
  !> Author: 
  !> Category: physical-observables, KKRhost
  !> Deprecated: False 
  !> Calculates the site-off diagonal XC-coupling parameters \(J_{ij}\). According to
  !> Lichtenstein et al. JMMM 67, 65 (1987)
  !-------------------------------------------------------------------------------
  !> @note 
  !> 
  !> * Adopted for TB-KKR code from Munich SPR-KKR package Sep 2004
  !> * For mpi-parallel version: moved energy loop from main1b into here. B. Zimmermann, Dez 2015
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine tbxccpljij(iftmat, ielast, ez, wez, nspin, ncpa, naez, natyp, noq, itoq, iqat, nshell, natomimp, atomimp, ratom, nofgij, nqcalc, iqcalc, ijtabcalc, ijtabsym, ijtabsh, &
    ish, jsh, dsymll, iprint, natypd, nsheld, lmmaxd, npol)

#ifdef CPP_MPI
    use :: mpi
#endif
    use :: mod_types, only: t_tgmat, t_mpi_c_grid, t_cpa
    use :: mod_runoptions, only: calc_exchange_couplings, calc_exchange_couplings_energy, disable_print_serialnumber
    use :: mod_mympi, only: myrank, master
    use :: mod_version_info, only: version_print_header
    use :: mod_md5sums
    use :: mod_cinit
    use :: mod_cmatmul
    use :: mod_initabjij
    use :: mod_cmatstr
    use :: mod_rotatespinframe, only: rotatematrix
    use :: mod_constants, only: pi, cone, czero, nsymaxd
    use :: mod_profiling, only: memocc

    implicit none
    ! .
    ! . Parameters
    integer :: nijmax
    parameter (nijmax=200)
    ! .
    ! . Scalar arguments
    integer :: ielast, nspin, iftmat, iprint, lmmaxd, naez, natomimp, natyp, natypd
    integer :: ncpa, nofgij, nqcalc, nsheld, npol
    complex (kind=dp) :: ez(*), wez(*)
    ! .
    ! . Array arguments
    integer :: atomimp(*), ijtabcalc(*), ijtabsh(*), ijtabsym(*), iqat(*), iqcalc(*)
    integer :: ish(nsheld, *), itoq(natypd, 2*nsymaxd), jsh(nsheld, 2*nsymaxd), noq(*), nshell(0:nsheld)
    complex (kind=dp) :: dsymll(lmmaxd, lmmaxd, *)
    real (kind=dp) :: ratom(3, *)
    ! .
    ! . Local scalars
    integer :: i1, ia, ifgmat, ifmcpa, iq, irec, ispin, isym, it, j1, ja, jq, jt, l1, i_all
    integer :: lm1, lm2, lstr, ns, nseff, nshcalc, nsmax, ntcalc, ie, ie_start, ie_end, ie_num
#ifdef CPP_MPI
    integer :: ierr
#endif
    complex (kind=dp) :: csum, wgtemp
#ifdef CPP_MPI
    complex (kind=dp) :: xintegdtmp
#endif
    character (len=8) :: fmt1
    character (len=22) :: fmt2
    character (len=8) :: jfbas
    character (len=13) :: jfnam
    character (len=26) :: jfnam2
    character (len=80) :: strbar, strtmp
    ! .
    ! . Local arrays
    integer, allocatable :: nijcalc(:), kijsh(:, :), jijdone(:, :, :)
    complex (kind=dp), allocatable :: jxcijint(:, :, :)
#ifndef CPP_MPI
    complex (kind=dp), allocatable :: xintegd(:, :, :)
#else
    complex (kind=dp), allocatable :: csum_store(:, :, :, :), csum_store2(:, :, :, :)
#endif
    complex (kind=dp) :: deltsst(lmmaxd, lmmaxd, natyp) 
    complex (kind=dp) :: dmatts(lmmaxd, lmmaxd, natyp, nspin)
    complex (kind=dp) :: dtilts(lmmaxd, lmmaxd, natyp, nspin), gmij(lmmaxd, lmmaxd)
    complex (kind=dp) :: gmji(lmmaxd, lmmaxd)
    complex (kind=dp) :: gs(lmmaxd, lmmaxd, nspin), tsst(lmmaxd, lmmaxd, natyp, 2)
    complex (kind=dp) :: w1(lmmaxd, lmmaxd), w2(lmmaxd, lmmaxd), w3(lmmaxd, lmmaxd)
    real (kind=dp) :: rsh(nsheld)
    integer :: jtaux(natyp)
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: max, sqrt
    ! ..
    ! .. External Subroutines ..
    ! ..
    ! .. Save statement
    save :: ifgmat, ifmcpa, jijdone, jxcijint, kijsh, nijcalc
#ifndef CPP_MPI
    save :: xintegd
#endif
    save :: nshcalc, nsmax
    ! ..
    ! .. Data statement
    data jfbas/'Jij.atom'/
    ! ..
    ! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    ! ==>                   -- initialisation step --                    <==
    ! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
    ! open(22,STATUS='unknown',FILE='integrand1.dat',
    ! &                         FORM='formatted')
    ! open(44,STATUS='unknown',FILE='integrand2.dat',
    ! &                         FORM='formatted')
    ! write(*,*) 'test brahim 2'

    ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    if (iprint>0) write (1337, 100)
#ifdef CPP_MPI
    ! this part of the code does not take much time, thus only the energy
    ! loop is parallel and not the second level over atoms, only the master
    ! in every energy loop does the calculation of the Jijs
    if (t_mpi_c_grid%myrank_ie==0) then

      ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
      ie_end = ielast
#endif
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

      ! open(22,STATUS='unknown',FILE='integrand.dat',
      ! &                         FORM='formatted')
      ifgmat = iftmat + 1
      ifmcpa = ifgmat + 1
      nsmax = max(naez, natyp)
      nshcalc = nshell(0) - nsmax
      ! ccc         IF ( NSHLOC.LT.NSHELD) STOP ' NSHLOC'
      ! ccc         IF ( NTLOC.LT.NATYP) STOP ' NTLOC'

      allocate (kijsh(nijmax,nshell(0)), stat=lm1)
      call memocc(lm1, product(shape(kijsh))*kind(kijsh), 'kijsh', 'tbxccpljij')
      if (lm1/=0) then
        write (6, 110) 'KIJSH'
        stop
      end if

      allocate (nijcalc(nshell(0)), jijdone(natyp,natyp,nshell(0)), stat=lm1)
      call memocc(lm1, product(shape(nijcalc))*kind(nijcalc), 'nijcalc', 'tbxccpljij')
      if (lm1/=0) then
        write (6, 110) 'JIJDONE/NIJCALC'
        stop
      end if

      ! The array CSUM_STORE could become quite large -> store to MPI-IO-file in future????
      ! Only allocate it for MPI-usage. If MPI is not used, two smaller array XINTEGD arrays is needed.
#ifdef CPP_MPI
      allocate (csum_store(natyp,natyp,nshell(0),ielast), stat=lm1)
      call memocc(lm1, product(shape(csum_store))*kind(csum_store), 'csum_store', 'tbxccpljij')
      if (lm1/=0) then
        write (6, 110) 'csum_store'
        stop
      end if
#else
      ! ALLOCATE (XINTEGD1(NATYP,NATYP,IELAST),STAT=LM1)
      allocate (xintegd(natyp,natyp,nshell(0)), stat=lm1)
      call memocc(lm1, product(shape(xintegd))*kind(xintegd), 'xintegd', 'tbxccpljij')
      if (lm1/=0) then
        write (6, 110) 'XINTEGD'
        stop
      end if
#endif
      allocate (jxcijint(natyp,natyp,nshell(0)), stat=lm1)
      call memocc(lm1, product(shape(jxcijint))*kind(jxcijint), 'jxcijint', 'tbxccpljij')
      if (lm1/=0) then
        write (6, 110) 'JXCIJINT'
        stop
      end if
      ! ..................
#ifdef CPP_MPI
      do ie = 1, ielast
        do ns = 1, nshell(0)
          do jt = 1, natyp
            do it = 1, natyp
              csum_store(it, jt, ns, ie) = czero
            end do
          end do
        end do
      end do                       ! IE
#else
      do ns = 1, nshell(0)
        do jt = 1, natyp
          do it = 1, natyp
            xintegd(it, jt, ns) = czero
            ! XINTEGD1(IT,JT,IE)= CZERO
          end do
        end do
      end do
#endif
      do ns = 1, nshell(0)
        nijcalc(ns) = 0
        do jt = 1, natyp
          do it = 1, natyp
            jijdone(it, jt, ns) = 0
            jxcijint(it, jt, ns) = czero
          end do
        end do
      end do

      call initabjij(iprint, naez, natyp, natomimp, nofgij, nqcalc, nsmax, nshell, iqcalc, atomimp, ish, jsh, ijtabcalc, ijtabsh, ijtabsym, nijcalc(1), kijsh(1,1), nijmax, &
        nshell(0), nsheld)
      ! --------------------------------------------------------------------
      do ns = nsmax + 1, nshell(0)
        nseff = ns - nsmax
        do i1 = 1, nijcalc(ns)
          ia = ish(ns, kijsh(i1,ns))
          ja = jsh(ns, kijsh(i1,ns))
          lm1 = (ia-1)*natomimp + ja
          iq = atomimp(ia)
          jq = atomimp(ja)

          do l1 = 1, noq(jq)
            jt = itoq(l1, jq)
            do j1 = 1, noq(iq)
              it = itoq(j1, iq)
              jijdone(it, jt, nseff) = 1
            end do
          end do
        end do
      end do
      ! ----------------------------------------------------------------------

      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      if (iprint>0) then
        do it = 1, natyp
          do jt = 1, natyp
            jtaux(jt) = 0
            do ns = nsmax + 1, nshell(0)
              nseff = ns - nsmax
              if (jijdone(it,jt,nseff)/=0) jtaux(jt) = jtaux(jt) + 1
            end do
          end do
          ! ccc               WRITE (6,99012) IT,(JTAUX(JT),JT=1,MIN(25,NATYP))
        end do
        write (1337, 230)
      end if
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

      ! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
      ! ==>                   INITIALISATION END                           <==
      ! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII


      ! ==>  read in the single site matrices (TSST) in the LOCAL frame

      ! ==>  set up the Delta t_i matrices DELTSST = TSST(UP) - TSST(DOWN)

      ! ==>  read in projection matrices DMAT and DTIL used to get
      ! ij            ij    _
      ! G     =  D  * G    * D
      ! ab       a    CPA    b
      ! with a/b the atom of type a/b sitting on site i/j
      ! for an atom having occupancy 1, DMAT/DTIL = unit matrix


      ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES
#ifdef CPP_MPI
      ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
      ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
      ie_start = 0
      ie_end = ielast
#endif

      do ie_num = 1, ie_end
        ie = ie_start + ie_num
        ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS SPIN
        do ispin = 1, nspin
          ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT NATYP
          do it = 1, natyp
            if (t_tgmat%tmat_to_file) then
              irec = ie + ielast*(ispin-1) + ielast*nspin*(it-1)
              read (iftmat, rec=irec) w1
            else
              irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(it-1)
              w1(:, :) = t_tgmat%tmat(:, :, irec)
            end if
            do j1 = 1, lmmaxd
              call zcopy(lmmaxd, w1(1,j1), 1, tsst(1,j1,it,ispin), 1)
            end do
            ! ----------------------------------------------------------------------
            if (ispin==2) then
              do j1 = 1, lmmaxd
                do i1 = 1, lmmaxd
                  deltsst(i1, j1, it) = tsst(i1, j1, it, 2) - tsst(i1, j1, it, 1)
                end do
              end do
            end if
            ! ----------------------------------------------------------------------
            if (ncpa==0) then
              call cinit(lmmaxd*lmmaxd, dmatts(1,1,it,ispin))
              do i1 = 1, lmmaxd
                dmatts(i1, i1, it, ispin) = cone
              end do
              call zcopy(lmmaxd*lmmaxd, dmatts(1,1,it,ispin), 1, dtilts(1,1,it,ispin), 1)
            else                 ! NCPA.EQ.0
              if (t_cpa%dmatproj_to_file) then
                write (*, *) 'test read proj'
                irec = ie + ielast*(ispin-1) + ielast*2*(it-1)
                read (ifmcpa, rec=irec) w1, w2
                do j1 = 1, lmmaxd
                  do i1 = 1, lmmaxd
                    dmatts(i1, j1, it, ispin) = w1(i1, j1)
                    dtilts(i1, j1, it, ispin) = w2(i1, j1)
                  end do
                end do
              else               ! t_cpa%dmatproj_to_file
                irec = ie_num + ie_end*(ispin-1)
                dmatts(:, :, it, ispin) = t_cpa%dmatts(:, :, it, irec)
                dtilts(:, :, it, ispin) = t_cpa%dtilts(:, :, it, irec)

              end if             ! t_cpa%dmatproj_to_file
            end if               ! NCPA.EQ.0
            ! ----------------------------------------------------------------------
          end do
          ! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
        end do
        ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        if (iprint>1) then
          write (1337,'(8X,60("-"),/,10X,  &
              "Single-site and projection matrices read in", " for IT=1,",I3,/)') natyp
          do ispin = 1, nspin
            write (1337, '(8X,60("+"),/,30X, " ISPIN = ",I1,/,8X,60("+"),/)') ispin
            do it = 1, natyp
              write (1337, '(12X," IE = ",I2," IT =",I3)') ie, it
              call cmatstr(' T MAT ', 7, tsst(1,1,it,ispin), lmmaxd, lmmaxd, 0, 0, 0, 1.0e-8_dp, 6)
              call cmatstr(' D MAT ', 7, dmatts(1,1,it,ispin), lmmaxd, lmmaxd, 0, 0, 0, 1.0e-8_dp, 6)
              call cmatstr(' D~ MAT', 7, dtilts(1,1,it,ispin), lmmaxd, lmmaxd, 0, 0, 0, 1.0e-8_dp, 6)
              if (it/=natyp) write (1337, '(8X,60("-"),/)')
            end do
            write (1337, '(8X,60("+"),/)')
          end do
          write (1337,'(8X,60("-"),/,10X,  &
              "Delta_t = t(it,DN) - t(it,UP) matrices for IT=1,", I3,/)') natyp
          do it = 1, natyp
            write (1337, '(12X," IE = ",I2," IT =",I3)') ie, it
            call cmatstr(' DEL T ', 7, deltsst(1,1,it), lmmaxd, lmmaxd, 0, 0, 0, 1.0e-8_dp, 6)
            if (it/=natyp) write (1337, '(8X,60("-"),/)')
          end do
          write (1337, '(8X,60("-"),/)')
        end if
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

        ! ***************************************************** loop over shells
        do ns = nsmax + 1, nshell(0)

          ! ----------------------------------------------------------------------

          ! ==>  get the off-diagonal Green function matrix Gij(UP) and Gji(DOWN)
          ! using the appropiate rotation (similar to ROTGLL routine)
          ! step one: read in GS for the representative shell
          ! ----------------------------------------------------------------------
          do ispin = 1, nspin
            if (t_tgmat%gmat_to_file) then
              irec = ie + ielast*(ispin-1) + ielast*nspin*(ns-1)
              read (ifgmat, rec=irec) w1
            else
              ie_end = t_mpi_c_grid%ntot2
              irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(ns-1)
              w1(:, :) = t_tgmat%gmat(:, :, irec)
            end if
            call zcopy(lmmaxd*lmmaxd, w1, 1, gs(1,1,ispin), 1)
          end do
          ! ----------------------------------------------------------------------
          ! step two: scan all I,J pairs needed out of this shell
          ! transform with the appropriate symmetry to get Gij

          ! ------------------------------------------------------ loop over pairs
          nseff = ns - nsmax
          do l1 = 1, nijcalc(ns)

            ia = ish(ns, kijsh(l1,ns))
            ja = jsh(ns, kijsh(l1,ns))

            lm1 = (ia-1)*natomimp + ja
            isym = ijtabsym(lm1)
            ! ----------------------------------------------------------------------
            do ispin = 1, nspin
              call zgemm('C', 'N', lmmaxd, lmmaxd, lmmaxd, cone, dsymll(1,1,isym), lmmaxd, gs(1,1,ispin), lmmaxd, czero, w2, lmmaxd)

              call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, w2, lmmaxd, dsymll(1,1,isym), lmmaxd, czero, w1, lmmaxd)

              if (ispin==1) then
                call zcopy(lmmaxd*lmmaxd, w1, 1, gmij, 1)
              else
                do lm2 = 1, lmmaxd
                  do lm1 = 1, lmmaxd

                    ! -> use Gji = Gij^T

                    gmji(lm1, lm2) = w1(lm2, lm1)
                  end do
                end do
              end if
            end do
            ! ----------------------------------------------------------------------

            ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
            if (iprint>2) then
              write (1337,'(8X,60("-"),/,10X,  &
                  " G_ij(DN) and G_ji(UP) matrices I =",I3," J =",I3,  &
                  " IE =",I3,/,8X,60("-"))') ia,ja,ie
              call cmatstr(' Gij DN', 7, gmij, lmmaxd, lmmaxd, 0, 0, 0, 1.0e-8_dp, 6)
              call cmatstr(' Gji UP', 7, gmji, lmmaxd, lmmaxd, 0, 0, 0, 1.0e-8_dp, 6)
            end if
            ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO


            ! ----------------------------------------------------------------------

            ! ==> calculate the exchange coupling constant J_ij via Eq. (19)
            ! modified for G instead of tau:
            ! J_ij ~ Trace [ (t_i(D)-t_i(U)) * Gij(U)
            ! * (t_j(D)-t_j(U)) * Gji(D)]
            ! in case of alloy system: perform projection on atom types

            ! ----------------------------------------------------------------------
            iq = atomimp(ia)
            jq = atomimp(ja)
            ! -------------------------------------------------- loop over occupants
            do j1 = 1, noq(jq)
              jt = itoq(j1, jq)
              do i1 = 1, noq(iq)
                it = itoq(i1, iq)

                ! --> Gjq,iq is projected on jt,it ==> Gjt,it

                call cmatmul(lmmaxd, lmmaxd, gmji, dtilts(1,1,it,2), w2)
                call cmatmul(lmmaxd, lmmaxd, dmatts(1,1,jt,2), w2, w1)

                ! --> Delta_j * Gjt,it

                call cmatmul(lmmaxd, lmmaxd, deltsst(1,1,jt), w1, w2)

                ! --> Giq,jq is projected on it,jt ==> Git,jt * Delta_j * Gjt,it

                call cmatmul(lmmaxd, lmmaxd, gmij, dtilts(1,1,jt,1), w3)
                call cmatmul(lmmaxd, lmmaxd, dmatts(1,1,it,1), w3, w1)

                ! --> Delta_i * Git,jt

                call cmatmul(lmmaxd, lmmaxd, deltsst(1,1,it), w1, w3)

                ! --> Delta_i * Git,jt * Delta_j * Gjt,it

                call cmatmul(lmmaxd, lmmaxd, w3, w2, w1)

                csum = czero
                do lm1 = 1, lmmaxd
                  csum = csum + w1(lm1, lm1)
                end do

#ifdef CPP_MPI
                ! BZ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! BZ! to ensure the correct order of energy points, the csum=values !!
                ! BZ! are stored for the MPI version, and the evaluation of JXCIJINT!!
                ! BZ! and XINTEGD is performed later                                !!
                ! BZ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                csum_store(it, jt, nseff, ie) = csum
#else

                jxcijint(it, jt, nseff) = jxcijint(it, jt, nseff) - wez(ie)*csum/real(nspin, kind=dp)

                xintegd(it, jt, nseff) = csum/(pi*4.d0)

                ! -------> perform substraction instead of addition
                ! because WGTE ~ -1/pi (WGTE = WEZ(IE)/NSPIN)
                ! Write out energy-resorved integrand and integral
                ! Phivos Mavropoulos 24.10.2012
                if (npol==0 .or. calc_exchange_couplings_energy) then
                  fmt2 = '(A,I5.5,A,I5.5,A,I5.5)'
                  write (jfnam2, fmt2) 'Jij_enrg.', it, '.', jt, '.', ns
                  if (ie==1) then
                    open (499, file=jfnam2, status='UNKNOWN')
                    call version_print_header(499, addition='; '//md5sum_potential//'; '//md5sum_shapefun, disable_print=disable_print_serialnumber)
                    write (499, fmt='(A)') '# Energy Re,Im ; j(E) Re,Im; J(E) Re,Im '
                    write (499, fmt='(3(A,I5))') '# IT=', it, ' JT=', jt, ' SHELL=', ns
                    write (499, fmt='(A,I6)') '#ENERGIES: ', ielast
                  else
                    open (499, file=jfnam2, status='OLD', position='APPEND')
                  end if
                  write (499, fmt='(6E12.4)') ez(ie), xintegd(it, jt, nseff), jxcijint(it, jt, nseff)/4.d0
                  close (499)
                end if           ! (npol==0 .or. calc_exchange_couplings_energy)
#endif

              end do             ! I1
            end do               ! J1, loop over occupants
            ! ----------------------------------------------------------------------
          end do                 ! L1 = 1,NIJCALC(NS)

          ! ----------------------------------------------------------------------
        end do                   ! loop over shells
        ! **********************************************************************
        ! write(22,*)DIMAG(XINTEGD(1,1,1)),DIMAG(JXCIJINT(1,1,1))
        ! write(44,*)DIMAG(XINTEGD(2,2,1)),DIMAG(XINTEGD(2,3,1))

      end do                     ! IE
      ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES


      ! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC MPI Communication
#ifdef CPP_MPI
      ! allocate
      allocate (csum_store2(natyp,natyp,nshell(0),ielast), stat=lm1)
      call memocc(lm1, product(shape(csum_store2))*kind(csum_store2), 'csum_store2', 'tbxccpljij')
      if (lm1/=0) then
        write (6, 110) 'csum_store2'
        stop
      end if

      ! initialize with zero
      do ie = 1, ielast
        do ns = 1, nshell(0)
          do jt = 1, natyp
            do it = 1, natyp
              csum_store2(it, jt, ns, ie) = czero
            end do
          end do
        end do
      end do                     ! IE

      ! perform communication and collect resutls in csum_store2
      lm1 = natyp*natyp*nshell(0)*ielast
      call mpi_reduce(csum_store, csum_store2, lm1, mpi_double_complex, mpi_sum, master, t_mpi_c_grid%mympi_comm_at, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Problem in MPI_REDUCE(csum_store)'
        stop 'TBXCCPLJIJ'
      end if                     ! IERR/=MPI_SUCCESS

      ! copy array back to original one
      if (myrank==master) then
        do ie = 1, ielast
          do ns = 1, nshell(0)
            do jt = 1, natyp
              do it = 1, natyp
                csum_store(it, jt, ns, ie) = csum_store2(it, jt, ns, ie)
              end do
            end do
          end do
        end do                   ! IE

        i_all = -product(shape(csum_store2))*kind(csum_store2)
        deallocate (csum_store2, stat=lm1)
        call memocc(lm1, i_all, 'csum_store2', 'tbxccpljij')
      end if                     ! myrank==master



      ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      ! ~~~~~~~ NEW WRITEOUT: integrand and energy-resolved integral ~~~~~~~~~~
      ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      if (myrank==master) then
        do ns = nsmax + 1, nshell(0)
          nseff = ns - nsmax
          do l1 = 1, nijcalc(ns)
            ia = ish(ns, kijsh(l1,ns))
            ja = jsh(ns, kijsh(l1,ns))
            iq = atomimp(ia)
            jq = atomimp(ja)
            do j1 = 1, noq(jq)
              jt = itoq(j1, jq)
              do i1 = 1, noq(iq)
                it = itoq(i1, iq)
                ! -------> perform substraction instead of addition
                ! because WGTE ~ -1/pi (WGTE = WEZ(IE)/NSPIN)
                ! Write out energy-resorved integrand and integral
                ! Phivos Mavropoulos 24.10.2012
                if (npol==0 .or. calc_exchange_couplings_energy) then
                  fmt2 = '(A,I5.5,A,I5.5,A,I5.5)'
                  write (jfnam2, fmt2) 'Jij_enrg.', it, '.', jt, '.', ns
                  open (499, file=jfnam2, status='UNKNOWN')
                  call version_print_header(499, addition='; '//md5sum_potential//'; '//md5sum_shapefun, disable_print=disable_print_serialnumber)
                  write (499, fmt='(A)') '# Energy Re,Im ; j(E) Re,Im; J(E) Re,Im '
                  write (499, fmt='(3(A,I5))') '# IT=', it, ' JT=', jt, ' SHELL=', ns
                  write (499, fmt='(A,I6)') '#ENERGIES: ', ielast
                end if           ! (NPOL==0 .OR. calc_exchange_couplings_energy)then

                do ie = 1, ielast
                  jxcijint(it, jt, nseff) = jxcijint(it, jt, nseff) - wez(ie)*csum_store(it, jt, nseff, ie)/real(nspin, kind=dp)
                  xintegdtmp = csum_store(it, jt, nseff, ie)/(pi*4.d0)
                  if (npol==0 .or. calc_exchange_couplings_energy) then
                    write (499, fmt='(6E12.4)') ez(ie), xintegdtmp, jxcijint(it, jt, nseff)/4.d0
                  end if         ! (NPOL==0 .OR. calc_exchange_couplings_energy)then

                end do           ! IE

                if (npol==0 .or. calc_exchange_couplings_energy) close (499)

              end do             ! I1
            end do               ! J1, loop over occupants
          end do                 ! L1 = 1,NIJCALC(NS)
        end do                   ! NS, loop over shells
      end if                     ! myrank==master
#endif

      if (myrank==master) then
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        wgtemp = cone/4d0
        ! -------> here factor 1/pi omitted since it is included in WGTE
        do ns = 1, nshcalc
          do jt = 1, natyp
            do it = 1, natyp
              jxcijint(it, jt, ns) = wgtemp*jxcijint(it, jt, ns)
            end do
          end do
          nseff = ns + nsmax
          rsh(ns) = 0d0
          do i1 = 1, 3
            rsh(ns) = rsh(ns) + ratom(i1, nseff)*ratom(i1, nseff)
          end do
          rsh(ns) = sqrt(rsh(ns))
        end do

        lm2 = 0
        ntcalc = 0
        do it = 1, natyp
          l1 = 0
          do ns = 1, nshcalc
            lm1 = 0
            do jt = 1, natyp
              if (jijdone(it,jt,ns)/=0) lm1 = lm1 + 1
            end do
            lm2 = max(lm1, lm2)
            l1 = max(l1, lm1)
          end do
          if (l1>0) then
            ntcalc = ntcalc + 1
            jtaux(ntcalc) = it
          end if
        end do
        write (strbar, '(19("-"))')
        lstr = 19
        do i1 = 1, lm2
          write (strtmp, '(A,15("-"))') strbar(1:lstr)
          lstr = lstr + 15
          strbar(1:lstr) = strtmp(1:lstr)
        end do

        write (1337, 120) strbar(1:lstr), strbar(1:lstr)
        do i1 = 1, ntcalc
          it = jtaux(i1)
          write (1337, 130, advance='no') it, iqat(it)
          l1 = 0
          do ns = 1, nshcalc
            lm1 = 0
            do jt = 1, natyp
              if (jijdone(it,jt,ns)/=0) lm1 = lm1 + 1
            end do
            if (lm1/=0) then
              lm2 = 0
              if (l1==0) then
                write (1337, 150, advance='no') rsh(ns)
                l1 = 1
              else
                write (1337, 160, advance='no') rsh(ns)
              end if
              do jt = 1, natyp
                if (jijdone(it,jt,ns)/=0) then
                  lm2 = lm2 + 1
                  if (lm2==1) then
                    write (1337, 170, advance='no') iqat(jt), aimag(jxcijint(it,jt,ns))*1d3, jt
                    write (1337, *) ns + nsmax, ' shell'
                  else
                    write (1337, 180, advance='no') aimag(jxcijint(it,jt,ns))*1d3, jt
                    write (1337, *) ns + nsmax, ' shell'
                  end if
                  if (lm2==lm1) write (1337, *)
                end if
              end do
            end if
          end do
          write (1337, 140) strbar(1:lstr)
        end do                   ! I1 = 1,NTCALC
        write (1337, *)
        ! ----------------------------------------------------------------------
        ! --> prepare output files

        do i1 = 1, ntcalc
          it = jtaux(i1)
          do ns = 1, nshcalc
            l1 = 1
            do jt = 1, natyp
              if (jijdone(it,jt,ns)/=0) then
                jijdone(it, l1, ns) = jt
                jxcijint(it, l1, ns) = jxcijint(it, jt, ns)
                l1 = l1 + 1
              end if
            end do
            do jt = l1, natyp
              jijdone(it, jt, ns) = 0
            end do
          end do
        end do

        do i1 = 1, ntcalc
          it = jtaux(i1)
          ! FMT1 = '(A,I1)'
          ! IF ( IT.GE.10 ) THEN
          ! FMT1 = '(A,I2)'
          ! IF ( IT.GE.100 ) FMT1 = '(A,I3)'
          ! endif
          fmt1 = '(A,I5.5)'
          write (jfnam, fmt1) jfbas, it

          ! write(JFNAME,FMT1)JFINTEG, IT
          open (49, file=jfnam)
          call version_print_header(49, addition='; '//md5sum_potential//'; '//md5sum_shapefun, disable_print=disable_print_serialnumber)
          write (49, 190) it, iqat(it)
          ! open(22,FILE=JFNAME)
          ! write(22,99014)IT,IQAT(IT)
          do l1 = 1, natyp
            lm1 = 0
            do ns = 1, nshcalc
              if (jijdone(it,l1,ns)/=0) then
                lm1 = lm1 + 1
                write (49, 200) rsh(ns), aimag(jxcijint(it,l1,ns)), jijdone(it, l1, ns), ns + nsmax ! fivos added NS+NSMAX
                ! do IE=1,IELAST
                ! XINTEGD(IT,JT,IE)= XINTEGD(IT,JT,IE)*(1.D0/(PI*4.D0))

                ! write(22,99015)XINTEGD(IT,JT,IE),JIJDONE(IT,L1,NS)
                ! end do
              end if
            end do
            ! end do
            if ((lm1/=0) .and. (l1/=natyp)) write (49, *) '&'
            ! IF ( (LM1.NE.0) .AND. (L1.NE.NATYP) ) WRITE (22,*) '&'
          end do                 ! l1=1,natyp
          close (49)
          ! Close(22)
        end do                   ! i1=1,ntcalc
        write (1337, 210)
        ! ----------------------------------------------------------------------
        ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
      end if                     ! myrank==master
#ifdef CPP_MPI
    end if                       ! t_mpi_c_grid%myrank_ie==0
#endif

    i_all = -product(shape(kijsh))*kind(kijsh)
    deallocate (kijsh, stat=lm1)
    call memocc(lm1, i_all, 'kijsh', 'tbxccpljij')

    i_all = -product(shape(nijcalc))*kind(nijcalc)
    deallocate (nijcalc, stat=lm1)
    call memocc(lm1, i_all, 'nijcalc', 'tbxccpljij')

    i_all = -product(shape(jijdone))*kind(jijdone)
    deallocate (jijdone, stat=lm1)
    call memocc(lm1, i_all, 'jijdone', 'tbxccpljij')

    i_all = -product(shape(jxcijint))*kind(jxcijint)
    deallocate (jxcijint, stat=lm1)
    call memocc(lm1, i_all, 'jxcijint', 'tbxccpljij')

100 format (79('='), /, 10x, 'TBXCCPLJIJ : Off-diagonal exchange coupling', ' constants J_ij', /, 79('='), /)
110 format (6x, 'ERROR: Cannot allocate array(s) :', a, /)
120 format (4x, a, 4('-'), /, 5x, ' IT/IQ ', 3x, 'R_IQ,JQ', 2x, 'JQ ', ' (J_IT,JT  JT)', /, 15x, ' [ ALAT ] ', 4x, '[ mRy ]', /, 4x, a, 4('-'))
130 format (5x, i3, 1x, i3)
140 format (4x, a, 4('-'))
150 format (f10.6)
160 format (12x, f10.6)
170 format (i4, f12.8, i3)
180 format (f12.8, i3)
190 format ('# off-diagonal exchange coupling constants ', /, '# for atom IT = ', i3, ' on site IQ = ', i3, /, '# R_IQ,JQ      J_IT,JT       JT', /, &
      '# ( ALAT )       ( Ry )', /, '#      ')
200 format (f12.8, 2x, e15.8, 2x, i3, i5)
210 format (6x, 'Output written into the files Jij.atomX', /, 6x, '  X = atom index in the input file')
220 format (10x, i3, 3x, 25(i3))
230 format (/, 8x, 60('-'), /)

240 format ('# Integrand  ', /, '# for atom IT = ', i3, ' on site IQ = ', i3, /, '#      J_IT,JT       JT', /, '#        ( Ry )', /, '#      ')
250 format (f12.8)
  end subroutine tbxccpljij

end module mod_tbxccpljij
