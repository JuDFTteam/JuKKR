!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_calctmat

  private
  public :: calctmat

contains

  !-------------------------------------------------------------------------------
  !> Summary: Computes single site t-matrix without SOC
  !> Author: 
  !> Category: KKRhost, single-site
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Computes single-site t-matrix starting from potential by calculating the
  !> scattering wavefunctions and from there the t-matrix
  !>
  !> For KREL = 1 (relativistic mode)                                 
  !>      and                                                                
  !>  NPOTD = 2 * NATYPD                                              
  !>  LMMAXD = 2 * (LMAXD+1)^2                                        
  !>  NSPIND = 1                                                      
  !>                                                                  
  !>  LDA+U implementation     Mar. 2002-Dec.2004                     
  !>                           ph.mavropoulos, h. ebert, v. popescu   
  !-------------------------------------------------------------------------------
  !> @notes:                                                           
  !>  average WLDAU for spherical wavefunctions:                      
  !>  The spherical part of the d or f wavefunction is found by adding
  !>  the average interaction potential `WLDAUAV` to the spherical      
  !>  potential. Then the non-spherical parts are found by using only 
  !>  the deviation of `WLDAU` from the average. This speeds up the     
  !>  convergence of the Born series. See also subroutines `regsol`, `pnstmat` and `pnsqns`
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine calctmat(icst, ins, ielast, nsra, ispin, nspin, i1, ez, drdi, rmesh, vins, visp, zat, irmin, ipan, & 
    ircut, cleb, loflm, icleb, iend, solver, soctl, ctl, vtrel, btrel, rmrel, drdirel, r2drdirel, zrel, jwsrel, idoldau, lopt, wldau, lly, deltae)
  
#ifdef CPP_MPI
    use :: mpi
    use :: mod_mympi, only: mpiadapt
    use :: mod_timing, only: timing_start, timing_stop, timings_1a
#endif
    use :: mod_mympi, only: myrank, nranks, master, distribute_work_energies
    use :: mod_runoptions, only: formatted_files, print_tmat
    use :: mod_types, only: t_tgmat, t_inc, t_mpi_c_grid, init_tgmat, t_lloyd, init_tlloyd
    use :: mod_datatypes, only: dp
    use :: global_variables, only: iemxd, lmmaxd, irmind, irmd, lmpotd, ncleb, krel, lm2d, lmaxd, ipand, mmaxd, nspind
    use :: mod_pnstmat, only: pnstmat
    use :: mod_cradwf, only: cradwf
    use :: mod_wfmesh, only: wfmesh
    use :: mod_regns, only: zgeinv1
    use :: mod_cmatstr, only: cmatstr
    use :: mod_drvreltmat, only: drvreltmat
    use :: mod_constants, only: cvlight, czero, cone

    implicit none
    ! ..
    ! .. Scalar Arguments ..
    real (kind=dp) :: zat
    integer :: i1, icst, ielast, iend, ins, ipan, ispin, nspin, nsra
    integer :: idoldau, jwsrel, lopt, zrel
    integer :: lly                 ! LLY <> 0 for applying Lloyds formula
    integer :: signde, ideriv
    complex (kind=dp) :: deltae    ! Energy difference for numerical derivative
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: ez(iemxd), tmat0(lmmaxd, lmmaxd), tmatll(lmmaxd, lmmaxd), dtmatll(lmmaxd, lmmaxd) ! LLY t-matrix derivative
    real (kind=dp) :: cleb(ncleb, 2), drdi(irmd), rmesh(irmd), vins(irmind:irmd, lmpotd), visp(irmd)
    character (len=10) :: solver
    real (kind=dp) :: soctl(krel*lmaxd+1)
    real (kind=dp) :: ctl(krel*lmaxd+1)
    integer :: icleb(ncleb, 4), ircut(0:ipand), loflm(lm2d)
    real (kind=dp) :: vtrel(irmd*krel+(1-krel))
    real (kind=dp) :: btrel(irmd*krel+(1-krel))
    real (kind=dp) :: drdirel(irmd*krel+(1-krel)), r2drdirel(irmd*krel+(1-krel)), rmrel(irmd*krel+(1-krel))
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: eryd, ek
    integer :: ie, irec, lm1, lm2, lmhi, lmlo, m1, mmax, irmin
    real (kind=dp) :: wldauav
    ! .. this routine does not need irregular wavefunctions
    logical, parameter :: lirrsol = .true.
    ! ..
    ! .. Local Arrays ..
    complex (kind=dp) :: alpha(0:lmaxd), fz(irmd, 0:lmaxd), pns(lmmaxd, lmmaxd, irmind:irmd, 2), pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd), tmat(0:lmaxd), &
      alphall(lmmaxd, lmmaxd), dalphall(lmmaxd, lmmaxd), & ! LLY alpha matrix and derivative
      alpha0(lmmaxd, lmmaxd), aux(lmmaxd, lmmaxd), & ! LLY alpha-matrix
      tralpha, tralpha1            ! LLY Trace of alpha^-1 * d alpha/dE for Lloyds formula
    real (kind=dp) :: wldau(mmaxd, mmaxd, nspind)
    real (kind=dp) :: rs(irmd, 0:lmaxd), sl(0:lmaxd)
    real (kind=dp) :: cutoff(irmd)
    integer :: ipiv(lmmaxd)        ! LLY
    character (len=9) :: txts(2)
#ifdef CPP_MPI
    integer :: ntot_pt(0:nranks-1), ioff_pt(0:nranks-1)
#endif
    integer :: ie_end, ie_num, ie_start, i11
    ! ..
    ! .. External Functions ..
    ! ..
    ! .. Data Statements
    data txts/'spin   UP', 'spin DOWN'/
    ! ..................................................................

    ! save timing information to array, needed to tackle load imbalance adaptively

    ! ==LDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAU

    if (idoldau==1) then
      wldauav = 0.d0
      lmlo = lopt*lopt + 1
      lmhi = (lopt+1)*(lopt+1)
      mmax = lmhi - lmlo + 1
      do m1 = 1, mmax
        wldauav = wldauav + wldau(m1, m1, ispin)
      end do
      wldauav = wldauav/dble(mmax)

      ! -> Note: Application if WLDAU makes the potential discontinuous.
      ! A cutoff can be used if necessary to make the potential continuous
      ! for example (array bounds should be adjusted):

      ! ccc            CUTOFF(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) *
      ! ccc     &                   ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
      ! ccc            CUTOFF(IR) = 1D0/CUTOFF(IR)

      do m1 = 1, irmd
        cutoff(m1) = 1.d0
      end do
    end if

    ! ==LDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAULDAU

    ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    if (myrank==master .and. t_inc%i_write>0) write (1337, *) 'atom: ', i1

    call distribute_work_energies(ielast, distribute_rest=.false.)
#ifdef CPP_MPI
    ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
    ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
    ie_start = 0
    ie_end = ielast
#endif

    ! now initialize arrays for tmat, gmat, and gref
    call init_tgmat(t_inc, t_tgmat, t_mpi_c_grid)
    if (lly/=0) call init_tlloyd(t_inc, t_lloyd, t_mpi_c_grid)


    do ie_num = 1, ie_end

      ie = ie_num + ie_start

#ifdef CPP_MPI
      ! start timing measurement for this pair of ie and i1, needed for MPIadapt
      if (mpiadapt>0) call timing_start('time_1a_ieiatom')
#endif

      if (t_inc%i_write>0) write (1337, *) 'CALCTMAT: IE=', ie, ' ATOM:', i1

      tmatll(:, :) = czero
      dtmatll(:, :) = czero

      alphall(:, :) = czero
      dalphall(:, :) = czero

      ! In case of Lloyds formula the derivative of t is needed.
      ! Then calculate t at E+dE, E-dE and average for t, subtract for dt/dE
      ideriv = 0                   ! LLY
      if (lly/=0) ideriv = 1       ! LLY
      do signde = -ideriv, ideriv, 2 ! LLY
        pz(:, :) = czero
        qz(:, :) = czero
        fz(:, :) = czero
        sz(:, :) = czero
        pns(:, :, irmind:irmd, :) = czero
        tmat0(:, :) = czero
        alpha0(:, :) = czero

        eryd = ez(ie) + signde*deltae/2.d0 ! LLY
        if (t_inc%i_write>0) write (1337, *) 'energy:', ie, '', eryd

        ! =======================================================================
        ! non/scalar-relativistic OR relativistic

        if (krel==0) then
          call wfmesh(eryd, ek, cvlight, nsra, zat, rmesh, sl, rs, ircut(ipan), irmd, lmaxd)
          call cradwf(eryd, ek, nsra, alpha, ipan, ircut, cvlight, rs, sl, pz, fz, qz, sz, tmat, visp, drdi, rmesh, zat, lirrsol, idoldau, lopt, wldauav, cutoff)

          do lm1 = 1, lmmaxd       ! LLY
            alpha0(lm1, lm1) = alpha(loflm(lm1)) ! LLY spherical (diag.) part of alpha
          end do                   ! LLY
          ! -----------------------------------------------------------------------
          ! spherical/non-spherical
          if (ins==0) then
            do lm1 = 1, lmmaxd
              tmat0(lm1, lm1) = tmat(loflm(lm1))
            end do
          else
            call pnstmat(drdi, ek, icst, pz, qz, fz, sz, pns, tmat0, vins, irmin, ipan, ircut, nsra, cleb, icleb, iend, loflm, & ! Added IRMIN 1.7.2014  &
              tmat, lmaxd, idoldau, lopt, lmlo, lmhi, wldau(1,1,ispin), wldauav, cutoff, alpha0) ! LLY In goes diag. alpha, out comes full alpha
          end if
          ! -----------------------------------------------------------------------
        else
          call drvreltmat(eryd, tmat0, vtrel, btrel, rmrel, drdirel, r2drdirel, zrel, jwsrel, solver, soctl, ctl, lmmaxd, lmaxd, irmd)
        end if

        ! non/scalar-relativistic OR relativistic
        ! =======================================================================

        ! In case of derivative calculation, average the t-matrix
        ! of the two energies and calculate the difference
        tmatll(:, :) = tmatll(:, :) + tmat0(:, :) ! LLY
        alphall(:, :) = alphall(:, :) + alpha0(:, :) ! LLY
        if (lly/=0) then
          dtmatll(:, :) = dtmatll(:, :) + signde*tmat0(:, :) ! LLY
          dalphall(:, :) = dalphall(:, :) + signde*alpha0(:, :) ! LLY
        end if

      end do                       ! SIGNDE = -IDERIV,IDERIV,2      ! LLY

      ! Average values of t-matrix and alpha at E+DE and E-DE
      tmatll(:, :) = tmatll(:, :)/real(1+ideriv, kind=dp) ! / 1 or 2 (IDERIV=1 or 2)! LLY
      alphall(:, :) = alphall(:, :)/real(1+ideriv, kind=dp) ! / 1 or 2 ! LLY

      if (lly/=0) then
        ! Construct derivative of t-matrix and alpha
        dtmatll(:, :) = dtmatll(:, :)/deltae
        dalphall(:, :) = dalphall(:, :)/deltae
        ! LLY Calculate Trace [alpha^-1 * d alpha/dE] for Lloyd's formula
        alpha0(:, :) = czero
        aux(:, :) = czero
        ! ALPHA0 = ALPHALL^-1
        call zgeinv1(alphall, alpha0, aux, ipiv, lmmaxd)
        ! AUX = ALPHALL^-1 * DALPHALL
        call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, alpha0, lmmaxd, dalphall, lmmaxd, czero, aux, lmmaxd)
        ! Trace of AUX is Trace of [alpha^-1 * d alpha/dE]
        tralpha = czero
        do lm1 = 1, lmmaxd
          tralpha = tralpha + aux(lm1, lm1)
        end do

        if (lly==-1) then          ! Calculate Tr(t^-1 dt/dE) instead of Tr(alpha^-1 * d alpha/de)
          alpha0(:, :) = tmatll(:, :)
          tmat0(:, :) = czero
          aux(:, :) = czero
          ! ALPHA0 = TMATLL^-1
          call zgeinv1(alpha0, tmat0, aux, ipiv, lmmaxd)
          ! Build trace
          tralpha1 = czero
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              tralpha1 = tralpha1 + tmat0(lm1, lm2)*dtmatll(lm2, lm1)
            end do
          end do
          tralpha1 = czero
          do lm1 = 1, lmmaxd
            tralpha1 = tralpha1 + dtmatll(lm1, lm1)/tmatll(lm1, lm1)
          end do

        end if

      end if                       ! (LLY.NE.0)


      tmat0(:, :) = tmatll(:, :)
      if (t_tgmat%tmat_to_file) then
        irec = ie + ielast*(ispin-1) + ielast*nspin*(i1-1)
        write (69, rec=irec) tmat0
        ! human readable writeout if test option is hit
        if (formatted_files) then
          write (696969, '(i9,20000F15.7)') irec, tmatll(:, :)
        end if
      else
#ifdef CPP_MPI
        i11 = i1-t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie)
#else
        i11 = i1
#endif
        irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(i11-1)
        t_tgmat%tmat(:, :, irec) = tmat0
      end if
      if (lly/=0) then             ! LLY
        if (t_lloyd%dtmat_to_file) then
          irec = ie + ielast*(ispin-1) + ielast*nspin*(i1-1)
          write (691, rec=irec) dtmatll ! LLY
        else
          irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(i1-1)

          t_lloyd%dtmat(:, :, irec) = dtmatll(:, :)
        end if                     ! t_lloyd%dtmat_to_file
        if (t_lloyd%tralpha_to_file) then
          irec = ie + ielast*(ispin-1) + ielast*nspin*(i1-1)
          write (692, rec=irec) tralpha ! LLY
        else
          irec = ie_num + ie_end*(ispin-1) + ie_end*nspin*(i1-1)
          t_lloyd%tralpha(irec) = tralpha
        end if
      end if                       ! LLY

      ! ----------------------------------------------------------------------
      if (print_tmat .and. (t_inc%i_write>0)) then
        write (1337, *)
        write (1337, 100, advance='no') '-----> t matrix for atom: ', i1
        if (krel==0) write (1337, 110, advance='no') txts(ispin)
        write (1337, 120) ', energy: ', eryd
        call cmatstr(' ', 1, tmat0, lmmaxd, lmmaxd, 2*krel+1, 2*krel+1, 0, 1.0e-8_dp, 6)
        write (1337, *)
      end if
      ! ----------------------------------------------------------------------

#ifdef CPP_MPI
      ! stop timing measurement for this pair of ie and i1, needed for MPIadapt
      if (mpiadapt>0) call timing_stop('time_1a_ieiatom', save_out=timings_1a(ie,i1))
#endif

    end do                         ! IE = 1,IELAST
    ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

100 format (a, i3)
110 format (', ', a)
120 format (a, 2f10.6)

    return
  end subroutine calctmat

end module mod_calctmat
