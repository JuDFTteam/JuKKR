module mod_decimaread
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine decimaread(ez, tk, nptp1, nptp2, nptp3, npol, ispin, lefttinvll, righttinvll, vacflag, ienergy, nlbasis, nrbasis, naez, kaoez, kmrot, ins, nspin, lmmax, ielast, &
    fileleft, fileright, krel, natypd, lmmaxd, nembd1)
    ! ,KORBIT)
    ! **********************************************************************
    ! *                                                                    *
    ! * This subroutine reads in the t-matrices of the left                *
    ! * and right host for the decimation method.                          *
    ! *                                                                    *
    ! * The t-matrices are writen out in kloopz1  (option 'deci-out')      *
    ! *                                                                    *
    ! * The host files contain the CMOMHOST data neeed to create the       *
    ! * interatomic potential in subroutine < vinterface >                 *
    ! * This is going to be read in < cmomsread >, the decimation files    *
    ! * are for this reason not rewinded.                                  *
    ! *                                                                    *
    ! * In case of 'vacuum' setting on one of the sides, no energy points  *
    ! * are read in and the VACLAG is set to TRUE.                         *
    ! *                                                                    *
    ! * A call of this routine with IENERGY = 0 means that we are not in   *
    ! * the energy loop - only the header of each decimation file read in  *
    ! *                                                                    *
    ! * IENERGY <> 0 reads in the matrices/energy at IENERGY -> returned   *
    ! *                                                                    *
    ! **********************************************************************
    use :: mod_version_info
    use :: mod_lngstring
    use :: mod_cinit
    implicit none
    ! ..
    integer :: krel, natypd, nembd1, lmmaxd ! ,KORBIT
    ! ..
    ! .. Scalar arguments
    integer :: nptp1, nptp2, nptp3, npol, ispin, ienergy
    integer :: nlbasis, nrbasis, naez, kmrot, ins, nspin, lmmax, ielast
    real (kind=dp) :: tk
    character (len=40) :: fileleft, fileright
    ! ..
    ! .. Array arguments
    integer :: kaoez(natypd, *)
    complex (kind=dp) :: ez(*)
    complex (kind=dp) :: lefttinvll(lmmaxd, lmmaxd, nembd1), righttinvll(lmmaxd, lmmaxd, nembd1)
    logical :: vacflag(2)
    ! ..
    ! .. Local variables
    integer :: ios, npt1l, npt2l, npt3l, npoll, insl
    integer :: nspinl, naezl, lmmaxl, krell, kmrotl
    integer :: iel, ih, lm1, lm2, i, idum, ih1, ihost, ilhost, nathost
    real (kind=dp) :: alatl, templ, e2l, qmt, qmp
    complex (kind=dp) :: el
    complex (kind=dp) :: w1(lmmaxd, lmmaxd)
    character (len=80) :: baner1
    character (len=40) :: filehost
    character (len=5) :: chhost(2), str5
    real (kind=dp) :: bravaisl(3, 3), rbasisl(3, nembd1)
    ! ..
    ! .. Data statements
    data chhost/'LEFT ', 'RIGHT'/
    ! ..
    ! ========================================================== IENERGY = 0
    if (ienergy<1) then
      write (1337, '(5X,A,/,8X,65("-"))') 'Reading in host Delta_t matrices'
      vacflag(1) = .false.
      vacflag(2) = .false.
      ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      do ihost = 1, 2
        filehost = fileleft
        nathost = nlbasis
        if (ihost==2) then
          filehost = fileright
          nathost = nrbasis
        end if
        ilhost = lngstring(filehost, 40)

        write (1337, '(8X,A5," side host: ")', advance='no') chhost(ihost)
        ! ----------------------------------------------------------------------
        if (filehost(1:7)=='vacuum') then
          write (1337, '(A)') 'VACUUM will be used'
          vacflag(ihost) = .true.
          ! ----------------------------------------------------------------------
        else
          write (1337, '(A,/)') filehost(1:ilhost)
          open (36+ihost, file=filehost, status='OLD', iostat=ios)
          call version_check_header(36+ihost)
          ! ......................................................................
          if (ios>0) then
            write (6, '(/,5X,2A)') 'ERROR: Can not open host file ', filehost(1:ilhost)
            write (6, '(12X,A,A)') 'vacuum    entry should be used to set', ' a vacuum-host on this side'
            stop '       < DECIMAREAD > '
          end if
          ! ......................................................................
          do ih = 1, 3
            read (36+ihost, 190) baner1
          end do
          read (36+ihost, 210) alatl, nspinl, naezl, lmmaxl, insl, krell, kmrotl
          ! ......................................................................
          if ((krell/=krel) .or. (kmrotl/=kmrot) .or. (insl/=ins) .or. (nspinl/=nspin) .or. (lmmaxl/=lmmax) .or. (naezl/=nathost)) then
            write (6, '(/,5X,2A)') 'ERROR: ', 'host not compatible with your input/sytem'
            write (6, '(14X,6(A6),/,8X,42("-"))') '  KREL', ' KMROT', '   INS', ' NSPIN', ' LMMAX', ' BASIS'
            write (6, '(8X,A6,6I6)') 'syst: ', krel, kmrot, ins, nspin, lmmax, nathost
            write (6, '(8X,A6,6I6,/)') 'host: ', krell, kmrotl, insl, nspinl, lmmaxl, naezl
            stop '       < DECIMAREAD > '
          end if
          ! ......................................................................
          write (1337, 130) alatl, nspinl, naezl, lmmaxl, insl, krell, kmrotl

          read (36+ihost, 190) baner1
          read (36+ihost, 150) bravaisl
          write (1337, 140) bravaisl
          read (36+ihost, 190) baner1
          ih = lngstring(baner1, 80)
          write (1337, 200) baner1(1:ih)
          ! ......................................................................
          if (krel==0) then
            do ih = 1, naezl
              read (36+ihost, 150)(rbasisl(i,ih), i=1, 3)
              write (1337, 150)(rbasisl(i,ih), i=1, 3)
            end do
          else
            do ih = 1, naezl
              read (36+ihost, 260)(rbasisl(i,ih), i=1, 3), qmt, qmp
              write (1337, 160)(rbasisl(i,ih), i=1, 3), qmt, qmp
            end do
          end if
          ! ......................................................................
          read (36+ihost, 220) e2l, templ
          read (36+ihost, 230) npt1l, npt2l, npt3l, npoll
          ! ......................................................................
          if ((npt1l/=nptp1) .or. (npt2l/=nptp2) .or. (npt3l/=nptp3) .or. (npoll/=npol) .or. (abs(templ-tk)>1.e-4_dp)) then
            write (6, '(/,5X,2A)') 'ERROR: Host energy mesh', ' not compatible with your input'
            write (6, '(14X,5(A6),/,8X,40("-"))') '   NE1', '   NE2', '   NE3', '  NPOL', '  TEMP'
            write (6, '(8X,A6,4I6,2X,F10.6)') 'input:', nptp1, nptp2, nptp3, npol, tk
            write (6, '(8X,A6,4I6,2X,F10.6,/)') ' host:', npt1l, npt2l, npt3l, npoll, templ
            stop '       < DECIMAREAD > '
          end if
          ! ......................................................................
          write (1337, 170) e2l, templ
          write (1337, 180) npt1l, npt2l, npt3l, npoll
        end if
        ! ----------------------------------------------------------------------
        write (1337, '(8X,65("-"))')
      end do
      ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      ! Headers read in
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      return
    end if
    ! ========================================================== IENERGY = 0


    ! ======================================================================
    ! READ IN Delta_t matrices for IENERGY > 0
    ! ======================================================================
    if ((ispin/=2) .or. (nspin/=1)) then
      ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
      do ihost = 1, 2
        nathost = nlbasis
        if (ihost==2) nathost = nrbasis
        ! ----------------------------------------------------------------------
        ! --> read in Delta_t if it is not a vacuum host
        ! attention: new format Dec. 2004

        if (.not. vacflag(ihost)) then
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ HOST ATOM
          ! LOOP
          do ih = 1, nathost
            call cinit(lmmaxd*lmmaxd, w1)

100         continue
            read (36+ihost, '(A5)', end=120) str5
            if (str5/='*****') go to 100

            read (36+ihost, 240) iel, el, idum
            ! ......................................................................
            if ((abs(el-ez(ienergy))>1.e-4_dp) .or. (iel/=ienergy)) then
              write (6, '(/,5X,2A)') 'ERROR: Host energy mesh', ' not compatible with your input'
              write (6, '(14X,2(A6),/,8X,43("-"))') '    IE', ' E(IE)'
              write (6, '(8X,A6,I6,1X,2(1X,F10.6))') 'input:', ienergy, ez(ienergy)
              write (6, '(8X,A6,I6,1X,2(1X,F10.6),/)') ' host:', iel, el
              stop '       < DECIMAREAD > '
            end if
            ! ......................................................................
            if (idum/=ih) then
              write (6, '(/,5X,2A,/)') 'ERROR reading host file', ' basis indexing wrong'
              stop '       < DECIMAREAD > '
            end if
            ! ......................................................................

            ! --> map the t-matrix corectly

            ih1 = kaoez(1, naez+(ihost-1)*nlbasis+ih)

            ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
110         continue
            read (36+ihost, 250) lm1, lm2, el
            if ((lm1+lm2)/=0) then
              w1(lm1, lm2) = el
              if ((lm1+lm2)<2*lmmaxd) go to 110
            end if

            ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            if (ihost==1) then
              do lm1 = 1, lmmaxd   ! (KREL+KORBIT+1)*LMMAX
                ! CALL ZCOPY((KREL+KORBIT+1)*LMMAX,W1(1,LM1),1,
                call zcopy(lmmaxd, w1(1,lm1), 1, lefttinvll(1,lm1,ih1), 1)
              end do
            else
              do lm1 = 1, lmmaxd   ! (KREL+KORBIT+1)*LMMAX
                ! CALL ZCOPY((KREL+KORBIT+1)*LMMAX,W1(1,LM1),1,
                call zcopy(lmmaxd, w1(1,lm1), 1, righttinvll(1,lm1,ih1), 1)
              end do
            end if
          end do
          ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        end if                     ! not vacuum
        ! ----------------------------------------------------------------------
      end do
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    end if                         ! ispin.ne.2 .or. .nspin.ne.1
    ! ======================================================================

    if (ienergy==ielast) write (1337, *)
    return
120 continue
    stop '        Error reading hostfile'

130 format (10x, 'ALAT=', f9.6, ' NSPIN=', i2, '  NAEZ=', i3, ' LMMAX=', i3, ' INS=', i1, ' KREL=', i1, ' KMROT=', i1)
140 format (10x, 'BRAVAIS ', /, 10x, 3f8.4, /, 10x, 3f8.4, /, 10x, 3f8.4)
150 format (10x, 3f8.4)
160 format (10x, 3f8.4, 2f9.4)
170 format (10x, 'EF=', f10.6, ' TEMP=', f10.4, ' Kelvin')
180 format (10x, 'N1=', i3, ' N2=', i3, ' N3=', i3, ' NPOL=', i3)

190 format (a80)
200 format (10x, a)
210 format (5x, f9.6, 7x, i2, 7x, i3, 7x, i3, 5x, i1, 6x, i1, 7x, i1)
220 format (3x, f10.6, 6x, f10.4)
230 format (3x, i3, 4x, i3, 4x, i3, 6x, i3)
240 format (7x, i5, 2d16.8, 6x, i3)
250 format (2i5, 1p, 2d22.14)
260 format (3f8.4, 2f9.4)

  end subroutine decimaread

end module mod_decimaread
