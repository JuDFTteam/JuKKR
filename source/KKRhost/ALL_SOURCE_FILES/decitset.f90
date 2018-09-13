module mod_decitset

contains

subroutine decitset(alat, bravsys, ez, ielast, nlbasis, nrbasis, fileleft, &
  fileright, ins, kvrel, krel, nspin, kmrot, vref, rmtref, nref, refpot, &
  lefttinv, righttinv, vacflag, nembd1, iemxd, irmd, ipand, lmaxd, lmgf0d, &
  lmmaxd, lm2d, nspind)
  use :: mod_datatypes, only: dp
  ! **********************************************************************
  ! *                                                                    *
  ! * This subroutine is thought as an alternative to the < decimaread > *
  ! * which requires an a priori calculated set of single-site matrices  *
  ! * over a fixed energy mesh. It is using the potential as written out *
  ! * in < outpothost > routine and determines the matrix                *
  ! *                                                                    *
  ! *            /         \-1   /             \-1                       *
  ! *            | Delta t |   = | t    - t    |                         *
  ! *            \         /     \  sys    ref /                         *
  ! *                                                                    *
  ! * for the left and the right host, using the energy mesh as read in  *
  ! * from the input-file                                                *
  ! *                                                                    *
  ! *                                        v.popescu - munich, Dec 04  *
  ! *                                                                    *
  ! * Notes: - no charge moments are calculated -- thus this option CAN  *
  ! *          NOT be used in SCF calculations                           *
  ! *        - non-spherical case not implemented, neither LDA+U (al-    *
  ! *          though the interface to regsol in decitmat is supplied)   *
  ! *        - CPA case not implemented - requires BZ integration        *
  ! *                                                                    *
  ! **********************************************************************
   use mod_calctref13
   use mod_decipotbas
   use mod_cmatstr
   use mod_decipothead
   use mod_decitmat
   use mod_lngstring
  implicit none
  ! ..
  ! .. Scalars arguments ..
  integer :: iemxd, nembd1, lmmaxd, ipand, nspind, irmd, lmaxd, lm2d, lmgf0d
  integer :: ielast, kmrot, nlbasis, nrbasis, nref, ins, kvrel, krel, nspin
  real (kind=dp) :: alat
  character (len=40) :: fileleft, fileright
  ! ..
  ! .. Array arguments ..
  integer :: refpot(nembd1)
  real (kind=dp) :: vref(*), rmtref(*), bravsys(3, 3)
  complex (kind=dp) :: ez(iemxd)
  complex (kind=dp) :: lefttinv(lmmaxd, lmmaxd, nembd1, nspind, iemxd), &
    righttinv(lmmaxd, lmmaxd, nembd1, nspind, iemxd)
  logical :: vacflag(2)
  ! ..
  ! .. Local scalars ..
  integer :: ihost, i, ll, mm, nqhost, ilhost
  integer :: nhost
  integer :: nq, nt, iqoff, itoff, ie, ih, iqh, ioq, info
  integer :: ipot, i1, ispin, nsra, lm1, lm2, irc1, iref
  integer :: ntleft, ntright, nthost
  ! .. LDA+U
  integer :: idoldau, lopt
  real (kind=dp) :: wldauav
  ! ..
  real (kind=dp) :: efermi, rirc
  complex (kind=dp) :: eryd, carg, cfctor
  character (len=40) :: filehost
  character (len=10) :: solver
  ! ..
  ! .. Local arrays
  integer :: krelh(2), nspinh(2), insh(2), ipvt(lmmaxd)
  integer :: noq(nembd1), kaoez(nembd1, nembd1), inhost(2)
  real (kind=dp) :: bravais(3, 3, 2), rbasis(3, nembd1), qmtet(nembd1), &
    qmphi(nembd1)
  character (len=5) :: chhost(2)
  character (len=9) :: txts(2)
  ! ..
  ! .. Allocatable local arrays
  integer :: ntmax
  real (kind=dp), allocatable :: zat(:), rws(:), rmt(:), conc(:)
  real (kind=dp), allocatable :: rr(:, :), drdi(:, :), visp(:, :), dror(:, :)
  real (kind=dp), allocatable :: socscl(:, :), cscl(:, :)
  integer, allocatable :: irws(:), ipan(:), iqat(:, :), ircut(:, :), loflm(:)
  complex (kind=dp), allocatable :: trefll(:, :, :), tmatll(:, :), dhmat(:, :, :)
  complex (kind=dp), allocatable :: dtrefll(:, :, :) ! LLY Lloyd
  complex (kind=dp), allocatable :: alpharef(:, :), dalpharef(:, :) ! LLY Lloyd Alpha
                                                       ! matrix and deriv.
  real (kind=dp), allocatable :: vtrel(:, :), btrel(:, :), r2drdirel(:, :)
  integer, allocatable :: zrel(:)
  ! ..
  ! .. External Functions
  logical, external :: test
  ! .. Data statements
  data chhost/'LEFT ', 'RIGHT'/
  data txts/'spin   UP', 'spin DOWN'/
  ! ......................................................................

  cfctor = alat/(8.e0_dp*atan(1.0e0_dp)) ! = ALAT/(2*PI)

  idoldau = 0
  lopt = -1
  wldauav = 0e0_dp
  allocate (loflm(lm2d), stat=i1)
  if (i1/=0) stop '    Allocate LOFLM'
  write (6, '(5X,A,/,8X,65("-"))') 'Reading in host potentials'
  vacflag(1) = .false.
  vacflag(2) = .false.
  nsra = 1
  if (kvrel>=1) nsra = 2
  i = 1
  do ll = 0, 2*lmaxd
    do mm = -ll, ll
      loflm(i) = ll
      i = i + 1
    end do
  end do
  ntleft = 0
  ntright = 0
  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
  nhost = 0
  do ihost = 1, 2
    filehost = fileleft
    nqhost = nlbasis
    if (ihost==2) then
      filehost = fileright
      nqhost = nrbasis
    end if
    ilhost = lngstring(filehost, 40)
    call decipothead(ihost, filehost, ilhost, nqhost, vacflag, alat, bravsys, &
      nq, nt, bravais(1,1,ihost), efermi, insh(ihost), krelh(ihost), &
      nspinh(ihost), ins, krel, nspin, kmrot)

    if (.not. vacflag(ihost)) then
      nhost = nhost + 1
      inhost(nhost) = ihost
      if (ihost==1) then
        ntleft = nt
      else
        ntright = nt
      end if
    end if


  end do

  if (ntleft+ntright<=0) then
    write (6, '(8X,"Vacuum will be considered on both sides",/, 8X,65("-"))')
    return
  end if

  ntmax = ntleft + ntright
  allocate (zat(ntmax), rws(ntmax), rmt(ntmax), conc(ntmax), stat=i1)
  if (i1/=0) stop '    Allocate ZAT/RWS/RMT/CONC'
  allocate (rr(irmd,ntmax), drdi(irmd,ntmax), stat=i1)
  if (i1/=0) stop '    Allocate RR/DRDI'
  allocate (visp(irmd,ntmax*nspind), stat=i1)
  if (i1/=0) stop '    Allocate VISP'
  allocate (irws(ntmax), ipan(ntmax), iqat(nembd1,ntmax), stat=i1)
  if (i1/=0) stop '    Allocate IRWS/IPAN/IQAT'
  allocate (ircut(0:ipand,ntmax), stat=i1)
  if (i1/=0) stop '    Allocate IRCUT'
  allocate (socscl(krel*lmaxd+1,krel*ntmax+(1-krel)), stat=i1)
  if (i1/=0) stop '    Allocate SOCSCL'
  allocate (cscl(krel*lmaxd+1,krel*ntmax+(1-krel)), stat=i1)
  if (i1/=0) stop '    Allocate CSCL'
  allocate (vtrel(irmd*krel+(1-krel),ntmax), stat=i1)
  if (i1/=0) stop '    Allocate VTREL'
  allocate (btrel(irmd*krel+(1-krel),ntmax), stat=i1)
  if (i1/=0) stop '    Allocate BTREL'

  ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
  do ihost = 1, 2
    write (6, '(8X,A5," side host: ")', advance='no') chhost(ihost)
    iqoff = 0
    itoff = 0
    nqhost = nlbasis
    nthost = ntleft
    filehost = fileleft
    if (ihost==2) then
      nqhost = nrbasis
      nthost = ntright
      iqoff = nlbasis
      itoff = ntleft
      filehost = fileright
    end if
    ilhost = lngstring(filehost, 40)

    if (filehost(1:7)=='vacuum') then
      write (6, '(A,/,8X,65("-"))') 'VACUUM will be used'
    else
      write (6, '(A,/)') filehost(1:ilhost)
      write (6, 130) krelh(ihost), nspinh(ihost), insh(ihost), kmrot, nqhost, &
        alat, efermi
      write (6, 140)((bravais(ll,mm,ihost),mm=1,3), ll=1, 3)
      call decipotbas(ihost, iqoff, itoff, nqhost, nthost, rbasis, qmtet, &
        qmphi, noq, kaoez, zat, iqat, conc, irws, ipan, ircut, rr, drdi, visp, &
        nspinh(ihost), krelh(ihost), solver, socscl, cscl, vtrel, btrel, irmd, &
        ipand, nembd1, ntmax, nspind, lmaxd)
      write (6, '(8X,65("-"))')
    end if
  end do
  ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  allocate (dror(irmd,ntmax), stat=i1)
  if (i1/=0) stop '    Allocate DROR'
  allocate (r2drdirel(irmd*krel+(1-krel),ntmax), stat=i1)
  if (i1/=0) stop '    Allocate R2DRDIREL'
  allocate (zrel(ntmax), stat=i1)
  if (i1/=0) stop '    Allocate ZREL'

  write (6, '(/,5X,A,/)') 'Calculating host (Delta_t)^(-1) matrices'
  if (krel==0) then
    do i = 1, ntleft + ntright
      irc1 = ircut(ipan(i), i)
      do i1 = 2, irc1
        dror(i1, i) = drdi(i1, i)/rr(i1, i)
      end do
    end do
  else
    do i = 1, ntleft + ntright
      irc1 = ircut(ipan(i), i)
      do i1 = 1, irc1
        r2drdirel(i1, i) = rr(i1, i)*rr(i1, i)*drdi(i1, i)
      end do
      zrel(i) = nint(zat(i))
    end do
  end if

  ! ******************************************************* energy loop IE
  allocate (trefll(lmmaxd,lmmaxd,nref), stat=i1)
  if (i1/=0) stop '    Allocate TREFLL'
  allocate (dtrefll(lmmaxd,lmmaxd,nref), stat=i1) ! LLY
  if (i1/=0) stop '    Allocate DTREFLL' ! LLY
  allocate (tmatll(lmmaxd,lmmaxd), dhmat(lmmaxd,lmmaxd,2), stat=i1)
  if (i1/=0) stop '    Allocate TMATLL/DHMAT'
  allocate (alpharef(0:lmaxd,nref), dalpharef(0:lmaxd,nref), stat=i1) ! LLY
                                                                      ! Lloyd
                                                                      ! Alpha
                                                                      ! matrix
                                                                      ! AND
                                                                      ! deriv.
  if (i1/=0) stop '    Allocate ALPHAREF/DALPHAREF'

  do ie = 1, ielast
    eryd = ez(ie)

    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    ! -> set up t matrices for the reference system

    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if (krel==0) then
      do i1 = 1, nref
        call calctref13(eryd, vref(i1), rmtref(i1), lmaxd, ih, trefll(1,1,i1), &
          dtrefll(1,1,i1), alpharef(0,i1), dalpharef(0,i1), lmaxd+1, lmgf0d)
      end do
    end if
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: HOST-LOOP
    do ilhost = 1, nhost
      ihost = inhost(ilhost)
      iqoff = 0
      itoff = 0
      nqhost = nlbasis
      if (ihost==2) then
        nqhost = nrbasis
        iqoff = nlbasis
        itoff = ntleft
      end if
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ sites in host
      do ih = 1, nqhost

        ! -> assign Delta_t = -t_ref

        ! Note: REFPOT(1) = REFPOT(NAEZ) (i.e., of the 2D system)

        iqh = iqoff + ih
        iref = refpot(iqh+1)
        do lm2 = 1, lmmaxd
          do lm1 = 1, lmmaxd
            dhmat(lm1, lm2, 1) = -trefll(lm1, lm2, iref)
          end do
        end do

        if (nspinh(ihost)>1) then
          do lm2 = 1, lmmaxd
            call zcopy(lmmaxd, dhmat(1,lm2,1), 1, dhmat(1,lm2,2), 1)
          end do
        end if
        ! ====================================================== spins and
        ! atoms
        do ispin = 1, nspinh(ihost)
          ! ----------------------------------------------------------------------
          do ioq = 1, noq(iqh)
            i1 = kaoez(ioq, iqh) + itoff
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! -> calculate t_sys for the atom I1 located on site IH

            ipot = (i1-1)*nspinh(ihost) + ispin
            irc1 = ircut(ipan(i1), i1)
            rirc = rr(irc1, i1)

            call decitmat(eryd, zat(i1), ipan(i1), rr(1,i1), dror(1,i1), &
              visp(1,ipot), ircut(0,i1), rirc, krel, nsra, ins, tmatll, loflm, &
              idoldau, lopt, wldauav, solver, socscl(1,krel*i1+(1- &
              krel)), cscl(1,krel*i1+(1-krel)), zrel(i1), vtrel(1,i1), btrel(1 &
              ,i1), drdi(1,i1), r2drdirel(1,i1), ipand, irmd, lmaxd, lmaxd+1, &
              lm2d, lmmaxd)

            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% tmat
            ! calculated

            ! -> Delta_t = Delta_t + CONC(I1)*t_mat(I1)

            carg = conc(i1)
            do lm2 = 1, lmmaxd
              call zaxpy(lmmaxd, carg, tmatll(1,lm2), 1, dhmat(1,lm2,ispin), &
                1)
            end do
          end do
          ! ----------------------------------------------------------------------
          ! tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt
          if (test('tmat    ')) then
            write (1337, *)
            write (1337, 100, advance='no') &
              '      ---> Delta_t  matrix for site: ', iqh
            if (krel==0) write (1337, 110, advance='no') txts(ispin)
            write (1337, 120) ', energy: ', eryd
            call cmatstr(' ', 1, dhmat(1,1,ispin), lmmaxd, lmmaxd, 2*krel+1, &
              2*krel+1, 0, 1e-8_dp, 6)
            write (1337, *)
          end if
          ! tttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttttt

          ! --> inversion

          call zgetrf(lmmaxd, lmmaxd, dhmat(1,1,ispin), lmmaxd, ipvt, info)
          call zgetri(lmmaxd, dhmat(1,1,ispin), lmmaxd, ipvt, tmatll, &
            lmmaxd*lmmaxd, info)
        end do
        ! ======================================================================

        ! --> scaling the host t-matrices to p.u.

        if (ihost==1) then
          do ispin = 1, nspinh(ihost)
            do lm2 = 1, lmmaxd
              do lm1 = 1, lmmaxd
                lefttinv(lm1, lm2, ih, ispin, ie) = cfctor* &
                  dhmat(lm1, lm2, ispin)
              end do
            end do
          end do
        else
          do ispin = 1, nspinh(ihost)
            do lm2 = 1, lmmaxd
              do lm1 = 1, lmmaxd
                righttinv(lm1, lm2, ih, ispin, ie) = cfctor* &
                  dhmat(lm1, lm2, ispin)
              end do
            end do
          end do
        end if
      end do
      ! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    end do
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do
  ! **********************************************************************
  deallocate (zat, rws, rmt, conc, rr, drdi, visp, stat=i1)
  if (i1/=0) stop '   Deallocate ZAT/RWS/RMT/.../VISP'
  deallocate (irws, ipan, iqat, ircut, loflm, stat=i1)
  if (i1/=0) stop '   Deallocate IRWS/IPAN/IQAT/IRCUT/LOFLM'
  deallocate (trefll, tmatll, dhmat, stat=i1)
  if (i1/=0) stop '   Deallocate TREFLL/TMATLL/DHMAT'
  deallocate (socscl, cscl, vtrel, btrel, stat=i1)
  if (i1/=0) stop '   Deallocate SOCSCL/CSCL/VTREL/BTREL'
  deallocate (alpharef, dalpharef, stat=i1)
  if (i1/=0) stop '   Deallocate ALPHAREF/DALPHAREF'
  if (krel==0) then
    deallocate (dror, stat=i1)
    if (i1/=0) stop '   Deallocate DROR'
  else
    deallocate (r2drdirel, zrel, stat=i1)
    if (i1/=0) stop '   Deallocate R2DRDIREL/ZREL'
  end if
100 format (a, i3)
110 format (', ', a)
120 format (a, 2f10.6)
130 format (10x, 'KREL= ', i1, ' NSPIN= ', i1, ' INS= ', i1, ' KMROT= ', i1, &
    /, 10x, 'NAEZ=', i3, ' ALAT= ', f9.6, ' EFERMI= ', f9.6)
140 format (10x, 'BRAVAIS ', /, 10x, 3f8.4, /, 10x, 3f8.4, /, 10x, 3f8.4, /, &
    10x, 'RBASIS')
end subroutine decitset

end module mod_decitset
