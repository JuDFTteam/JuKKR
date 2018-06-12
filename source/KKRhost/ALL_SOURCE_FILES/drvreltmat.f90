subroutine drvreltmat(eryd, tmatll, vt, bt, r, drdi, r2drdi, zat, jws, solver, &
  soctl, ctl, lmmaxd, lmaxd, irmd)
  use :: mod_datatypes, only: dp
  ! ********************************************************************
  ! *                                                                  *
  ! * driving routine to call relativistic < SSITE > routine           *
  ! * only to calculate the single-site t matrix                       *
  ! * v.popescu, munich, may 2004                                      *
  ! *                                                                  *
  ! ********************************************************************

  implicit none

  ! PARAMETER definitions
  integer :: nrmax
  parameter (nrmax=900)
  integer :: nlamax, nqmax, ntmax, nmmax
  parameter (nlamax=1, nqmax=1, ntmax=1, nmmax=1)
  integer :: nlmax, nkmmax, nmuemax, nkmpmax, nkmax, linmax
  parameter (nlmax=5)              ! this should be >= LMAXD + 1
  parameter (nkmmax=2*nlmax**2, nkmax=2*nlmax-1)
  parameter (nkmpmax=nkmmax+2*nlmax, nmuemax=2*nlmax)
  parameter (linmax=2*nlmax*(2*nlmax-1))

  ! Dummy arguments
  integer :: lmaxd, lmmaxd, irmd
  integer :: zat(ntmax), jws(nmmax)
  complex (kind=dp) :: tmatll(lmmaxd, lmmaxd)
  real (kind=dp) :: soctl(nlmax)
  real (kind=dp) :: ctl(nlmax)
  real (kind=dp) :: vt(nrmax), bt(nrmax)
  real (kind=dp) :: r(nrmax, nmmax), r2drdi(nrmax, nmmax)
  real (kind=dp) :: drdi(nrmax, nmmax)

  ! Local variables
  real (kind=dp) :: ameopo(nkmmax, nkmmax, nlamax, 3), &
    at(nrmax, nlamax, 3, ntmax)
  complex (kind=dp) :: bzj(linmax, ntmax), bzz(linmax, ntmax), &
    dzj(linmax, ntmax), dzz(linmax, ntmax), eryd, msst(nkmmax, nkmmax, ntmax), &
    ozj(linmax, ntmax), ozz(linmax, ntmax), p, qzj(linmax, ntmax), &
    qzz(linmax, ntmax), szj(linmax, ntmax), szz(linmax, ntmax)
  complex (kind=dp) :: tsst(nkmmax, nkmmax, ntmax), tsstlin(linmax, ntmax), &
    tzj(linmax, ntmax), tzz(linmax, ntmax)
  complex (kind=dp) :: ozzs(linmax, ntmax, 2), ozjs(linmax, ntmax, 2)
  logical :: calcint, getirrsol
  real (kind=dp) :: cgc(nkmpmax, 2)
  integer :: i, ihyper, ikm1lin(linmax), ikm2lin(linmax), il, imt(ntmax), &
    imue, iprint, iq, iqat(nqmax, ntmax), it, iwrirrwf, iwrregwf, j, &
    lopt(ntmax), mmax, nkm, nkmq(nqmax), nl, nlinq(nqmax), nlq(nqmax), nt, &
    nucleus
  integer :: nfilcbwf
  integer :: nsollm(nlmax, nmuemax), ltab(nmuemax), kaptab(nmuemax), &
    nmuetab(nmuemax)
  character (len=10) :: solver
  integer :: icall

  data icall/0/

  save :: icall, ikm1lin, ikm2lin, lopt, nlq, nkmq, iqat, imt, nkm, ihyper, &
    iprint, it, nt, nucleus, iwrregwf, iwrirrwf, calcint, getirrsol, nfilcbwf

  icall = icall + 1

  ! =======================================================================
  ! initialise relativistic and dummy variables and SAVE them
  ! =======================================================================
  if (icall==1) then

    if (lmaxd>nlmax-1) then
      write (6, *) ' LMAXD = ', lmaxd, ' > NLMAX-1 = ', nlmax - 1
      stop ' Increase NLMAX in < DRVRELTMAT > '
    end if

    if (irmd>nrmax) then
      write (6, *) ' IRMD = ', irmd, ' > NRMAX = ', nrmax
      write (6, *) ' Increase NRMAX in < sprkkr_rmesh.dim > '
      stop ' and in < DRVRELTMAT > '
    end if

    nl = lmaxd + 1                 ! no need to save, used only here once
    iprint = 0

    do i = 1, nmuemax
      ltab(i) = i/2
      if (2*ltab(i)==i) then
        kaptab(i) = ltab(i)
      else
        kaptab(i) = -ltab(i) - 1
      end if
      nmuetab(i) = 2*abs(kaptab(i))
    end do

    do il = 1, nlmax
      mmax = 2*il
      do imue = 1, mmax
        if ((imue==1) .or. (imue==mmax)) then
          nsollm(il, imue) = 1
        else
          nsollm(il, imue) = 2
        end if
      end do
    end do

    call ikmlin(iprint, nsollm, ikm1lin, ikm2lin, nlmax, nmuemax, linmax, &
      nlmax)

    call calccgc(ltab, kaptab, nmuetab, cgc, nkmax, nmuemax, nkmpmax)

    do it = 1, ntmax
      imt(it) = 1
      lopt(it) = -1                ! this should change for Brooks' OP
    end do

    do iq = 1, nqmax
      nlq(iq) = nl
      nkmq(iq) = lmmaxd
      nlinq(iq) = 2*nlq(iq)*(2*nlq(iq)-1)
      iqat(iq, 1) = 1
    end do

    nkm = lmmaxd
    ihyper = 0
    nt = 1
    it = 1
    nucleus = 0

    iwrregwf = 0
    iwrirrwf = 0
    calcint = .false.
    getirrsol = .false.
    nfilcbwf = 0
  end if                           ! ICALL.EQ.1
  ! =======================================================================

  call ssite(iwrregwf, iwrirrwf, nfilcbwf, calcint, getirrsol, soctl, ctl, &
    eryd, p, ihyper, iprint, ikm1lin, ikm2lin, nlq, nkmq, nlinq, nt, nkm, &
    iqat, tsst, msst, tsstlin, dzz, dzj, szz, szj, ozz, ozj, bzz, bzj, qzz, &
    qzj, tzz, tzj, vt, bt, at, zat, nucleus, r, drdi, r2drdi, jws, imt, &
    ameopo, lopt, solver, cgc, ozzs, ozjs, nlmax, nqmax, linmax, nrmax, nmmax, &
    ntmax, nkmmax, nkmpmax, nlamax)

  do j = 1, nkm
    call zcopy(nkm, tsst(1,j,it), 1, tmatll(1,j), 1)
  end do

  return
end subroutine drvreltmat
