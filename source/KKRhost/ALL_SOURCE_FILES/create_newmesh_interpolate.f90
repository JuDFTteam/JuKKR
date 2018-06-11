subroutine create_newmesh(nspin, r, irmin, irws, ipan, ircut, vins, visp, &
  r_log, npan_log, npan_eq, ncheb, npan_tot, rnew, rpan_intervall, &
  ipan_intervall, vinsnew, ntcell, thetas, thetasnew)

      Use mod_datatypes, Only: dp
  implicit none
  include 'inc.p'
  integer :: nspin, irmin(natypd), ipan(natypd), irws(natypd)
  integer :: lmmaxd
  parameter (lmmaxd=(lmaxd+1)**2)
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  integer :: irmind
  parameter (irmind=irmd-irnsd)
  integer :: npan_log, npan_eq, ncheb, npan_inst, npan_tot(natypd)
  integer :: ircut(0:ipand, natypd)
  real (kind=dp) :: r(irmd, natypd)
  real (kind=dp) :: fac
  parameter (fac=2d0)
  real (kind=dp) :: vins(irmind:irmd, lmpotd, nspotd), visp(irmd, nspotd), &
    vinsin(irmd, lmpotd, nspin)
  real (kind=dp) :: thetas(irid, nfund, ncelld), thetasin(irid, nfund, &
    ncelld), thetasnew(ntotd*(nchebd+1), nfund, ncelld)
  integer :: ntcell(natypd)
  integer :: i1, ipot, ipotm, ir, ispin, ir2, ip, icell, ishift, ilogpanshift, &
    ilinpanshift, npan_logtemp, imin, imax, iminnew, imaxnew, lm1
  real (kind=dp) :: r_log, rmin, rmax, rval
  real (kind=dp) :: rnew(ntotd*(nchebd+1), natypd), &
    rpan_intervall(0:ntotd, natypd)
  integer :: ipan_intervall(0:ntotd, natypd)
  real (kind=dp) :: vinsnew(ntotd*(nchebd+1), lmpotd, nspotd)

  vinsnew = 0d0
  thetasnew = 0d0
  ipotm = 0
! log panel
  do i1 = 1, natypd

    ipot = nspin*(i1-1) + 1
    npan_inst = ipan(i1) - 1
    npan_tot(i1) = npan_log + npan_eq + npan_inst

! NPAN_LOG

    rmin = r(2, i1)
    rmax = r_log
    rval = 0d0
    ishift = 0
    if (r_log>r(irmin(i1),i1)) then
      ilogpanshift = 1
      ilinpanshift = 0
    else
      ilogpanshift = 0
      ilinpanshift = 1
    end if
! equivalent panel
    if (ilinpanshift==1) then
      write (*, *) 'ERORR: non-spherical part of the potential needs'
      write (*, *) 'to be inside the log panel'
      write (*, *) 'atom (I1):', i1
      write (*, *) 'R_LOG', r_log
      write (*, *) 'R(IRMIN(I1), I1)', r(irmin(i1), i1)
      write (*, *) 'IRMIN(I1)', irmin(i1)
      stop 'Error creating newmesh'
    end if
! NPAN_EQ
    do ip = 0, npan_log - ilogpanshift
      rval = (fac**ip-1d0)/(fac**(npan_log-ilogpanshift)-1d0)
      rpan_intervall(ip+ishift, i1) = rmin + rval*(rmax-rmin)
      ipan_intervall(ip+ishift, i1) = (ip+ishift)*(ncheb+1)
      if (ishift==0 .and. rpan_intervall(ip,i1)>r(irmin(i1),i1)) then
        ishift = 1
        npan_logtemp = ip
        rpan_intervall(ip+1, i1) = rpan_intervall(ip, i1)
        ipan_intervall(ip+1, i1) = (ip+ishift)*(ncheb+1)
        rpan_intervall(ip, i1) = r(irmin(i1), i1)
        ipan_intervall(ip, i1) = ip*(ncheb+1)
      end if
    end do 
! intersection zone
! NPAN_INST
    ishift = 0
    rmin = r_log
    rmax = r(ircut(1,i1), i1)
    do ip = 0, npan_eq - ilinpanshift
      rpan_intervall(ip+ishift+npan_log, i1) = rmin + &
        ip*(rmax-rmin)/(npan_eq-ilinpanshift)
      ipan_intervall(ip+ishift+npan_log, i1) = (npan_log+ip+ishift)*(ncheb+1)
    end do 


    do ip = 1, npan_inst
      rpan_intervall(npan_log+npan_eq+ip, i1) = r(ircut(ip+1,i1), i1)
      ipan_intervall(npan_log+npan_eq+ip, i1) = (npan_log+npan_eq+ip)* &
        (ncheb+1)
    end do ! interpolate potential to new mesh
! save input potential to VINSIN
    npan_eq = npan_eq + npan_log - npan_logtemp
    npan_log = npan_logtemp

    call chebmesh(npan_tot, ncheb, rpan_intervall(0:,i1), rnew(1,i1))




    vinsin = 0d0
    do ir = 1, irws(i1)
      if (ir<irmin(i1)) then
        vinsin(ir, 1, 1) = visp(ir, ipot)
        vinsin(ir, 1, nspin) = visp(ir, ipot+nspin-1)
      else
        do lm1 = 1, lmpotd
          if (lm1==1) then
            vinsin(ir, lm1, 1) = visp(ir, ipot)
            vinsin(ir, lm1, nspin) = visp(ir, ipot+nspin-1)
          else
            vinsin(ir, lm1, 1) = vins(ir, lm1, ipot)
            vinsin(ir, lm1, nspin) = vins(ir, lm1, ipot+nspin-1)
          end if
        end do
      end if
    end do
! LM1
    do ispin = 1, nspin

      ipotm = ipotm + 1
! ISPIN
      do lm1 = 1, lmpotd
        imin = 1
        imax = irmin(i1)
        do ip = 1, npan_log
          iminnew = ipan_intervall(ip-1, i1) + 1
          imaxnew = ipan_intervall(ip, i1)
          call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm), &
            imax-imin+1, imaxnew-iminnew+1)
        end do
! interpolate shape function THETAS to new shape function THETASNEW
        imin = irmin(i1)
        imax = ircut(1, i1)
        do ip = npan_log + 1, npan_log + npan_eq
          iminnew = ipan_intervall(ip-1, i1) + 1
          imaxnew = ipan_intervall(ip, i1)
          call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm), &
            imax-imin+1, imaxnew-iminnew+1)
        end do
        ir2 = 0
        do ip = npan_log + npan_eq + 1, npan_tot(i1)
          ir2 = ir2 + 1
          imin = ircut(ir2, i1) + 1
          imax = ircut(ir2+1, i1)
          iminnew = ipan_intervall(ip-1, i1) + 1
          imaxnew = ipan_intervall(ip, i1)
          call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
            vinsin(imin:imax,lm1,ispin), vinsnew(iminnew:imaxnew,lm1,ipotm), &
            imax-imin+1, imaxnew-iminnew+1)
        end do
      end do ! save THETAS to THETASIN
! I1
    end do 
! set to 1 if NEWSOSOL under RUNOPT, otherwise 0
! SET ACCORDING TO lmax VALUE OF INPUTCARD
    icell = ntcell(i1)
    do lm1 = 1, nfund
      thetasin(:, lm1, icell) = thetas(:, lm1, icell)
      ir2 = 0
      do ip = npan_log + npan_eq + 1, npan_tot(i1)
        ir2 = ir2 + 1
        imin = ircut(ir2, i1) + 1
        imax = ircut(ir2+1, i1)
        iminnew = ipan_intervall(ip-1, i1) + 1
        imaxnew = ipan_intervall(ip, i1)
        call interpolspline(r(imin:imax,i1), rnew(iminnew:imaxnew,i1), &
          thetasin(imin-ircut(1,i1):imax-ircut(1,i1),lm1,icell), thetasnew( &
          iminnew:imaxnew,lm1,icell), imax-imin+1, imaxnew-iminnew+1)
      end do
    end do
  end do !      PARAMETER ( NRD = 20000, KPOIBZ = 32000 )
end subroutine



subroutine chebmesh(npan, ncheb, ri, ro)

  integer, intent (in) :: npan
  integer, intent (in) :: ncheb
  real (kind=dp), intent (in) :: ri(0:npan)
  real (kind=dp), intent (out) :: ro(npan*(ncheb+1))
  implicit none

!-----------------------------------------------

  integer :: i, k, ik
  real (kind=dp) :: tau, pi

  pi = 4d0*datan(1d0)
  do i = 1, npan
    do k = 0, ncheb
      ik = i*ncheb + i - k
      tau = dcos(((2*k+1)*pi)/(2*(ncheb+1)))
      tau = 0.5d0*((ri(i)-ri(i-1))*tau+ri(i)+ri(i-1))
      ro(ik) = tau
    end do
  end do
end subroutine

