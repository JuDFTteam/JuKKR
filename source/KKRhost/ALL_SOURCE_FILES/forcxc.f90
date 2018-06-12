subroutine forcxc(flm, flmc, lmax, nspin, nstart, nend, rhoc, v, r, alat, &
  drdi, irws, natref)
  ! >>>>>BEWARE!!! RM commented away!!! -->Dipole Tensor is useless
  ! SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,ALAT,
  ! +                  RM,NSHELL,DRDI,IRWS,NATREF)
  ! -----------------------------------------------------------------------
  ! calculates the force on nucleus m
  ! from a given non spherical charge density at the nucleus site r
  ! with core correction(exchange contribution)

  ! -----------------------------------------------------------------------
  use :: mod_types, only: t_inc
  use :: mod_datatypes, only: dp
  implicit none
  ! .. Parameters ..
  include 'inc.p'
  integer :: lmpotd
  parameter (lmpotd=(lpotd+1)**2)
  ! $     RM(3,*),
  ! INTEGER IRWS(*),NSHELL(*)
  real (kind=dp) :: alat
  integer :: lmax, natref, nend, nspin, nstart
  ! ..
  ! .. Local Scalars ..
  real (kind=dp) :: drdi(irmd, *), flm(-1:1, *), flmc(-1:1, *), r(irmd, *), &
    rhoc(irmd, *), &               ! ,DVOL
    v(irmd, lmpotd, *)
  ! ,J
  integer :: irws(*)
  ! ..
  ! .. Local Arrays ..
  real (kind=dp) :: dv, fac, pi, rws, trp, vint1, vol ! ..
  integer :: i, iatyp, iper, ipot, irep, irws1, ispin, lm, m ! .. External
                                                             ! Subroutines ..
  ! ..
  ! .. Save statement ..
  real (kind=dp) :: f(3, natypd), flmh(-1:1, natypd), flmxc(-1:1, natypd), &
    p(natypd), v1(irmd)
  ! ..
  ! .. Intrinsic Functions ..
  external :: simp3
  ! ..

  save :: pi


  intrinsic :: atan, dsqrt

  pi = 4.d0*atan(1.d0)
  fac = dsqrt((4.0d0*pi)/3.0d0)
  trp = 0.0d0
  if (lmax<1) then
    write (6, fmt=100)
    stop

  end if

  if (t_inc%i_write>0) write (1337, fmt=130)
  if (t_inc%i_write>0) write (1337, fmt=120)
  if (t_inc%i_write>0) write (1337, fmt=130)
  ! ---> determine the right potential numbers
  irep = 1
  do iatyp = nstart, nend

    iper = iatyp - natref
    p(iper) = 0.0d0
    if (t_inc%i_write>0) write (1337, fmt=140) iper

    irws1 = irws(iatyp)
    rws = r(irws1, iatyp)
    vol = 0.25*alat**3




    do m = -1, 1
      lm = 2 + m + 1

      do i = 1, irws1
        v1(i) = 0.0d0
      end do

      do ispin = 1, nspin

        ipot = nspin*(iatyp-1) + ispin


        dv = (-3.0d0*v(1,lm,ipot)-10.0d0*v(2,lm,ipot)+18.0d0*v(3,lm,ipot)- &
          6.0d0*v(4,lm,ipot)+v(5,lm,ipot))/(12.0d0*drdi(2,iatyp))

        v1(2) = rhoc(2, ipot)*(2.0d0*v(2,lm,ipot)/r(2,iatyp)+dv)/(4.0d0*pi) + &
          v1(2)

        do i = 3, irws1 - 2

          dv = (v(i-2,lm,ipot)-v(i+2,lm,ipot)+8.0d0*(v(i+1,lm,ipot)-v(i-1,lm, &
            ipot)))/(12.0d0*drdi(i,iatyp))

          v1(i) = rhoc(i, ipot)*(2.0d0*v(i,lm,ipot)/r(i,iatyp)+dv)/ &
            (4.0d0*pi) + v1(i)
        end do
        ! ---> integrate with simpson subroutine
        dv = (-v(irws1-4,lm,ipot)+6.0d0*v(irws1-3,lm,ipot)- &
          18.0d0*v(irws1-2,lm,ipot)+10.0d0*v(irws1-1,lm,ipot)+ &
          3.0d0*v(irws1,lm,ipot))/(12.0d0*drdi(irws1-1,iatyp))
        v1(irws1-1) = rhoc(irws1-1, ipot)*(2.0d0*v(irws1-1,lm,ipot)/r(irws1-1, &
          iatyp)+dv)/(4.0d0*pi) + v1(irws1-1)

        dv = (3.0d0*v(irws1-4,lm,ipot)-16.0d0*v(irws1-3,lm,ipot)+ &
          36.0d0*v(irws1-2,lm,ipot)-48.0d0*v(irws1-1,lm,ipot)+ &
          25.0d0*v(irws1,lm,ipot))/(12.0d0*drdi(irws1,iatyp))

        v1(irws1) = rhoc(irws1, ipot)*(2.0d0*v(irws1,lm,ipot)/r(irws1,iatyp)+ &
          dv)/(4.0d0*pi) + v1(irws1)
      end do



      call simp3(v1, vint1, 1, irws1, drdi(1,iatyp))

      flmh(m, iatyp) = flm(m, iatyp) - flmc(m, iatyp)
      flmxc(m, iatyp) = -fac*vint1 - flmc(m, iatyp)
      flm(m, iatyp) = flm(m, iatyp) + flmxc(m, iatyp)


    end do
    ! DO 60 J = 1,3
    if (t_inc%i_write>0) then
      write (1337, fmt=150) flmh(1, iatyp), flmc(1, iatyp), flmxc(1, iatyp), &
        flm(1, iatyp)
      write (1337, fmt=160) flmh(-1, iatyp), flmc(-1, iatyp), &
        flmxc(-1, iatyp), flm(-1, iatyp)
      write (1337, fmt=170) flmh(0, iatyp), flmc(0, iatyp), flmxc(0, iatyp), &
        flm(0, iatyp)
    end if
    ! P(IPER) = P(IPER) + RM(J,IREP)*NSHELL(IPER)*F(J,IATYP)*ALAT
    f(1, iatyp) = flm(1, iatyp)
    f(2, iatyp) = flm(-1, iatyp)
    f(3, iatyp) = flm(0, iatyp)
    ! 60    CONTINUE
    ! TRP = TRP + P(IPER)

    ! IREP = IREP + NSHELL(IPER)

    ! write (6,*) '-->Tensor is useless'
    ! WRITE (6,FMT=9700) P(IPER)


    ! DVOL = TRP/ (3.0D0*VOL)

    ! WRITE (6,FMT=9101)
  end do
  ! WRITE (6,FMT=9200)
  ! WRITE (6,FMT=9800) DVOL
  ! WRITE (6,FMT=9200)
  if (t_inc%i_write>0) then
    write (1337, fmt=130)

    ! 9101 FORMAT (1x,33 ('-'),' volume change ',33 ('-'),/,34x,
    ! +       ' in units Ry/(a(Bohr)**3 ')
    ! 9700 FORMAT (10x,'contribution to the trace of the dipol force tensor:'
    write (1337, fmt=110)
    write (1337, fmt=130)
  end if
  ! +       ,3x,e12.6,' Ry')
100 format (13x, 'error stop in subroutine force :', &
    ' the charge density has to contain non spherical', &
    ' contributions up to l=1 at least ')
  ! 9800 FORMAT (7x,' volume change dvol/vol=',2x,e12.6,' Ry/(a(Bohr))**3',
  ! +       /,7x,'( notice: has to be divided',
110 format (1x, 81('-'))
120 format (1x, 33('-'), ' force on the nucleus ', 33('-'), /, 34x, &
    ' in units Ry/(a(Bohr) ')
130 format (1x, '>')
140 format (3x, i5, 'th shell')
150 format (7x, 'fhx=', e13.6, 2x, 'fcx=', e13.6, 2x, 'fxcx=', e13.6, 2x, &
    'fx=', e13.6, ' Ry/(a(Bohr))')
160 format (7x, 'fhy=', e13.6, 2x, 'fcy=', e13.6, 2x, 'fxcy=', e13.6, 2x, &
    'fy=', e13.6, ' Ry/(a(Bohr))')
170 format (7x, 'fhz=', e13.6, 2x, 'fcz=', e13.6, 2x, 'fxcz=', e13.6, 2x, &
    'fz=', e13.6, ' Ry/(a(Bohr))')
  ! +       ' by the bulk modulus of the host)')

  ! >>>>>BEWARE!!! RM commented away!!! -->Dipole Tensor is useless
  ! SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,ALAT,
  ! +                  RM,NSHELL,DRDI,IRWS,NATREF)
  ! -----------------------------------------------------------------------
end subroutine forcxc
