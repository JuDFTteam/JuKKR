subroutine gradr(nspin, ist1, mesh, dx, drdi, drdi2, ro, zta, drr, ddrr, drru, &
  ddrru, rou, irmd)
!-----------------------------------------------------------------
!evaluates d(ro)/dr,d{d(ro)/dr}/dr.
!drr=d(ro)/dr, ddrr=d(drr)/dr.
!coded by T.Asada. Feb.1994.
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!------------------------------------------------------------------
  use :: mod_types, only: t_inc
      Use mod_datatypes, Only: dp
  implicit none
!.. Scalar Arguments ..
  real (kind=dp) :: dx
  integer :: irmd, ist1, mesh, nspin
!..
!.. Array Arguments ..
  real (kind=dp) :: ddrr(irmd), ddrru(irmd), drdi(irmd), drdi2(irmd), &
    drr(irmd), drru(irmd), ro(irmd), rou(irmd), zta(irmd)
!..
!.. Local Scalars ..
  real (kind=dp) :: d, drx, drx0, drx1, drx2, drx3, drxu, drxu0, drxu1, &
    drxu2, drxu3, drxx, drxx0, drxx1, drxx2, drxx3, drxxu, drxxu0, drxxu1, &
    drxxu2, drxxu3, f0, f1, f2, f3, f4, f5, g1, g2, g3, g4, g5, xlf
  integer :: i, i1, i2, i3, i4, i5, i6, ica, icg, iex, igd, igh, igl, ihb, &
    imj, ip9, ipg, ipw, ist, ivg, ivn, iwr, ixlf, j, ndvpt, nred
!..
!.. Statement Functions ..
  real (kind=dp) :: f131, f132, f133, f141, f142, f143, f144, f151, f152, &
    f153, f154, f155, f161, f162, f163, f164, f165, f166, f231, f232, f233, &
    f241, f242, f243, f244, f251, f252, f253, f254, f255, f261, f262, f263, &
    f264, f265, f266
!..
!.. Intrinsic Functions ..
  intrinsic :: dble
!..
!.. Save statement ..
  save :: ndvpt, igl, igh, imj, ica, icg, ivn, ipw, ipg, ivg, ip9, igd, ixlf, &
    iex, xlf, iwr
!..
!.. Data statements ..
!.....-----------------------------------------------------------------
!..uble function
!..ouble precision f131,f132,f133,f141,f142,f143,f144
!..ouble precision fl61,fl62,fl63,fl64,fl65,fl66
!..ouble precision fl51,fl52,fl53,fl54,fl55
!..ouble precision f231,f232,f233,f241,f242,f243,f244
!..ouble precision f251,f252,f253,f254,f255
!..ouble precision f261,f262,f263,f264,f265,f266

  data ndvpt/5/
  data igl, igh, imj, ica, icg, ivn, ipw, ipg, ivg, ip9, igd, ixlf, iex, &
    xlf/0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0.00d0/
  data iwr/0/
!     ..
!     .. Statement Function definitions ..
!  statement functions:

!.....three point formula for the 1st deriv.

!.....four point formula for the 1st deriv.

!.....five point formula for the 1st deriv.

!.....six point formula for the 1st deriv.

!.....three point formula for the 2nd deriv.

!.....four point formula for the 2nd deriv.

!.....five point formula for the 2nd deriv.

!.....six point formula for the 2nd deriv.
  f131(f0, f1, f2, d) = (-3*f0+4*f1-f2)/(2*d)
  f132(g1, f0, f1, d) = (-1*g1-0*f0+f1)/(2*d)
  f133(g2, g1, f0, d) = (g2-4*g1+3*f0)/(2*d)
  f141(f0, f1, f2, f3, d) = (-11*f0+18*f1-9*f2+2*f3)/(6*d)
  f142(g1, f0, f1, f2, d) = (-2*g1-3*f0+6*f1-f2)/(6*d)
  f143(g2, g1, f0, f1, d) = (g2-6*g1+3*f0+2*f1)/(6*d)
  f144(g3, g2, g1, f0, d) = (-2*g3+9*g2-18*g1+11*f0)/(6*d)
  f151(f0, f1, f2, f3, f4, d) = (-50*f0+96*f1-72*f2+32*f3-6*f4)/(24*d)
  f152(g1, f0, f1, f2, f3, d) = (-6*g1-20*f0+36*f1-12*f2+2*f3)/(24*d)
  f153(g2, g1, f0, f1, f2, d) = (2*g2-16*g1-0*f0+16*f1-2*f2)/(24*d)
  f154(g3, g2, g1, f0, f1, d) = (-2*g3+12*g2-36*g1+20*f0+6*f1)/(24*d)
  f155(g4, g3, g2, g1, f0, d) = (6*g4-32*g3+72*g2-96*g1+50*f0)/(24*d)
  f161(f0, f1, f2, f3, f4, f5, d) = (-274*f0+600*f1-600*f2+400*f3-150*f4+24*f5 &
    )/(120*d)
  f162(g1, f0, f1, f2, f3, f4, d) = (-24*g1-130*f0+240*f1-120*f2+40*f3-6*f4)/ &
    (120*d)
  f163(g2, g1, f0, f1, f2, f3, d) = (6*g2-60*g1-40*f0+120*f1-30*f2+4*f3)/ &
    (120*d)
  f164(g3, g2, g1, f0, f1, f2, d) = (-4*g3+30*g2-120*g1+40*f0+60*f1-6*f2)/ &
    (120*d)
  f165(g4, g3, g2, g1, f0, f1, d) = (6*g4-40*g3+120*g2-240*g1+130*f0+24*f1)/ &
    (120*d)
  f166(g5, g4, g3, g2, g1, f0, d) = (-24*g5+150*g4-400*g3+600*g2-600*g1+274*f0 &
    )/(120*d)
  f231(f0, f1, f2, d) = (f0-2*f1+f2)/(d*d)
  f232(g1, f0, f1, d) = (g1-2*f0+f1)/(d*d)
  f233(g2, g1, f0, d) = (g2-2*g1+f0)/(d*d)
  f241(f0, f1, f2, f3, d) = (6*f0-15*f1+12*f2-3*f3)/(3*d*d)
  f242(g1, f0, f1, f2, d) = (3*g1-6*f0+3*f1+0*f2)/(3*d*d)
  f243(g2, g1, f0, f1, d) = (0*g2+3*g1-6*f0+3*f1)/(3*d*d)
  f244(g3, g2, g1, f0, d) = (-3*g3+2*g2+15*g1+6*f0)/(3*d*d)
  f251(f0, f1, f2, f3, f4, d) = (35*f0-104*f1+114*f2-56*f3+11*f4)/(12*d*d)
  f252(g1, f0, f1, f2, f3, d) = (11*g1-20*f0+6*f1+4*f2-f3)/(12*d*d)
  f253(g2, g1, f0, f1, f2, d) = (-g2+16*g1-30*f0+16*f1-f2)/(12*d*d)
  f254(g3, g2, g1, f0, f1, d) = (-g3+4*g2+6*g1-20*f0+11*f1)/(12*d*d)
  f255(g4, g3, g2, g1, f0, d) = (11*g4-56*g3+114*g2-104*g1+35*f0)/(12*d*d)
  f261(f0, f1, f2, f3, f4, f5, d) = (225*f0-770*f1+1070*f2-780*f3+305*f4-50*f5 &
    )/(60*d*d)
  f262(g1, f0, f1, f2, f3, f4, d) = (50*g1-75*f0-20*f1+70*f2-30*f3+5*f4)/ &
    (60*d*d)
  f263(g2, g1, f0, f1, f2, f3, d) = (-5*g2+80*g1-150*f0+80*f1-5*f2+0*f3)/ &
    (60*d*d)
  f264(g3, g2, g1, f0, f1, f2, d) = (0*g3-5*g2+80*g1-150*f0+80*f1-5*f2)/ &
    (60*d*d)
  f265(g4, g3, g2, g1, f0, f1, d) = (5*g4-30*g3+70*g2-20*g1-75*f0+50*f1)/ &
    (60*d*d)
  f266(g5, g4, g3, g2, g1, f0, d) = (-50*g5+305*g4-780*g3+1070*g2-770*g1+225* &
    f0)/(60*d*d)
!     ..

!.....-----------------------------------------------------------------
  if ((iwr==1) .and. (t_inc%i_write>0)) then
    write (1337, "('igd, etc.', 14 i2, f10.4)") ipw, ipg, ivg, ip9, igd, ixlf, iex, xlf
  end if
  iwr = 0

  ist = ist1
!     write(6,*) 'ndvpt ist mesh dx drdi2' ,ndvpt,ist,mesh,dx,
!    &            drdi2(ist)

  if (ndvpt<3 .or. ndvpt>6) then
    write (6, fmt=120) ndvpt
    stop 18
  end if
!.....
!.....ro: total(core+val)(up+down) charge density.

  do i = ist, mesh
    rou(i) = ro(i)*(zta(i)+1.d0)/2.d0
  end do
!.....
  if (igd<=0) then

    do i = ist, mesh
      drr(i) = 0.d0
      ddrr(i) = 0.d0
      drru(i) = 0.d0
      ddrru(i) = 0.d0
    end do
    go to 110

  end if

  i1 = ist
  i2 = ist + 1
  i3 = ist + 2
  i4 = ist + 3
  i5 = ist + 4
  i6 = ist + 5

!.....drr:d(ro)/dr, ddrr=d(d(ro)/dr)/dr
!c.... drru,ddrru: for up   spin,
!.....

  if (nspin==1) go to 100
!.....
  if (ndvpt==3) then

    drx1 = f131(ro(i1), ro(i2), ro(i3), dx)
    drxu1 = f131(rou(i1), rou(i2), rou(i3), dx)
    drxx1 = f231(ro(i1), ro(i2), ro(i3), dx)
    drxxu1 = f231(rou(i1), rou(i2), rou(i3), dx)

  else if (ndvpt==4) then

    drx1 = f141(ro(i1), ro(i2), ro(i3), ro(i4), dx)
    drxu1 = f141(rou(i1), rou(i2), rou(i3), rou(i4), dx)
    drxx1 = f241(ro(i1), ro(i2), ro(i3), ro(i4), dx)
    drxxu1 = f241(rou(i1), rou(i2), rou(i3), rou(i4), dx)
    drx2 = f142(ro(i1), ro(i2), ro(i3), ro(i4), dx)
    drxu2 = f142(rou(i1), rou(i2), rou(i3), rou(i4), dx)
    drxx2 = f242(ro(i1), ro(i2), ro(i3), ro(i4), dx)
    drxxu2 = f242(rou(i1), rou(i2), rou(i3), rou(i4), dx)

  else if (ndvpt==5) then

    drx1 = f151(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), dx)
    drxu1 = f151(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), dx)
    drxx1 = f251(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), dx)
    drxxu1 = f251(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), dx)
    drx2 = f152(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), dx)
    drxu2 = f152(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), dx)
    drxx2 = f252(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), dx)
    drxxu2 = f252(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), dx)

  else if (ndvpt==6) then

    drx1 = f161(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxu1 = f161(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), rou(i6), dx)
    drxx1 = f261(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxxu1 = f261(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), rou(i6), dx)
    drx2 = f162(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxu2 = f162(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), rou(i6), dx)
    drxx2 = f262(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxxu2 = f262(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), rou(i6), dx)
    drx3 = f163(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxu3 = f163(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), rou(i6), dx)
    drxx3 = f263(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxxu3 = f263(rou(i1), rou(i2), rou(i3), rou(i4), rou(i5), rou(i6), dx)

  end if

  drr(i1) = drx1/drdi(i1)
  ddrr(i1) = (drxx1-drx1*drdi2(i1))/drdi(i1)**2
  drru(i1) = drxu1/drdi(i1)
  ddrru(i1) = (drxxu1-drxu1*drdi2(i1))/drdi(i1)**2

  if (ndvpt>3) then

    drr(i2) = drx2/drdi(i2)
    ddrr(i2) = (drxx2-drx2*drdi2(i2))/drdi(i2)**2
    drru(i2) = drxu2/drdi(i2)
    ddrru(i2) = (drxxu2-drxu2*drdi2(i2))/drdi(i2)**2

    if (ndvpt==6) then
      drr(i3) = drx3/drdi(i3)
      ddrr(i3) = (drxx3-drx3*drdi2(i3))/drdi(i3)**2
      drru(i3) = drxu3/drdi(i3)
      ddrru(i3) = (drxxu3-drxu3*drdi2(i3))/drdi(i3)**2
    end if

  end if

  nred = dble(ndvpt)/2 + .1d0

  do j = nred + ist, mesh - nred

    if (ndvpt==3) then

      drx = f132(ro(j-1), ro(j), ro(j+1), dx)
      drxu = f132(rou(j-1), rou(j), rou(j+1), dx)
      drxx = f232(ro(j-1), ro(j), ro(j+1), dx)
      drxxu = f232(rou(j-1), rou(j), rou(j+1), dx)

    else if (ndvpt==4) then

      drx = f142(ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxu = f142(rou(j-1), rou(j), rou(j+1), rou(j+2), dx)
      drxx = f242(ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxxu = f242(rou(j-1), rou(j), rou(j+1), rou(j+2), dx)

    else if (ndvpt==5) then

      drx = f153(ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxu = f153(rou(j-2), rou(j-1), rou(j), rou(j+1), rou(j+2), dx)
      drxx = f253(ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxxu = f253(rou(j-2), rou(j-1), rou(j), rou(j+1), rou(j+2), dx)

    else if (ndvpt==6) then

      drx = f164(ro(j-3), ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxu = f164(rou(j-3), rou(j-2), rou(j-1), rou(j), rou(j+1), rou(j+2), &
        dx)
      drxx = f264(ro(j-3), ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxxu = f264(rou(j-3), rou(j-2), rou(j-1), rou(j), rou(j+1), rou(j+2), &
        dx)

    end if

    drr(j) = drx/drdi(j)
    ddrr(j) = (drxx-drx*drdi2(j))/drdi(j)**2
    drru(j) = drxu/drdi(j)
    ddrru(j) = (drxxu-drxu*drdi2(j))/drdi(j)**2

  end do
!.....
  if (ndvpt==3) then

    drx0 = f133(ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxu0 = f133(rou(mesh-2), rou(mesh-1), rou(mesh), dx)
    drxx0 = f233(ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxxu0 = f233(rou(mesh-2), rou(mesh-1), rou(mesh), dx)

  else if (ndvpt==4) then

    drx1 = f143(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxu1 = f143(rou(mesh-3), rou(mesh-2), rou(mesh-1), rou(mesh), dx)
    drxx1 = f243(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxxu1 = f243(rou(mesh-3), rou(mesh-2), rou(mesh-1), rou(mesh), dx)
    drx0 = f144(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxu0 = f144(rou(mesh-3), rou(mesh-2), rou(mesh-1), rou(mesh), dx)
    drxx0 = f244(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxxu0 = f244(rou(mesh-3), rou(mesh-2), rou(mesh-1), rou(mesh), dx)

  else if (ndvpt==5) then

    drx1 = f154(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxu1 = f154(rou(mesh-4), rou(mesh-3), rou(mesh-2), rou(mesh-1), &
      rou(mesh), dx)
    drxx1 = f254(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxxu1 = f254(rou(mesh-4), rou(mesh-3), rou(mesh-2), rou(mesh-1), &
      rou(mesh), dx)
    drx0 = f155(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxu0 = f155(rou(mesh-4), rou(mesh-3), rou(mesh-2), rou(mesh-1), &
      rou(mesh), dx)
    drxx0 = f255(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxxu0 = f255(rou(mesh-4), rou(mesh-3), rou(mesh-2), rou(mesh-1), &
      rou(mesh), dx)

  else if (ndvpt==6) then

    drx2 = f164(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxu2 = f164(rou(mesh-5), rou(mesh-4), rou(mesh-3), rou(mesh-2), &
      rou(mesh-1), rou(mesh), dx)
    drxx2 = f264(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxxu2 = f264(rou(mesh-5), rou(mesh-4), rou(mesh-3), rou(mesh-2), &
      rou(mesh-1), rou(mesh), dx)

    drx1 = f165(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxu1 = f165(rou(mesh-5), rou(mesh-4), rou(mesh-3), rou(mesh-2), &
      rou(mesh-1), rou(mesh), dx)
    drxx1 = f265(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxxu1 = f265(rou(mesh-5), rou(mesh-4), rou(mesh-3), rou(mesh-2), &
      rou(mesh-1), rou(mesh), dx)

    drx0 = f166(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxu0 = f166(rou(mesh-5), rou(mesh-4), rou(mesh-3), rou(mesh-2), &
      rou(mesh-1), rou(mesh), dx)
    drxx0 = f266(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxxu0 = f266(rou(mesh-5), rou(mesh-4), rou(mesh-3), rou(mesh-2), &
      rou(mesh-1), rou(mesh), dx)


  end if

  if (ndvpt>3) then

    if (ndvpt==6) then
      drr(mesh-2) = drx2/drdi(mesh-2)
      drru(mesh-2) = drxu2/drdi(mesh-2)
      ddrr(mesh-2) = (drxx2-drx2*drdi2(mesh-2))/drdi(mesh-2)**2
      ddrru(mesh-2) = (drxxu2-drxu2*drdi2(mesh-2))/drdi(mesh-2)**2
    end if

    drr(mesh-1) = drx1/drdi(mesh-1)
    drru(mesh-1) = drxu1/drdi(mesh-1)
    ddrr(mesh-1) = (drxx1-drx1*drdi2(mesh-1))/drdi(mesh-1)**2
    ddrru(mesh-1) = (drxxu1-drxu1*drdi2(mesh-1))/drdi(mesh-1)**2

  end if

  drr(mesh) = drx0/drdi(mesh)
  drru(mesh) = drxu0/drdi(mesh)
  ddrr(mesh) = (drxx0-drx0*drdi2(mesh))/drdi(mesh)**2
  ddrru(mesh) = (drxxu0-drxu0*drdi2(mesh))/drdi(mesh)**2

  go to 110

100 continue

!.....
  if (ndvpt==3) then

    drx1 = f131(ro(i1), ro(i2), ro(i3), dx)
    drxx1 = f231(ro(i1), ro(i2), ro(i3), dx)

  else if (ndvpt==4) then

    drx1 = f141(ro(i1), ro(i2), ro(i3), ro(i4), dx)
    drxx1 = f241(ro(i1), ro(i2), ro(i3), ro(i4), dx)
    drx2 = f142(ro(i1), ro(i2), ro(i3), ro(i4), dx)
    drxx2 = f242(ro(i1), ro(i2), ro(i3), ro(i4), dx)

  else if (ndvpt==5) then

    drx1 = f151(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), dx)
    drxx1 = f251(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), dx)
    drx2 = f152(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), dx)
    drxx2 = f252(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), dx)

  else if (ndvpt==6) then

    drx1 = f161(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxx1 = f261(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drx2 = f162(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxx2 = f262(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drx3 = f163(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)
    drxx3 = f263(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6), dx)

  end if

  drr(i1) = drx1/drdi(i1)
  ddrr(i1) = (drxx1-drx1*drdi2(i1))/drdi(i1)**2

  if (ndvpt>3) then

    drr(i2) = drx2/drdi(i2)
    ddrr(i2) = (drxx2-drx2*drdi2(i2))/drdi(i2)**2

    if (ndvpt==6) then
      drr(i3) = drx3/drdi(i3)
      ddrr(i3) = (drxx3-drx3*drdi2(i3))/drdi(i3)**2
    end if

  end if

  nred = dble(ndvpt)/2 + .1d0

  if (mesh-nred<=ist) then
    write (6, fmt='(/'' MESH-NRED.LT.IST. MESH,NRED,IST='',3I4)') mesh, nred, &
      ist
    stop 13
  end if

  do j = nred + ist, mesh - nred

    if (ndvpt==3) then

      drx = f132(ro(j-1), ro(j), ro(j+1), dx)
      drxx = f232(ro(j-1), ro(j), ro(j+1), dx)

    else if (ndvpt==4) then

      drx = f142(ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxx = f242(ro(j-1), ro(j), ro(j+1), ro(j+2), dx)

    else if (ndvpt==5) then

      drx = f153(ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxx = f253(ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2), dx)

    else if (ndvpt==6) then

      drx = f164(ro(j-3), ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2), dx)
      drxx = f264(ro(j-3), ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2), dx)

    end if

    drr(j) = drx/drdi(j)
    ddrr(j) = (drxx-drx*drdi2(j))/drdi(j)**2
!           write(6,9000) j,drr(j)
!9000       format(1x,' j drr(j)',i5,e15.5)
  end do
!.....
  if (ndvpt==3) then

    drx0 = f133(ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxx0 = f233(ro(mesh-2), ro(mesh-1), ro(mesh), dx)

  else if (ndvpt==4) then

    drx1 = f143(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxx1 = f243(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drx0 = f144(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxx0 = f244(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)

  else if (ndvpt==5) then

    drx1 = f154(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxx1 = f254(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drx0 = f155(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)
    drxx0 = f255(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh), dx)

  else if (ndvpt==6) then

    drx2 = f164(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxx2 = f264(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)

    drx1 = f165(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxx1 = f265(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)

    drx0 = f166(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)
    drxx0 = f266(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1), &
      ro(mesh), dx)


  end if

  if (ndvpt>3) then

    if (ndvpt==6) then
      drr(mesh-2) = drx2/drdi(mesh-2)
      ddrr(mesh-2) = (drxx2-drx2*drdi2(mesh-2))/drdi(mesh-2)**2
    end if

    drr(mesh-1) = drx1/drdi(mesh-1)
    ddrr(mesh-1) = (drxx1-drx1*drdi2(mesh-1))/drdi(mesh-1)**2

  end if

  drr(mesh) = drx0/drdi(mesh)
  ddrr(mesh) = (drxx0-drx0*drdi2(mesh))/drdi(mesh)**2

110 continue



!      write(6,8000) nspin,ist1,mesh,dx
!8000 format(1x,' nspin ist1 mesh dx',3i5,2d20.10)
!     write(6,8001) (ro(kk),drr(kk),ddrr(kk),
!    &  drdi(kk),drdi2(kk), kk=ist1,mesh,20)
!8001 format(1x,' ro drr ddrr drdi drdi2',5f12.5)
  return
120 format (/, ' ndvpt should be ge.4 .or. le.6. ndvpt=', i3)
end subroutine
