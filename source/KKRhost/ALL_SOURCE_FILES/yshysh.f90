!-------------------------------------------------------------------------------
subroutine yshysh(x, y, z, r, yrealy)
  integer :: lmax
  parameter (lmax=lmaxd)
  integer :: lmax2
  parameter (lmax2=2*lmax)
  integer :: lmax2p, lmxp
  parameter (lmax2p=lmax2+1, lmxp=lmax2p*(lmax2p+1)/2)
!     ..
!     .. Scalar Arguments ..
  double precision :: r, x, y, z
!     ..
!     .. Array Arguments ..
  double precision :: yrealy(*)
!     ..
!     .. Local Scalars ..
  double precision :: a, arg, b, c, cosphi, costhe, eat, p, psq, rssq, save, &
    sinphi, sinthe, tave, tent, w, xa, ya
  integer :: i, ia, ib, ic, istopz, isuzy, j, jc, k, kov2, l, lave, lsuzy, &
    ltwoqp, m, mave, mp1, msuzy, n
!     ..
!     .. Local Arrays ..
  double precision :: cosmph(lmax2p), factor(50), plm(lmxp), sinmph(lmax2p)
!     ..
!     .. Intrinsic Functions ..
  intrinsic :: dabs, dsign, dsqrt
!     ..
  factor(1) = 1.0d00
  do i = 2, 50
    xa = i - 1
    factor(i) = xa*factor(i-1)
  end do
  psq = x*x + y*y
  rssq = psq + z*z
  if (rssq-1.0d-10) 100, 100, 110
100 sinthe = 0.0d00
  costhe = 1.0d00
  cosphi = 1.0d00
  sinphi = 0.0d00
  r = 0.0d00
  go to 140

110 if (psq-1.0d-10) 120, 120, 130
120 r = dsqrt(rssq)
  sinthe = 0.0d00
  costhe = dsign(1.0d00, z)
  sinphi = 0.0d00
  cosphi = 1.0d00
  go to 140

130 r = dsqrt(rssq)
  p = dsqrt(psq)
  sinthe = p/r
  costhe = z/r
  sinphi = y/p
  cosphi = x/p
!
140 xa = dabs(costhe)
  ya = dabs(sinthe)
!      write(6,*) costhe,sinthe
!
  if (xa-1.0d-08) 150, 150, 230
150 l = 0
  j = 0
  tave = 1.0d00
  lsuzy = 1
160 m = 0
  msuzy = 1
170 j = j + 1
  isuzy = lsuzy*msuzy
  if (isuzy>0) go to 180
  plm(j) = 0.0d00
  go to 190

180 k = l + m
  kov2 = k/2
  ia = k + 1
  ib = kov2 + 1
  jc = kov2 - m
  ic = jc + 1
  plm(j) = (((-1)**jc)*factor(ia))/(tave*factor(ib)*factor(ic))
190 if (m-l) 200, 210, 210
200 m = m + 1
  msuzy = -msuzy
  go to 170

210 if (l-lmax2) 220, 450, 450
220 l = l + 1
  lsuzy = -lsuzy
  tave = 2.0d00*tave
  go to 160
!
230 if (xa-0.99999999d00) 330, 240, 240
240 plm(1) = 1.0d00
  plm(2) = costhe
  l = 2
  j = 2
250 j = j + l
  a = l
  ltwoqp = l + l
  b = ltwoqp - 1
  c = l - 1
  k = j - l
  m = j - ltwoqp + 1
  plm(j) = (b*costhe*plm(k)-c*plm(m))/a
  if (l-lmax2) 260, 270, 270
260 l = l + 1
  go to 250
!
270 l = 1
  lave = 1
280 m = 1
  lave = lave + l
290 j = lave + m
  plm(j) = 0.0d00
  if (m-l) 300, 310, 310
300 m = m + 1
  go to 290

310 if (l-lmax2) 320, 450, 450
320 l = l + 1
  go to 280
  c
330 tent = (2.0d00*costhe)/ya
  plm(1) = 1.0d00
  plm(2) = costhe
  plm(3) = ya
  plm(5) = 3.0d00*ya*costhe
  l = 2
  j = 2
340 j = j + l
  a = l
  ltwoqp = l + l
  b = ltwoqp - 1
  c = l - 1
  k = j - l
  m = j - ltwoqp + 1
  plm(j) = (b*costhe*plm(k)-c*plm(m))/a
  if (l-lmax2) 350, 360, 360
350 l = l + 1
  go to 340

360 l = 3
  j = 5
370 j = j + l
  a = l - 1
  ltwoqp = l + l
  b = ltwoqp - 1
  c = l
  k = j - l
  m = j - ltwoqp + 1
  plm(j) = (b*costhe*plm(k)-c*plm(m))/a
  if (l-lmax2) 380, 390, 390
380 l = l + 1
  go to 370

390 l = 2
  lave = 3
400 lave = lave + l
  m = 1
410 j = lave + m
  k = j - 1
  n = k - 1
  eat = m
  a = tent*eat
  b = (m+l)*(l-m+1)
  plm(j) = a*plm(k) - b*plm(n)
  if (m+1-l) 420, 430, 430
420 m = m + 1
  go to 410

430 if (l-lmax2) 440, 450, 450
440 l = l + 1
  go to 400
  c
450 sinmph(1) = 0.0d00
  cosmph(1) = 1.0d00
  istopz = lmax2 + 1
  do i = 2, istopz
    j = i - 1
    sinmph(i) = sinphi*cosmph(j) + cosphi*sinmph(j)
    cosmph(i) = cosphi*cosmph(j) - sinphi*sinmph(j)
  end do
  c
  l = 0
460 m = 0
  lave = l*(l+1) + 1
  mave = ((l*(l+1))/2) + 1
  save = 2*l + 1
470 if (m/=0) go to 480
  arg = save/12.5663706144d00
  w = dsqrt(arg)
  yrealy(lave) = w*plm(mave)
  go to 490

480 ia = l - m + 1
  ib = l + m + 1
  arg = (save*factor(ia))/(6.28318530718d00*factor(ib))
  mp1 = m + 1
  w = dsqrt(arg)
  i = lave + m
  j = mave + m
  yrealy(i) = w*plm(j)*cosmph(mp1)
  i = lave - m
  yrealy(i) = w*plm(j)*sinmph(mp1)
490 if (m>=l) go to 500
  m = m + 1
  go to 470

500 if (l>=lmax2) go to 510
  l = l + 1
  go to 460
  c
510 return
  c
  c
  c
end subroutine
