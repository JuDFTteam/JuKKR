subroutine pnsqns(ar,cr,dr,drdi,ek,icst,pz,qz,fz,sz,pns,qns,nsra,fred, &
vins,ipan,ircut,cleb,icleb,iend,loflm,lkonv, &
ispin,ldau,nldau,lldau, &
wmldau,wmldauav,ldaucut, &
lmaxd, nspind, irmd, irnsd, ipand, ncleb)

  implicit none

  integer lmaxd
  integer nspind
  integer irmd
  integer irnsd
  integer ipand
  integer ncleb
  !     ..
  !     .. Scalar Arguments ..
  double complex     ek
  integer            icst,iend,ipan,lkonv,nsra,nldau,ispin,fred
  logical            ldau
  !     ..
  !     .. Array Arguments ..

  double complex     ar((lmaxd+1)**2,(lmaxd+1)**2)
  double complex     cr((lmaxd+1)**2,(lmaxd+1)**2)
  double complex     fz(irmd,0:lmaxd)

  !     DOUBLE COMPLEX     PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
  double complex     pns((lmaxd+1)**2,(lmaxd+1)**2, &
  (irmd-irnsd):irmd,2)

  double complex     pz(irmd,0:lmaxd)
  !     DOUBLE COMPLEX     QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
  double complex     qns((lmaxd+1)**2,(lmaxd+1)**2, &
  (irmd-irnsd):irmd,2)

  double complex     qz(irmd,0:lmaxd)
  double complex     sz(irmd,0:lmaxd)

  double precision   cleb(ncleb,2)
  double precision   drdi(irmd)
  !     DOUBLE PRECISION   VINS(IRMIND:IRMD,LMPOTD)
  double precision   vins(irmd-irnsd:irmd,(2*lmaxd+1)**2)
  !     DOUBLE PRECISION   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  double precision   wmldau(2*lmaxd+1,2*lmaxd+1,lmaxd+1,nspind)
  double precision   ldaucut(irmd)

  integer            icleb(ncleb,3),ircut(0:ipand),loflm(*)
  integer            lldau(lmaxd+1)
  !     ..
  !     .. Local Scalars ..
  integer            i,ir,irc1,lm1,lm2,lmmkonv,m1,m2, &
  lmlo,lmhi,ildau
  !     ..
  !     .. Local Arrays ..
  !     DOUBLE COMPLEX     CMAT(LMMAXD,LMMAXD,IRMIND:IRMD)
  double complex     cmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
  !     DOUBLE COMPLEX     DMAT(LMMAXD,LMMAXD,IRMIND:IRMD)
  double complex     dmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
  !     DOUBLE COMPLEX     DR(LMMAXD,LMMAXD)
  double complex     dr((lmaxd+1)**2, (lmaxd+1)**2)
  !     DOUBLE COMPLEX     EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2)
  double complex     efac((lmaxd+1)**2)
  double complex     pzekdr((lmaxd+1)**2,irmd-irnsd:irmd,2)
  !     DOUBLE COMPLEX     PZLM(LMMAXD,IRMIND:IRMD,2)
  double complex     pzlm((lmaxd+1)**2,(irmd-irnsd):irmd,2)
  !     DOUBLE COMPLEX     QZEKDR(LMMAXD,IRMIND:IRMD,2)
  double complex     qzekdr((lmaxd+1)**2,(irmd-irnsd):irmd,2)
  !     DOUBLE COMPLEX     QZLM(LMMAXD,IRMIND:IRMD,2)
  double complex     qzlm((lmaxd+1)**2, (irmd-irnsd):irmd,2)
  !     DOUBLE COMPLEX     TMATLL(LMMAXD,LMMAXD)
  double complex     tmatll((lmaxd+1)**2, (lmaxd+1)**2)
  !     DOUBLE PRECISION   VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
  double precision   vnspll((lmaxd+1)**2, (lmaxd+1)**2, &
  irmd-irnsd:irmd)

  double precision   wmldauav(lmaxd+1)
  !     ..
  !     .. External Subroutines ..
  external irwns,regns,vllns,wftsca

  integer             lmmaxd
  integer             irmind

  irmind= irmd-irnsd
  lmmaxd= (lmaxd+1)**2

  irc1 = ircut(ipan)
  !

  ! calculate V_LL' potential
  ! V_{LL'} = \sum_{L''} C_{L L' L''} V_{L''}
  call vllns(vnspll,vins,cleb,icleb,iend, &
  lmaxd, irmd, irnsd, ncleb)


  if (lkonv.ne.lmaxd) then
    lmmkonv = (lkonv+1)* (lkonv+1)
    do lm1 = 1,lmmaxd
      do lm2 = lmmkonv + 1,lmmaxd
        do i = irmind,irmd
          vnspll(lm2,lm1,i) = 0.0d0
          vnspll(lm1,lm2,i) = 0.0d0
        end do
      end do
    end do
  else
    lmmkonv = lmmaxd
  end if

  !-----------------------------------------------------------------------
  ! LDA+U
  ! Add WLDAU to non-spherical porential VINS in case of LDA+U
  ! Use the average wldau (=wldauav) and the deviation
  ! of wldau from this. Use the deviation in the Born series
  ! for the non-spherical wavefunction, while the average is
  ! used for the spherical wavefunction.
  !
  if (ldau) then
    do ildau=1,nldau
      if (lldau(ildau).ge.0) then
        !
        lmlo = lldau(ildau)*lldau(ildau) + 1
        lmhi = (lldau(ildau)+1)*(lldau(ildau)+1)
        !
        do ir = irmind,irmd
          !
          ! -> First add wldau to all elements ...
          !
          do lm2 = lmlo,lmhi
            m2 = lm2 - lmlo + 1
            do lm1 = lmlo,lmhi
              m1 = lm1 - lmlo + 1
              vnspll(lm1,lm2,ir) =vnspll(lm1,lm2,ir) &
              +wmldau(m1,m2,ildau,ispin)*ldaucut(ir)
            enddo
            !
            ! ... and then subtract average from diag. elements
            !
            vnspll(lm2,lm2,ir) =  vnspll(lm2,lm2,ir) &
            - wmldauav(ildau) * ldaucut(ir)
          enddo
        enddo
      endif
    enddo
  endif
  !
  ! LDA+U
  !-----------------------------------------------------------------------


  !
  !---> get wfts of same magnitude by scaling with efac
  !
  call wftsca(drdi,efac,pz,qz,fz,sz,nsra,pzlm,qzlm,pzekdr,qzekdr, &
  ek,loflm,irmind,irmd,lmaxd,lmmaxd)
  !
  !---> determine the irregular non sph. wavefunction
  !
  call irwns(cr,dr,efac,qns,vnspll,icst,ipan,ircut,nsra,pzlm,qzlm, &
  pzekdr,qzekdr,qns(1,1,irmind,1),cmat, &
  qns(1,1,irmind,2),dmat,irmind,irmd,ipand,lmmaxd)
  !
  !---> determine the regular non sph. wavefunction
  !
  call regns_old(ar,tmatll,efac,pns,vnspll,icst,ipan,ircut,pzlm,qzlm, &
  pzekdr,qzekdr,ek,pns(1,1,irmind,1),cmat, &
  pns(1,1,irmind,2),dmat,nsra,irmind,irmd,ipand,lmmaxd)
  !

  return

end
