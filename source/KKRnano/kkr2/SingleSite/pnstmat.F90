!>    construct the t-matrix from radial solutions and potential
      subroutine pnstmat(drdi,ek,icst,pz,qz,fz,sz,pns,tmatll,vins, &
                         ipan,ircut,nsra,cleb,icleb,iend,loflm,tmat, &
                         det,lkonv, &
                         ldau,nldau,lldau, &
                         wmldau_ispin,wmldauav,ldaucut, lmaxd, irmd, irnsd, ipand, ncleb)
      implicit none

      integer lmaxd
      integer irmd
      integer irnsd
      integer ipand
      integer ncleb
      double complex     ek,det
      integer            icst,iend,ipan,lkonv,nsra,nldau
      logical            ldau
      double complex    fz(irmd,0:lmaxd)
      double complex    pns((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd,2)
      double complex    pz(irmd,0:lmaxd)
      double complex    qz(irmd,0:lmaxd)
      double complex    sz(irmd,0:lmaxd)
      double complex    tmat(0:lmaxd)
      double complex    tmatll((lmaxd+1)**2, (lmaxd+1)**2)

      double precision   cleb(ncleb,2)
      double precision   drdi(irmd)
      double precision   vins(irmd-irnsd:irmd,(2*lmaxd+1)**2)

      double precision   wmldau_ispin(2*lmaxd+1,2*lmaxd+1,lmaxd+1)
      double precision   ldaucut(irmd)
      double precision   wmldauav(lmaxd+1)
      integer            icleb(ncleb,3),ircut(0:ipand),loflm(*),lldau(lmaxd+1)
     
      external :: regns, vllns, wftsca
      double complex, parameter :: zero=(0.d0,0.d0)
      integer            i,irc1,lm1,lm2,lmmkonv,ir,m1,m2,lmlo,lmhi,ildau

      double complex     ar((lmaxd+1)**2,(lmaxd+1)**2)
      double complex     cmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
      double complex     dmat((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)
      double complex     efac((lmaxd+1)**2)
      double complex     pzekdr((lmaxd+1)**2,irmd-irnsd:irmd,2)
      double complex     pzlm((lmaxd+1)**2, irmd-irnsd:irmd,2)
      double complex     qzekdr((lmaxd+1)**2,irmd-irnsd:irmd,2)
      double complex     qzlm((lmaxd+1)**2,irmd-irnsd:irmd,2)
      double precision   vnspll((lmaxd+1)**2,(lmaxd+1)**2,irmd-irnsd:irmd)

      integer            ipvt((lmaxd+1)**2)
      integer            info
      integer             irmind
      integer             lmmaxd

!     initialisation of some local arrays
      pzlm = zero
      qzlm = zero
      pzekdr = zero
      qzekdr = zero
      ar = zero
!     end initialisation

      irmind = irmd-irnsd
      lmmaxd = (lmaxd+1)**2

      irc1 = ircut(ipan)

      call vllns(vnspll,vins,cleb,icleb,iend, lmaxd, irmd, irnsd, ncleb)

      if (lkonv /= lmaxd) then
        lmmkonv = (lkonv+1)*(lkonv+1)
        do lm1 = 1, lmmaxd
          do lm2 = lmmkonv + 1, lmmaxd
            do i = irmind, irmd
              vnspll(lm2,lm1,i) = 0.d0
              vnspll(lm1,lm2,i) = 0.d0
            end do
          end do
        end do
      else
        lmmkonv = lmmaxd
      end if


!-----------------------------------------------------------------------
! lda+u
! add wldau to non-spherical porential vins in case of lda+u
! use the average wldau (=wldauav) and calculate the deviation
! of wldau from this. use the deviation in the born series
! for the non-spherical wavefunction, while the average is
! used for the spherical wavefunction.
!
      if (ldau) then
      do ildau=1,nldau
      if (lldau(ildau) >= 0) then

        lmlo = lldau(ildau)*lldau(ildau) + 1
        lmhi = (lldau(ildau)+1)*(lldau(ildau)+1)

        do ir = irmind,irmd
!
! -> first add wldau to all elements ...
!
          do lm2 = lmlo,lmhi
            m2 = lm2 - lmlo + 1
            do lm1 = lmlo,lmhi
              m1 = lm1 - lmlo + 1
              vnspll(lm1,lm2,ir) = vnspll(lm1,lm2,ir) + wmldau_ispin(m1,m2,ildau)*ldaucut(ir)
            enddo
!
! ... and then subtract average from diag. elements
!
            vnspll(lm2,lm2,ir) =  vnspll(lm2,lm2,ir) - wmldauav(ildau) * ldaucut(ir)
          enddo
        enddo
      endif
      enddo
      endif
!
! lda+u
!-----------------------------------------------------------------------

!
!---> get wfts of same magnitude by scaling with efac
!
      call wftsca(drdi,efac,pz,qz,fz,sz,nsra,pzlm,qzlm,pzekdr,qzekdr, ek,loflm,irmind,irmd,lmaxd,lmmaxd)
!
!---> determine the regular non sph. wavefunction
!
      call regns(ar,tmatll,efac,pns,vnspll,icst,ipan,ircut,pzlm,qzlm, &
                   pzekdr,qzekdr,ek,pns(1,1,irmind,1),cmat, &
                   pns(1,1,irmind,2),dmat,nsra,irmind,irmd,ipand,lmmaxd)

      do lm1 = 1,lmmkonv
        tmatll(lm1,lm1) = tmatll(lm1,lm1) + tmat(loflm(lm1))
      end do
        det = (1.d0,0.d0)
        call zgetrf(lmmaxd,lmmaxd,ar,lmmaxd,ipvt,info)
        do lm1 = 1,lmmaxd
          if (ipvt(lm1) /= lm1) det = -det
          det = ar(lm1,lm1)*det
        enddo

      end
