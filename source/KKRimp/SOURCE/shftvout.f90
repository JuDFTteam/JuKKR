MODULE MOD_SHFTVOUT

  CONTAINS

  !-------------------------------------------------------------------------------
  !> Summary: Shift outer potntial to new expansion around relaxed expansion center
  !> Author: 
  !> Category: KKRimp, potential, radial-grid
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Transform the 'outer' potential from the old expansion around
  !> the ideal positions to the new expansion around the shifted
  !> positions in case of lattice relaxations.
  !> The 'outer' potential is the madelung potential of the ideal
  !> host minus the intracell potential of the host for the perturbed
  !> cluster.
  !>
  !> @note
  !> gaunt2 has to be called bevor to set up the common block assleg
  !> @endnote
  !> @note PR: The previous note is probably not valid anymore. @endnote
  !>
  !> @warning ncleb is an empirical factor - it has to be optimized @endwarning
  !-------------------------------------------------------------------------------
  subroutine shftvout(vin, vout, sn, lpot, wg, yrg, lmmaxd, lassld, lpotd, lm3d, ncleb)
    use mod_ymy, only: ymy
    implicit none
    integer lpotd
    integer lmmaxd, lm3d
    integer lassld
    integer ncleb
    !     .. scalar arguments ..
    integer lpot
    !     .. array arguments ..
    real*8 sn(3),vin(*),vout(*)
    !     .. arrays in common ..
    real*8  wg(lassld),yrg(lassld,0:lassld,0:lassld)
    !     .. local scalars ..
    real*8 clecg,epi,factor,fpi,pi,r,r1,r2,r3,s,r0
    integer i,iend,j,l,l1,l2,l3,l3max,lm1,lm2,lm3,lmmax,lx,ly,m,m1,m1a,m1s,m2,m2a,m2s,m3,m3a,m3s,lm
    !     .. local arrays ..
    real*8 cleb(ncleb),dfac(0:lpotd,0:lpotd),y(lm3d)
    integer icleb(ncleb,3),loflm(lm3d)
    !     .. intrinsic functions ..
    intrinsic abs,real,sign


    pi = 4.d0*datan(1.d0)
  
    !---> determine the l-value for given lm
    i = 1
    do l = 0,2*lpot
      do m = -l,l
        loflm(i) = l
        i = i + 1
      end do
    end do
    do lm=1,lmmaxd
      vout(lm) = 0.d0
    end do
    
    fpi = 4.0d0*pi
    epi = 8.0d0*pi
    l3max = 2*lpot
    lmmax = (lpot+1)**2
    
    !--->calculate:                  (2*(l+l')+1)!!
    !                dfac(l,l')= ----------------------
    !                            (2*l+1)!! * (2*l'+1)!!
    dfac(0,0) = 1.d0
    do lx = 1,lpot
      dfac(lx,0) = dfac(lx-1,0)
      dfac(0,lx) = dfac(lx,0)
      do ly = 1,lx
        dfac(lx,ly) = dfac(lx,ly-1)*real(2* (lx+ly)+1)/real(2*ly+1)
        dfac(ly,lx) = dfac(lx,ly)
      end do
    end do

    !---> set up of the gaunt coefficients with an index field
    !     recognize that they are needed here only for l3=l1+l2
    i = 1
      do l1 = 0,lpot
        do l2 = 0,l1
        l3 = l1 - l2
        do m1 = -l1,l1
          do m2 = -l2,l2
            do m3 = -l3,l3
              m1s = sign(1,m1)
              m2s = sign(1,m2)
              m3s = sign(1,m3)
              
              if (m1s*m2s*m3s.ge.0) then
              
                m1a = abs(m1)
                m2a = abs(m2)
                m3a = abs(m3)
                
                factor = 0.0d0
                
                if (m1a+m2a.eq.m3a) factor = factor + real(3*m3s+sign(1,-m3))/8.0d0
                if (m1a-m2a.eq.m3a) factor = factor + real(m1s)/4.0d0
                if (m2a-m1a.eq.m3a) factor = factor + real(m2s)/4.0d0
                
                if (factor.ne.0.0d0) then
                
                  if (m1s*m2s.ne.1 .or. m2s*m3s.ne.1 .or. m1s*m3s.ne.1) factor = -factor
                  
                  s = 0.0d0
                  do j = 1,lassld
                    s = s + wg(j)*yrg(j,l1,m1a)*yrg(j,l2,m2a)*yrg(j,l3,m3a)
                  end do
                  clecg = s*factor
                  if (abs(clecg).gt.1.d-10) then
                    cleb(i) = clecg
                    icleb(i,1) = l1* (l1+1) + m1 + 1
                    icleb(i,2) = l2* (l2+1) + m2 + 1
                    icleb(i,3) = l3* (l3+1) + m3 + 1
                    i = i + 1
                  end if

                end if

              end if

            end do
          end do
        end do
      end do
    end do
    iend = i - 1


    if (ncleb.lt.iend) then

      stop 13

    else

      do lm1 = 1,lmmax
        vout(lm1)=0.0d0
      end do

      r1 = sn(1)
      r2 = sn(2)
      r3 = sn(3)
      
      r0 = r1**2+r2**2 +r3**2
      if (r0.gt.1.d-10) then
        call ymy(r1,r2,r3,r,y,l3max)
      else
        do lm=1,lmmax
          vout(lm) = vin(lm)
        end do
        return
      end if
      
      do i = 1,iend
        lm1 = icleb(i,1)
        lm2 = icleb(i,2)
        lm3 = icleb(i,3)
        l1 = loflm(lm1)
        l2 = loflm(lm2)
        l3 = loflm(lm3)

        vout(lm2) = vout(lm2) + fpi*(-1.d0)**l3*dfac(l2,l3)*cleb(i)*vin(lm1)*r**l3*y(lm3)
      end do

    end if

  end subroutine shftvout
end module mod_shftvout
