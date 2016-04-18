module AtomicCore_mod
  implicit none
  private
  public :: rhocore
  
  
  double precision, parameter, private :: r83sq=64.d0/9.d0, r1=1.d0/9.d0, r2=-5.d0*r1, r3=19.d0*r1
  double precision, parameter, private :: light_speed(1:2) = [1.d0, 274.0720442d0]
  
  contains
  
  
    double precision function rhocore(ebot, nsra, nspin, atom_id, drdi, r, visp, a, b, zat, ircut, rhoc, ecore, ncore, lcore, irmd) result(qc)
      integer, intent(in) :: irmd
      double precision, intent(in) :: a, b, zat, ebot
      integer, intent(in) :: atom_id, nspin, nsra, ncore(:) ! (nspin)
      double precision, intent(in) :: drdi(:), r(:), visp(:,:) ! (irmd,nspin)
      double precision, intent(inout) :: ecore(:,:) ! (20,nspin)
      double precision, intent(out) :: rhoc(:,:)! (irmd,nspin) 
      integer, intent(in) :: lcore(:,:), ircut(0:) ! (20,nspin) (0:ipand)

!     .. locals ..
      double precision :: rmax
      integer :: ipr, nr, ispin

! --------------------------------------------------------------
!     ipr=0 : do not write state dependent information
!     ipr=1 : write something
!     ipr=2 : write all (for debugging)
! --------------------------------------------------------------
      ipr = 0
      qc = 0.d0
 
      nr = ircut(1) ! end of core region
      rmax = r(nr)
      
      do ispin = 1, nspin

        qc = qc + corel(nsra,ipr,atom_id,rhoc(:,ispin),visp(:,ispin),ecore(:,ispin),lcore(:,ispin),ncore(ispin),drdi,zat,a,b,ispin,nspin,nr,rmax,irmd,ebot)

      enddo ! ispin
      if (ipr /= 0) write (6,fmt="(1x,5('*'),' core-relaxation for ',i3,'th cell was done ',5('*'))") atom_id
      
    endfunction rhocore

  
     double precision function corel(nsra, ipr, atom_id, rhoc, v, ecore, lcore, ncore, drdi, z, a, b, is, nspin, nr, rmax, irmd, ebot) result(qc)
!-----------------------------------------------------------------------
!     driver for core states
!-----------------------------------------------------------------------
!     lmxc = lmaxcore = (0,1,2,...), .e.g, argon core : lmxc = 1
!                                        krypton core : lmxc = 2
!     kfg = configuration of core, e.g., argon core: 3300=3s,3p,0d
!                                      krypton core: 4430=4s,4p,3d
!                                      xenon core: 5540=5s,5p,4d
!-----------------------------------------------------------------------
      use quadrature_mod, only: simpson

      double precision, intent(in) :: a,b,rmax,z,ebot,drdi(:),v(:)
      integer, intent(in) :: atom_id,ipr,irmd,is,ncore,nr,nspin,nsra,lcore(:)
      double precision, intent(inout) :: ecore(:)
      double precision, intent(out) :: rhoc(:)

!     .. locals ..
      integer, parameter :: nitmax=40, irnumx=10
      double precision :: e,e1,e2,ediff,ei,slope,sm,tol,value,wgt
      integer :: ic,in,inuc,ir,l,nc,nn,nre
      logical :: vlnc
      double precision :: f(irmd), g(irmd), rho(irmd)
      integer :: kfg(0:3)
      character(len=*), parameter :: spn(2)=['down','up  '], text(0:4)=['s   ','p   ','d   ','f   ','g   '], &
        F9="(1x,90 ('*'),/,'  n = ',i1,'  l = ',a4,'   nnode = ',i1,'  spin=',a4,i5,'th cell','    einput = ',1p,d16.8)"

      vlnc = .false.
      value =  1.d-8
      slope = -1.d-8
      e2 = 50.d0

      kfg(0:3) = 0
      do ic = 1, ncore
        do l = 0, 3
          if (lcore(ic) == l) kfg(l) = kfg(l) + 1
        enddo ! l
      enddo ! ic
      do l = 0, 3
        if (kfg(l) > 0) kfg(l) = kfg(l) + l
      enddo ! l

      tol = 1.d-12*(z*z + 1.d0)
      nc = 0
      inuc = -irnumx

      rhoc(1:irmd) = 0.d0
      rho(1:irmd) = 0.d0

      do l = 0, 3
        e1 = (-5.d0 - ((z + 1.d0)/(l + 1.d0))**2)*1.5d0 - 50.d0 ! guess value formula
        do in = l + 1, kfg(l)
          nn = in - l - 1
          nc = nc + 1
          inuc = inuc + irnumx
          e = ecore(nc)
          ei = e
          if (ipr /= 0) write(6,fmt=F9) in,text(l),nn,spn(is),atom_id,e

          sm = intcor(e1,e2,rho,g,f,v,value,slope,l,nn,e,nre,vlnc,a,b,z,rmax,nr,tol,irmd,ipr,nitmax,nsra)

          if (e > ebot) then
            write(6,'(''Error for L='',I1)') l
            write(6,*) 'E,EBOT',e,ebot
            write(6,*) 'The program found a core state above the bottom of the valence-band energy contour.'
            write(6,*) 'The results are very probably wrong.'
            write(6,*) 'The number of core states in the input potential should perhaps be decreased.'
            write(6,*) 'The program stops in corel.F90, atom_id=',atom_id
            stop 'Error 1 in corel.F90'
          endif

          ediff = e - ei
          ecore(nc) = e
          wgt = (l+1+l)*2.d0/(sm*nspin)
          if (ipr /= 0) write(6,fmt="(1x,'  einput =',1p,d16.8,'   eout - ein =',1p,d16.8,'   eoutput = ',1p,d16.8)") ei,ediff,e

          ! sum up contributions to total core charge
          do ir = 2, nre
            rhoc(ir) = rhoc(ir) + rho(ir)*wgt
            rho(ir) = 0.d0
          enddo ! ir

        enddo ! in

!     calculate the first valence state with (n_max+1, l) to verify
!     that it lies above the bottom of the energy contour
!     (treating it as an atomic state)
!     start energy: e(n_max, l) / 10.

        if (kfg(l) > 0) then
          in = kfg(l) + 1
          nn = in - l - 1
          e = ecore(nc)*0.1
          ei = e
          if (ipr /= 0) write(6,fmt=F9) in,text(l),nn,spn(is),atom_id,e

          sm = intcor(e1,e2,rho,g,f,v,value,slope,l,nn,e,nre,vlnc,a,b,z,rmax,nr,tol,irmd,ipr,nitmax,nsra)

          if (e < ebot) then
            write(6,'(''Error for L='',I1)') l
            write(6,*) 'E,EBOT',e,ebot
            write(6,*) 'The program found a core state below the bottom of the valence-band energy contour.'
            write(6,*) 'This state was not given in the input potential.'
            write(6,*) 'The results are very probably wrong.'
            write(6,*) 'Lower the bottom of the contour or increase the number of core states in the input potential.'
            write(6,*) 'The program stops in corel.F90, atom_id=',atom_id
            stop 'Error 2 in corel.F90'
          endif
        endif

      enddo ! l
      if (nc*irnumx > 150 .or. irnumx > 10) stop 'corel'

      qc = simpson(rhoc(1:nr), 1, nr, drdi) ! integrate core density to get core charge

    endfunction ! core_electrons


    double precision function intcor(f1, f2, rho, g, f, v, value, slope, l, nn, e, nre, vlnc, a, b, z, rn, nr, tol, irm, ipr, nitmax, nsra) result(sm)
      use Constants_mod, only: pi
    
      double precision, intent(in) :: a, b, f1, f2, rn, slope, tol, value, z
      double precision, intent(inout) :: e
      integer, intent(in) :: ipr, irm, l, nitmax, nn, nr, nsra
      integer, intent(out) :: nre
      logical, intent(in) :: vlnc
      double precision, intent(in) :: v(:)
      double precision, intent(out) :: f(:), g(:), rho(:)

!     .. locals ..
      double precision :: cvlight, de, dg1, dg2, dpsi1, dpsi2, drdikc, e1, e2, ea
      double precision :: gkc2, pkc1, pkc2, psi1, psi2, q, qkc1, qkc2, ratio
      double precision :: re, rkc, rpb, slop, tsrme, valu, vme, xxx
      integer :: k, k2, kc, iiter, nne, run
      double complex :: arg, cappai, dofe, hl(6)

      character(len=*), parameter :: F9020="(2i3, 2i4, 1p, 3d16.8, 1x, 2d9.2)", &
      F9000="(' l=', i3, '  nn=', i2, '  nr=', i4, '  f1/e/f2=', 3f10.3, /, ' tol=', 1p, d12.3, '  value/slope=', 2d12.3)", &
      F9030="(/, ' **** int: 0-pressure bcs not real')", F9050="(' *** int: stop after', i4, ' iterations')", &
      F9010="(13x, '  no boundary condition had to be used')", &
      F9040="(' state', i2, ', ', i1, ':', i4, 'x, ', i5, '/', i3, ',  bc=', 1p, 2d12.3, /, 14x, 'e=', d14.6, '   de=', d11.2, '   sm=', d12.4)"

      cvlight = light_speed(nsra)
      ea = exp(a)
      iiter = 0
      e1 = f1
      e2 = f2
      rho(1:irm) = 0.d0

      if (ipr == 2) write(6, fmt=F9000) l, nn, nr, f1, e, f2, tol, value, slope

      run = 1
      do while (run > 0)

        iiter = iiter + 1
        if (iiter > nitmax) then
          write(6, fmt=F9050) nitmax
          stop 'intcor'
        endif

        g(1:irm) = 0.d0
        f(1:irm) = 0.d0

        if (e <= e1 .or. e >= e2) e = 0.5d0*(e1 + e2)
        nre = nr
  !     write(6, *) e, e1, e2
        if (e > -1.d-8) return

        tsrme = 2.d0*sqrt(-e)
        re = (log(-tsrme*e/1.d-8)/tsrme - (z+z)/e)*2.d0

        nre = log(re/b + 1.d0)/a + 1.d0
        nre = (nre/2)*2 + 1
        nre = min0(nre, nr)
        nre = max0(nre, 35)

        xxx = 1.d0
        valu = 1.d-1
        slop = -1.d-1
        if (nre < nr .and. iiter == 1 .and. ipr /= 0) write(6, fmt=F9010)
        if (nre >= nr) then
          valu = value
          slop = slope
          if (.not. vlnc) then
  !--->     single site boundary condition
            vme = 0.d0 - e
            if (nsra == 1) then
              cappai = dcmplx(0.d0, dsqrt(vme))
            else
              cappai = dcmplx(0.d0, dsqrt((1.d0 - vme/cvlight/cvlight)*vme))
            endif
            arg = cappai*rn
            hl = hankel_core(l+2, arg)
            dofe = (l + 1.d0)/rn - cappai*hl(l+2)/hl(l+1)
            valu = 1.d-10
            slop = valu*dofe
            ! todo: valu and slop are real, maybe this part can be simplified ...
            ! since the Hankel functions are evaluated for purely imaginary arguments
          endif

        endif
        k2 = 30; if (nn == 0) k2 = nre/3
        nne = 0 ! init the number of nodes

        dg2 = intin(g, f, v, e, l, nne, valu, slop, nre, k2, kc, a, b, z, nsra)

        rkc = b*exp(a*kc - a) - b
        drdikc = a*(rkc + b)
        gkc2 = g(kc)
        psi2 = g(kc)
        dpsi2 = dg2/drdikc
        qkc2 = psi2*psi2 + dpsi2*dpsi2*rkc*rkc
        pkc2 = 0.5d0 - atan(rkc*dpsi2/psi2)/pi

        dg1 = intout(g, f, v, e, l, nne, kc, a, b, z, nsra)

        psi1 = g(kc)
        dpsi1 = dg1/drdikc
        qkc1 = psi1*psi1 + dpsi1*dpsi1*rkc*rkc
        pkc1 = 0.5d0 - atan(rkc*dpsi1/psi1)/pi
        if (nne == 9) nne = 0 ! why?
        if (nne == nn) then
  !         ratio1 = gkc2/g(kc)
  !         ratio = sqrt(qkc2/qkc1)
  !         if (ratio1 < 0.d0) ratio = -ratio
          ratio = sign(sqrt(qkc2/qkc1), gkc2*g(kc))
          g(1:kc) = g(1:kc)*ratio
          if (nsra == 1) then
            f(1:nre) = 0.d0
          else
            f(1:kc) = f(1:kc)*ratio
          endif
          sm = 0.d0
          rpb = b/ea
          q = ea*ea
          do k = 2, nre - 1, 2
            rpb = rpb*q
            sm = sm + rpb*(g(k)*g(k) + f(k)*f(k))
          enddo ! k
          rpb = b
          sm = sm + sm
          do k = 3, nre - 2, 2
            rpb = rpb*q
            sm = sm + rpb*(g(k)*g(k) + f(k)*f(k))
          enddo ! k
          sm = 2.d0*sm + rpb*q*(g(nre)*g(nre) + f(nre)*f(nre))
          sm = a*sm/3.d0
          de = pi*qkc2*(pkc2 - pkc1)/sm/rkc
          if (iiter >= nitmax-10 .or. ipr == 2) write(6, fmt=F9020) iiter, nne, nre, kc, e1, e, e2, de
          if (de > 0.d0) e1 = e
          if (de < 0.d0) e2 = e
          e = e + de
          if (abs(de) <= tol .or. iiter >= nitmax) run = 0 ! converged

        else
          if (iiter >= nitmax-10 .or. ipr == 2) write(6, fmt=F9020) iiter, nne, nre, kc, e1, e, e2
          if (nne > nn) e2 = e
          if (nne < nn) e1 = e
          e = 0.5d0*(e1 + e2)

        endif
      enddo ! while(run > 0)

      e = e - de
      do k = 1, nre
        rho(k) = g(k)*g(k) + f(k)*f(k)
      enddo ! k
      if (xxx <= 0.d0) write(6, fmt=F9030)
      if (iiter >= nitmax-10 .or. ipr >= 1 .or. xxx <= 0.d0) &
        write(6, fmt=F9040) l, nn, iiter, kc, nre, valu, slop, e, de, sm
    endfunction intcor
      
      
    double precision function intin(g, f, v, e, l, nne, valu, slop, k1, k2, kc, a, b, z, nsra) result(dg)
      integer, intent(in) :: k1, k2, l, nsra
      integer, intent(inout) :: nne
      integer, intent(out) :: kc
      double precision, intent(in) :: a, b, e, slop, valu, z
      double precision, intent(in) :: v(:)
      double precision, intent(out) :: f(:),  g(:)

!     .. locals ..
      double precision :: af1, af2, af3, ag1, ag2, ag3, b1, b2, cvlight, phiwgt, det, df1
      double precision :: df2, df3, dg1, dg2, dg3, dr, ea, ff, fllp1, gg, phi, q
      double precision :: r, rpb, sdg3, sg, sgp1, u, vb, x, y, zz
      double precision, parameter :: h83=-8.d0/3.d0 ! negative!!
      integer :: i, ir, jr, run
      double precision :: d(2,3)

      zz = z + z
      cvlight = light_speed(nsra); phiwgt = nsra - 1.d0
      fllp1 = l*(l + 1.d0)
      
      ea = exp(a)
      rpb = b*exp(a*k1 - a)
      r = rpb - b
      dr = a*rpb
      phi = (e + zz/r - v(k1))*dr/cvlight
      u = dr*cvlight + phi*phiwgt
      x = -dr/r
      y = -fllp1*x*x/u + phi
      g(k1) = valu
      f(k1) = (slop*dr + x*valu)/u
      q = 1.d0/sqrt(ea)
      ag1 = slop*dr
      af1 = x*f(k1) - y*g(k1)
      ir = k1
      dg3 = ag1
      if (k2 /= k1) then
        do i = 1, 3
          jr = ir
          ir = ir - 1
          rpb = rpb*q
          dr = rpb*a
          r = rpb - b
          gg = g(jr) - 0.5d0*ag1
          ff = f(jr) - 0.5d0*af1
          vb = (3.d0*v(jr) + 6.d0*v(ir) - v(ir-1))*.125d0
          phi = (e + zz/r - vb)*dr/cvlight
          u = dr*cvlight + phi*phiwgt
          x = -dr/r
          y = -fllp1*x*x/u + phi
          ag2 = u*ff - x*gg
          af2 = x*ff - y*gg
          gg = g(jr) - 0.5d0*ag2
          ff = f(jr) - 0.5d0*af2
          ag3 = u*ff - x*gg
          af3 = x*ff - y*gg
          rpb = rpb*q
          dr = a*rpb
          r = rpb - b
          phi = (e + zz/r - v(ir))*dr/cvlight
          u = dr*cvlight + phi
          if (nsra == 1) u = dr
          x = -dr/r
          y = -fllp1*x*x/u + phi
          gg = g(jr) - ag3
          ff = f(jr) - af3
          g(ir) = g(jr) - (ag1 + 2.d0*(ag2 + ag3) + u*ff - x*gg)/6.d0
          f(ir) = f(jr) - (af1 + 2.d0*(af2 + af3) + x*ff - y*gg)/6.d0
          sg   = dsign(1.d0,g(ir))
          sgp1 = dsign(1.d0,g(jr))
          if (sg*sgp1 < 0.d0) nne = nne + 1
          ag1 = u*f(ir) - x*g(ir)
          af1 = x*f(ir) - y*g(ir)
          if (ir == k2) then
            kc = ir
            dg = dg3
            return
          else
            d(1,i) = ag1
            d(2,i) = af1
          endif
        enddo
        q = 1.d0/ea
        dg1 = d(1,1)
        dg2 = d(1,2)
        dg3 = d(1,3)
        df1 = d(2,1)
        df2 = d(2,2)
        df3 = d(2,3)

        run = 1
        do while (run > 0)
          jr = ir
          ir = jr - 1
          rpb = rpb*q
          dr = a*rpb
          r = rpb - b
          phi = (e + zz/r - v(ir))*dr/cvlight
          u = dr*cvlight + phi*phiwgt
          x = -dr/r
          y = -fllp1*x*x/u + phi
          det = r83sq - x*x + u*y
          b1 = g(jr)*h83 + r1*dg1 + r2*dg2 + r3*dg3
          b2 = f(jr)*h83 + r1*df1 + r2*df2 + r3*df3
          g(ir) = (b1*(h83 - x) + b2*u)/det
          f(ir) = (b2*(h83 + x) - b1*y)/det
          sg   = dsign(1.d0,g(ir))
          sgp1 = dsign(1.d0,g(jr))
          if (sg*sgp1 < 0.d0) nne = nne + 1
          dg1 = dg2
          df1 = df2
          dg2 = dg3
          df2 = df3
          dg3 = u*f(ir) - x*g(ir)
          df3 = x*f(ir) - y*g(ir)
          sdg3 = dsign(1.d0,dg3)
          if (ir <= k2 .or. sg*sdg3 >= 0.d0) run = 0 ! stop
        enddo ! while
      endif
      kc = ir
      dg = dg3
    endfunction intin
      
      
    double precision function intout(g, f, v, e, l, nne, k2, a, b, z, nsra) result(dg)
      double precision, intent(in) :: a, b, e, z
      integer, intent(in) :: k2, l, nsra
      integer, intent(inout) :: nne
      double precision, intent(in) :: v(:)
      double precision, intent(out) :: f(:), g(:)

!     .. locals ..
      double precision :: aa,alfa,b1,b2,bb,beta,cvlight, phiwgt, det,df1,df2,df3
      double precision :: dg1,dg2,dg3,dr,ea,fllp1,p12,p21,phi,pp,qq,r
      double precision :: rpb,s,sg,sgm1,u,x,y,zz
      double precision, parameter :: h83=8.d0/3.d0 ! positive!!
      integer :: i,i1,k,km1,n
      double precision :: d(2,3), px(12), qx(12)

      zz = z + z
      cvlight = light_speed(nsra); phiwgt = nsra - 1.d0
      ea = exp(a)
      fllp1 = l*(l + 1.d0)
      aa = -zz/cvlight
      bb = fllp1 - aa*aa
      p21 = (v(1) - e)/cvlight
      p12 = cvlight - p21
      px(1) = 0.d0
      qx(1) = 0.d0
      if (z <= 20.d0 .or. nsra == 1) then

        s = 1.d0*l
        px(2) = 0.d0
        px(3) = 1.d0
        do k = 2, 9
          px(k+2) = ((v(1) - e)*px(k) - zz*px(k+1))/(k + 2.d0*l)/(k - 1.d0)
        enddo ! k
        do k = 2, 10
          qx(k) = px(k+1)*(l + k - 2.d0)/cvlight
        enddo ! k

      else
      
        s = sqrt(fllp1 + 1.d0 - aa*aa)
        px(2) = 1.d0
        qx(2) = (1.d0 - s)/aa
        do i = 3, 10
          n = i - 2
          alfa = p12*qx(i-1)
          beta = p21*aa*px(i-1)
          if (l /= 0) beta = beta - p12*aa*px(i-1) - p12*p21*px(i-2) + (n + s)*p12*qx(i-1)
          det = n*(n + 2.d0*s)*aa
          px(i) = (alfa*(n + s + 1.d0)*aa - aa*beta)/det
          qx(i) = (beta*(n + s - 1.d0) - bb*alfa)/det
        enddo ! i
        
      endif

      g(1) = 0.d0
      f(1) = 0.d0
      rpb = b
      do k = 2, 4
        rpb = rpb*ea
        r = rpb - b
        dr = a*rpb
        phi = (e + zz/r - v(k))*dr/cvlight
        u = dr*cvlight + phi*phiwgt
        x = -dr/r
        y = -fllp1*x*x/u + phi
        pp = px(10)
        qq = qx(10)
        do i1 = 3, 10
          i = 12 - i1
          pp = pp*r + px(i)
          qq = qq*r + qx(i)
        enddo ! i1
        g(k) = (r**s)*pp
        f(k) = (r**s)*qq
        sg   = dsign(1.d0,g(k))
        sgm1 = dsign(1.d0,g(k-1))
        if (sg*sgm1 < 0.d0) nne = nne + 1
        d(1,k-1) = u*f(k) - x*g(k)
        d(2,k-1) = x*f(k) - y*g(k)
      enddo ! k
      dg1 = d(1,1)
      dg2 = d(1,2)
      dg3 = d(1,3)
      df1 = d(2,1)
      df2 = d(2,2)
      df3 = d(2,3)
      do k = 5, k2
        km1 = k - 1
        rpb = rpb*ea
        r = rpb - b
        dr = a*rpb
        phi = (e + zz/r - v(k))*dr/cvlight
        u = dr*cvlight + phi*phiwgt
        x = -dr/r
        y = -fllp1*x*x/u + phi
        det = r83sq - x*x + u*y
        b1 = g(km1)*h83 + r1*dg1 + r2*dg2 + r3*dg3
        b2 = f(km1)*h83 + r1*df1 + r2*df2 + r3*df3
        g(k) = (b1*(h83 - x) + b2*u)/det
        f(k) = (b2*(h83 + x) - b1*y)/det
        sg   = dsign(1.d0,g(k))
        sgm1 = dsign(1.d0,g(km1))
        if (sg*sgm1 < 0.d0) nne = nne + 1
        dg1 = dg2
        dg2 = dg3
        dg3 = u*f(k) - x*g(k)
        df1 = df2
        df2 = df3
        df3 = x*f(k) - y*g(k)
      enddo ! k
      dg = dg3
    endfunction intout

      
    function hankel_core(l, arg) result(h)
!  this function uses the explicit formulas for the hankel
!  functions. for higher l-values these formulas may lead to
!  loss of significant figures. this function should be used
!  only for core states.
      integer, intent(in) :: l
      double complex, intent(in) :: arg
      double complex :: h(l) ! result array

!     .. locals ..
      double complex :: a1, a2, a3, a4
      double complex, parameter :: ci=(0.d0, 1.d0)

      if (l < 1) return ! nothing to do
      h(1) = -exp(arg*ci)/arg
      if (l == 1) return ! done
      a1 = 1.d0 - arg*ci
      h(2) = h(1)*a1/arg
      if (l == 2) return ! done
      a1 = 3.d0*a1
      a2 = arg*arg
      h(3) = h(1)*(a1 - a2)/a2
      if (l == 3) return ! done
      a1 = 5.d0*a1
      a3 = a2*arg*ci
      a4 = a2*arg
      a2 = 6.d0*a2
      h(4) = h(1)*(a1 - a2 + a3)/a4
      if (l == 4) return ! done
      a1 = 7.d0*a1
      a2 = 7.5d0*a2
      a3 = 10.d0*a3
      a4 = a4*arg
      h(5) = h(1)*(a1 - a2 + a3 + a4)/a4
      if (l == 5) return ! done
      h(6) = 9.d0*h(5)/arg - h(4)
      
      if (l /= 6) then
        write(6, fmt="(2x,' hankel :  l=',i2,' is too large')") l
        stop 'hankel.F90'
      endif ! 6

    endfunction hankel_core

endmodule AtomicCore_mod
