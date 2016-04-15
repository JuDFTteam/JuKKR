      subroutine corel(nsra, ipr, atom_id, rhoc, v, ecore, lcore, ncore, drdi, z, qc, a, b, is, nspin, nr, rmax, irmd, ebot)
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
      implicit none
      external :: intcor, simp3

      double precision, intent(in) :: a,b,rmax,z,ebot,drdi(*),v(*)
      integer, intent(in) :: atom_id,ipr,irmd,is,ncore,nr,nspin,nsra,lcore(*)
      double precision, intent(inout) :: ecore(*)
      double precision, intent(out) :: rhoc(*), qc ! todo: qc could be the return value

!     .. locals ..
      integer, parameter :: nitmax=40, irnumx=10
      double precision, parameter :: zero=0.d0
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

      rhoc(1:irmd) = zero
      rho(1:irmd) = zero

      do l = 0, 3
        e1 = (-5.d0 - ((z + 1.d0)/(l + 1.d0))**2)*1.5d0 - 50.d0 ! guess value formula
        do in = l + 1, kfg(l)
          nn = in - l - 1
          nc = nc + 1
          inuc = inuc + irnumx
          e = ecore(nc)
          ei = e
          if (ipr /= 0) write(6,fmt=F9) in,text(l),nn,spn(is),atom_id,e

          call intcor(e1,e2,rho,g,f,v,value,slope,l,nn,e,sm,nre,vlnc,a,b,z,rmax,nr,tol,irmd,ipr,nitmax,nsra)

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
            rho(ir) = zero
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

          call intcor(e1,e2,rho,g,f,v,value,slope,l,nn,e,sm,nre,vlnc,a,b,z,rmax,nr,tol,irmd,ipr,nitmax,nsra)

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

      endsubroutine ! core_electrons
