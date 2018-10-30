  subroutine dyn_susc_real2(omega,suscylm,suscy00,suscden,analytic,nonanalytic,enhanced)
! Dynamical KS susceptibility with real frequency
! Using density basis

  implicit none

  real(kind=c8b),    intent(in)  :: omega
  complex(kind=c8b), intent(out) :: suscylm(4,4,lmmax0,lmmax0,nasusc,nasusc)
  complex(kind=c8b), intent(out) :: suscy00(4,4,lmmax,lmmax,nasusc,nasusc)
  complex(kind=c8b), intent(out) :: suscden(4,lmmax,nasusc)
  logical,           intent(in)  :: analytic, nonanalytic, enhanced
! -----------------------------------------------------------------
!   i 2pi
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0)
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cminus = (-1.d0,0.d0)
! -----------------------------------------------------------------
  complex(kind=c8b) :: gfijw(nlmsb,nlmsb), gfjiw(nlmsb,nlmsb)
  complex(kind=c8b) :: gfij0(nlmsb,nlmsb), gfji0(nlmsb,nlmsb)
  complex(kind=c8b) :: norm, de, e
  integer(kind=i4b) :: q1, lm1, l1, m1, s1, p1
  integer(kind=i4b) :: q2, lm2, l2, m2, s2, p2
  integer(kind=i4b) :: q3, lm3, l3, m3, s3, p3
  integer(kind=i4b) :: q4, lm4, l4, m4, s4, p4
  integer(kind=i4b) :: ib, jb, is, js
  integer(kind=i4b) :: i2(2), i3(3), i4(4), ie, ia, ja, jlm, ilm, j, i, iq, jq, iq1, jq1
  integer(kind=i4b) :: ipiv(ngfsum), info
  real(kind=r8b)    :: gaunti, gauntj, doublegaunt, maxelem, efermi
  integer(kind=i4b) :: ne
!  real(kind=r8b)    :: r(nalmsb), c(nalmsb), rwork(2*nalmsb), rcond, ferr(nalmsb), berr(nalmsb)
!  complex(kind=c8b) :: work(2*nalmsb), temp(nalmsb,nalmsb)
!  complex(kind=c8b) :: x(nalmsb,nalmsb)
!  character*1       :: equed

  suscylm = 0.d0; suscy00 = 0.d0; maxelem = 0.d0; suscden = 0.d0
  kssusc0 = 0.d0
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  do ja=1,nasusc    ! atom j
  do ia=1,nasusc    ! atom i
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   ----------------------------------------------------------------
    if (lanalytic) then
!   Integration on usual box contour
    do ie=1,nescf     ! energy
!    write(*,*) "ie", ie
!   ----------------------------------------------------------------
      e  = escf(ie)
      de = descf(ie)
!     get the projected GFs
      call ratval_gf(ia,ja,e,gfij0)
!      gfij0 = 0.5d0*(gfij0 + transpose(gfij0))
      gfij0 = conjg(transpose(gfij0))
      call ratval_gf(ia,ja,e,gfji0)
!      gfji0 = 0.5d0*(gfji0 + transpose(gfji0))
      call ratval_gf(ia,ja,e+omega,gfijw)
!      gfijw = 0.5d0*(gfijw + transpose(gfijw))
      call ratval_gf(ia,ja,e-omega,gfjiw)
!      gfjiw = 0.5d0*(gfjiw + transpose(gfjiw))
      gfjiw = conjg(transpose(gfjiw))
      do jq=1+sum(nalmsbgf(1:ja-1)),sum(nalmsbgf(1:ja))
        i3 = i2almsbgf(:,jq)
        q3 = i3(1); q4 = i3(2)!; ja = i3(3)
        do iq=1+sum(nalmsbgf(1:ia-1)),sum(nalmsbgf(1:ia))
          i3 = i2almsbgf(:,iq)
          q1 = i3(1); q2 = i3(2)!; ia = i3(3)
          norm = (conjg(de)*gfij0(q2,q3)*gfjiw(q4,q1) - de*gfijw(q2,q3)*gfji0(q4,q1))/i2pi
          kssusc0(iq,jq) = kssusc0(iq,jq) + norm
        end do
      end do
!   ----------------------------------------------------------------
    end do            ! energy
    end if
!   ----------------------------------------------------------------
!   ----------------------------------------------------------------
    if (lnonanalytic) then
!   Integration along real axis
    efermi = real(escf(nescf))
    ne = max(int(abs(omega)/domega),11)
    do ie=1,ne        ! energy
!    write(*,*) "ie", ie
!     --------------------------------------------------------------
      e  = efermi + omega*((ie - 1)/(ne - 1.d0) - 1.d0)
      de = omega/(ne - 1) 
      if (ie == 1 .or. ie == ne) de = 0.5d0*de
!     get the projected GFs
      call ratval_gf(ia,ja,e,gfji0)
      gfji0 = conjg(transpose(gfji0))
!      call ratval_gf(ja,ia,e,gfji0)
      call ratval_gf(ia,ja,e+omega,gfijw)
!      call ratval_gf(ja,ia,e+omega,gfjiw)
!      gfjiw = conjg(transpose(gfjiw))
      do jq=1+sum(nalmsbgf(1:ja-1)),sum(nalmsbgf(1:ja))
        i3 = i2almsbgf(:,jq)
        q3 = i3(1); q4 = i3(2)!; ja = i3(3)
        do iq=1+sum(nalmsbgf(1:ia-1)),sum(nalmsbgf(1:ia))
          i3 = i2almsbgf(:,iq)
          q1 = i3(1); q2 = i3(2)!; ia = i3(3)
          norm = de*gfijw(q2,q3)*gfji0(q4,q1)/i2pi
          kssusc0(iq,jq) = kssusc0(iq,jq) + norm 
        end do
      end do
!   ----------------------------------------------------------------
    end do            ! energy
    end if
!   ----------------------------------------------------------------
!    maxelem = maxval(abs(kssusc0))
!    where (abs(kssusc0) < susctol*maxelem) kssusc0 = 0.d0
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  end do            ! atom i
  end do            ! atom j
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  call basis_change
  where (abs(kssusc0) < susctol) kssusc0 = 0.d0
! ------------------------------------------------------------------
  if (lenhanced) then
    call zgemm('N','N',ndensum,ndensum,ndensum,cminus,kssusc0,ngfsum,kernel,ngfsum,czero,denominator,ngfsum)
    do iq=1,ndensum
      denominator(iq,iq) = denominator(iq,iq) + 1.d0
    end do
    call zgesv(ndensum,ndensum,denominator,ngfsum,ipiv,kssusc0,ngfsum,info)
    if (info /= 0) stop 'dyn_susc_real2: failure in zgesv'
!    call zgesvx('N','N',nalmsb,nalmsb,denominator,nalmsb,temp,nalmsb,ipiv,equed,r,c,kssusc0,nalmsb,x,nalmsb,rcond,ferr,berr,work,rwork,info)
!    write(iodb,'("condition number of enhancement factor=",es16.3)') rcond
!    write(iodb,'("fwd & bkwd error in solution=",2es16.3)') maxval(abs(ferr)), maxval(abs(berr))
!    if (info /= 0) stop 'dyn_susc_real2: failure in zgesvx'
!    kssusc0 = x
  end if
! ------------------------------------------------------------------
! Symmetrize
!  call symmetrize(nalmsb,kssusc0,susctol)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! multipoles
  suscylm = 0.d0
  do jq=1,ndensum
    i4 = i2almsbden(:,jq)
    jb = i4(1); jlm = i4(2); js = i4(3); ja = i4(4)
    do iq=1,ndensum
      i4 = i2almsbden(:,iq)
      ib = i4(1); ilm = i4(2); is = i4(3); ia = i4(4)
      suscylm(is,js,ilm,jlm,ia,ja) = suscylm(is,js,ilm,jlm,ia,ja) + suscnorm(iq)*kssusc0(iq,jq)*suscnorm(jq)
    end do
  end do
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Final filtering (global)
!  maxelem = maxval(abs(suscylm))
!  where (abs(suscylm) < susctol*maxelem) suscylm = 0.d0
!  maxelem = maxval(abs(suscy00))
!  where (abs(suscy00) < susctol*maxelem) suscy00 = 0.d0
! All done!
  end subroutine dyn_susc_real2

