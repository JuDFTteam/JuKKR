  subroutine build_khartree(intrasite,intersite)
! Hartree kernel (Rydberd units)
! Using density basis
  use global

  implicit none

  logical,           intent(in)  :: intrasite, intersite
! -----------------------------------------------------------------
  real(kind=r8b),    parameter :: tol = 1.d-4
!   4pi
  real(kind=r8b),    parameter :: fourpi = 1.d0!/16.d0*atan(1.d0)
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ia2, ja, ja2, ndeni0, ndeni1, ndenj0, ndenj1, nri, nrj
  integer(kind=i4b) :: iq, ib, ilm, is, jq, jb, jlm, js, i, i4(4), il, jl, ir, klm
  real(kind=r8b)    :: dri(nrmax), drj(nrmax), rylm(lmmax4), rij, uij(3), prefac, rylmfac
  complex(kind=c8b) :: work1(nrmax), work2(nrmax)
  real(kind=r8b)    :: g0inter(lmmax2,lmmax2), qlmi, qlmj, start, finish
  real(kind=r8b)    :: khaintra(lmmax2,nasusc2), khainter(lmmax2,lmmax2,nasusc2,nasusc2)
  real(kind=r8b)    :: integrand(nrmax,nbmax2,nbmax2,nasusc2)
  complex(kind=c8b), external :: radint

!  write(iodb,'("in build_khartree")')
  call cpu_time(start)
!  write(iodb,'("starting ja2 loop")')
!  write(iodb,*) nasusc2, iasusc2
! Structure constants and basis multipoles
  khaintra = 0.d0; khainter = 0.d0
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
!   -------------------------------
!   no Hartree kernel for this atom
    if (ikha(ja) == 0) cycle
!   -------------------------------
!    write(iodb,*) ja, ja2
    nrj = nrpts(ja)
    drj(1:nrj) = drmesh(1:nrj,ja)!/rmesh(1:nrj,ja)**2
    ndenj0 = sum(nalmsbden(1:ja2-1))
    ndenj1 = sum(nalmsbden(1:ja2))
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
!     -------------------------------
!     no Hartree kernel for this atom
      if (ikha(ia) == 0) cycle
!     -------------------------------
      write(iodb,*) ia, ia2
      nri = nrpts(ia)
      dri(1:nri) = drmesh(1:nri,ia)!/rmesh(1:nrj,ja)**2
      ndeni0 = sum(nalmsbden(1:ia2-1))
      ndeni1 = sum(nalmsbden(1:ia2))
!     ***************
!     structural part
!     ***************
      if (ia /= ja) then
!       ++++++++ test: double factorials +++++
!        if (ia2 == 1 .and. ja2 == 2) then
!          write(iodb,'("(2m+2n-1)!!/(2m+1)!!(2n+1)!!")')
!          do il=0,nlmax2
!            do jl=0,nlmax2
!              prefac = 1.d0
!!             1/(2*il + 1)!!
!              do i=1,2*il+1,2
!                prefac = prefac/i
!              end do
!!             1/(2*jl + 1)!!
!              do i=1,2*jl+1,2
!                prefac = prefac/i
!              end do
!!             (2*(il+jl) - 1)!!
!              do i=1,2*(il+jl)-1,2
!                prefac = prefac*i
!              end do
!!             (-1)^jl
!              prefac = (-1)**jl*prefac
!              g0inter(il+1,jl+1) = prefac
!            end do
!            write(iodb,'(100f12.6)') g0inter(il+1,1:nlmax2+1)
!          end do
!          g0inter = 0.d0
!        end if
!       ++++++++++++++ end test ++++++++++++++
        write(iodb,'("structural part")')
!       connecting vector
        uij = ri(:,ja) - ri(:,ia)
        rij = sqrt(dot_product(uij,uij))
        if (abs(rij) < tol) stop 'build_khartree: rij = 0'
        uij = uij/rij
        call rymy(uij,nlmax4,lmmax4,rylm)
        do jlm=1,lmmax2
          do ilm=1,lmmax2
!           ------------------------------------------------------------
!           prefactor
            il = i2lm(2,ilm); jl = i2lm(2,jlm)
            prefac = 1.d0
!           1/(2*il + 1)!!
            do i=1,2*il+1,2
              prefac = prefac/i
            end do
!           1/(2*jl + 1)!!
            do i=1,2*jl+1,2
              prefac = prefac/i
            end do
!           (2*(il+jl) - 1)!!
            do i=1,2*(il+jl)-1,2
              prefac = prefac*i
            end do
!           8pi * (-1)^jl
            prefac = 2.d0*fourpi*(-1)**jl*prefac
!           ------------------------------------------------------------
!           structure dependent part
            rylmfac = 0.d0
            do i=1,lmmax4
!             the magnetic quantum number selection rule is built in the Gaunt coefficient
              if (i2lm(2,i) == il + jl) rylmfac = rylmfac + rgaunt(ilm,jlm,i)*rylm(i)/rij**(il+jl+1)
            end do
            g0inter(ilm,jlm) = prefac*rylmfac
!           ------------------------------------------------------------
          end do
        end do
!        write(iodb,'("g0inter for ia,ja=",2i4)') ia, ja
!        do ilm=1,lmmax2
!          write(iodb,'(1000es16.8)') g0inter(ilm,1:lmmax2)
!        end do
!     ******
      end if
!     ******
!     ------------------------------------------------------------------
      do jq=1+ndenj0,ndenj1
        i4 = i2almsbden(:,jq)
        jb = i4(1); jlm = i4(2); js = i4(3)!; ja = i4(4)
        jl = i2lm(2,jlm)
        do iq=1+ndeni0,ndeni1
          i4 = i2almsbden(:,iq)
          ib = i4(1); ilm = i4(2); is = i4(3)!; ia = i4(4)
          il = i2lm(2,ilm)
!         --------------------------------------------------------------
!         spin diagonal => WRONG!
!          if ((js == 3 .and. is == 3) .or. (js == 4 .and. is == 4)) then
          if (js > 2 .and. is > 2) then
!          if ((js == 4 .and. is == 4)) then
!           intrasite funny integral
            if (intrasite .and. ia2 == ja2 .and. ilm == jlm) then
!             integral over r' with switching between P and Q
!             ir <-> r; work2 stores data for r'
              do ir=1,nri
!               r' <= r --> phi(r') * r'**l / r**(l+1)
                work2(1:ir) = suscbasis(1:ir,jb,ilm,js,ia2)*rmesh(1:ir,ia)**il/rmesh(ir,ia)**(il+1)
!               r' > r  --> phi(r') * r**l / r'**(l+1)
                work2(ir+1:nri) = suscbasis(ir+1:nri,jb,ilm,js,ia2)*rmesh(ir,ia)**il/rmesh(ir+1:nri,ia)**(il+1)
!                if (ilm == 1 .and. is == 3 .and. ir == nri - 50) integrand(1:nri,ib,jb,ia2) = work2(1:nri)
!               integrate over r'
                work1(ir) = radint(nri,work2(1:nri),dri(1:nri),npanat(ia),ircutat(:,ia))
              end do
!             multiply by basis function and integrate over r
              work1(1:nri) = work1(1:nri)*suscbasis(1:nri,ib,ilm,is,ia2)
              qlmi = radint(nri,work1(1:nri),dri(1:nri),npanat(ia),ircutat(:,ia))
!              if (is == 3 .and. ib == 1 .and. jb == 1) write(*,'("qlmi=",5i4,2es16.8)') ia, ilm, is, ib, jb, qlmi
              khaintra(ilm,ia2) = khaintra(ilm,ia2) + 2.d0*fourpi*qlmi/(2.d0*il+1.d0)*suscnorm(iq)*suscnorm(jq)
              kernel(iq,jq) = kernel(iq,jq) + 2.d0*fourpi*qlmi/(2.d0*il+1.d0)
!             ++++ test: model radial functions ++++
!              if (ia2 == 1 .and. ib == 1 .and. is == 3) then
!!               integral over r' with switching between P and Q
!!               ir <-> r; work2 stores data for r'
!                do ir=1,nri
!!                 r' <= r --> r'**(l+2) * r'**l / r**(l+1)
!                  work2(1:ir) = rmesh(1:ir,ia)**(2*il+2)/rmesh(ir,ia)**(il+1)
!!                 r' > r  --> r'**(l+2) * r**l / r'**(l+1)
!                  work2(ir+1:nri) = rmesh(ir,ia)**il*rmesh(ir+1:nri,ia)
!!                 integrate over r'
!                  work1(ir) = radint(nri,work2(1:nri),dri(1:nri),npanat(ia),ircutat(:,ia))
!                end do
!!               multiply by r**(l+2) and integrate over r
!                work1(1:nri) = work1(1:nri)*rmesh(1:nri,ia)**(il+2)
!                qlmi = radint(nri,work1(1:nri),dri(1:nri),npanat(ia),ircutat(:,ia))
!                write(iodb,'("model qlmi=",3i4,2es16.8)') ia, ilm, is, rmesh(nri,ia), qlmi
!              end if
!             ++++++++++++++ end test ++++++++++++++
            end if
!           multipoles of the basis functions
            if (intersite .and. ia2 /= ja2) then
              work1(1:nri) = rmesh(1:nri,ia)**il*suscbasis(1:nri,ib,ilm,is,ia2)
              qlmi = radint(nri,work1,dri(1:nri),npanat(ia),ircutat(:,ia))
              work2(1:nrj) = rmesh(1:nrj,ja)**jl*suscbasis(1:nrj,jb,jlm,js,ja2)
              qlmj = radint(nrj,work2,drj(1:nrj),npanat(ja),ircutat(:,ja))
              khainter(ilm,jlm,ia2,ja2) = khainter(ilm,jlm,ia2,ja2) + qlmi*g0inter(ilm,jlm)*qlmj*suscnorm(iq)*suscnorm(jq)
              kernel(iq,jq) = kernel(iq,jq) + qlmi*g0inter(ilm,jlm)*qlmj
!              write(iodb,'("qlmi,qlmj=",8i4,2es16.8)') ia, ja, is, js, ilm, jlm, ib, jb, qlmi, qlmj
!             ++++ test: model radial functions ++++
!              if (ia2 == 1 .and. ja2 == 2 .and. ib == 1 .and. jb == 1 .and. ilm == jlm) then
!                work1(1:nri) = rmesh(1:nri,ia)**(2*il+2)
!                qlmi = radint(nri,work1,dri(1:nri),npanat(ia),ircutat(:,ia))
!                write(iodb,'("model qlmi=",3i4,2es16.8)') ia, ilm, is, rmesh(nri,ia), qlmi
!              end if
!             ++++++++++++++ end test ++++++++++++++
            end if
          end if
!         --------------------------------------------------------------
        end do
      end do
!     ------------------------------------------------------------------
      if (ia2 == ja2) then
        write(*,'("build_khartree: intra part")')
        write(*,'("ia=",i4,1000es16.8)') ia, khaintra(1:1,ia2)
!        do ir=1,nri
!          write(888,'(1000es16.8)') rmesh(ir,ia), ((real(integrand(ir,ib,jb,ia2)),ib=1,jb),jb=1,iwsusc2(1,3,ia2))
!        end do
        do il=0,nlmax
          do jl=0,nlmax2,2
            do ir=1,nri
!             r' <= r --> phi(r') * r'**l / r**(l+1)
              work2(1:ir) = phiref(1:ir,1,il,1,ia)**2*rmesh(1:ir,ia)**jl/rmesh(ir,ia)**(jl+1)
!             r' > r  --> phi(r') * r**l / r'**(l+1)
              work2(ir+1:nri) = phiref(ir+1:nri,1,il,1,ia)**2*rmesh(ir,ia)**jl/rmesh(ir+1:nri,ia)**(jl+1)
!              if (ilm == 1 .and. is == 3 .and. ir == nri - 50) integrand(1:nri,ib,jb,ia2) = work2(1:nri)
!             integrate over r'
              work1(ir) = radint(nri,work2(1:nri),dri(1:nri),npanat(ia),ircutat(:,ia))
            end do
!           multiply by basis function and integrate over r
            work1(1:nri) = work1(1:nri)*phiref(1:nri,1,il,1,ia)**2
            qlmi = radint(nri,work1(1:nri),dri(1:nri),npanat(ia),ircutat(:,ia))
!            write(*,'("F_k for ia,il,kl",3i4,es16.8)') ia, il, jl, qlmi
!           now check the double Gaunts
            do klm=1,lmmax2
              do jlm=1,lmmax
                do ilm=1,lmmax
                  if (i2lm(2,ilm) == il .and. i2lm(2,jlm) == il .and. i2lm(2,klm) == jl) then
                    qlmj = 2.d0*fourpi*rgaunt(ilm,jlm,klm)**2/(2.d0*jl + 1.d0)
!                    if (abs(qlmj) > 1.d-6) write(*,'("double gaunt for ia,klm,jlm,ilm",4i4,2es16.8)') ia, klm, jlm, ilm, qlmj, qlmi*qlmj
                  end if
                end do
              end do
            end do
          end do
        end do
      else
        write(*,'("build_khartree: inter part")')
        do ilm=1,1!lmmax2
          write(*,'("ia,ja,ilm=",3i4,1000es16.8)') ia, ja, ilm, khainter(ilm,1:1,ia2,ja2)
        end do
      end if
    end do
  end do
  call cpu_time(finish)
  write(*,'("time in build_khartree: ",f10.3," s")') finish - start
! All done!
  end subroutine build_khartree
