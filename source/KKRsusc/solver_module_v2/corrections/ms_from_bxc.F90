  subroutine ms_from_bxc
! check +- and -+ contributions to the z component of the spin moment
! ms = chi0*bxc
  use global

  implicit none

  complex(kind=c8b) :: mspm(lmmax,lmmax,nasusc2,nasusc2)
  complex(kind=c8b) :: msmp(lmmax,lmmax,nasusc2,nasusc2)
  complex(kind=c8b) :: msgf(lmmax,lmmax,nasusc2)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: twopi = 8.d0*atan(1.d0), tol = 1.d-7
  complex(kind=c8b), parameter :: i2pi = (0.d0,twopi), iu = (0.d0,1.d0)
  complex(kind=c8b) :: gfij(nlmsb,nlmsb), gfji(nlmsb,nlmsb)
  complex(kind=c8b) :: vdiff(nbmax,nbmax,lmmax,lmmax,nasusc2), msum(lmmax,lmmax,2), tmp, msph(nasusc2,nasusc2)
  complex(kind=c8b) :: msigmaz(nsmax,nsmax,nasusc2), msigmap(nsmax,nsmax,nasusc2), msigmam(nsmax,nsmax,nasusc2), spinrot(nsmax,nsmax)
  integer(kind=i4b) :: i3(3), ie, ia, ja, i1, i2, j1, j2, ib1, ilm1, is1, ib2, ilm2, is2, jb1, jlm1, js1, jb2, jlm2, js2, ia2, ja2, k
  real(kind=r8b)    :: re, im, u0(3), u1(3), orbrot(lmmax,lmmax), rot120p(lmmax,lmmax), rot120m(lmmax,lmmax)


! Rotation test
!  u0 = (/0.d0,0.d0,1.d0/)
!  u1 = (/0.d0,0.d0,1.d0/)
!  call orb_rotation(u0,u1,orbrot)
!  write(*,'("Rot: identity")')
!  do ilm1=1,lmmax
!    write(*,'(1000f8.4)') orbrot(ilm1,1:lmmax)
!  end do
!  u0 = (/0.d0,0.d0,1.d0/)
!  u1 = (/0.d0,0.d0,-1.d0/)
!  call orb_rotation(u0,u1,orbrot)
!  write(*,'("Rot: mirror z")')
!  do ilm1=1,lmmax
!    write(*,'(1000f8.4)') orbrot(ilm1,1:lmmax)
!  end do
!  u0 = (/1.d0,1.d0,0.d0/)
!  u1 = (/-1.d0,-1.d0,0.d0/)
!  call orb_rotation(u0,u1,orbrot)
!  write(*,'("Rot: mirror xy")')
!  do ilm1=1,lmmax
!    write(*,'(1000f8.4)') orbrot(ilm1,1:lmmax)
!  end do
!  u0 = (/0.d0,-1.d0,0.d0/)
!  u1 = (/sqrt(3.d0)/2.d0,0.5d0,0.d0/)
!  call orb_rotation(u0,u1,rot120p)
!  write(*,'("Rot: +2pi/3 about z")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') rot120p(ilm1,1:lmmax)
!  end do
!  u0 = (/0.d0,-1.d0,0.d0/)
!  u1 = (/-sqrt(3.d0)/2.d0,0.5d0,0.d0/)
!  call orb_rotation(u0,u1,rot120m)
!  write(*,'("Rot: -2pi/3 about z")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') rot120m(ilm1,1:lmmax)
!  end do
! GF blocks
!  write(*,'("Test diagonal blocks")')
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,1)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,1)
!      gfij(ilm1,ilm2) = gstruct(i1,i2,nescf)
!      re = real(gfij(ilm1,ilm2))
!      im = aimag(gfij(ilm1,ilm2))
!      if (abs(re) < tol) re = 0.d0
!      if (abs(im) < tol) im = 0.d0
!      gfij(ilm1,ilm2) = cmplx(re,im)
!    end do
!  end do
!  write(*,'(" dn for Fe Fe")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1-1,1:nlms:2)
!  end do
!  gfji(1:lmmax,1:lmmax) = matmul(matmul(rot120p,gfij(1:nlms:2,1:nlms:2)),rot120m)
!  write(*,'(" Diff")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1-1,1:nlms:2) - gfji(ilm1,1:lmmax)
!  end do
!  write(*,'(" up for Fe Fe")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1,2:nlms:2)
!  end do
!  gfji(1:lmmax,1:lmmax) = matmul(matmul(rot120p,gfij(2:nlms:2,2:nlms:2)),rot120m)
!  write(*,'(" Diff")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1,2:nlms:2) - gfji(ilm1,1:lmmax)
!  end do
!
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,2)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,2)
!      gfij(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,3)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,3)
!      gfji(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  write(*,'(" dn for Pt Pt 2 2")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1-1,1:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(1:nlms:2,1:nlms:2)),rot120p)
!  write(*,'(" Diff 3 3 - 2 2")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1-1,1:nlms:2) - work(ilm1,1:lmmax)
!  end do
!  write(*,'(" up for Pt Pt 2 2")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1,2:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(2:nlms:2,2:nlms:2)),rot120p)
!  write(*,'(" Diff 3 3 - 2 2")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1,2:nlms:2) - work(ilm1,1:lmmax)
!  end do
!
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,3)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,3)
!      gfij(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,4)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,4)
!      gfji(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  write(*,'(" dn for Pt Pt 3 3")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1-1,1:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(1:nlms:2,1:nlms:2)),rot120p)
!  write(*,'(" Diff 4 4 - 3 3")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1-1,1:nlms:2) - work(ilm1,1:lmmax)
!  end do
!  write(*,'(" up for Pt Pt 3 3")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1,2:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(2:nlms:2,2:nlms:2)),rot120p)
!  write(*,'(" Diff 4 4 - 3 3")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1,2:nlms:2) - work(ilm1,1:lmmax)
!  end do
!
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,4)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,4)
!      gfij(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,2)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,2)
!      gfji(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  write(*,'(" dn for Pt Pt 4 4")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1-1,1:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(1:nlms:2,1:nlms:2)),rot120p)
!  write(*,'(" Diff 2 2 - 4 4")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1-1,1:nlms:2) - work(ilm1,1:lmmax)
!  end do
!  write(*,'(" up for Pt Pt 4 4")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1,2:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(2:nlms:2,2:nlms:2)),rot120p)
!  write(*,'(" Diff 2 2 - 4 4")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1,2:nlms:2) - work(ilm1,1:lmmax)
!  end do
!
!  write(*,'("Test off-diagonal blocks")')
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,2)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,1)
!      gfij(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,3)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,1)
!      gfji(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  write(*,'(" dn for Fe Pt 1 2")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1-1,1:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(1:nlms:2,1:nlms:2)),rot120p)
!  write(*,'(" Diff 1 2 - 1 3")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1-1,1:nlms:2) - work(ilm1,1:lmmax)
!  end do
!  write(*,'(" up for Fe Pt 1 2")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1,2:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(2:nlms:2,2:nlms:2)),rot120p)
!  write(*,'(" Diff 1 2 - 1 3")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1,2:nlms:2) - work(ilm1,1:lmmax)
!  end do
!
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,3)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,1)
!      gfij(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,4)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,1)
!      gfji(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  write(*,'(" dn for Fe Pt 1 3")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1-1,1:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(1:nlms:2,1:nlms:2)),rot120p)
!  write(*,'(" Diff 1 3 - 1 4")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1-1,1:nlms:2) - work(ilm1,1:lmmax)
!  end do
!  write(*,'(" up for Fe Pt 1 3")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1,2:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(2:nlms:2,2:nlms:2)),rot120p)
!  write(*,'(" Diff 1 3 - 1 4")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1,2:nlms:2) - work(ilm1,1:lmmax)
!  end do
!
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,4)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,1)
!      gfij(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  do ilm2=1,nlms
!    i2 = alms2i(ilm2,2)
!    do ilm1=1,nlms
!      i1 = alms2i(ilm1,1)
!      gfji(ilm1,ilm2) = gstruct(i1,i2,nescf)
!    end do
!  end do
!  write(*,'(" dn for Fe Pt 1 4")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1-1,1:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(1:nlms:2,1:nlms:2)),rot120p)
!  write(*,'(" Diff 1 4 - 1 1")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1-1,1:nlms:2) - work(ilm1,1:lmmax)
!  end do
!  write(*,'(" up for Fe Pt 1 4")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfij(2*ilm1,2:nlms:2)
!  end do
!  work(1:lmmax,1:lmmax) = matmul(matmul(rot120m,gfij(2:nlms:2,2:nlms:2)),rot120p)
!  write(*,'(" Diff 1 4 - 1 2")')
!  do ilm1=1,lmmax
!    write(*,'(1000es12.4)') gfji(2*ilm1,2:nlms:2) - work(ilm1,1:lmmax)
!  end do
!
!
! compute Vup - Vdn = Tr magdir.sigma Vlmsb
  vdiff = 0.d0
  do ia2=1,nasusc2
    ia = iasusc2(ia2)
!   spin rotation
    u0 = (/0.d0,0.d0,1.d0/); u1 = magdir(:,ia)
    call spin_rotation(u0,u1,pauli,spinrot)
!   spin projection on local direction of magnetization: rotated \sigma_z
    msigmaz(:,:,ia2) = matmul(matmul(spinrot,pauli(:,:,3)),conjg(transpose(spinrot)))
    write(*,'("msigmaz for ia=",i4)') ia
    write(*,'(4f8.4)') ((msigmaz(is1,is2,ia2),is2=1,nsmax),is1=1,nsmax)
!   spin flip from dn to up: rotated \sigma_x + iu \sigma_y
    msigmap(:,:,ia2) = 0.5d0*matmul(matmul(spinrot,pauli(:,:,1) + iu*pauli(:,:,2)),conjg(transpose(spinrot)))
    write(*,'("msigmap for ia=",i4)') ia
    write(*,'(4f8.4)') ((msigmap(is1,is2,ia2),is2=1,nsmax),is1=1,nsmax)
!   spin flip from up to dn: rotated \sigma_x - iu \sigma_y
    msigmam(:,:,ia2) = 0.5d0*matmul(matmul(spinrot,pauli(:,:,1) - iu*pauli(:,:,2)),conjg(transpose(spinrot)))
    write(*,'("msigmam for ia=",i4)') ia
    write(*,'(4f8.4)') ((msigmam(is1,is2,ia2),is2=1,nsmax),is1=1,nsmax)
!   extract Vdiff using this projection
    do i2=1,nlmsba(ia)
      i3 = i2lmsb(:,i2,ia)
      ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
      do i1=1,nlmsba(ia)
        i3 = i2lmsb(:,i1,ia)
        ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
        vdiff(ib1,ib2,ilm1,ilm2,ia2) = vdiff(ib1,ib2,ilm1,ilm2,ia2) + vlmsbgf(i1,i2,ia)*msigmaz(is2,is1,ia2)
      end do
    end do
    write(iodb,'("vdiff for ia=",i4)') ia
    do ilm1=1,lmmax
    do ib1=1,nbmax
      write(iodb,'(1000es10.1)') vdiff(ib1,1:nbmax,ilm1,1:lmmax,ia2)
    end do
    end do
  end do
! add up KS +- and -+ susc multiplied by Vdup - Vdn in the basis
  mspm = 0.d0; msmp = 0.d0; msgf = 0.d0
  do ja2=1,nasusc2
    ja = iasusc2(ja2)
    do ia2=1,nasusc2
      ia = iasusc2(ia2)
      do ie=1,nescf
        call projected_gf(ie,ia,ja,gfij,.true.,.true.)
        call projected_gf(ie,ja,ia,gfji,.true.,.true.)
!       ----------------------------------------------------------------
        if (ia == ja) then
!         Magnetization projected on local spin quantization axis
          do i2=1,nlmsba(ia)
            i3 = i2lmsb(:,i2,ia)
            ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
            do i1=1,nlmsba(ia)
              i3 = i2lmsb(:,i1,ia)
              ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
              tmp = (conjg(descf(ie)*msigmaz(is1,is2,ia2)*gfji(i2,i1)) - descf(ie)*msigmaz(is2,is1,ia2)*gfij(i1,i2))/i2pi 
              msgf(ilm1,ilm2,ia2) = msgf(ilm1,ilm2,ia2) + overlap(i1,i2,ia)*tmp
            end do
          end do
        end if
!       ----------------------------------------------------------------
      do j2=1,nlmsba(ja)
      do j1=1,nlmsba(ja)
        i3 = i2lmsb(:,j1,ia)
        jb1 = i3(1); jlm1 = i3(2); js1 = i3(3)
        i3 = i2lmsb(:,j2,ia)
        jb2 = i3(1); jlm2 = i3(2); js2 = i3(3)
!       ----------------------------------------------------------------
!       spin of ja (dn,up)
!        if (js1 == 2 .and. js2 == 1) then
          do i2=1,nlmsba(ia)
          do i1=1,nlmsba(ia)
            i3 = i2lmsb(:,i1,ia)
            ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
            i3 = i2lmsb(:,i2,ia)
            ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
!           spin of ia (up,dn)
!            if (is1 == 2 .and. is2 == 1) then
!              write(*,'("+- ",7i4)') ie, ia, ja, i1, i2, j1, j2
!              tmp = (conjg(descf(ie)*gfji(j1,i1)*gfij(i2,j2)) - descf(ie)*gfij(i1,j1)*gfji(j2,i2))/i2pi
              tmp = (conjg(descf(ie)*msigmam(js2,js1,ja2)*gfji(j1,i1)*msigmap(is1,is2,ia2)*gfij(i2,j2)) &
                      - descf(ie)*msigmap(is2,is1,ia2)*gfij(i1,j1)*msigmam(js1,js2,ja2)*gfji(j2,i2))/i2pi
!              mspm(ilm1,ilm2,ia,ja) = mspm(ilm1,ilm2,ia,ja) + overlap(i1,i2,ia)*tmp*vdiff(jb1,jb2,jlm1,jlm2,ja)
              mspm(ilm1,ilm2,ia2,ja2) = mspm(ilm1,ilm2,ia2,ja2) + overlap(i1,i2,ia)*tmp*vdiff(jb1,jb2,jlm1,jlm2,ja2)
!            end if
!          end do
!          end do
!       ----------------------------------------------------------------
!       spin of ja (up,dn)
!        else if (js1 == 1 .and. js2 == 2) then
!          do i2=1,nlmsba(ia)
!          do i1=1,nlmsba(ia)
!            i3 = i2lmsb(:,i1,ia)
!            ib1 = i3(1); ilm1 = i3(2); is1 = i3(3)
!            i3 = i2lmsb(:,i2,ia)
!            ib2 = i3(1); ilm2 = i3(2); is2 = i3(3)
!           spin of ia (dn,up)
!            if (is1 == 1 .and. is2 == 2) then
!              write(*,'("+- ",7i4)') ie, ia, ja, i1, i2, j1, j2
!              tmp = (conjg(descf(ie)*gfji(j1,i1)*gfij(i2,j2)) - descf(ie)*gfij(i1,j1)*gfji(j2,i2))/i2pi
              tmp = (conjg(descf(ie)*msigmap(js2,js1,ja2)*gfji(j1,i1)*msigmam(is1,is2,ia2)*gfij(i2,j2)) &
                      - descf(ie)*msigmam(is2,is1,ia2)*gfij(i1,j1)*msigmap(js1,js2,ja2)*gfji(j2,i2))/i2pi
!              if (abs(tmp) < tol) tmp = 0.d0
              msmp(ilm1,ilm2,ia2,ja2) = msmp(ilm1,ilm2,ia2,ja2) + overlap(i1,i2,ia)*tmp*vdiff(jb1,jb2,jlm1,jlm2,ja2)
!            end if
          end do
          end do
!        end if
!       ----------------------------------------------------------------
      end do
      end do
    end do
  end do
  end do
! filter zeros
!  do ja=1,nasusc
!    do ia=1,nasusc
!      do ilm2=1,lmmax
!        do ilm1=1,lmmax
!          re = real(mspm(ilm1,ilm2,ia,ja))
!          im = aimag(mspm(ilm1,ilm2,ia,ja))
!          if (abs(re) < tol) re = 0.d0
!          if (abs(im) < tol) im = 0.d0
!          mspm(ilm1,ilm2,ia,ja) = cmplx(re,im)
!          re = real(msmp(ilm1,ilm2,ia,ja))
!          im = aimag(msmp(ilm1,ilm2,ia,ja))
!          if (abs(re) < tol) re = 0.d0
!          if (abs(im) < tol) im = 0.d0
!          msmp(ilm1,ilm2,ia,ja) = cmplx(re,im)
!        end do
!      end do
!      msmp(:,:,ia,ja) = (conjg(transpose(msmp(:,:,ia,ja))) - msmp(:,:,ia,ja))/i2pi
!    end do
!  end do
! print results
!  write(*,'("ms_from_bxc")')
!  do ia=1,nasusc
!    write(*,'("ia=",i4)') ia
!    write(*,'("mspm=")')
!    do ilm1=1,lmmax
!      write(*,'(1000es12.3)') sum(mspm(ilm1,1:lmmax,ia,:),dim=2)
!    end do
!    write(*,'("msmp=")')
!    do ilm1=1,lmmax
!      write(*,'(1000es12.3)') sum(msmp(ilm1,1:lmmax,ia,:),dim=2)
!    end do
!    do ja=1,nasusc
!      write(*,'("contribution from ja=",i4)') ja
!      write(*,'("mspm=")')
!      do ilm1=1,lmmax
!        write(*,'(1000es12.3)') mspm(ilm1,1:lmmax,ia,ja)
!      end do
!      write(*,'("msmp=")')
!      do ilm1=1,lmmax
!        write(*,'(1000es12.3)') msmp(ilm1,1:lmmax,ia,ja)
!      end do
!    end do
!  end do
! multipoles
  do jlm1=1,lmmax2
    msph = 0.d0
    do ja2=1,nasusc2
      do ia2=1,nasusc2
        do ilm2=1,lmmax
          do ilm1=1,lmmax
            msph(ia2,ja2) = msph(ia2,ja2) + rgaunt(ilm1,ilm2,jlm1)*mspm(ilm1,ilm2,ia2,ja2)
          end do
        end do
        re = real(msph(ia2,ja2)); im = aimag(msph(ia2,ja2))
        if (abs(re) < tol) re = 0.d0
        if (abs(im) < tol) im = 0.d0
        msph(ia2,ja2) = cmplx(re,im)
      end do
    end do
    write(*,'("mspm: multipole lm=",2i4)') i2lm(:,jlm1)
    do ia2=1,nasusc2
      tmp = 0.d0
      do ilm2=1,lmmax
        do ilm1=1,lmmax
          tmp = tmp + rgaunt(ilm1,ilm2,jlm1)*msgf(ilm1,ilm2,ia2)
        end do
      end do
      re = real(tmp); im = aimag(tmp)
      if (abs(re) < tol) re = 0.d0
      if (abs(im) < tol) im = 0.d0
      tmp = cmplx(re,im)
      write(*,'(i4,4es16.8," | ",1000es16.8)') iasusc2(ia2), tmp, sum(msph(ia2,:)), msph(ia2,:)
    end do
    msph = 0.d0
    do ja2=1,nasusc2
      do ia2=1,nasusc2
        do ilm2=1,lmmax
          do ilm1=1,lmmax
            msph(ia2,ja2) = msph(ia2,ja2) + rgaunt(ilm1,ilm2,jlm1)*msmp(ilm1,ilm2,ia2,ja2)
          end do
        end do
        re = real(msph(ia2,ja2)); im = aimag(msph(ia2,ja2))
        if (abs(re) < tol) re = 0.d0
        if (abs(im) < tol) im = 0.d0
        msph(ia2,ja2) = cmplx(re,im)
      end do
    end do
    where (abs(msph) < 1.d-7) msph = 0.d0
    write(*,'("msmp: multipole lm=",2i4)') i2lm(:,jlm1)
    do ia2=1,nasusc2
      tmp = 0.d0
      do ilm2=1,lmmax
        do ilm1=1,lmmax
          tmp = tmp + rgaunt(ilm1,ilm2,jlm1)*msgf(ilm1,ilm2,ia2)
        end do
      end do
      re = real(tmp); im = aimag(tmp)
      if (abs(re) < tol) re = 0.d0
      if (abs(im) < tol) im = 0.d0
      tmp = cmplx(re,im)
      write(*,'(i4,4es16.8," | ",1000es16.8)') ia2, tmp, sum(msph(ia2,:)), msph(ia2,:)
    end do
  end do
! All done!
  end subroutine ms_from_bxc
