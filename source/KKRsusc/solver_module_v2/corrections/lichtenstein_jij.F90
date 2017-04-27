  subroutine lichtenstein_jij()
! Magnetic exchange interactions
! Eq. 19 in JMMM 67, 65-74 (1987), in Juelich convention
! One factor of 1/2 is taken into the definition of the magnetic Hamiltonian:
! H = -1/2 \sum_{ij} e_i J_{ij} e_j

  use global

  implicit none

  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0)
  real(kind=r8b),    parameter :: pi = 4.d0*atan(1.d0)
! Which atoms for Jij's
  integer(kind=i4b) :: najij
  integer(kind=i4b) :: iajij(nasusc)
! Exchange interactions
  real(kind=r8b),    allocatable :: jijupdn(:,:), jijdnup(:,:)
! Storage for angular momentum matrices
  complex(kind=c8b), allocatable :: dti(:,:), dtj(:,:), work1(:,:), work2(:,:)
  complex(kind=c8b), allocatable :: gijup(:,:), gijdn(:,:), gjiup(:,:), gjidn(:,:)
! ----------------------------------------------------------------------
  integer(kind=i4b) :: ie, i2(2), ia2, ia, ja2, ja, jalms, ialms, jlms, jlm, js, ilms, ilm, is
  complex(kind=c8b) :: trace, element


! Use ijij to find which atoms to compute Jij for
  najij = 0
  do ia=1,nasusc
    if (ijij(ia) /= 0) then
      najij = najij + 1
      iajij(najij) = ia
    end if
  end do
  write(*,'("lichtenstein_jij: RAM=",f16.3," MB")') 16.d0*8*lmmax*lmmax/1024.d0**2
  allocate(dti(lmmax,lmmax),dtj(lmmax,lmmax),work1(lmmax,lmmax),work2(lmmax,lmmax),jijupdn(najij,najij),jijdnup(najij,najij))
  allocate(gijup(lmmax,lmmax),gjiup(lmmax,lmmax),gijdn(lmmax,lmmax),gjidn(lmmax,lmmax))
  jijupdn = 0.d0; jijdnup = 0.d0
! *************
! Energy loop
  do ie=1,nescf
! *************
! ++++++++++++++
  do ja2=1,najij
! ++++++++++++++
    ja = iajij(ja2)
!   ---------------------------
!   tup - tdn
    dtj = czero
    do jlms=1,nlms
      i2 = i2lms(:,jlms)
      jlm = i2(1); js = i2(2)
      do ilms=1,nlms
        i2 = i2lms(:,ilms)
        ilm = i2(1); is = i2(2)
        if (is == js) dtj(ilm,jlm) = dtj(ilm,jlm) + (2*is-3)*tmatrix(ilms,jlms,ja,ie)
      end do
    end do
!    trace = czero
!    do ilm=1,lmmax
!      trace = trace + dtj(ilm,ilm)
!    end do
!    write(iodb,'(" tr dtj=",2i4,2es16.8)') ie, ja, trace
!   ++++++++++++++
    do ia2=1,najij
!   ++++++++++++++
      ia = iajij(ia2)
!     ---------------------------
!     tup - tdn
      dti = czero
      do jlms=1,nlms
        i2 = i2lms(:,jlms)
        jlm = i2(1); js = i2(2)
        do ilms=1,nlms
          i2 = i2lms(:,ilms)
          ilm = i2(1); is = i2(2)
          if (is == js) dti(ilm,jlm) = dti(ilm,jlm) + (2*is-3)*tmatrix(ilms,jlms,ia,ie)
        end do
      end do
!      trace = czero
!      do ilm=1,lmmax
!        trace = trace + dti(ilm,ilm)
!      end do
!      write(iodb,'(" tr dti=",2i4,2es16.8)') ie, ia, trace
!     ---------------------------
      gijup = czero
      gijdn = czero
      gjiup = czero
      gjidn = czero
      do jalms=1,nalms
        jlms = i2alms(1,jalms)
        i2 = i2lms(:,jlms)
        jlm = i2(1); js = i2(2)
        do ialms=1,nalms
          ilms = i2alms(1,ialms)
          i2 = i2lms(:,ilms)
          ilm = i2(1); is = i2(2)
          element = gstruct(ialms,jalms,ie)
!         Gij up and dn
          if (i2alms(2,ialms) == ia .and. i2alms(2,jalms) == ja) then
            if (is == js .and. is == 1) gijdn(ilm,jlm) = element
            if (is == js .and. is == 2) gijup(ilm,jlm) = element
          end if
!         Gji up and dn
          if (i2alms(2,ialms) == ja .and. i2alms(2,jalms) == ia) then
            if (is == js .and. is == 1) gjidn(ilm,jlm) = element
            if (is == js .and. is == 2) gjiup(ilm,jlm) = element
          end if
        end do
      end do
!     ---------------------------
!      trace = czero
!      do ilm=1,lmmax
!        trace = trace + gijdn(ilm,ilm)
!      end do
!      write(iodb,'(" tr gijdn=",3i4,2es16.8)') ie, ia, ja, trace
!      trace = czero
!      do ilm=1,lmmax
!        trace = trace + gijup(ilm,ilm)
!      end do
!      write(iodb,'(" tr gijup=",3i4,2es16.8)') ie, ia, ja, trace
!      trace = czero
!      do ilm=1,lmmax
!        trace = trace + gjidn(ilm,ilm)
!      end do
!      write(iodb,'(" tr gjidn=",3i4,2es16.8)') ie, ja, ia, trace
!      trace = czero
!      do ilm=1,lmmax
!        trace = trace + gjiup(ilm,ilm)
!      end do
!      write(iodb,'(" tr gjiup=",3i4,2es16.8)') ie, ja, ia, trace
!     ---------------------------
!     dti * gij --> work1
      call zgemm('N','N',lmmax,lmmax,lmmax,cone,dti,lmmax,gijup,lmmax,czero,work1,lmmax)
!     (dti * gij) * dtj --> work2
      call zgemm('N','N',lmmax,lmmax,lmmax,cone,work1,lmmax,dtj,lmmax,czero,work2,lmmax)
!     (dti * gij * dtj) * gji --> work1
      call zgemm('N','N',lmmax,lmmax,lmmax,cone,work2,lmmax,gjidn,lmmax,czero,work1,lmmax)
      if (ljionsite .and. ia2 == ja2) then
!       dti * gii --> work1
        call zgemm('N','N',lmmax,lmmax,lmmax,-cone,dti,lmmax,gjiup,lmmax,cone,work1,lmmax)
        call zgemm('N','N',lmmax,lmmax,lmmax,+cone,dti,lmmax,gjidn,lmmax,cone,work1,lmmax)
      end if
      trace = czero
      do ilm=1,lmmax
        trace = trace + work1(ilm,ilm)
      end do
!      jijupdn(ia2,ja2) = jijupdn(ia2,ja2) + aimag(descf(ie)*trace)/(2.d0*pi)
      jijupdn(ia2,ja2) = jijupdn(ia2,ja2) + aimag(descf(ie)*trace)/pi
!     ---------------------------
!     dti * gij --> work1
      call zgemm('N','N',lmmax,lmmax,lmmax,cone,dti,lmmax,gijdn,lmmax,czero,work1,lmmax)
!     (dti * gij) * dtj --> work2
      call zgemm('N','N',lmmax,lmmax,lmmax,cone,work1,lmmax,dtj,lmmax,czero,work2,lmmax)
!     (dti * gij * dtj) * gji --> work1
      call zgemm('N','N',lmmax,lmmax,lmmax,cone,work2,lmmax,gjiup,lmmax,czero,work1,lmmax)
      if (ljionsite .and. ia2 == ja2) then
!       dti * gii --> work1
        call zgemm('N','N',lmmax,lmmax,lmmax,-cone,dti,lmmax,gjiup,lmmax,cone,work1,lmmax)
        call zgemm('N','N',lmmax,lmmax,lmmax,+cone,dti,lmmax,gjidn,lmmax,cone,work1,lmmax)
      end if
      trace = czero
      do ilm=1,lmmax
        trace = trace + work1(ilm,ilm)
      end do
!      jijdnup(ia2,ja2) = jijdnup(ia2,ja2) + aimag(descf(ie)*trace)/(2.d0*pi)
      jijdnup(ia2,ja2) = jijdnup(ia2,ja2) + aimag(descf(ie)*trace)/pi
!   ++++++
    end do
!   ++++++
! ++++++
  end do
! ++++++
! ******
  end do
! ******
  write(*,'(" Magnetic exchange interactions: Jupdn(ia,ja)")')
  write(*,'(" iajij=",100i4)') iajij(1:najij)
  do ia2=1,najij
    write(*,'(100es16.8)') jijupdn(ia2,1:najij)
  end do
  write(*,'(" Magnetic exchange interactions: Jdnup(ia,ja)")')
  write(*,'(" iajij=",100i4)') iajij(1:najij)
  do ia2=1,najij
    write(*,'(100es16.8)') jijdnup(ia2,1:najij)
  end do
  deallocate(jijupdn,jijdnup,dti,dtj,gijup,gijdn,gjiup,gjidn)
! All done!
  end subroutine lichtenstein_jij
