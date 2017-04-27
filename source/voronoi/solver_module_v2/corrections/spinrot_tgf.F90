  subroutine spinrot_tgf(magdir0,magdir1)
! Rotates the entire projected GF according to local spin rotations
  use global

  implicit none

! Initial spin axes
  real(kind=r8b), intent(in) :: magdir0(3,nasusc)
! Final spin axis
  real(kind=r8b), intent(in) :: magdir1(3,nasusc)
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
  real(kind=r8b),    parameter :: pi = 4.d0*atan(1.d0)
  complex(kind=c8b) :: spinrot(nsmax,nsmax,nasusc), trace, qe, ms(3), mo(3)
  complex(kind=c8b), allocatable :: pzsave(:,:), gfsave(:,:), tmsave(:,:), gsave(:,:)
  complex(kind=c8b), allocatable :: sr1(:,:,:), sr2(:,:,:), dtmat(:,:,:), spo(:,:,:)
  integer(kind=i4b) :: ie, i3(3), i2(2), ilms, jlms, ia, ja, jalms, ialms
  integer(kind=i4b) :: i, j, k, l, is, js, ib, jb, ilm, jlm
  real(kind=r8b)    :: start, finish, time1, time2


  allocate(pzsave(nlmsb,nlms),gfsave(nlmsb,nlmsb),tmsave(nlms,nlms),gsave(nlms,nlms),spo(nlms,nlms,nasusc))
  allocate(sr1(nlms,nlms,nasusc),sr2(nlmsb,nlmsb,nasusc),dtmat(nlms,nlms,nasusc))
! Spin rotation matrices
  sr1 = czero; sr2 = czero
  do ia=1,nasusc
    call spin_rotation(magdir0(:,ia),magdir1(:,ia),pauli,spinrot(:,:,ia))
!    write(*,'(i4,8f8.4)') ia, spinrot(:,:,ia)
!    write(*,'(" rot sx ",4f8.4)') matmul(spinrot(:,:,ia),matmul(pauli(:,:,1),conjg(transpose(spinrot(:,:,ia)))))
!    write(*,'(" rot sy ",4f8.4)') matmul(spinrot(:,:,ia),matmul(pauli(:,:,2),conjg(transpose(spinrot(:,:,ia)))))
!    write(*,'(" rot sz ",4f8.4)') matmul(spinrot(:,:,ia),matmul(pauli(:,:,3),conjg(transpose(spinrot(:,:,ia)))))
!    write(*,'(" rot s0 ",4f8.4)') matmul(spinrot(:,:,ia),matmul(pauli(:,:,4),conjg(transpose(spinrot(:,:,ia)))))
!   for the nlms basis
    do ilm=1,lmmax
      do js=1,nsmax
        j = lms2i(ilm,js)
        do is=1,nsmax
          i = lms2i(ilm,is)
          sr1(i,j,ia) = spinrot(is,js,ia)
        end do
      end do
    end do
!   for the nlmsb basis
    do j=1,nlmsba(ia)
      i3 = i2lmsb(:,j,ia)
      jb = i3(1); jlm = i3(2); js = i3(3)
      do i=1,nlmsba(ia)
        i3 = i2lmsb(:,i,ia)
        ib = i3(1); ilm = i3(2); is = i3(3)
        if (ilm == jlm .and. ib == jb) sr2(i,j,ia) = spinrot(is,js,ia)
      end do
    end do
  end do
! **************
  time1 = 0.d0; time2 = 0.d0
  qe = czero; ms = czero; mo = czero; spo = czero
  do ie=1,nesusc
    call cpu_time(start)
    dtmat = czero
    do ia=1,nasusc
      gfsave = czero; pzsave = czero; tmsave = czero
!     ------------------------------------------------------------------
!                   Rotation of the single-site GF
!     ------------------------------------------------------------------
!     gfpqsave = gfpq.U^\dagger
      call zgemm('N','C',nlmsba(ia),nlmsba(ia),nlmsba(ia),cone,gfpq(:,:,ia,ie),nlmsb,sr2(:,:,ia),nlmsb,czero,gfsave,nlmsb)
!     gfpq = U.gfpqsave
      call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ia),cone,sr2(:,:,ia),nlmsb,gfsave,nlmsb,czero,gfpq(:,:,ia,ie),nlmsb)
!     ------------------------------------------------------------------
!                 Rotation of the scattering solutions
!     ------------------------------------------------------------------
!     pzsave = pzr.U^\dagger
      call zgemm('N','C',nlmsba(ia),nlms,nlms,cone,pzr(:,:,ia,ie),nlmsb,sr1(:,:,ia),nlms,czero,pzsave,nlmsb)
!     pzr = U.pzsave
      call zgemm('N','N',nlmsba(ia),nlms,nlmsba(ia),cone,sr2(:,:,ia),nlmsb,pzsave,nlmsb,czero,pzr(:,:,ia,ie),nlmsb)
!     pzsave = pzl.U^T
      call zgemm('N','T',nlmsba(ia),nlms,nlms,cone,pzl(:,:,ia,ie),nlmsb,sr1(:,:,ia),nlms,czero,pzsave,nlmsb)
!     pzl = U^*.pzsave
      gfsave = conjg(sr2(:,:,ia))
      call zgemm('N','N',nlmsba(ia),nlms,nlmsba(ia),cone,gfsave,nlmsb,pzsave,nlmsb,czero,pzl(:,:,ia,ie),nlmsb)
!     ------------------------------------------------------------------
!                    Rotation of the t-matrices    
!     ------------------------------------------------------------------
!     tmsave = tmat.U^\dagger
      call zgemm('N','C',nlms,nlms,nlms,cone,tmatrix(:,:,ia,ie),nlms,sr1(:,:,ia),nlms,czero,tmsave,nlms)
!     dtmat = U.tmsave
      call zgemm('N','N',nlms,nlms,nlms,cone,sr1(:,:,ia),nlms,tmsave,nlms,czero,dtmat(:,:,ia),nlms)
!     dtmat = dtmat - tmat
      if (.not.lfreegf) call zaxpy(nlms*nlms,cminus,tmatrix(:,:,ia,ie),1,dtmat(:,:,ia),1)
!      trace = czero
!      do ilms=1,nlms
!        trace = trace + dtmat(ilms,ilms,ia)
!      end do
!      write(iodb,'("spinrot_tgf: dtmat trace for ie, ia=",2i4,2es16.8)') ie, ia, trace
!      if (ia < 3 .and. ie == nesusc) then
!        do ilms=1,nlms
!          write(iodb,'("dtmat=",100es10.2)') dtmat(ilms,1:nlms,ia)
!        end do
!      end if
    end do
    call cpu_time(finish)
    time1 = time1 + finish - start
!   --------------------------------------------------------------------
!                    Updating the structural GF    
!   --------------------------------------------------------------------
    call cpu_time(start)
!    write(iodb,'("spinrot_tgf: ie, dims=",20i8)') ie, nasusc, nlms, nalms, shape(alms2i), shape(dtmat), shape(gstruct(:,:,ie))
    call structural_gf(ie,nasusc,nlms,nalms,alms2i,dtmat,gstruct(:,:,ie))
!   Update total t-matrix: tmat = tmat + dtmat
    if (lfreegf) tmatrix(:,:,:,ie) = czero
    call zaxpy(nlms*nlms*nasusc,cone,dtmat,1,tmatrix(:,:,:,ie),1)
    call cpu_time(finish)
    time2 = time2 + finish - start
! ******
  end do
! ******
  write(*,'("spinrot_tgf: site diagonal time",f16.3," s")') time1
  write(*,'("spinrot_tgf: structural GF time",f16.3," s")') time2
! ----------------------------------------------------------------------
  deallocate(pzsave,gfsave,tmsave,sr1,sr2,dtmat,gsave,spo)
! All done!
  end subroutine spinrot_tgf
