  subroutine spinrot_tgf_sph(magdir0,magdir1)
! Rotates the entire projected GF according to local spin rotations
! Makes use of the assumption that the input GF is in collinear ASA form
  use global

  implicit none

! Initial spin axes
  real(kind=r8b), intent(in) :: magdir0(3,nasusc)
! Final spin axis
  real(kind=r8b), intent(in) :: magdir1(3,nasusc)
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0), cminus = (-1.d0,0.d0)
  complex(kind=c8b), allocatable :: dtmat(:,:,:)
  complex(kind=c8b) :: spinrot(nsmax,nsmax,nasusc)
  complex(kind=c8b) :: input(nsmax,nsmax), output(nsmax,nsmax), ul(nsmax,nsmax), ur(nsmax,nsmax)
  integer(kind=i4b) :: ie, i3(3), i2(2), ilms, jlms, ia, ja
  integer(kind=i4b) :: i, j, is, js, ib, jb, ilm, jlm
  real(kind=r8b)    :: start, finish, time1, time2


  allocate(dtmat(nlms,nlms,nasusc))
! Spin rotation matrices
  do ia=1,nasusc
    call spin_rotation(magdir0(:,ia),magdir1(:,ia),pauli,spinrot(:,:,ia))
!    write(*,'(i4,8f8.4)') ia, spinrot(:,:,ia)
!    write(*,'(" rot sx ",4f8.4)') matmul(spinrot(:,:,ia),matmul(pauli(:,:,1),conjg(transpose(spinrot(:,:,ia)))))
!    write(*,'(" rot sy ",4f8.4)') matmul(spinrot(:,:,ia),matmul(pauli(:,:,2),conjg(transpose(spinrot(:,:,ia)))))
!    write(*,'(" rot sz ",4f8.4)') matmul(spinrot(:,:,ia),matmul(pauli(:,:,3),conjg(transpose(spinrot(:,:,ia)))))
!    write(*,'(" rot s0 ",4f8.4)') matmul(spinrot(:,:,ia),matmul(pauli(:,:,4),conjg(transpose(spinrot(:,:,ia)))))
  end do
  time1 = 0.d0; time2 = 0.d0
! **************
  do ie=1,nesusc  ! energy
! **************
    call cpu_time(start)
    dtmat = czero
!   ++++++++++++++
    do ia=1,nasusc   ! atoms
!   ++++++++++++++
      ul = spinrot(:,:,ia)
      ur = conjg(transpose(spinrot(:,:,ia)))
!     ------------------------------------------------------------------
!                   Rotation of the single-site GF
!     ------------------------------------------------------------------
!     Onsite GF should be diagonal in lms, but not in basis functions
!     Do NOT use a spin-dependent basis!
      do ilm=1,lmmax
        do jb=1,iwsusc(i2lm(2,ilm),1,ia)
        do ib=1,iwsusc(i2lm(2,ilm),1,ia)
!         first copy collinear elements
          input(:,:) = czero
          do is=1,nsmax
            j = lmsb2i(jb,ilm,is,ia)
            i = lmsb2i(ib,ilm,is,ia)
            input(is,is) = gfpq(i,j,ia,ie)
          end do
!         now spin rotation
          output = matmul(matmul(ul,input),ur)
!         update the matrix elements of the onsite GF
          do js=1,nsmax
            j = lmsb2i(jb,ilm,js,ia)
            do is=1,nsmax
              i = lmsb2i(ib,ilm,is,ia)
              gfpq(i,j,ia,ie) = output(is,js)
            end do
          end do
        end do
        end do
      end do
!     ------------------------------------------------------------------
!                 Rotation of the scattering solutions
!     ------------------------------------------------------------------
!     Scattering wavefunctions should be diagonal: in-(l,m,s) --> out-(l,m,s)
!     Do NOT use a spin-dependent basis!
      do ilm=1,lmmax
        do ib=1,iwsusc(i2lm(2,ilm),1,ia)
!         ****   RIGHT SOLUTIONS   ****
!         first copy collinear elements
          input(:,:) = czero
          do is=1,nsmax
            ilms = lms2i(ilm,is)
            i = lmsb2i(ib,ilm,is,ia)
            input(is,is) = pzr(i,ilms,ia,ie)
          end do
!         now spin rotation
          output = matmul(matmul(ul,input),ur)
!         update the matrix elements of the right solution
          do js=1,nsmax
            jlms = lms2i(ilm,js)
            do is=1,nsmax
              i = lmsb2i(ib,ilm,is,ia)
              pzr(i,jlms,ia,ie) = output(is,js)
            end do
          end do
!         ****   LEFT SOLUTIONS   ****
!         first copy collinear elements
          input(:,:) = czero
          do is=1,nsmax
            ilms = lms2i(ilm,is)
            i = lmsb2i(ib,ilm,is,ia)
            input(is,is) = pzl(i,ilms,ia,ie)
          end do
!         now spin rotation
          output = matmul(matmul(ul,input),ur)
!         update the matrix elements of the left solution
          do is=1,nsmax
            ilms = lms2i(ilm,is)
            do js=1,nsmax
              j = lmsb2i(ib,ilm,js,ia)
              pzl(j,ilms,ia,ie) = output(is,js)
            end do
          end do
        end do
      end do
!     ------------------------------------------------------------------
!                    Rotation of the t-matrices    
!     ------------------------------------------------------------------
!     Onsite t-matrix should be diagonal in lms
      dtmat(:,:,ia) = czero
      do ilm=1,lmmax
!       first copy collinear elements
        input(:,:) = czero
        do is=1,nsmax
          ilms = lms2i(ilm,is)
          input(is,is) = tmatrix(ilms,ilms,ia,ie)
        end do
!       now spin rotation
        output = matmul(matmul(ul,input),ur)
!       update the matrix elements of the difference between t-matrices
        do js=1,nsmax
          jlms = lms2i(ilm,js)
          do is=1,nsmax
            ilms = lms2i(ilm,is)
            dtmat(ilms,jlms,ia) = output(is,js) - input(is,js)
          end do
        end do
      end do
!     ------------------------------------------------------------------
!   ++++++
    end do  ! atoms
!   ++++++
    call cpu_time(finish)
    time1 = time1 + finish - start
!   --------------------------------------------------------------------
!                    Updating the structural GF    
!   --------------------------------------------------------------------
    call cpu_time(start)
    call structural_gf(ie,nasusc,nlms,nalms,alms2i,dtmat,gstruct(:,:,ie))
!   Update total t-matrix: tmat = tmat + dtmat
    call zaxpy(nlms*nlms*nasusc,cone,dtmat,1,tmatrix(:,:,:,ie),1)
    call cpu_time(finish)
    time2 = time2 + finish - start
! ******
  end do  ! energy
! ******
  write(*,'("spinrot_tgf_sph: site diagonal time",f16.3," s")') time1
  write(*,'("spinrot_tgf_sph: structural GF time",f16.3," s")') time2
! ----------------------------------------------------------------------
  deallocate(dtmat)
! All done!
  end subroutine spinrot_tgf_sph
