  subroutine selfenergy_correction()
! Adding the self-energy to the GF

  implicit none

! complex parameters
  complex(kind=c8b), parameter :: cplus  = ( 1.d0, 0.d0)
  complex(kind=c8b), parameter :: cminus = (-1.d0, 0.d0)
  complex(kind=c8b), parameter :: czero  = ( 0.d0, 0.d0)
! change in t-matrix
  complex(kind=c8b) :: tmat(nlms,nlms,nasusc)
  complex(kind=c8b) :: tmatl(nlms,nlms), tmatr(nlms,nlms)
! SOC potential
  complex(kind=c8b) :: vsigmab(nlmsb,nlmsb)
! auxiliary
  complex(kind=c8b) :: work1(nlmsb,nlmsb), work2(nalms,nalms)
  complex(kind=c8b) :: gfpqsigma(nlmsb,nlmsb)
  complex(kind=c8b) :: pzlsigma(nlmsb,nlms), pzrsigma(nlmsb,nlms)
  complex(kind=c8b) :: edummy
  real(kind=r8b)    :: maxelem, start, finish
! -----------------------------------------------------------------
  integer(kind=i4b) :: i, j, ia, is, js, ja, k
  integer(kind=i4b) :: jlm, ilm, ilms, jlms, klms
  integer(kind=i4b) :: ib, jb, ie, i3(3), j3(3)
  integer(kind=i4b) :: ipiv(nlmsb), jpiv(nalms), info

  call cpu_time(start)
! loop over energies
  do ie=1,nesusc
    write(iodb,*) "Self-energy correction ie=", ie
! loop over atoms
    do ia=1,nasusc
!   construct potential
      edummy = esusc(ie)
!   compute the matrix elements in the basis
      call build_vsigmab(ia,ie,vsigmab)
!   ********************************************************************
!              solve Lippmann-Schwinger and Dyson equations
!   ********************************************************************
!     --> construct 1 - G0.Sigma
      work1 = 0.d0
      do i=1,nlmsba(ia)
        work1(i,i) = 1.d0
      end do
      call zgemm('N','N',nlmsba(ia),nlmsba(ia),nlmsba(ia),cminus,gfpq(:,:,ia,ie),nlmsb,vsigmab,nlmsb,cplus,work1,nlmsb)
!     --> LU decomposition of 1 - G0.Sigma
      call zgetrf(nlmsba(ia),nlmsba(ia),work1,nlmsb,ipiv,info)
      if (info /= 0) then
        write(*,*) "LU fail 1 - G0.Sigma at ie,ia=", ie, ia, info
        stop
      end if
!     --> new single-site GF: (1 - G0.Sigma).G = G0
      gfpqsigma = gfpq(:,:,ia,ie)
      call zgetrs('N',nlmsba(ia),nlmsba(ia),work1,nlmsb,ipiv,gfpqsigma,nlmsb,info)
      if (info /= 0) then
        write(*,*) "Self-energy single-site GF failed at ie,ia=", ie, ia
        stop
      end if
!     --> solve Lipp-Schw for rhs wfn: (1 - G0.Sigma).R^r = R0^r
!       nlms are the boundary conditions
      pzrsigma = pzr(:,:,ia,ie)  ! the initial rhs solution is needed in t^l
      call zgetrs('N',nlmsba(ia),nlms,work1,nlmsb,ipiv,pzrsigma,nlmsb,info)
      if (info /= 0) then
        write(*,*) "Self-energy rhs Lipp Schw failed at ie,ia=", ie, ia
        stop
      end if
!      write(iodb,*) "non-zero elements of rhs wfn"
!      write(iodb,'("  ib ilm  is jlm  js      pzr       pzrldau")')
      maxelem = maxval(abs(pzr(:,:,ia,ie)))
!      write(iodb,'("maxelem  ",es14.6)') maxelem
      do j=1,nlms
        do i=1,nlmsba(ia)
          if (abs(pzrsigma(i,j)) > sigmatol*maxelem) then
!        write(iodb,'(5i4,4es14.6)') i2lmsb(:,i), i2lms(:,j), pzr(i,j,ia,ie), pzrldau(i,j)
          else
            pzrsigma(i,j) = 0.d0
          end if
        end do
      end do
!     --> change in t-matrix: t^r = R0^l.Sigma.R^r
      tmatr = 0.d0; work1 = 0.d0
!     Sigma.R^r -> work
      call zgemm('N','N',nlmsba(ia),nlms,nlmsba(ia),cplus,vsigmab,nlmsb,pzrsigma,nlmsb,czero,work1,nlmsb)
!     R0^l.work -> t^r
      call zgemm('T','N',nlms,nlms,nlmsba(ia),cplus,pzl(:,:,ia,ie),nlmsb,work1,nlmsb,czero,tmatr,nlms)
!      write(iodb,*) "tmatr for ia=", ia
!      do j=1,nlms
!        write(iodb,'(100es10.2)') tmatr(:,j)
!      end do
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     --> construct 1 - G0^T.Sigma^T
      work1 = 0.d0
      do i=1,nlmsba(ia)
        work1(i,i) = 1.d0
      end do
      call zgemm('T','T',nlmsba(ia),nlmsba(ia),nlmsba(ia),cminus,gfpq(:,:,ia,ie),nlmsb,vsigmab,nlmsb,cplus,work1,nlmsb)
!     --> LU decomposition of 1 - G0^T.Sigma^T
      call zgetrf(nlmsba(ia),nlmsba(ia),work1,nlmsb,ipiv,info)
      if (info /= 0) then
        write(*,*) "LU fail 1 - G0^T.Sigma^T at ie,ia=", ie, ia, info
        stop
      end if
!     --> solve Lipp-Schw for left hand side wfn: (1 - G0^T.Sigma^T).R^l = R0^l
!       nlms are the boundary conditions
      pzlsigma = pzl(:,:,ia,ie)
      call zgetrs('N',nlmsba(ia),nlms,work1,nlmsb,ipiv,pzlsigma,nlmsb,info)
      if (info /= 0) then
        write(*,*) "LDA+U lhs Lipp Schw failed at ie,ia=", ie, ia
        stop
      end if
!      write(iodb,*) "non-zero elements of lhs wfn"
!      write(iodb,'("  ib ilm  is jlm  js      pzl       pzlldau")')
      maxelem = maxval(abs(pzl(:,:,ia,ie)))
!      write(iodb,'("maxelem  ",es14.6)') maxelem
      do j=1,nlms
        do i=1,nlmsba(ia)
          if (abs(pzlsigma(i,j)) > sigmatol*maxelem) then
!            write(iodb,'(5i4,4es14.6)') i2lmsb(:,i), i2lms(:,j), pzl(i,j,ia,ie), pzlldau(i,j)
          else
            pzlsigma(i,j) = 0.d0
          end if
        end do
      end do
!     --> change in t-matrix: t^l = R^l.Sigma.R0^r
      tmatl = 0.d0; work1 = 0.d0
!     R^l.Sigma -> work
      call zgemm('T','N',nlms,nlmsba(ia),nlmsba(ia),cplus,pzlsigma,nlmsb,vsigmab,nlmsb,czero,work1,nlmsb)
!     work.R0^r -> t^l
      call zgemm('N','N',nlms,nlms,nlmsba(ia),cplus,work1,nlmsb,pzr(:,:,ia,ie),nlmsb,czero,tmatl,nlms)
!      write(iodb,*) "tmatl for ia=", ia
!      do j=1,nlms
!        write(iodb,'(100es10.2)') tmatl(:,j)
!      end do
      maxelem = maxval(abs(tmatl - tmatr))
      write(iodb,'("|tmatl - tmatr| for ia=",i4," is ",es14.6)') ia, maxelem
!      do j=1,nlms
!        write(iodb,'(100es10.2)') abs(tmatl(:,j) - tmatr(:,j))
!      end do
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     --> update wfns
      pzl(:,:,ia,ie) = pzlsigma
      pzr(:,:,ia,ie) = pzrsigma
!     --> update t-matrix
!      if (ia == 1) tmat(:,:,ia) = tmatr
      tmat(:,:,ia) = tmatr
!     --> update single-site GF
      gfpq(:,:,ia,ie) = gfpqsigma
    end do
!   new structural GF: (1 - G0^s.t).G^s = G0^s
!   --> construct 1 - G0^s.t and G0^s
    do ja=1,nasusc
    do jlms=1,nlms
      j = alms2i(jlms,ja)
      do ia=1,nasusc
      do ilms=1,nlms
        i = alms2i(ilms,ia)
        work2(i,j) = 0.d0
        if (i == j) work2(i,j) = 1.d0
!     matrix multiplication
        do klms=1,nlms
          k = alms2i(klms,ja)
          work2(i,j) = work2(i,j) - gstruct(i,k,ie)*tmat(klms,jlms,ja)
        end do
      end do
      end do
    end do
    end do
!     --> solve Dyson equation
    call zgesv(nalms,nalms,work2,nalms,jpiv,gstruct(:,:,ie),nalms,info)
    if (info /= 0) then
      write(*,*) "Self-energy struct Dyson failed at ie=", ie
      stop
    end if
  end do
! Flag SOC as applied
!  ldau_applied = .true.
  call cpu_time(finish)
  write(*,'(/," Self-energy correction time=",f10.3," s",/)') finish - start
! All done!
  end subroutine selfenergy_correction

