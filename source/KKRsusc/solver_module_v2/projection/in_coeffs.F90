  subroutine in_coeffs(is,ie,e,ek,de)
! Get projection coefficients from outsusc.dat file
  use global

  implicit none

! --> spin channel
  integer(kind=i4b), intent(in)  :: is
! --> energy point label
  integer(kind=i4b), intent(in)  :: ie
! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(out) :: e, ek, de
! -----------------------------------------------------------------
  integer(kind=i4b) :: ia, ih, nb, nb1, ib, jb, il, il1, is1
  real(kind=r8b)    :: x(1000)
  character*60 :: header

  read(iomain,'(a)') header ! ie, e, ek, de
!  write(iodb,*) header
  read(iomain,*) x(1:6)
!  write(iodb,'(6es16.8)') x(1:6)
  e  = cmplx(x(1),x(2))
  ek = cmplx(x(3),x(4))
  de = cmplx(x(5),x(6))

  if (.not.lhdio) return

  read(iomain,'(a)') header ! proj coeffs
!  write(iodb,*) header
  do ia=1,nasusc
    read(iomain,'(a)') header ! ia, i1
!    write(iodb,*) header
    do il=0,nlmax
!   this should cross-check with read_wfns
      nb = iwsusc(il,is,ia)
      if (nb > 0) then
        read(iomain,*) il1, is1, nb
!        write(iodb,'(3i4)') il1, is1, nb
        read(iomain,'(a)') header ! pz coeffs
!        write(iodb,*) header
        do ib=1,nb
          read(iomain,*) x(1:2)
!          write(iodb,'(100es16.8)') x(1:2)
          pzc(ib,il,is,ia) = cmplx(x(1),x(2))
        end do
!     ---------------------------------------------------------------
        if (isra == 1) then
          read(iomain,'(a)') header ! fz coeffs
!          write(iodb,*) header
          do ib=1,nb
            read(iomain,*) x(1:2)
!            write(iodb,'(100es16.8)') x(1:2)
            fzc(ib,il,is,ia) = cmplx(x(1),x(2))
          end do
        end if
!     ---------------------------------------------------------------
      end if
    end do
  end do
  read(iomain,'(a)') header ! onsite coefficients
!  write(iodb,*) header
  do ia=1,nasusc
    read(iomain,'(a)') header ! ia, i1
!    write(iodb,*) header
    do il=0,nlmax
!   this was read in before
      nb = iwsusc(il,is,ia)
!      write(*,*) "nb=", nb
      if (nb > 0) then
        read(iomain,*) il1, is1, nb1
!        write(iodb,'(3i4)') il1, is1, nb
        read(iomain,'(a)') header ! pq coeff
!        write(iodb,*) header
        do ib=1,nb
          read(iomain,*) x(1:2*nb)
!          write(iodb,'(100es16.8)') x(1:2)
          do jb=1,nb
            pqc(jb,ib,il,is,ia) = cmplx(x(2*jb-1),x(2*jb))
          end do
        end do
!     -------------------------------------------------------------
        if (isra == 1) then
          read(iomain,'(a)') header ! ps coeff
!          write(iodb,*) header
          do ib=1,nb
            read(iomain,*) x(1:2*nb)
!            write(iodb,'(100es16.8)') x(1:2)
            do jb=1,nb
              psc(jb,ib,il,is,ia) = cmplx(x(2*jb-1),x(2*jb))
            end do
          end do
          read(iomain,'(a)') header ! fq coeff
!          write(iodb,*) header
          do ib=1,nb
            read(iomain,*) x(1:2*nb)
!            write(iodb,'(100es16.8)') x(1:2)
            do jb=1,nb
              fqc(jb,ib,il,is,ia) = cmplx(x(2*jb-1),x(2*jb))
            end do
          end do
          read(iomain,'(a)') ! fs coeff
!          write(iodb,*) header
          do ib=1,nb
            read(iomain,*) x(1:2*nb)
!            write(iodb,'(100es16.8)') x(1:2)
            do jb=1,nb
              fsc(jb,ib,il,is,ia) = cmplx(x(2*jb-1),x(2*jb))
            end do
          end do
        end if
!     -------------------------------------------------------------
      end if
    end do
  end do
! coefficients read in
  noregcoeffs = .false.
  noirrcoeffs = .false.
! All done
  end subroutine in_coeffs

