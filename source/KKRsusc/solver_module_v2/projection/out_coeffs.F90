  subroutine out_coeffs(is,ie,e,ek,de)
! Output projection coefficients to outsusc.dat file
  use global

  implicit none

! --> spin channel
  integer(kind=i4b), intent(in) :: is
! --> energy point label
  integer(kind=i4b), intent(in) :: ie
! --> energy point value, its square-root, integration weight
  complex(kind=c8b), intent(in) :: e, ek, de
! -----------------------------------------------------------------
  real(kind=r8b), parameter :: pi = 4.d0*atan(1.d0)
  integer(kind=i4b) :: ia, ih, nb, ib, il

  if (.not.lhdio) then
    esusc(ie) = e; eksusc(ie) = ek
    if (ikkr == 1 ) then  ! old impurity code
      desusc(ie) = -pi*de
    else                  ! KKRFLEX
      desusc(ie) = -pi*de/nsmax
    end if
    if (.not.lscfmesh) then
      escf(ie) = esusc(ie); descf(ie) = desusc(ie)
      if (ie == nesusc) then
        write(*,'(/," de sums to=",2f12.6,/)') sum(descf(1:nescf))
        efscf = real(escf(ie))  ! Fermi energy
      end if
    end if
    call save_coeffs(ie,is,.false.)
    return
  end if
  write(iomain,'(" ie, e, ek, de=",i4)') ie
  write(iomain,'(6es16.8)') e, ek, de
  write(iomain,'(" Projection coefficients:")')
  do ia=1,nasusc
    ih = iasusc(ia)
    write(iomain,'(" ia=",i4,"  i1=",i4)') ia, ih
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
      if (nb > 0) then
        write(iomain,'(4i4)') il, is, nb
        write(iomain,*) " pz coeff"
        do ib=1,nb
          write(iomain,'(100es16.8)') pzc(ib,il,is,ia)
        end do
!     ---------------------------------------------------------------
        if (isra == 1) then
          write(iomain,*) " fz coeff"
          do ib=1,nb
            write(iomain,'(100es16.8)') fzc(ib,il,is,ia)
          end do
        end if
!     ---------------------------------------------------------------
      end if
    end do
  end do
  write(iomain,'(" Onsite coefficients:")')
  do ia=1,nasusc
    ih = iasusc(ia)
    write(iomain,'(" ia=",i4,"  i1=",i4)') ia, ih
    do il=0,nlmax
      nb = iwsusc(il,is,ia)
      if (nb > 0) then
        write(iomain,'(4i4)') il, is, nb
        write(iomain,*) " pq coeff"
        do ib=1,nb
          write(iomain,'(100es16.8)') pqc(1:nb,ib,il,is,ia)
        end do
!     ---------------------------------------------------------------
        if (isra == 1) then
          write(iomain,*) " ps coeff"
          do ib=1,nb
            write(iomain,'(100es16.8)') psc(1:nb,ib,il,is,ia)
          end do
          write(iomain,*) " fq coeff"
          do ib=1,nb
            write(iomain,'(100es16.8)') fqc(1:nb,ib,il,is,ia)
          end do
          write(iomain,*) " fs coeff"
          do ib=1,nb
            write(iomain,'(100es16.8)') fsc(1:nb,ib,il,is,ia)
          end do
        end if
!     ---------------------------------------------------------------
      end if
    end do
  end do
! All done
  end subroutine out_coeffs
