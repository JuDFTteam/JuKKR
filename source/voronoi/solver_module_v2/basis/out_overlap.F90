  subroutine out_overlap(ia,ja)
! Mutual overlaps between basis functions for ia and ja
! Only gives right answer if the same radial mesh applies to everyone
  use global

  implicit none

! --> susc atom
  integer(kind=i4b), intent(in) :: ia, ja
! -----------------------------------------------------------------
  integer(kind=i4b) :: ib, jb, il, jl, is, js, nr, nb, mb
  complex(kind=c8b) :: work(nrmax)
  complex(kind=c8b) :: pp(nbmax)
  complex(kind=c8b), external :: radint

  nr = nrpts(ja)
  write(iodb,'("overlaps: ia, ja, is, js, il, jl")')
  do js=1,issusc(ja)
  do is=1,js
    do jl=0,nlmax
    do il=0,jl
      nb = iwsusc(jl,js,ja)
      mb = iwsusc(il,is,ia)
      if (nb > 0 .and. mb > 0) then
        write(iodb,'("block=",6i4)') ia, ja, is, js, il, jl
        do jb=1,nb
          pp = 0.d0
          do ib=1,mb
            work = phiref(:,jb,jl,js,ja)*phiref(:,ib,il,is,ia)
            pp(ib) = radint(nr,work(1:nr),drmesh(1:nr,ia),npanat(ia),ircutat(:,ia))
          end do
          where (abs(pp(1:mb)) < gstol) pp = 0.d0
          write(iodb,'(100es10.3)') real(pp(1:mb))
        end do
      end if
    end do
    end do
  end do
  end do
! All done!
  end subroutine out_overlap
