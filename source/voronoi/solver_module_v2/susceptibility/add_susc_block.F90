  subroutine add_susc_block(ia2,ja2,suscblock,work,chi0)
! Adds (ia2,ja2) block to KS suscepbilitity
! Using density basis
  use global

  implicit none

! Atom pair
  integer(kind=i4b), intent(in)    :: ia2, ja2
! KS susc block in GF basis
  complex(kind=c8b), intent(inout) :: suscblock(ngfmax,ngfmax)
! Work array
  complex(kind=c8b), intent(inout) :: work(ngfmax,ngfmax)
! Storage for the full susceptibility
  complex(kind=c8b), intent(inout) :: chi0(ndensum,ndensum)
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0)
  integer(kind=i4b) :: ia, ja
  integer(kind=i4b) :: iq, iq1, jq, jq1, ndeni0, ndenj0, ngfi0, ngfj0, ndeni1, ndenj1, ngfi1, ngfj1

  ja = iasusc2(ja2)
  ngfj0  = sum(nalmsbgf(1:ja-1))
  ngfj1  = sum(nalmsbgf(1:ja))
  ndenj0 = sum(nalmsbden(1:ja2-1))
  ndenj1 = sum(nalmsbden(1:ja2))
  ia = iasusc2(ia2)
  ngfi0  = sum(nalmsbgf(1:ia-1))
  ngfi1  = sum(nalmsbgf(1:ia))
  ndeni0 = sum(nalmsbden(1:ia2-1))
  ndeni1 = sum(nalmsbden(1:ia2))
! multiply the KS susc from the left by dengaunt
  work  = czero
  do jq1=1,ngfj1-ngfj0
    do iq=1,ndeni1-ndeni0
      do iq1=1,ngfi1-ngfi0
        if (abs(dengaunt(iq1,iq,ia2)) > ylmtol) then
          work(iq,jq1) = work(iq,jq1) + dengaunt(iq1,iq,ia2)*suscblock(iq1,jq1)
        end if
      end do
    end do
  end do
! multiply the KS susc from the right by dengaunt
  suscblock = czero
  do jq=1,ndenj1-ndenj0
    do jq1=1,ngfj1-ngfj0
      if (abs(dengaunt(jq1,jq,ja2)) > ylmtol) then
        do iq=1,ndeni1-ndeni0
          suscblock(iq,jq) = suscblock(iq,jq) + work(iq,jq1)*dengaunt(jq1,jq,ja2)
        end do
      end if
    end do
  end do
! copy block to storage
  do jq=1,ndenj1-ndenj0
    do iq=1,ndeni1-ndeni0
      chi0(iq+ndeni0,jq+ndenj0) = chi0(iq+ndeni0,jq+ndenj0) + suscblock(iq,jq)
    end do
  end do
! All done!
  end subroutine add_susc_block
