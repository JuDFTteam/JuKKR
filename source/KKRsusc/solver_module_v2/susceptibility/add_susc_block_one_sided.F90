  subroutine add_susc_block_one_sided(ia2,ja2,suscblock,work,chi0)
! Adds (ia2,ja2) block to KS correlation function
! Using gf basis (left) and  density basis (right)
  use global

  implicit none

! Atom pair
  integer(kind=i4b), intent(in)    :: ia2, ja2
! KS susc block in GF basis
  complex(kind=c8b), intent(inout) :: suscblock(ngfmax,ngfmax)
! Work array
  complex(kind=c8b), intent(inout) :: work(ngfmax,ngfmax)
! Storage for the full susceptibility
  complex(kind=c8b), intent(inout) :: chi0(ngradsum,ndensum)
! ----------------------------------------------------------------------
  complex(kind=c8b), parameter :: czero = (0.d0,0.d0), cone = (1.d0,0.d0)
  integer(kind=i4b) :: ia, ja, i4(4),i3(3), jlm1, jlm2, ilm1, ilm2, j1, i, q1, q2
  integer(kind=i4b) :: iq, iq1, jq, jq1, ndeni0, ndenj0, ngfi0, ngfj0, ndeni1, ndenj1, ngfi1, ngfj1, ngradi0, ngradi1, ngradj0, ngradj1

  ja = iasusc2(ja2)
  ngfj0  = sum(nalmsbgf(1:ja-1))
  ngfj1  = sum(nalmsbgf(1:ja))
  ndenj0 = sum(nalmsbden(1:ja2-1))
  ndenj1 = sum(nalmsbden(1:ja2))
  ngradj0 = sum(nalmsbgrad(1:ja2-1))
  ngradj1 = sum(nalmsbgrad(1:ja2))

  ia = iasusc2(ia2)
  ngfi0  = sum(nalmsbgf(1:ia-1))
  ngfi1  = sum(nalmsbgf(1:ia))
  ndeni0 = sum(nalmsbden(1:ia2-1))
  ndeni1 = sum(nalmsbden(1:ia2))
  ngradi0 = sum(nalmsbgrad(1:ia2-1))
  ngradi1 = sum(nalmsbgrad(1:ia2))

! multiply the KS correlation func from the right by dengaunt
!  work  = czero
!  do jq=1,ndenj1-ndenj0
!    do jq1=1,ngfj1-ngfj0
!      j1 = jq1 + ngfj0
!      i3 = i2almsbgf(:,j1)
!      q1 = i3(1); q2 = i3(2)
!      i3 = i2lmsb(:,q1,ja)
!      jlm1 = i3(2)
!      i3 = i2lmsb(:,q2,ja)
!      jlm2 = i3(2)
!      if((jlm1 == 5 .OR. jlm1 == 9) .AND. (jlm2 == 5 .OR. jlm2 ==9)) then
!      if (abs(dengaunt(jq1,jq,ja2)) > ylmtol) then
!        do iq=1,ngfi1-ngfi0
!          i = iq + ngfi0
!          i3 = i2almsbgf(:,i)
!          q1 = i3(1); q2 = i3(2)
!          i3 = i2lmsb(:,q1,ia)
!          ilm1 = i3(2)
!          i3 = i2lmsb(:,q2,ia)
!          ilm2 = i3(2)
!          if((ilm1 == 5 .OR. ilm1 == 9) .AND. (ilm2 == 5 .OR. ilm2 ==9)) then
!            work(iq,jq) = work(iq,jq) + suscblock(iq,jq1)*dengaunt(jq1,jq,ja2)
!          end if
!        end do
!      end if
!      end if
!    end do
!  end do
!-------------------------------------------------------------------------------
  work  = czero
  do jq=1,ndenj1-ndenj0
    do jq1=1,ngfj1-ngfj0
      if (abs(dengaunt(jq1,jq,ja2)) > ylmtol) then
        do iq=1,ngfi1-ngfi0
          work(iq,jq) = work(iq,jq) + suscblock(iq,jq1)*dengaunt(jq1,jq,ja2)
        end do
      end if
    end do
  end do

  suscblock(:,:)=czero
  suscblock(:,:)=work(:,:)

! copy block to storage
  do jq=1,ndenj1-ndenj0
    do iq=1,ngradi1-ngradi0
      chi0(iq+ngradi0,jq+ndenj0) = chi0(iq+ngradi0,jq+ndenj0) + suscblock(iq,jq)
    end do
  end do
! All done!

!--------------------------------------------------------------------------------
!  do jq=1,ngradj1-ngradj0
!    do iq=1,ngradi1-ngradi0
!      chi0(iq+ngradi0,jq+ngradj0) = chi0(iq+ngradi0,jq+ngradj0) + suscblock(iq,jq)
!    end do
!  end do

  end subroutine add_susc_block_one_sided
