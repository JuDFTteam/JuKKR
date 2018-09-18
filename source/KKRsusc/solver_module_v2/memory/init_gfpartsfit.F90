  subroutine init_gfpartsfit(my_rank)
! Allocates arrays for different parts of projected GF
  use global

  implicit none

  integer(kind=i4b) :: is, il, im, ib, i, lm, ia, ilms, jlms, my_rank
  real(kind=r8b)    :: ram


! -----------------------------------------------------------------------
!                             Pointers               
! -----------------------------------------------------------------------
! --> lmsb, new version: block sizes atom dependent
  allocate(lmsb2i(nbmax,lmmax,nsmax,nasusc))
  allocate(i2lmsb(3,nlmsb,nasusc))
  allocate(nlmsba(nasusc))
  do ia=1,nasusc
    i = 0
    do lm=1,lmmax
      do is=1,nsmax
        do ib=1,iwsusc(i2lm(2,lm),is,ia)
          i = i + 1
          lmsb2i(ib,lm,is,ia) = i
          i2lmsb(:,i,ia) = (/ib,lm,is/)
!          write(*,*) i, ib, is, lm, ia
        end do
      end do
    end do
    nlmsba(ia) = i
    if (my_rank == 0 .and. loutbasis) then
      write(*,'("nlmsb for ia=",i4," is ",i4)') ia, i
    end if ! my_rank
  end do
  if (my_rank == 0 .and. loutbasis) then
    write(*,'("nalmsb should be ",i8)') sum(nlmsba(1:nasusc))
  end if ! my_rank
! -----------------------------------------------------------------------
!   Storage for structural GF and coefficients in full projection basis
! -----------------------------------------------------------------------
! how much ram to be used in these arrays
  ram = 16.d0*2*nlms*nlmsb*nasusc*nesusc              ! size of reg coeffs arrays (lhs and rhs)
  ram = ram + 8.d0*nlmsb*nlmsb*nasusc                 ! overlap
  ram = ram + 16.d0*nlmsb*nlmsb*nasusc*nesusc         ! size of single-site GF array
  if (lsusc) ram = ram + 16.d0*nlmsb*nlmsb*nasusc     ! size of new potential and xc kernel arrays
! extra storage for small components
  if (isra == 1) then                                 
    ram = ram + 16.d0*2*nlmsb*nlms*nasusc*nesusc      ! fz left and right
    ram = ram + 16.d0*3*nlmsb*nlmsb*nasusc*nesusc     ! ps, fq, fs
  end if
  ram = ram/(1024.d0**2)
! time to cry 
  if (my_rank == 0) then 
    write(*,'("init_arrays: GF parts   RAM=",f16.3," MB")') ram
  end if ! my_rank
  allocate(overlap(nlmsb,nlmsb,nasusc))               ! overlaps
  overlap = 0.d0
  allocate(vlmsbgf(nlmsb,nlmsb,nasusc))    ! potential in the product basis
  vlmsbgf = 0.d0
! -----------------------------------------------------------------------
! Storage for coefficients of fitted GF
  if (lfit .and. ifit == 2) then
!   time to cry again
    ram = 16.d0*(1+numd+dend)*nlmsb*nlmsb*nasusc*nasusc
    ram = ram/(1024.d0**2)
    if (my_rank == 0) then 
      write(*,'("init_arrays: GF fit     RAM=",f16.3," MB")') ram
    end if ! my_rank
    allocate(gffit(numd+dend,nlmsb,nlmsb,nasusc,nasusc))
  end if
! -----------------------------------------------------------------------
! All done!
  end subroutine init_gfpartsfit
