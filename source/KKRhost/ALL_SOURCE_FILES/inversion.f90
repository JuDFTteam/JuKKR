! ************************************************************************
subroutine inversion(gllke, invmod, icheck)
! ************************************************************************
! This subroutine calculates the inversion of a matrix
! in 4 different ways depending on the form of the matrix

!     INVMOD = 0  ----> total inversion scheme
!     INVMOD = 1  ----> band matrix inversion scheme
!     INVMOD = 2  ----> corner band matrix inversion scheme
!     INVMOD = 3  ----> sparse matrix inversion scheme

! ------------------------------------------------------------------------
  implicit none

!     .. parameters ..
  include 'inc.p'
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************


! changed 3.11.99


  integer :: lmmaxd
  parameter (lmmaxd=(krel+korbit+1)*(lmaxd+1)**2)
  integer :: almd, ndim
  parameter (almd=naezd*lmmaxd, ndim=nprincd*lmmaxd)
  double complex :: ci, czero, cone
  parameter (ci=(0.d0,1.d0), czero=(0.d0,0.d0), cone=(1.d0,0.d0))

  double complex :: gllke(almd, almd), gdi(ndim, ndim, nlayerd), &
    gup(ndim, ndim, nlayerd), gdow(ndim, ndim, nlayerd)
  double complex, allocatable :: gtemp(:, :)
  integer :: i, i1, ip1, ii1, il1, ldi1, ip2, ii2, il2, ldi2, j, invmod
  integer :: lm1, lm2, info, ipvt(almd), nlayer
  integer :: icheck(naezd/nprincd, naezd/nprincd) 
! total matrix inversion
  external :: zgetrf, zgetrs, zcopy, invslab

  allocate (gtemp(almd,almd))

  nlayer = naezd/nprincd

  if (invmod==0) then 

!     write (6,*) '-------full inversion calculation--------'

    do i = 1, almd
      do j = 1, almd
        gtemp(i, j) = czero
        if (i==j) then
          gtemp(i, j) = cone
        end if
      end do
    end do




    call zgetrf(almd, almd, gllke, almd, ipvt, info)
    call zgetrs('N', almd, almd, gllke, almd, ipvt, gtemp, almd, info)
! slab or supercell
    call zcopy(almd*almd, gtemp, 1, gllke, 1)
! inversion


  else if ((invmod>=1) .and. (invmod<=2)) then 
!     this part now is correct also for    ! changes 20/10/99
!     supercell geometry : 20/10/99
!---> upper linear part
    do i1 = 1, nlayerd
      do ip1 = 1, nprincd
        do ip2 = 1, nprincd
          ii1 = (i1-1)*nprincd + ip1
          ii2 = (i1-1)*nprincd + ip2
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              ldi1 = lmmaxd*(ip1-1) + lm1
              il1 = lmmaxd*(ii1-1) + lm1
              ldi2 = lmmaxd*(ip2-1) + lm2
              il2 = lmmaxd*(ii2-1) + lm2
              gdi(ldi1, ldi2, i1) = gllke(il1, il2)
            end do
          end do
        end do
      end do
    end do


!---> lower linear part

    do i1 = 1, nlayerd
      do ip1 = 1, nprincd
        do ip2 = 1, nprincd
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              ldi1 = lmmaxd*(ip1-1) + lm1
              ldi2 = lmmaxd*(ip2-1) + lm2
              if (i1<=(nlayerd-1)) then
                ii1 = (i1-1)*nprincd + ip1
                ii2 = i1*nprincd + ip2
                il1 = lmmaxd*(ii1-1) + lm1
                il2 = lmmaxd*(ii2-1) + lm2
                gup(ldi1, ldi2, i1) = gllke(il1, il2)
              else
                ii1 = ip1
                ii2 = (nlayerd-1)*nprincd + ip2
                il1 = lmmaxd*(ii1-1) + lm1
                il2 = lmmaxd*(ii2-1) + lm2
                gdow(ldi1, ldi2, i1) = gllke(il1, il2)
              end if
            end do
          end do
        end do
      end do
    end do
!     end of the corrected part  20/10/99


    do i1 = 1, nlayerd
      do ip1 = 1, nprincd
        do ip2 = 1, nprincd
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              ldi1 = lmmaxd*(ip1-1) + lm1
              ldi2 = lmmaxd*(ip2-1) + lm2
              if (i1<=(nlayerd-1)) then
                ii1 = i1*nprincd + ip1
                ii2 = (i1-1)*nprincd + ip2
                il1 = lmmaxd*(ii1-1) + lm1
                il2 = lmmaxd*(ii2-1) + lm2
                gdow(ldi1, ldi2, i1) = gllke(il1, il2)
              else
                ii1 = (nlayerd-1)*nprincd + ip1
                ii2 = ip2
                il1 = lmmaxd*(ii1-1) + lm1
                il2 = lmmaxd*(ii2-1) + lm2
                gup(ldi1, ldi2, i1) = gllke(il1, il2)
              end if
            end do
          end do
        end do
      end do
    end do

!          write (6,*) '-------slab calculation--------'

    if (invmod==1) then
! supercell geometry inversion
      call invslab(gdi, gup, gdow, gllke, icheck)


!          write (6,*) '-------supercell calculation--------'
    else if (invmod==2) then 

      call invsupercell(gdi, gup, gdow, gllke, icheck)

! sparse matrix inversion

    end if
!     NOT YET IMPLEMENTED!!!!!!!!!

  else 




  end if

! ************************************************************************
  deallocate (gtemp)
! ************************************************************************
! This subroutine calculates the inversion of a matrix
  return
! in 4 different ways depending on the form of the matrix
end subroutine
