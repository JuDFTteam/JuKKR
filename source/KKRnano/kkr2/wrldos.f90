!----------------------------------------------------------------------
!>    Writes local density of states.
      subroutine write_ldos(den, ez, lmaxd1, iemxd, ititle, efermi, e1, e2, alatc, tK, nspin, i1)
      implicit none
      double precision, intent(in) :: alatc, e1, e2, efermi, tK
      integer, intent(in) :: iemxd, lmaxd1, nspin, i1
      double complex, intent(in) :: den(0:lmaxd1,iemxd,nspin), ez(iemxd)
      double precision :: dostot(0:lmaxd1,2), pdostot(0:lmaxd1,2)
      integer, intent(in) :: ititle(20,nspin)

!     .. locals ..
      double precision, parameter :: kB=0.6333659d-5
      double complex :: wez(iemxd) ! warning: wez is uninitialized !
      double precision :: dos, dossgn, efctor, pi
      integer :: ielast, ie, ipot, ispin, l
      character(len=16) :: fname

!     logical, external :: test

      pi = 4.d0*atan(1.d0)
      efctor = 1.d0
!     if (test('ev      ')) efctor = 13.6058d0
!
      ielast = iemxd
!
! initialize dostot
      dostot  = 0.d0
      pdostot = 0.d0

!=======================================================================
! write dos to file dos.i1.dat - begin
!=======================================================================
      write(unit=fname, fmt="(a,i4.4,a)") 'DOS.',i1,'.dat'
      open(48, file=fname, form='formatted', action='write')

      do ispin = 1, nspin
        ipot = nspin*(i1-1) + ispin
        dossgn = 1.d0 ; if (ispin /= nspin) dossgn = -1.d0

        write (48,fmt=9010) ititle(1:19,ispin)
        write (48,fmt=9020) i1
        write (48,fmt=9030) ispin, ielast, e1, e2, efermi, efctor
        write (48,fmt=9040) efermi
        write (48,fmt=9050) tK, pi*kB*tK, alatc

        do ie = 1, ielast
          dos = 0.d0
          do l = 0, lmaxd1
            dos = dos - 2.d0 * dimag(den(l,ie,ispin))/(pi*nspin)
            dostot(l,ispin) = dostot(l,ispin) + dimag(wez(ie)*den(l,ie,ispin))
          enddo
          write(48,fmt=9060) dble(ez(ie))*efctor, -2.d0*dimag(den(:,ie,ispin))*dossgn/(efctor*pi*nspin), dos*dossgn/efctor
        enddo ! ie

        write(48,fmt=9070) dostot(:,ispin)/(efctor*nspin)
        if (ispin /= nspin) write(48,fmt=9000)
      enddo ! ispin
      close(48)

!=======================================================================
! write dos to file dos.i1.dat - end
!=======================================================================
      return
 9000 format ('&')
 9010 format ('#',19a4)
 9020 format ('# I1    :',I8)
 9030 FORMAT ('# ISPIN :',I8,'   IELAST :',I5,/,'# E1,E2 :',2f12.5,' EFERMI :',f12.5,'   EFCTR',f10.6)
 9040 FORMAT ('# FERMI :',f12.5)
 9050 FORMAT ('# tK    =',f8.1,'   Kelvin =',3p,f8.3,' mRyd',0p,/,'# ALAT   :',f12.5)
 9060 FORMAT (1p,8e15.7)
 9065 FORMAT (1p,16d15.7)
 9070 FORMAT ('# Integrated DOS ',1p,d10.3,7d11.3)
 9080 format ('&')
      endsubroutine ! write_ldos

!> writes complex.dos file (complex density of states).
!     subroutine wrldos(den, ez, wez, lmaxd1, iemxd, npotd, ititle, efermi, e1, e2, alatc, tK, nspin, naez, ielast, i1, dostot) ! remove wez and npotd
      subroutine wrldos(den, ez, lmaxd1, iemxd, ititle, efermi, e1, e2, alatc, tK, nspin, naez, ielast, i1, dostot) ! remove wez and npotd
      implicit none
!     integer, intent(in) :: npotd ! todo: removed
!     double complex, intent(in) :: wez(iemxd) ! todo: removed
      double precision, intent(in) :: alatc, e1, e2, efermi, tK
      integer, intent(in) :: ielast, iemxd, lmaxd1, naez, nspin
      double complex, intent(in) :: den(0:lmaxd1,iemxd,nspin), ez(iemxd)
      integer, intent(in) :: ititle(20,nspin)
      double precision, intent(out) :: dostot(0:lmaxd1,2)

!     .. locals ..
      double precision, parameter :: kB=0.6333659d-5
      double complex   :: doscmplx
      double precision :: dossgn, pi, pispininv
      integer          :: i1, ie, ispin
      double precision, parameter :: efactor = 1.d0, efactorinv = 1.d0/efactor

      pi = 4.d0*atan(1.d0)
      pispininv = 1.d0/(pi*nspin)
      
      dostot = 0.d0

!
! open file complex.dos - kept for correspondence to complexdos3.f
!
      if(i1 == 1) then
        open(49,file='complex.dos', form='formatted', action='write')
        write(49,*) naez*nspin
        write(49,*) ielast
        write(49,*) lmaxd1
      endif

!=======================================================================
! write complex dos to file complex.dos - begin
!=======================================================================
!
        do ispin = 1,nspin
          dossgn = 1.d0
          if (ispin /= nspin) dossgn = -1.d0

          write(49,fmt=9010) ititle(1:19,ispin)
          write(49,fmt=9020) i1
          write(49,fmt=9030) ispin, ielast, e1, e2, efermi, efactor
          write(49,fmt=9040) efermi
          write(49,fmt=9050) tK, pi*kB*tK, alatc

          do ie = 1,ielast
            doscmplx = -2.d0*sum(den(:,ie,ispin))*pispininv
            write(49,fmt=9065) ez(ie)*efactor, -2.d0*den(:,ie,ispin)*dossgn*efactorinv*pispininv, doscmplx*dossgn*efactorinv
          enddo ! ie

          if (ispin /= nspin .or. i1 /= naez) write(49,fmt=9000)
        enddo ! ispin
        if(i1 == naez) close(49)
!
!=======================================================================
! write complex dos to file complex.dos - end
!=======================================================================
      return
!
 9000 format ('&')
 9010 format ('#',19a4)
 9020 FORMAT ('# I1    :',I8)
 9030 FORMAT ('# ISPIN :',I8,'   IELAST :',I5,/,'# E1,E2 :',2f12.5,' EFERMI :',f12.5,'   EFCTR',f10.6)
 9040 FORMAT ('# FERMI :',f12.5)
 9050 FORMAT ('# tK    =',f8.1,'   Kelvin =',3p,f8.3,' mRyd',0p,/,'# ALAT   :',f12.5)
 9060 FORMAT (1p,8e15.7)
 9065 FORMAT (1p,16d15.7)
 9070 FORMAT ('# Integrated DOS ',1p,d10.3,7d11.3)
 9080 FORMAT ('&')
      endsubroutine
