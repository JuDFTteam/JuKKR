
!================================================================
!> Collects contributions to the rms errors from all sites and
!> prints rms errors.
!
!> output of a) rms-error
!
! called by main2
!================================================================

subroutine RMSOUT_com(RMSAVQ,RMSAVM,ITER,NSPIN,NAEZ, &
MYLRANK, communicator)

  implicit none

  double precision RMSAVM,RMSAVQ

  integer ITER,NSPIN,NAEZ
  !     ..
  !     .. MPI variables ..
  !     .. L-MPI ..
  integer      MYLRANK, communicator

  double precision RMSQ, RMSM

  call allreduceRMS_com(RMSQ, RMSM, RMSAVQ,RMSAVM,NAEZ, communicator)

  ! ================== MYRANK.EQ.0 =======================================
  if(MYLRANK.eq.0) then

    call printRMSerror(ITER, NSPIN, RMSM, RMSQ)

  end if
  ! ============= MYRANK.EQ.0 ============================================

!  VFORM = .false.
!  if (ITER.eq.SCFSTEPS) then
!
!    inquire(file='VFORM',exist=VFORM)
!
!    if (VFORM) then
!
!      call writeFormattedPotential(fermiEnergy,VBC,NSPIN, &
!      KXC,LPOT,A,B,IRC, &
!      VINS,VISP,DRDI,IRNS,R,RWS,RMT,ALAT, &
!      ECORE,LCORE,NCORE,ZAT,ITITLE, &
!      MYLRANK, &
!      irmd, irnsd)
!
!    end if
!  end if

end

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Communicate the contributions of the rms error of charge and magnetisation
!> density and return results in RMSQ and RMSM.
subroutine allreduceRMS_com(RMSQ, RMSM, &  ! output
RMSAVQ_local,RMSAVM_local,NAEZ, & !in
communicator) !in

  implicit none

  INCLUDE 'mpif.h'

  double precision RMSAVM_local,RMSAVQ_local, &
  WORK1(2),WORK2(2)

  integer NAEZ


  !     .. local scalars ..
  double precision RMSQ,RMSM
  !     ..
  !     .. MPI variables ..
  integer communicator
  integer ierr

  !****************************************************** MPI COLLECT DATA

  WORK1(1) = RMSAVQ_local
  WORK1(2) = RMSAVM_local

  call MPI_ALLREDUCE(WORK1,WORK2,2, &
  MPI_DOUBLE_PRECISION,MPI_SUM,communicator, &
  IERR)

  RMSQ = SQRT(WORK2(1)/NAEZ)
  RMSM = SQRT(WORK2(2)/NAEZ)

  !****************************************************** MPI COLLECT DATA
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ======================================================================

end

!------------------------------------------------------------------------------
!> Output rms errors to screen.
subroutine printRMSerror(ITER, NSPIN, RMSM, RMSQ)
  implicit none
  integer :: ITER
  integer :: NSPIN
  double precision :: RMSM
  double precision :: RMSQ

  write(6,'(79(1H-),/)')
  if (NSPIN.eq.2) then
    write (6,fmt=9041) ITER,RMSQ,RMSM
  else
    write (6,fmt=9051) ITER,RMSQ
  end if
  write(6,'(79(1H-))')

9041 format ('      ITERATION',I4,' average rms-error : v+ + v- = ', &
  1p,d11.4,/,39x,' v+ - v- = ',1p,d11.4)
9051 format ('      ITERATION',I4,' average rms-error : v+ + v- = ', &
  1p,d11.4)
end subroutine


!================================================================
! Does formatted output of potential if file VFORM exists.
! and       b) (optional) potential in formatted format
!              files are named VPOT.X, X = atom number
!              >> output can be enable providing file VFORM

! called by main2
! subroutine called: rites

!------------------------------------------------------------------------------
!> Tests if file 'VFORM' exists
!> TODO: replace by better construct
logical function testVFORM()
  implicit none

  testVFORM = .false.
  inquire(file='VFORM',exist=testVFORM)

end function

!===========================================================================================
!> Output formatted potential files for each atom.
subroutine writeFormattedPotentialImpl(fermiEnergy, &
VBC,NSPIN, &
KXC,LPOT,A,B,IRC, &
VINS,VISP,DRDI,IRNS,R,RWS,RMT,ALAT, &
ECORE,LCORE,NCORE,ZAT,ITITLE, &
rank, irmd, irnsd) !  new input parameters after inc.p removal

  implicit none

  integer irmd
  integer irnsd

  double precision fermiEnergy
  double precision ALAT,VBC(2),A,B

  double precision VINS(IRMD-IRNSD:IRMD,(LPOT+1)**2,2), &
  VISP(IRMD,2), &
  DRDI(IRMD), &
  ECORE(20,2), &
  ZAT, &
  R(IRMD),RWS,RMT

  integer NSPIN,KXC,LPOT

  integer IRNS,IRC, &
  LCORE(20,NSPIN), &
  NCORE(NSPIN), &
  ITITLE(20,NSPIN)


  !     .. local scalars ..
  integer          I1,D1,D10,D100,D1000,OFF(3)
  character(12)    FNAME
  !     ..
  !     .. MPI variables ..
  !     .. L-MPI ..
  integer rank

  ! ......................................................................
  ! formatted output                       A.Thiess 09/09
  ! simplified: E.R.
  ! ......................................................................

  I1 = rank + 1

      D1 = mod(I1,10)
      D10 = int( (mod(I1,100) + 0.5)/10 )
      D100 = int( (mod(I1,1000) + 0.5)/100 )
      D1000 = int( (mod(I1,10000) + 0.5)/1000 )

      OFF(1) = iachar('1')-1
      OFF(2) = iachar('1')-1
      OFF(3) = iachar('1')-1

      if ( D10.ge.10 ) OFF(1) = iachar('7')
      if ( D100.ge.100 ) OFF(2) = iachar('7')
      if ( D1000.ge.1000 ) OFF(3) = iachar('7')
      FNAME='VPOT.' &
      //achar(D1000+OFF(3)) &
      //achar(D100+OFF(2)) &
      //achar(D10+OFF(1)) &
      //achar(D1+iachar('1')-1)

      open(11,file=FNAME,form='formatted')

      call RITES(11,NSPIN,ZAT,ALAT,RMT,RMT,RWS, &
      ITITLE,R,DRDI,VISP,A,B,KXC,IRNS,LPOT,VINS, &
      IRC,fermiEnergy,VBC,ECORE,LCORE,NCORE, &
      irmd, irnsd)

      close (11)
! ......................................................................
! ......................................................................

end
