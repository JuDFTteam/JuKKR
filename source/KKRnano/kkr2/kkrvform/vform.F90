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
  integer          I1,D1,D10,D100,D1000,D10000,OFF(4)
  character(12)    FNAME
  !     ..
  !     .. MPI variables ..
  !     .. L-MPI ..
  integer rank

  ! ......................................................................
  ! formatted output                       A.Thiess 09/09
  ! simplified: E.R.
  ! ......................................................................

  I1 = rank

      D1 = mod(I1,10)
      D10 = int( (mod(I1,100) + 0.5)/10 )
      D100 = int( (mod(I1,1000) + 0.5)/100 )
      D1000 = int( (mod(I1,10000) + 0.5)/1000 )
      D10000 = int( (mod(I1,100000) + 0.5)/10000 )

      OFF(1) = iachar('1')-1
      OFF(2) = iachar('1')-1
      OFF(3) = iachar('1')-1
      OFF(4) = iachar('1')-1

      if ( D10.ge.10 ) OFF(1) = iachar('7')
      if ( D100.ge.100 ) OFF(2) = iachar('7')
      if ( D1000.ge.1000 ) OFF(3) = iachar('7')
      if ( D10000.ge.10000 ) OFF(4) = iachar('7')

      FNAME='VPOT.' &
      //achar(D10000+OFF(4)) &
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
