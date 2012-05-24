
! process with MYLRANK(LMPIC).EQ.0 and LMPIC.EQ.1 writes results

subroutine RESULTS(LRECRES2,IELAST,ITSCF,LMAX,NAEZ,NPOL,NSPIN, &
  KPRE,KTE,LPOT,E1,E2,TK,EFERMI,ALAT,ITITLE,CHRGNT,ZAT,EZ,WEZ, &
  LDAU, &
  iemxd) ! new input parameter after inc.p removal

  implicit none

  integer iemxd

  !     .. Parameters ..
  !     ..
  !     .. Local Scalars ..
  integer IELAST,ITSCF,LMAX,NAEZ,NPOL,NSPIN
  integer KPRE,KTE
  integer I1,ISPIN,LPOT
  integer LRECRES1,LRECRES2
  !     ..
  double precision E1,E2,TK,EFERMI
  double precision CHRGNT,TOTSMOM,ALAT,PI
  !     ..
  logical TEST,LDAU
  !     ..
  !     .. Local Arrays ..
  double complex EZ(IEMXD),WEZ(IEMXD)
  double complex DEN(0:LMAX+1,IEMXD,NSPIN)
  double precision DOSTOT(0:LMAX+1,2)
  double precision ECOU(0:LPOT),EPOTIN,EULDAU,EDCLDAU, &
  ESPC(0:3,NSPIN),ESPV(0:LMAX+1,NSPIN), &
  EXC(0:LPOT)
  double precision ECORE(20,2)
  double precision ZAT(NAEZ)
  double precision CHARGE(0:LMAX+1,2)
  !     ..
  double precision CATOM(NSPIN),QC
  double precision VMAD

  !     INTEGER ITITLE(20,NPOTD)
  integer ITITLE(20,NSPIN*NAEZ)
  integer LCOREMAX

  integer NPOTD
  NPOTD = NSPIN*NAEZ

  PI = 4.0D0*ATAN(1.0D0)

  LRECRES1 = 8*43 + 16*(LMAX+2)


  if (NPOL.eq.0 .or. TEST('DOS     ')) then
    LRECRES1 = LRECRES1 + 32*(LMAX+2)*IEMXD
  end if


  open (71,access='direct',recl=LRECRES1,file='results1', &
  form='unformatted')


  do I1 = 1,NAEZ
    if (NPOL.eq.0 .or. TEST('DOS     ')) then
      read(71,rec=I1) QC,CATOM,CHARGE,ECORE,DEN
    else
      read(71,rec=I1) QC,CATOM,CHARGE,ECORE
    end if
    call WRMOMS(NAEZ,NSPIN,CHARGE,I1,LMAX,LMAX+1)
  end do


  if (NPOL.eq.0 .or. TEST('DOS     ')) then 
    do I1 = 1,NAEZ
      read(71,rec=I1) QC,CATOM,CHARGE,ECORE,DEN
      call WRLDOS(DEN,EZ,WEZ, &
      LMAX+1,IEMXD,NPOTD,ITITLE,EFERMI,E1,E2,ALAT,TK, &
      NSPIN,NAEZ,IELAST,I1,DOSTOT)
    end do
  end if


  TOTSMOM = 0.0D0
  do I1 = 1,NAEZ
    if (NPOL.eq.0 .or. TEST('DOS     ')) then
      read(71,rec=I1) QC,CATOM,CHARGE,ECORE,DEN
    else
      read(71,rec=I1) QC,CATOM,CHARGE,ECORE
    end if
    do ISPIN = 1,NSPIN
      if (ISPIN.ne.1) then
        write (6,fmt=9011) CATOM(ISPIN)
      else
        write (6,fmt=9001) I1,CATOM(ISPIN)
      end if
    end do
    write (6,fmt=9041) ZAT(I1),QC
    if (NSPIN.eq.2) TOTSMOM = TOTSMOM + CATOM(NSPIN)
  end do
  write(6,'(79(1H+))')
  write (6,fmt=9021) ITSCF,CHRGNT
  if (NSPIN.eq.2) write (6,fmt=9031) TOTSMOM
  write(6,'(79(1H+))')


9001 format ('  Atom ',I4,' charge in wigner seitz cell =',f10.6)
9011 format (7X,'spin moment in wigner seitz cell =',f10.6)
9021 format ('      ITERATION',I4, &
  ' charge neutrality in unit cell = ',f12.6)
9031 format ('                   ', &
  ' TOTAL mag. moment in unit cell = ',f12.6)
9041 format (4X,'nuclear charge  ',F10.6,9X,'core charge =   ',F10.6)
  !        WRITE(6,FMT=99001)
  !        WRITE(6,FMT=99002)
  !99001 FORMAT (79(1H=),/,18X,' MADELUNG POTENTIALS ',
  !     &        '(spherically averaged) ')
  !99002 FORMAT (/,25X,' ATOM ','  Delta_Q  ','     VMAD',/,25X,30(1H-))
  !99003 FORMAT (25X,I4,2X,F10.6,1X,F12.6)
  close(71)

  !=======================================================================
  ! output of information stored in 'results2'
  ! set KTE=1 in inputcard for output of energy contributions
  !=======================================================================

  open (72,access='direct',recl=LRECRES2,file='results2', &
  form='unformatted')

  do I1 = 1,NAEZ
    read(72,rec=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX, &
    EULDAU,EDCLDAU
  !        WRITE (6,FMT=99003) I1,(CATOM(1)-ZAT(I1)),VMAD
  end do
  !      WRITE(6,'(25X,30(1H-),/)')
  !      WRITE(6,'(79(1H=))')

  if (KTE.eq.1) then
    do I1 = 1,NAEZ
      read(72,rec=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX, &
      EULDAU,EDCLDAU
      call ETOTB1(ECOU,EPOTIN,ESPC,ESPV,EXC, &
      EULDAU,EDCLDAU,LDAU, &
      KPRE,LMAX,LPOT, &
      LCOREMAX,NSPIN,I1,NAEZ)
    end do
  end if

  close(72)

  !=======================================================================

  return

end
