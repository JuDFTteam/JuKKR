! ************************************************************************
SUBROUTINE calrmt(ipf,ipfe,ipe,imt,z,rmt,rws,rmtnew,alat,drdi,a,b,  &
    irws,r,ifile,kshape)
!***********************************************************************
!     this subroutine calculates imt and rmt(cal-rmt)
!                     and prints some informations about the used meshes
!        imtl = maximumnumber of meshpoints generating a radius
!               less or equal than rmt
!        imt  = number of meshpoint generating a new mt-radius closer th
!               mt-radius than every ather meshpoint
!***********************************************************************
!.. Scalar Arguments ..
      DOUBLE PRECISION A,ALAT,B,RMT,RMTNEW,RWS,Z
      INTEGER IFILE,IMT,IPE,IPF,IPFE,IRWS,KSHAPE
!..
!.. Array Arguments ..
      DOUBLE PRECISION DRDI(*),R(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION DRD1,DRDWS,RIMT,RIMTM1,RNUC
      INTEGER IDELTA,IH,IMTL,IRWSM2
!..
!.. Intrinsic Functions ..
      INTRINSIC EXP,LOG,MOD,REAL
!..
!.. External Subroutines ..
      EXTERNAL RCSTOP
!     ..
IF (kshape == 0) THEN
  rimt = LOG(rmt/b+1.d0)/a + 1.d0
  imtl = rimt
  irwsm2 = irws - 2
  idelta = (rimt-imtl)*2
  IF (idelta == 0) imt = imtl
  IF (idelta > 0) imt = imtl + 1
  rimtm1 = REAL(imt-1)
  rmtnew = b*EXP(a*rimtm1) - b
  
  IF (imt > irwsm2) THEN
    WRITE (ipf,FMT=9000)
    CALL rcstop('calrmt  ')
    
  END IF
  
ELSE
  
  IF (MOD(imt,2) == 0) THEN
    WRITE (ipf,FMT=*) ' error stop in calrmt - imt = ',imt,  &
        ' has to be odd to get proper core charge  '
    CALL rcstop('29      ')
    
  END IF
  
END IF

ih = irws/2
drd1 = drdi(1)
drdws = drdi(irws)
!----- nucleus radius rnuc in bohr's radii
rnuc = 2.2677022D-5* (2.d0*z)** (1.0D0/3.0D0)
!-----
IF (ifile /= 0) THEN
  WRITE (ipf,FMT=9010) z,a,b,rnuc,r(2),ih,r(ih),drd1,drdws
  WRITE (ipf,FMT=9020) irws,imt,rws,rmt,rmtnew,alat
  IF (ipe == 1) WRITE (ipfe,FMT=9010) z,a,b,rnuc,r(2),ih,r(ih), drd1,drdws
  IF (ipe == 1) WRITE (ipfe,FMT=9020) irws,imt,rws,rmt,rmtnew, alat
END IF



9000 FORMAT (1X,'potentials need more meshpoints',/,50 ('*'))
9010 FORMAT (' rmesh  z=',f5.2,'  a=',f7.4,'  b=',f9.6,'  rnuc=',f11.8,  &
    '  r(2)=',f11.8,/,' r(',i3,')=',f7.4,'   drdi(1)=',f11.8,  &
    '   drdi(irws)=',f9.6)
9020 FORMAT (' irws=',i6,' imt=',i6,/,' rws=',f12.8,' rmt=',f12.8,  &
    ' rmtnew=',f12.8,' alat=',f12.8)
END SUBROUTINE calrmt
