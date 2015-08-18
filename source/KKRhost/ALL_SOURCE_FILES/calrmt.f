c ************************************************************************
      SUBROUTINE CALRMT(IPF,IPFE,IPE,IMT,Z,RMT,RWS,RMTNEW,ALAT,DRDI,A,B,
     +                  IRWS,R,IFILE,KSHAPE)
c***********************************************************************
c     this subroutine calculates imt and rmt(cal-rmt)
c                     and prints some informations about the used meshes
c        imtl = maximumnumber of meshpoints generating a radius
c               less or equal than rmt
c        imt  = number of meshpoint generating a new mt-radius closer th
c               mt-radius than every ather meshpoint
c***********************************************************************
C     .. Scalar Arguments ..
      DOUBLE PRECISION A,ALAT,B,RMT,RMTNEW,RWS,Z
      INTEGER IFILE,IMT,IPE,IPF,IPFE,IRWS,KSHAPE
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(*),R(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DRD1,DRDWS,RIMT,RIMTM1,RNUC
      INTEGER IDELTA,IH,IMTL,IRWSM2
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP,LOG,MOD,REAL
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
      IF (KSHAPE.EQ.0) THEN
        RIMT = LOG(RMT/B+1.D0)/A + 1.D0
        IMTL = RIMT
        IRWSM2 = IRWS - 2
        IDELTA = (RIMT-IMTL)*2
        IF (IDELTA.EQ.0) IMT = IMTL
        IF (IDELTA.GT.0) IMT = IMTL + 1
        RIMTM1 = REAL(IMT-1)
        RMTNEW = B*EXP(A*RIMTM1) - B
c
        IF (IMT.GT.IRWSM2) THEN
          WRITE (IPF,FMT=9000)
          CALL RCSTOP('calrmt  ')

        END IF

      ELSE

        IF (MOD(IMT,2).EQ.0) THEN
          WRITE (IPF,FMT=*) ' error stop in calrmt - imt = ',IMT,
     +      ' has to be odd to get proper core charge  '
          CALL RCSTOP('29      ')

        END IF

      END IF
c
      IH = IRWS/2
      DRD1 = DRDI(1)
      DRDWS = DRDI(IRWS)
c----- nucleus radius rnuc in bohr's radii
      RNUC = 2.2677022D-5* (2.D0*Z)** (1.0D0/3.0D0)
c-----
      IF (IFILE.NE.0) THEN
        WRITE (IPF,FMT=9010) Z,A,B,RNUC,R(2),IH,R(IH),DRD1,DRDWS
        WRITE (IPF,FMT=9020) IRWS,IMT,RWS,RMT,RMTNEW,ALAT
        IF (IPE.EQ.1) WRITE (IPFE,FMT=9010) Z,A,B,RNUC,R(2),IH,R(IH),
     +      DRD1,DRDWS
        IF (IPE.EQ.1) WRITE (IPFE,FMT=9020) IRWS,IMT,RWS,RMT,RMTNEW,
     +      ALAT
      END IF



 9000 FORMAT (1x,'potentials need more meshpoints',/,50 ('*'))
 9010 FORMAT (' rmesh  z=',f5.2,'  a=',f6.4,'  b=',f8.6,'  rnuc=',f10.8,
     +       '  r(2)=',f10.8,/,' r(',i3,')=',f6.4,'   drdi(1)=',f10.8,
     +       '   drdi(irws)=',f8.6)
 9020 FORMAT (' irws=',i6,' imt=',i6,/,' rws=',f12.8,' rmt=',f12.8,
     +       ' rmtnew=',f12.8,' alat=',f12.8)
      END
