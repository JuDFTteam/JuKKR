C*==scfchrdns.f    processed by SPAG 6.05Rc at 14:43 on  7 Apr 2001
      SUBROUTINE SCFCHRDNS(NFILCBWF,R2DRDI,JWS,IMT,SHFTEF,TOTDOS,MUESPN,
     &                     MUEORB,IREL,IPRINT,NT,NL,NKM,ERYD,WE,EFERMI,
     &                     IECURR,NETAB,DOS,SMT,OMT,HFF,DOSI,SMTI,OMTI,
     &                     HFFI,DOSM,DOSL0,DOSINT,SMTM,SMTL0,SMTINT,
     &                     OMTM,OMTL0,OMTINT,HFFM,HFFL0,HFFINT,BCOR,
     &                     BCORS,DZZ,DZJ,SZZ,SZJ,OZZ,OZJ,BZZ,BZJ,
     &                     OZZS,OZJS,OMTLS0,TAUTLIN,NVALTOT,TXTT,
     &                     CONC,NAT,RHOCHR,RHOSPN,RHOORB,QEL,GDIA,GMDIA,
     &                     GOFF,NTMAX,NLMAX,NMUEMAX,LINMAX,NRMAX,NMMAX,
     &                     NKMMAX,EBAND,EBANDT)
C   ********************************************************************
C   *                                                                  *
C   * SUBROUTINE TO CALCULATE THE  CHARGE, SPIN  AND  ORBITAL DENSITY  *
C   *                  WITHIN AN ATOMIC CELL                           *
C   *                                                                  *
C   * 12/03/96 HE                                                      *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      INTEGER NTMAXCHK
      PARAMETER (NTMAXCHK=10)
      COMPLEX*16 C0
      PARAMETER (C0=(0.0D0,0.0D0))
      REAL*8 PI
      PARAMETER (PI=3.141592653589793238462643D0)
C
C Dummy arguments
C
      COMPLEX*16 EBAND,ERYD,WE
      REAL*8 EFERMI,MUEORB,MUESPN,SHFTEF,TOTDOS,NVALTOT
      INTEGER IECURR,IPRINT,IREL,LINMAX,NETAB,NFILCBWF,NKM,NKMMAX,NL,
     &        NLMAX,NMMAX,NMUEMAX,NRMAX,NT,NTMAX
      REAL*8 BCOR(NTMAX),BCORS(NTMAX),CONC(NTMAX),DOS(NTMAX),DOSI(NTMAX)
     &       ,GDIA(NKMMAX),GMDIA(NKMMAX),GOFF(NKMMAX),HFF(NTMAX),
     &       HFFI(NTMAX),OMT(NTMAX),OMTI(NTMAX),QEL(NTMAX),
     &       R2DRDI(NRMAX,NMMAX),RHOCHR(NRMAX,NTMAX),RHOORB(NRMAX,NTMAX)
     &       ,RHOSPN(NRMAX,NTMAX),SMT(NTMAX),SMTI(NTMAX)
      COMPLEX*16 BZJ(LINMAX,NTMAX),BZZ(LINMAX,NTMAX),DOSINT(NLMAX,NTMAX)
     &           ,DOSL0(NLMAX,NTMAX),DOSM(NMUEMAX),DZJ(LINMAX,NTMAX),
     &           DZZ(LINMAX,NTMAX),EBANDT(NTMAX),HFFINT(NLMAX,NTMAX),
     &           HFFL0(NLMAX,NTMAX),HFFM(NMUEMAX),OMTINT(NLMAX,NTMAX),
     &           OMTL0(NLMAX,NTMAX),OMTM(NMUEMAX),OZJ(LINMAX,NTMAX),
     &           OZZ(LINMAX,NTMAX),SMTINT(NLMAX,NTMAX),
     &           SMTL0(NLMAX,NTMAX),SMTM(NMUEMAX),SZJ(LINMAX,NTMAX),
     &           SZZ(LINMAX,NTMAX),TAUTLIN(LINMAX,NTMAX)
      COMPLEX*16 OZZS(LINMAX,NTMAX,2),OZJS(LINMAX,NTMAX,2),
     &           OMTLS0(NLMAX,NTMAX,2)
      INTEGER IMT(NTMAX),JWS(NMMAX),NAT(NTMAX)
      CHARACTER*4 TXTT(NTMAX)
C
C Local variables
C
      REAL*8 AUX,BDUM(3),CFF(2,2),CFG(2,2),CGF(2,2),CGG(2,2),
     &       CHKO(NTMAXCHK),CHKQ(NTMAXCHK),CHKS(NTMAXCHK),DEFERMI,DQ,MJ,
     &       MJMAX,MJMIN,R1M(2,2),RINT(NRMAX),TOTNOS
      DOUBLE PRECISION DBLE,DSQRT
      COMPLEX*16 DOSL,HFFL,JF(NRMAX,2,2),JG(NRMAX,2,2),OMTL,SMTL,WDS,
     &           WOF,WOG,WSF,WSG,WT,ZF(NRMAX,2,2),ZFJF,ZFZF,
     &           ZG(NRMAX,2,2),ZGJG,ZGZG
      COMPLEX*16 OMTLS(2),OMTMS(NMUEMAX,2)
      INTEGER I,IFLAG,IKM1,IKM2,IL,IM,IS,IT,JJ,JTOP,K1,K2,KA,
     &        KAP1,KAP2,KB,L,LIN,LMAX,MM,MUE,NSOL
      INTEGER IKAPMUE
      INTEGER NINT
C
      SAVE CHKO,CHKQ,CHKS
C
C*** End of declarations rewritten by SPAG
C
      DATA R1M/1.0D0,0.0D0,0.0D0,1.0D0/
C
      IF ( IECURR.EQ.1 ) THEN
C
         DO IT = 1,NT
            DO IL = 1,NL
               DOSINT(IL,IT) = C0
               SMTINT(IL,IT) = C0
               OMTINT(IL,IT) = C0
               HFFINT(IL,IT) = C0
            END DO
            EBANDT(IT) = C0
         END DO
C
C ----------------------------- account for spin degeneracy for IREL <=1
         IF ( IREL.LE.1 ) THEN
            DO IT = 1,NT
               IM = IMT(IT)
               JTOP = JWS(IM)
               DO I = 1,JTOP
                  RHOCHR(I,IT) = RHOCHR(I,IT)/2.0D0
                  RHOSPN(I,IT) = 0.0D0
                  RHOORB(I,IT) = 0.0D0
               END DO
            END DO
         END IF
C
         IF ( IPRINT.GT.0 ) THEN
            IF ( NT.GT.NTMAXCHK ) STOP '<SCFCHRDNS> NT > NTMAXCHK'
            DO IT = 1,NT
               IM = IMT(IT)
               JTOP = JWS(IM)
               DO I = 1,JTOP
                  RINT(I) = RHOCHR(I,IT)*R2DRDI(I,IM)
               END DO
               CALL RINTSIMP(RINT,JTOP,CHKQ(IT))
               DO I = 1,JTOP
                  RINT(I) = RHOSPN(I,IT)*R2DRDI(I,IM)
               END DO
               CALL RINTSIMP(RINT,JTOP,CHKS(IT))
               DO I = 1,JTOP
                  RINT(I) = RHOORB(I,IT)*R2DRDI(I,IM)
               END DO
               CALL RINTSIMP(RINT,JTOP,CHKO(IT))
            END DO
         END IF
      END IF
C
 100  CONTINUE
      TOTNOS = 0.0D0
      TOTDOS = 0.0D0
      MUESPN = 0.0D0
      MUEORB = 0.0D0
C
C
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
      DO IT = 1,NT
         IM = IMT(IT)
         JTOP = JWS(IM)
C
         LMAX = NL - 1
         LIN = 0
C
         DOS(IT) = 0.0D0
         SMT(IT) = 0.0D0
         OMT(IT) = 0.0D0
         HFF(IT) = 0.0D0
         DOSI(IT) = 0.0D0
         SMTI(IT) = 0.0D0
         OMTI(IT) = 0.0D0
         HFFI(IT) = 0.0D0
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
         DO L = 0,LMAX
            IL = L + 1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            KAP1 = -L - 1
            KAP2 = L
            IF ( L.EQ.0 ) KAP2 = KAP1
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
            DOSL = 0.0D0
            SMTL = 0.0D0
            OMTL = 0.0D0
            HFFL = 0.0D0
            OMTLS(1) = OMTL
            OMTLS(2) = OMTL
C
            IF ( IREL.GT.1 ) THEN
               MJMAX = DBLE(L) + 0.5D0
            ELSE
               MJMAX = DBLE(L)
            END IF
            MJMIN = -MJMAX
            MUE = 0
C
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            DO MJ = MJMIN,MJMAX,1.0D0
               MUE = MUE + 1
               DOSM(MUE) = 0.0D0
               SMTM(MUE) = 0.0D0
               OMTM(MUE) = 0.0D0
               HFFM(MUE) = 0.0D0
C 
               DO IS = 1,2
                  OMTMS(MUE,IS) = 0.0D0
               END DO
C
               IF ( IREL.LE.1 ) THEN
                  NSOL = 1
C           no coupling for:  abs(mue)= j   +  j=l+1/2 == kap=-l-1
               ELSE IF ( ABS(MJ).GT.DBLE(L) ) THEN
                  NSOL = 1
               ELSE
                  NSOL = 2
               END IF
C------------------------------------------------------------------------
               IKM1 = IKAPMUE(KAP1,NINT(MJ-0.5D0))
               IKM2 = IKAPMUE(KAP2,NINT(MJ-0.5D0))
C------------------------------------------------------------------------
               IF ( IREL.LE.1 ) THEN
                  IKM1 = IL
                  IKM2 = IL
                  IF ( NKM.NE.NL**2 ) WRITE (1337,99001) NKM
               END IF
C
C
C   COEFFICIENTS TO CALCULATE THE SPIN  MAGNETISATION
C
               CGG(1,1) = GDIA(IKM1)
               CGG(1,2) = GOFF(IKM1)
               CGG(2,1) = GOFF(IKM1)
               CGG(2,2) = GDIA(IKM2)
               CALL RINIT(4,CGF)
               CGF(1,1) = GMDIA(IKM1)
               CGF(2,2) = GMDIA(IKM2)
C
C   COEFFICIENTS TO CALCULATE THE ORBITAL MAGNETISATION
C
               CFG(1,1) = MJ*(KAP1+1.0D0)/(KAP1+0.5D0)
               CFG(2,2) = MJ*(KAP2+1.0D0)/(KAP2+0.5D0)
               CFG(1,2) = 0.5D0*DSQRT(1.0D0-(MJ/(KAP1+0.5D0))**2)
               CFG(2,1) = CFG(1,2)
               CALL RINIT(4,CFF)
               CFF(1,1) = MJ*(-KAP1+1.0D0)/(-KAP1+0.5D0)
               CFF(2,2) = MJ*(-KAP2+1.0D0)/(-KAP2+0.5D0)
C
C -------------------------------------------------- read wave functions
C
               CALL READWFUN(NFILCBWF,IT,L,MJ,NSOL,'REG','IRR',
     &                       IKM1,KAP1,IKM2,KAP2,NT,NKM,ZG,ZF,JG,JF,
     &                       JTOP,NRMAX)
C
               DO K1 = 1,NSOL
                  DO K2 = 1,NSOL
                     LIN = LIN + 1
                     WT = -TAUTLIN(LIN,IT)/PI
                     DOSM(MUE) = DOSM(MUE) + WT*DZZ(LIN,IT)
                     SMTM(MUE) = SMTM(MUE) + WT*SZZ(LIN,IT)
                     OMTM(MUE) = OMTM(MUE) + WT*OZZ(LIN,IT)
                     HFFM(MUE) = HFFM(MUE) + WT*BZZ(LIN,IT)
C
                     DO IS = 1,2
                        OMTMS(MUE,IS) = OMTMS(MUE,IS) 
     &                                   + WT*OZZS(LIN,IT,IS)
                     END DO
C
                     DO KA = 1,NSOL
                        DO KB = 1,NSOL
                           WDS = WE*WT*R1M(KA,KB)
                           WSG = WE*WT*CGG(KA,KB)
                           WSF = WE*WT*CGF(KA,KB)
                           WOG = WE*WT*CFG(KA,KB)
                           WOF = WE*WT*CFF(KA,KB)
                           DO I = 1,JTOP
                              ZGZG = ZG(I,KA,K1)*ZG(I,KB,K2)
                              ZFZF = ZF(I,KA,K1)*ZF(I,KB,K2)
                              RHOCHR(I,IT) = RHOCHR(I,IT)
     &                           + DIMAG(WDS*ZGZG+WDS*ZFZF)
                              RHOSPN(I,IT) = RHOSPN(I,IT)
     &                           + DIMAG(WSG*ZGZG-WSF*ZFZF)
                              RHOORB(I,IT) = RHOORB(I,IT)
     &                           + DIMAG(WOG*ZGZG-WOF*ZFZF)
                           END DO
                        END DO
                     END DO
C
C    NO IRREGULAR CONTRIBUTIONS TO THE BACKSCATTERING TERMS
                     IF ( K1.EQ.K2 ) THEN
                        DOSM(MUE) = DOSM(MUE) + DZJ(LIN,IT)/PI
                        SMTM(MUE) = SMTM(MUE) + SZJ(LIN,IT)/PI
                        OMTM(MUE) = OMTM(MUE) + OZJ(LIN,IT)/PI
                        HFFM(MUE) = HFFM(MUE) + BZJ(LIN,IT)/PI
C
                        DO IS = 1,2
                           OMTMS(MUE,IS) = OMTMS(MUE,IS) 
     &                                     + OZJS(LIN,IT,IS)/PI
                        END DO
C
                        DO KA = 1,NSOL
                           DO KB = 1,NSOL
                              WDS = WE*R1M(KA,KB)/PI
                              WSG = WE*CGG(KA,KB)/PI
                              WSF = WE*CGF(KA,KB)/PI
                              WOG = WE*CFG(KA,KB)/PI
                              WOF = WE*CFF(KA,KB)/PI
                              DO I = 1,JTOP
                                 ZGJG = ZG(I,KA,K1)*JG(I,KB,K2)
                                 ZFJF = ZF(I,KA,K1)*JF(I,KB,K2)
                                 RHOCHR(I,IT) = RHOCHR(I,IT)
     &                              + DIMAG(WDS*ZGJG+WDS*ZFJF)
                                 RHOSPN(I,IT) = RHOSPN(I,IT)
     &                              + DIMAG(WSG*ZGJG-WSF*ZFJF)
                                 RHOORB(I,IT) = RHOORB(I,IT)
     &                              + DIMAG(WOG*ZGJG-WOF*ZFJF)
                              END DO
                           END DO
                        END DO
                     END IF
                  END DO
               END DO
C
C
               DOSL = DOSL + DOSM(MUE)
               SMTL = SMTL + SMTM(MUE)
               OMTL = OMTL + OMTM(MUE)
               HFFL = HFFL + HFFM(MUE)
C
               DO IS = 1,2
                  OMTLS(IS) = OMTLS(IS) + OMTMS(MUE,IS)
               END DO
C
            END DO
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
            EBANDT(IT) = EBANDT(IT) + WE*DOSL*ERYD
C
            DOSINT(IL,IT) = DOSINT(IL,IT) + WE*DOSL
            SMTINT(IL,IT) = SMTINT(IL,IT) + WE*SMTL
            OMTINT(IL,IT) = OMTINT(IL,IT) + WE*OMTL
            HFFINT(IL,IT) = HFFINT(IL,IT) + WE*HFFL
C
            DOSL0(IL,IT) = DOSL
            SMTL0(IL,IT) = SMTL
            OMTL0(IL,IT) = OMTL
            HFFL0(IL,IT) = HFFL
C
            DO IS = 1,2
               OMTLS0(IL,IT,IS) = OMTLS(IS)
            END DO
C ----------------------------------------------------------------------
            IF ( (IPRINT.GT.0) .AND. (IECURR.EQ.NETAB) ) THEN
C
               IF ( IREL.GT.1 ) THEN
                  JJ = 2*L + 2
C
                  WRITE (1337,99005) IECURR,ERYD,L,IT,TXTT(IT),
     &                            'CRYSTAL TERMS       ',DOSINT(IL,IT),
     &                            SMTINT(IL,IT),OMTINT(IL,IT),
     &                            (HFFINT(IL,IT)*1D-3),DOSL,SMTL,OMTL,
     &                            (HFFL*1D-6),
     &                            ((-JJ-1+2*MUE),DOSM(MUE),SMTM(MUE),
     &                            OMTM(MUE),(HFFM(MUE)*1D-6),MUE=1,JJ)
C
               ELSE
                  JJ = 2*L + 1
C
                  WRITE (1337,99014) IECURR,ERYD,L,IT,TXTT(IT),
     &                            'CRYSTAL TERMS       ',DOSINT(IL,IT),
     &                            SMTINT(IL,IT),OMTINT(IL,IT),
     &                            (HFFINT(IL,IT)*1D-3),DOSL,SMTL,OMTL,
     &                            (HFFL*1D-6),
     &                            ((-L-1+MM),DOSM(MM),SMTM(MM),OMTM(MM),
     &                            (HFFM(MM)*1D-6),MM=1,JJ)
C
               END IF
            END IF
C ---------------------------------------------------------------------
C
C
            DOS(IT) = DOS(IT) + DIMAG(DOSL)
            SMT(IT) = SMT(IT) + DIMAG(SMTL)
            OMT(IT) = OMT(IT) + DIMAG(OMTL)
            HFF(IT) = HFF(IT) + DIMAG(HFFL)
            DOSI(IT) = DOSI(IT) + DIMAG(DOSINT(IL,IT))
            SMTI(IT) = SMTI(IT) + DIMAG(SMTINT(IL,IT))
            OMTI(IT) = OMTI(IT) + DIMAG(OMTINT(IL,IT))
            HFFI(IT) = HFFI(IT) + DIMAG(HFFINT(IL,IT))
         END DO
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C
C
         TOTNOS = TOTNOS + DOSI(IT)*CONC(IT)*NAT(IT)
         TOTDOS = TOTDOS + DOS(IT)*CONC(IT)*NAT(IT)
         MUESPN = MUESPN + SMTI(IT)*CONC(IT)*NAT(IT)
         MUEORB = MUEORB + OMTI(IT)*CONC(IT)*NAT(IT)
C
         IF ( IPRINT.GT.0 ) THEN
C
            WRITE (1337,99006) IECURR,ERYD,IT,TXTT(IT),DOSI(IT),
     &                      SMTI(IT),OMTI(IT),(HFFI(IT)*1D-3),DOS(IT),
     &                      SMT(IT),OMT(IT),(HFF(IT)*1D-3)
C
            IF ( IT.LT.NT ) THEN
               WRITE (1337,'(1X,79(''-''))')
            ELSE IF ( (IPRINT.GT.0) .OR. (IECURR.EQ.NETAB) ) THEN
               WRITE (1337,99010) TOTDOS,TOTNOS,MUESPN,MUEORB
            ELSE
               WRITE (1337,'('' '',79(''=''))')
            END IF
         END IF
C
      END DO
C TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
C
      IF ( IECURR.EQ.NETAB ) THEN
C
         EBAND = C0
         DO IT=1,NT
            EBAND = EBAND + EBANDT(IT)*CONC(IT)*NAT(IT)
         END DO

         IF ( IREL.GT.1 ) THEN
            DQ = NVALTOT - TOTNOS
         ELSE
            DQ = NVALTOT/2.0D0 - TOTNOS
         END IF
C
         DEFERMI = DQ/TOTDOS
C
         IF ( ABS(DQ).GT.1D-06 ) THEN
            WE = DEFERMI
            EFERMI = EFERMI + DEFERMI
            SHFTEF = DEFERMI
C
            WRITE (1337,'(/)')
            WRITE (1337,99012) (TXTT(IT),CONC(IT),IT=1,NT)
            WRITE (1337,99013) DQ,DEFERMI,EFERMI
C
            GOTO 100
         END IF
C
         IF ( IPRINT.GT.0 ) THEN
            IFLAG = 0
            DO IT = 1,NT
               IM = IMT(IT)
               JTOP = JWS(IM)
               DO I = 1,JTOP
                  RINT(I) = RHOCHR(I,IT)*R2DRDI(I,IM)
               END DO
               CALL RINTSIMP(RINT,JTOP,AUX)
               CHKQ(IT) = AUX - CHKQ(IT)
               IF ( ABS(CHKQ(IT)-DOSI(IT)).GT.1.0D-8 ) THEN
                  IFLAG = 1
                  WRITE (1337,99004) IT,'Q',DOSI(IT),CHKQ(IT),DOSI(IT)
     &                            /CHKQ(IT)
               END IF
               DO I = 1,JTOP
                  RINT(I) = RHOSPN(I,IT)*R2DRDI(I,IM)
               END DO
               CALL RINTSIMP(RINT,JTOP,AUX)
               CHKS(IT) = AUX - CHKS(IT)
               IF ( ABS(CHKS(IT)-SMTI(IT)).GT.1.0D-8 ) THEN
                  IFLAG = 1
                  WRITE (1337,99004) IT,'S',SMTI(IT),CHKS(IT),SMTI(IT)
     &                            /CHKS(IT)
               END IF
               DO I = 1,JTOP
                  RINT(I) = RHOORB(I,IT)*R2DRDI(I,IM)
               END DO
               CALL RINTSIMP(RINT,JTOP,AUX)
               CHKO(IT) = AUX - CHKO(IT)
               IF ( ABS(CHKO(IT)-OMTI(IT)).GT.1.0D-8 ) THEN
                  IFLAG = 1
                  WRITE (1337,99004) IT,'O',OMTI(IT),CHKO(IT),OMTI(IT)
     &                            /CHKO(IT)
               END IF
            END DO
C
            IF ( IFLAG.EQ.0 ) THEN
               WRITE (1337,99002)
            ELSE
               WRITE (1337,99003)
            END IF
         END IF
C
C
         DO IT = 1,NT
C
            WRITE (1337,99006) (IECURR+1),EFERMI,0.0D0,IT,TXTT(IT)
C
            BDUM(1) = BCORS(IT)*1D-3
            BDUM(2) = (BCOR(IT)-BCORS(IT))*1D-3
            BDUM(3) = BCOR(IT)*1D-3
C
            WRITE (1337,99007) (DIMAG(DOSL0(IL,IT)),DIMAG(DOSINT(IL,IT))
     &                     ,DIMAG(SMTL0(IL,IT)),DIMAG(SMTINT(IL,IT)),
     &                      DIMAG(OMTL0(IL,IT)),DIMAG(OMTINT(IL,IT)),
     &                      DIMAG(HFFINT(IL,IT))*1D-3,BDUM(IL),IL=1,
     &                      MIN(3,NL))
            IF ( NL.GT.3 ) WRITE (1337,99008)
     &                            (DIMAG(DOSL0(IL,IT)),DIMAG(DOSINT(IL,
     &                            IT)),DIMAG(SMTL0(IL,IT)),
     &                            DIMAG(SMTINT(IL,IT)),
     &                            DIMAG(OMTL0(IL,IT)),
     &                            DIMAG(OMTINT(IL,IT)),
     &                            DIMAG(HFFINT(IL,IT))*1D-3,IL=4,NL)
C
            WRITE (1337,99009) DOS(IT),DOSI(IT),SMT(IT),SMTI(IT),OMT(IT)
     &                     ,OMTI(IT),(HFFI(IT)*1D-3),
     &                      ((HFFI(IT)+BCOR(IT))*1D-3)
C
            IF ( IT.LT.NT ) THEN
               WRITE (1337,'(1X,79(''-''))')
            ELSE
               WRITE (1337,99011) TOTDOS,TOTNOS,MUESPN,MUEORB,
     &                            DIMAG(EBAND)
            END IF
C
            IM = IMT(IT)
            JTOP = JWS(IM)
C
            DO I = 1,JTOP
               RINT(I) = RHOCHR(I,IT)*R2DRDI(I,IM)
            END DO
C
            CALL RINTSIMP(RINT,JTOP,QEL(IT))
C
C ----------------------------- account for spin degeneracy for IREL <=1
            IF ( IREL.LE.1 ) THEN
               QEL(IT) = QEL(IT)*2.0D0
               DO I = 1,JTOP
                  RHOCHR(I,IT) = RHOCHR(I,IT)*2.0D0
                  RHOSPN(I,IT) = 0.0D0
                  RHOORB(I,IT) = 0.0D0
               END DO
            END IF
C
         END DO
C
      END IF
C
99001 FORMAT ('warning in <SCFCHRDNS>:  IREL<=1 and  NL**2 <> NKM=',I5)
99002 FORMAT (/,10X,'integrals in <SCFCHRDNS> agree within 1D-8',/)
99003 FORMAT (/,10X,'... integrals in <SCFCHRDNS>  NOT OK')
99004 FORMAT (' IT ',I3,2X,A,2X,F20.10,/,12X,4F20.10)
99005 FORMAT (/,I4,' E=',2F7.4,3X,'L=',I2,3X,'IT=',I2,2X,A,2X,A20,/,
     &        15X,'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',
     &        '   B_tot   [kG]',/,' INT(DE)  ',2F8.3,2X,2F8.3,2X,2F8.3,
     &        F10.1,F8.1,/,' SUM(MJ)  ',2F8.3,2X,2F8.3,2X,2F8.3,F10.1,
     &        F8.1,20(:,/,' MJ= ',I2,'/2 ',2F8.3,2X,2F8.3,2X,2F8.3,
     &        F10.1,F8.1))
99006 FORMAT (/,I4,' E=',2F7.4,10X,'IT=',I2,2X,A,:,/,15X,
     &        'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',
     &        '   B_tot   [kG]',/,' INT(DE) crystal  ',F8.3,10X,F8.3,
     &        10X,F8.3,10X,F8.1,/,' TOTAL   crystal  ',F8.3,10X,F8.3,
     &        10X,F8.3,10X,F8.1)
99007 FORMAT ('         DOS      NOS     P_spin   m_spin',
     &        '    P_orb    m_orb    B_val      B_core',/,'  s ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' s  ',F8.2,:,/,'  p ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' ns ',F8.2,:,/,'  d ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2,' cor',F8.2)
99008 FORMAT ('  f ',2F9.4,F10.4,F9.4,F10.5,F9.5,F8.2,:,/,'  g ',2F9.4,
     &        F10.4,F9.4,F10.5,F9.5,F8.2)
99009 FORMAT (' sum',2F9.4,F10.4,F9.4,F10.5,F9.5,F8.2,' v+c',F8.2)
99010 FORMAT (' ',79('-'),/,' TDOS/NOS ',2F8.3,' MUE-SPIN:',F8.3,
     &        '  MUE-ORB:',F8.3)
99011 FORMAT (' ',79('-'),/,' TOT',2F9.4,10X,F9.4,10X,F9.5,/,' E_band',
     &        F17.6,' [Ry]',/,' ',79('='))
99012 FORMAT ((' ',79('*'),/),/,' KKR-run for: ',15(A,F5.2))
99013 FORMAT (/,' results extrapolated to corrected FERMI - ENERGY:',/,
     &        ' CHARGE MISFIT     ',F9.5,/,' E_F CORRECTION    ',F9.5,/,
     &        ' NEW FERMI ENERGY  ',F9.5,/)
99014 FORMAT (/,I4,' E=',2F7.4,3X,'L=',I2,3X,'IT=',I2,2X,A,2X,A20,/,
     &        15X,'DOS  [1/Ry]  |  m_spin  [m_B]  |  m_orb   [m_B]  |',
     &        '   B_tot   [kG]',/,' INT(DE)  ',2F8.3,2X,2F8.3,2X,2F8.3,
     &        F10.1,F8.1,/,' SUM(ML)  ',2F8.3,2X,2F8.3,2X,2F8.3,F10.1,
     &        F8.1,20(:,/,' ML= ',I2,'   ',2F8.3,2X,2F8.3,2X,2F8.3,
     &        F10.1,F8.1))
      END
