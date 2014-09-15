C             New call for main2:
C             call ECOUB_NEW(CMOM,ECOU,LPOT,NSPIN,RHO2NS, &
C             VONS,ZAT(I1),R(:,I1), &
C             DRDI(:,I1),KVMAD,IRCUT(:,I1),IPAN(I1),IMAXSH,IFUNM(1,ICELL), &
C             ILM,GSH,THETAS(:,:,ICELL),LMSP(1,ICELL), &
C             irmd, irid, nfund, ipand, ngshd)


c 13.10.95 ***************************************************************
      SUBROUTINE ECOUB_NEW(CMOM,ECOU,LPOT,NSPIN,RHO2NS,VONS,Z,R,DRDI,
     +                 KVMAD,IRCUT,IPAN,IMAXSH,IFUNM,ILM,
     +                 GSH,THETAS,LMSP,
C                      new parameters after inc.p removal
     &                 irmd, irid, nfund, ipand, ngshd)
c ************************************************************************
c
c     attention : energy zero ---> electro static zero
c
c     calculate the electrostatic potential-energies without the
c     electron-nuclear interaction in the cell itself .
c     the energy of the representive atom i is given by
c
c                          rc
c      ecou(i) =  1/2 (  {  s dr' vons(r',lm,i)*rho2ns(r',lm,i,1) }
c                           0
c
c                                       -  z(i) * vmad ( ri )     )
c
c
c                                         ( {..} = summed over lm )
c             (see notes by b.drittler)
c     vons is the coulomb potential of the atom without the nuclear
c             potential of the atom
c     rho2ns(...,1) is the real charge density times r**2
c
c      both developed into spherical harmonics . (see deck rholm)
c
c     z    is the nuclear charge of the atom
c
c     vmad ( ri ) is a generalized madelung potential
c                 = 1/sqrt(4 pi) * vons(irws,1,is)
c                         - sqrt(4 pi) * 2 * cmom(1,ipot) / rws
c
c                                        ( <..> = spherical averaged )
c
c     attention : this subroutine has to be called before the
c                 exchange correlation potential is added to
c                 the potential vons .
c                 the energy calculated here is splitted into
c                 l-dependent parts to see the l -convergency .
c
c     attention : in case of shape corrections the contribution of
c                 the coulomb potential the of the nucleus is
c                 analytically cancelled only in the muffin tin sphere
c                 in the interstial region it has to be taken into
c                 account ! see deck epotins
c
c                 modified for band structure code
c                               b.drittler   jan. 1990
c-----------------------------------------------------------------------

      IMPLICIT NONE
C     .. Parameters ..
C     ..
C     ..

      INTEGER irmd
      INTEGER irid
      INTEGER nfund
      INTEGER ipand
      INTEGER ngshd

C     .. Scalar Arguments ..
      INTEGER KVMAD,LPOT,NSPIN
C     ..
C     .. Array Arguments ..
C     DOUBLE PRECISION CMOM(LMPOTD),DRDI(IRMD,*),ECOU(0:LPOTD),
C    +                 GSH(*),R(IRMD,*),RHO2NS(IRMD,LMPOTD,2),
C    +                 THETAS(IRID,NFUND,*),VONS(IRMD,LMPOTD,2),Z(*)
C     INTEGER IFUNM(*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
C    +        IRCUT(0:IPAND,*),IRWS(*),LMSP(*)

C     LPOTD = LPOT
C     LMPOTD = (LPOT+1)**2

      DOUBLE PRECISION CMOM((LPOT + 1)**2),DRDI(IRMD),
     &                 ECOU(0:LPOT),
     &                 GSH(ngshd),R(IRMD),
     &                 RHO2NS(IRMD,(LPOT + 1)**2,2),
     &                 THETAS(IRID,NFUND),
     &                 VONS(IRMD,(LPOT + 1)**2,2)

      DOUBLE PRECISION Z

C     TODO: dimension IFUNM?
      INTEGER IFUNM(*),ILM(NGSHD,3),IMAXSH(0:(LPOT + 1)**2),IPAN,
     &        IRCUT(0:IPAND),LMSP(*)

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RFPI,RHOSP,SIGN,VM,VMAD
      INTEGER I,IFUN,IPAN1,IPOT,IR,IRC1,IRH,IRS1,ISPIN,J,L,
     +        LM,LM2,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ER(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3,SIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C     ..
      RFPI = SQRT(16.D0*ATAN(1.0D0))
c
      IPAN1 = IPAN
      IRS1 = IRCUT(1)
      IRC1 = IRCUT(IPAN1)
c
c--->   determine the right potential numbers - the coulomb potential
c       is not spin dependend
c
      IPOT = NSPIN

      DO 80 L = 0,LPOT

        DO 10 I = 1,IRC1
          ER(I) = 0.0D0
   10   CONTINUE

        DO 70 ISPIN = 1,NSPIN

          IF (ISPIN.EQ.NSPIN) THEN
            SIGN = 1.0D0
          ELSE
            SIGN = -1.0D0
          END IF

          DO 60 M = -L,L
            LM = L*L + L + M + 1

            DO 20 I = 1,IRS1
              RHOSP = (RHO2NS(I,LM,1)+
     +                SIGN*RHO2NS(I,LM,NSPIN))/4.0D0
              ER(I) = ER(I) + RHOSP*VONS(I,LM,IPOT)
   20       CONTINUE

c
c--->           convolute with shape function
c
            DO 50 J = IMAXSH(LM-1) + 1,IMAXSH(LM)
              LM2 = ILM(J,2)
              IF (LMSP(ILM(J,3)).GT.0) THEN
                IFUN = IFUNM(ILM(J,3))

                IF (LM2.EQ.1) THEN
                  DO 30 IR = IRS1 + 1,IRC1
                    IRH = IR - IRS1
                    RHOSP = (RHO2NS(IR,LM,1)+
     +                      SIGN*RHO2NS(IR,LM,NSPIN))/2.0D0
c
c--->                 remember that in the interstial -2z/r has 
c                     to be taken into account
c
                    ER(IR) = ER(IR) + RHOSP*GSH(J)*
     +                       THETAS(IRH,IFUN)*
     +                       (VONS(IR,1,IPOT)/2.0D0-
     +                       Z/R(IR)*RFPI)
   30             CONTINUE

                ELSE              ! (LM2.EQ.1)

                  DO 40 IR = IRS1 + 1,IRC1
                    IRH = IR - IRS1
                    RHOSP = (RHO2NS(IR,LM,1)+
     +                      SIGN*RHO2NS(IR,LM,NSPIN))/2.0D0
                    ER(IR) = ER(IR) + RHOSP*GSH(J)*
     +                       THETAS(IRH,IFUN)*VONS(IR,LM2,IPOT)/
     +                       2.0D0
   40             CONTINUE

                END IF            ! (LM2.EQ.1)
              END IF   

   50       CONTINUE

   60     CONTINUE

   70   CONTINUE
c
c--->     now integrate
c
        CALL SIMPK(ER,ECOU(L),IPAN1,IRCUT,
     +             DRDI)

   80   CONTINUE                    ! L = 0,LPOT

C Terms calculated so far:
C     a) 1/2 * \int_0^{r_{MT}} \rho_L V_L r^2 dr +
C     b) 1/2 * \sum_{L'L''} \int_{r_{MT}}^{r_{BS}} C_{LL'L''} \tilde{\rho}_L V_{L'} \Theta_{L''} r^2 dr
C
C     a) muffin-tin b) interstitial
C     in interstitial use: \tilde{\rho} = \rho - \delta_{L,(0,0)} \sqrt{4 \pi} Z * r


c
c--->   calculate the madelung potential
c
        VMAD = VONS(IRS1,1,IPOT)/RFPI -
     +         RFPI*2.0D0*CMOM(1)/R(IRS1)
c
c--->   add to ecou
c
        ECOU(0) = ECOU(0) - Z*VMAD/2.0D0

C     Madelung term added: (why is it incomplete?)
C     \delta_{L,(0,0)} (-Z)/2 * \[ V_{(0,0)}(R)/sqrt{4 \pi} - \sqrt{4 \pi} * 2 * q_{0,0}(R) / R \]

c
c--->   option to calculate full generalized madelung potential
c                                    rc
c       vm(rn) = vmad +2*sqrt(4*pi)* s  dr*r*rho(lm=1,r)
c                                    0
        IF (KVMAD.EQ.1) THEN
          ER(1) = 0.0D0

          DO 90 I = 2,IRS1
            ER(I) = RHO2NS(I,1,1)/R(I)
   90     CONTINUE

          CALL SIMP3(ER,VM,1,IRS1,DRDI)
          VM = 2.0D0*RFPI*VM + VMAD
c
c         atom nr. iatyp is the iatyp-th atom on the potential cards
c         e. g., in binary alloys iatyp=1 and iatyp=2 refer to host
c
C         WRITE (6,FMT=9010) IATYP,VMAD
C         WRITE (6,FMT=9000) IATYP,VM
        END IF                      ! (KVMAD.EQ.1)

      RETURN

 9000 FORMAT (10x,'full generalized madelung pot. for atom',1x,i3,1x,
     +       ': ',1p,d14.6)
 9010 FORMAT (10x,'     generalized madelung pot. for atom',1x,i3,1x,
     +       ': ',1p,d14.6)
      END
