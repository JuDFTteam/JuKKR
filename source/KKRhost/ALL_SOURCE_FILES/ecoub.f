c 13.10.95 ***************************************************************
      SUBROUTINE ECOUB(CMOM,ECOU,LMAX,NSPIN,NATYP,RHO2NS,VM2Z,Z,R,DRDI,
     +                 IRWS,KVMAD,KSHAPE,IRCUT,IPAN,IMAXSH,IFUNM,ILM,
     +                 NTCELL,GSH,THETAS,LMSP)
c ************************************************************************
c
c     attention : energy zero ---> electro static zero
c
c     calculate the electrostatic potential-energies without the
c     electron-nuclear interaction in the cell itself .
c     the energy of the representive atom i is given by
c
c                          rc
c      ecou(i) =  1/2 (  {  s dr' vm2z(r',lm,i)*rho2ns(r',lm,i,1) }
c                           0
c
c                                       -  z(i) * vmad ( ri )     )
c
c
c                                         ( {..} = summed over lm )
c             (see notes by b.drittler)
c     vm2z is the coulomb potential of the atom without the nuclear
c             potential of the atom
c     rho2ns(...,1) is the real charge density times r**2
c
c      both developed into spherical harmonics . (see deck rholm)
c
c     z    is the nuclear charge of the atom
c
c     vmad ( ri ) is a generalized madelung potential
c                 = 1/sqrt(4 pi) * vm2z(irws,1,is)
c                         - sqrt(4 pi) * 2 * cmom(1,ipot) / rws
c
c                                        ( <..> = spherical averaged )
c
c     attention : this subroutine has to be called before the
c                 exchange correlation potential is added to
c                 the potential vm2z .
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
C     .. Parameters ..
      include 'inc.p'
C     ..
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHAPE,KVMAD,LMAX,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CMOM(LMPOTD,*),DRDI(IRMD,*),ECOU(0:LPOTD,*),
     +                 GSH(*),R(IRMD,*),RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 THETAS(IRID,NFUND,*),VM2Z(IRMD,LMPOTD,*),Z(*)
      INTEGER IFUNM(NATYPD,*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
     +        IRCUT(0:IPAND,*),IRWS(*),NTCELL(*),LMSP(NATYPD,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RFPI,RHOSP,SIGN,VM,VMAD
      INTEGER I,IATYP,ICELL,IFUN,IPAN1,IPOT,IR,IRC1,IRH,IRS1,ISPIN,J,L,
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
      DO 100 IATYP = 1,NATYP

        IF (KSHAPE.NE.0) THEN
          IPAN1 = IPAN(IATYP)
          ICELL = NTCELL(IATYP)
          IRS1 = IRCUT(1,IATYP)
          IRC1 = IRCUT(IPAN1,IATYP)
        ELSE
          IRS1 = IRWS(IATYP)
          IRC1 = IRS1
        END IF
c
c--->   determine the right potential numbers - the coulomb potential
c       is not spin dependend
c
        IPOT = IATYP*NSPIN

        DO 80 L = 0,LMAX

          DO 10 I = 1,IRC1
            ER(I) = 0.0D0
   10     CONTINUE

          DO 70 ISPIN = 1,NSPIN

            IF (ISPIN.EQ.NSPIN) THEN
              SIGN = 1.0D0
            ELSE
              SIGN = -1.0D0
            END IF

            DO 60 M = -L,L
              LM = L*L + L + M + 1

              DO 20 I = 1,IRS1
                RHOSP = (RHO2NS(I,LM,IATYP,1)+
     +                  SIGN*RHO2NS(I,LM,IATYP,NSPIN))/4.0D0
                ER(I) = ER(I) + RHOSP*VM2Z(I,LM,IPOT)
   20         CONTINUE

              IF (KSHAPE.NE.0) THEN
c
c--->           convolute with shape function
c
                DO 50 J = IMAXSH(LM-1) + 1,IMAXSH(LM)
                  LM2 = ILM(J,2)
                IF (LMSP(ICELL,ILM(J,3)).GT.0) THEN
                  IFUN = IFUNM(ICELL,ILM(J,3))

                  IF (LM2.EQ.1) THEN
                    DO 30 IR = IRS1 + 1,IRC1
                      IRH = IR - IRS1
                      RHOSP = (RHO2NS(IR,LM,IATYP,1)+
     +                        SIGN*RHO2NS(IR,LM,IATYP,NSPIN))/2.0D0
c
c--->                 remember that in the interstial -2z/r has 
c                     to be taken into account
c
                      ER(IR) = ER(IR) + RHOSP*GSH(J)*
     +                         THETAS(IRH,IFUN,ICELL)*
     +                         (VM2Z(IR,1,IPOT)/2.0D0-
     +                         Z(IATYP)/R(IR,IATYP)*RFPI)
   30               CONTINUE

                  ELSE              ! (LM2.EQ.1)

                    DO 40 IR = IRS1 + 1,IRC1
                      IRH = IR - IRS1
                      RHOSP = (RHO2NS(IR,LM,IATYP,1)+
     +                        SIGN*RHO2NS(IR,LM,IATYP,NSPIN))/2.0D0
                      ER(IR) = ER(IR) + RHOSP*GSH(J)*
     +                         THETAS(IRH,IFUN,ICELL)*VM2Z(IR,LM2,IPOT)/
     +                         2.0D0
   40               CONTINUE

                  END IF            ! (LM2.EQ.1)
                END IF   

   50           CONTINUE

              END IF                ! (KSHAPE.NE.0)

   60       CONTINUE

   70     CONTINUE
c
c--->     now integrate
c
          IF (KSHAPE.EQ.0) THEN
            CALL SIMP3(ER,ECOU(L,IATYP),1,IRS1,DRDI(1,IATYP))
          ELSE
            CALL SIMPK(ER,ECOU(L,IATYP),IPAN1,IRCUT(0,IATYP),
     +                 DRDI(1,IATYP))
          END IF

   80   CONTINUE                    ! L = 0,LMAX

c
c--->   calculate the madelung potential
c
        VMAD = VM2Z(IRS1,1,IPOT)/RFPI -
     +         RFPI*2.0D0*CMOM(1,IATYP)/R(IRS1,IATYP)
c
c--->   add to ecou
c
!         write(*,*) 'test',VM2Z(IRS1,1,IPOT),CMOM(1,IATYP),R(IRS1,IATYP)
!           write(*,*) 'test',ECOU(0,IATYP) , Z(IATYP),VMAD,2.0D0
        ECOU(0,IATYP) = ECOU(0,IATYP) - Z(IATYP)*VMAD/2.0D0
c
c--->   option to calculate full generalized madelung potential
c                                    rc
c       vm(rn) = vmad +2*sqrt(4*pi)* s  dr*r*rho(lm=1,r)
c                                    0
        IF (KVMAD.EQ.1) THEN
          ER(1) = 0.0D0

          DO 90 I = 2,IRS1
            ER(I) = RHO2NS(I,1,IATYP,1)/R(I,IATYP)
   90     CONTINUE

          CALL SIMP3(ER,VM,1,IRS1,DRDI(1,IATYP))
          VM = 2.0D0*RFPI*VM + VMAD
c
c         atom nr. iatyp is the iatyp-th atom on the potential cards
c         e. g., in binary alloys iatyp=1 and iatyp=2 refer to host
c
          WRITE (6,FMT=9010) IATYP,VMAD
          WRITE (6,FMT=9000) IATYP,VM
        END IF                      ! (KVMAD.EQ.1)

  100 CONTINUE                      ! IATYP = 1,NATYP

      RETURN

 9000 FORMAT (10x,'full generalized madelung pot. for atom',1x,i3,1x,
     +       ': ',1p,d14.6)
 9010 FORMAT (10x,'     generalized madelung pot. for atom',1x,i3,1x,
     +       ': ',1p,d14.6)
      END
