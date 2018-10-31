      SUBROUTINE GRADRL(NSPIN,MESH,L1MAX,DX,RHOL,RV,DRDI,IPAN,IPAND,
     +                  IRCUT,DRRL,DDRRL,DRRUL,DDRRUL,IRMD,LMPOTD)
c.....------------------------------------------------------------------
c     gradient of rl with rl defined by charge density=sum(rl*ylm).
c     mesh,l1max: max of mesh and l+1.
c     IRMD,LMPOTD: maxima of corresponding dimension parameters.
c     drrl=d(rl)/dr, ddrrl=d(drrl)/dr, drrul=d(rl-up)/dr,
c     ztal: zeta for each l-component necessary to get down-components.
c.....------------------------------------------------------------------
c.....------------------------------------------------------------------
      use mod_types, only: t_inc
      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION ZERO,ZERO1
      PARAMETER (ZERO=0.D0,ZERO1=1.D-12)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION DX
      INTEGER IPAN,IPAND,IRMD,L1MAX,LMPOTD,MESH,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DDRRL(IRMD,LMPOTD),DDRRUL(IRMD,LMPOTD),
     +                 DRDI(IRMD),DRRL(IRMD,LMPOTD),DRRUL(IRMD,LMPOTD),
     +                 RHOL(IRMD,2,LMPOTD),RV(IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CHGDEN,PI,R2,S4,SPIDEN
      INTEGER I1,IEN,IP,IR,IST,LLMAX
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRDI2(IRMD),RL1(IRMD),RL1UDM(IRMD),ZTAL(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL GRADR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ACOS,SQRT
C     ..
c.....------------------------------------------------------------------
      PI = ACOS(-1.d0)
      S4 = SQRT(4.D0*PI)
      LLMAX = L1MAX*L1MAX
c
      DO 30 IP = 1,IPAN
        IST = IRCUT(IP-1) + 1
        IEN = IRCUT(IP)
        if(t_inc%i_write>0) WRITE (1337,FMT=9010) IP,IST,IEN
        IF (IP.EQ.1) THEN
          DO 10 IR = IST,IEN
            DRDI2(IR) = DX
   10     CONTINUE
        ELSE
          DO 20 IR = IST,IEN
            DRDI2(IR) = ZERO
   20     CONTINUE
        END IF
   30 CONTINUE
c
c
      DO 110 I1 = 1,LLMAX


        IF (NSPIN.EQ.1) GO TO 50

        DO 40 IR = 2,MESH
          R2 = RV(IR)*RV(IR)
c        rl1(ir)=EXP(rv(ir))
c        if(nspin.eq.1) then
c        rl1(ir)=rhol(ir,1,I1)*2.e0
c        ztal(ir)=zero
c        else
c          rl1(ir)=EXP(rv(ir))
          CHGDEN = RHOL(IR,1,I1) + RHOL(IR,2,I1)
          SPIDEN = RHOL(IR,2,I1) - RHOL(IR,1,I1)
          IF (ABS(CHGDEN).GE.ZERO1) THEN
            RL1(IR) = CHGDEN
            ZTAL(IR) = SPIDEN/CHGDEN

          ELSE
            RL1(IR) = ZERO
            ZTAL(IR) = ZERO

          END IF
   40   CONTINUE

        GO TO 70


   50   CONTINUE
        DO 60 IR = 2,MESH
          R2 = RV(IR)*RV(IR)
          RL1(IR) = RHOL(IR,1,I1) + RHOL(IR,2,I1)
          ZTAL(IR) = ZERO
C     check
C
C        rl1(ir)=EXP(rv(ir))
C        ztal(ir)=zero
C
C
   60   CONTINUE

   70   CONTINUE
c
C        write(6,9011) I1,rv(ir),EXP(rv(ir)),rhol(ir,1,I1),
c    &                 rhol(ir,2,I1)
c9011    format(1x, ' I1 rv exp rhol',I5,4e15.6)
c
        RL1(1) = RL1(2)
        ZTAL(1) = ZTAL(2)
c
C       write(6,*) 'mesh',mesh
C       write(6,9000) (ir,rv(ir),rl1(ir),ztal(ir),
C    &     drdi(ir),drdi2(ir),ir=1,mesh,20)
C9000   format(1x,' ir rv rl1 ztal drdi drdi2',i5,5F15.5)
c

        DO 100 IP = 1,IPAN
          IST = IRCUT(IP-1) + 1
          IEN = IRCUT(IP)
c
C      write(6,9005) ip,ist,ien,I1
C9005  format(1x,/,' ip=',i5,' ist=',i5, ' ien=',i5,' I1=',i5)
C       write(6,*) ' before gradr '
          CALL GRADR(NSPIN,IST,IEN,1.d0,DRDI,DRDI2,RL1,ZTAL,DRRL(1,I1),
     +               DDRRL(1,I1),DRRUL(1,I1),DDRRUL(1,I1),RL1UDM,IRMD)
c
          IF (IP.EQ.1) THEN
            DO 80 IR = 1,4
              DRRL(IR,I1) = DRRL(5,I1)
              DDRRL(IR,I1) = DDRRL(5,I1)
              DRRUL(IR,I1) = DRRUL(5,I1)
              DDRRUL(IR,I1) = DDRRUL(5,I1)
   80       CONTINUE
          END IF
C
          IF (NSPIN.EQ.1) THEN
            DO 90 IR = IST,IEN
              DRRUL(IR,I1) = DRRL(IR,I1)/2.D0
              DDRRUL(IR,I1) = DDRRL(IR,I1)/2.D0
   90       CONTINUE
          END IF
C

c       if(I1.eq.1.or.I1.eq.4) then
c       write(6,9001) (ir,rv(ir),rl1(ir),drrl(ir,I1),ddrrl(ir,I1),
c    &    drrul(ir,I1),ddrrul(ir,I1),ir=ist,ien,20)
c9001 format(1x,' ir rv rl1 drrl ddrrl drrul ddrrul',i5,6F12.5)
c        end if

C     stop77
  100   CONTINUE

  110 CONTINUE
C
c     iir=10
c     llm=4
c     write(6,9996) iir,llm
c9996 format(1x,' iir llm',10i5)
c     write(6,9997) drrl(iir,llm),ddrrl(iir,llm),
c    &              drrul(iir,llm),ddrrul(iir,llm)
c9997 format(1x,' drrl ddrrl drrul ddrrul',5f10.5)
C
      RETURN
 9000 FORMAT (1x,' l1max=',i5,' mesh=',i5,'nspi=',i5,' ipan=',i5)
 9010 FORMAT (1x,'  ip ist ien',3i5)
      END
