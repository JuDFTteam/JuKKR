      subroutine rmtcalc(nt,nq,nm,alat,rws,brx,bry,brz,qx,qy,qz,rmt,imt)
      integer im,jm,nt,i1,i2,i3,iq,jq,nq,imt(20)
      double precision rmt(100),rws(100),brx(3),bry(3),brz(3),
     &     qx(100),qy(100),qz(100),drx,dry,drz,dqx,dqy,dqz,dij,alat,aux
C
C     nm is the number of atoms
C     subroutine from SPRKKR subroutine potfit.f
C
CCCCCCCCCCCC defining of equivalent meshes
      print*,'RMTCALC subroutine'
      DO IT = 1,NT
          print*,'imt(',it,')=',imt(it),'  which mesh type'
      enddo
      if (nm.lt.nt) then
          do i=1,nt-1
              do j=i+1,nt
                  print*,'i=',i,'j=',j
                  if (imt(i).eq.imt(j)) then
                      rws(j)=rws(i)
                      rmt(j)=rmt(i)
                      print*,'rmt(',j,')=',rmt(j)
                  endif
              enddo
          enddo
      endif
C
      DO IM = 1,nt
          RMT(IM) = RWS(IM)
          print*,'RMT(',IM,')=',RMT(IM)
      END DO
      print*,'     '
C     
      DO I1 = -1,1
          DO I2 = -1,1
              DO I3 = -1,1
                  DRX = I1*BRX(1) + I2*BRX(2) + I3*BRX(3)
                  DRY = I1*BRY(1) + I2*BRY(2) + I3*BRY(3)
                  DRZ = I1*BRZ(1) + I2*BRZ(2) + I3*BRZ(3)
C                  print*,'  i1=',i1,'  i2=',i2,'  i3=',i3
C                  print*,'  drx=',drx,'  dry=',dry,'  drz=',drz
C     
                  DO IQ = 1,NQ
C                      IM = IMT(IQ)
C                      IM = IMT(ITOQ(1,IQ))
C                      im=iq
                      DO JQ = 1,NQ
C                          JM = IMT(JQ)
C                          JM = IMT(ITOQ(1,JQ))
C                          jm=jq
                          DQX = DRX + QX(IQ) - QX(JQ)
                          DQY = DRY + QY(IQ) - QY(JQ)
                          DQZ = DRZ + QZ(IQ) - QZ(JQ)
C                          print*,'  dqx=',dqx,'  dqy=',dqy,'  dqz=',dqz
                          DIJ = SQRT(DQX**2+DQY**2+DQZ**2)*ALAT
C                          print*,'dij= ',DIJ
                          IF ( DIJ.GT.1D-8 ) THEN
C                              AUX = DIJ/(RWS(IM)+RWS(JM))
C                              RMT(IM) = MIN(RMT(IM),AUX*RWS(IM))
C                              RMT(JM) = MIN(RMT(JM),AUX*RWS(JM))
                              AUX = DIJ/(RWS(IQ)+RWS(JQ))
                              RMT(IQ) = MIN(RMT(IQ),AUX*RWS(IQ))
                              RMT(JQ) = MIN(RMT(JQ),AUX*RWS(JQ))
C                              print*,'  iq=',iq,'  jq=',jq, '  aux=',aux
C                              print*,' RMT(',IQ,')=',RMT(IQ),' RMT(',JQ,
C     &                             ')=',RMT(JQ)
                          END IF
                      END DO
                  END DO
              END DO
          END DO
      END DO
      end
