       SUBROUTINE INTERPOLATE_POTEN(LPOTD,IRMD,IRNSD,NATYPD,IPAND,
     +                              NSPOTD,NTOTD,NCHEBD,IRMDNEW,
     +                              NSPIN,R,IRMIN,IRWS,IRCUT,VINS,
     +                              VISP,NPAN_LOG,NPAN_EQ,
     +                              NPAN_TOT,RNEW,
     +                              IPAN_INTERVALL,VINSNEW)
       IMPLICIT NONE
       !INCLUDE 'inc.p'
       !include array dimensions in interface explicityly to get rid of
       !inc.p import and to be able to use routine for different number
       !of atoms
       integer, intent(in) :: LPOTD,IRMD,IRNSD,NATYPD,IPAND,
     +                        NSPOTD,NTOTD,NCHEBD
       INTEGER LMPOTD
       !PARAMETER (LMPOTD= (LPOTD+1)**2)
       INTEGER IRMIND
       !PARAMETER (IRMIND = IRMD-IRNSD)
       INTEGER NSPIN,IRMIN(NATYPD),IRWS(NATYPD),IRMDNEW
       INTEGER IRCUT(0:IPAND,NATYPD)
       INTEGER NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),NPAN_TOT(NATYPD)
       DOUBLE PRECISION R(IRMD,NATYPD)
       DOUBLE PRECISION 
     +   VINS((IRMD-IRNSD):IRMD,(LPOTD+1)**2,NSPOTD),
!     +   VINS(IRMIND:IRMD,LMPOTD,NSPOTD),
     +   VISP(IRMD,NSPOTD),
     +   VINSIN(IRMD,(LPOTD+1)**2,NSPIN)
!     +   VINSIN(IRMD,LMPOTD,NSPIN)
       DOUBLE PRECISION RNEW(IRMDNEW,NATYPD)
       !DOUBLE PRECISION RNEW(NTOTD*(NCHEBD+1),NATYPD)
       INTEGER IPAN_INTERVALL(0:NTOTD,NATYPD)
       DOUBLE PRECISION VINSNEW(IRMDNEW,(LPOTD+1)**2,NSPOTD)
!       DOUBLE PRECISION VINSNEW(NTOTD*(NCHEBD+1),LMPOTD,NSPOTD)
       INTEGER I1,IPOT,IPOTM,IMIN,IMAX,IP,IR,LM1,ISPIN,
     +         IMINNEW,IMAXNEW,IR2

       LMPOTD= (LPOTD+1)**2
       IRMIND = IRMD-IRNSD

       VINSNEW=0d0
       IPOTM=0

            
c interpolate potential to new mesh
       DO I1 = 1,NATYPD

        IPOT=NSPIN*(I1-1)+1

c save input potential to VINSIN
        VINSIN=0d0
        DO IR=1,IRWS(I1)
         IF (IR.LT.IRMIN(I1)) THEN
          VINSIN(IR,1,1)=VISP(IR,IPOT)
          VINSIN(IR,1,NSPIN)=VISP(IR,IPOT+NSPIN-1)
         ELSE  
          DO LM1=1,LMPOTD
           IF (LM1.EQ.1) THEN 
            VINSIN(IR,LM1,1)=VISP(IR,IPOT)
            VINSIN(IR,LM1,NSPIN)=VISP(IR,IPOT+NSPIN-1)
           ELSE
            VINSIN(IR,LM1,1)=VINS(IR,LM1,IPOT)
            VINSIN(IR,LM1,NSPIN)=VINS(IR,LM1,IPOT+NSPIN-1)
           ENDIF
          ENDDO
         ENDIF
        ENDDO

        DO ISPIN=1,NSPIN

         IPOTM=IPOTM+1

         DO LM1=1,LMPOTD

          IMIN=1
          IMAX=IRMIN(I1)
          DO IP=1,NPAN_LOG(I1)
           IMINNEW=IPAN_INTERVALL(IP-1,I1)+1
           IMAXNEW=IPAN_INTERVALL(IP,I1)
           CALL INTERPOLSPLINE(R(IMIN:IMAX,I1),RNEW(IMINNEW:IMAXNEW,I1),
     +           VINSIN(IMIN:IMAX,LM1,ISPIN),
     +           VINSNEW(IMINNEW:IMAXNEW,LM1,IPOTM),
     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
          ENDDO

          IMIN=IRMIN(I1)
          IMAX=IRCUT(1,I1)
          DO IP=NPAN_LOG(I1)+1,NPAN_LOG(I1)+NPAN_EQ(I1)
           IMINNEW=IPAN_INTERVALL(IP-1,I1)+1
           IMAXNEW=IPAN_INTERVALL(IP,I1)
           CALL INTERPOLSPLINE(R(IMIN:IMAX,I1),RNEW(IMINNEW:IMAXNEW,I1),
     +           VINSIN(IMIN:IMAX,LM1,ISPIN),
     +           VINSNEW(IMINNEW:IMAXNEW,LM1,IPOTM),
     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
          ENDDO

          IR2=0
          DO IP=NPAN_LOG(I1)+NPAN_EQ(I1)+1,NPAN_TOT(I1)
           IR2=IR2+1
           IMIN=IRCUT(IR2,I1)+1
           IMAX=IRCUT(IR2+1,I1)
           IMINNEW=IPAN_INTERVALL(IP-1,I1)+1
           IMAXNEW=IPAN_INTERVALL(IP,I1)
          CALL INTERPOLSPLINE(R(IMIN:IMAX,I1),RNEW(IMINNEW:IMAXNEW,I1),
     +           VINSIN(IMIN:IMAX,LM1,ISPIN),
     +           VINSNEW(IMINNEW:IMAXNEW,LM1,IPOTM),
     +           IMAX-IMIN+1,IMAXNEW-IMINNEW+1)
          ENDDO 
         ENDDO ! LM1
        ENDDO ! ISPIN
       ENDDO ! I1
       END
