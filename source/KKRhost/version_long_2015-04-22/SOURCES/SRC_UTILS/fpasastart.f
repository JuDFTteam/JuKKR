      program fpasastart
C-------------------------------------------------------------
C     fpasastart $old.pot $new.pot $shapefun
C-------------------------------------------------------------
C Reads in ASA potential and write out starting FP.
C On the input are needed: new number of r points (standart 484)
C                          IRNS1P,LMPOTP (standart 208, 25)
C                          new rws
C                          new rmt
C ! Until now only one critical radius allowed
      double precision alat,rws(200),rmt(200),rmtnew(200),z(100),efnew
     &     ,vbc(100),vm2z(1000,200),ecore(20,200),a(100),b(100),rcrit
     &     ,rmtfp,RIN(1000,200),ROUT(1000,200),anew(100),bnew(100)
     &     ,vm2znew(1000,200),U1(1000)
      integer irws(100),inew,lcore(20,200),icore,ncore(200),ir
      integer natm,i,ih,nspin,ispin
      DOUBLE PRECISION DRN(1000,100),SCALE(100),
     +                 XRN(1000,100)
      INTEGER MESHN(100),NM(100,100),NPAN(100),NFU(100)
      character*256 infile,outfile
      character*64 title(200)
      REAL*8 YLAG
C ------------------------------------------------------------
C     ----------INITIAL DATA--------------------------
      write(6,*) ' give the number of atoms: '
      read(5,*) natm
      write(6,*) ' give new rws (rcrit): '
      read(5,*) rcrit
      write(6,*) ' give new rmt: '
      read(5,*) rmtfp
      nrnew=484
      IRNS1P=208
      LMPOTP=25
      write(6,*) "give new nr,IRNS1P,LMPOTP "
      call flush(6)
      read(5,*) nrnew,IRNS1P,LMPOTP
      write(*,*) nrnew,IRNS1P,LMPOTP
C     ------OPEN files --------------
      call getarg(1,infile)
      open(unit=3,file=infile)
C     CYCLE - natm * 2spins
C     ih means which atom is considered
C     i  means ih spin-resolved
      nspin=2
C     --------------INPUT------------
      do ih = 1,natm
          do ispin = 1,nspin
              I = NSPIN* (IH-1) + ISPIN
              read(3,'(a)')title(i)
              print *,title(i)
              read(3,fmt=9030)rmt(ih),alat,rmtnew(ih)
              read(3,fmt=9040)z(ih),rws(ih),efnew,vbc(ispin)
              read(3,fmt=9050)irws(IH),a(IH),b(IH),ncore(I),inew
              IF (NCORE(i).GE.1) THEN
                  DO ICORE=1,NCORE(I)
                      READ (3,FMT=9070) LCORE(ICORE,i),
     &                     ECORE(ICORE,i)
                  END DO
              END IF
              READ (3,FMT=*) (VM2Z(IR,i),IR=1,irws(ih))
          enddo
      enddo
      close(3)

C     --------------SHAPEFUN------------
      print*,'reading in the shapefunfile'
      call getarg(3,outfile)
      open(unit=4,file=outfile)
      READ (4,FMT=9000) NCELL
      WRITE (6,FMT=*) '  ncell : ',NCELL
      READ (4,FMT=9010) (SCALE(ICELL),ICELL=1,NCELL)
      DO ICELL = 1,NCELL
          READ (4,FMT=9000) NPAN(ICELL),MESHN(ICELL)   
          WRITE (*,*) "NPAN,MESHN",NPAN(ICELL),MESHN(ICELL)   
          READ (4,FMT=9000) (NM(IPAN1,ICELL),IPAN1=2,NPAN(ICELL)+1)
          READ (4,FMT=9010) (XRN(IR,ICELL),DRN(IR,ICELL),IR=1,
     +         MESHN(ICELL))

          READ (4,FMT=9000) NFU(ICELL)
          NFUN = NFU(ICELL)
          WRITE (6,FMT=*) '  nfun  : ',NFUN,NFUND
          DO IFUN = 1,NFUN
              READ (4,FMT=9000) LM
              READ (4,FMT=9010) (rdum,N=1,MESHN(ICELL))
          END DO
      END DO
      close(4)
c
C SETTING UP INPUT AND OUUPUT MESH
      do ih = 1,natm
          do IR=1,irws(IH)
              rin(ir,ih)=b(ih)*(dexp(a(ih)*(dble(IR)-1.0d0))-1
     &             .0d0)
          end do
C     GET NEW IMT1 and B
          imt1=ANINT(LOG(rmtfp/b(1)+1.0D0)/a(1)) + 1
          imt1=ANINT(LOG(rmtfp/b(1)+1.0D0)/a(1)) + 1
          WRITE(*,*) IMT1,"              IMT1"
          anew(ih)=a(ih)
          IF (MOD(IMT1,2).EQ.0) THEN
              IMT1 = IMT1 + 1
          END IF
          bnew(ih)=rmtfp/(EXP(ANEW(IH)*REAL(IMT1-1))-1.0D0)
          write(*,*) bnew(ih),"           bnew"
          do IR=1,irws(IH)
              rout(ir,ih)=bnew(ih)*(dexp(anew(ih)*(dble(IR)-1.0d0))-1
     &             .0d0)
          end do
          IRNS=0
          DO ICELL=1,NCELL
              IRNS=IRNS+MESHN(ICELL)
          DO IRI = 1,MESHN(ICELL)
              IR = IRI + IMT1
              ROUT(IR,IH) = SCALE(ICELL)*ALAT*XRN(IRI,ICELL)
          END DO
          END DO
          IRMAX=IR
      end do
      do ih = 1,natm
          do ispin = 1,nspin
              I = NSPIN* (IH-1) + ISPIN
              VM2ZNEW(1,i)=VM2Z(1,i)
              DO IR=1,irws(IH)
                  U1(IR)=VM2Z(IR,i)
              END DO
              DO IR=2,IRMAX
                  VM2ZNEW(IR,i)=YLAG(ROUT(ir,ih),RIN(1,ih),U1,0,3
     &                 ,irws(IH))
              END DO
              DO IR=1,IRMAX
                  write(99,*) ir,rout(ir,ih),rin(ir,ih),VM2ZNEW(IR,i)
     &             ,VM2Z(IR,i)
              end do
          end do
      end do
      if (IRMAX.NE.nrnew) stop "IRMAX.NE.nrnew"
C     --------------OUTPUT------------
      print*,'writing out the potential'
      call getarg(2,outfile)
      open(unit=4,file=outfile)      
      do ih = 1,natm
          do ispin = nspin,1,-1
              I = NSPIN* (IH-1) + ISPIN
              print *,title(i)
              write(4,'(a)')title(i)
              write(4,fmt=9030)rmtfp,alat,rmtfp
              write(4,fmt=9040)z(ih),rcrit,efnew,vbc(ispin)
              write(4,fmt=9050)irmax,anew(IH),bnew(IH),ncore(i),inew
              IF (NCORE(i).GE.1) THEN
                  DO ICORE=1,NCORE(I)
                      write (4,FMT=9070) LCORE(ICORE,i),
     &                     ECORE(ICORE,i)
                  END DO
              END IF
              WRITE (4,FMT=9060) irmax,IRNS1P,LMPOTP,1
              write (4,fmt=9051) (VM2ZNEW(IR,i),IR=1,irmax)
              WRITE (4,FMT=9060) 1
          enddo
      enddo
      close(4)
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i3,/,2d15.8,/,2i2)
 9070 FORMAT (i5,1p,d20.11)
C 9051 FORMAT (1p,4d20.12)
 9051 FORMAT (1p,4D20.13)
 9000 FORMAT (16i5)
 9010 FORMAT (4d20.12)
 9060 FORMAT (10i5)
      end
      



      FUNCTION YLAG(XI,X,Y,IND1,N1,IMAX)
C   ********************************************************************
C   *                                                                  *
C   * lagrangian interpolation                                         *
C   * xi is interpolated entry into x-array                            *
C   * n is the order of lagrangran interpolation                       *
C   * y is array from which ylag is obtained by interpolation          *
C   * ind is the min-i for x(i).gt.xi                                  *
C   * if ind=0,x-array will be searched                                *
C   * imax is max index of x-and y-arrays                              *
C   *                                                                  *
C   * 07/12/94  HE  arg. IEX removed                                   *
C   ********************************************************************
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER IMAX,IND1,N1
      REAL*8 XI
      REAL*8 X(IMAX),Y(IMAX)
      REAL*8 YLAG
C
C Local variables
C
      REAL*8 D,P,S,XD
      INTEGER I,IND,INL,INU,J,N
      SAVE D,I,IND,INL,INU,J,N,P,S,XD
C
      IND = IND1
      N = N1
      IF ( N.GT.IMAX ) N = IMAX
      IF ( IND.GT.0 ) GOTO 200
      DO J = 1,IMAX
         IF ( ABS(XI-X(J)).LT.1.0D-12 ) GOTO 600
C      IF( ABS(XI-X(J)).LT.1.0e-06) GO TO 130
         IF ( XI.LT.X(J) ) GOTO 100
         IF ( XI.EQ.X(J) ) GOTO 600
      END DO
      GOTO 300
 100  CONTINUE
      IND = J
 200  CONTINUE
      IF ( IND.GT.1 ) THEN
      END IF
      INL = IND - (N+1)/2
      IF ( INL.LE.0 ) INL = 1
      INU = INL + N - 1
      IF ( INU.LE.IMAX ) GOTO 400
 300  CONTINUE
      INL = IMAX - N + 1
      INU = IMAX
 400  CONTINUE
      S = 0.0D0
      P = 1.0D0
      DO J = INL,INU
         P = P*(XI-X(J))
         D = 1.0D0
         DO I = INL,INU
            IF ( I.NE.J ) THEN
               XD = X(J)
            ELSE
               XD = XI
            END IF
            D = D*(XD-X(I))
         END DO
         S = S + Y(J)/D
      END DO
      YLAG = S*P
 500  CONTINUE
      RETURN
 600  CONTINUE
      YLAG = Y(J)
      GOTO 500
      END
C
C
