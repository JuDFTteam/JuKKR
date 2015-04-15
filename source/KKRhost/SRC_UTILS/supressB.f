      program supressB
      implicit none
C-------------------------------------------------------------
C     supressB $old.pot $new.pot
C     overwrites the UP-part by the DN-part for each atom
C-------------------------------------------------------------
C
      double precision
     &     alat1,rwsinp(20),rmtinp(20),rmtnewinp(20),z(20),efnew,vbc(40)
     &     ,vm2z(1000,40),ecore(20,20),a1(20),b1(20),magshift(600,20)
     &     ,xcsplit
      integer irws(40),inew,lcore(20,40),icore,ncore(40),ir
      integer natm,i,ih,nspin,ispin
      character*256 infile,outfile
      character*64 title
      xcsplit=0.05
C ------------------------------------------------------------
C     asks for initial data
      print 10
 10   format(' give the number of atoms: ',$)
      read(5,*) natm
C     ------OPEN files --------------
      call getarg(1,infile)
      open(unit=3,file=infile)
      call getarg(2,outfile)
      open(unit=4,file=outfile)      
C     CYCLE - natm * 2spins
C     ih means which atom is considered
C     i  means ih spin-resolved
      nspin=2
      do ih = 1,natm
          do ispin = 1,nspin
              I = NSPIN* (IH-1) + ISPIN
C     --------------INPUT------------
              read(3,'(a)')title
              print *,title
              read(3,fmt=9030)rmtinp(ih),alat1,rmtnewinp(ih)
              read(3,fmt=9040)z(ih),rwsinp(ih),efnew,vbc(ispin)
              read(3,fmt=9050)irws(IH),a1(IH),b1(IH),ncore(I),inew
              IF (NCORE(i).GE.1) THEN
                  DO ICORE=1,NCORE(I)
                      READ (3,FMT=9070) LCORE(ICORE,i),
     &                     ECORE(ICORE,i)
                  END DO
              END IF
              READ (3,FMT=*) (VM2Z(IR,i),IR=1,irws(ih))
C     --------------------------------------------------
C     if i = odd  -> WRITE
C     if i = even -> i=i-1
C     
              if (dble(i)/2.d0-i/2.eq.0.d0) i=i-1
C     
C     --------------OUTPUT------------
              write(4,'(a)')title
              write(4,fmt=9030)rmtinp(ih),alat1,rmtnewinp(ih)
              write(4,fmt=9040)z(ih),rwsinp(ih),efnew,vbc(ispin)
              write(4,fmt=9050)irws(IH),a1(IH),b1(IH),ncore(i),inew
C     
              IF (NCORE(i).GE.1) THEN
                  DO ICORE=1,NCORE(I)
                      write (4,FMT=9070) LCORE(ICORE,i),
     &                     ECORE(ICORE,i)
                  END DO
              END IF
              write (4,fmt=9051) (VM2Z(IR,i),IR=1,irws(ih))
          enddo
      enddo
      close(3)
      close(4)
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i3,/,2d15.8,/,2i2)
 9070 FORMAT (i5,1p,d20.11)
 9051 FORMAT (1p,4d20.12)
      end
