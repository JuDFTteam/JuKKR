      program number_potential
      implicit none
C
C     changes the order of atoms in potential
C     according to the order in `inputcard`
C
      double precision alat,rws(200),rmt(200),rmtnew(200),z(100),efnew
     &     ,vbc(100),vm2z(1000,200),ecore(100,200),a(100),b(100)
      integer irws(100),inew,lcore(100,200),icore,ncore(200),ir
      integer natm,natmnew,i,ih,nspin,ispin,order(200),iorder
      character*256 infile,outfile
      character*64 title(200)
C ------------------------------------------------------------
C     ----------INITIAL DATA--------------------------
      write(6,*) ' how many atoms in the input potential: '
      read(5,*) natm
      write(6,*) ' how many atoms should be in the output potential:'
      read(5,*) natmnew
      write(6,*) ' what is the order (DIMENSION = natmnew) '
      read(5,*) (order(i),i=1,natmnew)
      print*,(order(i),i=1,natmnew)
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
C     --------------OUTPUT------------
      print*,'writing out the potential'
      call getarg(2,outfile)
      open(unit=4,file=outfile)      
      do ih = 1,natmnew
          do ispin = 1,nspin
              iorder=order(ih)
              print*,'iorder=',iorder
              I = NSPIN* (iorder-1) + ISPIN
              print *,title(i)
              write(4,'(a)')title(i)
              write(4,fmt=9030)rmt(iorder),alat,rmtnew(iorder)
              write(4,fmt=9040)z(iorder),rws(iorder),efnew,vbc(ispin)
              write(4,fmt=9050)irws(iorder),a(iorder),b(iorder),ncore(i)
     &             ,inew
              IF (NCORE(i).GE.1) THEN
                  DO ICORE=1,NCORE(I)
                      write (4,FMT=9070) LCORE(ICORE,i),
     &                     ECORE(ICORE,i)
                  END DO
              END IF
              write (4,fmt=9051) (VM2Z(IR,i),IR=1,irws(iorder))
          enddo
      enddo
      close(4)
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i3,/,2d15.8,/,2i2)
 9070 FORMAT (i5,1p,d20.11)
 9051 FORMAT (1p,4d20.12)
      end
      
