      subroutine tbkkr_pformat(natm,rmtout,alat,z,rwsout,
     &     efnew,irws,a,b,VRUP,VRDN,title)
C     --------------------------------------------------------------
C     writing the potential for TBKKR Juelich
C     
      double precision alat,rwsout(100),rmtout(100),rmtnewout(100),
     &     z(100),efnew,vbc,vm2z(1000,100),ecore(20,200),a(100)
     &     ,b(100)
      integer irws(100),inew,lcore(20,200),icore,ncore(200),ir
      integer natm,jmt(200),i,ih,nspin,ispin,NQNTAB(15),LQNTAB(15)
     &     ,LCMAX
      character*80 title
      character*4 updn
      real*8 VRDN(1000,100),VRUP(1000,100)
      CHARACTER*2 TABCHSYM
      DATA NQNTAB/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6/
      DATA LQNTAB/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1/

C
      LCMAX=3
      inew=1
      nspin=2
      vbc=1.0d0
      print*,'writing out the TBKKR potential'
C
      do ih = 1,natm
          do ispin = 1,nspin
              I = NSPIN* (IH-1) + ISPIN        
C     
C     RMTOUT has been calculated by the subroutine rmtcalc 
C              rmtout(ih)=alat*sqrt(3.d0)/4.d0 ! only for bcc
              jmt(ih)=nint(1.d0/a(ih)*dlog((rmtout(ih)+b(ih))/b(ih))
     &             +1.d0)
              rmtnewout(ih)=b(ih)*(dexp(a(ih)*(dble(jmt(ih))-1.0d0))-
     &             1.0d0)
C              print*,'RMTNEW=',rmtnewout(ih)
C
              if (ispin.eq.1)  updn='DOWN'
              if (ispin.eq.2)  updn='UP  '
              print*,'iatm=',ih
              call coreshell(ih,i,z,ncore,lcore,ecore)
C     
              write(5,FMT=99002)TABCHSYM(int(z(ih))),title(1:10),updn
              write(5,fmt=9030)rmtout(ih),alat,rmtnewout(ih)
              write(5,fmt=9040)z(ih),rwsout(ih),efnew,vbc
              write(5,fmt=9050)irws(IH),a(IH),b(IH),ncore(i),inew
C     
              IF (NCORE(i).GE.1) THEN
                  DO LC=0,LCMAX
                      DO ICORE=1,NCORE(I)
                          if(lc.eq.LCORE(ICORE,i)) THEN
                              write (5,FMT=9070) LCORE(ICORE,i),
     &                             ECORE(ICORE,i)
                          end if
                      END DO
                  END DO
              END IF
              if (ispin.eq.1) then
                  do IR=1,irws(ih)
                      VM2Z(IR,i) = VRDN(IR,ih)
                  enddo
              else
                  do IR=1,irws(ih)
                      VM2Z(IR,i) = VRUP(IR,ih)
                  enddo
              endif
              write (5,fmt=9051) (VM2Z(IR,i),IR=1,irws(ih))
          enddo
      enddo
C     ADDITIONAL layers for a slab
      if (title(1:17).eq.' LMTO-PERLOV SLAB') then
          do ih=natm-5,1,-1
              do ispin = 1,nspin
                  I = NSPIN* (IH-1) + ISPIN        
C     
                  if (ispin.eq.1)  updn='DOWN'
                  if (ispin.eq.2)  updn='UP  '
                  print*,'ih=',ih
                  print*,'NCORE(',ih,')=',ncore(i)
                  write(5,FMT=99002)TABCHSYM(int(z(ih))),updn
                  write(5,fmt=9030)rmtout(ih),alat,rmtnewout(ih)
                  write(5,fmt=9040)z(ih),rwsout(ih),efnew,vbc
                  write(5,fmt=9050)irws(IH),a(IH),b(IH),ncore(i),inew
C     
                  IF (NCORE(i).GE.1) THEN
                      DO ICORE=1,NCORE(I)
                          write (5,FMT=9070) LCORE(ICORE,i),
     &                         ECORE(ICORE,i)
                      END DO
                  END IF
                  if (ispin.eq.1) then
                      do IR=1,irws(ih)
                          VM2Z(IR,i) = VRDN(IR,ih)
                      enddo
                  else
                      do IR=1,irws(ih)
                          VM2Z(IR,i) = VRUP(IR,ih)
                      enddo
                  endif
                  write (5,fmt=9051) (VM2Z(IR,i),IR=1,irws(ih))
              enddo
          enddo
      endif
C     
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i3,/,2d15.8,/,2i2)
 9070 FORMAT (i5,1p,d20.11)
 9051 FORMAT (1p,4d20.12)
99002 FORMAT (1X,a2,1x,a10,1x,a4)
      END
