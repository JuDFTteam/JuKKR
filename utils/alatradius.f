      program alatradius
      implicit none
C-------------------------------------------------------------
C     alatradius $old.pot $new.pot
C-------------------------------------------------------------
C used if one wants to estimate a new potential 
C for another lattice constant
C also the ratio between spheres can be changed
C-------------------------------------------------------------
C giving a lattice parameter ALAT 
C        the bravais lattice (bcc or fcc),
C    and the RATIOR between radii
C the program calculates the Wigner-Seitz radii RWS 
C and the Muffin-tin radii RMT
C  (program keeps the number IRWS unchanged)
C
C RMT lies between two points of the radial mesh
C program puts RMT to the nearer point of the radial mesh 
C and recalculating the value of RMT --> RMTNEW
C
C it also supplies the values of B such that the radial
C mesh R(I) is defined by
C
C       R(I) = B1* ( DEXP(A1*(I-1)) - 1 ) with A1 read from the potential
C
C structures: bcc, fcc, CsCl, CuAu, Cu3Au, ZnS
C !!! CuAu in the unit cell form (2 atoms, tetragonal)
C     c/a you have to specify
C
C     for ZnS specify the number of layers Vc, Fe, GaAs
C     Number of second half Vc2 = natm-vc-fe-asga-fe
C
      double precision alat1,alat2,COA,rwsinp(200),rmtinp(200)
     &     ,rwsout(200),rmtout(200),rmtnewinp(200),rmtnewout(200),z(200)
     &     ,efnew,vbc(200),vm2z(1000,200),ecore(20,200),pi,a1(200)
     &     ,b1(200),a2(200),b2(200),ratior,rwsideal,rmtideal,vshift(200)
      double precision volume
C      double precision vshift(1000,40),radiusinp(1000,40)
C     &     ,radiusout(1000,40)
      integer irws(200),inew,lcore(20,200),icore,ncore(200),ir
      integer natm,jmt(200),i,ih,nspin,ispin
      integer vc,fe,asga
      character brav*7
      parameter (pi=3.141592653589793238462643D0)
      character*256 infile,outfile
      character*64 title(200)
C ------------------------------------------------------------
      ratior=1.0d0
      COA=1.0d0
      volume=0.d0
      nspin=2
C     ----------INITIAL DATA--------------------------
      write(6,*) ' number of spins '
      read(5,*) nspin
      write(6,*) ' give the number of atoms: '
      read(5,*) natm
C     ----------Bravais lattice---------------------
 101  write(6,*) ' give the bravais lattice '
      if (natm.eq.1) 
     &     write (6,*) '(for natm=1, bcc and fcc implemeted)'
      if (natm.eq.2) 
     &     write (6,*) '(for natm=2, CsCl and CuAu implemeted)'
      if (natm.eq.4) 
     &     write (6,*) '(for natm=4, Cu3Au and ZnS implemeted)'
      if (natm.ge.10)
     &     write (6,*) '(only for Fe-GaAs systems)'
      read(5,*) brav
      if ((brav(1:3).ne.'bcc').and.(brav(1:3).ne.'fcc')
     &     .and.(brav(1:4).ne.'CsCl').and.(brav(1:4).ne.'CuAu')
     &     .and.(brav(1:5).ne.'Cu3Au').and.(brav(1:3).ne.'ZnS')
     &     .and.(brav(1:7).ne.'Fe-GaAs')) then 
          write(*,*)
          write(*,*)' only bcc, fcc, CsCl, CuAu, Cu3Au and ZnS
     &         and Fe/GaAs implemented so far' 
          write(*,*)
          goto 101
      end if
C     -----------------------------------------------
      if (brav(1:7).eq.'Fe-GaAs') then
          write(6,*) ' what is the number of Vc, Fe, AsGa layers'
          read(5,*) vc,fe,asga
          write(6,*) ' what is the ratio RWS(Ga)/RWS(As) '
          write(6,*)
     &     'if you want to keep RATIO according to the input, write "0"'
          read(5,*) ratior
      endif
      if ((natm.gt.1).and.(natm.lt.10)) then 
          write(6,*) ' what is the ratio RWS(1)/RWS(2) '
          write(6,*)
     &     'if you want to keep RATIO according to the input, write "0"'
          read(5,*) ratior
      endif
      write(6,*) ' give new lattice constant: '
      write(6,*) ' if you dont want to change it, write "0"'
      read(5,*) alat2
      if (brav(1:4).eq.'CuAu') then
          COA=1.41421356D0
          write(6,*) ' give the c/a '
          write(6,*) ' if you want c/a to be SQRT(2), write "0"'
          read(5,*) COA
          if (COA.eq.0) COA=1.41421356D0
      endif
C     ------OPEN files --------------
      call getarg(1,infile)
      open(unit=3,file=infile)
      call getarg(2,outfile)
      open(unit=4,file=outfile)      
C     CYCLE - natm * 2spins
C     ih means which atom is considered
C     i  means ih spin-resolved
      do ih = 1,natm
          do ispin = 1,nspin
              I = NSPIN* (IH-1) + ISPIN
C     --------------INPUT------------
              read(3,'(a)')title(i)
C              print *,'ih=',ih,' i=',i
              print *,title(i)
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
          enddo
      enddo
      close(3)
C     --------------CHANGE------------
      if ((natm.lt.10).and.(ratior.eq.0)) then 
          ratior=rwsinp(1)/rwsinp(2)
          print *,'rwsinp(1)=',rwsinp(1)
          print *,'rwsinp(2)=',rwsinp(2)
          print *,'ratior=',ratior
      elseif ((brav.eq.'Fe-GaAs').and.(ratior.eq.0)) then 
          ratior=rwsinp(2*vc+2*fe+4)/rwsinp(2*vc+2*fe+2)
          print *,'rwsinp(Ga)=',rwsinp(2*vc+2*fe+4)
          print *,'rwsinp(As)=',rwsinp(2*vc+2*fe+2)
          print *,'ratior=',ratior
      endif
      if (alat2.eq.0) alat2=alat1
      print *,'alat1=',alat1
      print *,'alat2=',alat2
C     -------------RMT and RWS--------------------------------------
C     ------------IDEAL --> RWS is calculated for RATIO = 1 ---------
C   
      if (brav.eq.'bcc'.or.brav.eq.'CsCl') then
          rwsideal=alat2/2.d0*(3.d0/pi)**(1.d0/3.d0)
          rmtideal=alat2*sqrt(3.d0)/4.d0
      elseif (brav.eq.'fcc'.or.brav.eq.'Cu3Au') then
          rwsideal=alat2/2.d0*(3.d0/2.d0/pi)**(1.d0/3.d0)
          rmtideal=alat2*sqrt(2.d0)/4.d0
      elseif (brav.eq.'CuAu') then
          rwsideal=alat2/2.d0*(3.d0/pi*COA)**(1.d0/3.d0)
C        c/a = SQRT(2) -- the unit cell from CuAu
C     if c/a > SQRT(2) atoms Cu-Cu nearest neighbors 
C     if c/a < SQRT(2) atoms Cu-Au nearest neighbors 
          if (COA.ge.1.41421356D0) then
              rmtideal=alat2*0.5d0
          else 
              rmtideal=alat2/4.d0*(2.d0+COA**2.d0)**0.5d0
          endif
      elseif ((brav.eq.'ZnS').or.(brav.eq.'Fe-GaAs')) then
          rwsideal = alat2/4.d0*(3.d0/pi)**(1.d0/3.d0)
          rmtideal = alat2/8.d0*sqrt(3.d0)
      endif
      print *,'rwsideal=',rwsideal
      print *,'rmtideal=',rmtideal
C      
      if (natm.eq.1) then 
          rwsout(1)=rwsideal
          rmtout(1)=rmtideal
      elseif (brav.eq.'CsCl'.or.brav.eq.'CuAu') then
          rmtout(2)=rmtideal*2.d0/(1.d0+ratior)
          rmtout(1)=rmtout(2)*ratior
          rwsout(2)=rwsideal*(2.d0/(1.d0+ratior**3.d0))**(1.d0/3.d0)
          rwsout(1)=rwsout(2)*ratior
      elseif (brav.eq.'Cu3Au') then
          rmtout(2)=rmtideal*2.d0/(1.d0+ratior)
          rmtout(1)=rmtout(2)*ratior
          rwsout(2)=rwsideal*(4.d0/(3.d0+ratior**3.d0))**(1.d0/3.d0)
          rwsout(1)=rwsout(2)*ratior
          rmtout(3)=rmtout(2)
          rmtout(4)=rmtout(2)
          rwsout(3)=rwsout(2)
          rwsout(4)=rwsout(2)
      elseif (brav.eq.'ZnS') then
          rmtout(2)=rmtideal*2.d0/(1.d0+ratior)
          rmtout(1)=rmtout(2)*ratior
          rwsout(2)=rwsideal*(2.d0/(1.d0+ratior**3.d0))**(1.d0/3.d0)
          rwsout(1)=rwsout(2)*ratior
          rmtout(3)=rmtout(2)
          rmtout(4)=rmtout(1)
          rwsout(3)=rwsout(2)
          rwsout(4)=rwsout(1)
      elseif (brav.eq.'Fe-GaAs') then
          do i=1,2*vc+2*fe
C              print *,'i=',i
              rwsout(i)=rwsideal
              rmtout(i)=rmtideal
          enddo
          do i=2*vc+2*fe,2*vc+2*fe+2*asga,4
C             
C     better to put always RATIOR = 1
C
              rwsout(i+2)=rwsideal*(2.d0/(1.d0+ratior**3.d0))
     &             **(1.d0/3.d0)
              rwsout(i+1)=rwsout(i+2)
              rwsout(i+3)=rwsout(i+2)*ratior
              rwsout(i+4)=rwsout(i+2)*ratior
C     
              rmtout(i+2)=rmtideal*2.d0/(1.d0+ratior)
              rmtout(i+1)=rmtout(i+2)
              rmtout(i+3)=rmtout(i+2)*ratior
              rmtout(i+4)=rmtout(i+2)*ratior
C     
              print *,i+1,i+4
          enddo
          do i=2*vc+2*fe+2*asga+1,natm
C              print *,'i=',i
              rwsout(i)=rwsideal
              rmtout(i)=rmtideal
          enddo
      endif
C
      print *,'brav=',brav
      do ih = 1,natm
          print *,'rmtinp(',ih,')=',rmtinp(ih)
          print *,'rmtout(',ih,')=',rmtout(ih)
          print *,'rwsinp(',ih,')=',rwsinp(ih)
          print *,'rwsout(',ih,')=',rwsout(ih)
      enddo
C     --------CHECK RADII-----------------------------
      if (brav.eq.'Fe-GaAs') then
          print*,'volume of the unit cell=', 
     &         alat2**3.d0/2.d0/4.d0*dble(natm/2)
      else 
          print*,'volume of the cube=',alat2**3.d0*COA
      endif
C
      if (brav.eq.'bcc') then
          print*,'volume radii=',8.d0*pi/3.d0*rwsout(1)**3.d0
      elseif (brav.eq.'fcc') then
          print*,'volume radii=',16.d0*pi/3.d0*rwsout(1)**3.d0
      elseif (brav.eq.'CsCl'.or.brav.eq.'CuAu'.or.brav.eq.'Cu3Au') then
          print*,'volume radii=',4.d0/3.d0*pi*(rwsout(1)**3.d0+rwsout(2)
     &         **3.d0)
      elseif (brav.eq.'ZnS') then  
          print*,'volume radii=',4.d0/3.d0*pi*8.d0*(rwsout(1)**3.d0
     &         +rwsout(2)**3.d0)
      elseif (brav.eq.'Fe-GaAs') then
          do i=1,natm
              volume=volume+4.d0/3.d0*pi*rwsout(i)**3.d0
          enddo
          print*,'volume radii=',volume
      endif
C     --------SECOND CYCLE------------------------------------
      do ih = 1,natm
          do ispin = 1,nspin
              I = NSPIN* (IH-1) + ISPIN        
              a2(ih) = a1(ih)
              b2(ih) = rwsout(ih)/(dexp(a2(ih)*dble(irws(ih)-1.d0))-1.d0
     &             )
C     calculates the position of RMT in the radial mesh
              jmt(ih)=nint(1.d0/a2(ih)*dlog((rmtout(ih)+b2(ih))/b2(ih))
     &             +1.d0)
C              print *, jmt(ih), irws(ih)
              rmtnewout(ih)=b2(ih)*(dexp(a2(ih)*dble(jmt(ih)-1.d0))-1.d0
     &             )
C     -------------------------------------------------------------------
C     shift the new potential at the WS-radius
C              vshift(ih) = 2.d0*z(ih)*(rwsinp(ih)-rwsout(ih))/
C     &             (rwsinp(ih)*rwsout(ih))
C     alternatively, shift the new potential at the MT-radius
              vshift(ih) = 2.d0*z(ih)*(rmtinp(ih)-rmtout(ih))
     &             /(rmtinp(ih)*rmtout(ih))
              print *,'z=',z(ih),' vshift=',vshift(ih)
              do ir=1,irws(ih)
                  vm2z(ir,i)=vm2z(ir,i)+vshift(ih)
              enddo
C     ----------- RECALCULATE VM2Z for each mesh point -------------
C              print *,'z=',z(ih)
C              do ir=20,irws(ih)
C                  radiusinp(ir,ih)=b1(ih)*(dexp(a1(ih)*dble(ir-1.d0))-1
C     &                 .d0)
C                  radiusout(ir,ih)=b2(ih)*(dexp(a2(ih)*dble(ir-1.d0))-1
C     &                 .d0)
C                  vshift(ir,ih)=(radiusinp(ir,ih)-radiusout(ir,ih))
C     &                 /radiusinp(ir,ih)/radiusout(ir,ih)
C                  vm2z(ir,i)=vm2z(ir,i)+vshift(ir,ih)
C              enddo
C     --------------OUTPUT------------
              write(4,'(a)')title(i)
              write(4,fmt=9030)rmtout(ih),alat2,rmtnewout(ih)
              write(4,fmt=9040)z(ih),rwsout(ih),efnew,vbc(ispin)
              write(4,fmt=9050)irws(IH),a2(IH),b2(IH),ncore(i),inew
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
      close(4)
C  19   format(i3,/,i3)
 9030 FORMAT (3f12.8)
 9040 FORMAT (f10.5,/,f10.5,2f15.10)
 9050 FORMAT (i3,/,2d15.8,/,2i2)
 9070 FORMAT (i5,1p,d20.11)
 9051 FORMAT (1p,4d20.12)
      end
