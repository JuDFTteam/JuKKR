      program potfit4tbkkr
      double precision RIN(1000,100),RHOUP(1000,100),VUP(1000,100)
     &     ,RHODN(1000,100),VDN(1000,100),QMPHI(20),QMTET(20)
      integer jws(100),jwsnewall,jwsnew(100),NMAX
      double precision yp1UP,yp1DN,ypnUP,ypnDN,rout(100),voutUP(100)
     &     ,voutDN(100),rws(100),rmt(100),aall,a(100),b(100),
     &     vrup(1000,100),vrdn(1000,100)
      character*256 infile,outfile
      parameter (NMAX=1000)
      double precision y2UP(NMAX),y2DN(NMAX)
C     
      double precision ALAT,BOA,COA,EFERMI,BRX(3),BRY(3),BRZ(3),QX(100)
     &     ,QY(100),QZ(100)
      integer IT,IR,NT,NQ,BRAVAIS,II
      character*80 title
C
C     ----------------------------------SPRKKR
C
      CHARACTER*80 LINE,system
      CHARACTER*10 SCFALG,SCFVXC,orbpol
      CHARACTER*12 RMESHTYP
      CHARACTER*4 STR4,txtt(20)
      logical scfstart,semicore,BREITINT,extfield,BLCOUPL,QAVAILABLE
      integer fmtvers,ichar,BRAVAISLAT(14)
      integer im,nm,ns,netab,irel,ibzint,nptbz,itrscf,KMROT
     &     ,ibas,idum,jmt(20),iqin,ITOQ(20,20),io,noq(20),jq,itp,nat(20)
     &     ,imt(20),iqat(20,20),ia
      double precision SCFMIX,SCFTOL,QMVEC(3),BEXT,dx(20)
     &     ,conc(20),smad(20,20),VT(750,20),BT(750,20)
      double precision z(100)
      logical yes
      integer i,j
C
      DATA BRAVAISLAT/12,13,14,11,8,9,10,4,5,6,7,2,3,1/
C
C
      external findreals,spline,splint,derivat,rmtcalc,tbkkr_pformat
C     
C     TBKKR   R(J) = B *  [exp(A(j-1)) -1]    
C
C     SPRKKR  R(J) = R(1)* exp(DX*(j-1))
C
C     LTMO quadratic mesh R(J) ~ B * A^2 
C
C     this program reads Perlov's and SPRKKR potential, 
C     inter(extra)polates to a
C     potential usable for TBKKR and writes 'misha' for plotting
C
C     !!!!!!!!! when using Perlov potential, in lmt.pot 
C     set the ALAT and all RWS up to 8 decimal numbers
C
      call getarg(1,infile)
      open(unit=4,file=infile)
C     ------------------------------
      READ (4,'(A)') LINE
C     --------------------SPRKKR
      IF ( LINE(1:34).EQ.'SPRKKR    selfconsistent potential' .OR. 
     &     LINE(1:24).EQ.'SPRKKR    SCF-start data' ) THEN
          IF ( LINE(11:19).EQ.'SCF-start' ) SCFSTART = .TRUE.
          READ (4,99003) TITLE
          IF ( TITLE(1:7).NE.'VERSION' ) THEN
              FMTVERS = 0
          ELSE
              FMTVERS = ICHAR(TITLE(10:10)) - ICHAR('1') + 1
              READ (4,99003) TITLE
          END IF
          READ (4,99003) SYSTEM
          READ (4,99006) NM       !number of r-meshes
          READ (4,99006) NQ       !number of atomic positions, NQ<NT =CPA
          READ (4,99006) NT       !number of atomic species
          READ (4,99006) NS
          READ (4,99006) NETAB
          READ (4,99006) IREL
          READ (4,99006) IBZINT
          READ (4,99006) NPTBZ
          READ (4,99005) STR4
          READ (4,99003) SCFVXC
          READ (4,99003) SCFALG
          READ (4,99006) ITRSCF
          READ (4,99008) SCFMIX
          READ (4,99008) SCFTOL
          IF ( FMTVERS.GE.4 ) READ (4,99009) SEMICORE
          READ (4,99009) BREITINT
          READ (4,99003) ORBPOL
          IF ( ORBPOL.EQ.'          ' ) ORBPOL = 'NONE'
          IF ( FMTVERS.GE.1 ) THEN
              READ (4,99009) EXTFIELD
              READ (4,99009) BLCOUPL
              READ (4,99008) BEXT
          END IF
          IF ( FMTVERS.GE.2 ) THEN
              READ (4,99006) KMROT
              READ (4,99008) (QMVEC(I),I=1,3)
          END IF
          READ (4,99008) EFERMI
          READ (4,99007) STR4,BRAVAIS
C-----------------------------------------------allow for old input LAT
          IF ( STR4.EQ.'LAT ' ) BRAVAIS = BRAVAISLAT(BRAVAIS)
          READ (4,99006) IBAS
          READ (4,99008) ALAT
          READ (4,99008) BOA
          READ (4,99008) COA
          call findreals(boa,1,yes)
          call findreals(coa,1,yes)
          print*,'alat=',alat,' boa=',boa,' coa=',coa
          IF ( FMTVERS.GE.3 ) THEN
C     
              print*,'translation vectors'
              DO I = 1,3
                  READ (4,99008) BRX(I),BRY(I),BRZ(I)
                  call findreals(brx(i),1,yes)
                  call findreals(bry(i),1,yes)
                  call findreals(brz(i),1,yes)
                  write(*,9420)BRX(I),BRY(I),BRZ(I)
              END DO
          END IF
C     ---------------------- read in radial mesh parameters
          READ (4,99003) STR4
          READ (4,99003) STR4
          READ (4,99003) RMESHTYP
          DO IM = 1,NM
              READ (4,99006) IDUM
              READ (4,99008) RWS(IM)
              READ (4,99008) RMT(IM)
              READ (4,99006) JWS(IM)
              READ (4,99006) JMT(IM)
              IF ( RMESHTYP.EQ.'EXPONENTIAL ' ) THEN
                  READ (4,99008) rin(1,IM)
                  READ (4,99008) DX(IM)
              ELSE
                  STOP ' check RMESHTYP in POTFILE '
              END IF
          END DO
C     
          IF ( FMTVERS.GE.1 ) THEN
              READ (4,'(/)')
              IF ( FMTVERS.EQ.1 ) THEN
                  DO IQ = 1,NQ
                      READ (4,99012) IQIN,QX(IQ),QY(IQ),QZ(IQ),NOQ(IQ),
     &                     (ITOQ(IO,IQ),CONC(ITOQ(IO,IQ)),IO=1,
     &                     NOQ(IQ))
                  END DO
              ELSE
                  DO IQ = 1,NQ
                      READ (4,99013) IQIN,QX(IQ),QY(IQ),QZ(IQ),NOQ(IQ),
     &                     (ITOQ(IO,IQ),CONC(ITOQ(IO,IQ)),IO=1,
     &                     NOQ(IQ))
                  END DO
                  print*,'atomic positions'
                  do IQ = 1,NQ
                      call findreals(qx(iq),1,yes)
                      call findreals(qy(iq),1,yes)
                      call findreals(qz(iq),1,yes)
                      write(*,9420)qx(iq),qy(iq),qz(iq)
                  enddo
                  DO IQ = 1,NQ
                      READ (4,99014) IQIN,QMTET(IQ),QMPHI(IQ)
                  END DO
              END IF
              QAVAILABLE = .TRUE.
C     --------------------------- MADELUNG MATRIX SMAD
              READ (4,'(/)')
              DO IQ = 1,NQ
                  READ (4,99010) (SMAD(IQ,JQ),JQ=1,NQ)
              END DO
C     
          END IF
C     
C     -------------------- read in atomic type information
          DO IT = 1,NT
              print*,'it= ',it,'nt= ',nt
              READ (4,'(/,10X,I10,/,A4)') ITP,TXTT(IT)
              READ (4,9001) Z(IT)
              READ (4,99006) NAT(IT)
              READ (4,99008) CONC(IT)
              READ (4,99006) IMT(IT)      !which mesh type
              READ (4,99006) (IQAT(IA,IT),IA=1,NAT(IT))
              IF ( FMTVERS.GE.4 ) READ (4,99006)
              IM = IMT(IT)
              IF ( SCFSTART ) THEN
                  DO I = 1,JWS(IM)
                      VT(I,IT) = 0.D0
                      BT(I,IT) = 0.D0
                  END DO
              ELSE
                  READ (4,99010) (VT(I,IT),I=1,JWS(IM))
                  READ (4,99010) (BT(I,IT),I=1,JWS(IM))
              END IF
              print*,'it= ',it,'nt= ',nt
          END DO
          WRITE (6,99003) ' SPRKKR input dataset read '
CCCCCCCCCCC
C     
          if (nm.lt.nt) then
              do i=1,nt-1
                  do j=i+1,nt
                      print*,'i=',i,'j=',j
                      if (imt(i).eq.imt(j)) then
                          rws(j)=rws(i)
                          rmt(j)=rmt(i)
                          jws(j)=jws(i)
                          jmt(j)=jmt(i)
                          rin(1,j)=rin(1,i)
                          dx(j)=dx(i)
                          print*,'rmt(',j,')=',rmt(j)
                      endif
                  enddo
              enddo
          endif
          DO IT = 1,NT
              DO I = 1,JWS(IM)
C-----------------exponential mesh SPRKKR 
                  rin(i,it)=rin(1,it)*dexp(dx(it)*(i-1))
C                  print*,rin(i,it),VT(I,IT),BT(I,IT)
                  vUP(i,it)= VT(I,IT) + BT(I,IT)
                  vUP(i,it)=  vUP(i,it) + 2.D0*Z(IT)/rin(i,it)
                  vDN(i,it)= VT(I,IT) - BT(I,IT)
                  vDN(i,it)= vDN(i,it) + 2.D0*Z(IT)/rin(i,it)
C                  print*,rin(i,it),vUP(I,IT),vDN(I,IT)
              enddo
              print*,'VUP(jws',it,')=',VUP(jws(it),it),'  VDN(jws',it
     &             ,')=',VDN(jws(it),it)
              
          END DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     ---------------PERLOV LMTO
      elseif ( (LINE(1:12).EQ.' LMTO-PERLOV') .OR. 
     &         (LINE(1:12).EQ.' lmto-perlov') ) THEN
          print*,'!!!!!! check the ALAT and all RWS'
          print*,'complete up to 8 decimal numbers!!!!!!!'
          print*,'    '
          READ (4,*) NT
          READ (4,*) NQ
          READ (4,*) BRAVAIS
          READ (4,*) jws
          READ (4,*) EFERMI
          READ (4,*) ALAT,BOA,COA
          call findreals(boa,1,yes)
          call findreals(coa,1,yes)
          print*,'alat=',alat,' boa=',boa,' coa=',coa
C     
          print*,'translation vectors'
          DO IQ = 1,3
              READ (4,*) BRX(IQ),BRY(IQ),BRZ(IQ),II
              call findreals(brx(iq),1,yes)
              call findreals(bry(iq),1,yes)
              call findreals(brz(iq),1,yes)
              if (abs(brx(iq)).lt.1.d-15) brx(iq)=0.d0
              if (abs(bry(iq)).lt.1.d-15) bry(iq)=0.d0
              if (abs(brz(iq)).lt.1.d-15) brz(iq)=0.d0
              write(*,9420)BRX(iq),BRY(iq),BRZ(iq)
          enddo
          print*,'atomic positions'
          DO IQ = 1,NQ
              READ (4,*) QX(IQ),QY(IQ),QZ(IQ),II
              call findreals(qx(iq),1,yes)
              call findreals(qy(iq),1,yes)
              call findreals(qz(iq),1,yes)
              if (abs(qx(iq)).lt.1.d-15) qx(iq)=0.d0
              if (abs(qy(iq)).lt.1.d-15) qy(iq)=0.d0
              if (abs(qz(iq)).lt.1.d-15) qz(iq)=0.d0
              write(*,9420)qx(iq),qy(iq),qz(iq)
          enddo
          READ (4,*) KMROT
          IF ( KMROT.NE.0 ) THEN
              DO IQ = 1,NQ
                  READ (4,*) QMPHI(IQ),QMTET(IQ)
              END DO
          END IF
          DO IT = 1,NT
              READ (4,*) II,RWS(II),z(it)
              jws(IT) = jws(1)
              RMT(it) = 0.0D0
C     reading the radius, chargeUP, VUP, chargeDN, VDN
              DO IR = 1,jws(it)
                  READ (4,*) RIN(IR,it),RHOUP(IR,it),VUP(IR,it),RHODN(IR
     &                 ,it),VDN(IR,it)
              END DO
              print*,'VUP(jws',it,')=',VUP(jws(it),it),'  VDN(jws',it
     &             ,')=',VDN(jws(it),it)
C     seting the imt(it) variable for rmtcalc subroutine
              imt(it)=it
          END DO
          WRITE (6,99003) ' PERLOV input dataset read '
      endif
C     
      close(4)
C     
      title=line
      print*,'title=',title
C     ------------------------------------------------------------
C     INTERPOLATION
C     using spline and splint interpolation from Numerical Recipes
      jwsnewall=484
      aall=0.025d0
      ball=-9999.0
      print*,'how many mesh points do you want (353)?'
      read(*,*)jwsnewall
      print*,'what A for b*[exp(A(ir-1)-1] (0.025d0)?'
      read(*,*)aall
      print*,'what B for b*[exp(A(ir-1)-1] (-9999.0)?'
      print*,' IF YOU ARE USING ASA POTENTIAL THEN USE (-9999.0)?'
      read(*,*) ball
      
      do it=1,nt
          jwsnew(it)=jwsnewall
          a(it)=aall
          b(it)=ball
      enddo
C
      open(unit=10,file='misha')
      do it=1,nt
          write(10,*) "################## it=",it
          rout(it)=0.0d0
          yp1UP=0.0d0
          yp1DN=0.0d0
          ypnUP=0.0d0
          ypnDN=0.0d0
          do i=1,1000
              y2UP(i)=0.0d0
              y2DN(i)=0.0d0
          end do
C          print*,'calculating the first derivative for boundaries'
          call derivat(rin(1,it),vUP(1,it),jws,yp1UP,ypnUP)
          call derivat(rin(1,it),vDN(1,it),jws,yp1DN,ypnDN)
C          print*,'calculating second derivatives y2'
          call spline(rin(1,it),vUP(1,it),jws,yp1UP,ypnUP,y2UP)
          call spline(rin(1,it),vDN(1,it),jws,yp1DN,ypnDN,y2DN)
C     
C     set the number of r-mesh points
          IF(b(it).eq.-9999.0) then
           b(it) = rws(it)/(dexp(a(it)*(dble(jwsnew(it))-1.0d0))-1.0d0)
          end if
          print*,'z(',it,')=',z(it)
C          print*,'jwsnew(',it,')=',jwsnew(it),' a(',it,')=',a(it),' b('
C     &         ,it,')=',b(it),' rws(',it,')=',rws(it)
          do IR=1,jwsnew(it)
              rout(it)=b(it)*(dexp(a(it)*(dble(IR)-1.0d0))-1.0d0)
              call splint(rin(1,it),vUP(1,it),y2UP,jws,rout(it)
     &             ,voutUP(it))
              call splint(rin(1,it),vDN(1,it),y2DN,jws,rout(it)
     &             ,voutDN(it))
              write (10,fmt=9051)rout(it),voutUP(it),voutDN(it)
C     VRUP(IR,IT) important to store the potential 
C     for the whole r-mesh
              vrup(ir,it)=voutUP(it)
              vrdn(ir,it)=voutDN(it)
          enddo
C          print*,'rout(',it,')=',rout(it),' voutUP(',jwsnewall,',it,')='
C     &         ,voutUP(it),' voutDN(',jwsnewall,',it,')=',voutDN(it)
C          print*,' '
          write(10,*) "          "
      enddo
      close(10)
C     ------------ TBKKR potential format
C     
      print*,'nm=',nm,'  nq=',nq,'  nt=',nt
      do it = 1,nt
          print*,'rmt(',it,')=',rmt(it)        
      enddo
      call rmtcalc(nt,nq,nm,alat,rws,brx,bry,brz,qx,qy,qz,rmt,imt)
      do it = 1,nt
          print*,'rmtcalc(',it,')=',rmt(it)        
      enddo
C     
      call getarg(2,outfile)
      open(unit=5,file=outfile)      
      call tbkkr_pformat(nt,rmt,alat,z,rws,efermi,jwsnew,a,b,vrup
     &     ,vrdn,title)
      close(5)
C     
C     
C
 9001 format (10x,d15.0)
 9051 FORMAT (3d20.12)
 9420 format (f12.8,2x,f12.8,2x,f12.8)
99003 FORMAT (10X,A)  
99005 FORMAT (10X,A,I5)
99006 FORMAT (10X,20I10)
99007 FORMAT (A4,6X,20I10)
99008 FORMAT (10X,5F20.10)
99009 FORMAT (10X,L3)
99010 FORMAT (5E16.9)
99012 FORMAT (I4,4X,2(F10.6,1X),F10.6,8X,I3,10X,20(I3,F6.3))
99013 FORMAT (I4,4X,2(F14.10,1X),F14.10,8X,I3,10X,20(I3,F6.3))
99014 FORMAT (I4,9X,F15.10,9X,F15.10)
      
      end
      
