      program magdip
      double precision  alat,spinmom(1000,5),orbmom(1000,5),mom(1000,5)
     &     ,br(2,2),R(3),qinp(3,1000),qi(3,1000),magdir(3,5),normfac
     &     ,betrag,nomin(5),sum1(5),sumR(5),M12(1000,1000,5),Msum(5)
     &     ,energy(5),clightryd
      integer i,j,iatom,jatom,natom,itranslat,idir,ndir,xyz
      logical cartesian
C
C     E_d(magdir) = SUM_qi1,qi2 [(mom1*mom2)/c^2 * M12]
C
C     M12 = SUM_br [1/|br+qi1+qi2|^3 *
C                (1-3*[(br+qi1+qi2)*magdir]^2/|br+qi1+qi2|^2)]
C     M12 = SUM_br [sum1]
C     Msum = SUM_qi1,qi2 [M12]
C
C     E_d (energy)-- magnetic dipole-dipole interaction energy [in Rydbergs]
C
C     ndir -- number of magnetisation directions
C     magdir-- magnetisation direction (normalised to 1)
C     br  -- bravais vectors 
C     CARTESIAN -- (F/T) .false. or .true.
C     qinp -- input atomic positions (in CARTESIAN or CRYSTALOGR.)
C     qi -- atomic positions (in CARTESIAN coordinates)
C     mom -- magnetic moment
C
      parameter (clightryd = 2.d0*137.036d0)
      read(5,*) ndir
      do idir=1,ndir
          read(5,*) (magdir(xyz,idir),xyz=1,3)
      enddo
      read(5,*) alat
      read(5,*)
      read(5,*) br(1,1),br(2,1)
      read(5,*) br(1,2),br(2,2)
      do j=1,2
          do i=1,2
              br(i,j)=br(i,j)*alat
          enddo
          print*,(br(i,j),i=1,2)
      enddo
      read(5,*)
      read(5,*) natom
      read(5,*)
      cartesian=.false.  ! default means crystallograpic atom. positions
      read(5,*) cartesian
      read(5,*)
      do i = 1,natom
          read(5,*) qinp(1,i),qinp(2,i),qinp(3,i),
     &         (spinmom(i,idir),orbmom(i,idir),idir=1,ndir)
          qi(3,i)=qinp(3,i)*alat
          do idir=1,ndir
              mom(i,idir)=spinmom(i,idir)+orbmom(i,idir)
         enddo
      enddo
      read(5,'(i8)')itranslat
      print*,'itranslat=',itranslat
      IF (.NOT.CARTESIAN) then  ! Rescale lattice
          print*,'cartesian false'
          print*,'atom. pos. will be recalculated to the CARTESIAN'
C     qi_vec =  qinp_vec * sum_i [br_vec_i]
          DO i=1,natom
              DO j=1,2
                  qi(j,i)= qinp(1,i)*br(j,1)+qinp(2,i)*br(j,2)
              enddo
          enddo
      else 
          print*,'cartesian true'
          do i = 1,natom
              do j=1,3
                  qi(j,i)= qinp(j,i)*alat
              enddo
          enddo
      endif
C     initializing
      do idir=1,ndir
          Msum(idir)=0.0d0
      enddo
CC
CC     normalise magdir
CC
      print*,'  '
      do idir=1,ndir
          print*,'magdir(',idir,')= ',(magdir(xyz,idir),xyz=1,3)
          normfac=(magdir(1,idir)**2.d0+magdir(2,idir)**2.d0+
     &         magdir(3,idir)**2.d0)**0.5d0
          print*,'norm factor= ',normfac
          do xyz=1,3
              magdir(xyz,idir)=magdir(xyz,idir)/normfac
              print*,idir,magdir(xyz,idir)
          enddo
      enddo
      print*,'  '
CC
CC     Msum = SUM_qi1,qi2 [M_q1,q2] = M11 + M12 + M21 + M22 
CC     M12  = sumR * mom1 * mom2
CC     sumR = SUM_br (1/betrag^(3/2) * (1-3*nomin/betrag))
CC
      do iatom=1,natom
          do jatom=iatom,natom
              do idir=1,ndir
                  sumR(idir)=0.0d0
              enddo
              print*,'   '
              print*,'iatom=',iatom,' jatom=',jatom
              print*,qi(1,iatom),qi(2,iatom),qi(3,iatom)
              print*,qi(1,jatom),qi(2,jatom),qi(3,jatom)
C     to get good results, i=(-100,100) and j=(-100,100) is already enough
      do i=-itranslat,itranslat
          do j=-itranslat,itranslat
              if (.NOT.(i.eq.0.and.j.eq.0.and.iatom.eq.jatom)) then
                  R(1)=br(1,1)*i+br(1,2)*j
                  R(2)=br(2,1)*i+br(2,2)*j
                  R(3)=0.d0
C     
                  betrag=(R(1)+qi(1,iatom)-qi(1,jatom))**2.d0+
     &                 (R(2)+qi(2,iatom)-qi(2,jatom))**2.d0+
     &                 (R(3)+qi(3,iatom)-qi(3,jatom))**2.d0
                  do idir=1,ndir
                      nomin(idir)=
     &                   ((R(1)+qi(1,iatom)-qi(1,jatom))*magdir(1,idir)
     &                   +(R(2)+qi(2,iatom)-qi(2,jatom))*magdir(2,idir)
     &                   +(R(3)+qi(3,iatom)-qi(3,jatom))*magdir(3,idir)
     &                     )**2.d0
                      sum1(idir)=1.d0/betrag**(3.d0/2.d0)*
     &                     (1.d0-3.d0*nomin(idir)/betrag)
                      sumR(idir)=sumR(idir)+sum1(idir)
                  enddo
CC
CC     this 'if' is for the output to see only interesting parts
C                  if (i.ge.0.and.j.ge.0.and.i.le.10.and.j.le.10.or.
C     &                 i.ge.0.and.i.le.10.and.j.ge.0.and.mod(j,100).eq.0
C     &                 .or.i.ge.0.and.j.ge.0.and.mod(i,100).eq.0.and
C     &                 .mod(j,100).eq.0) then
C                      print*,i,j,'betrag=',betrag,' nomin=',nomin
C     &                     ,' sum1=',sum1
C                  endif
              endif
          enddo 
      enddo
      do idir=1,ndir
          print*,'sumR(',idir,')=',sumR(idir)
          M12(iatom,jatom,idir)=sumR(idir)*mom(iatom,idir)*mom(jatom
     &         ,idir)
      if(iatom.ne.jatom) then 
          Msum(idir)=Msum(idir)+M12(iatom,jatom,idir)*2.d0
          print*,'  M12(',iatom,jatom,idir,') = M12(',jatom,iatom,idir
     &         ,')=',M12(iatom,jatom,idir),'  Msum(',idir,')=',Msum(idir
     &         )
      else
          Msum(idir)=Msum(idir)+M12(iatom,jatom,idir)
          print*,'  M12(',iatom,jatom,idir,') =',M12(iatom,jatom,idir)
     &         ,'  Msum(',idir,')=',Msum(idir)
      endif
      enddo
      enddo
      enddo
      do idir=1,ndir
          energy(idir)=Msum(idir)/(clightryd**2.d0)
          print*,'energy(',(magdir(xyz,idir),xyz=1,3),')=',energy(idir)
      enddo
      print*,'   '
      do idir=1,ndir
          do jdir=idir+1,ndir
              print*,idir,jdir,'energy(',(magdir(xyz,idir),xyz=1,3)
     &             ,')-energy(',(magdir(xyz,jdir),xyz=1,3),')=',
     &             energy(idir)-energy(jdir)
          enddo
      enddo
C     
      end
