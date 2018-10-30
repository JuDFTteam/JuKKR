      subroutine pdbplot(rbasis,natyp,rr,nr,bbox,zatom,alat)
      implicit none
c#@# KKRtags: VORONOI visualization
      integer nbig
      parameter (nbig=10000)
      include 'inc.geometry'
c
c
c
      integer ntot,natyp,nr,kaoez
      real*8 RBASIS(3,*),RR(3,0:NRD),bbox(3),zatom(*),alat
c
      integer jatom,irevsort(nbig),isort(nbig),natoms
      integer pos,nnn,ixx,nplot,ia,i,itot,ib,i0,n,iatom
      integer ide(nbig)
      real*8 rg(nbig),rsort(nbig),atomx(nbig),atomy(nbig),
     *       atomz(nbig)
      real*8 one,ocuup,temp,plotx(nbig),ploty(nbig),plotz(nbig)
      real*8 plotxa,plotya,plotza,color(3),radious,z,y,x,test,
     &  test1,occup
      character*1 xx1
      character*2 xx2
      character*3 xx3,elem_name
      character*4 xx4,name(nbig)
      character*5 name5
      character*24 name24
      logical whichatoms(0:113)
c--------------------------------------------------------
      open(27,FILE='lattice.pdb',status='unknown')
      xx1=' '
      xx2='  '
      xx3='   '
      xx4='    '
      one = 1.0d0
      occup = 1.00
      temp =1.0d0
      ixx = 1
      do i=0,113
      whichatoms(i) =.false.
      end do
c  ========== init finished ================

      itot = 0

      do JATOM=1,natyp
         CALL ELEMENTDATABASE(Zatom(JATOM),elem_name,color,radious)
         I0 = Zatom(JATOM)
         do n=0,nr              ! change it to see more atoms
            PLOTXa = RR(1,N)+RBASIS(1,JATOM)
            PLOTYa = RR(2,N)+RBASIS(2,JATOM)
            PLOTZa = RR(3,N)+RBASIS(3,JATOM)
            if ((plotxa.lt.bbox(1)).and.(plotxa.gt.-0.01d0).and.
     &           (plotya.lt.bbox(2)).and.(plotya.gt.-0.01d0).and.
     &           (plotza.lt.bbox(3)).and.(plotza.gt.-0.01d0)) then
              ! name5 = elem_name+xx2
               name24= elem_name // '                     '            
c                                 123456789012345678901           
               x = RR(1,N)*alat*0.529+RBASIS(1,JATOM)*alat*0.529d0
               y = RR(2,N)*alat*0.529+RBASIS(2,JATOM)*alat*0.529d0
               z = RR(3,N)*alat*0.529+RBASIS(3,JATOM)*alat*0.529d0
c               write(6,800) x,y,z,rr(1,n),rr(2,n),rr(3,n),
c     &       rbasis(1,jatom),rbasis(2,jatom),rbasis(3,jatom)
c 800           format(10F8.3)
               itot = itot + 1
               if (itot.gt.nbig) stop 'increase NBIG, in pdbplot' 
               plotx(itot) = x
               ploty(itot) = y
               plotz(itot) = z
               name(itot) = xx1 // elem_name
                                ! write out in Angstrom
            end if
         end do
       end do ! jatom

       NPLOT = itot

       write(6,*) 'bbox',bbox(1),bbox(2),bbox(3)
       write(6,*) 'Includes ',nplot,' atoms'

       do ia=1,nplot
         x = plotx(ia)
         y = ploty(ia)  
         z = plotz(ia)
      !   name5 = name(ia) 
         write(27,5000) 'ATOM  ',ia,name(ia),xx1,'XXX','1',1,xx1,
     &       x,y,z,occup,temp,ixx,xx4,xx2,xx2
       end do ! ia plot atoms
c
c Now calculate conectivity take out vacancies
c
       natoms = 0
       do ia=1,nplot
!         write(6,*) ia,name(ia)
         if (name(ia).ne.' H1 ') THEN
!             write(6,*) 'chosen',ia,name(ia)
             natoms = natoms  + 1
             atomx(natoms) = plotx(ia)
             atomy(natoms) = ploty(ia)
             atomz(natoms) = plotz(ia)
             ide(natoms) = ia
         END IF
       end do 
       write(6,*) ' No of atoms ',natoms
       if (nplot.gt.255) then
       write(6,*) 'Reduce atoms to get connectivity'
       else
        
cc -------------------------------------------------
c  Now calculate connectivity and plot
cc--------------------------------------------------
          do iatom =1,natoms
             do ia=1,natoms
                rsort(ia) = SQRT((atomx(ia)-atomx(iatom))**2+
     &                           (atomy(ia)-atomy(iatom))**2+
     &                           (atomz(ia)-atomz(iatom))**2)
           !  write(6,*) 'atoms',rsort(ia)
             end do    
             
c     
             CALL DSORT(RSORT,ISORT,Natoms,POS)
!             write(6,*) ' Ordering ',pos
c     Now use correct order, first atom is always the zero!
             do IA =1,natoms
                IB = ISORT(IA)
                irevsort(ia) = ide(ib)
                rg(ia) = rsort(ib)
               ! write(6,*) 'lalala',ia,ib,rg(ia)
              END DO
              
c     
         test = rg(1)
         if (abs(test).gt.1.d-5) STOP ' Error connectivity'
c now find indeces to first neigbours only
         nnn = 0
         do ia=2,natoms-1
            test = rg(ia)
            test1 = rg(ia+1)
            if (abs(test-test1).gt.0.001d0) then
             nnn = ia-1
             EXIT
            END IF
         end do   
    !     write(6,*) 'atom ',iatom,ide(iatom),' has ',nnn,' neighbors'
c     
        if (nnn.gt.0.and.nnn.le.4) 
     &       write(27,4000) (irevsort(i),i=1,nnn+1)
        if (nnn.gt.4.and.nnn.le.8) then
            write(27,4000) (irevsort(i),i=1,5)
            write(27,4000) irevsort(1),(irevsort(i),i=6,nnn+1)
        end if
        if (nnn.gt.8.and.nnn.le.12) then
            write(27,4000) (irevsort(i),i=1,5)
            write(27,4000) irevsort(1),(irevsort(i),i=6,9)
            write(27,4000) irevsort(1),(irevsort(i),i=10,nnn+1)
        end if
        if (nnn.gt.12.and.nnn.le.16) then
            write(27,4000) (irevsort(i),i=1,5)
            write(27,4000) irevsort(1),(irevsort(i),i=6,9)
            write(27,4000) irevsort(1),(irevsort(i),i=10,13)
            write(27,4000) irevsort(1),(irevsort(i),i=14,nnn+1)
        end if
        if (nnn.gt.16.and.nnn.le.20) then
            write(27,4000) (irevsort(i),i=1,5)
            write(27,4000) irevsort(1),(irevsort(i),i=6,9)
            write(27,4000) irevsort(1),(irevsort(i),i=10,13)
            write(27,4000) irevsort(1),(irevsort(i),i=14,17)
            write(27,4000) irevsort(1),(irevsort(i),i=18,nnn+1)
        end if

        end do ! iatom

        end if ! nplot.gt.255 
c
c  Now write necesary colors
c
         do JATOM=1,natyp
         CALL ELEMENTDATABASE(Zatom(JATOM),elem_name,color,radious)
         I0 = Zatom(JATOM)
         
         name24= '#####' // elem_name // '                '            
c                 12345         678 9012345678901234   
         if (.not.whichatoms(i0)) then
         write(27,6000) name24,(color(i),i=1,3),radious
         whichatoms(i0) = .true.
         end if
        end do 
      write(27,*)'END'
      close(27)

c      write(27,5000) 'ATOM  ',ind,name5,xx1,xx3,1,xx1,x,y,z,occup,temp,
c     &                ixx,xx4,xx2,xx2)
c      write(27,6000)name24,(color(i),i=1,3),radious 
c             
 4000 format('CONECT',5I5)
 5000 FORMAT  
     & (A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3,2F6.2,1X,I3,2X,A4,2A2)
 6000 FORMAT('COLOR ',A24,3F9.3,F6.2,A10)
      END 





