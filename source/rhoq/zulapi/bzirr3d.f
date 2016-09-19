!info 
!info   find irreducible BZ and create mesh in it.
!info           original version Arthur Ernst, Daresbury, 27/03/98
!                                fcc lattice bug corrected
!                                modified 26/5/99, 2/06/99
C Modified on 20.01.2000 To use the 2d inversion!
      subroutine bzirr3d(nkp,nkxyz,kpoibz,kp,recbv,bravais,
     &                   rfctor,wtkp,VLBZ,rsymat,nsymat,ISYMINDEX,irr)
!===========================================================================
!  Input:
!         nkxyz : original k-mesh net in the 3 directions of the reciprocal
!                 lattice vectors (not xyz directions).
!         recbv : reciprocal unit cell 
!                 (normalized: recbv(i,j)=G(i,j)*alat/(2*pi), i=x,y,z j=1,2,3)
!                 if G is the "normal" rec. lattice vector.
!      bravais  : direct unit cell in a.u.
!      rfctor   : normalization factor = alat/(2*pi)
!      rsymat   : symmetry rotation matrices of real lattice
!      nsymat   : number of rotation matrices of real lattice
!         irr   : if .true. then irreducible BZ else take all BZ
!  Output:
!         nkp   : number of points in irreducible BZ
!         kp    : k-point mesh
!         wtkp  : weights for k-points
!        VOLBZ  : volume of the BZ
!  Inside:
!         ibk   : Flag showing if the mesh point jx,jy,jz has already been 
!                 taken. Could also be a logical variable.
!==========================================================================
      implicit none
      integer MAXK1,MAXK2,MAXK3,NSYMAXD
      PARAMETER (MAXK1=500,MAXK2=500,MAXK3=500,NSYMAXD=48)
!i/o
      integer nkp,nkxyz(3),ibk(0:MAXK1,0:MAXK2,0:MAXK3),kpoibz,nsymat
      INTEGER ISYMINDEX(*),nkxyz1(3)
      double precision kp(3,*),wtkp(kpoibz),rfctor,VLBZ
      double precision recbv(3,3),bravais(3,3),RSYMAT(64,3,3)
      logical irr
!local
      integer i,j,jx,jy,jz,nsym,iws,iwt,is,ja,jb,jc,ix,iy,iz,k,nk,n,isym
      integer i0,n0,INVERS2D(64)
      double precision u(3,3,48),rq(3,3),gq(3,3),x,y,z,a,b,c,ox,oy,oz, 
     &     w0,pi,alat,v1
      logical linver,lsurf,test
      parameter(pi=3.1415926535897932384626d0)
      double precision ddet33
      external ddet33
!==========================================================================
C
C --->  ALAT=2*PI*RFCTOR
C
c Check if we are in surface mode
      lsurf = .false.
      if (bravais(1,3).eq.0.d0.and.bravais(2,3).eq.0.d0.and.
     &     bravais(3,3).eq.0.d0) lsurf = .true. 

      if (lsurf) nkxyz(3) = 1
        ALAT = RFCTOR*2.d0*pi
      write(6,*) 'bzirr3d:',nkp,kpoibz,rfctor,'lsurf:',lsurf
c     write(6,*) (nkxyz(i),i=1,3)
      if (nkxyz(1).GT.MAXK1.or.nkxyz(2).GT.MAXK2.or.nkxyz(3).GT.MAXK3)
     &    then
          write(6,*) 'Sub BZIRR3d: Increase MAXK ',(nkxyz(i),i=1,2,3),
     &                                              MAXK1,MAXK2,MAXK3
          STOP
      end if
      do i=1,3
      nkxyz1(i) = nkxyz(i)
      end do 
c-------------------
c     if (mod(nkxyz(1),2).ne.0) nkxyz1(1) = nkxyz(1)+1
c     if (mod(nkxyz(2),2).ne.0) nkxyz1(2) = nkxyz(2)+1
c     if (.not.lsurf) then 
c     if (mod(nkxyz(3),2).ne.0) nkxyz1(3) = nkxyz(3)+1    
c     end if
c-------------------
c
c Create small unit cell for the integration in the reciprocal space (gq),
c and its reciprocal cell (big!) in the real space normalized by 2*pi (rq).
c
      do i = 1, 3
         do j = 1, 3
!             ! test
!             if(i==2) recbv(2,j) = (2.0d0*recbv(2,j) + recbv(1,j))/2.0d0
!             if(i==2) bravais(2,j) = (2.0d0*bravais(2,j) + 
!      &               bravais(1,j))/2.0d0
!             ! test end
            gq(i,j)=recbv(i,j)/nkxyz1(j)
            rq(i,j)=bravais(i,j)*nkxyz1(j) !/alat
         enddo 
      enddo
!===========================================================================
      if(irr) then
         if (.not.lsurf) THEN
            write(6,*) "RECBV",(RECBV(1,j),j=1,3)
            write(6,*) "RECBV",(RECBV(2,j),j=1,3)
            write(6,*) "RECBV",(RECBV(3,j),j=1,3)
c            write(6,*) "nsym",nsym
            call symlat(nsym,recbv,u)
            write(6,191) nsym
 191        FORMAT('Reciprocal lattice has ',I4,' symmetries')
         END IF                 ! (.not.lsurf)

         nsym = nsymat
         write(6,192) nsym
 192     FORMAT('The real lattice ',I4,'  symmetries will be used')
         LINVER=.true.
         if (TEST('noinvers').OR.TEST('fullBZ  ')) LINVER=.false.
c     
c     In case of surface the Inversion  matrix cannot be in the symmetry  
c     operations (sub findgroup)
c     
         IF (LSURF) THEN
c     In case of surface look if the in plane inversion is included            
c     if not add it explicitly!
            do i=1,nsym             
               IF (ISYMINDEX(i).eq.12) LINVER=.false.
            end do   
            
            IF (LINVER) THEN
c     
c     Since I am in here the 2d-inversion is not a symmetry opperation
c     so I add it in any case
c     
               write(6,*) 
     &      '** In Plane Inversion is not a real space symmetry   **'
             write(6,*)
     &      '** it will be used in the 2D reciprocal space        **'
             write(6,*) 
     &      '##   This will double the symmetry operations   ##     '
c
c Now include all the symmetries by hand!
c            
             DO I=1,64
                INVERS2d(i) = 100
             end do
             INVERS2d(1)  = 12  !E 
             INVERS2d(12) = 1   !C2z
             INVERS2d(15) = 18  !C4z
             INVERS2d(18) = 15  !C4z-1
             INVERS2d(34) = 35  !IC2x
             INVERS2d(35) = 34  !IC2y
             INVERS2d(43) = 44  !IC2a
             INVERS2d(44) = 43  !IC2b
             
             INVERS2d(49) = 52  ! C3z 
             INVERS2d(52) = 49  ! C6z-1
             INVERS2d(50) = 51  ! C3z-1
             INVERS2d(51) = 50  ! C6z
             INVERS2d(61) = 63  ! IC3z
             INVERS2d(63) = 61  ! IC6z
             INVERS2d(62) = 64  ! IC3z-1
             INVERS2d(64) = 62  ! IC6z-1
c     Test --------------------------
             do I=1,nsym
                i0 = ISYMINDEX(i)
                if (INVERS2d(i0).eq.100) then
                   write(6,*) 'ERROR in 2D BZIRR3d',i,i0
c     If this occurs the symmetry op. does not leave the z-axis unchanged!
c     look at sub pointgrp.f
                   STOP
                end if 
             end do
c--------------------------------
             N0 = 0
             do i=1,nsym
                i0 = ISYMINDEX(i) 
                N0 = N0 + 1 
                ISYMINDEX(NSYM+N0) = INVERS2d(I0)
             end do 
             NSYM = NSYM + N0
             write(6,*) 
     &            '     Number of symmetries used in the 2D  BZ :',NSYM
c
          END IF                ! LINVER
c     
c     
c     
c     
       ELSE                     ! (LSURF)
c
c
          do i=1,nsym
             if (ISYMINDEX(i).eq.25) LINVER=.false. 
c     matrix no.25 is the inversion 
          end do            
          
          
          IF (LINVER) THEN
             write(6,*) '** Inversion is not a real space symmetry   **'
             write(6,*) '** it will be used in reciprocal space      **'
             write(6,*) '## This will double the symmetry operations ##'
             write(6,*) 
     &                  '   Number of symmetries used for BZ :',2*NSYM
             do i=1,NSYM
                i0 = ISYMINDEX(i)
                if (i0.le.24)              ISYMINDEX(NSYM+i) = I0 + 24
                IF (I0.gt.24.and.I0.le.48) ISYMINDEX(NSYM+i) = I0 - 24
                IF (I0.gt.48.and.I0.le.56) ISYMINDEX(NSYM+i) = I0 + 8
                IF (I0.gt.56.and.I0.le.64) ISYMINDEX(NSYM+i) = I0 - 8
             end do
             NSYM = 2*NSYM
          END IF
          
       END IF                   ! (LSURF)   
       
         
         
         if (nsym.gt.48) then
            write(6,*) 'Dim problem in bzirr3d.'
            STOP
         end if
         if (lsurf) then
            if (nsym.gt.12) then
               write(6,*) 'Dim problem in bzirr3d, surf mode.'
               STOP
            end if
         endif
         do n=1,nsym
             isym = isymindex(n)
             do i=1,3
                do j=1,3
                u(i,j,n) = rsymat(isym,i,j)
                end do
             end do
         end do
         
      else
         nsym=1
         u(1,1,1)=1.d0
         u(1,2,1)=0.d0
         u(1,3,1)=0.d0
         u(2,1,1)=0.d0
         u(2,2,1)=1.d0
         u(2,3,1)=0.d0
         u(3,1,1)=0.d0
         u(3,2,1)=0.d0
         u(3,3,1)=1.d0
      endif
!==========================================================================
      do i=0,nkxyz1(1)
         do j=0,nkxyz1(2)
            do k=0,nkxyz1(3)
               ibk(i,j,k)=0
            enddo
         enddo
      enddo
      nk=nkxyz1(1)*nkxyz1(2)*nkxyz1(3)
      ix=nkxyz1(1)
      iy=nkxyz1(2)
      iz=nkxyz1(3)
      nkp=0
      iws=0
      do jx = 0, ix-1
         do jy = 0, iy-1
            do jz = 0, iz-1

               if(ibk(jx,jy,jz).eq.0) then
                  nkp=nkp+1
                  iwt=0
                  x=jx*gq(1,1)+jy*gq(1,2)+jz*gq(1,3)
                  y=jx*gq(2,1)+jy*gq(2,2)+jz*gq(2,3)
                  z=jx*gq(3,1)+jy*gq(3,2)+jz*gq(3,3)
                  do is=1,nsym
                     ox=u(1,1,is)*x+u(1,2,is)*y+u(1,3,is)*z
                     oy=u(2,1,is)*x+u(2,2,is)*y+u(2,3,is)*z
                     oz=u(3,1,is)*x+u(3,2,is)*y+u(3,3,is)*z

c Check if the rotated k-point belongs to the mesh.
                     a=ox*rq(1,1)+oy*rq(2,1)+oz*rq(3,1)
                     b=ox*rq(1,2)+oy*rq(2,2)+oz*rq(3,2)
                     c=ox*rq(1,3)+oy*rq(2,3)+oz*rq(3,3)

                     ja=dnint(a)
                     jb=dnint(b)
                     jc=dnint(c)
                     if(abs(ja-a)+abs(jb-b)+abs(jc-c).gt.1.e-3) then
                        write(*,*)"ERROR in bzirr3d!"
                     endif
c Check passed (or not if "ERROR in bzirr3d!" occured).
c
c Now find the mesh point in the unit cell which is equivalent by translation
c to the rotated one.
                     ja=mod(ja,ix)
                     if(ja.lt.0) ja=ja+ix
                     jb=mod(jb,iy)
                     if(jb.lt.0) jb=jb+iy
                     jc=mod(jc,iz)
                     if(jc.lt.0) jc=jc+iz
                     if(ibk(ja,jb,jc).eq.0) then
                        ibk(ja,jb,jc)=nkp
                        iwt=iwt+1
                     endif
c Mesh point in the unit cell found.
                  enddo
                  
                  kp(1,nkp)=x
                  kp(2,nkp)=y
                  kp(3,nkp)=z
                  wtkp(nkp)=dble(iwt)/nk
                  if (nkp.gt.kpoibz) then
                     write(6,*) 'Bzirr3d increase kpoibz '
                     STOP
                  end if
                  iws=iws+iwt
               endif
            enddo
         enddo
      enddo
      if (lsurf) then
         v1 = dabs(recbv(1,1)*recbv(2,2)-recbv(1,2)*recbv(2,1))
      else
         v1 = ddet33(recbv) 
      endif
      VLBZ = 0d0
      do i=1,nkp
         wtkp(i) = wtkp(i)*v1/dfloat(nsym)
         VLBZ = VLBZ + wtkp(i)*dfloat(nsym)
      end do   
            write(6,*) 'bzirr3d:',nkp,kpoibz,rfctor,'lsurf:',lsurf
            write(6,*) (nkxyz1(i),i=1,3)    
      return
      end






