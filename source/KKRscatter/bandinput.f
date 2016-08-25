       SUBROUTINE BANDINPUT(IE,ksum,BZKP,NKPOID)
       implicit none
c *****************************************************
c This subroutine reads info from inputcard and sets up
c the k-points along lines for band-structure plots
c Complex band structure is also possible, extra 
c OPT('COMPLEX') is neaded for this
c                                           6.4.2001
c *****************************************************
       INTEGER IE,DIRECNO,ier,ksum,ins,NKPOID
       DOUBLE PRECISION KPOINTA(6),KPOINTB(6),BZKP(6,*)
       double precision test,dxyz(6),dist
       character*80 UIO
       integer ixyz,i,iline,j,i1,NKMESH
       LOGICAL OPT
c
       
      
      CALL IoInput('DIRECNO   ',UIO,1,7,IER)
      READ (UNIT=UIO,FMT=*)  direcno
      
      IF (IE.eq.1) write (6,*) 'direcno  =',direcno
      
      KSUM = 0
      DO I=1,direcno
         ILINE = (I-1)*3 + 1
         CALL IoInput('DIRECDEF  ',UIO,ILINE,7,IER)
         READ (UNIT=UIO,FMT=*) NKMESH
         IF (IE.eq.1) write(6,*) NKMESH
         CALL IoInput('DIRECDEF  ',UIO,ILINE+1,7,IER)
         READ (UNIT=UIO,FMT=*)  (KPOINTA(J),J=1,6)
         IF (IE.eq.1) write(6,*)  (KPOINTA(J),J=1,6)
         CALL IoInput('DIRECDEF  ',UIO,ILINE+2,7,IER)
         READ (UNIT=UIO,FMT=*)  (KPOINTB(J),J=1,6)
         IF (IE.eq.1) write(6,*)  (KPOINTb(J),J=1,6)
c     ---
         test = abs(KPOINTA(4))+abs(KPOINTA(5))+abs(KPOINTA(6))+
     &        abs(KPOINTb(4))+abs(KPOINTb(5))+abs(KPOINTb(6))
         IF (test.gt.1.d-6) THEN            
            WRITE(6,*) '** COMPLEX   BAND   CALCULATION **'
            
            IF (.not.opt('COMPLEX   ')) THEN
               WRITE(6,*) 'USE OPTION : COMPLEX '
               WRITE(6,*) ' Program is STOPING!!'
               STOP 'BANDINPUT'
            END IF
            
            test = abs(KPOINTA(1)-KPOINTb(1))+
     &           abs(KPOINTA(2)-KPOINTB(2))+
     &           abs(KPOINTA(3)-KPOINTB(3))
            IF (TEST.GT.1.d-6) THEN
               WRITE(6,*) ' REAL PART IS CHANGING '
               STOP 'BANDINPUT'
            END IF  
            DIST = SQRT((kpointa(4)-kpointb(4))**2 +
     &           (kpointa(5)-kpointb(5))**2 +
     &           (kpointa(6)-kpointb(6))**2)
            
            IF (DIST.GT.0.D0.AND.NKMESH.GT.0) THEN
               
               do ixyz=4,6
                  dxyz(ixyz) = (kpointb(ixyz)-kpointa(ixyz))/
     &                 (nkmesh-1)
               end do
               
               do i1 = 1,nkmesh
                  ksum = ksum + 1
                  do ixyz=1,3
                     bzkp(ixyz,ksum) = kpointa(ixyz)
                  end do
                  
                  do ixyz=4,6
                     bzkp(ixyz,ksum) = kpointa(ixyz) + dxyz(ixyz)*(i1-1)
                  end do        ! ixyz
               end do           ! i1
               
            END IF

         ELSE                   ! real band structure
c     -------------------
            
            DIST = SQRT((kpointa(1)-kpointb(1))**2 +
     &           (kpointa(2)-kpointb(2))**2 +
     &           (kpointa(3)-kpointb(3))**2)
            
            IF (DIST.GT.0.D0.AND.NKMESH.GT.0) THEN
         
               do ixyz=1,3
                  dxyz(ixyz) = (kpointb(ixyz)-kpointa(ixyz))/
     &                 (nkmesh-1)
               end do
               
               do i1 = 1,nkmesh
                  ksum = ksum + 1
                  do ixyz=1,3
                     bzkp(ixyz,ksum) = kpointa(ixyz) + 
     &                    dxyz(ixyz)*(i1-1)    
                  end do        ! ixyz
                  do ixyz=4,6
                     bzkp(ixyz,ksum) = kpointa(ixyz)
                  end do
               end do           ! i1
               
            ELSE
               WRITE(6,*) ' POINTS ARE THE SAME OR NO POINTS !!!'
               STOP 'BANDINPUT'      
            END IF 
         END IF
      ENDDO
        
      IF (IE.EQ.1) THEN
      WRITE(6,*) "BANDINPUT",KSUM
      DO i1=1,KSUM
        WRITE(6,"(3e17.9)") (BZKP(J,I1),J=1,3)
      END DO
      IF (KSUM.GT.NKPOID) STOP ' INCREASE NKPOID in banstr_SO.f'    
        WRITE(6,*) 'K POINTS FOR BAND CALC. ',KSUM
      END IF
      
      RETURN 
      END


