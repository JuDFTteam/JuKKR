      SUBROUTINE GFMASK(ICHECK,ICC,INVMOD,NSH1,NSH2,NLAYER,naez,IE,
     &                     NSHELL,IATCONDL,IATCONDR,NCONDPAIR,IEGFOUT)
c     *************************************************************
c     * This subroutine prepares the ICHECK matrix that is used
c     * for calculating the proper GF matrix elements for impurity
c     * or transport calculations
c     *                                            29.02.2000
c     ************************************************************ 
      implicit none
      include 'inc.p'
      integer icc,invmod,nlayer,naez,ie,nshell,IEGFOUT
      INTEGER IATCONDL(*),IATCONDR(*),NCONDPAIR
      integer ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)    
      INTEGER NSH1(*),NSH2(*),ICOUPLE(NAEZD,NAEZD)
      integer i,j,ii,istep1,ilt1,istep2,ilt2,il2,il1
      integer k,i1,i2,j1,j2
      LOGICAL OPT

        IF (INVMOD.EQ.1) THEN
           
           DO I=1,NLAYER
              DO J=1,NLAYER
                 if (i.eq.j) then
                    ICHECK(I,J)=1
                 else
                    ICHECK(I,J)=0
                 endif
              ENDDO
           ENDDO
           
        ENDIF
        
        IF (INVMOD.EQ.2) THEN
           
           DO I=1,NLAYER
              DO J=1,NLAYER
                 if ((i.eq.j).OR.(J.EQ.NLAYER).OR.((I.EQ.NLAYER))) then
                    ICHECK(I,J)=1
                 else
                    ICHECK(I,J)=0
                 endif
                 
              ENDDO
           ENDDO
           
        ENDIF          ! end of changes look below also


      IF ((ICC.NE.0).AND.(INVMOD.EQ.1)) THEN
c
c Prepare the matrix icouple which has 1 in all nn' (atom) blocks
c that we nead .
c
         DO I=1,NAEZ
            DO J=1,NAEZ
               ICOUPLE(I,J) = 0
               DO II=1,NSHELL
                  IF (((NSH1(II).EQ.I).AND.(NSH2(II).EQ.J)).OR.
     +                 ((NSH1(II).EQ.J).AND.(NSH2(II).EQ.I))) THEN 
                     ICOUPLE(I,J) = 1 
                  ENDIF                       
               ENDDO
            ENDDO
         ENDDO

         IF (OPT('CONDUCT ')) THEN ! Conductivity calculation
            write(6,*) '::::::::::::::::::::::::::::::::::::::::::::::'
            WRITE(6,*) 'ICHECK MATRIX SET FOR CONDUCTIVITY CALCULATION'
             DO I=1,NLAYER
              DO J=1,NLAYER
                    ICHECK(I,J)=0
              ENDDO
           ENDDO

            do I=1,NAEZ
               do J=1,NAEZ
                  ICOUPLE(I,J) = 0
               end do
            end do
                  DO II=1,NCONDPAIR
                        I = IATCONDL(II)
                        J = IATCONDR(II)
                        WRITE(6,*) 'THE ',I,J,'SET TO ONE'

                           ICOUPLE(I,J) = 1
                        if (i.gt.j) then 
                        write(6,*) ' SUB GFMASK '
                        write(6,*) ' Put layers in different order '
                        STOP 'I am stopping'    
                        END IF                          
                          ! ICOUPLE(J,I) = 1
                        
                  ENDDO
c
c Now given the matrix icouple prepare the matric icheck
c which has 1 in all principal-layer blocks that we nead
c this will be used in the matrix inversion
c
         do i=1,naez
            i1 = mod(i,NPRINCD) + 1
            do j=i,naez ! upper triagle 
               j1 = mod(j,NPRINCD) + 1
               
c The element i,j -> i1, j1
c 
               if (ICOUPLE(I,J).eq.1) THEN
                 ICHECK(I1,J1) = 1                
                 
c      
c Now go down in the colunm until i1=j1 
c             
                 do ii=i1,j1
                     ICHECK(ii,j1) = 1  ! column i1
                 end do
                 do ii=j1,nlayer
                     ICHECK(ii,ii) = 1  ! diagonal 
                 end do                  
               end if 

            end do
         end do 
        
         ELSE! Conductivity calculation
c
c Now given the matrix icouple prepare the matric icheck
c which has 1 in all principal-layer blocks that we nead
c this will be used in the matrix inversion
c

         istep1=0
         ilt1=1
         do 3 il1=1,naez
            istep1=istep1+1
            if (istep1.gt.nprincd) then
               ilt1=ilt1+1
               istep1=1
            endif

            ilt2=1
            istep2=0
            do 4 il2=1,naez 
               istep2=istep2+1
               if (istep2.gt.nprincd) then
                  ilt2=ilt2+1
                  istep2=1
               endif
               if (icouple(il1,il2).eq.1) then
                  icheck(ilt1,ilt2)=1 
               endif
 5             continue
 4          continue
 3       continue      
         
c     added 24.2.2000
c     in the case of calculation of single blocks 
c     it has to put the correct value to icheck in order
c     calculate all the elements also necessary to calculate
c     that single block.

         do j=1,nlayer
c     loop over the element icheck(i,j) with fixed j and i<j
            if (j.ne.1) then
               do i=1,j-1
                  if (icheck(i,j).eq.1) then
                     do k=i+1,j
                        icheck(k,j)=1
                     enddo
                  endif
               enddo
            endif
c     loop over the element icheck(i,j) with fixed j and i>j
         
            if (j.ne.nlayer) then
               do i=nlayer,j+1,-1
                  if (icheck(i,j).eq.1) then
                     do k=i-1,j,-1
                        icheck(k,j)=1
                     enddo
                  endif
               enddo
            endif
          
         enddo

         END IF ! finished OPT('CONDUCT' ) LOOP  

         if ((ie.eq.1).OR.(ie.eq.iegfout)) then
c           write (6,*) 'this is the icheck matrix'
            do i=1,nlayer
c              write (6,FMT='(20i3)') (icheck(i,j),j=1,nlayer)
            enddo
         endif

      ENDIF      
c
      END



