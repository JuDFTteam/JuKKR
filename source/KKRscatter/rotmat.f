      subroutine rotmat(iopt,li,nrot,symopm,vecg) 
!- Converts rotation/rotoinversion matrix <-> (nrot,vecg,li)            
! ----------------------------------------------------------------------
!i Inputs:                                                              
!i   iopt  := -1 to convert (nrot,vecg,li) in symopm                    
!i          =  1 to convert symopm in (nrot,vecg,li)                    
!i Inputs/Outputs:                                                      
!io  li    :if T: inversion or rotoinversion                            
!io  nrot  :rotation angle = 2*pi/nrot                                  
!io  symopm:symmetry operation matrix                                   
!io  vecg  :rotation axis                                               
! ----------------------------------------------------------------------
      implicit none 
! Passed parameters:                                                    
      integer iopt,nrot 
      double precision vecg(3),symopm(3,3) 
      logical li 
! Local parameters:                                                     
      integer i,idamax,in,j,k 
      double precision costbn,detop,ddet33,dnrm2,omcos,                 
     &                 sintbn,sinpb3,tiny,twopi,vfac                    
      character*144 messg 
      parameter(twopi=6.28318530717958648d0) 
      parameter(tiny=1.0d-3) 
! External calls:                                                       
      external daxpy,dcopy,ddet33,dinit,dnrm2,dscal,errmsg,idamax,nrmliz 
! Intrinsic functions:                                                  
      intrinsic  dabs,dacos,dcos,dmax1,dsign,dsin,dsqrt,iabs,idnint 
                                                                        
      if (iopt.eq.-1) then 
        call dinit(symopm,9) 
        in = iabs(nrot) 
        if (in.eq.1) then 
          call dcopy(3,1.d0,0,symopm,4) 
        elseif (in.eq.2.or.in.eq.3.or.in.eq.4.or.in.eq.6) then 
          sintbn = dsin(twopi/nrot) 
          costbn = dcos(twopi/nrot) 
          omcos  = 1.d0-costbn 
c          if (dnrm2(vecg).lt.tiny)                                      
c     &      call errmsg(' ROTMAT: zero rotation vector.$',4)            
          call nrmliz(1,vecg,vecg) 
          symopm(1,1)=omcos*vecg(1)*vecg(1)+costbn 
          symopm(1,2)=omcos*vecg(1)*vecg(2)-sintbn*vecg(3) 
          symopm(1,3)=omcos*vecg(1)*vecg(3)+sintbn*vecg(2) 
          symopm(2,1)=omcos*vecg(2)*vecg(1)+sintbn*vecg(3) 
          symopm(2,2)=omcos*vecg(2)*vecg(2)+costbn 
          symopm(2,3)=omcos*vecg(2)*vecg(3)-sintbn*vecg(1) 
          symopm(3,1)=omcos*vecg(3)*vecg(1)-sintbn*vecg(2) 
          symopm(3,2)=omcos*vecg(3)*vecg(2)+sintbn*vecg(1) 
          symopm(3,3)=omcos*vecg(3)*vecg(3)+costbn 
        else 
          call errmsg(' ROTMAT: bad nrot.$',3) 
        endif 
        if (li) call dscal(9,-1.d0,symopm(1,1),1) 
                                                                        
      elseif (iopt.eq.1) then 
! ----- First calculate determinant.                                    
        detop=ddet33(symopm) 
        if (dabs(dabs(detop)-1.0d0).gt.tiny)                            
     &    call errmsg(' ROTMAT: determinant is not +/- 1$',4)           
        detop=dsign(1.d0,detop) 
        li=detop.lt.0.d0 
! ----- multiply operation symopm with detop                            
        call dscal(9,detop,symopm(1,1),1) 
! ----- For the rotation angle we have due to the normalization of v:   
! ----- sum_i symopm(i,i) = sum_i (1-cos) v_i*v_i+3*cos = 1 + 2 * cos,  
        costbn=-0.5d0 
        call daxpy(3,0.5d0,symopm(1,1),4,costbn,0) 
      write(*,*) 'rotmat'                   
        if (dabs(costbn-1.d0).lt.tiny) then 
          nrot=1 
          call dinit(vecg,3) 
        else 
          nrot=idnint(twopi/dacos(dmax1(-1.d0,costbn))) 
! ------- for nrot > 2 the matrix is non-symmetric and the rotation     
! ------- axis can be calculated from the antisymmetric part.           
! ------- for nrot = 2 this not possible. However, the squared vector   
! ------- components are given by:  mat(i,i) = 2 v_i * v_i - 1.         
! ------- This is used for the largest component. The others are taken  
! ------- from: mat(i,j) = 2 v_i * v_j for i ne j. This way we also     
! ------- get the right phases between the components.                  
          if (nrot.eq.2) then 
            do i=1,3 
              vecg(i)=0.5d0*(symopm(i,i)+1.0d0) 
            enddo 
            j=idamax(3,vecg,1) 
            if (vecg(j).lt.0.0d0) then 
              write(messg,402)j,symopm(j,j) 
              call errmsg(messg,4) 
            endif 
            vecg(j)=dsqrt(vecg(j)) 
            vfac=0.5d0/vecg(j) 
            do i=1,3 
              if (i.ne.j) vecg(i)=vfac*symopm(i,j) 
            enddo 
          else 
            vecg(1)=symopm(3,2)-symopm(2,3) 
            vecg(2)=symopm(1,3)-symopm(3,1) 
            vecg(3)=symopm(2,1)-symopm(1,2) 
          endif 
! ------- next renormalize at least one component to 1 in order to      
! ------- allow for abbreviations as 'D', 'X', 'Y' or 'Z'               
          sinpb3=dsqrt(.75d0) 
          if (dabs((sinpb3-dabs(vecg(1)))*(sinpb3-dabs(vecg(2)))*       
     &            (sinpb3-dabs(vecg(3)))).gt.tiny) then                 
           do j=3,1,-1 
             vfac=dabs(vecg(j)) 
             if(vfac.gt.tiny) call dscal(3,1.d0/vfac,vecg,1) 
           enddo 
          endif 
        endif 
        call dscal(9,detop,symopm(1,1),1) 
      endif 
                                                                        
  402 format(' ROTMAT: Bad component ',i1,' of operation ',             
     &       '. Diagonal element =',f9.5,'$')                           
      END                                           



