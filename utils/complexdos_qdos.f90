program complexdos
   implicit none
   ! principle of dos here: two-point contour integration
   ! for dos in the middle of the two points. the input dos
   ! and energy must be complex. parameter deltae should be
   ! of the order of magnitude of eim. 
   !
   !    
   !      <-2*deltae->   _
   !           /\        |     dos=(n(1)+n(2))/2 + (n(1)-n(2))*eim/deltae
   !          /  \       |
   !        (1)  (2)   2*i*eim=2*i*pi*kb*tk
   !        /      \     |
   !       /        \    |
   !------------------------ (real e axis)
   integer*4 iemaxd,lmaxd
   parameter(iemaxd=1000,lmaxd=10)
   integer*4 npot,iemax,lmax
   real*8 eim,deltae,tk,kb,pi
   ! total dos stored in dos(lmax+1,ie)
   complex*16 dos(0:lmaxd+1,iemaxd),ez(iemaxd)
   complex*16 dosnew(0:lmaxd+1), dentot
   real*8 temp(2*lmaxd+6)
   real*8 enew,ef,ev
   integer*4 ie,ii,ll,iheader, ispin, iq, ierr, nqdos, nspin, natyp, dumint
   double precision, allocatable :: qvec(:,:)
   double complex, allocatable :: den(:,:,:)
   character(len=80) :: text
   character(len=45) :: text1
   character(len=20) :: text2
   character(len=80) :: tmpname

   ev = 13.6058

   ! if only tk is known, use bolzmann constant, pi to find eim:
   ! kb=0.6333659d-5
   ! pi=3.14159265358979312d0
   ! eim=pi*kb*tk
   
   ! find parameters npot, nspin, and nqdos
   tmpname = "cqdos.01.1.dat"
   open (49,file=tmpname, form='formatted')
   ! read header
   read(49,9002) text
   if(text(1:9)/='# serial:') rewind(49)
   read(49,9002) text
   read(49,*) lmax, npot, nspin, nqdos, iemax
   read(49,9002) text
   allocate(qvec(3,nqdos), den(0:lmaxd+1,iemaxd, nqdos), stat=ierr)
   if(ierr/=0) stop 'error qvec allocation'
   

   do ii=1,npot
      do ispin=1,nspin
      
         ! open cqdos file
         if(npot>=100) then
            tmpname = "cqdos."//char(48+ii/100)//char(48+mod(ii/10,10))//char(48+mod(ii,10))//"."//char(48+ispin)//".dat"
         else
            tmpname = "cqdos."//char(48+mod(ii/10,10))//char(48+mod(ii,10))//"."//char(48+ispin)//".dat"
         end if
         open (49,file=tmpname, form='formatted')
         rewind(49)
         
         
         ! open output file
         if(npot>=100) then
            tmpname = "qdos.interpol."//char(48+ii/100)//char(48+mod(ii/10,10))//char(48+mod(ii,10))//"."//char(48+ispin)//".dat"
         else
            tmpname = "qdos.interpol."//char(48+mod(ii/10,10))//char(48+mod(ii,10))//"."//char(48+ispin)//".dat"
         end if
         open (50,file=tmpname, form='formatted')
         rewind(50)
         
         
         write(*,*) 'reading potential',ii, ispin
         ! read header:
         read(49,9002) text
         if(text(1:9)/='# serial:') then
            write(*,*) ' serial number not found in file: rewind'
            rewind(49)
         else
            write(50,9002) text
         end if
         read(49,9002) text
         read(49,*) lmax, natyp, dumint, nqdos, iemax
         read(49,9002) text
         write(50,9002) text
         
         ! read dos: (total dos stored at dos(lmax+1,ie))
         do ie=1,iemax
           do iq=1,nqdos
            read(49,9001) ez(ie),qvec(1,iq),qvec(2,iq),qvec(3,iq), dentot, (den(ll,ie,iq),ll=0,lmax+1)
           end do ! iq=1,nqdos
         end do ! ie=1,ielast
         
         ! compute and write out corrected dos at new (middle) energy points:
         do iq=1,nqdos
           do ie=2,iemax-1
            deltae = dreal(ez(ie+1) - ez(ie))
            eim = dimag(ez(ie))
            enew = dreal(ez(ie)) ! real quantity
            if(iq==1) WRITE(*,'(i3,A,i3,A,3es14.7)') ie-1,' out of ',iemax-2, '; ', deltae, eim, enew

            do ll=0,lmax+1
               dosnew(ll) = den(ll,ie,iq) + 0.5d0*(den(ll,ie-1,iq)-den(ll,ie+1,iq))*dcmplx(0.d0,eim)/deltae
            enddo
            
            
            dentot = dcmplx(0.d0,0.d0)
            do ll = 0,lmax+1
              dentot = dentot + dosnew(ll)
            enddo
            !               Re(e) Im(e)   qx       qy         qz            tot            s,p,d,...,ns    
            write(50,9001) enew,0.0d0,qvec(1,iq),qvec(2,iq),qvec(3,iq),dimag(dentot),(dimag(dosnew(ll)),ll=0,lmax+1)

           end do ! ie=2,iemax-1
         end do ! iq=1,nqdos
         
         close(49)
         close(50)
         
      end do ! ispin
   
   enddo ! ii=1,npot
   
   deallocate(qvec)
   
   stop 'telos'
   9001 format(5f10.6,80e16.8)
   9002 format(a80)
   
end
   
                                                                           