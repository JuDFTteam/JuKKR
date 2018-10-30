      SUBROUTINE TETRAEDER(nkpoid,nkpoid_n,mesh_n,k,k_compl,n_max,
     +   lines,rows,n_trian,trian)

c     in this subroutine the volume of 1/48 of the first brillouinzone is
c     calculated. The triangles calculated in TRIANGLES are connected to
c     the origin to form tetraeder

c     the file triangles is constructed which allows to visualize the
c     triangles

      implicit none
      
      integer          ::  n,j,rows,lines,n_trian,
     +                     nkpoid,nkpoid_n,n_max
      double precision ::  a(3),b(3),c(3),pi,
     +                     volume,volume_one,
     +                     k(3,nkpoid+nkpoid_n),
     +                     k_compl(3,(nkpoid+nkpoid_n)*48)
      integer          ::  trian(n_max**2+2*lines*rows+rows*3,3),
     +                     mesh_n(rows,lines),nb,nges

c      write(6,*) "lines:", lines
c      write(6,*) "rows:", rows
c      write(6,*) "n_trian:", n_trian

      pi=3.141592653589793238

      open(unit=56,file="triangles", form="formatted")
      open(unit=55,file="triangles_for_tetraeder", form="formatted")
      volume=0.d0
      WRITE(6,*) "n_trian in tetraeder", n_trian
      do n=1,n_trian
c        write(6,"((A5),(4(I5,2x)))") "n:",n,
c     +                  (trian(n,j),j=1,3)
        do j=1,3
          a(j)=k(j,trian(n,1))
          b(j)=k(j,trian(n,2))
          c(j)=k(j,trian(n,3))
        end do
        write(55,"((I5),(3(e16.8,2X)))") n,
     +                  (a(j),j=1,3)
        write(55,"((I5),(3(e16.8,2X)))") n,
     +                  (b(j),j=1,3)
        write(55,"((I5),(3(e16.8,2X)))") n,
     +                  (c(j),j=1,3)
        write(55,"((I5),(3(e16.8,2X)))") n,
     +                  (a(j),j=1,3)
        write(55,"((I5),(3(e16.8,2X)))") 
        write(55,"((I5),(3(e16.8,2X)))") 
 
        write(56,"((I5),(3(e16.8,2X)))") n,
     +                  (a(j),j=1,3)
        write(56,"((I5),(3(e16.8,2X)))") n,
     +                  (b(j),j=1,3)
        write(56,"((I5),(3(e16.8,2X)))") n,
     +                  (c(j),j=1,3)
        write(56,"((I5),(3(e16.8,2X)))") n,
     +                  (a(j),j=1,3)
        write(56,"((I5),(3(e16.8,2X)))") 
        write(56,"((I5),(3(e16.8,2X)))") 

        CALL VOLUME_TETR(a,b,c,volume_one)
c        write(6,"((A15,X),(I5,X),(e16.8))") "n,volume_one:",n,volume_one
        volume=volume+volume_one
      end do
      write(6,*) "volume WO neck:",48*volume

      do n=1,rows-1
c        write(6,"((A5,X),(I5))") "row",nkpoid+mesh_n(n,lines)
        do j=1,3
          a(j)=k(j,nkpoid+mesh_n(n,lines))
          b(j)=k(j,nkpoid+mesh_n(n+1,lines))
          c(j)=0.5d0
        end do
        write(55,"((I5),(3(e16.8,2X)))") n+n_trian,
     +                  (a(j),j=1,3)
        write(55,"((I5),(3(e16.8,2X)))") n+n_trian,
     +                  (b(j),j=1,3)
        write(55,"((I5),(3(e16.8,2X)))") n+n_trian,
     +                  (c(j),j=1,3)
        write(55,"((I5),(3(e16.8,2X)))") n+n_trian,
     +                  (a(j),j=1,3)
        write(55,"((I5),(3(e16.8,2X)))") 
        write(55,"((I5),(3(e16.8,2X)))") 
 
        CALL VOLUME_TETR(a,b,c,volume_one)
c        write(6,"((A15,X),(I5,X),(e16.8))") "n,volume_one:",n,volume_one
        volume=volume+volume_one
      end do

      close(55)
      close(56)


      write(6,*) "volume of the first BZ:", 48.d0*volume

      open(unit=56,file="triangles_complete", form="formatted")
      WRITE(6,*) "n_trian in tetraeder", n_trian
      nges=0
      do nb=1,48
        do n=1,n_trian
        nges=nges+1
c        write(6,"((A5),(4(I5,2x)))") "n:",n,
c     +                  (trian(n,j),j=1,3)
        do j=1,3
          a(j)=k_compl(j,trian(n,1)+(nb-1)*(nkpoid+nkpoid_n))
          b(j)=k_compl(j,trian(n,2)+(nb-1)*(nkpoid+nkpoid_n))
          c(j)=k_compl(j,trian(n,3)+(nb-1)*(nkpoid+nkpoid_n))
        end do
 
        write(56,"((I8),(3(e16.8,2X)))") nges,
     +                  (a(j),j=1,3)
        write(56,"((I8),(3(e16.8,2X)))") nges,
     +                  (b(j),j=1,3)
        write(56,"((I8),(3(e16.8,2X)))") nges,
     +                  (c(j),j=1,3)
        write(56,"((I8),(3(e16.8,2X)))") nges,
     +                  (a(j),j=1,3)
        write(56,"((I8),(3(e16.8,2X)))") 
        write(56,"((I8),(3(e16.8,2X)))") 

        end do
      end do
      close(56)

      END SUBROUTINE TETRAEDER
