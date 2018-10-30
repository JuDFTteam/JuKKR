      SUBROUTINE CREATE_GRID_FS_3D(KF_IRR,NKF,KF_COMPL,DIREC_VEL_FS)

      IMPLICIT NONE 
c *****************************************************
c This subroutine integrates the function f over the Fermisufarce
c needs the file kpoints_FS (contains the kpoints on the Fermisurface)
c  important variables:
c      
c      nkpoid_total            total number of kpoints
c      nkpoid_n                number of kpoints on the neck (for copper)
c      nkpoid_b                number of kpoints in the main part
c      k(nkpoid_total,3)       kpoints  
c      mesh(lines,rows)        numbers consecutively the kpoints     
c      area                    contains the area of the triangles 
c      border(nkpoid_total)    gives the number of triangles adjoining the
c                              kpoint 
c      border(nkpoid_total,20) contains the numbers of the triangles
c                              adjoining one kpoint 
c      trian(num of triangles) stores the kpoints adjoining one triangle                  

c     input variables

      INTEGER,INTENT(IN)            ::   NKF
      DOUBLE PRECISION,INTENT(IN)   ::   KF_IRR(3,*),KF_COMPL(3,*)
      DOUBLE PRECISION,INTENT(OUT)  ::   DIREC_VEL_FS(3,NKF)
c      DOUBLE PRECISION,INTENT(OUT)  ::   DIREC_VEL_FS(3,NKF,48)

c     local variables

      integer                       ::   NK,NB,I,NKPOID_B,NKPOID_N,
     +                                   K_MAX,N_TRIAN,
     +                                   n_max,
     +                                   num_rows,num_lines(800),N,J

      DOUBLE PRECISION              ::   K_STORE 

      DOUBLE PRECISION,ALLOCATABLE  ::   AREA(:),vel_FS(:)

      INTEGER,ALLOCATABLE           ::   border(:),mesh(:,:),i_max(:),
     +                                   border_tri(:,:),trian(:,:),
     +                                   mesh_n(:,:)
      CHARACTER*80 UIO
      LOGICAL                       ::   OPT

c evaluate the kpoints on the fermi surface, check first how many points are in the neck 

c      nkpoid_total=NKF
      nkpoid_n=0
      nkpoid_b=0
        
      WRITE(6,*) "NK main - neck"
      K_STORE=KF_IRR(1,1)

      DO NK=1,NKF
        NKPOID_N=NKPOID_N+1
c 
c   for Rb or other almost spherical fermisurface the next lines must be
c   commented out - they are necessary for copper, since the neck is
c   treated in a more accurate way
c
c -------------------------------------------------------------------
        IF (OPT('NECK    ')) THEN 
          WRITE(6,"((I5),(6e17.9))") NK,(KF_IRR(J,NK),J=1,3),K_STORE,
     +                          abs(K_STORE-KF_IRR(1,NK))
          IF (NK>1 .AND. NKPOID_B==0 .AND.
     +             abs(K_STORE-KF_IRR(1,NK)) > 5.D-2) then
            NKPOID_B=NK-1
            NKPOID_N=1
          END IF
          K_STORE=KF_IRR(1,NK)
        END IF
c -------------------------------------------------------------------

      END DO

      if (nkpoid_b==0) then
        nkpoid_n=0
        nkpoid_b=NKF           
      end if

c       nkpoid_total=nkpoid_total-1
c      NKPOID_N=NKPOID_N
c      NKPOID_B=NKPOID_B-1

      write (6,*) "Number of kpoints:", NKF          
      write (6,*) "Number of kpoints main:", nkpoid_b
      write (6,*) "Number of kpoints neck:", nkpoid_n


c read in the kpoints on the fermi surface (for copper only the main
c part)

      k_max=0
      do n=1,nkpoid_b
        do i=1,n
          k_max=k_max+1
          if ((i>1) .AND. (k_max>1)) THEN
            if ((abs(KF_IRR(3,k_max)-KF_IRR(3,k_max-1)) > 1.d-1)) then 
              k_max=k_max-1 
              exit
            end if
          end if
          if (k_max==nkpoid_b) exit
        end do
        n_max=n
        if (k_max==nkpoid_b) exit
      end do

      write (6,*) "n_max"          
 
c  i_max() contains the number of kpoints in one column

      allocate(i_max(n_max))

      k_max=0
      do n=1,nkpoid_b
        do i=1,n
c          write(6,*) "n,i:", n,i
          k_max=k_max+1
c          write(6,*) "kmax,k(k_max,3):", k_max, k(k_max,3) 
          if ((i>1) .AND. (k_max>1)) THEN
            if ((abs(KF_IRR(3,k_max)-KF_IRR(3,k_max-1)) > 1.d-1)) then 
c          if ((i>1) .AND. (k_max>1) .AND. 
c     +        (abs(KF_IRR(3,k_max)-KF_IRR(3,k_max-1)) > 1.d-1)) then 
              k_max=k_max-1 
              exit
            end if
          end if
          i_max(n)=i
          if (k_max==nkpoid_b) exit
        end do
c        write(6,*) "n,i_max(n):",n,i_max(n)
        if (k_max==nkpoid_b) exit
      end do

c      write(6,*) "n:",n
      write(6,*) "number of columns in the main part (n_max):",n_max
      write(6,*) "maximal value of the num of lines in the main part:"
     +                        ,maxval(i_max)

      allocate(mesh(n_max,n_max))
      allocate(border(nkpoid_b+nkpoid_n))
      allocate(border_tri(nkpoid_b+nkpoid_n,20))
      border=0
      border_tri=0

c for copper: read in the kpoints close to the neck   

      num_lines=0
      num_rows=0

      if(nkpoid_n>0 ) then

        num_lines=0
        k_max=nkpoid_b

        do n=1,nkpoid_n
          do i=1,nkpoid_n
            k_max=k_max+1
c            write(6,*) "n,i,k(k_max,3):", n,i, k(k_max,3) 
            IF ((i>1) .AND. 
     +         (abs(KF_IRR(3,k_max)-KF_IRR(3,k_max-1)) > 2.d-2)) THEN 
c            if ((i>1) .AND. (k_max>1) .AND. 
c     +        (abs(KF_IRR(3,k_max)-KF_IRR(3,k_max-1)) > 2.d-2)) then 
              k_max=k_max-1 
              exit
            end if
            num_lines(n)=i
            if (k_max==nkpoid_b+nkpoid_n) exit
          end do
          num_rows=n
c          write(6,*) "n,num_lines(n):",n,num_lines(n)
          if (k_max==nkpoid_b+nkpoid_n) exit
        end do
      end if

      write(6,*) "num of columns in the neck:", num_rows 
      write(6,*) "maximal value of the num of lines in the neck:",
     +                           maxval(num_lines)


       allocate(mesh_n(num_rows,maxval(num_lines)))
       allocate(trian(n_max**2+2*maxval(num_lines)*num_rows+num_rows*3,
     +                                                             3))
       allocate(area(n_max**2+2*maxval(num_lines)*num_rows+num_rows*3))


c construction of the triangles using the kpoints on the fermi surface


       WRITE(6,*) "before triangles" 

      CALL TRIANGLES(nkpoid_b,n_max,i_max,KF_IRR,border,border_tri,area,
     +   nkpoid_n,maxval(num_lines),num_rows,n_trian,mesh_n,trian)

      WRITE(6,*) "n_trian after triangles", n_trian
c calculation of the volume of the BZ, using the constructed triangles by conc necting them with the origin to form tetraeder 

c      DO n=1,NKF*48
c        WRITE(6,"((I5),(3e17.9))") n,(KF_COMPL(J,N),J=1,3)
c      END DO

      CALL TETRAEDER(nkpoid_b,nkpoid_n,mesh_n,KF_IRR,KF_COMPL,n_max,
     +   maxval(num_lines),num_rows,n_trian,trian)

c calculation of the area of the whole fermisurface
c at each kpoint the adjoining the triangles contribute with 1/3 of
c their area

      CALL AREA_FS(nkpoid_b,nkpoid_n,border,border_tri,area,n_max,
     +     maxval(num_lines),num_rows)


 9090 FORMAT(I5,3F17.9)

c construction of the kpoints above and below the fermi surface; they
c are situated on the surface normal at the distance delta_k_FS 
c important for the calculation of the Fermi velocity

c the surface normal is just calculated as the mean of the surface
c normal of the neighboring triangles


      write (6,*) "calculate the direction of the fermi velocity"

      CALL DIREC_VEL_FERMI(NKF,border,border_tri,KF_IRR,trian,
     +      n_max,num_rows,maxval(num_lines),direc_vel_FS(:,:))

c      CALL CREATE_CLUST_3D(NKF,direc_vel_FS(:,:,1),
c     +                     direc_vel_FS) 


c     DO NK=1,NKF
c       write(80, "((I5),(3e17.9))") NK,(DIREC_VEL_FS(J,NK),J=1,3)
c     END DO

c     write(80,*) " "
c     write(80,*) " "

      END SUBROUTINE CREATE_GRID_FS_3D
