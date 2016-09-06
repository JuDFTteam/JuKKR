subroutine ADDVIRATOMS14(LINTERFACE,NVIRT,naez, naezd, natypd,NEMB,NEMBD,RBASIS, &        
                       LCARTESIAN,BRAVAIS,NCLS,NINEQ,REFPOT,KAOEZ,NOQ,NREF,RMTREFAT,I25) 

      implicit none
!interface variables
      logical             :: LINTERFACE,LCARTESIAN,LABSCORD
      integer             :: naez,NCLS,NINEQ
      integer             :: NAEZD,NEMBD
      integer             :: NAtypd,NREF,NVIRT
      integer             :: NEMB,NACLSD
      integer             :: IVIR(1000)
      parameter (NACLSD=1000)

      integer             :: I,I1,J
      real*8              :: RBASIS(3,*),RBASISOLD(3,NEMB+NAEZD),RBASISSAVE(3,NEMB+NAEZD)
      real*8              :: RMTREFAT(NAEZD+NEMBD)
      integer             :: REFPOT(*),REFPOTOLD(NAEZD+NEMB)
      integer             :: NOQ(*)
      INTEGER             :: KAOEZ(NATYPD,*),KAOEZOLD(1,NEMB+NAEZD)
      real*8              :: DIFF,RMAXCLUS,VEC1(3),VEC2(3,NACLSD)
      integer             :: NBR(3),NMAX,NMAXZ,N1,N2,N3,IQ
      external            :: GETCLUSNXYZ
    


!local variables
      CHARACTER(len=40)   :: I25
      integer             :: NREFOLD
      integer             :: natomimp
      real*8,allocatable  :: ratomimp(:,:)
!      real*8,allocatable  :: rbasislist(:,:)
      integer,allocatable :: ATOMIMP(:)
      integer             :: iatom,ibasis
      integer             :: ierr
!       real*8              :: ratomvtest(3)
      real*8              :: rbasisnew(3)
      real*8,allocatable  :: RCLSNEW(:,:),RBASISNEW1(:,:)

      real*8              :: bravais(3,3)
      real*8,allocatable  :: bravaisinv(:,:)
      real*8              :: tol
      integer             :: ndim
      integer             :: naeznew

      tol = 1.d-5


      write(1337,*) 'LINTERFACE',LINTERFACE
      write(1337,*) 'NAEZ',NAEZ
      write(1337,*) 'NAEZD',NAEZD
      write(1337,*) 'NEMB',NEMB
      write(1337,*) 'RBASISOLD'
      do ibasis=1,naez+nemb
        write(1337,*) ibasis,rbasis(:,ibasis)
        REFPOTOLD(IBASIS)=REFPOT(IBASIS)
        KAOEZOLD(1,IBASIS)=KAOEZ(1,IBASIS)
        RBASISSAVE(:,IBASIS)=RBASIS(:,IBASIS)
      end do



!  -----------------------------------------------------------------
!                      read the scoef file
!  -----------------------------------------------------------------

      open(unit=32452345,file=I25,iostat=ierr)
      write(1337,*) '*',I25,'*'
      write(1337,*) '*',Ierr,'*'
      if (ierr/=0) stop '[addvirtual] file not found'
      read(32452345,*) natomimp
      write(1337,*) 'natomimp',natomimp
      allocate(ratomimp(3,natomimp))
      allocate(ATOMIMP(natomimp))
      do iatom=1,natomimp
        read(32452345,*) ratomimp(:,iatom),ATOMIMP(iatom)
        write(1337,'(A,I,A,3F)') 'IMPATOM ',iatom,' :',ratomimp(:,iatom)
      end do

!  -----------------------------------------------------------------
!                      set bulk/surface
!  -----------------------------------------------------------------
      IF ( LINTERFACE ) THEN
         NDIM = 2
         WRITE (1337,'(23X,A)') 'ADDVIRTUAL : surface geometry mode'
      ELSE
         NDIM = 3
         WRITE (1337,'(23X,A)') 'ADDVIRTUAL : bulk geometry mode'
      END IF
!  -----------------------------------------------------------------
!                      read bravais vectors
!  -----------------------------------------------------------------
      
! invert bravais vectors
      ALLOCATE(BRAVAISINV(NDIM,NDIM))
      BRAVAISINV = BRAVAIS(1:NDIM,1:NDIM)
      call INVERSE_D1(BRAVAISINV)
  

      NREFOLD=0
      DO IBASIS=1,NAEZ+NEMB
        NREFOLD = MAX(NREFOLD,REFPOTOLD(IBASIS)) 
      ENDDO
      write(1337,*) 'Number of reference potentials is currently', nrefold

! Change basis vectors to cartesian coordinates in order to calculate distances.
      IF (.NOT.LCARTESIAN) THEN

         IF (LINTERFACE) THEN
            DO I=1,NAEZ+NEMB
               DO J=1,NDIM
                  RBASISOLD(J,I)= (RBASIS(1,I)*BRAVAIS(J,1)+RBASIS(2,I)*BRAVAIS(J,2))
               ENDDO
               RBASISOLD(3,I)=RBASIS(3,I)
            ENDDO
         ELSE
            DO I=1,NAEZ+NEMB
               DO J=1,NDIM
                  RBASISOLD(J,I)= (RBASIS(1,I)*BRAVAIS(J,1)+RBASIS(2,I)*BRAVAIS(J,2)+RBASIS(3,I)*BRAVAIS(J,3))
               ENDDO
            ENDDO
         ENDIF

      ELSE

         DO I=1,NAEZ+NEMB
            RBASISOLD(:,I)=RBASIS(:,I)
         ENDDO

      ENDIF


      ! If the 1st imp. atom in the list is at (0,0,0) then all coordinates are assumed
      ! relative to the 1st imp atom, otherwise relative to the lattice coords (absolute coords).
      LABSCORD = .FALSE.
      J = 0
      DO J = 1,3
         IF ( ABS(RATOMIMP(J,1)).GT.1D-8 ) LABSCORD = .TRUE.
      END DO

      ALLOCATE(RCLSNEW(3,NATOMIMP))
      ALLOCATE(RBASISNEW1(3,NATOMIMP))
      DO I = 1,NATOMIMP
         CALL DCOPY(3,RATOMIMP(1,I),1,RCLSNEW(1,I),1)
      END DO
      IF ( .NOT.LABSCORD ) THEN
         IQ = ATOMIMP(1)
         DO I = 1,NATOMIMP
            CALL DAXPY(3,1D0,RBASISOLD(1,IQ),1,RCLSNEW(1,I),1)
         END DO
      END IF
      RMAXCLUS = 0D0
      DO I = 2,NATOMIMP
         DIFF = 0D0
         DO J = 1,3
            DIFF = DIFF + (RCLSNEW(J,I)-RCLSNEW(J,1))**2
         END DO
         DIFF = SQRT(DIFF)
         RMAXCLUS = MAX(RMAXCLUS,DIFF)
      END DO

      NBR(1:3) = 0
      CALL GETCLUSNXYZ(RMAXCLUS,BRAVAIS,NDIM,DIFF,NBR)
      NMAX = MAX(NBR(1),NBR(2),NBR(3))
      NMAXZ = NMAX 
      IF ( NDIM.EQ.2 ) NMAXZ = 0
      IQ=0
      DO N1 = -NMAX,NMAX
         DO N2 = -NMAX,NMAX
            DO N3 = -NMAXZ,NMAXZ
               
               VEC1(1:3) = DBLE(N1)*BRAVAIS(1:3,1)+DBLE(N2)*BRAVAIS(1:3,2) + DBLE(N3)*BRAVAIS(1:3,3)

               DO I1 = 1,NAEZ
                  IQ=IQ+1
                  DIFF = 0D0
                  VEC2(1:3,IQ) = VEC1(1:3) + RBASISOLD(1:3,I1)
               END DO

            END DO
         END DO
      END DO

      IBASIS=0
      DO 100 I=1,NATOMIMP
         DO I1=1,IQ
            DIFF = SQRT((RCLSNEW(1,I)-VEC2(1,I1))**2+(RCLSNEW(2,I)-VEC2(2,I1))**2+(RCLSNEW(3,I)-VEC2(3,I1))**2)
            IF ( DIFF.LE.(TOL) ) GOTO 100 ! Position is on lattice, do not set as virtual atom          
         ENDDO
         CALL RTOBASIS(BRAVAIS,RCLSNEW(:,I),RBASISNEW,NDIM)
         IF (LINTERFACE) THEN
            DO J=1,2
               RBASISNEW1(J,I)=RBASISNEW(1)*BRAVAISINV(J,1)+RBASISNEW(2)*BRAVAISINV(J,2)
            ENDDO
            RBASISNEW1(3,I)=RBASISNEW(3)
         ELSE
            RBASISNEW1(1:3,I)=RBASISNEW(1)*BRAVAISINV(1:3,1)+RBASISNEW(2)*BRAVAISINV(1:3,2)+RBASISNEW(3)*BRAVAISINV(1:3,3)
         ENDIF
         WRITE(1337,*) 'rnew',RBASISNEW1(:,I)
         IF (I.GT.1) THEN
            DO I1=1,I-1
               DIFF=SQRT((RBASISNEW1(1,I)-RBASISNEW1(1,I1))**2+(RBASISNEW1(2,I)-RBASISNEW1(2,I1))**2+(RBASISNEW1(3,I)-RBASISNEW1(3,I1))**2)
               IF (DIFF.LE.1D-05) GOTO 100
            ENDDO
         ENDIF
         IBASIS=IBASIS+1
         IVIR(IBASIS)=I
100   CONTINUE
      ! IBASIS is the number of virtual atoms
      WRITE(1337,*) 'ibasis',IBASIS,(IVIR(J),J=1,IBASIS)

      if (ibasis+naez>naezd) then
         write(*,*) '[addvirtual] naez increased to ',ibasis
         write(*,*) '[addvirtual] naezd is',naezd
         write(*,*) '[addvirtual] naeznew > naezd please change naezd'
         stop 'addvistual'
      else
         write(1337,*) 'NAEZ will soon be increased to : ',ibasis+naez
      end if


      NVIRT = IBASIS
      NINEQ=NINEQ+NVIRT
      NCLS=NCLS+NVIRT
      naeznew = NVIRT+naez
      if (naeznew>naezd) then
         write(*,*) '[addvirtual] naez increased to ',naeznew
         write(*,*) '[addvirtual] naezd is',naezd
         write(*,*) '[addvirtual] naeznew > naezd please change naezd'
         stop 'addvirtual'
      end if

      do i=naeznew+1,naeznew+nemb ! Added +1 : Phivos 25.7.2014
         RBASIS(:,I)=RBASISSAVE(:,I-IBASIS)
         REFPOT(I)=REFPOTOLD(I-IBASIS)
         KAOEZ(1,I)=KAOEZOLD(1,I-IBASIS)
      ENDDO
      DO I = NAEZNEW + NEMB , NAEZNEW + 1, -1
         RMTREFAT(I) = RMTREFAT(I-NVIRT)  ! Shift values of embedded positions in array
      ENDDO
      DO I = NAEZNEW - NEMB, NAEZNEW
         RMTREFAT(I) = 1.D-20
      ENDDO
      

      DO I=NAEZ+1,NAEZNEW
         RBASIS(:,I)=RBASISNEW1(:,IVIR(I-NAEZ))
         REFPOT(I)=NREFOLD+1
         NOQ(I)=0
         KAOEZ(1,I)=-1
      ENDDO
      DO I=1,NAEZNEW+NEMB
         NREF=MAX(NREF,REFPOT(I))
      ENDDO

!  -----------------------------------------------------------------
!     write out stuff
!  -----------------------------------------------------------------
    write(1337,*) 'addvirtual: List of new basis atoms including virtual atoms'
    do j = 1, naeznew+nemb
      write(1337,*) rbasis(:,j)
    end do
    write(1337,*) '-------------------------------------------------'

     write(1337,*) 'naeznew is now ',naeznew
     write(1337,*) 'setting naez to naeznew'
     naez=naeznew
     write(1337,*) 'updating rbasis array with virtual basis sites'



    do ibasis = 1, naeznew+nemb
       write(1337,*) 'REFPOT',REFPOT(ibasis)
       write(1337,*) 'NOQ',   NOQ(ibasis)
       write(1337,*) 'KAOEZ', KAOEZ(1,ibasis)
    end do









!    stop 'end of ADDVIRTUAL'
!   deallocate(rbasislist,ratomimp,atomimp)
   deallocate(ratomimp,atomimp)
   DEALLOCATE(BRAVAISINV)
   DEALLOCATE(RCLSNEW)
   DEALLOCATE(RBASISNEW1)
 contains

logical function vec_in_list(vec,veclist,bound)
  ! --------------------------
  ! checks if the vector vec is in the vector list veclist
  ! in the range of (1,bound)
  ! --------------------------
  integer        :: bound
  real*8         :: vec(3)
  real*8         :: veclist(3,bound)
  integer        :: ilist
  real*8         :: tempvec(3),diff
  vec_in_list=.false.
  do ilist = 1, bound
    tempvec = vec-veclist(:,ilist)
    diff = sqrt(tempvec(1)**2+tempvec(2)**2+tempvec(3)**2)
    if (diff<10e-5) vec_in_list=.true.
  end do
end function vec_in_list


subroutine rtobasis(bravais,rpos,rbasis,ndim)
  ! --------------------------
  ! converts a spacial vector rpos to a basis vector rbasis
  ! such that rbasis = bravais * n with n in [0,1]^ndim
  ! --------------------------
  implicit none
  real*8              :: bravais(3,3)
  real*8              :: bravais_inv(ndim,ndim)
  integer             :: ndim
  real*8              :: rpos(3)
  real*8              :: rbasis(3)
  
  real*8              :: ncoeffreal(ndim)
  integer             :: ncoeffint(ndim)
  integer             :: idim

  ! --------------------------
  ! first invert the bravais matrix => bravais_inv
  ! --------------------------
     bravais_inv = bravais(1:ndim,1:ndim)
     call inverse_d1(bravais_inv)
  ! --------------------------
  ! then do n = Bravais* rpos
  ! --------------------------
     ncoeffreal=0.0D0
     do idim=1,ndim
     ncoeffreal = ncoeffreal + bravais_inv(1:ndim,idim)*rpos(idim)
     end do

  ! ------old method ---------
  ! take the smaller integer value of n => ncoeffint
  ! --------------------------
        do idim = 1, ndim
  !       ncoeffint(idim) = floor(ncoeffreal(idim))
        end do

  ! ------new method---------
  ! take the smaller integer value of n + 0.5 => ncoeffint
  ! --------------------------
     do idim = 1, ndim
       ncoeffint(idim) = floor(ncoeffreal(idim)+0.5)
     end do
  ! --------------------------
  ! do rbasis = rpos - bravias * ncoeffint
  ! --------------------------
     rbasis = 0.0D0
     do idim = 1, ndim
     rbasis = rbasis + bravais(:,idim)*ncoeffint(idim)
     end do
     rbasis = rpos - rbasis
end subroutine rtobasis

subroutine inverse_d1(mat)
  implicit none
  real(8), intent(inout) :: mat(:,:)
  real(8), allocatable   :: work(:)
  integer, allocatable   :: ipiv(:)
  integer                :: n,info
  n = size(mat,1)
  if(size(mat,2).ne.n) stop 'inverse_d1: array dimensions differ.'
  allocate ( ipiv(n),work(n) )
  call dgetrf(n,n,mat,n,ipiv,info)      ; if(info.ne.0) stop 'inverse_d1: dpotrf failed.'
  call dgetri(n,mat,n,ipiv,work,n,info) ; if(info.ne.0) stop 'inverse_d1: dpotri failed.'
end subroutine inverse_d1


end subroutine ADDVIRATOMS14
