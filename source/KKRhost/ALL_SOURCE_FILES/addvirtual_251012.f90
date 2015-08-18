subroutine ADDVIRATOMS(LINTERFACE, naez, naezd, natypd,NEMB,RBASIS, &        
                        REFPOT, KAOEZ,NOQ,NREF) 

      implicit none
!interface variables
      logical             :: LINTERFACE,LCARTESIAN,LABSCORD
      integer             :: naez
      integer             :: NAEZD
      integer             :: NAtypd,NREF
      integer             :: NEMB,NACLSD
      parameter (NACLSD=1000)

      integer             :: I,I1,J
      real*8              :: RBASIS(3,*)
      real*8              :: DIFF,RMAXCLUS,VEC1(3),VEC2(3,NACLSD)
      integer             :: NBR(3),NMAX,NMAXZ,N1,N2,N3,IQ
      external            :: GETCLUSNXYZ
    

      integer             :: REFPOT(*)
!       integer             :: ntcell(*)
!       double precision    :: mtfac(*)
!       integer             :: irns(*) 
!       double precision    :: rmtref(*)
! 
!       INTEGER             :: IQAT(NAEZD,NATYPD)
!       integer             :: cls(*)
!       real*8              :: CONC(NATYPD)
      integer             :: NOQ(*)
!       integer             :: nat(natypd)
      INTEGER             :: KAOEZ(NATYPD,*) 

!local variables
      integer             :: first
      CHARACTER(len=40)   :: I25
      integer             :: NREFOLD
      integer             :: natomimp
      real*8,allocatable  :: ratomimp(:,:)
!      real*8,allocatable  :: rbasislist(:,:)
      integer,allocatable :: ATOMIMP(:)
      integer             :: iatom,ibasis,idim1,idim2
      integer             :: ierr
!       real*8              :: ratomvtest(3)
      real*8              :: rbasisnew(3)
      real*8,allocatable  :: RCLSNEW(:,:)

      real*8              :: bravais(3,3)
      integer             :: ndim
      integer             :: naeznew

      CHARACTER(len=80)   :: UIO
      write(*,*) 'LINTERFACE',LINTERFACE
      write(*,*) 'naez',naez
      write(*,*) 'naezd',naezd
      do iatom=1,naez+nemb
        write(*,'(A,I,A,3F)') 'rbasisold',iatom,':',rbasis(:,iatom)
      end do

    do ibasis = 1, naezd+nemb
write(*,*) '---------------------------------------------'
write(*,*) ibasis
write(*,*) '---------------------------------------------'
write(*,*) 'REFPOT',         REFPOT(ibasis)
write(*,*) 'KAOEZ',        KAOEZ(1,ibasis)
    end do
    do ibasis = 1,naezd
write(*,*) 'NOQ',         NOQ(ibasis)





!  -----------------------------------------------------------------
!                      read the scoef file
!  -----------------------------------------------------------------
      CALL IoInput('FILES     ',UIO,1+4,7,IERR)
                      READ (UNIT=UIO,FMT='(A40)')  I25


      open(unit=32452345,file=I25,iostat=ierr)
      write(*,*) '*',I25,'*'
      write(*,*) '*',Ierr,'*'
      if (ierr/=0) stop '[readvirtual] file not found'
      read(32452345,*) natomimp
      write(*,*) 'natomimp',natomimp
      allocate(ratomimp(3,natomimp))
      allocate(ATOMIMP(natomimp))
      do iatom=1,natomimp
        read(32452345,*) ratomimp(:,iatom),ATOMIMP(iatom)
        write(*,'(A,I,A,3F)') 'IMPATOM ',iatom,' :',ratomimp(:,iatom)
      end do

!  -----------------------------------------------------------------
!                      set bulk/surface
!  -----------------------------------------------------------------
      IF ( LINTERFACE ) THEN
         NDIM = 2
         WRITE (6,'(23X,A)') 'ADDVIRTUAL : surface geometry mode'
      ELSE
         NDIM = 3
         WRITE (6,'(23X,A)') '  ADDVIRTUAL : bulk geometry mode'
      END IF
!  -----------------------------------------------------------------
!                      read bravais vectors
!  -----------------------------------------------------------------
      DO IDIM1 = 1,NDIM
         CALL IOINPUT('BRAVAIS   ',UIO,IDIM1,7,IERR)
         READ (UNIT=UIO,FMT=*) (BRAVAIS(IDIM2,IDIM1),IDIM2=1,NDIM)
         WRITE(*,'(A,I,A,3F)') 'BRAVAIS',IDIM1,': ',BRAVAIS(:,IDIM1)
      END DO

!      bravais(1,:)=(/ 1.0D0, 0.0D0, 0.0D0 /)
!      bravais(2,:)=(/ 0.0D0, 1.0D0, 0.0D0 /)
!      bravais(3,:)=(/ 0.0D0, 0.0D0, 1.0D0 /)


!      rbasisold (:,1)= (/ -0.9D0, -0.9D0, -0.9D0 /)
!      rbasisold (:,2)= (/ 0.15D0, 0.15D0, 0.15D0 /)
!      rbasisold (:,3)= (/ 0.2D0, 0.2D0, 0.2D0 /)
!      naez=3



!      if (bravais(1,3)**2+bravais(2,3)**2+bravais(3,3)**2<10e-10) then
!        ndim = 2
!      else 
!        ndim = 3
!      end if


      NREFOLD=0
      DO IBASIS=1,NAEZ
        NREFOLD = MAX(NREFOLD,REFPOT(IBASIS)) 
      ENDDO
      write(*,*) 'Number of reference potentials is currently', nrefold

      LABSCORD = .FALSE.
      J = 0
      DO WHILE ( J.LT.3.AND..NOT.LABSCORD )
         J = J + 1
         IF ( ABS(RATOMIMP(J,1)).GT.1D-8 ) LABSCORD = .TRUE.
      END DO
      ALLOCATE(RCLSNEW(3,NATOMIMP))
      DO I = 1,NATOMIMP
         CALL DCOPY(3,RATOMIMP(1,I),1,RCLSNEW(1,I),1)
      END DO
      IF ( .NOT.LABSCORD ) THEN
         IQ = ATOMIMP(1)
         DO I = 1,NATOMIMP
            CALL DAXPY(3,1D0,RBASIS(1,IQ),1,RCLSNEW(1,I),1)
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

      DO I = 1,3
         NBR(I) = 0
      END DO
      CALL GETCLUSNXYZ(RMAXCLUS,BRAVAIS,NDIM,DIFF,NBR)
      NMAX = MAX(NBR(1),NBR(2),NBR(3))
      NMAXZ = NMAX 
      IF ( NDIM.EQ.2 ) NMAXZ = 0
      IQ=0
         DO N1 = -NMAX,NMAX
            DO N2 = -NMAX,NMAX
               DO N3 = -NMAXZ,NMAXZ

                  DO J = 1,3
                     VEC1(J) = DBLE(N1)*BRAVAIS(J,1)+DBLE(N2)*BRAVAIS(J,2) + DBLE(N3)*BRAVAIS(J,3)
                  END DO

                  DO I1 = 1,NAEZ
                     IQ=IQ+1
                     DIFF = 0D0
                     DO J = 1,3
                        VEC2(J,IQ) = VEC1(J) + RBASIS(J,I1)
                     END DO
                  END DO
               END DO
            END DO
         END DO
       IBASIS=0
       DO 100 I=1,NATOMIMP
        DO I1=1,IQ
         DIFF = SQRT((RCLSNEW(1,I)-VEC2(1,I1))**2+(RCLSNEW(2,I)-VEC2(2,I1))**2+(RCLSNEW(3,I)-VEC2(3,I1))**2)
         IF ( DIFF.LE.(1D-5) ) THEN
          GOTO 100
         ENDIF
        ENDDO
        IBASIS=IBASIS+1
       RBASIS(:,NAEZ+IBASIS)=RCLSNEW(:,I)
 100  CONTINUE
      WRITE(6,*) 'ibasis',IBASIS
!  -----------------------------------------------------------------
!             redefine existing basis vectors such that
!             they satisfy:
!  -----------------------------------------------------------------
!             rbasis = sum_i n(i)* BRAVAIS(i) with n in [0,1]
!  -----------------------------------------------------------------
!     allocate(rbasislist(3,naez+natomimp))
!     do ibasis = 1, naez
!       call rtobasis(bravais,rbasisold(:,ibasis),rbasisnew,ndim)
!       rbasislist(:,ibasis) = rbasisnew
!     end do!iatom

!  -----------------------------------------------------------------
!     .... (to be commented)
!  -----------------------------------------------------------------
!     ibasis=naez
!     do iatom=1,natomimp
!       call rtobasis(bravais,ratomimp(:,iatom),rbasisnew,ndim)
!       write(*,'(A,I4,A,3F14.8,A,3F14.8)') 'atom',iatom,'(',ratomimp(:,iatom),') transformed to (',rbasisnew
!       if (.not. vec_in_list(rbasisnew,rbasislist,ibasis)) then
!         ibasis = ibasis + 1
!       end if
!     end do !natomimp
!     write(*,*) 'rbasis_test'
!     do iatom = 1,naez
!     write(*,FMT='(i4,3f16.8)') iatom,rbasislist(1,iatom),rbasislist(2,iatom),rbasislist(3,iatom)
!     write(*,FMT='(i4,3f16.8)') iatom,rbasisold(1,iatom),rbasisold(2,iatom),rbasisold(3,iatom)
!     enddo
     if (ibasis+naez>naezd) then
       write(*,*) '[readvirtual] naez increased to ',ibasis
       write(*,*) '[readvirtual] naezd is',naezd
       write(*,*) '[readvirtual] naeznew > naezd please change naezd'
       stop
     else
       write(*,*) 'NAEZ will soon be increased to : ',ibasis+naez
     end if


!     ibasis=naez+1
!     do iatom=1,natomimp
!        write(*,*) '*'
!       call rtobasis(bravais,ratomimp(:,iatom),rbasisnew,ndim)
!       if (.not. vec_in_list(rbasisnew,rbasislist,ibasis)) then
!         rbasislist(:,ibasis) = rbasisnew

!         REFPOT(ibasis)=NREFOLD+1
!         NOQ(ibasis) = 0
!         KAOEZ(1,ibasis) = -1 !ibasis
!         ibasis = ibasis + 1
!         write(*,*) iatom,' added'
!       else
!       write(*,*) iatom,'already in basis'
!       end if
!     end do !natomimp
     naeznew = ibasis+naez
     do i=naez+1,naeznew
       refpot(i)=NREFOLD+1
       NOQ(i)=0
       KAOEZ(1,i)=-1
     enddo 
     DO i=1,naeznew
       NREF=MAX(NREF,REFPOT(I))
     ENDDO
     if (naeznew>naezd) then
       write(*,*) '[readvirtual] naez increased to ',naeznew
       write(*,*) '[readvirtual] naezd is',naezd
       write(*,*) '[readvirtual] naeznew > naezd please change naezd'
       stop
     end if
        

!  -----------------------------------------------------------------
!     write out stuff
!  -----------------------------------------------------------------
    write(*,*) 'List of new basis atoms'
    do j = 1, naeznew
      write(*,*) rbasis(:,j)
    end do
    write(*,*) '-------------------------------------------------'

     write(*,*) 'naeznew is now ',naeznew
     write(*,*) 'setting naez to naeznew'
     naez=naeznew
     write(*,*) 'updating rbasis array with virtual basis sites'
!     rbasisold=0.0D0
!     rbasisold = rbasislist(:,1:naeznew)

!      ratomvtest = (/ 1.5D0,1.5D0,1.5D0/)

!      call rtobasis(bravais,ratomvtest,rbasis,ndim)
!      write(*,*) 'RBASIS is ',rbasis



    do ibasis = 1, naeznew
write(*,*) '---------------------------------------------'
write(*,*) ibasis
write(*,*) '---------------------------------------------'
! write(*,*) 'LMXC',LMXC(ibasis)
! write(*,*) 'KFG',         KFG(:,ibasis)
write(*,*) 'REFPOT',         REFPOT(ibasis)
! write(*,*) 'NTCELL',         NTCELL(ibasis)
! write(*,*) 'MTFAC',         MTFAC(ibasis)
! write(*,*) 'IRNS',         IRNS(ibasis)
! write(*,*) 'RMTREF',         RMTREF(REFPOT(ibasis))
! 
! write(*,*) 'IQAT',         IQAT(1,ibasis) 
! write(*,*) 'CLS',         CLS(ibasis)
! write(*,*) 'CONC',         CONC(ibasis)
write(*,*) 'NOQ',         NOQ(ibasis)
! write(*,*) 'NAT',        NAT(ibasis)
write(*,*) 'KAOEZ',        KAOEZ(1,ibasis)
! write(*,*) '---------------------------------------------'
    end do









!    stop 'end of ADDVIRTUAL'
!   deallocate(rbasislist,ratomimp,atomimp)
   deallocate(ratomimp,atomimp)
   DEALLOCATE(RCLSNEW)
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
  integer                :: n,info,i,j
  n = size(mat,1)
  if(size(mat,2).ne.n) stop 'inverse_d1: array dimensions differ.'
  allocate ( ipiv(n),work(n) )
  call dgetrf(n,n,mat,n,ipiv,info)      ; if(info.ne.0) stop 'inverse_d1: dpotrf failed.'
  call dgetri(n,mat,n,ipiv,work,n,info) ; if(info.ne.0) stop 'inverse_d1: dpotri failed.'
end subroutine inverse_d1


end subroutine ADDVIRATOMS
