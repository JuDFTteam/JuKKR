!**********************************************************************
subroutine GLL95(E,CLEB,ICLEB,LOFLM,IEND,TREFLL,DTREFLL,ATOM, &
                 REFPOT,RATOM,NATOM,ALAT,GREF0,DGDEOUT, &
                 LLY_G0TR,ICLS,CLS,I3, &
!                new input parameters after inc.p removal
                 naezd, lmaxd, naclsd, ncleb, nrefd, LLY)

! **********************************************************************
!
!     solution of the DYSON equation for a cluster of potentials
!     (TREFLL) centered at positions RATOM in free space,
!
! ----------------------------------------------------------------------

  implicit none

  integer :: naezd
  integer :: lmaxd
  integer :: naclsd
  integer :: ncleb
  integer :: nrefd
  integer :: LLY

  !
  !     INTEGER LMGF0D,NGD
  !     PARAMETER (LMGF0D= (LMAXD+1)**2,NGD=LMGF0D*NACLSD)
  !     INTEGER LLYNGD
  !     PARAMETER (LLYNGD=LLY*(LMGF0D*NACLSD-1)+1)

  double complex :: CONE
  double complex :: CZERO
  parameter (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
  !     ..
  !     .. Scalar Arguments ..
  double complex :: E
  double complex :: LLY_G0TR
  double precision :: ALAT
  integer :: I3
  integer :: IEND
  integer :: NATOM
  integer :: ICLS
  integer :: INFO
  !     ..
  !     .. Array Arguments ..
  integer :: IPVT(NACLSD*(LMAXD+1)**2)

  !     DOUBLE COMPLEX GREF0(NGD,LMGF0D)
  double complex :: GREF0(NACLSD*(LMAXD+1)**2,(LMAXD+1)**2)

  !     DOUBLE COMPLEX TREFLL(LMGF0D,LMGF0D,NREFD)
  double complex :: TREFLL((LMAXD+1)**2,(LMAXD+1)**2,NREFD)

  !     DOUBLE COMPLEX DTREFLL(LMGF0D,LMGF0D,NREFD)
  double complex :: DTREFLL((LMAXD+1)**2,(LMAXD+1)**2,NREFD)


  !                    DGDEOUT(LLYNGD,LMGF0D)
  double complex :: DGDEOUT(LLY*(NACLSD*(LMAXD+1)**2-1)+1, &
                            (LMAXD+1)**2)


  double precision :: CLEB(*)
  double precision :: RATOM(3,*)
  integer :: ATOM(*)
  integer :: ICLEB(NCLEB,3)
  integer :: LOFLM(*)
  integer :: REFPOT(*)
  integer :: CLS(NAEZD)
  !     ..
  !     .. Local Scalars ..
  integer :: I
  integer :: LM1
  integer :: LM2
  integer :: N1
  integer :: N2
  integer :: NDIM
  integer :: NLM1
  integer :: NLM2
  !     ..
  !     .. Local Arrays ..
  !     DOUBLE COMPLEX DGLLDE(LMGF0D,LMGF0D)
  !     DOUBLE COMPLEX GLL(LMGF0D,LMGF0D),GREF(NGD,NGD),GTREF(NGD,LMGF0D)

  double complex :: DGLLDE((LMAXD+1)**2,(LMAXD+1)**2)
  double complex :: GLL((LMAXD+1)**2,(LMAXD+1)**2)

! ---------------------------------------------------------------------
!     The following arrays can be very large (> 60 MB in typical cases)
!     therefore they are allocated on the heap
! ---------------------------------------------------------------------

  !     DOUBLE COMPLEX GREF (NACLSD*(LMAXD+1)**2,NACLSD*(LMAXD+1)**2)
  !     DOUBLE COMPLEX GTREF(NACLSD*(LMAXD+1)**2,(LMAXD+1)**2)

  double complex, allocatable, dimension(:,:) :: GREF
  double complex, allocatable, dimension(:,:) :: GTREF

  !     DOUBLE COMPLEX DGTDE(LLYNGD,LMGF0D)
  !     DOUBLE COMPLEX DGTDE(LLY*(NACLSD*(LMAXD+1)**2-1)+1,(LMAXD+1)**2)
  double complex, allocatable, dimension(:,:) ::  DGTDE

  !     DOUBLE COMPLEX DGTDE0(LLYNGD,LLYNGD)
  !     DOUBLE COMPLEX DGTDE0(LLY*(NACLSD*(LMAXD+1)**2-1)+1,
  !    &                      LLY*(NACLSD*(LMAXD+1)**2-1)+1)

  double complex, allocatable, dimension(:,:) :: DGTDE0

  !     DOUBLE COMPLEX DGDE(LLYNGD,LLYNGD)
  !     DOUBLE COMPLEX DGDE(LLY*(NACLSD*(LMAXD+1)**2-1)+1,
  !    &                    LLY*(NACLSD*(LMAXD+1)**2-1)+1)

  double complex, allocatable, dimension(:,:) :: DGDE


  double precision :: RDIFF(3)
  !     ..
  !     .. External Subroutines ..
  external GFREE,GREFSY,ZCOPY,ZGEMM

  !     .. External Functions ..
  logical :: TEST
  external TEST
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic ABS,DBLE
  !     ..

  integer :: memory_stat
  logical :: memory_fail

  integer :: LMGF0D
  integer :: NGD
  integer :: LLYNGD


  LMGF0D = (LMAXD+1)**2
  NGD = LMGF0D*NACLSD
  LLYNGD = LLY*(LMGF0D*NACLSD-1) + 1

  !     Allocate arrays
  memory_stat = 0
  memory_fail = .false.

  allocate(GREF(NGD,NGD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  allocate(GTREF(NGD,LMGF0D), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  allocate(DGTDE(LLYNGD,LMGF0D), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  allocate(DGTDE0(LLYNGD,LLYNGD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  allocate(DGDE(LLYNGD,LLYNGD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (memory_fail == .true.) then
    write(*,*) "GLL95: FATAL Error, failure to allocate memory."
    write(*,*) "       Probably out of memory."
    stop
  end if

  if (TEST('flow    ')) write (6,fmt=*) '>>> GLL95'

  NDIM = LMGF0D*NATOM

  !
  ! ---> construct free Green's function
  !
  do N1 = 1,NATOM
    do N2 = 1,NATOM
      do I = 1,3
        !            RDIFF(I) = (RATOM(I,N1) - RATOM(I,N2))*ALAT
        !           changed P.Z. 4.7.97
        RDIFF(I) = - (RATOM(I,N1)-RATOM(I,N2))*ALAT
      end do

      if (N1/=N2) then

        call GFREE(RDIFF,E,GLL,CLEB,ICLEB,LOFLM,IEND, lmaxd, ncleb)

        do LM2 = 1,LMGF0D
          NLM2 = (N2-1)*LMGF0D + LM2
          do LM1 = 1,LMGF0D
            NLM1 = (N1-1)*LMGF0D + LM1
            GREF(NLM1,NLM2) = GLL(LM1,LM2)
          end do
        end do
      else
        do LM2 = 1,LMGF0D
          NLM2 = (N2-1)*LMGF0D + LM2
          do LM1 = 1,LMGF0D
            NLM1 = (N1-1)*LMGF0D + LM1
            GREF(NLM1,NLM2) = CZERO
          end do
        end do

      end if  ! (N1 /= N2)
    end do
   end do

   if (TEST('flow    ')) write (6,fmt=*) 'GFREE o.k.'
   ! ----------------------------------------------------------------------

   if (LLY==1) then
     !
     ! ---> construct derivative of free Green's function
     !
     do N1 = 1,NATOM
       do N2 = 1,NATOM
         do I = 1,3
           RDIFF(I) = - (RATOM(I,N1)-RATOM(I,N2))*ALAT
         end do

         if (N1/=N2) then

           call DGFREE(RDIFF,E,DGLLDE,CLEB,ICLEB,LOFLM,IEND, &
           lmaxd, ncleb)

           do LM2 = 1,LMGF0D
             NLM2 = (N2-1)*LMGF0D + LM2
             do LM1 = 1,LMGF0D
               NLM1 = (N1-1)*LMGF0D + LM1
               DGDE(NLM1,NLM2) = DGLLDE(LM1,LM2)
             end do
           end do
         else
           do LM2 = 1,LMGF0D
             NLM2 = (N2-1)*LMGF0D + LM2
             do LM1 = 1,LMGF0D
               NLM1 = (N1-1)*LMGF0D + LM1
               DGDE(NLM1,NLM2) = CZERO
             end do
           end do

         end if  ! (N1 /= N2)

       end do
     end do

   endif


   call ZCOPY(NGD*LMGF0D,GREF,1,GREF0,1)

   if (LLY==1) then

     do N2 = 1,NATOM
       NLM2 = (N2-1)*LMGF0D + 1

       call ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,DGDE(1,NLM2),NGD, &
                  TREFLL(1,1,REFPOT(ABS(ATOM(N2)))),LMGF0D, &
                  CZERO,GTREF,NGD)

       call ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,GREF(1,NLM2),NGD, &
                  DTREFLL(1,1,REFPOT(ABS(ATOM(N2)))),LMGF0D, &
                  CONE,GTREF,NGD)

       call ZCOPY(NGD*LMGF0D,GTREF,1,DGTDE0(1,NLM2),1)
     end do
   end if

   do N2 = 1,NATOM
     NLM2 = (N2-1)*LMGF0D + 1

     call ZGEMM('N','N',NDIM,LMGF0D,LMGF0D,-CONE,GREF(1,NLM2),NGD, &
     TREFLL(1,1,REFPOT(ABS(ATOM(N2)))),LMGF0D, &
     CZERO,GTREF,NGD)

     call ZCOPY(NGD*LMGF0D,GTREF,1,GREF(1,NLM2),1)
   end do

   if (LLY==1) then
     do N2 = 1,LMGF0D
       do N1 = 1,NGD
         DGTDE(N1,N2) = DGTDE0(N1,N2)
       enddo
     enddo
   endif


   call GREFSY(GREF,GREF0,IPVT,NDIM,DGTDE, &
               LLY_G0TR, &
               naezd, lmaxd, naclsd, LLY)

   if (LLY==1) then

     call ZGEMM('N','N',NDIM,LMGF0D,NDIM,-CONE,DGTDE0,NGD, &
                GREF0,NGD,CONE,DGDE,NGD)

     call ZGETRS('N',NDIM,LMGF0D,GREF,NGD,IPVT,DGDE,NGD,INFO)

     do N2 = 1,LMGF0D
       do N1 = 1,NGD
         DGDEOUT(N1,N2) = DGDE(N1,N2)
       enddo
     enddo

   endif

   !     Deallocate arrays

   deallocate(GREF)
   deallocate(GTREF)
   deallocate(DGTDE)
   deallocate(DGTDE0)
   deallocate(DGDE)

   if (TEST('flow    ')) write (6,fmt=*) 'GREFSY o.k.'

 end subroutine GLL95

