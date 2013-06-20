module ShapeStandardMesh_mod
  implicit none

  public :: MESH
  private :: MESH0

contains

  !----------------------------------------------------------------------------
  !> Radial mesh generation for non-spherical part.
  !> @brief
  !>
  !-----------------------------------------------------------------------
  !>     THIS ROUTINE DEFINES A UNIQUE SUITABLE RADIAL MESH 'XRN,DRN' OF
  !>     'MESHN' POINTS,DISTRIBUTED INTO 'NPAN' PANNELS  DEFINED  BY THE
  !>     CRITICAL POINTS 'CRT'
  !-----------------------------------------------------------------------
  !> @param[in,out] CRT    array of critical points, dimension NPAND, unsorted on input, sorted on output
  !> @param         NPAN   number of critical points
  !> @param[in,out] NM     number of points in panels/critical intervals
  !> @param         XRN    radial mesh points
  !> @param         DRN    weights for radial mesh
  !> @param[out]    MESHN  number of mesh points generated
  !> @param         NPOI   [in]? number of mesh points requested ???
  !> @param[in]     KEYPAN set to 0 to generate mesh
  !> @param[in]     NMIN   minimal number of points for panels/critical intervals
  !> @param[in]     MESHND dimension of CRT, XRN and DRN arrays
  !> @param[in]     NPAND  dimension of NM array
  !> @param[in]     VERBOSITY  0=no screen output, 1=output info about radial mesh
  subroutine MESH(CRT,NPAN,NM,XRN,DRN,MESHN,NPOI,KEYPAN,NMIN, &
  MESHND,NPAND,VERBOSITY) ! new input parameters

    use shape_constants_mod, only: DP
    implicit none

    !     .. PARAMETER STATEMENTS ..

    integer, intent(in) :: MESHND
    integer, intent(in) :: NPAND
    integer, intent(in) :: VERBOSITY
    !     .. SCALAR ARGUMENTS ..

    integer ::   NPAN,MESHN,NPOI,KEYPAN,NMIN

    !     .. ARRAY ARGUMENTS ..

    integer ::   NM(NPAND)
    real(kind=DP) ::    CRT(NPAND),XRN(MESHND),DRN(MESHND)

    !     .. LOCAL SCALARS ..

    integer ::   IORD,IPAN,N1,N2,K
    real(kind=DP) ::    C,D

    !-----------------------------------------------------------------------

    ! sort critical points
    do IORD=1,NPAN
      C=CRT(IORD)
      do IPAN=NPAN,IORD,-1
        if(CRT(IPAN) > C)   goto 20
        C=CRT(IPAN)
        CRT(IPAN)=CRT(IORD)
        CRT(IORD)=C
20    continue
      end do
    end do

    !     CALCULATE AN APPROPRIATE MESH

    if (KEYPAN == 0) then
      call MESH0(CRT,NM,NPAN,NPOI,NMIN,NPAND)
    end if

    if (VERBOSITY > 0) then
      write(6,103)
    end if

    N2=0
    do IPAN = 1,NPAN-1

      if (VERBOSITY > 0) then
        write(6,104) IPAN,CRT(IPAN),CRT(IPAN+1),NM(IPAN)
      end if

      N1 = N2 + 1
      N2 = N2 + NM(IPAN)
      if (MESHND >= N2) then
        C = (CRT(IPAN+1)-CRT(IPAN))/DBLE(N2-N1)
        D = CRT(IPAN) - C*DBLE(N1)
        do K = N1,N2
          XRN(K) = C*DBLE(K) + D
          DRN(K) = C
        end do
      else
        goto 70
      end if
    end do

    if (VERBOSITY > 0) then
      write(6,105)
    end if

    MESHN = N2
    !      WRITE(6,101)(K,DRN(K),XRN(K),K=1,MESHN)
    return

70  write (6,102) MESHND
    stop
101 format (/'    NEW MESH  K,DRN,XRN'/(1H ,I3,2F12.7))
102 format ('   *** FROM MESH  :    NXR=',I4,' IS TOO SMALL')
103 format(/50('-')/'I',13X,'SUITABLE RADIAL MESH',15X,'I'/'I',13X, &
    20('*'),15X,'I'/'I',3X,'IPAN',7X,'FROM',7X,'TO',13X,'POINTS  I'/ &
    'I',48X,'I')
104 format('I',2X,I5,2e14.7,I10,'   I')
105 format(50('-'))
  end subroutine MESH

!------------------------------------------------------------------------------
  subroutine MESH0(CRT,NM,NPAN,NAPROX,NMIN,NPAND)
    use shape_constants_mod, only: DP
    implicit none
    ! ***********************************************************
    ! *  THIS SUBROUTINE CALCULATES AN APPROPRIATE MESH FOR
    ! *  THE SHAPE FUNCTIONS. MORE THAN NMIN POINTS BETWEEN TWO
    ! *  CRITICAL POINTS
    ! *  In case of more dense mesh increase NMIN
    ! *
    ! ***********************************************************

    integer, intent(in) :: NPAND

    real(kind=DP) :: CRT(NPAND)
    real(kind=DP) :: DIST,D1
    integer ::   NM(NPAND)
    integer ::   NAPROX,NPAN,NMIN,N,NTOT,I,NA
    intrinsic DFLOAT
    !     DATA NMIN/3/  ! 7

    if ((NPAN-1)*NMIN > NAPROX) then
      write(6,*) NPAN,NMIN,NAPROX
      stop ' INCREASE NUMBER OF POINTS'
    end if

    DIST=ABS(CRT(1)-CRT(NPAN))
    do I=1,NPAN-1
      D1=ABS(CRT(I)-CRT(I+1))
      NM(I)=NAPROX*D1/DIST
      if (NM(I) < NMIN) then
        NM(I)=NMIN
      end if
    end do
    N = NMIN*(NPAN-1)
    N = NAPROX -N
    if (N <= 0)  stop '*** INCREASE NUMBER OF MESH POINTS ***'
    D1=DBLE(N)/DBLE(NAPROX)
    NTOT=N
    do I=1,NPAN-1
      NA = NINT(D1*FLOAT(NM(I)))  ! DBLE instead of FLOAT ?
      if (NM(I) > NMIN .and. (NTOT-NA) > 0) then
        NM(I)=NMIN+NA
        NTOT=NTOT-NA
      end if
    end do
    NM(1) = NM(1) + NTOT
    return
  end subroutine MESH0

end module ShapeStandardMesh_mod
