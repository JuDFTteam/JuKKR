SUBROUTINE rotmat(iopt,li,nrot,symopm,vecg)
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
      integer i,idamax,in,j
      double precision costbn,detop,ddet33,dnrm2,omcos, &                 
                       sintbn,sinpb3,tiny,twopi,vfac                    
      character*144 messg 
      parameter(twopi=6.28318530717958648d0) 
      parameter(tiny=1.0d-3) 
! External calls:                                                       
      external daxpy,dcopy,ddet33,rinit,dnrm2,dscal,errmsg,idamax,nrmliz 
! Intrinsic functions:                                                  
      intrinsic  dabs,dacos,dcos,dmax1,dsign,dsin,dsqrt,iabs,idnint 

IF (iopt == -1) THEN
  CALL rinit(symopm,9)
  in = IABS(nrot)
  IF (in == 1) THEN
    CALL dcopy(3,1.d0,0,symopm,4)
  ELSE IF (in == 2.OR.in == 3.OR.in == 4.OR.in == 6) THEN
    sintbn = DSIN(twopi/nrot)
    costbn = DCOS(twopi/nrot)
    omcos  = 1.d0-costbn
    IF (dnrm2(3,vecg,1) < tiny)  &
        CALL errmsg(' ROTMAT: zero rotation vector.$',4)
    CALL nrmliz(1,vecg,vecg)
    symopm(1,1)=omcos*vecg(1)*vecg(1)+costbn
    symopm(1,2)=omcos*vecg(1)*vecg(2)-sintbn*vecg(3)
    symopm(1,3)=omcos*vecg(1)*vecg(3)+sintbn*vecg(2)
    symopm(2,1)=omcos*vecg(2)*vecg(1)+sintbn*vecg(3)
    symopm(2,2)=omcos*vecg(2)*vecg(2)+costbn
    symopm(2,3)=omcos*vecg(2)*vecg(3)-sintbn*vecg(1)
    symopm(3,1)=omcos*vecg(3)*vecg(1)-sintbn*vecg(2)
    symopm(3,2)=omcos*vecg(3)*vecg(2)+sintbn*vecg(1)
    symopm(3,3)=omcos*vecg(3)*vecg(3)+costbn
  ELSE
    CALL errmsg(' ROTMAT: bad nrot.$',3)
  END IF
  IF (li) CALL dscal(9,-1.d0,symopm(1,1),1)
  
ELSE IF (iopt == 1) THEN
! ----- First calculate determinant.
  detop=ddet33(symopm)
  IF (DABS(DABS(detop)-1.0D0) > tiny)  &
      CALL errmsg(' ROTMAT: determinant is not +/- 1$',4)
  detop=DSIGN(1.d0,detop)
  li=detop < 0.d0
! ----- multiply operation symopm with detop
  CALL dscal(9,detop,symopm(1,1),1)
! ----- For the rotation angle we have due to the normalization of v:
! ----- sum_i symopm(i,i) = sum_i (1-cos) v_i*v_i+3*cos = 1 + 2 * cos,
  costbn=-0.5D0
  CALL daxpy(3,0.5D0,symopm(1,1),4,costbn,0)
  IF (DABS(costbn-1.d0) < tiny) THEN
    nrot=1
    CALL rinit(vecg,3)
  ELSE
    nrot=IDNINT(twopi/DACOS(DMAX1(-1.d0,costbn)))
! ------- for nrot > 2 the matrix is non-symmetric and the rotation
! ------- axis can be calculated from the antisymmetric part.
! ------- for nrot = 2 this not possible. However, the squared vector
! ------- components are given by:  mat(i,i) = 2 v_i * v_i - 1.
! ------- This is used for the largest component. The others are taken
! ------- from: mat(i,j) = 2 v_i * v_j for i ne j. This way we also
! ------- get the right phases between the components.
    IF (nrot == 2) THEN
      DO i=1,3
        vecg(i)=0.5D0*(symopm(i,i)+1.0D0)
      END DO
      j=idamax(3,vecg,1)
      IF (vecg(j) < 0.0D0) THEN
        WRITE(messg,402)j,symopm(j,j)
        CALL errmsg(messg,4)
      END IF
      vecg(j)=DSQRT(vecg(j))
      vfac=0.5D0/vecg(j)
      DO i=1,3
        IF (i /= j) vecg(i)=vfac*symopm(i,j)
      END DO
    ELSE
      vecg(1)=symopm(3,2)-symopm(2,3)
      vecg(2)=symopm(1,3)-symopm(3,1)
      vecg(3)=symopm(2,1)-symopm(1,2)
    END IF
! ------- next renormalize at least one component to 1 in order to
! ------- allow for abbreviations as 'D', 'X', 'Y' or 'Z'
    sinpb3=DSQRT(.75D0)
    IF (DABS((sinpb3-DABS(vecg(1)))*(sinpb3-DABS(vecg(2)))*  &
          (sinpb3-DABS(vecg(3)))) > tiny) THEN
      DO j=3,1,-1
        vfac=DABS(vecg(j))
        IF(vfac > tiny) CALL dscal(3,1.d0/vfac,vecg,1)
      END DO
    END IF
  END IF
  CALL dscal(9,detop,symopm(1,1),1)
END IF

402 FORMAT(' ROTMAT: Bad component ',i1,' of operation ',  &
    '. Diagonal element =',f9.5,'$')
END SUBROUTINE rotmat
