! ************************************************************************
SUBROUTINE invslab(gdi,gup,gdow,gin,icheck)
! ************************************************************************
! ************************************************************************
!
! ---> ALGORITM FOR SLAB GEOMETRY
!
! ------------------------------------------------------------------------
!
! ---> factorization D ^-1 = (prod L) * M * (prod U)
!
!      see notes R. Zeller
!
! ------------------------------------------------------------------------
IMPLICIT NONE

!     .. parameters ..
INCLUDE 'inc.p'

! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************

INTEGER, PARAMETER :: LMMAXD = (KREL+KORBIT+1) * (LMAXD+1)**2
INTEGER, parameter :: NDIM= max(NPRINCD*LMMAXD, 1)
INTEGER, PARAMETER :: ALMD= NAEZD*LMMAXD
DOUBLE COMPLEX, PARAMETER :: CI=(0.D0,1.D0)
DOUBLE COMPLEX, PARAMETER :: CZERO=(0.D0,0.D0)
DOUBLE COMPLEX, PARAMETER :: CONE=(1.D0,0.D0)
INTEGER IPVT(NDIM)
COMPLEX*16     GUP(NDIM,NDIM,NLAYERD), &
               GDOW(NDIM,NDIM,NLAYERD), &
               GDI(NDIM,NDIM,NLAYERD), &
               DMAT(NDIM,NDIM,NLAYERD), &
               DINVER(NDIM,NDIM,NLAYERD), &
               F(NDIM,NDIM),E(NDIM,NDIM), &
               G(NDIM,NDIM), &
               CUNIT(NDIM,NDIM),GIN(ALMD,ALMD), &
               GDIOLD(NDIM,NDIM,NLAYERD)
!.. local scalars ..
INTEGER N,LM,INFO, &
        I,J, &
        IROW
!.. data statements ..
INTEGER ICHECK(NLAYERD,NLAYERD)

!.. external subroutines ..
EXTERNAL CINIT,ZCOPY,ZGEMM,ZGETRF,ZGETRS,BTOM

!---> to be changed: set up the triangular matrix.
!     first version:


CALL cinit(ndim*ndim,e)
CALL cinit(ndim*ndim,f)
CALL cinit(ndim*ndim,g)
CALL cinit(ndim*ndim,cunit)
CALL cinit(ndim*ndim*nlayerd,dmat)

DO  n = 1,ndim
  cunit(n,n) = cone
END DO

DO n=1,nlayerd
  CALL zcopy(ndim*ndim,gdi(1,1,n),1,gdiold(1,1,n),1)
END DO


!---> calculate D_1 = (M_11)**(-1)


CALL zcopy(ndim*ndim,gdi(1,1,1),1,e(1,1),1)
CALL zcopy(ndim*ndim,cunit,1,dmat(1,1,1),1)

CALL zgetrf(ndim,ndim,e(1,1),ndim,ipvt,info)
CALL zgetrs('N',ndim,ndim,e(1,1),ndim,ipvt,dmat(1,1,1), ndim,info)

!---> claculate D_N (2 <= N <= NLAYERD)

DO  n = 2,nlayerd
  
!---> F = D(N-1) * M1(N-1)
  
  
  CALL zgemm('N','N',ndim,ndim,ndim,cone,dmat(1,1,n-1),  &
      ndim,gup(1,1,n-1),ndim,czero,f(1,1),ndim)
  
  
!---> D(N) = [MDI(N) - MUP(N-1)*DMAT(N-1)*DOWM(N-1) ]^(-1)
  
  CALL zcopy(ndim*ndim,gdi(1,1,n),1,e(1,1),1)
  CALL zcopy(ndim*ndim,cunit(1,1),1,dmat(1,1,n),1)
  
  CALL zgemm('N','N',ndim,ndim,ndim,-cone,gdow(1,1,n-1),  &
      ndim,f(1,1),ndim,cone,e(1,1),ndim)
  
  CALL zcopy(ndim*ndim,e(1,1),1,dinver(1,1,n),1)
  
  CALL zgetrf(ndim,ndim,e(1,1),ndim,ipvt,info)
  CALL zgetrs('N',ndim,ndim,e(1,1),ndim,ipvt, dmat(1,1,n),ndim,info)
  
END DO

!     At this point the matrix DMAT(ndim,ndim,nlayerd) contains the
!     matrices [of dimension (ndim,ndim)]  D^n, n=1,..,nlayerd

!---> calculate Z_n for 1 =< n <= n-1

DO  n = nlayerd,1,(-1)
  
  IF (n == nlayerd) THEN
    
    CALL zcopy(ndim*ndim,dmat(1,1,nlayerd),1, e(1,1),1)
    
    CALL btom(nlayerd,nlayerd,e,ndim,gin,almd,.false.)
    
    CALL zcopy(ndim*ndim,dmat(1,1,nlayerd),1, gdi(1,1,nlayerd),1)
  ELSE
    
    CALL zgemm('N','N',ndim,ndim,ndim,cone,gdow(1,1,n),  &
        ndim,dmat(1,1,n),ndim,czero,f(1,1),ndim)
    
    CALL bofm(n+1,n+1,e,ndim,gin,almd)
    
    CALL zgemm('N','N',ndim,ndim,ndim,cone,e(1,1),  &
        ndim,f(1,1),ndim,czero,g(1,1),ndim)
    CALL zgemm('N','N',ndim,ndim,ndim,cone,gup(1,1,n),  &
        ndim,g(1,1),ndim,czero,f(1,1),ndim)
    
    DO  lm = 1,ndim
      f(lm,lm) = cone + f(lm,lm)
    END DO
    
    CALL zgemm('N','N',ndim,ndim,ndim,cone,dmat(1,1,n),  &
        ndim,f(1,1),ndim,czero,e(1,1),ndim)
    
    CALL btom(n,n,e,ndim,gin,almd,.false.)
    
    CALL zcopy(ndim*ndim,e(1,1),1,gdi(1,1,n),1)
    
  END IF
  
!     here start the two loops on the row index,
!     in order to span all the matrix
!     and to calculate just the blocks that are needed
!     for the construction of the cluster of green's function
  
  IF (icheck(n,n) == 0) GO TO 200
  
  IF (n == 1) GO TO 100
  
!     this is the loop for element G_ij with i<j
  
  DO  irow=(n-1),1,(-1)
    
    IF (icheck(irow,n) == 1) THEN
      
      CALL bofm(irow+1,n,e,ndim,gin,almd)
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,gup(1,1,irow),  &
          ndim,e(1,1),ndim,czero,f(1,1),ndim)
      
      CALL zgemm('N','N',ndim,ndim,ndim,-cone,dmat(1,1,irow),  &
          ndim,f(1,1),ndim,czero,e(1,1),ndim)
      
      CALL btom(irow,n,e,ndim,gin,almd,.false.)
      
    END IF
    
    
    enddo
  
  100    CONTINUE
  
  IF (n == nlayerd) GO TO 200
  
!     this is the loop for element G_ij with i>j
  
  DO  irow=n+1,nlayerd,1
    
    IF (icheck(irow,n) == 1) THEN
      
      CALL zcopy(ndim*ndim,cunit(1,1),1,e(1,1),1)
      
      CALL bofm(irow,irow,f,ndim,gin,almd)
      
      CALL zgetrf(ndim,ndim,f(1,1),ndim,ipvt,info)
      CALL zgetrs('N',ndim,ndim,f(1,1),ndim,ipvt, e(1,1),ndim,info)
      
      DO i=1,ndim
        DO j=1,ndim
          f(i,j)=gdiold(i,j,irow)-(dinver(i,j,irow)-e(i,j))
        END DO
      END DO
      
      CALL zcopy(ndim*ndim,cunit(1,1),1,e(1,1),1)
      
      CALL zgetrf(ndim,ndim,f(1,1),ndim,ipvt,info)
      CALL zgetrs('N',ndim,ndim,f(1,1),ndim,ipvt, e(1,1),ndim,info)
      
      CALL zgemm('N','N',ndim,ndim,ndim,-cone,e(1,1),  &
          ndim,gdow(1,1,irow-1),ndim,czero,f(1,1),ndim)
      
      CALL bofm(irow-1,n,e,ndim,gin,almd)
      
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,f(1,1),  &
          ndim,e(1,1),ndim,czero,g(1,1),ndim)
      
!     corrected 15.3.2000
      
      CALL btom(irow,n,g,ndim,gin,almd,.false.)
      
      
    END IF
    
    
    enddo
  
  200    CONTINUE
  
END DO





RETURN

END SUBROUTINE invslab
