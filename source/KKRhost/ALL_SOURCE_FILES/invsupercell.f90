! ************************************************************************
SUBROUTINE invsupercell(m2,m1,m3,gin,icheck)
! ************************************************************************
!
! ---> ALGORITM FOR SUPERCELL GEOMETRY
!
! ------------------------------------------------------------------------
!
! ---> factorization D ^-1 = (prod L) * M * (prod U)
!
!      see notes R. Zeller
!
! ------------------------------------------------------------------------

implicit none

!.. Parameters ..
include 'inc.p'

! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
!

INTEGER, PARAMETER :: LMMAXD = (KREL+KORBIT+1) * (LMAXD+1)**2
INTEGER, parameter :: NDIM= max(NPRINCD*LMMAXD, 1)
INTEGER, PARAMETER :: ALMD= NAEZD*LMMAXD
DOUBLE COMPLEX, PARAMETER :: CZERO=(0.D0,0.D0)
DOUBLE COMPLEX, PARAMETER :: CONE=(1.D0,0.D0)

!.. array arguments
DOUBLE COMPLEX &
     M1(NDIM,NDIM,NLAYERD), &
     M2(NDIM,NDIM,NLAYERD), &
     M3(NDIM,NDIM,NLAYERD), &
     GIN(ALMD,ALMD)

!.. local scalars
INTEGER N,LM,INFO, &
        IROW,ICOL,NL

!.. local arrays
INTEGER IPVT(NDIM),ICHECK(NLAYERD,NLAYERD)
DOUBLE COMPLEX A(NDIM,NDIM), &
               B(NDIM,NDIM,NLAYERD), &
               C(NDIM,NDIM,NLAYERD), &
               D(NDIM,NDIM,NLAYERD), &
               E(NDIM,NDIM), &
               F(NDIM,NDIM), &
               G(NDIM,NDIM), &
               CUNIT(NDIM,NDIM)
! --------------------------------------------------------------------
!.. External Subroutines ..
EXTERNAL CINIT,ZCOPY,ZGEMM,ZGETRF,ZGETRS,BTOM
INTRINSIC ABS,DIMAG


! ---> START OF THE FACTORIZATION L * M * U


! ---> N =1

!     initialize all the matricex

CALL cinit(ndim*ndim,a)
CALL cinit(ndim*ndim*nlayerd,b)
CALL cinit(ndim*ndim*nlayerd,c)
CALL cinit(ndim*ndim*nlayerd,d)
CALL cinit(ndim*ndim,e)
CALL cinit(ndim*ndim,f)
CALL cinit(ndim*ndim,g)
CALL cinit(ndim*ndim,cunit)

! ------------------------------------------------------------------------

! ---> cunit = complex unity matrix of order NDIM

DO  n = 1,ndim
  cunit(n,n) = cone
END DO


CALL zcopy(ndim*ndim,m2(1,1,1),1,e(1,1),1)
CALL zcopy(ndim*ndim,cunit,1,d(1,1,1),1)
CALL zgetrf(ndim,ndim,e(1,1),ndim,ipvt,info)
CALL zgetrs('N',ndim,ndim,e(1,1),ndim,ipvt,d(1,1,1),ndim,info)


nl=nlayerd


IF (nl == 1) GO TO 980

CALL zcopy(ndim*ndim,m2(1,1,nl),1,a(1,1),1)
CALL zcopy(ndim*ndim,m1(1,1,nl),1,b(1,1,1),1)
CALL zcopy(ndim*ndim,m3(1,1,nl),1,c(1,1,1),1)


! ------------------------------------------------------------------------

! ---> 2 <= N < NL-1

IF (nl == 2) GO TO 970

IF (nl == 3) GO TO 960

DO  n = 2,nl-2
  
  
! ---> E = D(N-1) * C(N-1)
  
  
  CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,n-1),ndim,  &
      c(1,1,n-1),ndim,czero,e(1,1),ndim)
  
  
! ---> F = D(N-1) * M1(N-1)
  
  
  CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,n-1),ndim,  &
      m1(1,1,n-1),ndim,czero,f(1,1),ndim)
  
  
! ---> A = A - B(N-1)*D(N-1)*C(N-1)
  
  
  CALL zgemm('N','N',ndim,ndim,ndim,-cone,b(1,1,n-1),ndim,  &
      e(1,1),ndim,cone,a(1,1),ndim)
  
  
! ---> B(N) = - B(N-1)*D(N-1)*M1(N-1)
  
  CALL zgemm('N','N',ndim,ndim,ndim,-cone,b(1,1,n-1),ndim,  &
      f(1,1),ndim,czero,b(1,1,n),ndim)
  
  
! ---> C(N) = - M3(N-1)*D(N-1)*C(N-1)
  
  
  CALL zgemm('N','N',ndim,ndim,ndim,-cone,m3(1,1,n-1),ndim,  &
      e(1,1),ndim,czero,c(1,1,n),ndim)
  
  
! ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
  
  
  CALL zcopy(ndim*ndim,m2(1,1,n),1,e(1,1),1)
  CALL zgemm('N','N',ndim,ndim,ndim,-cone,m3(1,1,n-1),ndim,  &
      f(1,1),ndim,cone,e(1,1),ndim)
  CALL zcopy(ndim*ndim,cunit(1,1),1,d(1,1,n),1)
  CALL zgetrf(ndim,ndim,e(1,1),ndim,ipvt,info)
  CALL zgetrs('N',ndim,ndim,e(1,1),ndim,ipvt,d(1,1,n),ndim,info)
  
  
END DO

! ------------------------------------------------------------------------

! ---> N = NL - 1



960  n = nl - 1


! ---> E = D(N-1) * C(N-1)


CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,n-1),ndim,  &
    c(1,1,n-1),ndim,czero,e(1,1),ndim)


! ---> F = D(N-1) * M1(N-1)


CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,n-1),ndim,  &
    m1(1,1,n-1),ndim,czero,f(1,1),ndim)


! ---> A = A - B(N-1)*D(N-1)*C(N-1)


CALL zgemm('N','N',ndim,ndim,ndim,-cone,b(1,1,n-1),ndim,  &
    e(1,1),ndim,cone,a(1,1),ndim)


! ---> B(N) = - B(N-1)*D(N-1)*M1(N-1) + M3(N)


CALL zcopy(ndim*ndim,m3(1,1,n),1,b(1,1,n),1)
CALL zgemm('N','N',ndim,ndim,ndim,-cone,b(1,1,n-1),ndim,  &
    f(1,1),ndim,cone,b(1,1,n),ndim)


! ---> C(N) = - M3(N-1)*D(N-1)*C(N-1) + M1(N)


CALL zcopy(ndim*ndim,m1(1,1,n),1,c(1,1,n),1)
CALL zgemm('N','N',ndim,ndim,ndim,-cone,m3(1,1,n-1),ndim,  &
    e(1,1),ndim,cone,c(1,1,n),ndim)


! ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1


CALL zcopy(ndim*ndim,m2(1,1,n),1,e(1,1),1)
CALL zgemm('N','N',ndim,ndim,ndim,-cone,m3(1,1,n-1),ndim,  &
    f(1,1),ndim,cone,e(1,1),ndim)
CALL zcopy(ndim*ndim,cunit(1,1),1,d(1,1,n),1)
CALL zgetrf(ndim,ndim,e(1,1),ndim,ipvt,info)
CALL zgetrs('N',ndim,ndim,e(1,1),ndim,ipvt,d(1,1,n),ndim,info)


! ------------------------------------------------------------------------

! ---> N = NL


970  n = nl


! ---> D(NL) = (A - B(NL-1)*D(NL-1)*C(NL-1))^-1


! ---> E = D(NL-1) * C(NL-1)


CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,nl-1),ndim,  &
    c(1,1,nl-1),ndim,czero,e(1,1),ndim)


! ---> A = A - B(NL-1) * E


CALL zgemm('N','N',ndim,ndim,ndim,-cone,b(1,1,nl-1),ndim,  &
    e(1,1),ndim,cone,a(1,1),ndim)


! ---> D(NL) = (A)^-1


CALL zcopy(ndim*ndim,cunit(1,1),1,d(1,1,nl),1)
CALL zgetrf(ndim,ndim,a(1,1),ndim,ipvt,info)
CALL zgetrs('N',ndim,ndim,a(1,1),ndim,ipvt,d(1,1,nl),ndim,info)


980  CONTINUE                  ! jump label for NL=1

! ------------------------------------------------------------------------

! --->  END OF FACTORIZATION

! ------------------------------------------------------------------------


! ---> HERE IT STARTS LOOP OVER THE DIAGONAL ELEMENT.
! ---> THE PROGRAM CHECKS ID ICHECK(N,N) = 1 AND, IF SO,
! ---> IT CALCULATES THE DIAGONAL BLOCK (N,N)
! ---> THEN IT MAKES TWO LOOPS, ONE OVER A ROW INDEX `IROW`
! ---> AND THE OTHER OVER A COLUMN INDEX `ICOL`, AND CHECKS
! ---> WHICH ARE THE ELEMENTS THAT HAS TO BE CALCULATED.
! ---> (THE ONES FOR WHICH ICHECK = 1)


! ---> IT STARTS THE LOOP OVER N


DO  n=nl,1,(-1)        ! START OF THE LOOP OVER THE DIAGONAL
  
  IF (n == nl) THEN
    
!     write (6,*) 'it calculates the element ','(',nl,',',nl,')'
    
    
! ---> GTOT(NL,NL) = D(NL)
    
    
    CALL zcopy(ndim*ndim,d(1,1,nl),1,e(1,1),1)
    
    CALL btom(nl,nl,e,ndim,gin,almd,.false.)
    
    
  ELSE IF (n == nl-1) THEN
    
    IF (icheck(nl-1,nl-1) == 1) THEN
      
! ---> GTOT(NL-1,NL-1) = D(NL-1) + D(NL-1)*C(NL-1)*D(NL)*B(NL-1)*D(NL-1)
      
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,b(1,1,nl-1),ndim,  &
          d(1,1,nl-1),ndim,czero,e(1,1),ndim)
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,nl),ndim,  &
          e(1,1),ndim,czero,f(1,1),ndim)
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,c(1,1,nl-1),ndim,  &
          f(1,1),ndim,czero,e(1,1),ndim)
      
      DO lm = 1,ndim
        e(lm,lm) = cone + e(lm,lm)
      END DO
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,nl-1),ndim,  &
          e(1,1),ndim,czero,f(1,1),ndim)
      
      CALL btom(nl-1,nl-1,f,ndim,gin,almd,.false.)
      
!     write (6,*) 'it calculates the element ','(',nl-1,',',nl-1,')'
      
    END IF
    
  ELSE
    
    IF (icheck(n,n) == 1) THEN
      
      
! ---> GTOT(N,N) = D(N) + D(N)*( M(N,N+1)*GTOT(N+1,N+1) +
!                  + C(N)*Z(NL,N+1) )*M(N+1,N)*D(N) -
!                  - Z(N,NL)*B(N)*D(N)
      
      CALL bofm(nl,n+1,f,ndim,gin,almd)
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,c(1,1,n),ndim,  &
          f(1,1),ndim,czero,e(1,1),ndim)
      
      CALL bofm(n+1,n+1,f,ndim,gin,almd)
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,m1(1,1,n),ndim,  &
          f(1,1),ndim,cone,e(1,1),ndim)
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,m3(1,1,n),ndim,  &
          d(1,1,n),ndim,czero,f(1,1),ndim)
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,e(1,1),ndim,  &
          f(1,1),ndim,czero,g(1,1),ndim)
      
      DO lm = 1,ndim
        g(lm,lm) = cone + g(lm,lm)
      END DO
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,n),ndim,  &
          g(1,1),ndim,czero,e(1,1),ndim)
      
      CALL zgemm('N','N',ndim,ndim,ndim,cone,b(1,1,n),ndim,  &
          d(1,1,n),ndim,czero,f(1,1),ndim)
      
      CALL bofm(n,nl,g,ndim,gin,almd)
      
      CALL zgemm('N','N',ndim,ndim,ndim,-cone,g(1,1),ndim,  &
          f(1,1),ndim,cone,e(1,1),ndim)
      
      CALL btom(n,n,e,ndim,gin,almd,.false.)
      
!     write (6,*) 'it calculates the element ','(',n,',',n,')'
      
    END IF
    
  END IF
  
  
  DO irow = (n-1),1,(-1) ! LOOP OVER THE ROW FOR THE COLUMN N
    
    IF (icheck(irow,n) == 1) THEN
      
      IF ((n == nl).AND.(irow == (nl-1))) THEN
        
        CALL zgemm('N','N',ndim,ndim,ndim,cone,d(1,1,nl-1),ndim,  &
            c(1,1,nl-1),ndim,czero,f(1,1),ndim)
        
        CALL zgemm('N','N',ndim,ndim,ndim,-cone,f(1,1),ndim,  &
            d(1,1,nl),ndim,czero,g(1,1),ndim)
        
        CALL btom(nl-1,nl,g,ndim,gin,almd,.false.)
        
      ELSE
        
!     M(I,I+1) * Z(I+1,J)
        
        CALL bofm(irow+1,n,e,ndim,gin,almd)
        
        CALL zgemm('N','N',ndim,ndim,ndim,cone,m1(1,1,irow),ndim,  &
            e(1,1),ndim,czero,g(1,1),ndim)
        
!     M(I,I+1) * Z(I+1,J) + C(I) * Z(N,J)
        
        CALL bofm(nl,n,e,ndim,gin,almd)
        
        CALL zgemm('N','N',ndim,ndim,ndim,cone,c(1,1,irow),ndim,  &
            e(1,1),ndim,cone,g(1,1),ndim)
        
!     -D(I) * ( M(I,I+1)*Z(I+1,J)+C(I)*Z(N,J) )
        
        CALL zgemm('N','N',ndim,ndim,ndim,-cone,d(1,1,irow),ndim,  &
            g(1,1),ndim,czero,e(1,1),ndim)
        
        CALL btom(irow,n,e,ndim,gin,almd,.false.)
        
      END IF
      
!     write (6,*) 'it calculates the element ','(',irow,',',n,')'
      
    END IF
    
  END DO                     ! LOOP OVER THE ROW FOR THE COLUMN N
  
  
  DO icol = (n-1),1,(-1)    ! LOOP OVER THE COLUMN FOR THE ROW N
    
    IF (icheck(n,icol) == 1) THEN
      
      IF ((n == nl).AND.(icol == (nl-1))) THEN
        
        CALL zgemm('N','N',ndim,ndim,ndim,cone,b(1,1,nl-1),ndim,  &
            d(1,1,nl-1),ndim,czero,e(1,1),ndim)
        
        CALL zgemm('N','N',ndim,ndim,ndim,-cone,d(1,1,nl),ndim,  &
            e(1,1),ndim,czero,g(1,1),ndim)
        
        CALL btom(nl,nl-1,g,ndim,gin,almd,.false.)
        
      ELSE
        
!     Z(I,J+1) * M(J+1,J)
        
        CALL bofm(n,icol+1,e,ndim,gin,almd)
        
        CALL zgemm('N','N',ndim,ndim,ndim,cone,e(1,1),ndim,  &
            m3(1,1,icol),ndim,czero,g(1,1),ndim)
        
!     Z(I,J+1) * M(J+1,J) + Z(I,N) * B(J)
        
        CALL bofm(n,nl,e,ndim,gin,almd)
        
        CALL zgemm('N','N',ndim,ndim,ndim,cone,e(1,1),ndim,  &
            b(1,1,icol),ndim,cone,g(1,1),ndim)
        
!     -( Z(I,J+1) * M(J+1,J)+Z(I,N) * B(J) ) * D(J)
        
        CALL zgemm('N','N',ndim,ndim,ndim,-cone,g(1,1),ndim,  &
            d(1,1,icol),ndim,czero,e(1,1),ndim)
        
        CALL btom(n,icol,e,ndim,gin,almd,.false.)
        
      END IF
      
!     write (6,*) 'it calculates the element ','(',n,',',icol,')'
      
    END IF
    
  END DO                     ! LOOP OVER THE COLUMN FOR THE ROW N
  
  
END DO



RETURN

END SUBROUTINE invsupercell
