c ************************************************************************
      SUBROUTINE INVSLAB(GDI,GUP,GDOW,GIN,ICHECK)
c ************************************************************************
c ************************************************************************
c
c ---> ALGORITM FOR SLAB GEOMETRY
c
c ------------------------------------------------------------------------
c
c ---> factorization D ^-1 = (prod L) * M * (prod U)
c
c      see notes R. Zeller
c
c ------------------------------------------------------------------------
      implicit none

c     .. parameters ..
      include 'inc.p'

      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER NDIM
      PARAMETER (NDIM= NPRINCD*LMAXSQ)
      INTEGER ALMD
      PARAMETER (ALMD= NAEZD*LMAXSQ)
      DOUBLE COMPLEX CI,CZERO,CONE
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))
      INTEGER IPVT(NDIM)
      COMPLEX*16     GUP(NDIM,NDIM,NLAYERD),
     +               GDOW(NDIM,NDIM,NLAYERD),
     +               GDI(NDIM,NDIM,NLAYERD),
     +               DMAT(NDIM,NDIM,NLAYERD),
     +               DINVER(NDIM,NDIM,NLAYERD),
     +               F(NDIM,NDIM),E(NDIM,NDIM),
     +               G(NDIM,NDIM),
     +               CUNIT(NDIM,NDIM),GIN(ALMD,ALMD),
     +               GDIOLD(NDIM,NDIM,NLAYERD)
c    .. local scalars ..
      INTEGER N,LM,INFO,
     +        I,J,
     +        IROW
c     .. data statements ..
      INTEGER ICHECK(NLAYERD,NLAYERD)
c
c     .. external subroutines ..
      EXTERNAL CINIT,ZCOPY,ZGEMM,ZGETRF,ZGETRS,BTOM
c
c
c---> to be changed: set up the triangular matrix.
c     first version:
c

      CALL CINIT(NDIM*NDIM,E)
      CALL CINIT(NDIM*NDIM,F)
      CALL CINIT(NDIM*NDIM,G)
      CALL CINIT(NDIM*NDIM,CUNIT)
      CALL CINIT(NDIM*NDIM*NLAYERD,DMAT)

      DO 20 N = 1,NDIM
         CUNIT(N,N) = CONE
 20   CONTINUE

      DO N=1,NLAYERD
         CALL ZCOPY(NDIM*NDIM,GDI(1,1,N),1,GDIOLD(1,1,N),1)
      ENDDO

c
c---> calculate D_1 = (M_11)**(-1)
c

      CALL ZCOPY(NDIM*NDIM,GDI(1,1,1),1,E(1,1),1)
      CALL ZCOPY(NDIM*NDIM,CUNIT,1,DMAT(1,1,1),1)

      CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,DMAT(1,1,1),
     +     NDIM,INFO)
c
c---> claculate D_N (2 <= N <= NLAYERD)
c
      DO 30 N = 2,NLAYERD
c
c---> F = D(N-1) * M1(N-1)
c

         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,DMAT(1,1,N-1),
     +        NDIM,GUP(1,1,N-1),NDIM,CZERO,F(1,1),NDIM)

c
c---> D(N) = [MDI(N) - MUP(N-1)*DMAT(N-1)*DOWM(N-1) ]^(-1)
c
         CALL ZCOPY(NDIM*NDIM,GDI(1,1,N),1,E(1,1),1)
         CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,DMAT(1,1,N),1)

         CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,GDOW(1,1,N-1),
     +        NDIM,F(1,1),NDIM,CONE,E(1,1),NDIM)

         CALL ZCOPY(NDIM*NDIM,E(1,1),1,DINVER(1,1,N),1)

         CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
         CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,
     +        DMAT(1,1,N),NDIM,INFO)

 30   CONTINUE

c     At this point the matrix DMAT(ndim,ndim,nlayerd) contains the
c     matrices [of dimension (ndim,ndim)]  D^n, n=1,..,nlayerd
c
c---> calculate Z_n for 1 =< n <= n-1
c
      DO 40 N = NLAYERD,1,(-1)

        if (n.eq.nlayerd) then

           CALL ZCOPY(NDIM*NDIM,DMAT(1,1,NLAYERD),1,
     +          E(1,1),1)

           CALL BTOM(NLAYERD,NLAYERD,E,NDIM,GIN,ALMD,.FALSE.)

           CALL ZCOPY(NDIM*NDIM,DMAT(1,1,NLAYERD),1,
     +          GDI(1,1,NLAYERD),1)
        else

           CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,GDOW(1,1,N),
     +          NDIM,DMAT(1,1,N),NDIM,CZERO,F(1,1),NDIM)

           CALL BOFM(N+1,N+1,E,NDIM,GIN,ALMD)

           CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,E(1,1),
     +          NDIM,F(1,1),NDIM,CZERO,G(1,1),NDIM)
           CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,GUP(1,1,N),
     +          NDIM,G(1,1),NDIM,CZERO,F(1,1),NDIM)

           DO 45 LM = 1,NDIM
              F(LM,LM) = CONE + F(LM,LM)
 45        CONTINUE

           CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,DMAT(1,1,N),
     +          NDIM,F(1,1),NDIM,CZERO,E(1,1),NDIM)

           CALL BTOM(N,N,E,NDIM,GIN,ALMD,.FALSE.)

           CALL ZCOPY(NDIM*NDIM,E(1,1),1,GDI(1,1,N),1)

        endif

c     here start the two loops on the row index, in order to span all the matrix
c     and to calculates just the blocks that are needed for the construction of
c     the cluster of green's function

        if (icheck(n,n).eq.0) go to 200

        if (n.eq.1) go to 100

c     this is the loop for element G_ij with i<j

        DO 50 IROW=(N-1),1,(-1)

           if (icheck(irow,n).eq.1) then

             CALL BOFM(IROW+1,N,E,NDIM,GIN,ALMD)

             CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,GUP(1,1,IROW),
     +            NDIM,E(1,1),NDIM,CZERO,F(1,1),NDIM)

             CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,DMAT(1,1,IROW),
     +            NDIM,F(1,1),NDIM,CZERO,E(1,1),NDIM)

             CALL BTOM(IROW,N,E,NDIM,GIN,ALMD,.FALSE.)

          endif


 50    ENDDO

 100    CONTINUE

        if (n.eq.nlayerd) go to 200

c     this is the loop for element G_ij with i>j

        DO 60 IROW=N+1,NLAYERD,1

           if (icheck(irow,n).eq.1) then

              CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,E(1,1),1)

              CALL BOFM(IROW,IROW,F,NDIM,GIN,ALMD)

              CALL ZGETRF(NDIM,NDIM,F(1,1),NDIM,IPVT,INFO)
              CALL ZGETRS('N',NDIM,NDIM,F(1,1),NDIM,IPVT,
     +             E(1,1),NDIM,INFO)

              DO I=1,NDIM
                 DO J=1,NDIM
                    F(I,J)=GDIOLD(I,J,IROW)-(DINVER(I,J,IROW)-E(I,J))
                 ENDDO
              ENDDO

              CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,E(1,1),1)

              CALL ZGETRF(NDIM,NDIM,F(1,1),NDIM,IPVT,INFO)
              CALL ZGETRS('N',NDIM,NDIM,F(1,1),NDIM,IPVT,
     +             E(1,1),NDIM,INFO)

              CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,E(1,1),
     +             NDIM,GDOW(1,1,IROW-1),NDIM,CZERO,F(1,1),NDIM)

              CALL BOFM(IROW-1,N,E,NDIM,GIN,ALMD)


              CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,F(1,1),
     +             NDIM,E(1,1),NDIM,CZERO,G(1,1),NDIM)

c     corrected 15.3.2000

              CALL BTOM(IROW,N,G,NDIM,GIN,ALMD,.FALSE.)

c
           endif


 60     ENDDO

 200    CONTINUE

 40   CONTINUE




C
      RETURN
C
      END
