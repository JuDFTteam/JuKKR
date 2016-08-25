c ************************************************************************
      SUBROUTINE INVSLAB_SO(GDI,GUP,GDOW,GIN,ICHECK)
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
      INTEGER LMMAXSO
      PARAMETER (LMMAXSO= NSPOD*(LMAXD+1)**2)
      INTEGER NDIM
      PARAMETER (NDIM= NPRINCD*LMAXSQ)
      INTEGER NDIMSO
      PARAMETER (NDIMSO=NSPOD*NPRINCD*LMAXSQ)
      INTEGER ALMD
      PARAMETER (ALMD= NAEZD*LMAXSQ)
      INTEGER ALMDSO
      PARAMETER (ALMDSO= NSPOD*NAEZD*LMAXSQ)
      DOUBLE COMPLEX CI,CZERO,CONE
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))
      INTEGER IPVT(NDIMSO)
      COMPLEX*16     GUP(NDIMSO,NDIMSO,NLAYERD),
     +               GDOW(NDIMSO,NDIMSO,NLAYERD),
     +               GDI(NDIMSO,NDIMSO,NLAYERD),
     +               DMAT(NDIMSO,NDIMSO,NLAYERD),
     +               DINVER(NDIMSO,NDIMSO,NLAYERD),
     +               F(NDIMSO,NDIMSO),E(NDIMSO,NDIMSO),
     +               G(NDIMSO,NDIMSO),
     +               CUNIT(NDIMSO,NDIMSO),GIN(ALMDSO,ALMDSO),
     +               GDIOLD(NDIMSO,NDIMSO,NLAYERD)
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

      CALL CINIT(NDIMSO*NDIMSO,E)
      CALL CINIT(NDIMSO*NDIMSO,F)
      CALL CINIT(NDIMSO*NDIMSO,G)
      CALL CINIT(NDIMSO*NDIMSO,CUNIT)
      CALL CINIT(NDIMSO*NDIMSO*NLAYERD,DMAT)

      DO 20 N = 1,NDIMSO
         CUNIT(N,N) = CONE
 20   CONTINUE

      DO N=1,NLAYERD
         CALL ZCOPY(NDIMSO*NDIMSO,GDI(1,1,N),1,GDIOLD(1,1,N),1)
      ENDDO

c
c---> calculate D_1 = (M_11)**(-1)
c

      CALL ZCOPY(NDIMSO*NDIMSO,GDI(1,1,1),1,E(1,1),1)
      CALL ZCOPY(NDIMSO*NDIMSO,CUNIT,1,DMAT(1,1,1),1)

      CALL ZGETRF(NDIMSO,NDIMSO,E(1,1),NDIMSO,IPVT,INFO)
      CALL ZGETRS('N',NDIMSO,NDIMSO,E(1,1),NDIMSO,IPVT,DMAT(1,1,1),
     +     NDIMSO,INFO)
c
c---> claculate D_N (2 <= N <= NLAYERD)
c
      DO 30 N = 2,NLAYERD
c
c---> F = D(N-1) * M1(N-1)
c

         CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,CONE,DMAT(1,1,N-1),
     +        NDIMSO,GUP(1,1,N-1),NDIMSO,CZERO,F(1,1),NDIMSO)

c
c---> D(N) = [MDI(N) - MUP(N-1)*DMAT(N-1)*DOWM(N-1) ]^(-1)
c
         CALL ZCOPY(NDIMSO*NDIMSO,GDI(1,1,N),1,E(1,1),1)
         CALL ZCOPY(NDIMSO*NDIMSO,CUNIT(1,1),1,DMAT(1,1,N),1)

         CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,-CONE,GDOW(1,1,N-1),
     +        NDIMSO,F(1,1),NDIMSO,CONE,E(1,1),NDIMSO)

         CALL ZCOPY(NDIMSO*NDIMSO,E(1,1),1,DINVER(1,1,N),1)

         CALL ZGETRF(NDIMSO,NDIMSO,E(1,1),NDIMSO,IPVT,INFO)
         CALL ZGETRS('N',NDIMSO,NDIMSO,E(1,1),NDIMSO,IPVT,
     +        DMAT(1,1,N),NDIMSO,INFO)

 30   CONTINUE

c     At this point the matrix DMAT(ndim,ndim,nlayerd) contains the
c     matrices [of dimension (ndim,ndim)]  D^n, n=1,..,nlayerd
c
c---> calculate Z_n for 1 =< n <= n-1
c
      DO 40 N = NLAYERD,1,(-1)

        if (n.eq.nlayerd) then

           CALL ZCOPY(NDIMSO*NDIMSO,DMAT(1,1,NLAYERD),1,
     +          E(1,1),1)

           CALL BTOM(NLAYERD,NLAYERD,E,NDIMSO,GIN,ALMDSO,.FALSE.)

           CALL ZCOPY(NDIMSO*NDIMSO,DMAT(1,1,NLAYERD),1,
     +          GDI(1,1,NLAYERD),1)
        else

           CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,CONE,GDOW(1,1,N),
     +          NDIMSO,DMAT(1,1,N),NDIMSO,CZERO,F(1,1),NDIMSO)

           CALL BOFM(N+1,N+1,E,NDIMSO,GIN,ALMDSO)

           CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,CONE,E(1,1),
     +          NDIMSO,F(1,1),NDIMSO,CZERO,G(1,1),NDIMSO)
           CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,CONE,GUP(1,1,N),
     +          NDIMSO,G(1,1),NDIMSO,CZERO,F(1,1),NDIMSO)

           DO 45 LM = 1,NDIMSO
              F(LM,LM) = CONE + F(LM,LM)
 45        CONTINUE

           CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,CONE,DMAT(1,1,N),
     +          NDIMSO,F(1,1),NDIMSO,CZERO,E(1,1),NDIMSO)

           CALL BTOM(N,N,E,NDIMSO,GIN,ALMDSO,.FALSE.)

           CALL ZCOPY(NDIMSO*NDIMSO,E(1,1),1,GDI(1,1,N),1)

        endif

c     here start the two loops on the row index, in order to span all the matrix
c     and to calculates just the blocks that are needed for the construction of
c     the cluster of green's function

        if (icheck(n,n).eq.0) go to 200

        if (n.eq.1) go to 100

c     this is the loop for element G_ij with i<j

        DO 50 IROW=(N-1),1,(-1)

           if (icheck(irow,n).eq.1) then

             CALL BOFM(IROW+1,N,E,NDIMSO,GIN,ALMDSO)

             CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,CONE,GUP(1,1,IROW),
     +            NDIMSO,E(1,1),NDIMSO,CZERO,F(1,1),NDIMSO)

             CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,-CONE,
     +            DMAT(1,1,IROW),NDIMSO,F(1,1),NDIMSO,CZERO,E(1,1),
     +            NDIMSO)

             CALL BTOM(IROW,N,E,NDIMSO,GIN,ALMDSO,.FALSE.)

          endif


 50    ENDDO

 100    CONTINUE

        if (n.eq.nlayerd) go to 200

c     this is the loop for element G_ij with i>j

        DO 60 IROW=N+1,NLAYERD,1

           if (icheck(irow,n).eq.1) then

              CALL ZCOPY(NDIMSO*NDIMSO,CUNIT(1,1),1,E(1,1),1)

              CALL BOFM(IROW,IROW,F,NDIMSO,GIN,ALMDSO)

              CALL ZGETRF(NDIMSO,NDIMSO,F(1,1),NDIMSO,IPVT,INFO)
              CALL ZGETRS('N',NDIMSO,NDIMSO,F(1,1),NDIMSO,IPVT,
     +             E(1,1),NDIMSO,INFO)

              DO I=1,NDIMSO
                 DO J=1,NDIMSO
                    F(I,J)=GDIOLD(I,J,IROW)-(DINVER(I,J,IROW)-E(I,J))
                 ENDDO
              ENDDO

              CALL ZCOPY(NDIMSO*NDIMSO,CUNIT(1,1),1,E(1,1),1)

              CALL ZGETRF(NDIMSO,NDIMSO,F(1,1),NDIMSO,IPVT,INFO)
              CALL ZGETRS('N',NDIMSO,NDIMSO,F(1,1),NDIMSO,IPVT,
     +             E(1,1),NDIMSO,INFO)

              CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,-CONE,E(1,1),
     +             NDIMSO,GDOW(1,1,IROW-1),NDIMSO,CZERO,F(1,1),NDIMSO)

              CALL BOFM(IROW-1,N,E,NDIMSO,GIN,ALMDSO)


              CALL ZGEMM('N','N',NDIMSO,NDIMSO,NDIMSO,CONE,F(1,1),
     +             NDIMSO,E(1,1),NDIMSO,CZERO,G(1,1),NDIMSO)

c     corrected 15.3.2000

              CALL BTOM(IROW,N,G,NDIMSO,GIN,ALMDSO,.FALSE.)

c
           endif


 60     ENDDO

 200    CONTINUE

 40   CONTINUE




C
      RETURN
C
      END
