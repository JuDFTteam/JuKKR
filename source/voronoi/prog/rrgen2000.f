C 02.08.95 ***************************************************************
      SUBROUTINE RRGEN (BV1,AA,LSURF,RR,NR)
C ************************************************************************
c     generates a number of real space vectors to construct the
c     clusters representing the local surrounding of the atoms in
c     routine CLSGEN
c ------------------------------------------------------------------------
      implicit none
      include 'inc.geometry'
c
      LOGICAL LSURF
c     .. scalar arguments
      INTEGER LATT,                 ! lattice type
     +        NR                    ! number of real space vectors
      REAL*8          AA           ! lattice constant alata
c
C    .. array arguments
      REAL*8         
     +     BV1(3,3),                ! true basis vectors
     +     RR(3,0:NRD)              ! real space vectors
C
c    .. local scalars     
c
      REAL*8         
     +     R,RS,
     +     R1,R2,R3,RR2,epsshl,rmax
C
      INTEGER
     +     POS,
     +     N1,N2,N3,
     +     i,j,k
c
c     .. local arrays
      REAL*8         
     +     RABS(NRD),RR1(3,NRD),
     +     V(3),
     +     VX(3),VY(3),VZ(3),
     +     VX0(3),VY0(3),VZ0(3)
      INTEGER IND(NRD)
c
c    
      DATA  epsshl /1.0d-5/
c
      INTRINSIC ABS,DFLOAT,MIN,SQRT
      LOGICAL TEST,OPT
      EXTERNAL DSORT,VADD,VEQ,VSUB,SCALPR,TEST,OPT
c
c ------------------------------------------------------------------------
      write(6,*) '>>> RRGEN: generation of real space mesh rr(nr)'
      CALL SCALPR(BV1(1,1),BV1(1,1),R1)
      CALL SCALPR(BV1(1,2),BV1(1,2),R2)
      CALL SCALPR(BV1(1,3),BV1(1,3),R3)
      RMAX = 8.d0
c
      R1=SQRT(R1) !/AA
      R2=SQRT(R2) !/AA
      R3=SQRT(R3) !/AA
      R=1.5D0*RMAX+SQRT(R1*R1+R2*R2+R3*R3) + EPSSHL
      RS=R*R
c
                        N1 = (R/R1) ! +1  
                        N2 = (R/R2) ! +1
      IF (.NOT.(LSURF)) N3 = (R/R3) ! +1
                        N1 = MIN(12,N1)
                        N2 = MIN(12,N2)
      IF (.NOT.(LSURF)) N3 = MIN(12,N3)
                        N1 = MAX(2,N1)
                        N2 = MAX(2,N2)
      IF (.NOT.(LSURF)) N3 = MAX(2,N3)

      IF (LSURF) N3=0
      IF (OPT('SLAB    ').AND. .NOT.TEST('CONT    ')) N3 = 0
      IF (OPT('WIRE    ').AND. .NOT.TEST('CONT    ')) THEN 
        N1 = 0
        N2 = 0
      END IF
c
      WRITE(6,1100) R
 1100 FORMAT(' r        :',F10.2)
      write(6,1101) rs
 1101 FORMAT(' r**2     :',F10.2)
      write(6,1102) n1,n2,n3
 1102 FORMAT(' n1,n2,n3 :',3I5)
      nr=0
      rr(1,0)=0.0d0
      rr(2,0)=0.0d0
      rr(3,0)=0.0d0
      if (test('RR      ')) 
     +     write (6,1003) nr,rr(1,nr),rr(2,nr),rr(3,nr),0.0
      call vmul(bv1(1,1),dfloat(-n1-1),vx0(1))
      call vmul(bv1(1,2),dfloat(-n2-1),vy0(1))
      call vmul(bv1(1,3),dfloat(-n3-1),vz0(1))
      call veq(vx0,vx)
      do 10 i=-n1,n1
        call vadd(vx,bv1(1,1),vx)
        call veq(vy0,vy)
        do 20 j= -n2,n2
          call vadd(vy,bv1(1,2),vy)
          call veq(vz0,vz)
          do 30 k= -n3,n3
            call vadd(vz,bv1(1,3),vz)
            call vadd(vx,vy,v)
            call vadd(v,vz,v)
            call scalpr(v,v,rr2)
c           if ( ((rr2.le.rs*aa*aa).or.(abs(i)+abs(j)+abs(k).le.6))
            if ( ((rr2.le.rs).or.(abs(i)+abs(j)+abs(k).le.6))
     +            .and.(rr2.gt.epsshl) ) then
              nr=nr+1
              rr1(1,nr) = v(1) !/aa
              rr1(2,nr) = v(2) !/aa
              rr1(3,nr) = v(3) !/aa
              rabs(nr)= sqrt(rr2) !/AA
            end if  
 30       continue
 20     continue  
 10   continue
c ------------------------------------------------------------------------
      if (nr.gt.nrd) then
        write(6,*) 'RRGEN: Please, change the parameter nrd in',
     *   ' inc.p to ',nr
        stop
      end if  
c ------------------------------------------------------------------------
      CALL DSORT(RABS,IND,NR,POS)
      DO 105 I = 1,NR
        POS = IND(I)
        RR(1,I) = RR1(1,POS)
        RR(2,I) = RR1(2,POS)
        RR(3,I) = RR1(3,POS)
        IF (TEST('RR      ')) WRITE (6,1003) 
     +       I,RR(1,I),RR(2,I),RR(3,I),RABS(POS)
 105  END DO
      
      write(6,1001) nr

      return
 1000 FORMAT(I6,3F12.5)
 1001 FORMAT(I6 ,' real space vectors created.')
 1003 FORMAT(I6,3F12.3,F15.4)
      END                           ! RRGEN
