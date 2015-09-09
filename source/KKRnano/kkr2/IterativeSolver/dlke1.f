C     >> Input parameters
C>    @param     ALAT    lattice constant a
C>    @param     NACLS   number of atoms in cluster
C>    @param     RR      array of real space vectors
C>    @param     EZOA
C>    @param     BZKP
C>    @param     nrd     There are nrd+1 real space vectors in RR
C>    @param     naclsd. maximal number of atoms in cluster

C     << Output parameters
C>    @param     EIKRM   Fourier exponential factor with minus sign
C>    @param     EIKRP   Fourier exponential factor with plus sign

      subroutine dlke1(alat,nacls,rr,ezoa,
     &                 bzkp,eikrm,eikrp,
     &                 nrd, naclsd)
      implicit none
c ----------------------------------------------------------------------
c
c     Fourier transformation of the cluster Greens function
c     Prepares the calculation (calculates Fourier factors) for dlke0
c ----------------------------------------------------------------------

      integer nrd
      integer naclsd

      double complex ci,cone
      parameter (ci= (0.0d0,1.0d0),cone=(1.d0,0.d0))
c     ..
c     .. scalar arguments ..
      double precision alat
c     ..
c     .. array arguments ..
      integer ezoa(*),nacls
      double complex eikrp(naclsd),eikrm(naclsd)
      double precision bzkp(*),rr(3,0:nrd)
c     ..
c     .. local scalars ..
      double precision convpu,tpi
      integer m
      double complex tt
c     ..
c     .. local arrays ..
      double complex arg(3)
c     ..
c     .. external subroutines ..
      external cinit,test,opt,zaxpy
c     ..
c     .. intrinsic functions ..
      intrinsic atan,exp
c     ..
c     .. save statement ..
      save
c     ..
c
      tpi = 8.0d0*atan(1.0d0)         
      convpu = alat/tpi

      do 90 m = 1,nacls

c
c     
c     Here we do   --                  nn'
c                  \                   ii'          ii'
c                  /  exp(+ik(x  -x ))G   (E)  =   G   (k,E)
c                  --          n'  n   LL'          LL'
c                  n'
c  Be careful a minus sign must be included here. RR is not
c  symmetric around each atom. The minus comes from the fact that
c  the repulsive potential GF is calculated for 0n and not n0!                   
c  and that is why we need a minus sign extra!
c  

           arg(1) = -ci*tpi*rr(1,ezoa(m))
           arg(2) = -ci*tpi*rr(2,ezoa(m))
           arg(3) = -ci*tpi*rr(3,ezoa(m))
c
        tt = bzkp(1)*arg(1)+bzkp(2)*arg(2)+bzkp(3)*arg(3)
c
c  convert to p.u. and multiply with 1/2
        eikrp(m) = exp( tt) * convpu * 0.5d0
        eikrm(m) = exp(-tt) * convpu * 0.5d0
c
 90   continue                    

      return
      end
