! 13.10.95 ***************************************************************
SUBROUTINE ecoub(cmom,ecou,lmax,nspin,natyp,rho2ns,vm2z,z,r,drdi,  &
    irws,kvmad,kshape,ircut,ipan,imaxsh,ifunm, ilm_map,ntcell,gsh,thetas,lmsp)
! ************************************************************************

!     attention : energy zero ---> electro static zero

!     calculate the electrostatic potential-energies without the
!     electron-nuclear interaction in the cell itself .
!     the energy of the representive atom i is given by

!                          rc
!      ecou(i) =  1/2 (  {  s dr' vm2z(r',lm,i)*rho2ns(r',lm,i,1) }
!                           0

!                                       -  z(i) * vmad ( ri )     )


!                                         ( {..} = summed over lm )
!             (see notes by b.drittler)
!     vm2z is the coulomb potential of the atom without the nuclear
!             potential of the atom
!     rho2ns(...,1) is the real charge density times r**2

!      both developed into spherical harmonics . (see deck rholm)

!     z    is the nuclear charge of the atom

!     vmad ( ri ) is a generalized madelung potential
!                 = 1/sqrt(4 pi) * vm2z(irws,1,is)
!                         - sqrt(4 pi) * 2 * cmom(1,ipot) / rws

!                                        ( <..> = spherical averaged )

!     attention : this subroutine has to be called before the
!                 exchange correlation potential is added to
!                 the potential vm2z .
!                 the energy calculated here is splitted into
!                 l-dependent parts to see the l -convergency .

!     attention : in case of shape corrections the contribution of
!                 the coulomb potential the of the nucleus is
!                 analytically cancelled only in the muffin tin sphere
!                 in the interstial region it has to be taken into
!                 account ! see deck epotins

!                 modified for band structure code
!                               b.drittler   jan. 1990
!-----------------------------------------------------------------------
      use mod_types, only: t_inc
      IMPLICIT NONE
!.. Parameters ..
      include 'inc.p'
!..
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
!..
!.. Scalar Arguments ..
      INTEGER KSHAPE,KVMAD,LMAX,NATYP,NSPIN
!..
!.. Array Arguments ..
DOUBLE PRECISION CMOM(LMPOTD,*),DRDI(IRMD,*),ECOU(0:LPOTD,*), &
                 GSH(*),R(IRMD,*),RHO2NS(IRMD,LMPOTD,NATYPD,*), &
                 THETAS(IRID,NFUND,*),VM2Z(IRMD,LMPOTD,*),Z(*)
INTEGER IFUNM(NATYPD,*),ILM_MAP(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*), &
        IRCUT(0:IPAND,*),IRWS(*),NTCELL(*),LMSP(NATYPD,*)
!..
!.. Local Scalars ..
DOUBLE PRECISION RFPI,RHOSP,SIGN,VM,VMAD
INTEGER I,IATYP,ICELL,IFUN,IPAN1,IPOT,IR,IRC1,IRH,IRS1,ISPIN,J,L, &
        LM,LM2,M
!..
!.. Local Arrays ..
      DOUBLE PRECISION ER(IRMD)
!..
!.. External Subroutines ..
      EXTERNAL SIMP3,SIMPK
!..
!.. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
!..
rfpi = SQRT(16.d0*ATAN(1.0D0))

DO  iatyp = 1,natyp
  
  IF (kshape /= 0) THEN
    ipan1 = ipan(iatyp)
    icell = ntcell(iatyp)
    irs1 = ircut(1,iatyp)
    irc1 = ircut(ipan1,iatyp)
  ELSE
    irs1 = irws(iatyp)
    irc1 = irs1
  END IF
  
!--->   determine the right potential numbers - the coulomb potential
!       is not spin dependend
  
  ipot = iatyp*nspin
  
  DO  l = 0,lmax
    
    DO  i = 1,irc1
      er(i) = 0.0D0
    END DO
    
    DO  ispin = 1,nspin
      
      IF (ispin == nspin) THEN
        SIGN = 1.0D0
      ELSE
        SIGN = -1.0D0
      END IF
      
      DO  m = -l,l
        lm = l*l + l + m + 1
        
        DO  i = 1,irs1
          rhosp = (rho2ns(i,lm,iatyp,1)+ SIGN*rho2ns(i,lm,iatyp,nspin))/4.0D0
          er(i) = er(i) + rhosp*vm2z(i,lm,ipot)
        END DO
        
        IF (kshape /= 0) THEN
          
!--->           convolute with shape function
          
          DO  j = imaxsh(lm-1) + 1,imaxsh(lm)
            lm2 = ilm_map(j,2)
            IF (lmsp(icell,ilm_map(j,3)) > 0) THEN
              ifun = ifunm(icell,ilm_map(j,3))
              
              IF (lm2 == 1) THEN
                DO  ir = irs1 + 1,irc1
                  irh = ir - irs1
                  rhosp = (rho2ns(ir,lm,iatyp,1)+  &
                      SIGN*rho2ns(ir,lm,iatyp,nspin))/2.0D0
                  
!--->                 remember that in the interstial -2z/r has
!                     to be taken into account
                  
                  er(ir) = er(ir) + rhosp*gsh(j)* thetas(irh,ifun,icell)*  &
                      (vm2z(ir,1,ipot)/2.0D0- z(iatyp)/r(ir,iatyp)*rfpi)
                END DO
                
              ELSE              ! (LM2.EQ.1)
                
                DO  ir = irs1 + 1,irc1
                  irh = ir - irs1
                  rhosp = (rho2ns(ir,lm,iatyp,1)+  &
                      SIGN*rho2ns(ir,lm,iatyp,nspin))/2.0D0
                  er(ir) = er(ir) + rhosp*gsh(j)*  &
                      thetas(irh,ifun,icell)*vm2z(ir,lm2,ipot)/ 2.0D0
                END DO
                
              END IF            ! (LM2.EQ.1)
            END IF
            
          END DO
          
        END IF                ! (KSHAPE.NE.0)
        
      END DO
      
    END DO
    
!--->     now integrate
    
    IF (kshape == 0) THEN
      CALL simp3(er,ecou(l,iatyp),1,irs1,drdi(1,iatyp))
    ELSE
      CALL simpk(er,ecou(l,iatyp),ipan1,ircut(0,iatyp), drdi(1,iatyp))
    END IF
    
  END DO
  
  
!--->   calculate the madelung potential
  
  vmad = vm2z(irs1,1,ipot)/rfpi - rfpi*2.0D0*cmom(1,iatyp)/r(irs1,iatyp)
  
!--->   add to ecou
  
  ecou(0,iatyp) = ecou(0,iatyp) - z(iatyp)*vmad/2.0D0
  
!--->   option to calculate full generalized madelung potential
!                                    rc
!       vm(rn) = vmad +2*sqrt(4*pi)* s  dr*r*rho(lm=1,r)
!                                    0
  IF (kvmad == 1) THEN
    er(1) = 0.0D0
    
    DO  i = 2,irs1
      er(i) = rho2ns(i,1,iatyp,1)/r(i,iatyp)
    END DO
    
    CALL simp3(er,vm,1,irs1,drdi(1,iatyp))
    vm = 2.0D0*rfpi*vm + vmad
    
!         atom nr. iatyp is the iatyp-th atom on the potential cards
!         e. g., in binary alloys iatyp=1 and iatyp=2 refer to host
    
    IF(t_inc%i_write>0) WRITE (1337,FMT=9010) iatyp,vmad
    IF(t_inc%i_write>0) WRITE (1337,FMT=9000) iatyp,vm
  END IF                      ! (KVMAD.EQ.1)
  
END DO

RETURN

9000 FORMAT (10X,'full generalized madelung pot. for atom',1X,i3,1X,  &
    ': ',1P,d14.6)
9010 FORMAT (10X,'     generalized madelung pot. for atom',1X,i3,1X,  &
    ': ',1P,d14.6)
END SUBROUTINE ecoub
