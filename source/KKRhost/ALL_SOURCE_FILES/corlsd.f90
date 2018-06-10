SUBROUTINE corlsd(rs,zta,ec,vcup,vcdn,ecrs,eczta,alfc)
!.....-----------------------------------------------------------------
!     uniform-gas correlation of perdew and wang 1991
!.....-----------------------------------------------------------------
!     input: seitz radius (rs), relative spin polarization (zta)
!     output: correlation energy per electron (ec),
!             up- and down-spin potentials (vcup,vcdn),
!             derivatives of ec wrt rs (ecrs) &zta (eczta).
!     output: correlation contribution (alfc) to the spin stiffness
!.....-----------------------------------------------------------------
!.. Scalar Arguments ..
      DOUBLE PRECISION ALFC,EC,ECRS,ECZTA,RS,VCDN,VCUP,ZTA
!..
!.. Local Scalars ..
      DOUBLE PRECISION ALFM,ALFRSM,COMM,EP,EPRS,EU,EURS,F,FZ,FZZ,GAM, &
                       THRD,THRD4,Z4
!..
!.. External Subroutines ..
      EXTERNAL GCOR91
!..
!.. Save statement ..
      SAVE GAM,FZZ,THRD,THRD4
!..
!.. Data statements ..
!.....-----------------------------------------------------------------
      DATA GAM,FZZ/0.5198421d0,1.709921d0/
      DATA THRD,THRD4/0.333333333333d0,1.333333333333d0/
!..
!.....-----------------------------------------------------------------
f = ((1.d0+zta)**thrd4+ (1.d0-zta)**thrd4-2.d0)/gam
CALL gcor91(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,  &
    0.49294D0,1.00D0,rs,eu,eurs)
CALL gcor91(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,  &
    0.62517D0,1.00D0,rs,ep,eprs)
CALL gcor91(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,  &
    0.49671D0,1.00D0,rs,alfm,alfrsm)
!  alfm is minus the spin stiffness alfc
alfc = -alfm
z4 = zta**4
ec = eu* (1.d0-f*z4) + ep*f*z4 - alfm*f* (1.d0-z4)/fzz
!  energy done. now the potential:
ecrs = eurs* (1.d0-f*z4) + eprs*f*z4 - alfrsm*f* (1.d0-z4)/fzz
fz = thrd4* ((1.d0+zta)**thrd- (1.d0-zta)**thrd)/gam
eczta = 4.d0* (zta**3)*f* (ep-eu+alfm/fzz) +  &
    fz* (z4*ep-z4*eu- (1.d0-z4)*alfm/fzz)
comm = ec - rs*ecrs/3.d0 - zta*eczta
vcup = comm + eczta
vcdn = comm - eczta

RETURN
END SUBROUTINE corlsd
