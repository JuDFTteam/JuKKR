      SUBROUTINE CORLSD(RS,ZTA,EC,VCUP,VCDN,ECRS,ECZTA,ALFC)
c.....-----------------------------------------------------------------
c     uniform-gas correlation of perdew and wang 1991
c.....-----------------------------------------------------------------
c     input: seitz radius (rs), relative spin polarization (zta)
c     output: correlation energy per electron (ec),
c             up- and down-spin potentials (vcup,vcdn),
c             derivatives of ec wrt rs (ecrs) &zta (eczta).
c     output: correlation contribution (alfc) to the spin stiffness
c.....-----------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALFC,EC,ECRS,ECZTA,RS,VCDN,VCUP,ZTA
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALFM,ALFRSM,COMM,EP,EPRS,EU,EURS,F,FZ,FZZ,GAM,
     +                 THRD,THRD4,Z4
C     ..
C     .. External Subroutines ..
      EXTERNAL GCOR91
C     ..
C     .. Save statement ..
      SAVE GAM,FZZ,THRD,THRD4
C     ..
C     .. Data statements ..
c.....-----------------------------------------------------------------
      DATA GAM,FZZ/0.5198421d0,1.709921d0/
      DATA THRD,THRD4/0.333333333333d0,1.333333333333d0/
C     ..
c.....-----------------------------------------------------------------
      F = ((1.d0+ZTA)**THRD4+ (1.d0-ZTA)**THRD4-2.d0)/GAM
      CALL GCOR91(0.0310907d0,0.21370d0,7.5957d0,3.5876d0,1.6382d0,
     +            0.49294d0,1.00d0,RS,EU,EURS)
      CALL GCOR91(0.01554535d0,0.20548d0,14.1189d0,6.1977d0,3.3662d0,
     +            0.62517d0,1.00d0,RS,EP,EPRS)
      CALL GCOR91(0.0168869d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,
     +            0.49671d0,1.00d0,RS,ALFM,ALFRSM)
c  alfm is minus the spin stiffness alfc
      ALFC = -ALFM
      Z4 = ZTA**4
      EC = EU* (1.d0-F*Z4) + EP*F*Z4 - ALFM*F* (1.d0-Z4)/FZZ
c  energy done. now the potential:
      ECRS = EURS* (1.d0-F*Z4) + EPRS*F*Z4 - ALFRSM*F* (1.d0-Z4)/FZZ
      FZ = THRD4* ((1.d0+ZTA)**THRD- (1.d0-ZTA)**THRD)/GAM
      ECZTA = 4.d0* (ZTA**3)*F* (EP-EU+ALFM/FZZ) +
     +        FZ* (Z4*EP-Z4*EU- (1.d0-Z4)*ALFM/FZZ)
      COMM = EC - RS*ECRS/3.d0 - ZTA*ECZTA
      VCUP = COMM + ECZTA
      VCDN = COMM - ECZTA

      RETURN
      END
