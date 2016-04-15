      subroutine rhocore(ebot, nsra, ispin, nspin, atom_id, drdi, r, visp, a, b, zat, ircut, rhoc, qc, ecore, ncore, lcore, irmd, ipand)
      implicit none
      external :: corel

      integer, intent(in) :: irmd, ipand
      double precision, intent(in) :: a, b, zat, ebot
      integer, intent(in) :: atom_id, ispin, ncore, nspin, nsra
      double precision, intent(in) :: drdi(irmd), r(irmd), visp(irmd)
      double precision, intent(inout) :: ecore(20)
      double precision, intent(out) :: rhoc(irmd,2)
      integer, intent(in) :: ircut(0:ipand), lcore(20)

!     .. locals ..
      double precision :: qc, qc1, rmax
      integer :: ipr, nr
      
! --------------------------------------------------------------
!     ipr=0 : do not write state dependent information
!     ipr=1 : write something
!     ipr=2 : write all (for debugging)
! --------------------------------------------------------------
        ipr = 0
      
        if (ispin == 1) qc = 0.d0
  !     nr end of core region
        nr = ircut(1)
        rmax = r(nr)

        call corel(nsra,ipr,atom_id,rhoc(1,ispin),visp,ecore,lcore,ncore,drdi,zat,qc1,a,b,ispin,nspin,nr,rmax,irmd,ebot)

        if (ipr /= 0) write (6,fmt="(1x,5('*'),' core-relaxation for ',i3,'th cell was done ',5('*'))") atom_id
        qc = qc + qc1

      endsubroutine rhocore
