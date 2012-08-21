! bottleneck for lots of k-points
module DLKE0_smat_mod

contains
! Modification of original DLKE0 using sparse matrix
! argument list changed: removed IC, NACLS -> num_cluster_atoms (scalar)
! **********************************************************************
!> @param smat        sparse block matrix, it has to be properly constructed
!>                    and all non-zero blocks have to be marked
!> @param rowstarts   number of first block in each row
subroutine DLKE0_smat(atom_index,smat, ia, ka, kvstr, EIKRM,EIKRP,num_cluster_atoms, &
                 ATOM,NUMN0,INDN0,GINP, &
                 naez, lmmaxd, naclsd)
! **********************************************************************
  implicit none

  double complex, dimension(:), intent(inout) :: smat
  integer, dimension(:), intent(in) :: ia
  integer, dimension(:), intent(in) :: kvstr
  integer, dimension(:), intent(in) :: ka

  integer, intent(in) :: naez
  integer, intent(in) :: lmmaxd
  integer, intent(in) :: naclsd

  integer, intent(in) :: atom_index

  double complex EIKRM(NACLSD),EIKRP(NACLSD)
  integer, intent(in) :: num_cluster_atoms
  integer ATOM(NACLSD),INDN0(NAEZ,naclsd),NUMN0(naez)
  double complex, intent(in) :: GINP(lmmaxd, lmmaxd, NACLSD)

  !     ..
  integer J,LM1,LM2,M,N1,N2,IND1,IND2
  integer :: start
  integer :: gllh_index
  integer :: lmmax1, lmmax2
  integer, parameter :: CZERO = (0.0d0, 0.0d0)
  !     ..
  ! ----------------------------------------------------------------------

  smat = CZERO

  do M = 1,num_cluster_atoms
     
     do N1 = 1,NUMN0(atom_index)
        IND1 = INDN0(atom_index,N1)
        if(ATOM(M).eq.IND1) then

           lmmax1 = kvstr(atom_index + 1) - kvstr(atom_index)
           lmmax2 = kvstr(IND1 + 1) - kvstr(IND1)

           do LM1 = 1,lmmax1
              do LM2 = 1,lmmax2
              
                 gllh_index = ka(ia(atom_index + N1 - 1)) - 1 + (LM2 - 1) * lmmaxd + LM1
                 !GLLH(LM1,AM+LM2,I) = GLLH(LM1,AM+LM2,I) &
                 !                  + EIKRM(M)*GINP(LM2,LM1,M)

                 smat(gllh_index) = smat(gllh_index) + EIKRM(M)*GINP(LM2,LM1,M)
              enddo
           enddo

        endif
     enddo

     J = ATOM(M)

     start = ka(ia(J)) - 1
     do N2 = 1,NUMN0(J)
        IND2 = INDN0(J,N2)
        if(atom_index.eq.IND2) then

           lmmax1 = kvstr(J + 1) - kvstr(J)
           lmmax2 = kvstr(IND2 + 1) - kvstr(IND2)

           do LM1 = 1,lmmax1
              do LM2 = 1,lmmax2
              
                 gllh_index = ka(ia(J + N2 - 1)) - 1 + (LM2 - 1) * lmmaxd + LM1

                 !gllh_index = (J - 1) * NACLSD*lmmaxd*naez + (AN + LM2 - 1) * naez + LM1
                 !GLLH(LM1,AN+LM2,J) = GLLH(LM1,AN+LM2,J) &
                 !                  + EIKRP(M)*GINP(LM1,LM2,M)
                 !GLLH(gllh_index) = GLLH(gllh_index) &
                 !     + EIKRP(M)*GINP(LM1,LM2,M)
                 
                 smat(gllh_index) = smat(gllh_index) + EIKRP(M)*GINP(LM1,LM2,M)
              enddo
           enddo

        endif
     enddo
  enddo

end subroutine DLKE0_smat

end module
