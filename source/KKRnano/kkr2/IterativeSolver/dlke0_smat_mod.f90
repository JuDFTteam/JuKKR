module DLKE0_smat_mod

contains

  ! **********************************************************************
  subroutine DLKE0_smat(site_index,smat,ia,ka,kvstr,EIKRM,EIKRP,NACLS, &
  ATOM,NUMN0,INDN0,GINP, &
  naez, lmmaxd, naclsd)
  ! **********************************************************************
    implicit none

    double complex, dimension(:), intent(inout) :: smat
    integer, dimension(:), intent(in) :: ia
    integer, dimension(:), intent(in) :: ka
    integer, dimension(:), intent(in) :: kvstr

    integer naez
    integer lmmaxd
    integer naclsd

    integer :: site_index
    double complex, intent(in) :: GINP(lmmaxd, lmmaxd, NACLSD)
    !double complex, dimension(:) :: GLLH(lmmaxd, NACLSD * lmmaxd, *)
    double complex EIKRM(NACLSD),EIKRP(NACLSD)
    integer ATOM(NACLSD),NACLS,INDN0(NAEZ,naclsd),NUMN0(naez)

    integer J,LM1,LM2,M,N1,N2,IND1,IND2
    integer lmmax1, lmmax2
    integer ind

    ! ----------------------------------------------------------------------

    do M = 1,NACLS

      do N1 = 1,NUMN0(site_index)
        IND1 = INDN0(site_index,N1)
        if(ATOM(M).eq.IND1) then

          lmmax1 = kvstr(site_index + 1) - kvstr(site_index)
          lmmax2 = kvstr(IND1 + 1) - kvstr(IND1)

          do LM1 = 1,lmmax1
            do LM2 = 1,lmmax2

              ind = ka(ia(site_index) + N1 - 1) + (LM2 - 1) * lmmax1 + LM1 - 1

              smat(ind) = smat(ind) + EIKRM(M) * GINP(LM2,LM1,M)

            enddo
          enddo

        endif
      enddo

      J = ATOM(M)

      do N2 = 1,NUMN0(J)
        IND2 = INDN0(J,N2)
        if(site_index.eq.IND2) then

          lmmax1 = kvstr(J + 1) - kvstr(J)
          lmmax2 = kvstr(IND2 + 1) - kvstr(IND2)

          do LM2 = 1,lmmax2
            do LM1 = 1,lmmax1

              ind = ka(ia(J) + N2 - 1) + (LM2 - 1) * lmmax1 + LM1 - 1

              smat(ind) = smat(ind) + EIKRP(M) * GINP(LM1,LM2,M)

            enddo
          enddo

        endif
      enddo

    enddo

  end subroutine

end module
