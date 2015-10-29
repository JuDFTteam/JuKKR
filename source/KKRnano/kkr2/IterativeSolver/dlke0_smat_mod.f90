module dlke0_smat_mod
  implicit none
  private
  public :: dlke0_smat ! todo: merge into module kkrmat01_new

  contains

  subroutine dlke0_smat(site_index, smat, ia, ka, kvstr, eikrm, eikrp, nacls, atom, numn0, indn0, Ginp, naez, lmmaxd, naclsd)
    integer, intent(in) :: site_index
    double complex, intent(inout) :: smat(:)
    integer, intent(in) :: ia(:)
    integer, intent(in) :: ka(:)
    integer, intent(in) :: kvstr(:)
    double complex, intent(in) :: eikrm(naclsd), eikrp(naclsd)
    integer, intent(in) :: nacls, atom(naclsd) 
    integer, intent(in) :: numn0(naez)
    integer, intent(in) :: indn0(naez,naclsd)
    integer, intent(in) :: naez, lmmaxd, naclsd
    
    double complex, intent(in) :: Ginp(lmmaxd,lmmaxd,naclsd)
    
    integer :: j, lm1, lm2, m, n1, n2, ind1, ind2, lmmax1, lmmax2, ind

    do m = 1, nacls

      do n1 = 1, numn0(site_index)
        ind1 = indn0(site_index,n1)
        if (atom(m) == ind1 .and. atom(m) > 0) then

          lmmax1 = kvstr(site_index + 1) - kvstr(site_index)
          lmmax2 = kvstr(ind1 + 1) - kvstr(ind1)

          do lm1 = 1, lmmax1
            do lm2 = 1, lmmax2

              ind = ka(ia(site_index) + n1-1) + lmmax1*(lm2-1) + (lm1-1)

              smat(ind) = smat(ind) + eikrm(m) * Ginp(lm2,lm1,m)

            enddo ! lm2
          enddo ! lm1

        endif
      enddo ! n1

      j = atom(m)
      if (j < 1) cycle

      do n2 = 1, numn0(j)
        ind2 = indn0(j,n2)
        if (site_index == ind2) then

          lmmax1 = kvstr(j + 1) - kvstr(j)
          lmmax2 = kvstr(ind2 + 1) - kvstr(ind2)

          do lm2 = 1, lmmax2
            do lm1 = 1, lmmax1

              ind = ka(ia(j) + n2 - 1) + (lm2 - 1) * lmmax1 + lm1 - 1

              smat(ind) = smat(ind) + eikrp(m) * Ginp(lm1,lm2,m)

            enddo ! lm1
          enddo ! lm2

        endif
      enddo ! n2

    enddo ! m

  endsubroutine ! dlke0_smat

endmodule ! dlke0_smat_mod
