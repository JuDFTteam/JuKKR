module mod_rotgll

contains

  subroutine rotgll(gmatll, natomimp, ijtabsym, ijtabsh, dsymll, symunitary, igf, rc, crel, rrel, krel, lmmaxd, irec)
    ! **********************************************************************
    ! *                                                                    *
    ! *   it calculates all the elements of the Green Function of          *
    ! *   the impurity cluster using the GF calculated for the             *
    ! *   representative pairs.                                            *
    ! *   the representative pair and the symmetry operation D are given   *
    ! *   through the arrays IJTABSH and IJTABSYM set up in < SHELLGEN2K > *
    ! *                                                                    *
    ! *     _ _                                                            *
    ! *     n n'                      n n'                                 *
    ! *     m m'               T      m m'                                 *
    ! *    G    (E) = SUM    D     * G    (E) * D                          *
    ! *     L L'      L1 L2   L L1    L1L2      L2 L'                      *
    ! *                                                                    *
    ! *                   _                 _                              *
    ! *              n    n            n'   n'                             *
    ! *   where   D R  = R     and  D R  = R                               *
    ! *              m    m            m    m                              *
    ! *                                                                    *
    ! **********************************************************************
    use :: mod_mympi, only: myrank, master
    use :: mod_datatypes, only: dp, sp
    use :: mod_changerep
    use :: mod_cmatstr

    implicit none
    ! ..
    ! .. Parameter definitions
    complex (kind=dp) :: czero, cone
    parameter (czero=(0.0d0,0.0d0), cone=(1.d0,0.d0))
    ! ..
    ! .. Scalar arguments
    integer :: ngclus, lmmaxd, irec
    integer :: igf, krel, natomimp
    ! ..
    ! .. Array arguments
    integer :: ijtabsym(*), ijtabsh(*)
    complex (kind=dp) :: gmatll(lmmaxd, lmmaxd, *), dsymll(lmmaxd, lmmaxd, *)
    complex (kind=dp) :: crel(lmmaxd, lmmaxd), rc(lmmaxd, lmmaxd), rrel(lmmaxd, lmmaxd)
    logical :: symunitary(*)
    ! ..
    ! .. Local arrays
    complex (kind=dp), allocatable :: gll(:, :, :, :), tpg(:, :)
    complex (kind=dp), allocatable :: gclust(:)
    ! ..
    ! .. Local scalars
    integer :: ilin, iq, icall, ish, isym, jq
    integer :: lm1, lm2, nlin, ilm, jlm
    character (len=1) :: cnt
    character (len=4) :: str4i, str4j
    character (len=18) :: str18
    ! ..
    ! .. External Functions
    logical :: test, opt
    external :: test
    ! ..
    ! .. Data statement
    data icall/1/
    ! ..
    ! .. Save statement
    save :: icall
    ! ..
    allocate (gll(lmmaxd,lmmaxd,natomimp,natomimp), tpg(lmmaxd,lmmaxd), stat=lm1)
    if (lm1/=0) then
      write (6, 100) ' GLL/TPG'
      stop '           < ROTGLL > '
    end if
100 format (6x, 'ERROR: failed to allocate array(s) :', a, /)
    ! **********************************************************************
    if (icall==1) then
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      if (myrank==master) write (1337, '(79("="))')
      if (myrank==master) write (1337, '(6X,2A)') 'ROTGLL : Expand GF for all pairs by rotation', ' and write out (all E-points)'
      if (myrank==master) write (1337, '(79("="))')
      if (myrank==master) write (1337, *)
      ! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
    end if
    ! ***********************************************************************
    do iq = 1, natomimp
      do jq = 1, natomimp
        ! -----------------------------------------------------------------------
        ilin = (iq-1)*natomimp + jq
        ish = ijtabsh(ilin)
        isym = ijtabsym(ilin)
        ! -----------------------------------------------------------------------
        ! for REL CASE look if it is a unitary / ANTI - unitary rotation
        ! -----------------------------------------------------------------------
        cnt = 'N'
        if (.not. symunitary(isym)) cnt = 'T'

        call zgemm('C', cnt, lmmaxd, lmmaxd, lmmaxd, cone, dsymll(1,1,isym), lmmaxd, gmatll(1,1,ish), lmmaxd, czero, tpg, lmmaxd)

        call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, tpg, lmmaxd, dsymll(1,1,isym), lmmaxd, czero, gll(1,1,iq,jq), lmmaxd)
        ! -----------------------------------------------------------------------
      end do
    end do
    ! ***********************************************************************

    ! visualise Gij
    if (test('Gmatij  ')) then
      write (1337, '(/,4X,70("+"),/,4X,A,I4)') 'cluster G_ij matrices for i,j = 1,', natomimp

      do iq = 1, natomimp
        write (1337, '(/,8X,66("="))')
        do jq = 1, natomimp

          write (str4i, '(I4)') iq
          write (str4j, '(I4)') jq
          str18 = '   i =' // str4i(1:4) // ' j =' // str4j(1:4)
          if (krel==0) then
            call cmatstr(str18, 18, gll(1,1,iq,jq), lmmaxd, lmmaxd, 0, 0, 0, 1.0e-8_dp, 6)
          else
            call changerep(gll(1,1,iq,jq), 'REL>RLM', tpg, lmmaxd, lmmaxd, rc, crel, rrel, str18, 18)
          end if
          if (jq<natomimp) write (1337, '(/,9X,65("-"))')
        end do
        if (iq==natomimp) write (1337, '(/,8X,66("="),/)')
      end do
      write (1337, '(4X,70("+"))')
    end if

    ! ***********************************************************************

    ! --> output of GF

    ngclus = lmmaxd*natomimp
    if (igf/=0) then
      icall = icall + 1
      nlin = 0

      allocate (gclust(ngclus*ngclus), stat=lm1)
      if (lm1/=0) then
        write (6, 100) ' GCLUST'
        stop '           < ROTGLL > '
      end if
      do jq = 1, natomimp
        do lm2 = 1, lmmaxd
          jlm = (jq-1)*lmmaxd + lm2
          do iq = 1, natomimp
            do lm1 = 1, lmmaxd
              nlin = nlin + 1
              if (nlin>ngclus*ngclus) stop '<ROTGLL>: NLIN.GT.(NATOMIMP*LMMAXD)**2'
              gclust(nlin) = gll(lm1, lm2, iq, jq)
              ! test
              ! WRITE(214321,'(4i,2E)') LM1,LM2,IQ,JQ,GCLUST(NLIN)
              ! writeout of green_host for WRTGREEN option
              if (opt('WRTGREEN') .and. myrank==master) then
                ilm = (iq-1)*lmmaxd + lm1
                write (58, '((2I5),(2e17.9))') jlm, ilm, gll(lm1, lm2, iq, jq)
              end if
            end do
          end do
        end do
      end do




      if ((opt('KKRFLEX '))) then
#ifdef CPP_MPI
        irec = irec
#else
        irec = icall
#endif
        ! force single precision complex writeout to minimize file size etc.
        ! maybe this can be removed in the future
        write (888, rec=irec) cmplx(gclust, kind=sp)
        if ((opt('GPLAIN  '))) then
          write (8888, '(50000E25.16)') gclust
        end if
      end if

      ! ==== the following write-out has been disabled, because it was assumed to be     !no-green
      ! ====  obsolete with the implementation of the MPI-communicated arrays. If I am   !no-green
      ! ====  wrong and the write-out is needed in subsequent parts, construct a         !no-green
      ! ====  test-option around it so that it is only written out in this case.         !no-green
      ! IF ( .not. OPT('KKRFLEX ') ) THEN                                        !no-green
      ! WRITE(88,REC=ICALL) GCLUST                                             !no-green
      ! endif                                                                   !no-green

      deallocate (gclust)
    end if                         ! IGF/=0
    ! ***********************************************************************
    deallocate (gll, tpg)
    return
  end subroutine rotgll

end module mod_rotgll
