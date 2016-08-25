      SUBROUTINE WRITETBKKR_DATA(ie,lmaxd, lmmaxd, lmax, nspd, nspo,
     +                         nspoh, nrd, nembd, nemb, nclsd, ncls,
     +                         natypd, natyp, naezd, naez, naclsd,
     +                         ielast, ins, ALAT, BRAVAIS, RECBV,
     +                         RBASIS, CLS, NACLS, EQINV, EZOA, ATOM,
     +                         KAOEZ, RCLS, RR, energy, TMATLL, GINP)
      IMPLICIT NONE     

c-----------------------------------------------------------------
c write down parameters for external program PKKR, B.Zimmermann, Juelich, 05.2013
c-----------------------------------------------------------------

      INTEGER :: ie,lmaxd, lmmaxd, lmax, nspd, nspo, nspoh,
     +           nrd, nembd, nemb, nclsd, ncls, natypd, natyp,
     +           naezd, naez, naclsd, ielast, ins
      DOUBLE PRECISION :: ALAT,
     +                    BRAVAIS(3,3),
     +                    RECBV(3,3),
     +                    RBASIS(3,NAEZD+NEMBD)
      INTEGER :: CLS(NATYPD),
     +           NACLS(NCLSD),
     +           EQINV(NAEZD),
     +           EZOA(NACLSD,NAEZD),
     +           ATOM(NACLSD,NAEZD),
     +           KAOEZ(NAEZD+NEMBD)
      DOUBLE PRECISION :: RCLS(3,NACLSD,NCLSD),
     +                    RR(3,0:NRD)
      DOUBLE COMPLEX :: energy,
     +                  GINP(NACLSD*LMMAXD,LMMAXD,NCLSD),
     +                  TMATLL(NSPD*LMMAXD,NSPD*LMMAXD,NAEZD)

c...locals...
      INTEGER :: I1, I2, J, LM1, LM2, ICLS
      WRITE(6,*) 'in WRITETBKKR_DATA'

      if(ie==1)then
        open(unit=11,file='TBkkr_params.txt',
     +       form='formatted', action='write')
              write(11,'(2I8)') lmaxd, lmax
              write(11,'(3I8)') nspd, nspo, nspoh
              write(11,'(1I8)') nrd
              write(11,'(2I8)') nembd, nemb
              write(11,'(2I8)') nclsd, ncls
              write(11,'(2I8)') natypd, natyp
              write(11,'(2I8)') naezd, naez
              write(11,'(1I8)') naclsd
              write(11,'(1I8)') IELAST
              write(11,'(1I8)') INS
        close(11) 

         
        !open file for rest (double scalars + arrays)
        open(unit=12,file='TBkkr_container.txt',
     +           form='formatted',action='write')

              !write out lattice information
              write(12,'(A)') 'alat:'
              write(12,'(ES25.16)') ALAT
              write(12,'(A)') 'bravais:'
              write(12,'(3ES25.16)') ((BRAVAIS(I1,I2),I1=1,3),I2=1,3)
              write(12,'(A)') 'recbv:'
              write(12,'(3ES25.16)') ((RECBV(I1,I2),I1=1,3),I2=1,3)
              write(12,'(A)') 'RBASIS:'
              write(12,'(3ES25.16)') ((RBASIS(J,I1), J=1,3),
     +                                  I1=1,NAEZD+NEMBD)

              !write out cluster information
              write(12,'(A)') 'CLS:'
              write(12,'(1I8)') (CLS(I1),I1=1,NATYPD)
              write(12,'(A)') 'NACLS:'
              write(12,'(1I8)') (NACLS(I1), I1=1,NCLSD)
              write(12,'(A)') 'RCLS:'
              do I2=1,NCLSD
               do I1=1,NACLSD
                 write(12,'(3ES25.16)') RCLS(:,I1,I2)
               end do
              end do
              write(12,'(A)') 'EQINV:'
              write(12,'(1I8)') (EQINV(I1),I1=1,NAEZD)
              write(12,'(A)') 'EZOA:'
              write(12,'(1I8)') ((EZOA(I1,I2),I1=1,NACLSD),I2=1,NAEZD)
              write(12,'(A)') 'ATOM:'
              write(12,'(1I8)') ((ATOM(I1,I2),I1=1,NACLSD),I2=1,NAEZD)
              write(12,'(A)') 'KAOEZ:'
              write(12,'(1I8)') (KAOEZ(I1),   I1=1,NAEZD+NEMBD)
              write(12,'(A)') 'RR:'
              do I1=0,NRD
                write(12,'(3ES25.16)') RR(:,I1)
              end do

      end if!ie==1

      !====================================!
      !=== write out energy information ===!

      write(12,'(A)') 'energy(ie):'
      write(12,'(2ES25.16)') energy

      !write out (energydependent) T-matrix
      write(12,'(A)') 'TMATLL(ie):'
      do I1=1,NAEZD
        do LM2=1,NSPD*LMMAXD
          do LM1=1,NSPD*LMMAXD
            write(12,'(2ES25.16)') TMATLL(LM1,LM2,I1)
          end do
        end do
      end do

      !write out (energydependent) Green's function
      write(12,'(A)') 'GINP(ie):'
      do ICLS=1,NCLSD
        do LM2=1,LMMAXD
          do LM1=1,NACLSD*LMMAXD
            write(12,'(2ES25.16)')  GINP(LM1,LM2,ICLS)
          end do
        end do
      end do

      if(ie==ielast)then
        close(12)
        write(*,*) "Containerfile written"

        if(ielast==1) THEN
          write(*,*) "File is !NOT! capable of calculating FERMI-VL"
        elseif(ielast==3) THEN
          write(*,*) "File is capable of also calculating FERMI-VL"
        else
          write(*,*) "IELAST is neither 1 nor 2: smthng went wrong!"
        end if!ielast==1

      end if!ie==ielast

      END
