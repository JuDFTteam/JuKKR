 program test_line_order
   use mod_calconfs,    only: order_lines
   integer     :: k(7008)
   integer     :: n_k
   integer, allocatable   :: k_ord(:),band(:)

   n_k = 7008
   allocate(k_ord(n_k), band(n_k), STAT=ierr)
     if(ierr/=0) stop 'Problem allocating kpt2irr_ord and band_indices'

   open (unit = 2, file = "temp.kpts")
   read (2,'(10I8)') k 
   close (2)
   !k(1) = 7
   !k(2) = 2
   !k(3) = 5
   !k(4) = 7
   !k(5) = 3
   !k(6) = 4
   !k(7) = 6
   !k(8) = 8
   !k(9) = 4
   !k(10) = 6
   !k(11) = 1
   !k(12) = 5


   call order_lines(k, n_k, k_ord, band)

   WRITE(*,*) k_ord(:)
   WRITE(*,*) band(:)

 end program test_line_order






