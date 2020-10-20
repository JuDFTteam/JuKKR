module mod_broyden

  use mod_datatypes, only: dp

  private
  public :: broyden

contains
  
  !-------------------------------------------------------------------------------
  !> Summary: Broyden mixing
  !> Author:
  !> Category: KKRimp, potential
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine broyden (vector, vlen, alpha, rms, iter, n_init, mbroylen, mvlen)
    !
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     History: Original code written by D.D. Johnson (see PRB 38, 12807)
    !                  note:  there are a few typos in that paper but 
    !                  the code is working!
    !              Rewritten by W. A. Shelton for LSMS code 6/21/94
    !                  this version is easy to read (no goto!!!! more comments ...)
    !                  and is setup for MPP machines (not tested)
    !              Rewritten by T. C. Schulthess, ORNL, March 97
    !                  this version should work for any code (see comments below)
    !
    !     Bug fixes:   TCS, 8/5/97 see comments below 
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    !     further comments on how to use this subroutine:
    !     (if nothing useful stands here I had no time yet to update these
    !     comments, please consult usage in lkkr code version 0.6.3 or later,
    !     or call Thomas Schulthess (423) 5768067)
    !
    !     vector(r,i) -> i=1: old vector (input), scratch (ouput)
    !                 -> i=2: new vector (input), mixed vector (output)
    !     vlen        -> length of vector
    !     alpha       -> linear mixing factor
    !     rms         -> RMS difference between old and new vector
    !     iter        -> iteration number (if 1, linear mixing, broyden reset)
    !     broylen     -> number of iterations that are used for mixing (<=mbroylen)
    !     u, vt, f, df, vold, and w are working arrays that need to be saved
    !                   between call to this subroutine
    !     a, b, d, cm, and ipiv are working arrays that need not be saved
    !     mbroylen    -> maximum number of iterations that can be saved
    !     mvlen       -> maximum length of vectors
    !
    !     See declaration for exact dimentsions and types
    !
    !     There are two options for matrix inversions, a Gaussian
    !     elimination routine called invert1 and calls to lapack routines
    !     with pivoting (see comments "using invert1" and "using lapack").
    !     Obviously only one can be used, comment out the other one.
    !
    !     When using this subroutine in a parallel code in which only parts
    !     of the vectors are known on every node, make sure that the calls
    !     to gldsum (global sum) are correct (LKKR and LSMS codes have
    !     different calls).
    !
    !     In a serial code, either comment out the calls to glbsum or
    !     provide a dummy subroutine
    !
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
    implicit none 
    !
    integer :: mbroylen, mvlen ! used for dimensioning 
    real (kind=dp) :: vector(mvlen,2)
    integer :: vlen
    real (kind=dp) :: alpha
    real (kind=dp) :: rms
    integer :: iter
    integer :: n_init,broylen
    real (kind=dp), allocatable, save:: u(:,:)
    real (kind=dp), allocatable, save:: vt(:,:)
    real (kind=dp), allocatable, save:: f(:)
    real (kind=dp), allocatable, save:: df(:)
    real (kind=dp), allocatable, save:: vold(:)
    real (kind=dp) :: a(mbroylen,mbroylen)
    real (kind=dp) :: b(mbroylen,mbroylen)
    real (kind=dp) :: d(mbroylen,mbroylen)
    real (kind=dp) :: cm(mbroylen)
    real (kind=dp), allocatable, save:: w(:)
    integer :: ipiv(mbroylen) 
    !
    integer :: i,j,k,info,ntasks
    real (kind=dp) :: fac1,fac2,fnorm,dfnorm,w0,work,aij,gmi,cmj,wtmp 
    !
    integer :: lastit,lastm1,nn
    real (kind=dp), parameter :: zero=0.d0, one=1.d0
    real (kind=dp) :: amix
    save :: lastit,amix 
    
    broylen=mbroylen
    
    if (.not. allocated(u)) then
    allocate(u(mvlen,mbroylen))
    allocate(vt(mvlen,mbroylen))
    allocate(f(mvlen))
    allocate(df(mvlen))
    allocate(vold(mvlen))
    allocate(w(mbroylen))
    end if
    
    !
    if (iter<=n_init) then
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     first n_init iteration: preform linear mixing, load f and vold, set
    !                      different pointers and variables
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    lastit = iter-1          ! initialize pointers
    lastm1 = lastit-1 
    
    amix = alpha           ! for safety reasons
    
    do k = 1,vlen
        f(k) = vector(k,2) - vector(k,1)
        vold(k) = vector(k,1)
    enddo 
    
    do k = 1,vlen          ! this is the linear mixing
        vector(k,2) = vector(k,1) + amix * f(k)
    enddo 
    
    else ! (iter<=n_init)
    
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !     iter > n_init: this is where the non-linear mixing is done
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    lastit = lastit+1      ! update pointers
    lastm1 = lastit-1 
    
    if( iter > broylen) then ! set current lenght of broyden cycle
        nn = broylen
    else
        nn = lastit         !lastm1
    endif 
    
    w0=.01d0               ! set weighting factor for the zeroth iteration 
    
    !---- find: f[i] := vector(2)[i] - vector(1)[i]
    do k = 1,vlen
        df(k) = vector(k,2) - vector(k,1) - f(k)
    enddo
  
    do k = 1,vlen
        f(k) = vector(k,2) - vector(k,1)
    enddo 
    
    !---- find: fnorm  := |f|
    dfnorm = zero
    fnorm = zero
    
    do k = 1,vlen
        dfnorm = dfnorm + df(k)*df(k)
        fnorm  = fnorm  + f(k)*f(k)
    enddo 
    
    dfnorm = sqrt( dfnorm )
    fnorm  = sqrt( fnorm ) 
    
    !---- set: vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|
    fac2 = one/dfnorm
    fac1 = amix*fac2
    
    do k = 1,vlen
        vector(k,2) = fac1*df(k) + fac2*(vector(k,1) - vold(k))
        vold(k) = vector(k,1)
        vector(k,1) = fac2*df(k)
    enddo 
    
    !---- store vector(1) and vector(2) in the stacks u and vt restpectively
    call broy_sav(u,vt,vector,iter-1,broylen,vlen,mvlen) 
    
    !---- calculate coefficient matrices, a(i,j), and sum cm(i) for corrections:
    do j=1,nn - 1          ! off diagonal elements of a(i,j)
        do i = j+1,nn
        aij = zero
        do k = 1,vlen
            aij = aij + vt(k,j)*vt(k,i)
        enddo 
        a(i,j) = aij
        a(j,i) = aij
        enddo
    enddo 
  
    do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
        aij = zero
        cmj = zero
        do k=1,vlen
        cmj = cmj + vt(k,i)*f(k)
        aij = aij + vt(k,i)*vt(k,i)
        enddo 
        a(i,i) = aij
        cm(i) = cmj
    enddo 
    
    !---- shift down weights in stack
    if(iter-1 > broylen)then
        do i=1,broylen-1
        w(i)=w(i+1)
        enddo
    endif
    
    wtmp = zero
    if( rms > 1.0d-09 ) wtmp=2.0*sqrt(0.010d+00/rms)
    if( wtmp < one )    wtmp=1.00d+00
    
    if(iter > broylen)then
        w(broylen)=wtmp
    else
        w(lastit)=wtmp       !w(lastm1)=wtmp
    endif 
    
    !---- now calculate the b-matrix:
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> uses lapack
    do i=1,nn
        do j=1,nn
        b(j,i)= a(j,i)*w(j)*w(i)
        enddo
        b(i,i)= w0**2 + a(i,i)*w(i)*w(i)
    enddo
    
    call dgetrf(nn,nn,b,mbroylen,ipiv,info)
    call dgetri(nn,b,mbroylen, ipiv, d, nn, info )
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< uses lapack
  
    !---- mix vectors: 
    do k=1,vlen
        vector(k,2)= vold(k) + amix * f(k)
    enddo
    do i=1,nn
        gmi = zero
        do j=1,nn
        gmi = gmi + cm(j)*b(j,i)*w(j)
        enddo
    
        do k=1,vlen
        vector(k,2) = vector(k,2) - gmi*u(k,i)*w(i)
        enddo
    enddo 
    
    endif  ! (iter<=n_init)
  
  end subroutine broyden
  

  !-------------------------------------------------------------------------------
  !> Summary: Broyden mixing
  !> Author:
  !> Category: KKRimp, potential
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !-------------------------------------------------------------------------------
  subroutine broy_sav(fins,fots,vector,itscf,istore,ivsiz,mivsiz) 
  
    integer :: itscf
    integer :: ivsiz
    integer :: mivsiz
    integer :: istore 
    real (kind=dp) ::  vector(mivsiz,2)
    real (kind=dp) ::  fins(mivsiz,istore)
    real (kind=dp) ::  fots(mivsiz,istore) 
    
    if( itscf <= istore ) then 
    
      !     ==================================================================
      !     Load the first istore iterations in increasing iteration count
      !     ==================================================================
      do i = 1,ivsiz
        fins(i,itscf) = vector(i,2)
      enddo 
      
      do i = 1,ivsiz
        fots(i,itscf) = vector(i,1)
      enddo 
      
    else ! itscf<=istore 
    
      !     ==================================================================
      !     Re-load so that the ordering is in increasing iteration count
      !     ==================================================================
      do j = 1,istore - 1 
      
        do i = 1,ivsiz
          fins(i,j) = fins(i,j+1)
        enddo 
        
        do i = 1,ivsiz
          fots(i,j) = fots(i,j+1)
        enddo 
      
      enddo 
      
      !     ==================================================================
      !     Load current charge densities in the last storage location
      !     ==================================================================
      do i = 1,ivsiz
        fins(i,istore) = vector(i,2)
      enddo 
      
      do i = 1,ivsiz
        fots(i,istore) = vector(i,1)
      enddo 
      
    endif ! itscf<=istore
  
  end subroutine broy_sav
  

end module mod_broyden
