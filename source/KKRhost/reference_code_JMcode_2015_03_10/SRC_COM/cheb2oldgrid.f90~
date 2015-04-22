subroutine cheb2oldgrid(nrmax,nrmaxnew,lmmaxpot,rmesh,ncheb, &
                        npan_tot,rpan_intervall,ipan_intervall, &
                        arrayin,arrayout,irmd)
implicit none
integer  :: ncheb,npan_tot,nrmax,nrmaxnew,lmmaxpot,irmd
double precision :: rmesh(irmd)
double precision :: rpan_intervall(0:npan_tot)
integer :: ipan_intervall(0:npan_tot)
double complex :: arrayin(nrmaxnew,lmmaxpot)
double complex :: arrayout(nrmax,lmmaxpot)
!local
integer :: in,ir,ir2,ir3,it,ilm
integer :: intsub(2,nrmax)
double precision :: Cinvmatrix(0:ncheb,0:ncheb)
double precision,allocatable :: CCmatrix(:,:)
double complex :: alphaparams(0:ncheb,lmmaxpot)
double precision  :: rmeshnorm(nrmax)

! divide the mesh into subintervals
intsub(1,:)=0
intsub(2,:)=-1
in=1

ir=1
it=1
do while (ir.LE.nrmax .and. in.LE.npan_tot)
  if (abs(rmesh(ir)-rpan_intervall(in)).LT.10d-14) then
    intsub(it,in)=ir
    intsub(2,in)=ir
    it=1
    ir=ir+1
    in=in+1
  elseif (((rmesh(ir)-rpan_intervall(in).LT.-1d-10) ) .and. &
          ((rmesh(ir)-rpan_intervall(in-1)).GT.-1d-10) ) then
  intsub(it,in)=ir
  intsub(2,in)=ir
  it=2
  ir=ir+1
  elseif ((rmesh(ir)-rpan_intervall(in-1)).LT.-1d-10)  then
    intsub(1,in)=0
    intsub(2,in)=-1
    it=1
    ir=ir+1
  else
    it=1
    in=in+1
  end if
end do
if (abs(rmesh(nrmax)- rpan_intervall(npan_tot))>1d-15) then
   write(*,*) rmesh(nrmax),rpan_intervall(npan_tot),rmesh(nrmax)-rpan_intervall(npan_tot)
  stop 'error with the new and old mesh'
else
  
end if


call getCinvmatrix(Ncheb,Cinvmatrix)

in=0
do while (in<=npan_tot)
  in=in+1
  if (intsub(2,in)<intsub(1,in)) cycle
  alphaparams=0.0D0
  do ilm=1,lmmaxpot
  do ir2=0,ncheb
    do ir=0,ncheb
      alphaparams(ir2,ilm)=alphaparams(ir2,ilm)+ Cinvmatrix(ir2,ir)* arrayin(ipan_intervall(in-1)+1+ir,ilm)
    end do
  end do
  end do

  do ir=intsub(1,in),intsub(2,in)
    ir2=ir+1-intsub(1,in)
    rmeshnorm(ir2)=(2*rmesh(ir)-(rpan_intervall(in)+rpan_intervall(in-1))) &
                                    /(rpan_intervall(in)-rpan_intervall(in-1))
    if (abs(rmeshnorm(ir2))>1.0D0 .and. abs(rmeshnorm(ir2))-1.0D0<10e-14) then
      rmeshnorm(ir2)=sign(1.0D0,rmeshnorm(ir2))
    end if
  end do

  allocate(CCmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

  call getCCmatrix(Ncheb,rmeshnorm(1:intsub(2,in)-intsub(1,in)+1),intsub(2,in)-intsub(1,in)+1,CCmatrix)
  do ilm=1,lmmaxpot
  do ir2=intsub(1,in),intsub(2,in)
    do ir=0,ncheb
      ir3=ir2-intsub(1,in)+1
      arrayout(ir2,ilm)=arrayout(ir2,ilm)+ CCmatrix(ir3,ir)* alphaparams(ir,ilm)
    end do
  end do
  end do !ilm=1,lmmaxpot
  deallocate(CCmatrix) !(CCmatrix(intsub(2,in)-intsub(1,in)+1,0:ncheb))

end do !in

end subroutine cheb2oldgrid

