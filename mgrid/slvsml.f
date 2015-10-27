*deck slvsml 
c***begin prologue     slvsml
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            poisson equation on a coarse grid is solved
c***                   exactly using red-black gauss-seidel iteration.
c***description        the coarse grid finite difference equations
c***                   are solved exactly using gauss-seidel.  what
c***                   that means is that we iterate the gs cycle to
c***                   convergence.  we assume we can do that because the
c***                   grid is coarse enough to enable relaxation to
c***                   proceed smoothly.
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       rhoin
      subroutine slvsml(u,rhs,h,k,tol,n,maxit,type,prnt)
c
      implicit integer (a-z)
      logical prnt
      character*(*) type  
      real*8 u, rhs, h, k, hsq, hsqi, fac, rms, tol
      dimension u(n,n), rhs(n,n)
      common/io/inp, iout
      hsq=h*h
      hsqi=1.d0/hsq
      fac=4.d0-k*k*hsq
      if( prnt ) then
         write(iout,1)
         write(iout,2)
      endif   
      do 10 i=1,maxit
         call relax(u,rhs,h,k,n,type)
         call resid(u,rhs,rms,rms,h,k,n,'only-rms')
         if (prnt) then
             write(iout,3) i, rms
         endif
         if (rms.le.tol) go to 20
   10 continue
      write(iout,4) rms, tol   
      call lnkerr('no convergence in coarse solution')
   20 write(iout,5) i, rms
      return
    1 format(/,1x,'convergence information on coarse grid solution')
    2 format(/,15x,'iteration',13x,'rms error')
    3 format(18x,i3,11x,e15.8)   
    4 format(/,1x,'no convergence:  delta = ',e15.8,' tol = ',e15.8)
    5 format(/,1x,'converged in ',i3,' iterations',/1x,
     1            'rms weighted residual = ',e15.8)
      end





