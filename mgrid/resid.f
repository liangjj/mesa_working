*deck resid
c***begin prologue     resid
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            residual computation for poisson equation
c***                   exactly using red-black gauss-seidel iteration.
c***description 
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       resid
      subroutine resid(u,rhs,res,rms,h,k,n,dir)
c
      implicit integer (a-z)
      logical prnt
      character*(*) dir  
      real*8 u, rhs, res, rms, h, k, hsq, hsqi, fac, r 
      dimension u(n,n), rhs(n,n), res(n,n)
      common/io/inp, iout
      hsq=h*h
      hsqi=1.d0/hsq
      fac=4.d0-k*k*hsq
      rms=0.d0
      if (dir.eq.'only-rms') then
          do 10 j=2,n-1
             do 20 l=2,n-1
                r=( rhs(j,l) + u(j+1,l) + u(j-1,l) + u(j,l+1)
     1                       + u(j,l-1) - fac*u(j,l) )
                rms=rms+r*r
   20        continue
   10     continue
      else
         do 30 j=2,n-1
            do 40 l=2,n-1
               res(j,l)=( rhs(j,l) + u(j+1,l) + u(j-1,l) + u(j,l+1)
     1                             + u(j,l-1) - fac*u(j,l) )
               rms=rms+res(j,l)*res(j,l)
 40         continue
 30      continue   
      endif
      rms=sqrt(h*rms)               
      return
      end





