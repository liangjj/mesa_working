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
      subroutine resid(res,u,rhs,h,hsqi,fac,n,rms)
c
      implicit integer (a-z)
      logical prnt
      character*(*) dir  
      real*8 res, u, rhs, h, hsqi, fac, rms 
      dimension u(n,n), rhs(n,n), res(n,n)
      rms=0.d0
      do 30 j=2,n-1
         do 40 l=2,n-1
            res(j,l) =  rhs(j,l) - hsqi*( u(j+1,l) + u(j-1,l) +
     1                                    u(j,l+1) + u(j,l-1) ) 
     2                           + fac*u(j,l) 
            rms=rms+res(j,l)*res(j,l)
 40      continue
 30   continue   
      call rzero(res(1,1),n)
      call rzero(res(1,n),n)
      do 50 i=1,n
         res(1,i)=0.d0
         res(n,i)=0.d0
 50   continue   
      rms=h*h*rms
      rms=sqrt(rms)
      return
      end





