*deck resid
c***begin prologue     resid
c***date written       941208   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m1200, link 1200, multigrid
c***author             schneider, b. i.(nsf)
c***source             m1200
c***purpose            residual computation for poisson equation
c***description 
c***references        
c
c***routines called    iosys, util and mdutil
c***end prologue       resid
      subroutine resid(res,u,rhs,h,k,hsqi,fac,n,rms)
c
      implicit integer (a-z)
      logical prnt
      character*(*) dir  
      real*8 res, u, rhs, h, k, hsqi, fac, rms 
      dimension u(n,n,n), rhs(n,n,n), res(n,n,n)
      call rzero(res,n*n*n)
      rms=0.d0
      do 30 i=2,n-1
         do 40 j=2,n-1
            do 50 k=2,n-1
               res(i,j,k) =  rhs(i,j,k) - hsqi*( u(i+1,j,k) + 
     1                                           u(i-1,j,k) +
     2                                           u(i,j+1,k) + 
     3                                           u(i,j-1,k) +
     4                                           u(i,j,k+1) +
     5                                           u(i,j,k-1) ) 
     2                                  + fac*u(i,j,k) 
               rms=rms+res(i,j,k)*res(i,j,k)
   50       continue            
 40      continue
 30   continue   
      rms=h*h*h*rms
      rms=sqrt(rms)
      return
      end





