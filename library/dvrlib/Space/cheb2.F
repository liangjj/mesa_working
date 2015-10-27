*deck cheb2.f 
c***begin prologue     cheb2
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           chebyshev weights
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            weight functions and their first and second
c***                   derivatives
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       cheb2
      subroutine cheb2(wt,dwt,ddwt,x,n)
c
      implicit integer (a-z)
      real*8 x, wt, dwt, ddwt, inv, inv2
      dimension x(n), wt(n), dwt(n), ddwt(n)
      common/io/inp, iout      
      do 10 i=1,n
         inv = 1.d0 / ( 1.d0 - x(i) * x(i) )
         fac = ( 1.d0 - x(i) * x(i) ) ** .25d0 
         wt(i) = fac
         inv2 = inv * inv
         dwt(i) = - .5d0 * x(i) * fac * inv
         ddwt(i) = -.5d0 * fac * inv * x(i) *x(i) * fac *inv2
 10   continue   
      return
      end















