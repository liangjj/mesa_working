*deck jacobi.f 
c***begin prologue     jacobi
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           jacobi weights
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            weight functions and their first and second
c***                   derivatives
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       jacobi
      subroutine jacobi(wt,dwt,ddwt,x,alf,bet,n)
c
      implicit integer (a-z)
      real*8 x, wt, dwt, ddwt, alf, bet, inv, inv1, fac, fac1
      dimension x(n), wt(n), dwt(n), ddwt(n)
      common/io/inp, iout      
      do 10 i=1,n
         inv = 1.d0 / ( 1.d0 - x(i) )
         inv1 = 1.d0 / ( 1.d0 + x(i) )
         fac = ( 1.d0 - x(i) )**(.5*alf )
         fac1 = ( 1.d0 + x(i) )**(.5*bet )
         wt(i) = fac * fac1
         dwt(i) = - .5d0 * alf * fac * inv * fac1 + fac * .5d0 * bet *
     1               fac1 * inv1
         ddwt(i) = -.5d0 * alf * ( .5d0 * alf - 1.d0 ) * fac * fac1 * 
     1               inv * inv + .5d0 * bet * ( .5d0 * bet - 1.d0 ) *
     2               fac * fac1 *inv1 * inv1
 10   continue   
      return
      end















