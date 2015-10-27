*deck lagure.f 
c***begin prologue     lagure
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           laguree weights
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            weight functions and their first and second
c***                   derivatives
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       lagure
      subroutine lagure(wt,dwt,ddwt,x,alf,n)
c
      implicit integer (a-z)
      real*8 x, wt, dwt, ddwt, alf, fac, inv, inv2, exfac
      dimension x(n), wt(n), dwt(n), ddwt(n)
      common/io/inp, iout      
      do 10 i=1,n
         fac = x(i)**(.5d0*alf)
         inv = 1.d0 / x(i)
         inv2 = inv * inv 
         exfac = exp(-.5d0*x(i))
         wt(i) = fac*exfac
         dwt(i) = .5d0 * alf * fac * inv * exfac -
     1            .5d0 * fac * exfac   
         ddwt(i) = ( .5d0 * alf * ( .5d0 * alf - 1.d0 ) * fac * inv2
     1             - .5d0 * alf * fac * inv + .25d0 * fac ) * exfac
 10   continue   
      return
      end















