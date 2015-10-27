*deck cpy3.f 
c***begin prologue     cpy3
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           tridiagonal matrix fill
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            fill a tridiagonal matrix
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       cpy3
      subroutine cpy3(band,e,d,n)
c
      implicit integer (a-z)
      real*8 band, e, d
      dimension band(2,n), e(n), d(n)
      common/io/inp, iout      
      do 10 i=1,n
         e(i) = band(1,i)
         d(i) = band(2,i)
 10   continue   
      return
      end
