*deck %W%  %G%
      subroutine blprt(x,n)
      implicit integer(a-z)
      real*8 x(*)
c
      ix=0
      do 10 i=1,n
c1       write(6,20)(x(ix+j),j=1,i)
  20     format(/,20(5(1x,f12.8)))
         ix=ix+i
  10  continue
c
      return
      end
