*deck @(#)norm.f	5.1  11/6/94
      subroutine norm(x,y)
      implicit none
c
      real*8 x(3),y(3)
c
      integer i
      integer inp,iout
      real*8 zero,half,tenth,xm
c
      common/io/inp,iout
c
      parameter (zero=0.0d+00,half=0.5d+00,tenth=0.1d+00)      
c
 1000 format(' m620: norm of vector is less than 0.1')
c
      xm=zero
      do 10 i=1,3
         xm=xm+x(i)*x(i)
   10 continue
      xm=sqrt(xm)
      if (xm.lt.tenth) then
         write (iout,1000)
      endif
      do 20 i=1,3
         y(i)=x(i)/xm
   20 continue
c
c
      return
      end
