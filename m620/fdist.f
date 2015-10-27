*deck @(#)fdist.f	5.1  11/6/94
      subroutine fdist(x,y,dist)
      implicit none
c
      real*8 x(3),y(3),dist
c
      integer i
      real*8 z(3)
c
      do 10 i=1,3
         z(i)=y(i)-x(i)
   10 continue
      dist=sqrt(z(1)**2+z(2)**2+z(3)**2)
c
c
      return
      end
