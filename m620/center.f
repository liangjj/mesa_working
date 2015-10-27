*deck @(#)center.f	5.1  11/6/94
      subroutine center(x,y,z,zl)
      implicit none
c
      real*8 x(3),y(3),z(3),zl(3)
c
      integer i
      real*8 half
c
      parameter (half=0.5d+00)
c
c
      do 10 i=1,3
         zl(i)=y(i)-x(i)
         z(i)=half*(x(i)+y(i))
   10 continue
      call norm(zl,zl)
c
c
      return
      end
