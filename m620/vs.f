*deck @(#)vs.f	5.1  11/6/94
      subroutine vs(a,b,c)
      implicit none
      integer i
      real*8 a(3),b(3),c(3)
c
c
      do 1 i=1,3
         c(i)=a(i)-b(i)
    1 continue
c
c
      return
      end
