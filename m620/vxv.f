*deck @(#)vxv.f	5.1  11/6/94
      subroutine vxv(a,b,c)
      implicit none
      real*8 a(3),b(3),c(3)
c
c
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
c
c
      return
      end
