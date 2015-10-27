*deck @(#)vsamul.f	5.1  11/6/94
      subroutine vsamul(a,b,c,xx,n)
      implicit real*8(a-h,o-z)
      real*8 a(n),b(n),c(n)
      real*8 xx
      integer n
c
      do 10 i=1,n
         a(i)=a(i)+b(i)*c(i)*xx
  10  continue
c
      return
      end
