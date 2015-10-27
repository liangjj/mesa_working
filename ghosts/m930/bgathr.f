*deck  @(#)bgathr.f	4.1 7/7/93
      subroutine bgathr(n,a,b,index)
      real*8 a(n),b(*)
      dimension index(n)
c
      do 1 i=1,n
         a(i)=b(index(i))
   1  continue
c
      return
      end
