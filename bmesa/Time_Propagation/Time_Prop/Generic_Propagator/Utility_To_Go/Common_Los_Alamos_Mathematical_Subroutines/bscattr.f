*deck scatter   @(#)clams.f     3.2   11/9/87
      subroutine bscattr(n,a,index,b)
c
      implicit integer (a-z)
c
      real*8 a(*),b(n)
      integer index(n)
c
      do 1 i=1,n
         a(index(i))=b(i)
    1 continue
c
      return
      end
