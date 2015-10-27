*deck @(#)bscattr.f	5.1  11/6/94
      subroutine bscattr(n,a,index,b)
      real*8 a(*),b(n)
      dimension index(n)
c
      do 1 i=1,n
         a(index(i))=a(index(i))+b(i)
   1  continue
c
      return
      end
