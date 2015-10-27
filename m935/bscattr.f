*deck @(#)bscattr.f	5.1  11/6/94
      subroutine bscattr(n,xb,index,xm)
      real*8 xm(n),xb(*)
      integer index(n)
c
      do 1 i=1,n
         xb(index(i))=xb(index(i))+xm(i)
  1   continue
c
      return
      end
