*deck @(#)blfold.f	1.1  11/30/90
      subroutine blfold(tr,sq,n)
      implicit integer(a-z)
      real*8 tr(*),sq(n,n)
c
      ix=0
c
      do 1 i=1,n
         do 2 j=1,i
            ix=ix+1
            tr(ix)=sq(i,j)+sq(j,i)
    2    continue
    1 continue
c
      return
      end
