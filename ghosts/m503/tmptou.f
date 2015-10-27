*deck @(#)tmptou.f	3.1  11/17/92
      subroutine tmptou(u,temp,ind,big,small)
      implicit integer (a-z)
      real*8 u(big,big), temp(small,small)
      integer ind(small)
c
c
      do 10 i=1,small
         do 20 j=1,small
            u(ind(i),ind(j))=temp(i,j)
   20    continue
   10 continue
c
c
      return
      end
