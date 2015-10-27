*deck @(#)symdm.f	5.1  11/6/94
      subroutine symdm(g,n)
c
      real*8 g(n,n)
c
      do 10 i=1,n
         do 5 j=1,i
            g(i,j)=(g(i,j)+g(j,i))*0.5d+00
            g(j,i)=g(i,j)
 5       continue
 10   continue
c
c
      return
      end
