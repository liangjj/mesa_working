*deck @(#)fold.f	5.1  11/6/94
      subroutine fold(xout,xin,n)
      real*8 xout(*),xin(n,n)
c
      ix=0
      do 1 i=1,n
         do 2 j=1,i
            ix=ix+1
            xout(ix)=xin(i,j)+xin(j,i)
  2      continue
         xout(ix)=xout(ix)*.5d0
  1   continue
c
      return
      end
