*deck @(#)yfold.f	1.1  11/30/90
      subroutine yfold(ymat,scr,nob)
      implicit integer(a-z)
      real*8 ymat(nob,nob),scr(nob,nob)
c
      do 10 i=1,nob
         do 20 j=1,nob
            scr(i,j)=ymat(i,j)+ymat(j,i)
 20      continue
 10   continue
c
      do 30 i=1,nob
         do 40 j=1,nob
            ymat(j,i)=scr(j,i)*.5d+00
 40      continue
 30   continue
c
      return
      end
