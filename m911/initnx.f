*deck @(#)initnx.f	5.1  11/6/94
      subroutine initnx(m,n)
      implicit integer(a-z)
      dimension m(n)
c
c
      ix=0
      do 10 i=1,n
         m(i)=ix
         ix=ix+n
   10 continue
c
c
      return
      end
