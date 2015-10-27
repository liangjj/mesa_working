*deck @(#)mkx.f	1.1 9/8/91
      subroutine mkx(x,rmin,rdel,n)
      implicit integer (a-z)
      real *8 x, rmin, rdel
      dimension x(n)
      x(1)=rmin
      do 10 i=2,n
         x(i)=x(i-1) + rdel
   10 continue
      return
      end
