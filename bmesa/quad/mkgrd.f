*deck @(#)mkgrd.f	1.1 9/8/91
      subroutine mkgrd(x,rmin,rmax,rdel,n)
      implicit integer (a-z)
      real *8 x, rmin, rmax, rdel
      dimension x(n)
      common /io/ inp,iout
      rdel=(rmax-rmin)/(n-1)
      x(1)=rmin
      do 10 i=2,n
         x(i)=x(i-1) + rdel
   10 continue
      return
      end
