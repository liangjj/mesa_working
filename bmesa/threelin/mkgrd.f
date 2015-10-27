c $Header: mkgrd.f,v 1.2 92/12/12 09:34:46 bis Exp $
*deck mkgrd.f
      subroutine mkgrd(x,rmin,rmax,rdel,n)
      implicit integer (a-z)
      real *8 x, rmin, rmax, rdel
      dimension x(0:n)
      common /io/ inp, iout
      x(0)=rmin
      do 10 i=1,n
         x(i)=x(i-1) + rdel
   10 continue
      return
      end
