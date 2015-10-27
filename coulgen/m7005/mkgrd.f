*deck @(#)mkgrd.f	1.1 9/8/91
      subroutine mkgrd(x,xinv,xs,rmin,rdel,n,toskp)
      implicit integer (a-z)
      real *8 x, xinv, xs, rmin, rdel, one
      dimension x(n), xinv(n), xs(*)
      common /io/ inp, iout
      data one/1.d+00/
      x(1)=rmin
      xinv(1)=one/x(1)
      do 10 i=2,n
         x(i)=x(i-1) + rdel
         xinv(i)=one/x(i)
   10 continue
      nspln=0
      do 20 i=1,n,toskp
         nspln=nspln+1
         xs(nspln)=x(i)
   20 continue
      return
      end
