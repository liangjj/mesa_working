      subroutine mkx(x,xinv,rmin,rdel,n)
      implicit integer (a-z)
      real *8 x, xinv, rmin, rdel, one
      dimension x(n), xinv(n)
      data one / 1.0d+00 /
      do 10 i=1,n
         x(i)=rmin+(i-1)*rdel
         xinv(i)=one/x(i)
   10 continue
      return
      end
