*deck @(#)mkxexp.f	1.1 9/8/91
      subroutine mkxexp(x,xinv,expfn,alpha,rmin,rdel,n)
      implicit integer (a-z)
      real *8 x, xinv, expfn, rmin, rdel, alpha, one
      dimension x(n), expfn(n), xinv(n)
      data one/1.d+00/
      x(1)=rmin
      xinv(1)=one/x(1)
      expfn(1)=exp(-alpha*rmin)
      do 10 i=2,n
         x(i)=x(i-1) + rdel
         xinv(i)=one/x(i)
         expfn(i)=exp(-alpha*x(i))
   10 continue
      return
      end
