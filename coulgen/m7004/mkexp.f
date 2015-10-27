*deck @(#)mkexp.f	1.1 9/8/91
      subroutine mkexp(x,expfn,alpha,n)
      implicit integer (a-z)
      real *8 x, expfn, alpha
      dimension x(n), expfn(n)
      do 10 i=1,n
         expfn(i)=exp(-alpha*x(i))
   10 continue
      return
      end
