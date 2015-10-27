*deck @(#)mkfun.f	1.1 9/8/91
      subroutine mkfun(x,fun,n,mmax,type)
      implicit integer (a-z)
      real *8 x, fun
      character*(*) type
      dimension x(n), fun(n,0:mmax,2)
      common/io/inp,iout
      if (type.eq.'exponential') then
          do 10 i=1,n
             fun(i,0,1)=exp(-x(i))
   10     continue
      else
          do 20 m=0,mmax
             do 30 i=1,n
                fun(i,m,1)=sin(m*x(i))
                fun(i,m,2)=cos(m*x(i))
   30        continue
   20     continue
      endif               
      return
      end
