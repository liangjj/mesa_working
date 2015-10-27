*deck @(#)fulint.f	1.1 9/8/91
      subroutine fulint(fun,wt,ints,n,mmax,type)
      implicit integer (a-z)
      real *8 fun, wt, ints, sum
      character*(*) type
      dimension fun(n,0:mmax,2), wt(n,n-1), ints(0:mmax,2)
      common/io/inp,iout
      do 10 i=1,n
         sum=0.d0
         do 20 j=1,n-1
            sum=sum+wt(i,j)
   20    continue
         wt(i,1)=sum
   10 continue
      if (type.eq.'exponential') then
          do 30 i=1,n
             ints(0,1)=ints(0,1)+wt(i,1)*fun(i,0,1)
   30     continue
      else
          do 40 m=0,mmax
             do 50 i=1,n
                ints(m,1)=ints(m,1)+fun(i,m,1)*wt(i,1)
                ints(m,2)=ints(m,2)+fun(i,m,2)*wt(i,1)
   50        continue
   40     continue
      endif               
      return
      end


