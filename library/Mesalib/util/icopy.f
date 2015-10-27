*deck icopy
      subroutine icopy(ia,ib,n)
      implicit integer (a-z)
      dimension ia(n), ib(n)
      do 10 i=1,n
   10 ib(i)=ia(i)
      return
      end
