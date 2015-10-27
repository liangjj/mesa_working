*deck icopy
      subroutine icopy(ia,ib,n)
      implicit integer (a-z)
      dimension ia(n), ib(n)
      call vec_$icopy(ia,ib,n)
      return
      end
