*deck copy
      subroutine copy(a,b,n)
      implicit integer (a-z)
      real*8 a,b
      dimension a(n), b(n)
      call vec_$dcopy(a,b,n)
      return
      end
