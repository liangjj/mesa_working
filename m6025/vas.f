      subroutine vas ( n, a, ia, s )
      implicit real *8 (a-h,o-z)
      dimension a(ia,n)
c
      if ( ia .lt. 1 ) return
      if ( n .lt. 1 ) return
c
      do 1 i = 1, n
         a(1,i) = a(1,i) + s
 1    continue
c
      return
      end
