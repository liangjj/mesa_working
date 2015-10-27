      subroutine v2fit1 ( x, a, n, m, nmax )
c
c     v2fit1 -- a(i,j) = x(i)**(j-1)
c
      implicit real ( a - h, o - z )
      dimension x(n), a(nmax,m)
      call vsets ( n, a, 1, 1e0 )
      do 2 j = 2, m
      call vpv ( n, a(1,j-1),1, x,1, a(1,j),1 )
 2    continue
      return
      end
