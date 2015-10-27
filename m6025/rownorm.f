      subroutine rownorm(a,mdepth,b,m,n)
      implicit real*8 ( a - h, o - z )
c max normalization
      dimension a(mdepth,1),b(1),c(60)
      common /ab/ pivot( 300),irange( 600)
      do 2 i=1,m
      do 1 k=1,n
    1 c(k)=a(i,k)
      c(n+1) = b(i)
      call mmax(c,n+1,pivot(i))
      call unpak(pivot(i),irange(i))
    2 continue
      do 4 i=1,m
      b(i)=expad(b(i),-irange(i))
      do 3 k=1,n
    3 a(i,k)=expad(a(i,k),-irange(i))
    4 continue
      return
      end
