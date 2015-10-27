      subroutine normal(a,mdepth,m,n,z,istage)
      implicit real*8 ( a - h, o - z )
c based on max element
c irange contains exp of a matrix cols. 1-n
c istage is 1 to normalize a    is 2  to renormalize solution z
      common/irg/irang(30)
      dimension pivot(30),irange(30)
      equivalence(irang,irange)
      dimension a(mdepth,1),z(n)
      go to (1,5) istage
    1 do 2 k=1,n
      call mmax(a(1,k),m,pivot(k))
      call unpak (pivot(k),irange(k))
    2 continue
      do 4  k=1,n
      do 4  i=1,m
    4 a(i,k) = expad(a(i,k),-irange(k))
      go to 7
    5 do 6 k=1,n
    6 z(k)= expad(z(k),-irange(k))
    7 continue
      return
      end
