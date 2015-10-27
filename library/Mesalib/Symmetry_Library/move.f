*deck @(#)move.f	5.1  11/6/94
      subroutine move(maxap3,a,b,n)
      implicit real*8(a-h,o-z)
c
c     move n sets of coordinates from a to b.
c
      dimension a(maxap3,3), b(maxap3,3)
c
c     call rtrace(6hmove  ,1)
      do 100 i=1,n
         do 100 j=1,3
 100        b(i,j) = a(i,j)
            return
            end
