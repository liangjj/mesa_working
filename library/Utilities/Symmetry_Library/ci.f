*deck @(#)ci.f	5.1  11/6/94
      subroutine ci(t,n,told)
c
      implicit integer (a-z)
c
      real*8 t(3,3,n),told(3,3,n)
c
c     apply the inversion operator to the representation matrices.
c
      do 3 k=1,n
         do 2 j=1,3
            do 1 i=1,3
               t(i,j,k)=-told(i,j,k)
    1       continue
    2    continue
    3 continue
c
c
      return
      end
