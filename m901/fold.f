*deck @(#)fold.f	5.1  11/6/94
      subroutine fold (sa,s,n)
c
c
c
      implicit real*8 (a-h,o-z)
c
cvax  extended dummy sa,s
c
c
c     dimension sa(1),s(n,n)
c
c     ij=0
c     do 2 i=2,n
c     do 1 j=1,i-1
c     ij=ij+1
c     sa(ij)=sa(ij)+s(i,j)-s(j,i)
c   1 continue
c   2 continue
c     return
c
      dimension sa(1),s(1)
c
      k=1
      jisv=1
      do 2 i=1,n
      ji=jisv
      ij=i
      do 1 j=1,i-1
      sa(k)=sa(k)+s(ij)-s(ji)
      k=k+1
      ij=ij+n
      ji=ji+1
    1 continue
      jisv=jisv+n
    2 continue
      return
      end
