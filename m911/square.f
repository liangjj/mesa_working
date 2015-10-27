*deck @(#)square.f	5.1  11/6/94
      subroutine square(ca,c,n)
c
c
c
      implicit real*8 (a-h,o-z)
c
cvax  extended dummy ca,c
c
c
c     dimension ca(1),c(n,n)
c
c     ij=0
c     c(1,1)=0.0
c     do 2 i=2,n
c     do 1 j=1,i-1
c     ij=ij+1
c     t=ca(ij)
c     c(i,j)=t
c     c(j,i)=-t
c   1 continue
c     c(i,i)=0.0
c   2 continue
c     return
c
      dimension ca(1),c(1)
c
      k=1
      jisv=1
      do 2 i=1,n
         ij=i
         ji=jisv
cdir$ ivdep
         do 1 j=1,i-1
            t=ca(k)
            c(ij)=t
            c(ji)=-t
            k=k+1
            ij=ij+n
            ji=ji+1
    1    continue
         c(ij)=0.0d+00
         jisv=jisv+n
    2 continue
c
c
      return
      end
