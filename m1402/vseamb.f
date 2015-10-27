*deck @(#)vseamb.f	5.1  11/6/94
      subroutine vseamb(a,b,n,m)
      implicit real*8(a-h,o-z)
      dimension a(m,m),b(m,m)
c
c     dimension a(n),b(n)
c     do 10 i=1,n
c     a(i)=a(i)-b(i)
c 10  continue
c
       do 20 i=1,m
          do 30 j=1,m
           a(i,j)=a(i,j)-b(i,j)
  30      continue
  20   continue
      return
      end
