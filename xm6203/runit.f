*deck @(#)runit.f	1.2  10/27/94
      subroutine runit(a,n)
c
      real*8 a(n,n)
c
      call rzero(a,n**2)
      do 1 i=1,n
         a(i,i)=1.0d+00
    1 continue
c
c
      return
      end
