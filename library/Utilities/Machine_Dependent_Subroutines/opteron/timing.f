*deck @(#)timing.f	5.1  11/6/94
      subroutine timing(a,b,c)
      real*8 a,b,c
      real ta,tb,etime,time
      real t(2)
c
c     etime returns the elapsed cpu and system time. 
c     this routine is called at the beginning(lxopen) and end(chainx)
c     of each link and the code(tsumry) does its own differencing.
c     mesa passes real*8 numbers and etime returns real*4 precision
c     so a conversion(dble) is needed.
c
c      time=etime(ta,tb)
      time=etime(t)
c      a=dble(ta)
c      b=dble(tb)
      a=dble(t(1))
      b=dble(t(2))
      c=0.0d+00
c
c
      return
      end
