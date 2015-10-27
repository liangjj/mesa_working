*deck @(#)timing.f	1.1  11/30/90
      subroutine timing(a,b,c)
      real*8 a,b,c
      real ta,tb,tc
c
c     secnds returns the elapsed cpu and system time. 
c     this routine is called at the beginning(lxopen) and end(chainx)
c     of each link and the code(tsumry) does its own differencing.
c     mesa passes real*8 numbers and etime returns real*4 precision
c     so a conversion(dble) is needed.
c
      ta=secnds(0.0)
      a=ta
      b=0.d0
      c=0.d0
c
c
      return
      end
