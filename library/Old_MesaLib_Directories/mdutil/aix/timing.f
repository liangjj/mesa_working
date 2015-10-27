*deck @(#)timing.f	5.1  11/6/94
      subroutine timing(a,b,c)
      real*8 second
      real*8 a,b,c
c
c     mclock returns the elapsed user time. in units of 1/60 sec. 
c     this routine is called at the beginning(lxopen) and end(chainx)
c     of each link and the code(tsumry) does its own differencing.
c     mesa passes real*8 numbers and mclock returns integer precision
c     so a conversion is needed.
c
      a=second(a)
      b=0.0d+00
      c=0.0d+00
c
c
      return
      end
